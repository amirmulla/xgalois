#!/usr/bin/env python3
"""Live differential oracle runner for xgalois v0.4.

Executes randomized and boundary differential tests over NDJSON against a
compiled `xgalois_oracle_cli` binary, using SageMath when invoked under
`sage -python` and exact pure-Python reference arithmetic under `python3`.

Usage:
  python3 test/oracle/sage/live_oracle.py \\
      --cli=/tmp/xgalois_oracle_cli --seed=2026 --trials=25 --suites=all
"""

from __future__ import annotations

import argparse
import math
import pathlib
import random
import sys
import time
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

_SAGE_DIR = pathlib.Path(__file__).resolve().parent
if str(_SAGE_DIR) not in sys.path:
  sys.path.insert(0, str(_SAGE_DIR))

from oracles import HAS_SAGE, SuiteResult, sage  # noqa: E402
from oracles.gs_oracle import run_gs_oracle  # noqa: E402
from oracles.irs_oracle import run_irs_oracle  # noqa: E402
from oracles.lrc_oracle import run_lrc_oracle  # noqa: E402
from transport import OracleClient, OracleRpcError  # noqa: E402


# ============================================================================
# Live Algebra (Groups & Rings) Oracle Suite
# ============================================================================


def _extended_gcd_int(a: int, b: int) -> Tuple[int, int, int]:
  old_r, r = a, b
  old_s, s = 1, 0
  old_t, t = 0, 1
  while r != 0:
    q = old_r // r
    old_r, r = r, old_r - q * r
    old_s, s = s, old_s - q * s
    old_t, t = t, old_t - q * t
  if old_r < 0:
    old_r, old_s, old_t = -old_r, -old_s, -old_t
  return old_r, old_s, old_t


def _generalized_crt_int(
    residues: Sequence[int], moduli: Sequence[int]
) -> Optional[Tuple[int, int]]:
  cur_mod = abs(moduli[0])
  cur_rem = residues[0] % cur_mod
  for r_i, m_i_raw in zip(residues[1:], moduli[1:]):
    m_i = abs(m_i_raw)
    r_i %= m_i
    g, s, _ = _extended_gcd_int(cur_mod, m_i)
    diff = r_i - cur_rem
    if diff % g != 0:
      return None
    k = diff // g
    m_div_g = m_i // g
    new_mod = cur_mod * m_div_g
    step = (k * s) % m_div_g
    cur_rem = (cur_rem + cur_mod * step) % new_mod
    cur_mod = new_mod
  return cur_rem, cur_mod


def _ref_negacyclic_mul(
    a: Sequence[int], b: Sequence[int], q: int, n: int
) -> List[int]:
  if HAS_SAGE:
    r_poly = sage.PolynomialRing(sage.Zmod(q), "x")
    x_var = r_poly.gen()
    q_ring = r_poly.quotient(x_var**n + 1)
    prod_q = q_ring(list(a)) * q_ring(list(b))
    return [int(prod_q.lift()[i]) % q for i in range(n)]
  prod = [0] * n
  for i in range(n):
    for j in range(n):
      s = i + j
      if s < n:
        prod[s] = (prod[s] + a[i] * b[j]) % q
      else:
        prod[s - n] = (prod[s - n] - a[i] * b[j]) % q
  return prod


def run_algebra_oracle(
    client: OracleClient, *, seed: int = 2026, trials: int = 25
) -> SuiteResult:
  """Verifies Symmetric/Alternating/Dihedral groups, NegacyclicRing, and CRT live."""
  rng = random.Random(seed + 404)
  result = SuiteResult(name="algebra")

  # 1. SymmetricGroup orders & solvability
  for deg in (3, 4, 5, 6, 7, 12, 20):
    resp = client.call(
        "groups_and_rings_check", {"sub_op": "symmetric_group", "degree": deg}
    )
    if HAS_SAGE:
      g_sage = sage.SymmetricGroup(deg)
      exp_order = str(int(g_sage.order()))
      exp_solvable = bool(g_sage.is_solvable()) if deg <= 7 else None
    else:
      exp_order = str(math.factorial(deg))
      exp_solvable = (deg <= 4) if deg <= 7 else None

    if resp["order"] != exp_order or resp["bsgs_order"] != exp_order:
      result.record_fail(
          f"sym_{deg}: order={resp['order']} bsgs={resp['bsgs_order']}"
          f" != {exp_order}"
      )
    elif exp_solvable is not None and resp.get("is_solvable") != exp_solvable:
      result.record_fail(
          f"sym_{deg}: is_solvable={resp.get('is_solvable')} != {exp_solvable}"
      )
    else:
      result.record_pass()

  # 2. AlternatingGroup orders & solvability
  for deg in (3, 4, 5, 6, 7, 15):
    resp = client.call(
        "groups_and_rings_check", {"sub_op": "alternating_group", "degree": deg}
    )
    if HAS_SAGE:
      g_sage = sage.AlternatingGroup(deg)
      exp_order = str(int(g_sage.order()))
      exp_solvable = bool(g_sage.is_solvable()) if deg <= 7 else None
    else:
      exp_order = str(math.factorial(deg) // 2)
      exp_solvable = (deg <= 4) if deg <= 7 else None

    if resp["order"] != exp_order or resp["bsgs_order"] != exp_order:
      result.record_fail(
          f"alt_{deg}: order={resp['order']} bsgs={resp['bsgs_order']}"
          f" != {exp_order}"
      )
    elif exp_solvable is not None and resp.get("is_solvable") != exp_solvable:
      result.record_fail(
          f"alt_{deg}: is_solvable={resp.get('is_solvable')} != {exp_solvable}"
      )
    else:
      result.record_pass()

  # 3. DihedralGroup orders & centers
  for n in (3, 4, 5, 6, 8, 9, 12, 16):
    resp = client.call(
        "groups_and_rings_check", {"sub_op": "dihedral_group", "n": n}
    )
    if HAS_SAGE:
      d_sage = sage.DihedralGroup(n)
      exp_order = str(int(d_sage.order()))
      exp_center = int(d_sage.center().order())
    else:
      exp_order = str(2 * n)
      exp_center = 2 if (n % 2 == 0) else 1

    if resp["order"] != exp_order or resp["center_order"] != exp_center:
      result.record_fail(
          f"dih_{n}: order={resp['order']} center={resp['center_order']}"
          f" != ({exp_order}, {exp_center})"
      )
    else:
      result.record_pass()

  # 4. Randomized NegacyclicRing Z_q[x]/(x^n + 1) multiplications
  neg_trials = max(6, trials // 2)
  for idx in range(neg_trials):
    q = rng.choice([17, 97, 257, 3329, 12289])
    n = rng.choice([4, 8, 16])
    a = [rng.randrange(q) for _ in range(n)]
    b = [rng.randrange(q) for _ in range(n)]
    expected = _ref_negacyclic_mul(a, b, q, n)
    resp = client.call(
        "groups_and_rings_check",
        {
            "sub_op": "negacyclic_mul",
            "modulus_q": q,
            "degree": n,
            "a": a,
            "b": b,
        },
    )
    if resp["product"] != expected:
      result.record_fail(f"negacyclic_q{q}_n{n}#t{idx}: product mismatch")
    else:
      result.record_pass()

  # 5. Generalized CRT over Z (coprime, non-coprime consistent, inconsistent)
  crt_trials = max(6, trials // 2)
  for idx in range(crt_trials):
    kind = idx % 3
    if kind == 0:
      moduli = rng.sample([3, 5, 7, 11, 13, 17, 19], 3)
      residues = [rng.randrange(m) for m in moduli]
    elif kind == 1:
      # Consistent non-coprime moduli: sample target x, take x mod m_i
      moduli = [6, 10, 14, 15]
      x_secret = rng.randrange(210)
      residues = [x_secret % m for m in moduli]
    else:
      # Deliberately inconsistent non-coprime pair
      moduli = [6, 10, 15]
      residues = [1, 2, 0]  # 1 mod 2 != 2 mod 2

    ref = _generalized_crt_int(residues, moduli)
    if HAS_SAGE and ref is not None:
      assert int(sage.CRT_list(residues, moduli)) == ref[0]
      assert int(sage.lcm(moduli)) == ref[1]

    resp = client.call(
        "groups_and_rings_check",
        {"sub_op": "crt_int", "residues": residues, "moduli": moduli},
    )
    if ref is None:
      if resp["consistent"]:
        result.record_fail(
            f"crt_inconsistent#t{idx}: expected consistent=False"
        )
      else:
        result.record_pass()
    else:
      if (
          not resp["consistent"]
          or resp.get("solution") != str(ref[0])
          or resp.get("lcm") != str(ref[1])
      ):
        result.record_fail(
            f"crt_consistent#t{idx}: got {resp} != expected ({ref[0]}, {ref[1]})"
        )
      else:
        result.record_pass()

  return result


# ============================================================================
# CLI Driver
# ============================================================================

SUITE_REGISTRY: Dict[str, Callable[..., SuiteResult]] = {
    "lrc": run_lrc_oracle,
    "irs": run_irs_oracle,
    "gs": run_gs_oracle,
    "algebra": run_algebra_oracle,
}


def parse_suites(raw: str) -> List[str]:
  raw = raw.strip().lower()
  if raw == "all":
    return list(SUITE_REGISTRY.keys())
  chosen = [part.strip() for part in raw.split(",") if part.strip()]
  if not chosen:
    raise ValueError("--suites cannot be empty")
  for name in chosen:
    if name not in SUITE_REGISTRY:
      raise ValueError(
          f"Unknown suite '{name}'. Valid options: all,"
          f" {', '.join(SUITE_REGISTRY.keys())}"
      )
  return chosen


def main(argv: Optional[Sequence[str]] = None) -> int:
  parser = argparse.ArgumentParser(
      description="Run live differential Sage/Python oracle suites against xgalois_oracle_cli."
  )
  parser.add_argument(
      "--cli",
      required=True,
      type=pathlib.Path,
      help="Path to compiled xgalois_oracle_cli executable.",
  )
  parser.add_argument(
      "--seed",
      type=int,
      default=2026,
      help="Deterministic PRNG seed (default: 2026).",
  )
  parser.add_argument(
      "--trials",
      type=int,
      default=50,
      help="Number of trials per suite (default: 50).",
  )
  parser.add_argument(
      "--suites",
      type=str,
      default="all",
      help="Suites to run: all or comma-separated subset of lrc,irs,gs,algebra.",
  )
  parser.add_argument(
      "--timeout",
      type=float,
      default=30.0,
      help="Per-RPC timeout in seconds (default: 30.0).",
  )
  args = parser.parse_args(argv)

  try:
    selected_suites = parse_suites(args.suites)
  except ValueError as exc:
    print(f"ERROR: {exc}", file=sys.stderr)
    return 2

  backend = "sagemath" if HAS_SAGE else "python3_exact_reference"
  print(
      f"XGalois Live Differential Oracle | backend={backend} | "
      f"seed={args.seed} | trials={args.trials} | suites={','.join(selected_suites)}"
  )

  total_passed = 0
  total_failed = 0
  suite_results: List[SuiteResult] = []
  t0 = time.monotonic()

  try:
    with OracleClient(args.cli, timeout=args.timeout) as client:
      pong = client.ping()
      if pong.get("status") != "ok":
        print(f"ERROR: ping failed: {pong}", file=sys.stderr)
        return 1
      print(f"Connected to xgalois_oracle_cli v{pong.get('version', '?')}")

      for suite_name in selected_suites:
        st0 = time.monotonic()
        runner = SUITE_REGISTRY[suite_name]
        sr = runner(client, seed=args.seed, trials=args.trials)
        elapsed = time.monotonic() - st0
        suite_results.append(sr)
        total_passed += sr.passed
        total_failed += sr.failed
        status_str = "PASS" if sr.ok else "FAIL"
        print(
            f"  [{status_str}] {sr.name:<8s}: "
            f"{sr.passed} passed, {sr.failed} failed ({elapsed:.2f}s)"
        )
        for fail_msg in sr.failures[:10]:
          print(f"    - {fail_msg}", file=sys.stderr)
        if len(sr.failures) > 10:
          print(
              f"    ... and {len(sr.failures) - 10} more failures",
              file=sys.stderr,
          )
  except (FileNotFoundError, OracleRpcError) as exc:
    print(f"FATAL ORACLE ERROR: {exc}", file=sys.stderr)
    return 1

  total_elapsed = time.monotonic() - t0
  all_ok = total_failed == 0 and all(sr.ok for sr in suite_results)
  overall = "PASS" if all_ok else "FAIL"
  print(
      f"Summary [{overall}]: {total_passed} passed, {total_failed} failed "
      f"across {len(suite_results)} suite(s) in {total_elapsed:.2f}s"
  )
  return 0 if all_ok else 1


if __name__ == "__main__":
  sys.exit(main())
