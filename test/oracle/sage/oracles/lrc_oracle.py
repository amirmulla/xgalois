#!/usr/bin/env python3
"""Live differential oracle for Locally Recoverable Codes (Tamo-Barg & Availability).

Verifies XGalois LocallyRecoverableCode against independent exact reference math
(SageMath when running under `sage -python`, exact Python finite-field math
under standard `python3`).
"""

from __future__ import annotations

import math
import random
from typing import Any, Dict, List, Sequence, Tuple

from . import (
    FiniteFieldRef,
    GF16_MOD,
    GF64_MOD,
    HAS_SAGE,
    SuiteResult,
    sage,
)


def pklk_bound(n: int, k: int, r: int, delta: int = 2) -> int:
  """Prakash-Kamath-Lalitha-Kumar Singleton-like upper bound for (r, delta)-LRC."""
  ceil_k_over_r = (k + r - 1) // r
  return n - k + 1 - (ceil_k_over_r - 1) * (delta - 1)


def _reference_tamo_barg_encode(
    ff: FiniteFieldRef,
    eval_points: Sequence[int],
    good_poly: Sequence[int],
    message: Sequence[int],
    r: int,
) -> List[int]:
  """Independently encodes `message` via Tamo-Barg evaluation polynomial f_a(x)."""
  k = len(message)
  q_k, r_k = divmod(k, r)
  if HAS_SAGE:
    f_sage = ff.sage_field
    pr = sage.PolynomialRing(f_sage, "X")
    x_var = pr.gen()
    g_poly = pr([ff.to_sage(c) for c in good_poly])
    f_poly = pr(0)
    idx = 0
    for j in range(q_k):
      g_pow = g_poly**j
      for i in range(r):
        f_poly += ff.to_sage(message[idx]) * g_pow * (x_var**i)
        idx += 1
    if r_k > 0:
      g_pow = g_poly**q_k
      for i in range(r_k):
        f_poly += ff.to_sage(message[idx]) * g_pow * (x_var**i)
        idx += 1
    return [ff.from_sage(f_poly(ff.to_sage(pt))) for pt in eval_points]

  out: List[int] = []
  for pt in eval_points:
    g_val = ff.eval_poly(good_poly, pt)
    acc = 0
    idx = 0
    g_pow = 1
    for _ in range(q_k):
      x_pow = 1
      for _ in range(r):
        term = ff.mul(message[idx], ff.mul(g_pow, x_pow))
        acc = ff.add(acc, term)
        x_pow = ff.mul(x_pow, pt)
        idx += 1
      g_pow = ff.mul(g_pow, g_val)
    if r_k > 0:
      x_pow = 1
      for _ in range(r_k):
        term = ff.mul(message[idx], ff.mul(g_pow, x_pow))
        acc = ff.add(acc, term)
        x_pow = ff.mul(x_pow, pt)
        idx += 1
    out.append(acc)
  return out


def _verify_single_set_case(
    client: Any,
    ff: FiniteFieldRef,
    n: int,
    k: int,
    r: int,
    delta: int,
    partition_type: str,
    rng: random.Random,
    result: SuiteResult,
    tag: str,
) -> None:
  s = r + delta - 1
  expected_d = pklk_bound(n, k, r, delta)
  message = [ff.random_element(rng) for _ in range(k)]

  # 1. Single-erasure local repair check
  erased_pos = rng.randrange(n)
  resp = client.call(
      "lrc_encode_and_repair",
      {
          "field": ff.spec(),
          "n": n,
          "k": k,
          "r": r,
          "delta": delta,
          "partition_type": partition_type,
          "message": message,
          "erased_positions": [erased_pos],
          "allow_global_fallback": False,
      },
  )

  eval_pts: List[int] = resp["evaluation_points"]
  good_poly: List[int] = resp["good_poly"]
  fibers: List[List[int]] = resp["fiber_partition"]
  repair_sets: List[List[int]] = resp["repair_sets"]
  codeword: List[int] = resp["codeword"]

  if len(eval_pts) != n or len(set(eval_pts)) != n:
    result.record_fail(f"{tag}: evaluation points not n distinct elements")
    return
  if len(good_poly) - 1 != s:
    result.record_fail(
        f"{tag}: good_poly degree {len(good_poly) - 1} != s={s}"
    )
    return

  # Verify good_poly is constant on every fiber and distinct across fibers
  fiber_vals = set()
  seen_coords = set()
  for f_idx, fiber in enumerate(fibers):
    if len(fiber) != s:
      result.record_fail(f"{tag}: fiber {f_idx} size {len(fiber)} != s={s}")
      return
    vals = {ff.eval_poly(good_poly, eval_pts[c]) for c in fiber}
    if len(vals) != 1:
      result.record_fail(f"{tag}: good_poly not constant on fiber {f_idx}")
      return
    fv = next(iter(vals))
    if fv in fiber_vals:
      result.record_fail(f"{tag}: duplicate good_poly fiber value across fibers")
      return
    fiber_vals.add(fv)
    seen_coords.update(fiber)
  if seen_coords != set(range(n)):
    result.record_fail(f"{tag}: fibers do not partition [0, {n})")
    return

  # Verify PKLK optimal distance
  if resp["minimum_distance"] != expected_d:
    result.record_fail(
        f"{tag}: minimum_distance={resp['minimum_distance']} != PKLK={expected_d}"
    )
    return

  # Independent codeword verification
  expected_cw = _reference_tamo_barg_encode(ff, eval_pts, good_poly, message, r)
  if codeword != expected_cw:
    result.record_fail(f"{tag}: encoded codeword mismatch vs reference")
    return
  if resp["unencoded_message"] != message:
    result.record_fail(f"{tag}: unencoded_message mismatch")
    return

  # Verify strict single-erasure local repair
  if not resp["repair_success"]:
    result.record_fail(f"{tag}: single-erasure local repair failed")
    return
  if resp["symbols_read"] != r:
    result.record_fail(
        f"{tag}: single-erasure symbols_read={resp['symbols_read']} != r={r}"
    )
    return
  if resp["repaired_codeword"] != codeword:
    result.record_fail(f"{tag}: repaired_codeword != original codeword")
    return

  # Independently verify Lagrange repair from first r indices of repair set
  rset = repair_sets[erased_pos]
  if erased_pos in rset or len(rset) != s - 1:
    result.record_fail(f"{tag}: invalid repair set for coord {erased_pos}")
    return
  helpers = rset[:r]
  ref_recovered = ff.lagrange_eval(
      [eval_pts[h] for h in helpers],
      [codeword[h] for h in helpers],
      eval_pts[erased_pos],
  )
  if ref_recovered != codeword[erased_pos]:
    result.record_fail(f"{tag}: independent Lagrange repair mismatch")
    return
  result.record_pass()

  # 2. If delta >= 3, test multi-erasure local-only repair of (delta - 1)
  #    erasures inside a single fiber (strictly r symbols read, no global fallback)
  if delta >= 3:
    fiber = fibers[rng.randrange(len(fibers))]
    erased_in_fiber = sorted(rng.sample(fiber, delta - 1))
    m_resp = client.call(
        "lrc_encode_and_repair",
        {
            "field": ff.spec(),
            "n": n,
            "k": k,
            "r": r,
            "delta": delta,
            "partition_type": partition_type,
            "message": message,
            "erased_positions": erased_in_fiber,
            "allow_global_fallback": False,
        },
    )
    if not m_resp["repair_success"]:
      result.record_fail(
          f"{tag}(multi-local): failed repairing {delta - 1} erasures in fiber"
      )
      return
    if m_resp["used_global_fallback"]:
      result.record_fail(f"{tag}(multi-local): unexpectedly used global fallback")
      return
    if m_resp["symbols_read"] != r:
      result.record_fail(
          f"{tag}(multi-local): symbols_read={m_resp['symbols_read']} != r={r}"
      )
      return
    if m_resp["repaired_codeword"] != codeword:
      result.record_fail(f"{tag}(multi-local): repaired codeword mismatch")
      return
    result.record_pass()


def _verify_availability_case(
    client: Any,
    ff: FiniteFieldRef,
    construction: str,
    params: Dict[str, Any],
    localities: Sequence[int],
    local_distances: Sequence[int],
    expected_designed_d: int,
    rng: random.Random,
    result: SuiteResult,
    tag: str,
) -> None:
  k = params["k_expected"]
  message = [ff.random_element(rng) for _ in range(k)]
  test_coord = rng.randrange(params["n_expected"])

  req: Dict[str, Any] = {
      "field": ff.spec(),
      "construction": construction,
      "localities": list(localities),
      "message": message,
      "test_coord": test_coord,
  }
  if local_distances:
    req["local_distances"] = list(local_distances)
  if construction == "multiplicative":
    req["n"] = params["n_expected"]
    req["k"] = k
  else:
    req["num_superblocks"] = params["num_superblocks"]
    req["message_superblocks"] = params["message_superblocks"]

  resp = client.call("lrc_availability_verify", req)

  n = resp["n"]
  avail = resp["availability"]
  if n != params["n_expected"] or resp["k"] != k or avail != len(localities):
    result.record_fail(
        f"{tag}: shape mismatch got (n={n}, k={resp['k']}, t={avail})"
    )
    return
  if resp["designed_distance"] != expected_designed_d:
    result.record_fail(
        f"{tag}: designed_distance={resp['designed_distance']} !="
        f" expected={expected_designed_d}"
    )
    return

  deltas = list(local_distances) if local_distances else [2] * len(localities)
  all_repair_sets: List[List[List[int]]] = resp["all_repair_sets"]
  if len(all_repair_sets) != n:
    result.record_fail(f"{tag}: all_repair_sets length {len(all_repair_sets)} != {n}")
    return

  # Verify pairwise disjointness of repair sets at EVERY coordinate
  for i in range(n):
    sets_i = all_repair_sets[i]
    if len(sets_i) != avail:
      result.record_fail(f"{tag}: coord {i} has {len(sets_i)} != {avail} sets")
      return
    seen_union = set()
    for u, rset in enumerate(sets_i):
      expected_rset_size = localities[u] + deltas[u] - 2
      if i in rset:
        result.record_fail(f"{tag}: coord {i} appears in its own repair set {u}")
        return
      if len(rset) != expected_rset_size:
        result.record_fail(
            f"{tag}: coord {i} set {u} size {len(rset)} != {expected_rset_size}"
        )
        return
      rset_set = set(rset)
      if len(rset_set) != len(rset):
        result.record_fail(f"{tag}: duplicate index in repair set {u} of coord {i}")
        return
      if seen_union & rset_set:
        result.record_fail(
            f"{tag}: repair sets for coord {i} are NOT pairwise disjoint!"
        )
        return
      seen_union |= rset_set

  # Verify independent recovery from EVERY disjoint repair set
  repaired_entries = resp["repaired_via_each_set"]
  if len(repaired_entries) != avail:
    result.record_fail(f"{tag}: repaired_via_each_set len != {avail}")
    return
  for u, entry in enumerate(repaired_entries):
    if (
        not entry["success"]
        or not entry["matches_original"]
        or entry["symbols_read"] != localities[u]
    ):
      result.record_fail(f"{tag}: independent repair failed on set {u}: {entry}")
      return
  if not resp.get("all_coords_all_sets_ok", False):
    result.record_fail(f"{tag}: exhaustive all-coords all-sets check failed")
    return

  # For multiplicative CRT availability, independently cross-check Lagrange
  # interpolation on every repair set at test_coord
  if construction == "multiplicative":
    eval_pts = resp["evaluation_points"]
    codeword = resp["codeword"]
    for u in range(avail):
      helpers = all_repair_sets[test_coord][u][: localities[u]]
      ref_val = ff.lagrange_eval(
          [eval_pts[h] for h in helpers],
          [codeword[h] for h in helpers],
          eval_pts[test_coord],
      )
      if ref_val != codeword[test_coord]:
        result.record_fail(
            f"{tag}: independent Lagrange check failed for set {u} at"
            f" coord {test_coord}"
        )
        return

  result.record_pass()


def run_lrc_oracle(
    client: Any, *, seed: int = 2026, trials: int = 25
) -> SuiteResult:
  """Runs the full live differential LRC oracle suite."""
  rng = random.Random(seed + 101)
  result = SuiteResult(name="lrc")

  f13 = FiniteFieldRef(13, 1)
  f16 = FiniteFieldRef(2, 4, GF16_MOD)
  f17 = FiniteFieldRef(17, 1)
  f64 = FiniteFieldRef(2, 6, GF64_MOD)

  # Single-set configurations covering multiplicative, additive, r|k, r!|k, and (r, delta)
  single_set_configs: List[Tuple[FiniteFieldRef, int, int, int, int, str, str]] = [
      (f13, 12, 6, 3, 2, "multiplicative", "mult_f13_12_6_r3"),
      (f13, 12, 5, 3, 2, "multiplicative", "mult_f13_12_5_r3_nondiv"),
      (f16, 15, 6, 2, 2, "multiplicative", "mult_f16_15_6_r2"),
      (f16, 15, 5, 2, 2, "multiplicative", "mult_f16_15_5_r2_nondiv"),
      (f17, 16, 8, 3, 2, "multiplicative", "mult_f17_16_8_r3_nondiv"),
      (f16, 16, 6, 3, 2, "additive", "add_f16_16_6_r3"),
      (f16, 16, 7, 3, 2, "additive", "add_f16_16_7_r3_nondiv"),
      (f16, 15, 6, 3, 3, "multiplicative", "rdelta_f16_15_6_r3_d3"),
      (f16, 15, 5, 3, 3, "multiplicative", "rdelta_f16_15_5_r3_d3_nondiv"),
      (f16, 16, 4, 2, 3, "additive", "rdelta_add_f16_16_4_r2_d3"),
      (f64, 63, 18, 6, 4, "multiplicative", "rdelta_f64_63_18_r6_d4"),
  ]

  for t_idx in range(trials):
    cfg = single_set_configs[t_idx % len(single_set_configs)]
    ff, n, k, r, delta, ptype, name = cfg
    _verify_single_set_case(
        client,
        ff,
        n,
        k,
        r,
        delta,
        ptype,
        rng,
        result,
        f"{name}#t{t_idx}",
    )

  # Ensure every single_set_config ran at least once even if trials < len(configs)
  if trials < len(single_set_configs):
    for cfg in single_set_configs[trials:]:
      ff, n, k, r, delta, ptype, name = cfg
      _verify_single_set_case(
          client, ff, n, k, r, delta, ptype, rng, result, f"{name}#extra"
      )

  # Availability configurations (CRT multiplicative & direct-sum additive)
  def crt_designed_d(
      n: int, k: int, locs: Sequence[int], deltas: Sequence[int]
  ) -> int:
    bs = [r + d - 1 for r, d in zip(locs, deltas)]
    exps = [
        e for e in range(n) if all((e % s_u) < r_u for s_u, r_u in zip(bs, locs))
    ]
    return n - exps[k - 1]

  avail_trials = max(4, trials // 3)
  avail_configs = [
      (
          f13,
          "multiplicative",
          {"n_expected": 12, "k_expected": 4},
          [2, 3],
          [2, 2],
          crt_designed_d(12, 4, [2, 3], [2, 2]),
          "avail_mult_f13_12_4_r2_3",
      ),
      (
          f16,
          "multiplicative",
          {"n_expected": 15, "k_expected": 6},
          [2, 4],
          [2, 2],
          crt_designed_d(15, 6, [2, 4], [2, 2]),
          "avail_mult_f16_15_6_r2_4",
      ),
      (
          f16,
          "additive",
          {
              "num_superblocks": 2,
              "message_superblocks": 1,
              "n_expected": 16,
              "k_expected": 3,
          },
          [1, 3],
          [2, 2],
          (2 - 1 + 1) * 2 * 2,
          "avail_add_f16_16_3_r1_3",
      ),
      (
          f64,
          "additive",
          {
              "num_superblocks": 2,
              "message_superblocks": 1,
              "n_expected": 32,
              "k_expected": 9,
          },
          [3, 3],
          [2, 2],
          (2 - 1 + 1) * 2 * 2,
          "avail_add_f64_32_9_r3_3",
      ),
  ]

  for idx in range(avail_trials):
    ff, constr, pdict, locs, deltas, exp_d, name = avail_configs[
        idx % len(avail_configs)
    ]
    _verify_availability_case(
        client,
        ff,
        constr,
        pdict,
        locs,
        deltas,
        exp_d,
        rng,
        result,
        f"{name}#a{idx}",
    )

  return result
