#!/usr/bin/env python3
"""Live differential oracle for Guruswami-Sudan list decoding of GRS codes.

Verifies list completeness, soundness, and exact Johnson radius against
SageMath `GRSGuruswamiSudanDecoder` (when running under `sage -python`) and
exact brute-force Hamming-ball enumeration on small codes over GF(8), GF(11),
and GF(16).
"""

from __future__ import annotations

import itertools
import math
import random
from typing import Any, Dict, List, Sequence, Set, Tuple

from . import (
    FiniteFieldRef,
    GF8_MOD,
    GF16_MOD,
    HAS_SAGE,
    SuiteResult,
    sage,
)
from .irs_oracle import reference_grs_encode


def johnson_radius(n: int, k: int) -> int:
  return n - 1 - math.isqrt(n * (k - 1))


def hamming_distance(a: Sequence[int], b: Sequence[int]) -> int:
  return sum(1 for x, y in zip(a, b) if x != y)


def _enumerate_all_codewords(
    ff: FiniteFieldRef,
    eval_points: Sequence[int],
    col_mults: Sequence[int],
    k: int,
) -> List[Tuple[int, ...]]:
  """Enumerates every codeword of GRS(eval_points, col_mults, k) (for q^k <= 4096)."""
  n = len(eval_points)
  # Precompute generator rows: row i is v_j * alpha_j^i
  gen_rows: List[List[int]] = []
  pows = list(col_mults)
  for _ in range(k):
    gen_rows.append(list(pows))
    pows = [ff.mul(pows[j], eval_points[j]) for j in range(n)]

  codewords: List[Tuple[int, ...]] = []
  for msg in itertools.product(range(ff.q), repeat=k):
    cw = [0] * n
    for i, coeff in enumerate(msg):
      if coeff == 0:
        continue
      row = gen_rows[i]
      for j in range(n):
        cw[j] = ff.add(cw[j], ff.mul(coeff, row[j]))
    codewords.append(tuple(cw))
  return codewords


def _sage_gs_decode(
    ff: FiniteFieldRef,
    eval_points: Sequence[int],
    col_mults: Sequence[int],
    k: int,
    received: Sequence[int],
    tau: int,
) -> Set[Tuple[int, ...]]:
  """Runs SageMath GRSGuruswamiSudanDecoder at radius `tau`."""
  f_sage = ff.sage_field
  s_evals = [ff.to_sage(a) for a in eval_points]
  s_mults = [ff.to_sage(v) for v in col_mults]
  code = sage.codes.GeneralizedReedSolomonCode(s_evals, k, s_mults)
  dec = code.decoder("GuruswamiSudan", tau=tau)
  rx_vec = sage.vector(f_sage, [ff.to_sage(v) for v in received])
  decoded_list = dec.decode_to_code(rx_vec)
  return {tuple(ff.from_sage(x) for x in cw) for cw in decoded_list}


def _make_two_codeword_midpoint(
    ff: FiniteFieldRef,
    eval_points: Sequence[int],
    col_mults: Sequence[int],
    k: int,
    tau: int,
    rng: random.Random,
) -> Tuple[List[int], Tuple[int, ...], Tuple[int, ...]]:
  """Constructs a received word r with TWO distinct codewords c1, c2 within distance <= tau.

  Requires d = n - k + 1 <= 2 * tau. Pick two messages whose difference vanishes
  at k - 1 evaluation points so d_H(c1, c2) = n - (k - 1) = d, then set r to
  agree with c1 on ceil(d / 2) differing positions and with c2 on floor(d / 2).
  """
  n = len(eval_points)
  d = n - k + 1
  assert d <= 2 * tau

  # Build nonzero difference polynomial delta(x) = scale * prod_{z in Z}(x - z)
  # with |Z| = k - 1 so delta has degree k - 1 and vanishes on Z.
  zero_positions = rng.sample(range(n), k - 1)
  poly = [ff.random_nonzero(rng)]
  for z_idx in zero_positions:
    z = eval_points[z_idx]
    neg_z = ff.neg(z)
    nxt = [0] * (len(poly) + 1)
    for i, c in enumerate(poly):
      nxt[i] = ff.add(nxt[i], ff.mul(c, neg_z))
      nxt[i + 1] = ff.add(nxt[i + 1], c)
    poly = nxt
  assert len(poly) == k

  m1 = [ff.random_element(rng) for _ in range(k)]
  m2 = [ff.add(a, b) for a, b in zip(m1, poly)]
  c1 = reference_grs_encode(ff, eval_points, col_mults, m1)
  c2 = reference_grs_encode(ff, eval_points, col_mults, m2)
  diff_positions = [j for j in range(n) if c1[j] != c2[j]]
  assert len(diff_positions) == d

  # Split differing positions between c1 and c2 so d_H(r, c1) <= tau and d_H(r, c2) <= tau
  take_from_c2 = d // 2
  rx = list(c1)
  for pos in diff_positions[:take_from_c2]:
    rx[pos] = c2[pos]
  assert hamming_distance(rx, c1) <= tau and hamming_distance(rx, c2) <= tau
  return rx, tuple(c1), tuple(c2)


def _check_gs_response(
    ff: FiniteFieldRef,
    eval_points: Sequence[int],
    col_mults: Sequence[int],
    k: int,
    received: Sequence[int],
    tau: int,
    resp: Dict[str, Any],
    all_codewords: Sequence[Tuple[int, ...]] | None,
    must_contain: Sequence[Tuple[int, ...]],
    result: SuiteResult,
    tag: str,
) -> None:
  n = len(eval_points)
  exp_j = johnson_radius(n, k)
  if resp["johnson_radius"] != exp_j:
    result.record_fail(
        f"{tag}: johnson_radius={resp['johnson_radius']} != {exp_j}"
    )
    return
  actual_radius = resp["decoding_radius"]
  if actual_radius < tau:
    result.record_fail(
        f"{tag}: decoding_radius={actual_radius} < requested tau={tau}"
    )
    return

  got_list = [tuple(cw) for cw in resp["codewords"]]
  got_set = set(got_list)
  if len(got_set) != len(got_list):
    result.record_fail(f"{tag}: duplicate codewords in GS output list")
    return

  # Soundness: every returned codeword must be a valid GRS codeword within actual_radius
  for cw, msg in zip(got_list, resp["messages"]):
    re_enc = tuple(reference_grs_encode(ff, eval_points, col_mults, msg))
    if re_enc != cw:
      result.record_fail(f"{tag}: returned codeword failed GRS re-encoding check")
      return
    if hamming_distance(cw, received) > actual_radius:
      result.record_fail(
          f"{tag}: returned codeword distance {hamming_distance(cw, received)}"
          f" > actual_radius={actual_radius}"
      )
      return

  for req_cw in must_contain:
    if req_cw not in got_set:
      result.record_fail(
          f"{tag}: expected codeword missing from GS list (dist="
          f"{hamming_distance(req_cw, received)}, radius={actual_radius})"
      )
      return

  # Exact completeness via brute-force Hamming-ball enumeration on small codes
  if all_codewords is not None:
    expected_ball = {
        cw
        for cw in all_codewords
        if hamming_distance(cw, received) <= actual_radius
    }
    if got_set != expected_ball:
      missing = len(expected_ball - got_set)
      extra = len(got_set - expected_ball)
      result.record_fail(
          f"{tag}: Hamming-ball mismatch (missing={missing}, extra={extra},"
          f" actual_radius={actual_radius})"
      )
      return

  # Cross-check against SageMath GRSGuruswamiSudanDecoder when available
  if HAS_SAGE:
    sage_set = _sage_gs_decode(
        ff, eval_points, col_mults, k, received, actual_radius
    )
    if got_set != sage_set:
      result.record_fail(
          f"{tag}: mismatch vs Sage GRSGuruswamiSudanDecoder "
          f"(xg={len(got_set)}, sage={len(sage_set)})"
      )
      return

  result.record_pass()


def run_gs_oracle(
    client: Any, *, seed: int = 2026, trials: int = 25
) -> SuiteResult:
  """Runs the full live differential Guruswami-Sudan list-decoding oracle suite."""
  rng = random.Random(seed + 303)
  result = SuiteResult(name="gs")

  f8 = FiniteFieldRef(2, 3, GF8_MOD)
  f11 = FiniteFieldRef(11, 1)
  f16 = FiniteFieldRef(2, 4, GF16_MOD)

  eval7 = f8.primitive_eval_points(7)
  eval10 = f11.primitive_eval_points(10)
  eval15 = f16.primitive_eval_points(15)

  mults7 = [((i * 2) % 7) + 1 for i in range(7)]
  mults10 = [((i * 3) % 10) + 1 for i in range(10)]
  mults15 = [((i * 3) % 15) + 1 for i in range(15)]

  # Small codes with exact brute-force Hamming-ball enumeration:
  # - GRS(7, 2) over GF(8): d=6, t_u=2, tau_J=4 (64 codewords)
  # - GRS(7, 3) over GF(8): d=5, t_u=2, tau_J=3 (512 codewords)
  # - GRS(10, 2) over GF(11): d=9, t_u=4, tau_J=6 (121 codewords)
  # - GRS(15, 2) over GF(16): d=14, t_u=6, tau_J=11 (256 codewords)
  # - GRS(15, 3) over GF(16): d=13, t_u=6, tau_J=9 (4096 codewords)
  small_configs = [
      (f8, eval7, [1] * 7, 2, 4, "gs_f8_7_2_tau4"),
      (f8, eval7, mults7, 3, 3, "gs_f8_7_3_tau3"),
      (f11, eval10, mults10, 2, 6, "gs_f11_10_2_tau6"),
      (f16, eval15, [1] * 15, 2, 10, "gs_f16_15_2_tau10"),
      (f16, eval15, mults15, 3, 9, "gs_f16_15_3_tau9"),
  ]

  cached_tables: Dict[str, List[Tuple[int, ...]]] = {}
  for ff, evals, mults, k, _tau, name in small_configs:
    cached_tables[name] = _enumerate_all_codewords(ff, evals, mults, k)

  # 1. Exact Hamming-ball completeness & soundness on random & multi-codeword cases
  for idx in range(trials):
    ff, evals, mults, k, tau, name = small_configs[idx % len(small_configs)]
    all_cws = cached_tables[name]
    n = len(evals)
    d = n - k + 1

    # Alternate between: (a) deliberate 2-codeword list midpoint, (b) random tau errors,
    # and (c) uniformly random received word (tests empty/multi-codeword balls!)
    kind = idx % 3
    must_contain: List[Tuple[int, ...]] = []
    if kind == 0 and d <= 2 * tau:
      rx, c1, c2 = _make_two_codeword_midpoint(ff, evals, mults, k, tau, rng)
      must_contain = [c1, c2]
      sub = "midpoint2"
    elif kind == 1:
      msg = [ff.random_element(rng) for _ in range(k)]
      cw = reference_grs_encode(ff, evals, mults, msg)
      rx = list(cw)
      for pos in rng.sample(range(n), tau):
        rx[pos] = ff.add(rx[pos], ff.random_nonzero(rng))
      must_contain = [tuple(cw)]
      sub = "randerr"
    else:
      rx = [ff.random_element(rng) for _ in range(n)]
      sub = "uniform"

    resp = client.call(
        "grs_guruswami_sudan",
        {
            "field": ff.spec(),
            "evaluation_points": evals,
            "column_multipliers": mults,
            "k": k,
            "received": rx,
            "tau": tau,
        },
    )
    _check_gs_response(
        ff,
        evals,
        mults,
        k,
        rx,
        tau,
        resp,
        all_cws,
        must_contain,
        result,
        f"{name}_{sub}#g{idx}",
    )

  # 2. Larger GRS(15, 5) over GF(16): d=11, t_u=5, tau_J=7 (test tau=6 and tau=7)
  large_trials = max(4, trials // 4)
  for idx in range(large_trials):
    tau = 6 if (idx % 2 == 0) else 7
    if idx % 2 == 0:
      rx, c1, c2 = _make_two_codeword_midpoint(
          f16, eval15, mults15, 5, tau, rng
      )
      must_contain = [c1, c2]
      sub = "midpoint2"
    else:
      msg = [f16.random_element(rng) for _ in range(5)]
      cw = reference_grs_encode(f16, eval15, mults15, msg)
      rx = list(cw)
      for pos in rng.sample(range(15), tau):
        rx[pos] = f16.add(rx[pos], f16.random_nonzero(rng))
      must_contain = [tuple(cw)]
      sub = "randerr"

    resp = client.call(
        "grs_guruswami_sudan",
        {
            "field": f16.spec(),
            "evaluation_points": eval15,
            "column_multipliers": mults15,
            "k": 5,
            "received": rx,
            "tau": tau,
        },
    )
    _check_gs_response(
        f16,
        eval15,
        mults15,
        5,
        rx,
        tau,
        resp,
        None,
        must_contain,
        result,
        f"gs_f16_15_5_tau{tau}_{sub}#L{idx}",
    )

  return result
