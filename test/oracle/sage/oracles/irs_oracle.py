#!/usr/bin/env python3
"""Live differential oracle for Interleaved RS (IRS) and Power Decoding.

Verifies homogeneous and heterogeneous InterleavedGRSCode (`ssb`, `ssb_exact`,
`bky_hankel`) strictly beyond unique decoding radius up to t_max against an
independent reference BKY Hankel solver (using SageMath matrices when available
and exact pure-Python RREF otherwise), plus PowerDecodingDecoder on RS(63, 10)
over GF(64) at t = 30, 31, 32.
"""

from __future__ import annotations

import random
from typing import Any, Dict, List, Optional, Sequence, Tuple

from . import (
    FiniteFieldRef,
    GF16_MOD,
    GF32_MOD,
    GF64_MOD,
    HAS_SAGE,
    SuiteResult,
    sage,
)


def compute_collaborative_radius(redundancies: Sequence[int]) -> int:
  """Exact max t >= 0 such that sum_{r_m > t} (r_m - t) >= t."""
  if not redundancies:
    return 0
  max_r = max(redundancies)
  best_t = 0
  for t in range(1, max_r + 1):
    eqs = sum(r - t for r in redundancies if r > t)
    if eqs >= t:
      best_t = t
    else:
      break
  return best_t


def compute_dual_multipliers(
    ff: FiniteFieldRef,
    eval_points: Sequence[int],
    col_multipliers: Sequence[int],
) -> List[int]:
  """Computes exact GRS dual multipliers w_j = (v_j * prod_{i != j}(a_j - a_i))^-1."""
  n = len(eval_points)
  duals: List[int] = []
  for j in range(n):
    prod = col_multipliers[j]
    aj = eval_points[j]
    for i in range(n):
      if i != j:
        prod = ff.mul(prod, ff.sub(aj, eval_points[i]))
    duals.append(ff.inv(prod))
  return duals


def reference_grs_encode(
    ff: FiniteFieldRef,
    eval_points: Sequence[int],
    col_multipliers: Sequence[int],
    message: Sequence[int],
) -> List[int]:
  if HAS_SAGE:
    f_sage = ff.sage_field
    s_evals = [ff.to_sage(a) for a in eval_points]
    s_mults = [ff.to_sage(v) for v in col_multipliers]
    code = sage.codes.GeneralizedReedSolomonCode(s_evals, len(message), s_mults)
    enc = code.encoder("EvaluationPolynomial")
    pr = sage.PolynomialRing(f_sage, "X")
    msg_poly = pr([ff.to_sage(c) for c in message])
    return [ff.from_sage(c) for c in enc.encode(msg_poly)]
  return [
      ff.mul(v, ff.eval_poly(message, a))
      for a, v in zip(eval_points, col_multipliers)
  ]


def reference_irs_syndromes(
    ff: FiniteFieldRef,
    eval_points: Sequence[int],
    dual_mults: Sequence[int],
    received_matrix: Sequence[Sequence[int]],
    redundancies: Sequence[int],
) -> List[List[int]]:
  """Computes syndrome sequences S^(m)_i = sum_j y_{j, m} * w_j * alpha_j^i."""
  n = len(eval_points)
  ell = len(redundancies)
  max_r = max(redundancies) if redundancies else 0
  # Precompute w_j * alpha_j^i
  h_rows: List[List[int]] = []
  pows = list(dual_mults)
  for _ in range(max_r):
    h_rows.append(list(pows))
    pows = [ff.mul(pows[j], eval_points[j]) for j in range(n)]

  seqs: List[List[int]] = []
  for m in range(ell):
    r_m = redundancies[m]
    seq_m: List[int] = []
    for i in range(r_m):
      acc = 0
      row = h_rows[i]
      for j in range(n):
        if received_matrix[j][m] != 0:
          acc = ff.add(acc, ff.mul(row[j], received_matrix[j][m]))
      seq_m.append(acc)
    seqs.append(seq_m)
  return seqs


def reference_bky_hankel_locator(
    ff: FiniteFieldRef,
    seqs: Sequence[Sequence[int]],
    max_errors: int,
) -> Optional[List[int]]:
  """Exact reference BKY stacked Hankel solver returning Lambda(z) = [1, L1, ..., L_tau]."""
  if all(all(v == 0 for v in s) for s in seqs):
    return [1]
  if max_errors == 0:
    return None

  for tau in range(max_errors, 0, -1):
    rows: List[List[int]] = []
    for seq_m in seqs:
      if len(seq_m) <= tau:
        continue
      for i in range(len(seq_m) - tau):
        row = [seq_m[i + tau - 1 - c] for c in range(tau)]
        row.append(ff.neg(seq_m[i + tau]))
        rows.append(row)
    if len(rows) < tau:
      continue

    if HAS_SAGE:
      f_sage = ff.sage_field
      mat = sage.matrix(
          f_sage,
          len(rows),
          tau + 1,
          [[ff.to_sage(v) for v in r] for r in rows],
      )
      rref = mat.rref()
      pivots = list(rref.pivots())
      if len(pivots) == tau and pivots[-1] < tau:
        coeffs = [1] + [ff.from_sage(rref[c, tau]) for c in range(tau)]
        if coeffs[-1] != 0:
          return coeffs
    else:
      rref, pivots = ff.rref_with_pivots(rows)
      if len(pivots) == tau and pivots[-1] < tau:
        coeffs = [1] + [rref[c][tau] for c in range(tau)]
        if coeffs[-1] != 0:
          return coeffs

  return None


def reference_bky_error_positions(
    ff: FiniteFieldRef,
    eval_points: Sequence[int],
    dual_mults: Sequence[int],
    received_matrix: Sequence[Sequence[int]],
    redundancies: Sequence[int],
) -> Optional[List[int]]:
  """Runs reference BKY Hankel solver and Chien search for non-zero eval points."""
  max_errors = compute_collaborative_radius(redundancies)
  seqs = reference_irs_syndromes(
      ff, eval_points, dual_mults, received_matrix, redundancies
  )
  lam = reference_bky_hankel_locator(ff, seqs, max_errors)
  if lam is None:
    return None
  deg = len(lam) - 1
  if deg == 0:
    return []
  positions: List[int] = []
  for j, a in enumerate(eval_points):
    if a == 0:
      continue
    if ff.eval_poly(lam, ff.inv(a)) == 0:
      positions.append(j)
  if len(positions) != deg:
    return None
  return positions


def _sample_full_rank_burst_instance(
    ff: FiniteFieldRef,
    eval_points: Sequence[int],
    col_mults: Sequence[int],
    dual_mults: Sequence[int],
    dims: Sequence[int],
    t_errors: int,
    rng: random.Random,
    max_attempts: int = 30,
) -> Tuple[List[List[int]], List[List[int]], List[List[int]], List[int]]:
  """Generates a random IRS codeword + burst error pattern verified uniquely solvable by ref BKY."""
  n = len(eval_points)
  ell = len(dims)
  redundancies = [n - k for k in dims]

  messages = [[ff.random_element(rng) for _ in range(k)] for k in dims]
  cols = [
      reference_grs_encode(ff, eval_points, col_mults, msg) for msg in messages
  ]
  cw_matrix = [[cols[m][j] for m in range(ell)] for j in range(n)]

  for _ in range(max_attempts):
    err_pos = sorted(rng.sample(range(n), t_errors))
    rx_matrix = [list(row) for row in cw_matrix]
    for pos in err_pos:
      # Ensure burst error row is non-zero
      row_err = [ff.random_element(rng) for _ in range(ell)]
      if all(v == 0 for v in row_err):
        row_err[rng.randrange(ell)] = ff.random_nonzero(rng)
      for m in range(ell):
        rx_matrix[pos][m] = ff.add(rx_matrix[pos][m], row_err[m])

    ref_pos = reference_bky_error_positions(
        ff, eval_points, dual_mults, rx_matrix, redundancies
    )
    if ref_pos == err_pos:
      return messages, cw_matrix, rx_matrix, err_pos

  raise RuntimeError(
      f"Failed to sample uniquely decodable IRS burst pattern after"
      f" {max_attempts} attempts"
  )


def _verify_irs_config(
    client: Any,
    ff: FiniteFieldRef,
    eval_points: Sequence[int],
    col_mults: Sequence[int],
    dims: Sequence[int],
    t_errors: int,
    rng: random.Random,
    result: SuiteResult,
    tag: str,
) -> None:
  n = len(eval_points)
  redundancies = [n - k for k in dims]
  expected_t_max = compute_collaborative_radius(redundancies)
  expected_t_u = min(redundancies) // 2
  assert t_errors <= expected_t_max

  dual_mults = compute_dual_multipliers(ff, eval_points, col_mults)
  messages, cw_matrix, rx_matrix, err_pos = _sample_full_rank_burst_instance(
      ff, eval_points, col_mults, dual_mults, dims, t_errors, rng
  )

  for mode in ("ssb", "ssb_exact", "bky_hankel"):
    resp = client.call(
        "irs_decode",
        {
            "field": ff.spec(),
            "evaluation_points": list(eval_points),
            "column_multipliers": list(col_mults),
            "dims": list(dims),
            "received_matrix": rx_matrix,
            "mode": mode,
            "fallback_to_independent": False,
        },
    )
    if resp["collaborative_radius"] != expected_t_max:
      result.record_fail(
          f"{tag}[{mode}]: collaborative_radius={resp['collaborative_radius']}"
          f" != {expected_t_max}"
      )
      return
    if resp["unique_radius"] != expected_t_u:
      result.record_fail(
          f"{tag}[{mode}]: unique_radius={resp['unique_radius']} !="
          f" {expected_t_u}"
      )
      return
    if not resp["success"]:
      result.record_fail(
          f"{tag}[{mode}]: failed on reference-solvable t={t_errors} pattern:"
          f" {resp.get('error_message', '')}"
      )
      return
    if resp["used_fallback"]:
      result.record_fail(
          f"{tag}[{mode}]: unexpectedly reported used_fallback=True"
      )
      return
    if resp["codeword_matrix"] != cw_matrix:
      result.record_fail(f"{tag}[{mode}]: decoded codeword_matrix mismatch")
      return
    if resp["burst_error_positions"] != err_pos:
      result.record_fail(
          f"{tag}[{mode}]: burst_error_positions={resp['burst_error_positions']}"
          f" != {err_pos}"
      )
      return
    if resp["messages"] != messages:
      result.record_fail(f"{tag}[{mode}]: decoded messages mismatch")
      return
    result.record_pass()


def _verify_power_decoding_rs63_10(
    client: Any,
    ff64: FiniteFieldRef,
    rng: random.Random,
    result: SuiteResult,
    trials_per_t: int,
) -> None:
  """Verifies PowerDecodingDecoder on low-rate RS(63, 10) over GF(64) at t = 30, 31, 32."""
  n, k = 63, 10
  eval_pts = ff64.primitive_eval_points(n)
  col_mults = [1] * n
  dual_mults = compute_dual_multipliers(ff64, eval_pts, col_mults)
  # Virtual IRS dims for optimal ell_p = 3: k_s = s*(k - 1) + 1 = [10, 19, 28]
  virtual_dims = [10, 19, 28]
  virtual_reds = [n - km for km in virtual_dims]
  expected_t_p = compute_collaborative_radius(virtual_reds)  # 33

  for t_err in (30, 31, 32):
    for trial in range(trials_per_t):
      msg = [ff64.random_element(rng) for _ in range(k)]
      cw = reference_grs_encode(ff64, eval_pts, col_mults, msg)

      # Sample an error pattern that the reference virtual-IRS BKY solver uniquely
      # resolves so we have deterministic ground truth at every t in {30, 31, 32}.
      rx: Optional[List[int]] = None
      for _ in range(25):
        err_pos = sorted(rng.sample(range(n), t_err))
        cand_rx = list(cw)
        for pos in err_pos:
          cand_rx[pos] = ff64.add(cand_rx[pos], ff64.random_nonzero(rng))
        # Build virtual powered matrix Y[j][m] = cand_rx[j]^(m+1)
        v_mat = [
            [y, ff64.mul(y, y), ff64.mul(y, ff64.mul(y, y))] for y in cand_rx
        ]
        ref_pos = reference_bky_error_positions(
            ff64, eval_pts, dual_mults, v_mat, virtual_reds
        )
        if ref_pos == err_pos:
          rx = cand_rx
          break
      if rx is None:
        result.record_fail(
            f"power_rs63_10_t{t_err}#t{trial}: failed to sample reference-solvable pattern"
        )
        continue

      resp = client.call(
          "grs_power_decode",
          {
              "field": ff64.spec(),
              "evaluation_points": eval_pts,
              "k": k,
              "received": rx,
              "mode": "ssb",
              "fallback_to_bm": False,
          },
      )
      if resp["powering_degree"] != 3:
        result.record_fail(
            f"power_rs63_10_t{t_err}#t{trial}: powering_degree={resp['powering_degree']} != 3"
        )
        continue
      if resp["decoding_radius"] != expected_t_p:
        result.record_fail(
            f"power_rs63_10_t{t_err}#t{trial}: decoding_radius={resp['decoding_radius']} != {expected_t_p}"
        )
        continue
      if not resp["success"]:
        result.record_fail(
            f"power_rs63_10_t{t_err}#t{trial}: PowerDecoding failed on"
            f" reference-solvable t={t_err}: {resp.get('error_message', '')}"
        )
        continue
      if resp["codeword"] != cw or resp["message"] != msg:
        result.record_fail(
            f"power_rs63_10_t{t_err}#t{trial}: decoded codeword/message mismatch"
        )
        continue
      result.record_pass()


def run_irs_oracle(
    client: Any, *, seed: int = 2026, trials: int = 25
) -> SuiteResult:
  """Runs the full live differential IRS + Power Decoding oracle suite."""
  rng = random.Random(seed + 202)
  result = SuiteResult(name="irs")

  f16 = FiniteFieldRef(2, 4, GF16_MOD)
  f32 = FiniteFieldRef(2, 5, GF32_MOD)
  f64 = FiniteFieldRef(2, 6, GF64_MOD)
  f17 = FiniteFieldRef(17, 1)

  eval15 = f16.primitive_eval_points(15)
  eval16_f17 = f17.primitive_eval_points(16)
  eval31 = f32.primitive_eval_points(31)
  eval63 = f64.primitive_eval_points(63)

  ones15 = [1] * 15
  mults15 = [((i * 3) % 15) + 1 for i in range(15)]
  ones16 = [1] * 16
  ones31 = [1] * 31
  mults31 = [((i * 5) % 31) + 1 for i in range(31)]
  ones63 = [1] * 63

  # Homogeneous & heterogeneous IRS configs strictly beyond unique radius t_u
  irs_configs = [
      # Homogeneous RS(15, 7; l=3) over GF(16): r=8, t_u=4, t_max=6
      (f16, eval15, ones15, [7, 7, 7], 5, "hom_f16_15_7_l3_t5"),
      (f16, eval15, ones15, [7, 7, 7], 6, "hom_f16_15_7_l3_t6"),
      # Homogeneous GRS(15, 6; l=4) with non-trivial multipliers: r=9, t_u=4, t_max=7
      (f16, eval15, mults15, [6, 6, 6, 6], 7, "hom_grs_f16_15_6_l4_t7"),
      # Prime field F_17: RS(16, 6; l=3): r=10, t_u=5, t_max=7
      (f17, eval16_f17, ones16, [6, 6, 6], 7, "hom_f17_16_6_l3_t7"),
      # Heterogeneous IRS(15; [5, 7, 9]) over GF(16): r=[10, 8, 6], t_u=3, t_max=6
      (f16, eval15, ones15, [5, 7, 9], 6, "het_f16_15_5_7_9_t6"),
      # Heterogeneous GRS(31; [9, 13, 17]) over GF(32): r=[22, 18, 14], t_u=7, t_max=13
      (f32, eval31, mults31, [9, 13, 17], 12, "het_grs_f32_31_9_13_17_t12"),
      (f32, eval31, mults31, [9, 13, 17], 13, "het_grs_f32_31_9_13_17_t13"),
      # Homogeneous RS(63, 27; l=3) over GF(64): r=36, t_u=18, t_max=27
      (f64, eval63, ones63, [27, 27, 27], 27, "hom_f64_63_27_l3_t27"),
  ]

  num_irs_runs = max(len(irs_configs), trials)
  for idx in range(num_irs_runs):
    ff, evals, mults, dims, t_err, name = irs_configs[idx % len(irs_configs)]
    _verify_irs_config(
        client, ff, evals, mults, dims, t_err, rng, result, f"{name}#i{idx}"
    )

  # Power Decoding on low-rate RS(63, 10) over GF(64) at t = 30, 31, 32
  power_trials_per_t = max(2, trials // 8)
  _verify_power_decoding_rs63_10(client, f64, rng, result, power_trials_per_t)

  return result
