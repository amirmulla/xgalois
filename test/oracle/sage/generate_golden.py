#!/usr/bin/env python3
"""Deterministic Tier-1 Oracle Golden Vector Generator for xgalois v0.3.

Pure Python 3 (non-preparsed): executes under BOTH `sage -python` (using
SageMath 10.5+ directly when `sage.all` is importable) AND plain `python3`
(using exact arbitrary-precision reference math).

Usage:
  python3 test/oracle/sage/generate_golden.py \\
      --seed=2026 --out=test/oracle/golden/v0.3_algebra_fields_poly.json

  python3 test/oracle/sage/generate_golden.py \\
      --seed=2026 --check-dir=test/oracle/golden
"""

from __future__ import annotations

import argparse
import json
import math
import pathlib
import random
import sys
from typing import Any, Dict, List, Optional, Sequence, Tuple

try:
  import sage.all as sage  # type: ignore[import-not-found]

  HAS_SAGE = True
except ImportError:
  sage = None
  HAS_SAGE = False


GOLDEN_FILENAME = "v0.3_algebra_fields_poly.json"
SCHEMA_VERSION = "1.0"


# ============================================================================
# Canonical binary / prime / polynomial helpers (pure Python exact math)
# ============================================================================


def poly_ascending_to_mask(coeffs: Sequence[int]) -> int:
  mask = 0
  for deg, c in enumerate(coeffs):
    if c & 1:
      mask |= 1 << deg
  return mask


def exponents_to_ascending(m: int, exponents: Sequence[int]) -> List[int]:
  coeffs = [0] * (m + 1)
  for e in exponents:
    coeffs[e] = 1
  return coeffs


def gf2m_mul(a: int, b: int, mod_mask: int, m: int) -> int:
  res = 0
  top_bit = 1 << m
  while b:
    if b & 1:
      res ^= a
    a <<= 1
    if a & top_bit:
      a ^= mod_mask
    b >>= 1
  return res


def gf2m_pow(base: int, exp: int, mod_mask: int, m: int) -> int:
  res = 1
  b = base
  e = exp
  while e > 0:
    if e & 1:
      res = gf2m_mul(res, b, mod_mask, m)
    b = gf2m_mul(b, b, mod_mask, m)
    e >>= 1
  return res


def gf2m_inv(a: int, mod_mask: int, m: int) -> int:
  if a == 0:
    raise ZeroDivisionError("Zero has no multiplicative inverse in GF(2^m)")
  return gf2m_pow(a, (1 << m) - 2, mod_mask, m)


def gf2m_trace(a: int, mod_mask: int, m: int) -> int:
  acc = a
  cur = a
  for _ in range(1, m):
    cur = gf2m_mul(cur, cur, mod_mask, m)
    acc ^= cur
  assert acc in (0, 1), f"Trace must lie in GF(2), got {acc}"
  return acc


def gf2m_half_trace(a: int, mod_mask: int, m: int) -> int:
  assert m % 2 == 1, "HalfTrace is only defined for odd degree m"
  acc = a
  cur = a
  for _ in range((m - 1) // 2):
    cur = gf2m_mul(cur, cur, mod_mask, m)
    cur = gf2m_mul(cur, cur, mod_mask, m)
    acc ^= cur
  return acc


def tonelli_shanks_canonical(a: int, p: int) -> Optional[int]:
  """Returns canonical sqrt min(r, p - r) mod odd prime p, or None."""
  a %= p
  if a == 0:
    return 0
  if pow(a, (p - 1) // 2, p) != 1:
    return None
  if p % 4 == 3:
    r = pow(a, (p + 1) // 4, p)
  else:
    q = p - 1
    s = 0
    while q % 2 == 0:
      q //= 2
      s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1:
      z += 1
    m_val = s
    c = pow(z, q, p)
    t = pow(a, q, p)
    r = pow(a, (q + 1) // 2, p)
    while t != 1:
      i = 0
      tmp = t
      while tmp != 1:
        tmp = (tmp * tmp) % p
        i += 1
      b = pow(c, 1 << (m_val - i - 1), p)
      m_val = i
      c = (b * b) % p
      t = (t * c) % p
      r = (r * b) % p
  assert (r * r) % p == a
  return min(r, p - r)


# ============================================================================
# Quadratic integer rings: Z[i] and Z[omega]
# ============================================================================


def round_div(num: int, den: int) -> int:
  q, r = divmod(num, den)
  if 2 * r > den:
    q += 1
  return q


def gaussian_mul(x: Tuple[int, int], y: Tuple[int, int]) -> Tuple[int, int]:
  return (x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0])


def gaussian_norm(x: Tuple[int, int]) -> int:
  return x[0] * x[0] + x[1] * x[1]


def gaussian_canonical(x: Tuple[int, int]) -> Tuple[int, int]:
  a, b = x
  if a == 0 and b == 0:
    return (0, 0)
  for u in ((1, 0), (0, -1), (-1, 0), (0, 1)):
    c = gaussian_mul(u, x)
    if c[0] > 0 and c[1] >= 0:
      return c
  return x


def gaussian_gcd(a: Tuple[int, int], b: Tuple[int, int]) -> Tuple[int, int]:
  x, y = a, b
  while y != (0, 0):
    ny = gaussian_norm(y)
    q = (
        round_div(x[0] * y[0] + x[1] * y[1], ny),
        round_div(x[1] * y[0] - x[0] * y[1], ny),
    )
    qy = gaussian_mul(q, y)
    r = (x[0] - qy[0], x[1] - qy[1])
    x, y = y, r
  return gaussian_canonical(x)


def eisenstein_mul(x: Tuple[int, int], y: Tuple[int, int]) -> Tuple[int, int]:
  bd = x[1] * y[1]
  return (x[0] * y[0] - bd, x[0] * y[1] + x[1] * y[0] - bd)


def eisenstein_conj(x: Tuple[int, int]) -> Tuple[int, int]:
  return (x[0] - x[1], -x[1])


def eisenstein_norm(x: Tuple[int, int]) -> int:
  return x[0] * x[0] - x[0] * x[1] + x[1] * x[1]


def eisenstein_canonical(x: Tuple[int, int]) -> Tuple[int, int]:
  if x == (0, 0):
    return (0, 0)
  units = ((1, 0), (1, 1), (0, 1), (-1, 0), (-1, -1), (0, -1))
  for u in units:
    c = eisenstein_mul(u, x)
    if c[1] >= 0 and c[0] > c[1]:
      return c
  return x


def eisenstein_divrem(
    x: Tuple[int, int], y: Tuple[int, int]
) -> Tuple[Tuple[int, int], Tuple[int, int]]:
  ny = eisenstein_norm(y)
  prod = eisenstein_mul(x, eisenstein_conj(y))
  q0 = (round_div(prod[0], ny), round_div(prod[1], ny))
  best_q = q0
  qy = eisenstein_mul(best_q, y)
  best_r = (x[0] - qy[0], x[1] - qy[1])
  best_n = eisenstein_norm(best_r)
  for da, db in ((1, 0), (-1, 0), (0, 1), (0, -1), (1, 1), (-1, -1)):
    cq = (q0[0] + da, q0[1] + db)
    cqy = eisenstein_mul(cq, y)
    cr = (x[0] - cqy[0], x[1] - cqy[1])
    cn = eisenstein_norm(cr)
    if cn < best_n:
      best_q, best_r, best_n = cq, cr, cn
  return best_q, best_r


def eisenstein_gcd(a: Tuple[int, int], b: Tuple[int, int]) -> Tuple[int, int]:
  x, y = a, b
  while y != (0, 0):
    _, r = eisenstein_divrem(x, y)
    x, y = y, r
  return eisenstein_canonical(x)


# ============================================================================
# Generalized CRT, Matrix Berkowitz, Skew Polynomials, Field Towers
# ============================================================================


def extended_gcd_int(a: int, b: int) -> Tuple[int, int, int]:
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


def generalized_crt_int(
    residues: Sequence[int], moduli: Sequence[int]
) -> Optional[Tuple[int, int]]:
  cur_mod = abs(moduli[0])
  cur_rem = residues[0] % cur_mod
  for r_i, m_i_raw in zip(residues[1:], moduli[1:]):
    m_i = abs(m_i_raw)
    r_i = r_i % m_i
    g, s, _ = extended_gcd_int(cur_mod, m_i)
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


def matrix_det_and_charpoly_mod(
    mat: Sequence[int], n: int, mod: int
) -> Tuple[int, List[int]]:
  """Exact Berkowitz charpoly det(lambda*I - A) (descending) & det mod `mod`."""
  if HAS_SAGE:
    ring = sage.Zmod(mod)
    m_sage = sage.matrix(ring, n, n, list(mat))
    cp = m_sage.charpoly()
    # Sage charpoly() is monic degree n: x^n + c_1 x^{n-1} + ... + c_n
    coeffs_desc = [int(cp[n - i]) % mod for i in range(n + 1)]
    det_val = int(m_sage.determinant()) % mod
    return det_val, coeffs_desc

  poly = [1 % mod]
  for k in range(n - 1, -1, -1):
    m = n - 1 - k
    toeplitz = [0] * (m + 2)
    toeplitz[0] = 1 % mod
    toeplitz[1] = (-mat[k * n + k]) % mod
    if m > 0:
      v = [mat[(k + 1 + i) * n + k] % mod for i in range(m)]
      for j in range(2, m + 2):
        dot = 0
        for i in range(m):
          dot = (dot + mat[k * n + (k + 1 + i)] * v[i]) % mod
        toeplitz[j] = (-dot) % mod
        if j < m + 1:
          next_v = [0] * m
          for r in range(m):
            acc = 0
            for c in range(m):
              acc = (acc + mat[(k + 1 + r) * n + (k + 1 + c)] * v[c]) % mod
            next_v[r] = acc
          v = next_v
    next_poly = [0] * (m + 2)
    for r in range(m + 2):
      acc = 0
      for c in range(min(r, m) + 1):
        acc = (acc + toeplitz[r - c] * poly[c]) % mod
      next_poly[r] = acc
    poly = next_poly
  det_val = (-poly[-1]) % mod if (n % 2 == 1) else (poly[-1] % mod)
  return det_val, poly


def skew_mul_gf2m(
    f: Sequence[int], g: Sequence[int], mod_mask: int, m: int
) -> List[int]:
  """Skew polynomial multiplication over GF(2^m)[x; sigma(a)=a^2]."""
  if not f or not g:
    return []
  out = [0] * (len(f) + len(g) - 1)
  h = list(g)
  for fk in f:
    if fk != 0:
      for j, hj in enumerate(h):
        if hj != 0:
          out[j] ^= gf2m_mul(fk, hj, mod_mask, m)
    # h <- x * h using x * c = c^2 * x
    h = [0] + [gf2m_mul(c, c, mod_mask, m) for c in h]
  while out and out[-1] == 0:
    out.pop()
  return out


def skew_eval_operator_gf2m(
    f: Sequence[int], x: int, mod_mask: int, m: int
) -> int:
  acc = 0
  cur = x
  for fi in f:
    if fi != 0 and cur != 0:
      acc ^= gf2m_mul(fi, cur, mod_mask, m)
    cur = gf2m_mul(cur, cur, mod_mask, m)
  return acc


def skew_minimal_poly_gf2m(
    points: Sequence[int], mod_mask: int, m: int
) -> List[int]:
  poly = [1]
  for pt in points:
    beta = skew_eval_operator_gf2m(poly, pt, mod_mask, m)
    if beta == 0:
      continue
    # c = sigma(beta) / beta = beta^2 / beta = beta in char 2 with sigma(a)=a^2
    c = gf2m_mul(gf2m_mul(beta, beta, mod_mask, m), gf2m_inv(beta, mod_mask, m), mod_mask, m)
    poly = skew_mul_gf2m([c, 1], poly, mod_mask, m)
  return poly


def tower_f7_2_mul(
    a: Tuple[int, int], b: Tuple[int, int], xi: int, p: int
) -> Tuple[int, int]:
  """Multiplication in F_p[u]/(u^2 - xi)."""
  return (
      (a[0] * b[0] + xi * a[1] * b[1]) % p,
      (a[0] * b[1] + a[1] * b[0]) % p,
  )


def tower_f7_2_inv(a: Tuple[int, int], xi: int, p: int) -> Tuple[int, int]:
  norm = (a[0] * a[0] - xi * a[1] * a[1]) % p
  inv_norm = pow(norm, p - 2, p)
  return ((a[0] * inv_norm) % p, ((-a[1]) * inv_norm) % p)


# ============================================================================
# Vector Generation
# ============================================================================


def build_group_vectors() -> List[Dict[str, Any]]:
  vectors: List[Dict[str, Any]] = []

  if HAS_SAGE:
    s5_order = str(int(sage.SymmetricGroup(5).order()))
    a5 = sage.AlternatingGroup(5)
    a5_order = str(int(a5.order()))
    a5_solvable = bool(a5.is_solvable())
    d8 = sage.DihedralGroup(8)
    d8_order = str(int(d8.order()))
    d8_center = int(d8.center().order())
    d4 = sage.DihedralGroup(4)
    d4_order = str(int(d4.order()))
    d4_center = int(d4.center().order())
    s20_order = str(int(sage.SymmetricGroup(20).order()))
    s4 = sage.SymmetricGroup(4)
    s4_order = str(int(s4.order()))
    s4_comm_order = int(s4.commutator().order())
    s4_syl2_order = int(s4.sylow_subgroup(2).order())
    s4_syl3_order = int(s4.sylow_subgroup(3).order())
    s4_solvable = bool(s4.is_solvable())
    snf_mat = sage.matrix(sage.ZZ, [[2, 0, 0], [0, 4, 0], [0, 0, 6]])
    snf_diag = snf_mat.smith_form()[0].diagonal()
    snf_factors = [str(int(d)) for d in snf_diag if int(d) > 1]
  else:
    s5_order = str(math.factorial(5))
    a5_order = str(math.factorial(5) // 2)
    a5_solvable = False
    d8_order = str(2 * 8)
    d8_center = 2  # Z(D_n) = {e, r^(n/2)} for even n >= 2
    d4_order = str(2 * 4)
    d4_center = 2
    s20_order = str(math.factorial(20))
    s4_order = str(math.factorial(4))
    s4_comm_order = 12  # [S_4, S_4] = A_4
    s4_syl2_order = 8  # 24 = 2^3 * 3
    s4_syl3_order = 3
    s4_solvable = True
    # diag(2, 4, 6): d1 = gcd(2,4,6)=2, d1*d2 = gcd(8,12,24)=4 => d2=2, d3=12
    snf_factors = ["2", "2", "12"]

  vectors.append({
      "id": "grp_s5_order",
      "category": "groups",
      "operation": "symmetric_group_order",
      "status": "required",
      "degree": 5,
      "expected_order": s5_order,
  })
  vectors.append({
      "id": "grp_a5_order_nonsolvable",
      "category": "groups",
      "operation": "alternating_group_properties",
      "status": "required",
      "degree": 5,
      "expected_order": a5_order,
      "expected_is_solvable": a5_solvable,
  })
  vectors.append({
      "id": "grp_d8_order_center",
      "category": "groups",
      "operation": "dihedral_group_properties",
      "status": "required",
      "n": 8,
      "expected_order": d8_order,
      "expected_center_order": d8_center,
  })
  vectors.append({
      "id": "grp_d4_order_center",
      "category": "groups",
      "operation": "dihedral_group_properties",
      "status": "required",
      "n": 4,
      "expected_order": d4_order,
      "expected_center_order": d4_center,
  })
  vectors.append({
      "id": "grp_s20_order",
      "category": "groups",
      "operation": "symmetric_group_order",
      "status": "required",
      "degree": 20,
      "expected_order": s20_order,
  })
  vectors.append({
      "id": "grp_s4_structure",
      "category": "groups",
      "operation": "symmetric_group_structure",
      "status": "required",
      "degree": 4,
      "expected_order": s4_order,
      "expected_commutator_order": s4_comm_order,
      "expected_sylow_2_order": s4_syl2_order,
      "expected_sylow_3_order": s4_syl3_order,
      "expected_is_solvable": s4_solvable,
  })
  vectors.append({
      "id": "grp_snf_z3_2_4_6",
      "category": "groups",
      "operation": "snf_abelian_decomposition",
      "status": "required",
      "relations": [[2, 0, 0], [0, 4, 0], [0, 0, 6]],
      "expected_free_rank": 0,
      "expected_invariant_factors": snf_factors,
  })
  return vectors


def build_ring_vectors(rng: random.Random) -> List[Dict[str, Any]]:
  vectors: List[Dict[str, Any]] = []

  # 1. Z/12Z unit group & totient
  if HAS_SAGE:
    r12 = sage.Zmod(12)
    z12_units = sorted(int(u) for u in r12 if u.is_unit())
    z12_totient = int(sage.euler_phi(12))
  else:
    z12_units = [u for u in range(1, 12) if math.gcd(u, 12) == 1]
    z12_totient = len(z12_units)
  vectors.append({
      "id": "ring_z12_units",
      "category": "rings",
      "operation": "modular_ring_units",
      "status": "required",
      "modulus": 12,
      "expected_units": z12_units,
      "expected_unit_count": len(z12_units),
      "expected_euler_totient": z12_totient,
  })

  # 2. Z[i] Euclidean GCD & norm
  # Choose deterministic non-trivial pair sharing Gaussian prime factor (2 + i) * (1 + i)
  g_a = (32, 9)  # 32 + 9i
  g_b = (4, 11)  # 4 + 11i
  # Also construct a pair with nontrivial shared factor (2 + i) * (3 + 2i) = 4 + 7i
  g_c = gaussian_mul((2, 1), (7, 4))  # (10, 15)
  g_d = gaussian_mul((2, 1), (3, -5))  # (11, -7)
  if HAS_SAGE:
    zi = sage.QuadraticField(-1, "i").ring_of_integers()
    i_gen = zi.gen(1)
    g1_sage = sage.gcd(zi(g_a[0] + g_a[1] * i_gen), zi(g_b[0] + g_b[1] * i_gen))
    g2_sage = sage.gcd(zi(g_c[0] + g_c[1] * i_gen), zi(g_d[0] + g_d[1] * i_gen))
    g1_norm = int(g1_sage.norm())
    g2_norm = int(g2_sage.norm())
    g1_canon = gaussian_canonical(gaussian_gcd(g_a, g_b))
    g2_canon = gaussian_canonical(gaussian_gcd(g_c, g_d))
    assert gaussian_norm(g1_canon) == g1_norm
    assert gaussian_norm(g2_canon) == g2_norm
  else:
    g1_canon = gaussian_gcd(g_a, g_b)
    g1_norm = gaussian_norm(g1_canon)
    g2_canon = gaussian_gcd(g_c, g_d)
    g2_norm = gaussian_norm(g2_canon)

  vectors.append({
      "id": "ring_gaussian_gcd",
      "category": "rings",
      "operation": "gaussian_integer_gcd",
      "status": "required",
      "cases": [
          {
              "a": [g_a[0], g_a[1]],
              "b": [g_b[0], g_b[1]],
              "expected_gcd_canonical": [g1_canon[0], g1_canon[1]],
              "expected_gcd_norm": g1_norm,
          },
          {
              "a": [g_c[0], g_c[1]],
              "b": [g_d[0], g_d[1]],
              "expected_gcd_canonical": [g2_canon[0], g2_canon[1]],
              "expected_gcd_norm": g2_norm,
          },
      ],
  })

  # 3. Z[omega] Euclidean GCD & norm (omega^2 + omega + 1 = 0)
  e_shared = (3, 1)  # norm = 9 - 3 + 1 = 7
  e_a = eisenstein_mul(e_shared, (4, 1))  # norm = 7 * 13
  e_b = eisenstein_mul(e_shared, (5, 2))  # norm = 7 * 19
  e_c = (7, 2)
  e_d = (3, 5)
  if HAS_SAGE:
    k3 = sage.QuadraticField(-3, "s")
    ok = k3.ring_of_integers()
    w = (-1 + k3.gen()) / 2
    assert w * w + w + 1 == 0
    eg1 = sage.gcd(ok(e_a[0] + e_a[1] * w), ok(e_b[0] + e_b[1] * w))
    eg2 = sage.gcd(ok(e_c[0] + e_c[1] * w), ok(e_d[0] + e_d[1] * w))
    e1_norm = int(eg1.norm())
    e2_norm = int(eg2.norm())
    e1_canon = eisenstein_gcd(e_a, e_b)
    e2_canon = eisenstein_gcd(e_c, e_d)
    assert eisenstein_norm(e1_canon) == e1_norm
    assert eisenstein_norm(e2_canon) == e2_norm
  else:
    e1_canon = eisenstein_gcd(e_a, e_b)
    e1_norm = eisenstein_norm(e1_canon)
    e2_canon = eisenstein_gcd(e_c, e_d)
    e2_norm = eisenstein_norm(e2_canon)

  vectors.append({
      "id": "ring_eisenstein_gcd",
      "category": "rings",
      "operation": "eisenstein_integer_gcd",
      "status": "required",
      "cases": [
          {
              "a": [e_a[0], e_a[1]],
              "b": [e_b[0], e_b[1]],
              "expected_gcd_canonical": [e1_canon[0], e1_canon[1]],
              "expected_gcd_norm": e1_norm,
          },
          {
              "a": [e_c[0], e_c[1]],
              "b": [e_d[0], e_d[1]],
              "expected_gcd_canonical": [e2_canon[0], e2_canon[1]],
              "expected_gcd_norm": e2_norm,
          },
      ],
  })

  # 4. Z_17[x]/(x^4 + 1) negacyclic product
  q_neg = 17
  n_neg = 4
  neg_a = [rng.randrange(q_neg) for _ in range(n_neg)]
  neg_b = [rng.randrange(q_neg) for _ in range(n_neg)]
  if HAS_SAGE:
    r_poly = sage.PolynomialRing(sage.Zmod(q_neg), "x")
    x_var = r_poly.gen()
    q_ring = r_poly.quotient(x_var**n_neg + 1)
    prod_q = q_ring(neg_a) * q_ring(neg_b)
    neg_prod = [int(prod_q.lift()[i]) % q_neg for i in range(n_neg)]
  else:
    neg_prod = [0] * n_neg
    for i in range(n_neg):
      for j in range(n_neg):
        s = i + j
        if s < n_neg:
          neg_prod[s] = (neg_prod[s] + neg_a[i] * neg_b[j]) % q_neg
        else:
          neg_prod[s - n_neg] = (neg_prod[s - n_neg] - neg_a[i] * neg_b[j]) % q_neg

  vectors.append({
      "id": "ring_negacyclic_z17_x4",
      "category": "rings",
      "operation": "negacyclic_ring_mul",
      "status": "required",
      "modulus_q": q_neg,
      "degree": n_neg,
      "a": neg_a,
      "b": neg_b,
      "expected_product": neg_prod,
  })

  # 5. M_3(Z/6Z) Berkowitz determinant & characteristic polynomial
  m3_mod = 6
  m3_mat = [2, 3, 1, 4, 5, 2, 1, 2, 3]
  m3_det, m3_charpoly = matrix_det_and_charpoly_mod(m3_mat, 3, m3_mod)
  vectors.append({
      "id": "ring_matrix_m3_z6_det",
      "category": "rings",
      "operation": "matrix_ring_berkowitz_det",
      "status": "required",
      "modulus": m3_mod,
      "n": 3,
      "matrix": m3_mat,
      "expected_det": m3_det,
      "expected_char_poly_desc": m3_charpoly,
  })

  # 6. Generalized CRT over Z (coprime, non-coprime consistent, non-coprime inconsistent)
  crt_coprime_res = [2, 3, 2]
  crt_coprime_mod = [3, 5, 7]
  if HAS_SAGE:
    sol_cp = int(sage.crt(crt_coprime_res, crt_coprime_mod))
    lcm_cp = int(sage.lcm(crt_coprime_mod))
  else:
    res_cp = generalized_crt_int(crt_coprime_res, crt_coprime_mod)
    assert res_cp is not None
    sol_cp, lcm_cp = res_cp

  crt_nc_res = [3, 7, 13]
  crt_nc_mod = [6, 10, 14]
  # Check consistency: 3=1 mod 2, 7=1 mod 2, 13=1 mod 2 => consistent!
  res_nc = generalized_crt_int(crt_nc_res, crt_nc_mod)
  assert res_nc is not None
  sol_nc, lcm_nc = res_nc
  if HAS_SAGE:
    assert int(sage.CRT_list(crt_nc_res, crt_nc_mod)) == sol_nc
    assert int(sage.lcm(crt_nc_mod)) == lcm_nc

  crt_inc_res = [1, 2]
  crt_inc_mod = [4, 6]
  assert generalized_crt_int(crt_inc_res, crt_inc_mod) is None

  vectors.append({
      "id": "ring_crt_integer",
      "category": "rings",
      "operation": "generalized_crt_int",
      "status": "required",
      "cases": [
          {
              "name": "coprime",
              "residues": crt_coprime_res,
              "moduli": crt_coprime_mod,
              "expected_consistent": True,
              "expected_solution": sol_cp,
              "expected_lcm": lcm_cp,
          },
          {
              "name": "non_coprime_consistent",
              "residues": crt_nc_res,
              "moduli": crt_nc_mod,
              "expected_consistent": True,
              "expected_solution": sol_nc,
              "expected_lcm": lcm_nc,
          },
          {
              "name": "non_coprime_inconsistent",
              "residues": crt_inc_res,
              "moduli": crt_inc_mod,
              "expected_consistent": False,
          },
      ],
  })

  # 7. FrobeniusSkewRing over GF(2^4, x^4 + x + 1)
  gf16_mod = [1, 1, 0, 0, 1]
  gf16_mask = poly_ascending_to_mask(gf16_mod)
  skew_f = [3, 5, 11]  # 3 + 5*x + 11*x^2 in GF(2^4)[x; sigma]
  skew_x = 7  # 1 + alpha + alpha^2
  skew_eval = skew_eval_operator_gf2m(skew_f, skew_x, gf16_mask, 4)
  indep_pts = [1, 2, 4]  # 1, alpha, alpha^2 (F_2-dim = 3)
  dep_pts = [1, 2, 3]  # 1, alpha, 1 + alpha (F_2-dim = 2)
  min_poly_indep = skew_minimal_poly_gf2m(indep_pts, gf16_mask, 4)
  min_poly_dep = skew_minimal_poly_gf2m(dep_pts, gf16_mask, 4)
  for pt in indep_pts:
    assert skew_eval_operator_gf2m(min_poly_indep, pt, gf16_mask, 4) == 0
  for pt in dep_pts:
    assert skew_eval_operator_gf2m(min_poly_dep, pt, gf16_mask, 4) == 0

  vectors.append({
      "id": "ring_skew_frobenius_gf16",
      "category": "rings",
      "operation": "frobenius_skew_ring",
      "status": "required",
      "field": {"p": 2, "m": 4, "modulus": gf16_mod},
      "frobenius_step": 1,
      "eval_poly": skew_f,
      "eval_point": skew_x,
      "expected_eval": skew_eval,
      "minimal_poly_cases": [
          {
              "points": indep_pts,
              "expected_degree": len(min_poly_indep) - 1,
              "expected_coeffs": min_poly_indep,
          },
          {
              "points": dep_pts,
              "expected_degree": len(min_poly_dep) - 1,
              "expected_coeffs": min_poly_dep,
          },
      ],
  })

  return vectors


def build_field_vectors(rng: random.Random) -> List[Dict[str, Any]]:
  vectors: List[Dict[str, Any]] = []

  # 1. GF(2^8) with 0x11D and AES 0x11B
  def make_gf256_vector(
      vec_id: str,
      mod_coeffs: List[int],
      fixed_cases: List[Tuple[int, int, int]],
  ) -> Dict[str, Any]:
    mask = poly_ascending_to_mask(mod_coeffs)
    cases = []
    if HAS_SAGE:
      r2 = sage.PolynomialRing(sage.GF(2), "x")
      x = r2.gen()
      mod_poly = sum(c * x**i for i, c in enumerate(mod_coeffs))
      f256 = sage.GF(2**8, name="a", modulus=mod_poly)
      for a_int, b_int, exp in fixed_cases:
        a_el = f256.from_integer(a_int)
        b_el = f256.from_integer(b_int)
        mul_val = int((a_el * b_el).to_integer())
        inv_val = int((a_el ** (-1)).to_integer())
        pow_val = int((a_el**exp).to_integer())
        cases.append({
            "a": a_int,
            "b": b_int,
            "exp": exp,
            "expected_mul": mul_val,
            "expected_inv": inv_val,
            "expected_pow": pow_val,
        })
    else:
      for a_int, b_int, exp in fixed_cases:
        cases.append({
            "a": a_int,
            "b": b_int,
            "exp": exp,
            "expected_mul": gf2m_mul(a_int, b_int, mask, 8),
            "expected_inv": gf2m_inv(a_int, mask, 8),
            "expected_pow": gf2m_pow(a_int, exp, mask, 8),
        })
    return {
        "id": vec_id,
        "category": "fields",
        "operation": "gf2_extension_arith",
        "status": "required",
        "field": {"p": 2, "m": 8, "modulus": mod_coeffs},
        "cases": cases,
    }

  cases_11d = [
      (0x57, 0x83, 17),
      (rng.randint(1, 255), rng.randint(1, 255), rng.randint(2, 250)),
      (rng.randint(1, 255), rng.randint(1, 255), rng.randint(2, 250)),
  ]
  vectors.append(
      make_gf256_vector(
          "field_gf256_0x11d", [1, 0, 1, 1, 1, 0, 0, 0, 1], cases_11d
      )
  )

  cases_11b = [
      (0x57, 0x83, 13),  # FIPS-197: 0x57 * 0x83 = 0xc1
      (0x53, 0xCA, 254),  # FIPS-197: 0x53 * 0xca = 1, 0x53^-1 = 0xca
      (rng.randint(1, 255), rng.randint(1, 255), rng.randint(2, 250)),
  ]
  vectors.append(
      make_gf256_vector(
          "field_gf256_aes_0x11b", [1, 1, 0, 1, 1, 0, 0, 0, 1], cases_11b
      )
  )

  # 2. FieldTower F_{7^2} = F_7[x]/(x^2 - 3) -> modulus [4, 0, 1]
  p_tower = 7
  xi_tower = 3
  mod_tower = [(-xi_tower) % p_tower, 0, 1]  # [4, 0, 1]
  ta = (2, 5)
  tb = (4, 3)
  if HAS_SAGE:
    r7 = sage.PolynomialRing(sage.GF(p_tower), "u")
    u = r7.gen()
    f49 = sage.GF(p_tower**2, name="w", modulus=u**2 - xi_tower)
    w = f49.gen()
    ea = ta[0] + ta[1] * w
    eb = tb[0] + tb[1] * w
    emul = ea * eb
    einv = ea ** (-1)
    t_mul = [int(emul.polynomial()[0]), int(emul.polynomial()[1])]
    t_inv = [int(einv.polynomial()[0]), int(einv.polynomial()[1])]
    t_trace = int(ea.trace())
    t_norm = int(ea.norm())
  else:
    t_mul = list(tower_f7_2_mul(ta, tb, xi_tower, p_tower))
    t_inv = list(tower_f7_2_inv(ta, xi_tower, p_tower))
    t_trace = (2 * ta[0]) % p_tower
    t_norm = (ta[0] * ta[0] - xi_tower * ta[1] * ta[1]) % p_tower

  vectors.append({
      "id": "field_tower_f7_2",
      "category": "fields",
      "operation": "field_tower_relative_ops",
      "status": "required",
      "field": {"p": p_tower, "m": 2, "modulus": mod_tower},
      "a": list(ta),
      "b": list(tb),
      "expected_mul": t_mul,
      "expected_inv": t_inv,
      "expected_relative_trace": t_trace,
      "expected_relative_norm": t_norm,
  })

  # 3. CyclotomicCosets(15, 2)
  n_cyc, q_cyc = 15, 2
  visited = [False] * n_cyc
  cosets: List[List[int]] = []
  for s in range(n_cyc):
    if visited[s]:
      continue
    coset = []
    cur = s
    while not visited[cur]:
      visited[cur] = True
      coset.append(cur)
      cur = (cur * q_cyc) % n_cyc
    cosets.append(coset)
  vectors.append({
      "id": "field_cyclotomic_cosets_15_2",
      "category": "fields",
      "operation": "cyclotomic_cosets",
      "status": "required",
      "n": n_cyc,
      "q": q_cyc,
      "expected_cosets": cosets,
  })

  # 4. GaloisFieldPrimeBig<256>: Curve25519 and secp256k1
  def make_prime_big_vector(
      vec_id: str, curve_name: str, prime: int
  ) -> Dict[str, Any]:
    a = rng.randrange(2, prime - 1)
    b = rng.randrange(2, prime - 1)
    root_seed = rng.randrange(2, prime - 1)
    qr = (root_seed * root_seed) % prime
    # Find deterministic quadratic non-residue
    qnr = rng.randrange(2, prime - 1)
    while pow(qnr, (prime - 1) // 2, prime) != prime - 1:
      qnr = (qnr + 1) % prime
      if qnr < 2:
        qnr = 2
    if HAS_SAGE:
      fp = sage.GF(prime)
      mul_val = int(fp(a) * fp(b))
      inv_val = int(fp(a) ** (-1))
      r_sage = int(fp(qr).sqrt())
      sqrt_val = min(r_sage, prime - r_sage)
      assert not fp(qnr).is_square()
    else:
      mul_val = (a * b) % prime
      inv_val = pow(a, prime - 2, prime)
      sqrt_val = tonelli_shanks_canonical(qr, prime)
      assert sqrt_val is not None
      assert tonelli_shanks_canonical(qnr, prime) is None
    return {
        "id": vec_id,
        "category": "fields",
        "operation": "prime_big_256_arith",
        "status": "required",
        "curve": curve_name,
        "p_hex": hex(prime),
        "a_hex": hex(a),
        "b_hex": hex(b),
        "expected_mul_hex": hex(mul_val),
        "expected_inv_hex": hex(inv_val),
        "qr_hex": hex(qr),
        "expected_sqrt_hex": hex(sqrt_val),
        "qnr_hex": hex(qnr),
    }

  p_25519 = (1 << 255) - 19
  p_secp256k1 = (1 << 256) - (1 << 32) - 977
  vectors.append(
      make_prime_big_vector("field_prime_big_curve25519", "curve25519", p_25519)
  )
  vectors.append(
      make_prime_big_vector("field_prime_big_secp256k1", "secp256k1", p_secp256k1)
  )

  # 5. GF(2^163) NIST B-163: x^163 + x^7 + x^6 + x^3 + 1
  m_163 = 163
  exp_163 = [163, 7, 6, 3, 0]
  mod_163_asc = exponents_to_ascending(m_163, exp_163)
  mask_163 = poly_ascending_to_mask(mod_163_asc)
  a_163 = rng.randrange(1, 1 << m_163)
  b_163 = rng.randrange(1, 1 << m_163)

  # Ensure one trace-0 element and one trace-1 element
  tr0_elem = rng.randrange(1, 1 << m_163)
  if gf2m_trace(tr0_elem, mask_163, m_163) != 0:
    tr0_elem ^= 1  # Tr(1) = 163 mod 2 = 1 in GF(2^163)
  assert gf2m_trace(tr0_elem, mask_163, m_163) == 0

  tr1_elem = tr0_elem ^ 1
  assert gf2m_trace(tr1_elem, mask_163, m_163) == 1

  if HAS_SAGE:
    r2 = sage.PolynomialRing(sage.GF(2), "x")
    x = r2.gen()
    mod_p163 = sum(x**e for e in exp_163)
    f163 = sage.GF(2**m_163, name="t", modulus=mod_p163)
    ea = f163.from_integer(a_163)
    eb = f163.from_integer(b_163)
    mul_163 = int((ea * eb).to_integer())
    inv_163 = int((ea ** (-1)).to_integer())
    tr_a_163 = int(ea.trace())
    e_tr0 = f163.from_integer(tr0_elem)
    assert int(e_tr0.trace()) == 0
    ht_acc = e_tr0
    ht_cur = e_tr0
    for _ in range((m_163 - 1) // 2):
      ht_cur = ht_cur**4
      ht_acc += ht_cur
    assert ht_acc * ht_acc + ht_acc == e_tr0
    half_trace_163 = int(ht_acc.to_integer())
  else:
    mul_163 = gf2m_mul(a_163, b_163, mask_163, m_163)
    inv_163 = gf2m_inv(a_163, mask_163, m_163)
    tr_a_163 = gf2m_trace(a_163, mask_163, m_163)
    half_trace_163 = gf2m_half_trace(tr0_elem, mask_163, m_163)
    ht_sq = gf2m_mul(half_trace_163, half_trace_163, mask_163, m_163)
    assert (ht_sq ^ half_trace_163) == tr0_elem

  vectors.append({
      "id": "field_gf2_163_trace_halftrace",
      "category": "fields",
      "operation": "gf2_binary_large_arith",
      "status": "required",
      "field": {
          "p": 2,
          "m": m_163,
          "modulus": mod_163_asc,
          "modulus_exponents": exp_163,
      },
      "a_hex": hex(a_163),
      "b_hex": hex(b_163),
      "expected_mul_hex": hex(mul_163),
      "expected_inv_hex": hex(inv_163),
      "expected_trace_a": tr_a_163,
      "trace0_elem_hex": hex(tr0_elem),
      "expected_half_trace_hex": hex(half_trace_163),
      "trace1_elem_hex": hex(tr1_elem),
  })

  return vectors


def build_coding_vectors(rng: random.Random) -> List[Dict[str, Any]]:
  vectors: List[Dict[str, Any]] = []
  gf16_mod = [1, 1, 0, 0, 1]  # x^4 + x + 1
  gf16_mask = poly_ascending_to_mask(gf16_mod)

  # Evaluation points alpha^0, alpha^1, ..., alpha^14 with primitive alpha = x = 2
  eval_points = []
  cur = 1
  for _ in range(15):
    eval_points.append(cur)
    cur = gf2m_mul(cur, 2, gf16_mask, 4)

  def eval_poly_gf16(coeffs: Sequence[int], x_val: int) -> int:
    acc = 0
    for c in reversed(coeffs):
      acc = gf2m_mul(acc, x_val, gf16_mask, 4) ^ c
    return acc

  # 1. RS(15, 9) over GF(2^4, x^4 + x + 1) with 3 errors (BM unique decoding)
  rs_msg = [rng.randrange(16) for _ in range(9)]
  if HAS_SAGE:
    r2 = sage.PolynomialRing(sage.GF(2), "x")
    x = r2.gen()
    f16 = sage.GF(16, name="a", modulus=x**4 + x + 1)
    s_evals = [f16.from_integer(e) for e in eval_points]
    grs_code = sage.codes.GeneralizedReedSolomonCode(s_evals, 9)
    enc = grs_code.encoder("EvaluationPolynomial")
    poly_ring = enc.message_space()
    msg_poly = poly_ring([f16.from_integer(c) for c in rs_msg])
    rs_cw = [int(c.to_integer()) for c in enc.encode(msg_poly)]
  else:
    rs_cw = [eval_poly_gf16(rs_msg, pt) for pt in eval_points]

  rs_err_positions = [1, 7, 13]
  rs_err = [0] * 15
  for pos in rs_err_positions:
    rs_err[pos] = rng.randint(1, 15)
  rs_received = [c ^ e for c, e in zip(rs_cw, rs_err)]

  vectors.append({
      "id": "coding_rs_15_9_bm_3err",
      "category": "coding",
      "operation": "rs_encode_bm_decode",
      "status": "required",
      "field": {"p": 2, "m": 4, "modulus": gf16_mod},
      "n": 15,
      "k": 9,
      "eval_points": eval_points,
      "message": rs_msg,
      "expected_codeword": rs_cw,
      "error_positions": rs_err_positions,
      "error_vector": rs_err,
      "received": rs_received,
  })

  # 2. GRS(15, 5) over GF(2^4, x^4 + x + 1) with 6 errors (Guruswami-Sudan list decoding)
  # d = 11 => unique radius = 5, Johnson radius = 15 - 1 - floor(sqrt(60)) = 7.
  grs_col_mults = [((i * 3) % 15) + 1 for i in range(15)]
  grs_msg = [rng.randrange(16) for _ in range(5)]
  if HAS_SAGE:
    r2 = sage.PolynomialRing(sage.GF(2), "x")
    x = r2.gen()
    f16 = sage.GF(16, name="a", modulus=x**4 + x + 1)
    s_evals = [f16.from_integer(e) for e in eval_points]
    s_mults = [f16.from_integer(v) for v in grs_col_mults]
    grs_code = sage.codes.GeneralizedReedSolomonCode(s_evals, 5, s_mults)
    enc = grs_code.encoder("EvaluationPolynomial")
    poly_ring = enc.message_space()
    msg_poly = poly_ring([f16.from_integer(c) for c in grs_msg])
    grs_cw = [int(c.to_integer()) for c in enc.encode(msg_poly)]
  else:
    grs_cw = [
        gf2m_mul(v, eval_poly_gf16(grs_msg, pt), gf16_mask, 4)
        for v, pt in zip(grs_col_mults, eval_points)
    ]

  grs_err_positions = [0, 2, 5, 8, 11, 14]
  grs_err = [0] * 15
  for pos in grs_err_positions:
    grs_err[pos] = rng.randint(1, 15)
  grs_received = [c ^ e for c, e in zip(grs_cw, grs_err)]

  vectors.append({
      "id": "coding_grs_15_5_gs_6err",
      "category": "coding",
      "operation": "grs_guruswami_sudan_list_decode",
      "status": "required",
      "field": {"p": 2, "m": 4, "modulus": gf16_mod},
      "n": 15,
      "k": 5,
      "eval_points": eval_points,
      "column_multipliers": grs_col_mults,
      "message": grs_msg,
      "expected_codeword": grs_cw,
      "error_positions": grs_err_positions,
      "error_vector": grs_err,
      "received": grs_received,
      "target_errors": 6,
      "expected_johnson_radius": 7,
  })

  # 3. Promoted v0.4 LRC vectors (required)

  # 3a. Tamo-Barg optimal LRC (12, 6, r=3) over F_13 with g(x) = x^4 (d = 6)
  p13 = 13
  gamma13 = 2  # smallest primitive root mod 13
  omega13_4 = pow(gamma13, (p13 - 1) // 4, p13)  # 8
  h13_4 = [pow(omega13_4, t, p13) for t in range(4)]
  lrc12_evals: List[int] = []
  for ell in range(3):
    rep = pow(gamma13, ell, p13)
    for h_el in h13_4:
      lrc12_evals.append((rep * h_el) % p13)
  lrc12_msg = [rng.randrange(p13) for _ in range(6)]
  if HAS_SAGE:
    f13 = sage.GF(p13)
    pr13 = sage.PolynomialRing(f13, "X")
    x13 = pr13.gen()
    g13 = x13**4
    f_lrc12 = sum(
        f13(lrc12_msg[j * 3 + i]) * (g13**j) * (x13**i)
        for j in range(2)
        for i in range(3)
    )
    lrc12_cw = [int(f_lrc12(f13(pt))) for pt in lrc12_evals]
  else:
    lrc12_cw = []
    for pt in lrc12_evals:
      g_val = pow(pt, 4, p13)
      acc = 0
      for j in range(2):
        gp = pow(g_val, j, p13)
        for i in range(3):
          acc = (acc + lrc12_msg[j * 3 + i] * gp * pow(pt, i, p13)) % p13
      lrc12_cw.append(acc)

  vectors.append({
      "id": "coding_tamo_barg_lrc_12_6_3",
      "category": "coding",
      "operation": "tamo_barg_lrc",
      "status": "required",
      "field": {"p": 13, "m": 1, "modulus": [0, 1]},
      "n": 12,
      "k": 6,
      "locality_r": 3,
      "local_distance": 2,
      "partition_type": "multiplicative",
      "good_poly_degree": 4,
      "expected_min_distance": 6,
      "eval_points": lrc12_evals,
      "message": lrc12_msg,
      "expected_codeword": lrc12_cw,
      "erased_position": 5,
      "expected_symbols_read": 3,
      "description": "Tamo-Barg optimal LRC (12, 6, r=3) over F_13 with g(x) = x^4",
  })

  # 3b. (r, delta)-LRC (15, 6, r=3, delta=3) over GF(16) with g(x) = x^5 (d = 8)
  omega16_5 = gf2m_pow(2, 15 // 5, gf16_mask, 4)  # alpha^3 = 8
  h16_5 = [gf2m_pow(omega16_5, t, gf16_mask, 4) for t in range(5)]
  lrc15_evals: List[int] = []
  for ell in range(3):
    rep = gf2m_pow(2, ell, gf16_mask, 4)
    for h_el in h16_5:
      lrc15_evals.append(gf2m_mul(rep, h_el, gf16_mask, 4))
  lrc15_msg = [rng.randrange(16) for _ in range(6)]
  if HAS_SAGE:
    r2 = sage.PolynomialRing(sage.GF(2), "x")
    x = r2.gen()
    f16 = sage.GF(16, name="a", modulus=x**4 + x + 1)
    pr16 = sage.PolynomialRing(f16, "X")
    x16 = pr16.gen()
    g16 = x16**5
    f_lrc15 = sum(
        f16.from_integer(lrc15_msg[j * 3 + i]) * (g16**j) * (x16**i)
        for j in range(2)
        for i in range(3)
    )
    lrc15_cw = [
        int(f_lrc15(f16.from_integer(pt)).to_integer()) for pt in lrc15_evals
    ]
  else:
    lrc15_cw = []
    for pt in lrc15_evals:
      g_val = gf2m_pow(pt, 5, gf16_mask, 4)
      acc = 0
      for j in range(2):
        gp = gf2m_pow(g_val, j, gf16_mask, 4)
        for i in range(3):
          xp = gf2m_pow(pt, i, gf16_mask, 4)
          term = gf2m_mul(
              lrc15_msg[j * 3 + i], gf2m_mul(gp, xp, gf16_mask, 4), gf16_mask, 4
          )
          acc ^= term
      lrc15_cw.append(acc)

  vectors.append({
      "id": "coding_lrc_rdelta_15_6_r3_d3",
      "category": "coding",
      "operation": "tamo_barg_lrc",
      "status": "required",
      "field": {"p": 2, "m": 4, "modulus": gf16_mod},
      "n": 15,
      "k": 6,
      "locality_r": 3,
      "local_distance": 3,
      "partition_type": "multiplicative",
      "good_poly_degree": 5,
      "expected_min_distance": 8,
      "eval_points": lrc15_evals,
      "message": lrc15_msg,
      "expected_codeword": lrc15_cw,
      "erased_position": 2,
      "multi_erased_positions": [5, 8],
      "expected_symbols_read": 3,
      "description": "Tamo-Barg (r, delta)-LRC (15, 6, r=3, delta=3) over GF(16) with g(x) = x^5",
  })

  # 3c. Multiplicative CRT Availability LRC (12, 4, r=[2, 3]) over F_13
  avail_evals = [pow(gamma13, j, p13) for j in range(12)]
  avail_exps = [
      e for e in range(12) if (e % 3) < 2 and (e % 4) < 3
  ][:4]  # [0, 1, 4, 6]
  assert avail_exps == [0, 1, 4, 6]
  avail_msg = [rng.randrange(p13) for _ in range(4)]
  if HAS_SAGE:
    f13 = sage.GF(p13)
    avail_cw = [
        int(
            sum(
                f13(avail_msg[i]) * (f13(pt) ** avail_exps[i]) for i in range(4)
            )
        )
        for pt in avail_evals
    ]
  else:
    avail_cw = [
        sum(avail_msg[i] * pow(pt, avail_exps[i], p13) for i in range(4)) % p13
        for pt in avail_evals
    ]

  vectors.append({
      "id": "coding_lrc_availability_12_4_r2_3",
      "category": "coding",
      "operation": "lrc_availability",
      "status": "required",
      "field": {"p": 13, "m": 1, "modulus": [0, 1]},
      "construction": "multiplicative",
      "n": 12,
      "k": 4,
      "localities": [2, 3],
      "local_distances": [2, 2],
      "expected_availability": 2,
      "expected_designed_distance": 12 - avail_exps[-1],
      "crt_exponents": avail_exps,
      "eval_points": avail_evals,
      "message": avail_msg,
      "expected_codeword": avail_cw,
      "test_coord": 3,
      "description": "Multiplicative CRT Availability LRC (12, 4, r=[2, 3]) over F_13",
  })

  # 4. Promoted v0.4 Interleaved RS (15, 7; l=3) over GF(16) with t=6 > t_u=4
  irs_n, irs_k, irs_ell = 15, 7, 3
  irs_red = irs_n - irs_k  # 8
  irs_msgs = [[rng.randrange(16) for _ in range(irs_k)] for _ in range(irs_ell)]
  if HAS_SAGE:
    r2 = sage.PolynomialRing(sage.GF(2), "x")
    x = r2.gen()
    f16 = sage.GF(16, name="a", modulus=x**4 + x + 1)
    s_evals = [f16.from_integer(e) for e in eval_points]
    grs7 = sage.codes.GeneralizedReedSolomonCode(s_evals, irs_k)
    enc7 = grs7.encoder("EvaluationPolynomial")
    pr16 = enc7.message_space()
    irs_cols = [
        [
            int(c.to_integer())
            for c in enc7.encode(pr16([f16.from_integer(v) for v in m]))
        ]
        for m in irs_msgs
    ]
  else:
    irs_cols = [
        [eval_poly_gf16(m, pt) for pt in eval_points] for m in irs_msgs
    ]
  irs_cw_mat = [[irs_cols[m][j] for m in range(irs_ell)] for j in range(irs_n)]

  # Dual multipliers w_j = (prod_{i != j}(alpha_j - alpha_i))^-1 in GF(16)
  irs_duals = []
  for j in range(irs_n):
    prod = 1
    for i in range(irs_n):
      if i != j:
        prod = gf2m_mul(prod, eval_points[j] ^ eval_points[i], gf16_mask, 4)
    irs_duals.append(gf2m_inv(prod, gf16_mask, 4))

  def _irs_hankel_full_rank(rx_mat: Sequence[Sequence[int]], tau: int) -> bool:
    seqs: List[List[int]] = []
    for m in range(irs_ell):
      seq_m: List[int] = []
      for i in range(irs_red):
        acc = 0
        for j in range(irs_n):
          w_ap = gf2m_mul(
              irs_duals[j], gf2m_pow(eval_points[j], i, gf16_mask, 4),
              gf16_mask, 4
          )
          acc ^= gf2m_mul(w_ap, rx_mat[j][m], gf16_mask, 4)
        seq_m.append(acc)
      seqs.append(seq_m)
    rows: List[List[int]] = []
    for seq_m in seqs:
      for i in range(irs_red - tau):
        row = [seq_m[i + tau - 1 - c] for c in range(tau)]
        row.append(seq_m[i + tau])
        rows.append(row)
    # RREF over GF(16)
    m_rows = len(rows)
    n_cols = tau + 1
    pivots: List[int] = []
    cur = 0
    for col in range(n_cols):
      if cur >= m_rows:
        break
      piv = -1
      for r_i in range(cur, m_rows):
        if rows[r_i][col] != 0:
          piv = r_i
          break
      if piv == -1:
        continue
      rows[cur], rows[piv] = rows[piv], rows[cur]
      inv_p = gf2m_inv(rows[cur][col], gf16_mask, 4)
      for c in range(col, n_cols):
        rows[cur][c] = gf2m_mul(rows[cur][c], inv_p, gf16_mask, 4)
      for r_i in range(m_rows):
        if r_i != cur and rows[r_i][col] != 0:
          fac = rows[r_i][col]
          for c in range(col, n_cols):
            rows[r_i][c] ^= gf2m_mul(fac, rows[cur][c], gf16_mask, 4)
      pivots.append(col)
      cur += 1
    return len(pivots) == tau and pivots[-1] < tau and rows[tau - 1][tau] != 0

  irs_err_positions = [0, 3, 5, 8, 11, 14]
  irs_received: List[List[int]] = []
  while True:
    cand_rx = [list(row) for row in irs_cw_mat]
    for pos in irs_err_positions:
      row_err = [rng.randrange(16) for _ in range(irs_ell)]
      if all(v == 0 for v in row_err):
        row_err[0] = 1
      for m in range(irs_ell):
        cand_rx[pos][m] ^= row_err[m]
    if _irs_hankel_full_rank(cand_rx, len(irs_err_positions)):
      irs_received = cand_rx
      break

  vectors.append({
      "id": "coding_interleaved_rs_15_7_l3",
      "category": "coding",
      "operation": "interleaved_rs_collaborative_decode",
      "status": "required",
      "field": {"p": 2, "m": 4, "modulus": gf16_mod},
      "n": irs_n,
      "k": irs_k,
      "interleaving_order": irs_ell,
      "unique_radius": 4,
      "collaborative_radius": 6,
      "eval_points": eval_points,
      "messages": irs_msgs,
      "expected_codeword_matrix": irs_cw_mat,
      "burst_error_positions": irs_err_positions,
      "received_matrix": irs_received,
      "description": "Interleaved RS (15, 7; l=3) over GF(16) collaborative burst decoding at t=6 > t_u=4",
  })

  # 5. Planned v0.4 placeholder (Gabidulin kept planned for Gatekeeper promotion)
  vectors.append({
      "id": "planned_v04_gabidulin_8_4_gf256",
      "category": "coding",
      "operation": "gabidulin_rank_decode",
      "status": "planned",
      "target_version": "v0.4",
      "field": {"p": 2, "m": 8, "modulus": [1, 0, 1, 1, 1, 0, 0, 0, 1]},
      "base_field": {"p": 2, "m": 1, "modulus": [0, 1]},
      "n": 8,
      "k": 4,
      "expected_min_rank_distance": 5,
      "unique_rank_radius": 2,
      "description": "Gabidulin MRD code G(8, 4) over GF(2^8)/GF(2) with Frobenius skew ring",
  })

  return vectors


def generate_golden_payload(seed: int) -> Dict[str, Any]:
  rng = random.Random(seed)
  vectors: List[Dict[str, Any]] = []
  vectors.extend(build_group_vectors())
  vectors.extend(build_ring_vectors(rng))
  vectors.extend(build_field_vectors(rng))
  vectors.extend(build_coding_vectors(rng))

  required_count = sum(1 for v in vectors if v["status"] == "required")
  planned_count = sum(1 for v in vectors if v["status"] == "planned")

  return {
      "schema_version": SCHEMA_VERSION,
      "generator": "sagemath" if HAS_SAGE else "python3_exact_reference",
      "seed": seed,
      "summary": {
          "total_vectors": len(vectors),
          "required_vectors": required_count,
          "planned_vectors": planned_count,
      },
      "vectors": vectors,
  }


def serialize_golden(payload: Dict[str, Any]) -> str:
  return json.dumps(payload, indent=2, sort_keys=True) + "\n"


def check_golden_dir(check_dir: pathlib.Path, seed: int) -> int:
  target = check_dir / GOLDEN_FILENAME
  if not target.is_file():
    print(f"DRIFT ERROR: missing golden file: {target}", file=sys.stderr)
    return 1

  expected_payload = generate_golden_payload(seed)
  expected_text = serialize_golden(expected_payload)
  actual_text = target.read_text(encoding="utf-8")

  if actual_text == expected_text:
    print(
        f"OK: {target} matches freshly generated golden vectors "
        f"(seed={seed}, generator={expected_payload['generator']})."
    )
    return 0

  # Allow cross-backend verification (e.g., SageMath checking a Python-generated
  # golden file or vice versa) ONLY if every single field other than `generator`
  # matches identically.
  try:
    actual_payload = json.loads(actual_text)
  except json.JSONDecodeError as exc:
    print(f"DRIFT ERROR: invalid JSON in {target}: {exc}", file=sys.stderr)
    return 1

  if actual_payload.get("generator") != expected_payload.get("generator"):
    actual_norm = dict(actual_payload)
    expected_norm = dict(expected_payload)
    actual_norm.pop("generator", None)
    expected_norm.pop("generator", None)
    if serialize_golden(actual_norm) == serialize_golden(expected_norm):
      print(
          f"OK (cross-backend exact match): {target} "
          f"(disk={actual_payload.get('generator')}, "
          f"runtime={expected_payload.get('generator')}, seed={seed})."
      )
      return 0

  print(
      f"DRIFT ERROR: {target} differs from generated golden output "
      f"(seed={seed}). Re-run with --out={target} to update.",
      file=sys.stderr,
  )
  return 1


def main(argv: Optional[Sequence[str]] = None) -> int:
  parser = argparse.ArgumentParser(
      description="Generate or drift-check xgalois v0.3 golden oracle vectors."
  )
  parser.add_argument(
      "--seed",
      type=int,
      default=2026,
      help="Deterministic PRNG seed (default: 2026).",
  )
  parser.add_argument(
      "--out",
      type=pathlib.Path,
      default=None,
      help="Output JSON file path.",
  )
  parser.add_argument(
      "--check-dir",
      type=pathlib.Path,
      default=None,
      help="Directory containing checked-in golden JSON files to drift-check.",
  )
  args = parser.parse_args(argv)

  if args.check_dir is not None:
    return check_golden_dir(args.check_dir, args.seed)

  payload = generate_golden_payload(args.seed)
  text = serialize_golden(payload)

  if args.out is not None:
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(text, encoding="utf-8")
    print(
        f"Wrote {args.out} "
        f"({payload['summary']['required_vectors']} required, "
        f"{payload['summary']['planned_vectors']} planned, "
        f"generator={payload['generator']}, seed={args.seed})."
    )
  else:
    sys.stdout.write(text)
  return 0


if __name__ == "__main__":
  sys.exit(main())
