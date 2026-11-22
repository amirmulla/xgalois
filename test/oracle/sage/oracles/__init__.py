#!/usr/bin/env python3
"""Shared exact finite-field reference arithmetic and suite helpers for live oracles.

Pure Python 3 (non-preparsed): automatically uses SageMath (`sage.all`) when
running under `sage -python`, and falls back to exact pure-Python GF(p)/GF(2^m)
arithmetic when running under standard `python3`.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import random
from typing import Any, Dict, List, Optional, Sequence, Tuple

try:
  import sage.all as sage  # type: ignore[import-not-found]

  HAS_SAGE = True
except ImportError:
  sage = None
  HAS_SAGE = False


@dataclass
class SuiteResult:
  """Pass/fail summary for a single differential oracle suite."""

  name: str
  passed: int = 0
  failed: int = 0
  failures: List[str] = field(default_factory=list)

  def record_pass(self) -> None:
    self.passed += 1

  def record_fail(self, msg: str) -> None:
    self.failed += 1
    self.failures.append(msg)

  @property
  def ok(self) -> bool:
    return self.failed == 0 and self.passed > 0


def poly_ascending_to_mask(coeffs: Sequence[int]) -> int:
  mask = 0
  for deg, c in enumerate(coeffs):
    if c & 1:
      mask |= 1 << deg
  return mask


class FiniteFieldRef:
  """Exact finite field GF(p) or GF(2^m) with explicit modulus.

  Uses SageMath when available and exact pure-Python arithmetic otherwise.
  Elements are canonically represented as Python `int`s in [0, q - 1].
  """

  def __init__(
      self,
      p: int,
      m: int = 1,
      modulus: Optional[Sequence[int]] = None,
  ) -> None:
    self.p = int(p)
    self.m = int(m)
    if self.m == 1:
      self.modulus = [0, 1] if modulus is None else [int(c) for c in modulus]
      self.q = self.p
      self._mask = 0
    elif self.p == 2 and 2 <= self.m <= 16:
      if modulus is None or len(modulus) != self.m + 1:
        raise ValueError(
            f"GF(2^{self.m}) requires explicit ascending modulus of length"
            f" {self.m + 1}"
        )
      self.modulus = [int(c) & 1 for c in modulus]
      if self.modulus[self.m] != 1:
        raise ValueError("Binary extension modulus must be monic of degree m")
      self.q = 1 << self.m
      self._mask = poly_ascending_to_mask(self.modulus)
    else:
      raise ValueError(f"Unsupported field GF({self.p}^{self.m})")

    self._sage_field: Any = None
    if HAS_SAGE:
      if self.m == 1:
        self._sage_field = sage.GF(self.p)
      else:
        r2 = sage.PolynomialRing(sage.GF(2), "x")
        x = r2.gen()
        mod_poly = sum(c * x**i for i, c in enumerate(self.modulus))
        self._sage_field = sage.GF(self.q, name="a", modulus=mod_poly)

  def spec(self) -> Dict[str, Any]:
    return {"p": self.p, "m": self.m, "modulus": list(self.modulus)}

  @property
  def sage_field(self) -> Any:
    return self._sage_field

  def to_sage(self, a: int) -> Any:
    if not HAS_SAGE:
      raise RuntimeError("SageMath is not available")
    if self.m == 1:
      return self._sage_field(a % self.p)
    if hasattr(self._sage_field, "from_integer"):
      return self._sage_field.from_integer(int(a))
    return self._sage_field.fetch_int(int(a))

  def from_sage(self, elem: Any) -> int:
    if self.m == 1:
      return int(elem) % self.p
    if hasattr(elem, "to_integer"):
      return int(elem.to_integer())
    return int(elem.integer_representation())

  def add(self, a: int, b: int) -> int:
    if self.m == 1:
      return (a + b) % self.p
    return a ^ b

  def sub(self, a: int, b: int) -> int:
    if self.m == 1:
      return (a - b) % self.p
    return a ^ b

  def neg(self, a: int) -> int:
    if self.m == 1:
      return (-a) % self.p
    return a

  def mul(self, a: int, b: int) -> int:
    if self.m == 1:
      return (a * b) % self.p
    res = 0
    aa = a
    bb = b
    top = 1 << self.m
    mask = self._mask
    while bb:
      if bb & 1:
        res ^= aa
      aa <<= 1
      if aa & top:
        aa ^= mask
      bb >>= 1
    return res

  def pow(self, base: int, exp: int) -> int:
    if self.m == 1:
      return pow(base, exp, self.p)
    res = 1
    b = base
    e = exp
    while e > 0:
      if e & 1:
        res = self.mul(res, b)
      b = self.mul(b, b)
      e >>= 1
    return res

  def inv(self, a: int) -> int:
    if a == 0:
      raise ZeroDivisionError("Zero has no inverse")
    if self.m == 1:
      return pow(a, self.p - 2, self.p)
    return self.pow(a, self.q - 2)

  def div(self, a: int, b: int) -> int:
    return self.mul(a, self.inv(b))

  def eval_poly(self, coeffs: Sequence[int], x: int) -> int:
    acc = 0
    for c in reversed(coeffs):
      acc = self.add(self.mul(acc, x), c)
    return acc

  def lagrange_eval(
      self, xs: Sequence[int], ys: Sequence[int], target_x: int
  ) -> int:
    r = len(xs)
    assert len(ys) == r
    acc = 0
    for a in range(r):
      num = 1
      den = 1
      for b in range(r):
        if b != a:
          num = self.mul(num, self.sub(target_x, xs[b]))
          den = self.mul(den, self.sub(xs[a], xs[b]))
      acc = self.add(acc, self.mul(ys[a], self.div(num, den)))
    return acc

  def rref_with_pivots(
      self, rows: Sequence[Sequence[int]]
  ) -> Tuple[List[List[int]], List[int]]:
    if not rows:
      return [], []
    m_rows = len(rows)
    n_cols = len(rows[0])
    mat = [list(r) for r in rows]
    pivots: List[int] = []
    cur_row = 0
    for col in range(n_cols):
      if cur_row >= m_rows:
        break
      piv = -1
      for r in range(cur_row, m_rows):
        if mat[r][col] != 0:
          piv = r
          break
      if piv == -1:
        continue
      if piv != cur_row:
        mat[cur_row], mat[piv] = mat[piv], mat[cur_row]
      inv_p = self.inv(mat[cur_row][col])
      if inv_p != 1:
        for c in range(col, n_cols):
          mat[cur_row][c] = self.mul(mat[cur_row][c], inv_p)
      for r in range(m_rows):
        if r == cur_row:
          continue
        factor = mat[r][col]
        if factor == 0:
          continue
        for c in range(col, n_cols):
          mat[r][c] = self.sub(mat[r][c], self.mul(factor, mat[cur_row][c]))
      pivots.append(col)
      cur_row += 1
    return mat, pivots

  def primitive_element(self) -> int:
    order = self.q - 1
    factors: List[int] = []
    tmp = order
    d = 2
    while d * d <= tmp:
      if tmp % d == 0:
        factors.append(d)
        while tmp % d == 0:
          tmp //= d
      d += 1
    if tmp > 1:
      factors.append(tmp)
    for cand in range(2, self.q):
      if all(self.pow(cand, order // f) != 1 for f in factors):
        return cand
    return 1

  def primitive_eval_points(self, n: int) -> List[int]:
    if n > self.q - 1:
      raise ValueError(f"n={n} exceeds q-1={self.q - 1}")
    g = self.primitive_element()
    pts = []
    cur = 1
    for _ in range(n):
      pts.append(cur)
      cur = self.mul(cur, g)
    return pts

  def random_element(self, rng: random.Random) -> int:
    return rng.randrange(self.q)

  def random_nonzero(self, rng: random.Random) -> int:
    return rng.randrange(1, self.q)


# Standard explicit-modulus fields shared across live oracle suites
GF8_MOD = [1, 1, 0, 1]  # x^3 + x + 1
GF16_MOD = [1, 1, 0, 0, 1]  # x^4 + x + 1
GF32_MOD = [1, 0, 1, 0, 0, 1]  # x^5 + x^2 + 1
GF64_MOD = [1, 1, 0, 0, 0, 0, 1]  # x^6 + x + 1
GF256_MOD = [1, 0, 1, 1, 1, 0, 0, 0, 1]  # x^8 + x^4 + x^3 + x^2 + 1 (0x11D)

__all__ = [
    "FiniteFieldRef",
    "GF8_MOD",
    "GF16_MOD",
    "GF32_MOD",
    "GF64_MOD",
    "GF256_MOD",
    "HAS_SAGE",
    "SuiteResult",
    "poly_ascending_to_mask",
    "sage",
]
