"""Modulus-explicit interoperability helpers between XGalois and SageMath.

`sage.all` is imported lazily inside functions that construct SageMath objects,
so this module can be imported and tested in plain Python without Sage installed.
"""

from __future__ import annotations

from typing import Any, Iterable, List, Optional, Sequence


def modulus_mask_to_coeffs(mask: int, m: int) -> List[int]:
    """Convert a degree-m binary polynomial bitmask to ascending [c_0, ..., c_m].

    Bit i of `mask` is coefficient c_i in {0, 1}. Requires m >= 1 and bit m == 1
    (monic of degree m) with no higher bits set.
    """
    if not isinstance(m, int) or isinstance(m, bool) or m < 1:
        raise ValueError(f"Extension degree m must be an integer >= 1, got {m!r}")
    if not isinstance(mask, int) or isinstance(mask, bool) or mask < 0:
        raise ValueError(f"Modulus mask must be a non-negative integer, got {mask!r}")
    if mask >> (m + 1) != 0:
        raise ValueError(
            f"Modulus mask {mask:#x} has bits above degree {m}"
        )
    if ((mask >> m) & 1) == 0:
        raise ValueError(
            f"Modulus mask {mask:#x} is not monic of degree {m}"
        )
    return [(mask >> i) & 1 for i in range(m + 1)]


def coeffs_to_modulus_mask(coeffs: Sequence[int]) -> int:
    """Convert ascending binary polynomial coefficients [c_0, ..., c_m] to bitmask."""
    if len(coeffs) < 2:
        raise ValueError(
            "Binary modulus coefficient list must have length >= 2 (degree >= 1)"
        )
    mask = 0
    for deg, c in enumerate(coeffs):
        if isinstance(c, bool) or not isinstance(c, int) or c not in (0, 1):
            raise ValueError(
                f"Binary modulus coefficient at degree {deg} must be 0 or 1, got {c!r}"
            )
        if c == 1:
            mask |= 1 << deg
    if coeffs[-1] != 1:
        raise ValueError(
            f"Binary modulus polynomial must be monic (leading coeff == 1), got {coeffs[-1]!r}"
        )
    return mask


def _import_sage_gf() -> Any:
    """Lazily import `GF` from `sage.all`."""
    try:
        from sage.all import GF  # type: ignore[import-not-found]
    except ImportError as exc:
        raise ImportError(
            "SageMath ('sage.all') is required for SageMath conversion functions."
        ) from exc
    return GF


def to_sage_field(gf: Any) -> Any:
    """Construct the exact SageMath finite field matching an XGalois field.

    Preserves the exact defining irreducible polynomial for binary extension
    fields GF(2^m) so element bitmasks map identically to polynomial-basis
    representations in SageMath.
    """
    GF = _import_sage_gf()
    p = int(gf.characteristic)
    m = int(gf.degree)
    if m == 1:
        return GF(p)
    if p == 2:
        mask = int(gf.modulus)
        coeffs = modulus_mask_to_coeffs(mask, m)
        poly_ring = GF(2)["x"]
        mod_poly = poly_ring(coeffs)
        return GF(2**m, name="a", modulus=mod_poly)
    raise ValueError(
        f"Unsupported XGalois field for SageMath conversion: GF({p}^{m})"
    )


def to_sage_element(gf: Any, val: int, sage_field: Optional[Any] = None) -> Any:
    """Convert an integer XGalois field element to a SageMath field element."""
    field = sage_field if sage_field is not None else to_sage_field(gf)
    p = int(gf.characteristic)
    m = int(gf.degree)
    ival = int(val)
    if m == 1:
        return field(ival % p)
    if p == 2:
        if ival < 0 or ival >= (1 << m):
            raise ValueError(
                f"Element {ival} out of range for GF(2^{m}) (expected 0 <= val < {1 << m})"
            )
        bits = [(ival >> i) & 1 for i in range(m)]
        return field(bits)
    raise ValueError(
        f"Unsupported XGalois field for element conversion: GF({p}^{m})"
    )


def from_sage_element(elem: Any) -> int:
    """Convert a SageMath finite field element back to an XGalois integer value.

    For GF(2^m), returns the polynomial-basis bitmask `sum(c_i * 2^i)`.
    For GF(p), returns the canonical integer representative in `[0, p - 1]`.
    """
    if hasattr(elem, "to_integer"):
        return int(elem.to_integer())
    if hasattr(elem, "polynomial"):
        poly = elem.polynomial()
        coeffs = poly.list()
        mask = 0
        for i, c in enumerate(coeffs):
            if int(c) & 1:
                mask |= 1 << i
        return mask
    return int(elem)


def to_sage_poly(
    gf: Any, coeffs: Iterable[int], sage_field: Optional[Any] = None
) -> Any:
    """Convert ascending XGalois polynomial coefficients to a SageMath polynomial."""
    field = sage_field if sage_field is not None else to_sage_field(gf)
    poly_ring = field["x"]
    sage_coeffs = [to_sage_element(gf, c, sage_field=field) for c in coeffs]
    return poly_ring(sage_coeffs)


def from_sage_poly(poly: Any) -> List[int]:
    """Convert a SageMath univariate polynomial to ascending XGalois integer coeffs."""
    return [from_sage_element(c) for c in poly.list()]


__all__ = [
    "modulus_mask_to_coeffs",
    "coeffs_to_modulus_mask",
    "to_sage_field",
    "to_sage_element",
    "from_sage_element",
    "to_sage_poly",
    "from_sage_poly",
]
