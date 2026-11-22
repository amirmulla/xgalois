"""XGalois Python package (nanobind C++20 extension + SageMath interop & oracle)."""

from __future__ import annotations

from typing import Any

__version__ = "1.0.0"

_NATIVE_SYMBOLS = (
    "GFPrime",
    "GF2X",
    "TamoBargLRC",
    "InterleavedRS",
    "GRSCode",
    "GabidulinBinaryCode",
    "create_tamo_barg_lrc_prime",
    "create_tamo_barg_lrc_binary",
    "create_availability_lrc_prime",
    "create_availability_lrc_binary",
    "symmetric_group_order",
    "alternating_group_order",
    "dihedral_group_order",
    "negacyclic_mul",
    "crt_int",
)

_BUILD_HELP = (
    "The compiled 'xgalois._xgalois' C++20 extension is not available. "
    "Build and install it with `pip install .` (scikit-build-core) or "
    "`cmake -B build -DXGALOIS_BUILD_PYTHON=ON && cmake --build build`. "
    "Pure-Python modules (`xgalois.sage_interop`, `xgalois.oracle`) remain "
    "importable independently."
)

_NATIVE_IMPORT_ERROR: ImportError | None = None

try:
    from . import _xgalois as _xgalois  # type: ignore[attr-defined]
except ImportError as _exc:
    _xgalois = None
    _NATIVE_IMPORT_ERROR = _exc
    _HAS_NATIVE_EXTENSION = False
else:
    _HAS_NATIVE_EXTENSION = True
    GFPrime = _xgalois.GFPrime
    GF2X = _xgalois.GF2X
    TamoBargLRC = _xgalois.TamoBargLRC
    InterleavedRS = _xgalois.InterleavedRS
    GRSCode = _xgalois.GRSCode
    GabidulinBinaryCode = _xgalois.GabidulinBinaryCode
    create_tamo_barg_lrc_prime = _xgalois.create_tamo_barg_lrc_prime
    create_tamo_barg_lrc_binary = _xgalois.create_tamo_barg_lrc_binary
    create_availability_lrc_prime = _xgalois.create_availability_lrc_prime
    create_availability_lrc_binary = _xgalois.create_availability_lrc_binary
    symmetric_group_order = _xgalois.symmetric_group_order
    alternating_group_order = _xgalois.alternating_group_order
    dihedral_group_order = _xgalois.dihedral_group_order
    negacyclic_mul = _xgalois.negacyclic_mul
    crt_int = _xgalois.crt_int


def has_native_extension() -> bool:
    """Return True if the compiled C++20 `_xgalois` extension is loaded."""
    return _HAS_NATIVE_EXTENSION


def __getattr__(name: str) -> Any:
    if name in _NATIVE_SYMBOLS:
        raise ImportError(
            f"Cannot access 'xgalois.{name}': {_BUILD_HELP}"
        ) from _NATIVE_IMPORT_ERROR
    raise AttributeError(f"module 'xgalois' has no attribute {name!r}")


__all__ = [
    "__version__",
    "has_native_extension",
    *_NATIVE_SYMBOLS,
]
