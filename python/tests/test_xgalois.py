"""Unit tests for XGalois Python bindings, SageMath interop, and oracle layer."""

from __future__ import annotations

import os
import pathlib
import stat
import sys
import tempfile
import types
import unittest
from typing import Any, Dict, List, Optional, Tuple
from unittest import mock

# Ensure `python/` is on sys.path when invoked via
# `python3 -m unittest discover -s python/tests -v` from repo root.
_PYTHON_DIR = pathlib.Path(__file__).resolve().parents[1]
if str(_PYTHON_DIR) not in sys.path:
    sys.path.insert(0, str(_PYTHON_DIR))

import xgalois  # noqa: E402
from xgalois import oracle, sage_interop  # noqa: E402

try:
    import numpy as np  # type: ignore[import-not-found]

    _HAS_NUMPY = True
except ImportError:
    np = None
    _HAS_NUMPY = False


# ============================================================================
# Helpers / Mocks for pure-Python testing
# ============================================================================


class _MockPrimeField:
    def __init__(self, p: int) -> None:
        self.characteristic = p
        self.degree = 1
        self.order = p
        self.modulus = p

    def add(self, a: int, b: int) -> int:
        return (a + b) % self.characteristic


class _MockBinaryExtField:
    def __init__(self, m: int, modulus_mask: int) -> None:
        self.characteristic = 2
        self.degree = m
        self.order = 1 << m
        self.modulus = modulus_mask

    def add(self, a: int, b: int) -> int:
        return a ^ b


class _MockSagePoly:
    def __init__(self, coeffs: List[Any]) -> None:
        self._coeffs = list(coeffs)

    def list(self) -> List[Any]:
        return list(self._coeffs)


class _MockSagePolyRing:
    def __init__(self, base_field: Any) -> None:
        self.base_field = base_field

    def __call__(self, coeffs: List[Any]) -> _MockSagePoly:
        return _MockSagePoly(coeffs)


class _MockSageElement:
    def __init__(self, val: int, is_binary_ext: bool) -> None:
        self._val = val
        self._is_binary_ext = is_binary_ext

    def to_integer(self) -> int:
        return self._val

    def __int__(self) -> int:
        return self._val


class _MockSageField:
    def __init__(
        self, order: int, *, name: Optional[str] = None, modulus: Any = None
    ) -> None:
        self.order = order
        self.name = name
        self.modulus = modulus
        self._is_binary_ext = modulus is not None

    def __getitem__(self, var: str) -> _MockSagePolyRing:
        del var
        return _MockSagePolyRing(self)

    def __call__(self, val: Any) -> _MockSageElement:
        if isinstance(val, list):
            mask = 0
            for i, bit in enumerate(val):
                if int(bit) & 1:
                    mask |= 1 << i
            return _MockSageElement(mask, True)
        return _MockSageElement(int(val) % self.order, self._is_binary_ext)


def _mock_sage_gf(order: int, name: Optional[str] = None, modulus: Any = None) -> _MockSageField:
    return _MockSageField(order, name=name, modulus=modulus)


class _MockNativeModule:
    """In-memory stub implementing the `_xgalois` C++ extension surface."""

    class GFPrime:
        def __init__(self, p: int) -> None:
            self.characteristic = p
            self.degree = 1
            self.order = p

        def add(self, a: int, b: int) -> int:
            return (a + b) % self.characteristic

    class GF2X:
        def __init__(self, m: int, modulus: Any = "", check: bool = True) -> None:
            del check
            self.characteristic = 2
            self.degree = m
            self.order = 1 << m
            self.modulus = (
                sage_interop.coeffs_to_modulus_mask(modulus)
                if isinstance(modulus, list)
                else 0
            )

        def add(self, a: int, b: int) -> int:
            return a ^ b

    class _MockLRC:
        def __init__(self, n: int, k: int, r: int, delta: int) -> None:
            self.length = n
            self.dimension = k
            self.locality = r
            self.local_distance = delta
            self.availability = 1
            self.minimum_distance = n - k
            self.designed_distance = n - k

        def encode(self, msg: List[int]) -> List[int]:
            return list(msg) + [0] * (self.length - len(msg))

        def unencode(self, cw: List[int]) -> List[int]:
            return list(cw[: self.dimension])

        def evaluation_points(self) -> List[int]:
            return list(range(1, self.length + 1))

        def repair_sets(self, coord: int) -> List[List[int]]:
            return [[(coord + 1) % self.length, (coord + 2) % self.length]]

        def fiber_partition(self, idx: int = 0) -> List[List[int]]:
            del idx
            return [list(range(self.length))]

        def good_poly(self) -> List[int]:
            return [0, 0, 0, 1]

        def local_repair(
            self, received: List[int], erased_index: int, repair_set_idx: int = 0
        ) -> Dict[str, Any]:
            del received
            return {
                "success": True,
                "recovered_symbol": 42,
                "erased_index": erased_index,
                "repair_set_index": repair_set_idx,
                "symbols_read": self.locality,
                "error_message": "",
            }

        def local_repair_multiple(
            self,
            received: List[int],
            erased_positions: List[int],
            allow_global_fallback: bool = True,
        ) -> Dict[str, Any]:
            del erased_positions, allow_global_fallback
            return {
                "success": True,
                "codeword": list(received),
                "message": list(received[: self.dimension]),
                "num_local_repairs": 2,
                "num_global_repairs": 0,
                "total_symbols_read": 2 * self.locality,
                "used_global_fallback": False,
                "error_message": "",
            }

    @classmethod
    def create_tamo_barg_lrc_prime(
        cls,
        gf: Any,
        n: int,
        k: int,
        r: int,
        delta: int = 2,
        partition_type: Optional[str] = None,
    ) -> Any:
        del gf, partition_type
        return cls._MockLRC(n, k, r, delta)

    @classmethod
    def create_tamo_barg_lrc_binary(
        cls,
        gf: Any,
        n: int,
        k: int,
        r: int,
        delta: int = 2,
        partition_type: Optional[str] = None,
    ) -> Any:
        del gf, partition_type
        return cls._MockLRC(n, k, r, delta)

    class InterleavedRS:
        def __init__(self, gf: Any, eval_pts: List[int], *args: Any) -> None:
            del gf, eval_pts, args

        def decode(
            self,
            received: List[List[int]],
            mode: str = "ssb",
            fallback: bool = True,
            erasure_positions: Optional[List[int]] = None,
        ) -> Dict[str, Any]:
            del mode, fallback, erasure_positions
            return {
                "success": True,
                "codeword_matrix": received,
                "messages": [[1, 2]],
                "burst_error_positions": [],
                "collaborative_radius": 2,
                "unique_radius": 1,
                "used_fallback": False,
                "error_message": "",
            }

    class GRSCode:
        def __init__(
            self,
            gf: Any,
            eval_pts: List[int],
            k: int,
            col_mults: Optional[List[int]] = None,
        ) -> None:
            del gf, eval_pts, col_mults
            self.k = k

        def decode_power(
            self,
            received: List[int],
            powering_degree: Optional[int] = None,
            mode: str = "ssb",
            fallback_to_bm: bool = True,
            erasure_positions: Optional[List[int]] = None,
        ) -> Dict[str, Any]:
            del mode, fallback_to_bm, erasure_positions
            return {
                "success": True,
                "codeword": received,
                "message": received[: self.k],
                "num_errors": 0,
                "decoding_radius": 3,
                "powering_degree": powering_degree or 2,
                "error_message": "",
            }

        def decode_guruswami_sudan(
            self, received: List[int], tau: int
        ) -> Dict[str, Any]:
            return {
                "codewords": [received],
                "messages": [received[: self.k]],
                "decoding_radius": tau,
                "johnson_radius": tau,
                "multiplicity": 1,
                "list_size": 2,
            }

    class GabidulinBinaryCode:
        def __init__(self, gf: Any, n: int, k: int) -> None:
            del gf
            self.n = n
            self.k = k

        def decode(
            self, received: List[int], decoder: str = "GabidulinSkewEEA"
        ) -> Dict[str, Any]:
            del decoder
            return {
                "success": True,
                "codeword": received,
                "message": received[: self.k],
                "minimum_distance": self.n - self.k + 1,
                "error_rank": 0,
                "error_message": "",
            }

    @staticmethod
    def symmetric_group_order(n: int) -> str:
        acc = 1
        for i in range(2, n + 1):
            acc *= i
        return str(acc)

    @staticmethod
    def alternating_group_order(n: int) -> str:
        if n < 2:
            return "1"
        acc = 1
        for i in range(2, n + 1):
            acc *= i
        return str(acc // 2)

    @staticmethod
    def dihedral_group_order(n: int) -> str:
        return str(2 * n)

    @staticmethod
    def negacyclic_mul(q: int, n: int, a: List[int], b: List[int]) -> List[int]:
        out = [0] * n
        for i, av in enumerate(a):
            for j, bv in enumerate(b):
                idx = i + j
                if idx < n:
                    out[idx] = (out[idx] + av * bv) % q
                else:
                    out[idx - n] = (out[idx - n] - av * bv) % q
        return out

    @staticmethod
    def crt_int(
        residues: List[int], moduli: List[int]
    ) -> Optional[Tuple[str, str]]:
        if residues == [1, 2] and moduli == [2, 2]:
            return None
        return ("23", "105")


# ============================================================================
# Tests: Package __init__
# ============================================================================


class TestPackageInit(unittest.TestCase):
    def test_version_and_pure_python_imports(self) -> None:
        self.assertEqual(xgalois.__version__, "1.0.0")
        self.assertTrue(callable(xgalois.has_native_extension))
        self.assertTrue(hasattr(sage_interop, "modulus_mask_to_coeffs"))
        self.assertTrue(hasattr(oracle, "NanobindTransport"))

    def test_missing_native_extension_raises_helpful_import_error(self) -> None:
        if xgalois.has_native_extension():
            self.skipTest("Native extension is built")
        with self.assertRaises(ImportError) as ctx:
            _ = xgalois.GFPrime
        msg = str(ctx.exception)
        self.assertIn("pip install .", msg)
        self.assertIn("XGALOIS_BUILD_PYTHON=ON", msg)


# ============================================================================
# Tests: sage_interop
# ============================================================================


class TestSageInterop(unittest.TestCase):
    def test_modulus_mask_coeffs_roundtrip(self) -> None:
        cases = [
            (2, 0b111, [1, 1, 1]),  # x^2 + x + 1
            (3, 0b1011, [1, 1, 0, 1]),  # x^3 + x + 1
            (4, 0b10011, [1, 1, 0, 0, 1]),  # x^4 + x + 1
            (8, 0x11B, [1, 1, 0, 1, 1, 0, 0, 0, 1]),  # AES polynomial
            (16, 0x1100B, [1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1]),
        ]
        for m, mask, expected_coeffs in cases:
            coeffs = sage_interop.modulus_mask_to_coeffs(mask, m)
            self.assertEqual(coeffs, expected_coeffs)
            self.assertEqual(sage_interop.coeffs_to_modulus_mask(coeffs), mask)

    def test_modulus_validation_errors(self) -> None:
        with self.assertRaises(ValueError):
            sage_interop.modulus_mask_to_coeffs(0b011, 3)  # non-monic
        with self.assertRaises(ValueError):
            sage_interop.modulus_mask_to_coeffs(0b11011, 3)  # bits above m
        with self.assertRaises(ValueError):
            sage_interop.modulus_mask_to_coeffs(0b11, 0)  # m < 1
        with self.assertRaises(ValueError):
            sage_interop.coeffs_to_modulus_mask([1])  # too short
        with self.assertRaises(ValueError):
            sage_interop.coeffs_to_modulus_mask([1, 2, 1])  # non-binary
        with self.assertRaises(ValueError):
            sage_interop.coeffs_to_modulus_mask([1, 1, 0])  # non-monic

    def test_sage_field_element_poly_roundtrip_with_mock_sage(self) -> None:
        fake_sage_all = types.SimpleNamespace(GF=_mock_sage_gf)
        with mock.patch.dict(sys.modules, {"sage": types.SimpleNamespace(), "sage.all": fake_sage_all}):
            # Prime field GF(17)
            gf_p = _MockPrimeField(17)
            sf_p = sage_interop.to_sage_field(gf_p)
            self.assertEqual(sf_p.order, 17)
            el_p = sage_interop.to_sage_element(gf_p, 13, sage_field=sf_p)
            self.assertEqual(sage_interop.from_sage_element(el_p), 13)
            poly_p = sage_interop.to_sage_poly(gf_p, [3, 0, 16], sage_field=sf_p)
            self.assertEqual(sage_interop.from_sage_poly(poly_p), [3, 0, 16])

            # Binary extension field GF(2^4) with x^4 + x + 1 (0b10011)
            gf_2 = _MockBinaryExtField(4, 0b10011)
            sf_2 = sage_interop.to_sage_field(gf_2)
            self.assertEqual(sf_2.order, 16)
            self.assertEqual(sf_2.modulus.list(), [1, 1, 0, 0, 1])
            for val in range(16):
                el_2 = sage_interop.to_sage_element(gf_2, val, sage_field=sf_2)
                self.assertEqual(sage_interop.from_sage_element(el_2), val)
            poly_2 = sage_interop.to_sage_poly(gf_2, [5, 11, 0, 15], sage_field=sf_2)
            self.assertEqual(sage_interop.from_sage_poly(poly_2), [5, 11, 0, 15])


# ============================================================================
# Tests: oracle.py (SubprocessNdjsonTransport & NanobindTransport)
# ============================================================================


class TestOracleTransports(unittest.TestCase):
    def test_subprocess_ndjson_transport(self) -> None:
        script_body = (
            "#!/usr/bin/env python3\n"
            "import json, sys\n"
            "for line in sys.stdin:\n"
            "    req = json.loads(line)\n"
            "    method = req['method']\n"
            "    if method == 'ping':\n"
            "        resp = {'id': req['id'], 'ok': True, 'result': {'status': 'ok', 'version': 'mock'}}\n"
            "    elif method == 'fail_me':\n"
            "        resp = {'id': req['id'], 'ok': False, 'error': 'boom'}\n"
            "    else:\n"
            "        resp = {'id': req['id'], 'ok': True, 'result': {'echo': req['params']}}\n"
            "    sys.stdout.write(json.dumps(resp) + '\\n')\n"
            "    sys.stdout.flush()\n"
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            cli = pathlib.Path(tmpdir) / "mock_oracle_cli.py"
            cli.write_text(script_body, encoding="utf-8")
            os.chmod(cli, os.stat(cli).st_mode | stat.S_IXUSR)

            with oracle.OracleClient.from_cli(cli, timeout=5.0) as client:
                self.assertEqual(client.ping()["status"], "ok")
                echo = client.call("echo", {"x": 123})
                self.assertEqual(echo["echo"], {"x": 123})
                with self.assertRaises(oracle.OracleRpcError) as ctx:
                    client.call("fail_me", {})
                self.assertIn("boom", str(ctx.exception))

        with self.assertRaises(FileNotFoundError):
            oracle.SubprocessNdjsonTransport("/nonexistent/xgalois_oracle_cli")

    def test_nanobind_transport_all_methods_with_mock_module(self) -> None:
        with oracle.OracleClient.from_nanobind(module=_MockNativeModule) as client:
            ping = client.ping()
            self.assertEqual(ping["status"], "ok")
            self.assertEqual(ping["transport"], "nanobind")

            lrc_res = client.call(
                "lrc_encode_and_repair",
                {
                    "field": {"p": 13, "m": 1},
                    "n": 6,
                    "k": 4,
                    "r": 2,
                    "delta": 2,
                    "message": [1, 2, 3, 4],
                    "erased_positions": [0],
                },
            )
            self.assertTrue(lrc_res["repair_success"])
            self.assertEqual(lrc_res["repaired_codeword"][0], 42)
            self.assertEqual(lrc_res["symbols_read"], 2)

            irs_res = client.call(
                "irs_decode",
                {
                    "field": {"p": 17, "m": 1},
                    "evaluation_points": [1, 2, 3, 4, 5, 6],
                    "k": 2,
                    "ell": 2,
                    "received_matrix": [[1, 2]] * 6,
                },
            )
            self.assertTrue(irs_res["success"])

            pwr_res = client.call(
                "grs_power_decode",
                {
                    "field": {"p": 17, "m": 1},
                    "evaluation_points": [1, 2, 3, 4, 5, 6],
                    "k": 2,
                    "received": [1, 2, 3, 4, 5, 6],
                },
            )
            self.assertTrue(pwr_res["success"])

            gs_res = client.call(
                "grs_guruswami_sudan",
                {
                    "field": {"p": 17, "m": 1},
                    "evaluation_points": [1, 2, 3, 4, 5, 6],
                    "k": 2,
                    "received": [1, 2, 3, 4, 5, 6],
                    "tau": 3,
                },
            )
            self.assertEqual(gs_res["decoding_radius"], 3)

            gab_res = client.call(
                "gabidulin_decode",
                {
                    "field": {"p": 2, "m": 4, "modulus": [1, 1, 0, 0, 1]},
                    "n": 4,
                    "k": 2,
                    "received": [1, 2, 4, 8],
                },
            )
            self.assertTrue(gab_res["success"])
            self.assertEqual(gab_res["minimum_distance"], 3)

            self.assertEqual(
                client.call(
                    "groups_and_rings_check",
                    {"sub_op": "symmetric_group", "degree": 5},
                )["order"],
                "120",
            )
            self.assertEqual(
                client.call(
                    "groups_and_rings_check",
                    {"sub_op": "alternating_group", "degree": 5},
                )["order"],
                "60",
            )
            self.assertEqual(
                client.call(
                    "groups_and_rings_check",
                    {"sub_op": "dihedral_group", "n": 6},
                )["order"],
                "12",
            )
            self.assertEqual(
                client.call(
                    "groups_and_rings_check",
                    {
                        "sub_op": "negacyclic_mul",
                        "modulus_q": 17,
                        "degree": 2,
                        "a": [0, 1],
                        "b": [0, 1],
                    },
                )["product"],
                [16, 0],
            )
            crt_ok = client.call(
                "groups_and_rings_check",
                {"sub_op": "crt_int", "residues": [2, 3, 2], "moduli": [3, 5, 7]},
            )
            self.assertTrue(crt_ok["consistent"])
            self.assertEqual(crt_ok["solution"], "23")
            crt_bad = client.call(
                "groups_and_rings_check",
                {"sub_op": "crt_int", "residues": [1, 2], "moduli": [2, 2]},
            )
            self.assertFalse(crt_bad["consistent"])

            with self.assertRaises(oracle.OracleRpcError):
                client.call("unknown_method", {})


# ============================================================================
# Tests: Compiled C++20 `_xgalois` Nanobind Bindings
# ============================================================================


@unittest.skipUnless(
    xgalois.has_native_extension(),
    "Compiled _xgalois C++20 extension not built; skipping native binding tests",
)
class TestNativeBindings(unittest.TestCase):
    def test_gf_prime_scalar_and_numpy(self) -> None:
        gf = xgalois.GFPrime(17)
        self.assertEqual(gf.order, 17)
        self.assertEqual(gf.characteristic, 17)
        self.assertEqual(gf.degree, 1)
        self.assertEqual(gf.add(12, 9), 4)
        self.assertEqual(gf.sub(3, 10), 10)
        self.assertEqual(gf.mul(6, 7), 8)
        self.assertEqual(gf.mul(5, gf.inv(5)), 1)
        self.assertEqual(gf.div(8, 6), 7)
        self.assertEqual(gf.pow(3, 4), 13)

        if _HAS_NUMPY:
            a = np.array([1, 5, 10, 16], dtype=np.uint32)
            b = np.array([2, 13, 8, 2], dtype=np.uint32)
            dst = np.zeros(4, dtype=np.uint32)
            gf.add_array(dst, a, b)
            self.assertEqual(dst.tolist(), [3, 1, 1, 1])
            gf.mul_array(dst, a, b)
            self.assertEqual(dst.tolist(), [2, 14, 12, 15])
            gf.fma_array(dst, a, 3)
            self.assertEqual(dst.tolist(), [5, 12, 8, 12])

    def test_gf2x_scalar_and_numpy(self) -> None:
        gf = xgalois.GF2X(4, [1, 1, 0, 0, 1])  # x^4 + x + 1 = 0b10011 = 19
        self.assertEqual(gf.order, 16)
        self.assertEqual(gf.characteristic, 2)
        self.assertEqual(gf.degree, 4)
        self.assertEqual(gf.modulus, 0b10011)
        self.assertEqual(gf.add(5, 3), 6)
        self.assertEqual(gf.sub(5, 3), 6)
        self.assertEqual(gf.mul(7, gf.inv(7)), 1)
        self.assertEqual(gf.div(gf.mul(5, 9), 9), 5)
        self.assertEqual(gf.pow(2, 4), 3)  # x^4 == x + 1

        gf_str = xgalois.GF2X(8, "x^8 + x^4 + x^3 + x + 1")
        self.assertEqual(gf_str.modulus, 0x11B)

        if _HAS_NUMPY:
            a = np.array([1, 2, 4, 8], dtype=np.uint32)
            b = np.array([3, 5, 7, 9], dtype=np.uint32)
            dst = np.zeros(4, dtype=np.uint32)
            gf.add_array(dst, a, b)
            self.assertEqual(dst.tolist(), [2, 7, 3, 1])
            gf.mul_array(dst, a, b)
            self.assertEqual(dst.tolist(), [gf.mul(int(x), int(y)) for x, y in zip(a, b)])
            prev = dst.copy()
            gf.fma_array(dst, a, 5)
            self.assertEqual(
                dst.tolist(),
                [gf.add(int(p), gf.mul(int(x), 5)) for p, x in zip(prev, a)],
            )

    def test_tamo_barg_lrc(self) -> None:
        gf = xgalois.GFPrime(13)
        lrc = xgalois.create_tamo_barg_lrc_prime(gf, 6, 4, 2, 2)
        self.assertEqual(lrc.length, 6)
        self.assertEqual(lrc.dimension, 4)
        self.assertEqual(lrc.locality, 2)
        self.assertEqual(lrc.local_distance, 2)
        self.assertEqual(lrc.availability, 1)
        self.assertGreaterEqual(lrc.minimum_distance, 2)
        self.assertEqual(len(lrc.repair_sets(0)[0]), 2)

        msg = [1, 4, 7, 11]
        cw = lrc.encode(msg)
        corrupted = list(cw)
        corrupted[0] = gf.add(corrupted[0], 1)
        rep = lrc.local_repair(corrupted, 0, 0)
        self.assertTrue(rep["success"])
        self.assertEqual(rep["recovered_symbol"], cw[0])
        self.assertEqual(rep["symbols_read"], 2)

        corrupted2 = list(cw)
        corrupted2[0] = gf.add(corrupted2[0], 1)
        corrupted2[3] = gf.add(corrupted2[3], 1)
        rep_multi = lrc.local_repair_multiple(corrupted2, [0, 3], True)
        self.assertTrue(rep_multi["success"])
        self.assertEqual(rep_multi["codeword"], cw)

        # Binary LRC
        gf2 = xgalois.GF2X(4, [1, 1, 0, 0, 1])
        lrc_bin = xgalois.create_tamo_barg_lrc_binary(gf2, 6, 4, 2, 2)
        cw_bin = lrc_bin.encode([1, 2, 3, 4])
        cor_bin = list(cw_bin)
        cor_bin[1] = gf2.add(cor_bin[1], 1)
        rep_bin = lrc_bin.local_repair(cor_bin, 1, 0)
        self.assertTrue(rep_bin["success"])
        self.assertEqual(rep_bin["recovered_symbol"], cw_bin[1])

    def test_interleaved_rs_and_grs(self) -> None:
        gf = xgalois.GFPrime(17)
        eval_pts = list(range(1, 13))  # n = 12
        irs = xgalois.InterleavedRS(gf, eval_pts, 4, 3)  # k=4, ell=3 -> r=8, t_collab=6
        self.assertEqual(irs.length, 12)
        self.assertEqual(irs.interleaving_order, 3)
        self.assertEqual(irs.unique_radius, 4)
        self.assertEqual(irs.collaborative_radius, 6)

        msgs = [[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12]]
        cw_mat = irs.encode(msgs)
        recv = [list(row) for row in cw_mat]
        # Inject 5 burst errors (beyond unique radius 4, within collaborative radius 6)
        for pos, err in zip([0, 2, 5, 8, 11], [(1, 2, 3), (4, 1, 5), (2, 7, 1), (3, 3, 8), (6, 2, 4)]):
            for col in range(3):
                recv[pos][col] = gf.add(recv[pos][col], err[col])
        dec = irs.decode(recv, mode="ssb", fallback=True)
        self.assertTrue(dec["success"])
        self.assertEqual(dec["codeword_matrix"], cw_mat)
        self.assertEqual(dec["messages"], msgs)

        # GRS BM, Power, Guruswami-Sudan
        grs = xgalois.GRSCode(gf, eval_pts, 4)
        cw = grs.encode([1, 3, 5, 7])
        cor = list(cw)
        for pos in (0, 3, 7):
            cor[pos] = gf.add(cor[pos], 1)
        bm_res = grs.decode_bm(cor)
        self.assertTrue(bm_res["success"])
        self.assertEqual(bm_res["codeword"], cw)

        pwr_res = grs.decode_power(cor)
        self.assertTrue(pwr_res["success"])
        self.assertEqual(pwr_res["codeword"], cw)

        gs_res = grs.decode_guruswami_sudan(cor, 5)
        self.assertIn(cw, gs_res["codewords"])

    def test_gabidulin_binary_code(self) -> None:
        gf2 = xgalois.GF2X(4, [1, 1, 0, 0, 1])
        gab = xgalois.GabidulinBinaryCode(gf2, 4, 2)
        self.assertEqual(gab.length, 4)
        self.assertEqual(gab.dimension, 2)
        self.assertEqual(gab.minimum_distance, 3)
        cw = gab.encode([3, 5])
        # Inject rank-1 error: [1, 1, 0, 0]
        cor = [gf2.add(cw[0], 1), gf2.add(cw[1], 1), cw[2], cw[3]]
        self.assertEqual(gab.rank_weight([1, 1, 0, 0]), 1)
        dec = gab.decode(cor, decoder="GabidulinSkewEEA")
        self.assertTrue(dec["success"])
        self.assertEqual(dec["codeword"], cw)
        self.assertEqual(dec["error_rank"], 1)

    def test_algebra_and_rings(self) -> None:
        self.assertEqual(xgalois.symmetric_group_order(5), "120")
        self.assertEqual(xgalois.alternating_group_order(5), "60")
        self.assertEqual(xgalois.dihedral_group_order(6), "12")
        self.assertEqual(xgalois.negacyclic_mul(17, 2, [0, 1], [0, 1]), [16, 0])
        self.assertEqual(xgalois.crt_int([2, 3, 2], [3, 5, 7]), ("23", "105"))
        self.assertIsNone(xgalois.crt_int([1, 2], [2, 2]))

    def test_live_nanobind_oracle_transport(self) -> None:
        with oracle.OracleClient.from_nanobind() as client:
            self.assertEqual(client.ping()["status"], "ok")
            res = client.call(
                "groups_and_rings_check",
                {"sub_op": "symmetric_group", "degree": 6},
            )
            self.assertEqual(res["order"], "720")


if __name__ == "__main__":
    unittest.main()
