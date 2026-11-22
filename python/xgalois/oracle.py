"""Unified dual-mode oracle client for XGalois.

Supports both:
- `SubprocessNdjsonTransport`: persistent NDJSON subprocess speaking to
  `xgalois_oracle_cli` (zero Python extension dependency), and
- `NanobindTransport`: direct in-process calls into `xgalois._xgalois` with
  zero IPC overhead.
"""

from __future__ import annotations

import collections
import json
import os
import pathlib
import selectors
import subprocess
import threading
import time
from typing import Any, Deque, Dict, List, Optional, Protocol, Union


class OracleRpcError(RuntimeError):
    """Raised when an oracle RPC call fails or times out."""

    def __init__(
        self,
        message: str,
        *,
        method: Optional[str] = None,
        request_id: Optional[int] = None,
        stderr_tail: str = "",
    ) -> None:
        self.method = method
        self.request_id = request_id
        self.stderr_tail = stderr_tail
        parts = [message]
        if method is not None:
            parts.append(f"method={method}")
        if request_id is not None:
            parts.append(f"id={request_id}")
        if stderr_tail:
            parts.append(f"stderr_tail={stderr_tail!r}")
        super().__init__(" | ".join(parts))


class OracleTransport(Protocol):
    """Protocol implemented by both subprocess and in-process transports."""

    def send_request(
        self,
        method: str,
        params: Optional[Dict[str, Any]] = None,
        *,
        timeout: Optional[float] = None,
    ) -> Dict[str, Any]: ...

    def stderr_tail(self, max_lines: int = 20) -> str: ...

    def close(self) -> None: ...


class SubprocessNdjsonTransport:
    """Manages a persistent xgalois_oracle_cli subprocess speaking NDJSON."""

    def __init__(
        self,
        cli_path: Union[str, pathlib.Path],
        *,
        timeout: float = 30.0,
        env: Optional[Dict[str, str]] = None,
    ) -> None:
        self.cli_path = pathlib.Path(cli_path)
        if not self.cli_path.is_file():
            raise FileNotFoundError(
                f"xgalois_oracle_cli binary not found: {self.cli_path}"
            )
        self.timeout = float(timeout)
        self._env = env
        self._proc: Optional[subprocess.Popen[bytes]] = None
        self._stderr_lines: Deque[str] = collections.deque(maxlen=200)
        self._stderr_thread: Optional[threading.Thread] = None
        self._stdout_buf = bytearray()
        self._next_id = 0
        self._start()

    def _start(self) -> None:
        self._proc = subprocess.Popen(
            [str(self.cli_path)],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            bufsize=0,
            env=self._env,
        )
        self._stderr_thread = threading.Thread(
            target=self._drain_stderr, daemon=True
        )
        self._stderr_thread.start()

    def _drain_stderr(self) -> None:
        assert self._proc is not None and self._proc.stderr is not None
        try:
            while True:
                raw = self._proc.stderr.readline()
                if not raw:
                    break
                self._stderr_lines.append(
                    raw.decode("utf-8", errors="replace").rstrip("\r\n")
                )
        except Exception:
            pass

    def stderr_tail(self, max_lines: int = 20) -> str:
        if not self._stderr_lines:
            return ""
        return "\n".join(list(self._stderr_lines)[-max_lines:])

    def _read_line_with_timeout(self, timeout: float) -> bytes:
        assert self._proc is not None and self._proc.stdout is not None
        fd = self._proc.stdout.fileno()
        deadline = time.monotonic() + timeout

        nl_idx = self._stdout_buf.find(b"\n")
        if nl_idx != -1:
            line = bytes(self._stdout_buf[: nl_idx + 1])
            del self._stdout_buf[: nl_idx + 1]
            return line

        sel = selectors.DefaultSelector()
        try:
            sel.register(fd, selectors.EVENT_READ)
            while True:
                remaining = deadline - time.monotonic()
                if remaining <= 0.0:
                    raise TimeoutError(
                        f"Timed out after {timeout:.2f}s waiting for oracle CLI response"
                    )
                events = sel.select(timeout=remaining)
                if not events:
                    raise TimeoutError(
                        f"Timed out after {timeout:.2f}s waiting for oracle CLI response"
                    )
                chunk = os.read(fd, 65536)
                if not chunk:
                    if self._stdout_buf:
                        line = bytes(self._stdout_buf)
                        self._stdout_buf.clear()
                        return line
                    rc = self._proc.poll()
                    raise EOFError(
                        f"xgalois_oracle_cli closed stdout unexpectedly (returncode={rc})"
                    )
                self._stdout_buf.extend(chunk)
                nl_idx = self._stdout_buf.find(b"\n")
                if nl_idx != -1:
                    line = bytes(self._stdout_buf[: nl_idx + 1])
                    del self._stdout_buf[: nl_idx + 1]
                    return line
        finally:
            sel.close()

    def send_request(
        self,
        method: str,
        params: Optional[Dict[str, Any]] = None,
        *,
        timeout: Optional[float] = None,
    ) -> Dict[str, Any]:
        if self._proc is None or self._proc.stdin is None:
            raise OracleRpcError("Transport is closed", method=method)

        self._next_id += 1
        req_id = self._next_id
        payload = {
            "id": req_id,
            "method": method,
            "params": params if params is not None else {},
        }
        wire = (json.dumps(payload, separators=(",", ":")) + "\n").encode("utf-8")

        try:
            self._proc.stdin.write(wire)
            self._proc.stdin.flush()
        except BrokenPipeError as exc:
            rc = self._proc.poll()
            raise OracleRpcError(
                f"Broken pipe sending request (returncode={rc}): {exc}",
                method=method,
                request_id=req_id,
                stderr_tail=self.stderr_tail(),
            ) from exc

        eff_timeout = self.timeout if timeout is None else float(timeout)
        try:
            raw_line = self._read_line_with_timeout(eff_timeout)
        except (TimeoutError, EOFError) as exc:
            raise OracleRpcError(
                str(exc),
                method=method,
                request_id=req_id,
                stderr_tail=self.stderr_tail(),
            ) from exc

        try:
            resp = json.loads(raw_line.decode("utf-8"))
        except json.JSONDecodeError as exc:
            raise OracleRpcError(
                f"Invalid JSON from oracle CLI: {exc} (raw={raw_line[:200]!r})",
                method=method,
                request_id=req_id,
                stderr_tail=self.stderr_tail(),
            ) from exc

        if not isinstance(resp, dict):
            raise OracleRpcError(
                f"Expected JSON object response, got {type(resp).__name__}",
                method=method,
                request_id=req_id,
                stderr_tail=self.stderr_tail(),
            )
        if resp.get("id") != req_id:
            raise OracleRpcError(
                f"Response id mismatch: expected {req_id}, got {resp.get('id')}",
                method=method,
                request_id=req_id,
                stderr_tail=self.stderr_tail(),
            )
        return resp

    def close(self) -> None:
        proc = self._proc
        if proc is None:
            return
        self._proc = None
        try:
            if proc.stdin is not None:
                proc.stdin.close()
        except Exception:
            pass
        try:
            proc.wait(timeout=2.0)
        except subprocess.TimeoutExpired:
            proc.terminate()
            try:
                proc.wait(timeout=2.0)
            except subprocess.TimeoutExpired:
                proc.kill()
                proc.wait(timeout=2.0)
        finally:
            for stream in (proc.stdout, proc.stderr):
                if stream is not None:
                    try:
                        stream.close()
                    except Exception:
                        pass
            if self._stderr_thread is not None:
                self._stderr_thread.join(timeout=1.0)
                self._stderr_thread = None

    def __enter__(self) -> "SubprocessNdjsonTransport":
        return self

    def __exit__(self, exc_type: Any, exc: Any, tb: Any) -> None:
        self.close()


class NanobindTransport:
    """In-process oracle transport invoking `xgalois._xgalois` directly."""

    def __init__(self, module: Optional[Any] = None) -> None:
        if module is None:
            from . import _xgalois as native_mod

            if native_mod is None:
                raise ImportError(
                    "NanobindTransport requires the compiled 'xgalois._xgalois' "
                    "extension. Build with `pip install .` or pass a CLI path to "
                    "SubprocessNdjsonTransport."
                )
            self._mod = native_mod
        else:
            self._mod = module
        self._next_id = 0
        self._closed = False

    def stderr_tail(self, max_lines: int = 20) -> str:
        return ""

    def close(self) -> None:
        self._closed = True

    def __enter__(self) -> "NanobindTransport":
        return self

    def __exit__(self, exc_type: Any, exc: Any, tb: Any) -> None:
        self.close()

    def _build_field(self, field_spec: Dict[str, Any]) -> tuple[Any, bool]:
        """Return (field_obj, is_binary)."""
        p = int(field_spec["p"])
        m = int(field_spec.get("m", 1))
        if m == 1:
            return self._mod.GFPrime(p), False
        if p == 2 and 2 <= m <= 16:
            if "modulus" not in field_spec:
                raise ValueError(
                    "Binary extension field requires explicit ascending 'modulus' array"
                )
            modulus = field_spec["modulus"]
            return self._mod.GF2X(m, modulus, True), True
        raise ValueError(
            "Unsupported field specification: expected prime field (m == 1) or "
            "binary extension field (p == 2, 2 <= m <= 16)"
        )

    def _handle_ping(self) -> Dict[str, Any]:
        return {
            "version": "1.0.0",
            "status": "ok",
            "transport": "nanobind",
        }

    def _handle_lrc_encode_and_repair(self, params: Dict[str, Any]) -> Dict[str, Any]:
        gf, is_binary = self._build_field(params["field"])
        n = int(params["n"])
        k = int(params["k"])
        r = int(params["r"])
        delta = int(params.get("delta", 2))
        partition_type = params.get("partition_type", None)

        if is_binary:
            code = self._mod.create_tamo_barg_lrc_binary(
                gf, n, k, r, delta, partition_type
            )
        else:
            code = self._mod.create_tamo_barg_lrc_prime(
                gf, n, k, r, delta, partition_type
            )

        msg = [int(x) for x in params["message"]]
        codeword = list(code.encode(msg))
        unencoded = list(code.unencode(codeword))

        erased_positions = [int(i) for i in params.get("erased_positions", [])]
        allow_global_fallback = bool(params.get("allow_global_fallback", True))

        corrupted = list(codeword)
        for pos in erased_positions:
            if 0 <= pos < n:
                corrupted[pos] = int(gf.add(corrupted[pos], 1))

        repaired = list(corrupted)
        repair_success = True
        symbols_read = 0
        num_local_repairs = 0
        num_global_repairs = 0
        used_global_fallback = False
        repair_error = ""

        if len(erased_positions) == 1:
            pos = erased_positions[0]
            res = code.local_repair(corrupted, pos, 0)
            repair_success = bool(res["success"])
            symbols_read = int(res["symbols_read"])
            repair_error = str(res.get("error_message", ""))
            if repair_success:
                repaired[pos] = int(res["recovered_symbol"])
                num_local_repairs = 1
        elif len(erased_positions) > 1:
            res = code.local_repair_multiple(
                corrupted, erased_positions, allow_global_fallback
            )
            repair_success = bool(res["success"])
            symbols_read = int(res["total_symbols_read"])
            num_local_repairs = int(res["num_local_repairs"])
            num_global_repairs = int(res["num_global_repairs"])
            used_global_fallback = bool(res["used_global_fallback"])
            repair_error = str(res.get("error_message", ""))
            if repair_success:
                repaired = list(res["codeword"])

        repair_sets = [list(code.repair_sets(i)[0]) for i in range(n)]
        out: Dict[str, Any] = {
            "codeword": codeword,
            "repaired_codeword": repaired,
            "minimum_distance": int(code.minimum_distance),
            "designed_distance": int(code.designed_distance),
            "symbols_read": symbols_read,
            "evaluation_points": list(code.evaluation_points()),
            "repair_sets": repair_sets,
            "fiber_partition": [list(f) for f in code.fiber_partition(0)],
            "good_poly": list(code.good_poly()),
            "repair_success": repair_success,
            "used_global_fallback": used_global_fallback,
            "num_local_repairs": num_local_repairs,
            "num_global_repairs": num_global_repairs,
            "unencoded_message": unencoded,
        }
        if repair_error:
            out["repair_error"] = repair_error
        return out

    def _handle_irs_decode(self, params: Dict[str, Any]) -> Dict[str, Any]:
        gf, _ = self._build_field(params["field"])
        eval_pts = [int(x) for x in params["evaluation_points"]]
        col_mults = params.get("column_multipliers", None)
        if col_mults is not None:
            col_mults = [int(x) for x in col_mults]

        if "dims" in params:
            dims = [int(d) for d in params["dims"]]
            irs = self._mod.InterleavedRS(gf, eval_pts, dims, col_mults)
        else:
            k = int(params["k"])
            ell = int(params["ell"])
            irs = self._mod.InterleavedRS(gf, eval_pts, k, ell, col_mults)

        received = [[int(v) for v in row] for row in params["received_matrix"]]
        mode = str(params.get("mode", "ssb"))
        fallback = bool(params.get("fallback_to_independent", True))
        erasure_positions = [int(i) for i in params.get("erasure_positions", [])]
        return dict(irs.decode(received, mode, fallback, erasure_positions))

    def _handle_grs_power_decode(self, params: Dict[str, Any]) -> Dict[str, Any]:
        gf, _ = self._build_field(params["field"])
        eval_pts = [int(x) for x in params["evaluation_points"]]
        col_mults = params.get("column_multipliers", None)
        if col_mults is not None:
            col_mults = [int(x) for x in col_mults]
        k = int(params["k"])
        grs = self._mod.GRSCode(gf, eval_pts, k, col_mults)

        received = [int(x) for x in params["received"]]
        powering_degree = params.get("powering_degree", None)
        if powering_degree is not None:
            powering_degree = int(powering_degree)
        mode = str(params.get("mode", "ssb"))
        fallback_to_bm = bool(params.get("fallback_to_bm", True))
        erasure_positions = [int(i) for i in params.get("erasure_positions", [])]
        return dict(
            grs.decode_power(
                received, powering_degree, mode, fallback_to_bm, erasure_positions
            )
        )

    def _handle_grs_guruswami_sudan(self, params: Dict[str, Any]) -> Dict[str, Any]:
        gf, _ = self._build_field(params["field"])
        eval_pts = [int(x) for x in params["evaluation_points"]]
        col_mults = params.get("column_multipliers", None)
        if col_mults is not None:
            col_mults = [int(x) for x in col_mults]
        k = int(params["k"])
        grs = self._mod.GRSCode(gf, eval_pts, k, col_mults)

        received = [int(x) for x in params["received"]]
        tau = int(params["tau"])
        return dict(grs.decode_guruswami_sudan(received, tau))

    def _handle_gabidulin_decode(self, params: Dict[str, Any]) -> Dict[str, Any]:
        gf, is_binary = self._build_field(params["field"])
        if not is_binary:
            raise ValueError(
                "gabidulin_decode requires a binary extension field GF(2^m)"
            )
        n = int(params["n"])
        k = int(params["k"])
        code = self._mod.GabidulinBinaryCode(gf, n, k)
        received = [int(x) for x in params["received"]]
        decoder = str(params.get("decoder", "GabidulinSkewEEA"))
        return dict(code.decode(received, decoder))

    def _handle_groups_and_rings_check(self, params: Dict[str, Any]) -> Dict[str, Any]:
        sub_op = str(params["sub_op"])
        if sub_op == "symmetric_group":
            deg = int(params["degree"])
            order = str(self._mod.symmetric_group_order(deg))
            return {"order": order, "bsgs_order": order}
        if sub_op == "alternating_group":
            deg = int(params["degree"])
            order = str(self._mod.alternating_group_order(deg))
            return {"order": order, "bsgs_order": order}
        if sub_op == "dihedral_group":
            n = int(params["n"])
            order = str(self._mod.dihedral_group_order(n))
            return {"order": order}
        if sub_op == "negacyclic_mul":
            q = int(params["modulus_q"])
            deg = int(params["degree"])
            a = [int(x) for x in params["a"]]
            b = [int(x) for x in params["b"]]
            return {"product": list(self._mod.negacyclic_mul(q, deg, a, b))}
        if sub_op == "crt_int":
            residues = [int(x) for x in params["residues"]]
            moduli = [int(x) for x in params["moduli"]]
            ans = self._mod.crt_int(residues, moduli)
            if ans is None:
                return {"consistent": False}
            return {
                "consistent": True,
                "solution": str(ans[0]),
                "lcm": str(ans[1]),
            }
        raise ValueError(f"Unknown groups_and_rings_check sub_op: {sub_op!r}")

    def _dispatch(self, method: str, params: Dict[str, Any]) -> Dict[str, Any]:
        if method == "ping":
            return self._handle_ping()
        if method == "lrc_encode_and_repair":
            return self._handle_lrc_encode_and_repair(params)
        if method == "irs_decode":
            return self._handle_irs_decode(params)
        if method == "grs_power_decode":
            return self._handle_grs_power_decode(params)
        if method == "grs_guruswami_sudan":
            return self._handle_grs_guruswami_sudan(params)
        if method == "gabidulin_decode":
            return self._handle_gabidulin_decode(params)
        if method == "groups_and_rings_check":
            return self._handle_groups_and_rings_check(params)
        raise ValueError(f"Unknown RPC method: {method!r}")

    def send_request(
        self,
        method: str,
        params: Optional[Dict[str, Any]] = None,
        *,
        timeout: Optional[float] = None,
    ) -> Dict[str, Any]:
        del timeout
        if self._closed:
            raise OracleRpcError("Transport is closed", method=method)
        self._next_id += 1
        req_id = self._next_id
        try:
            result = self._dispatch(method, params if params is not None else {})
            return {"id": req_id, "ok": True, "result": result}
        except Exception as exc:
            return {"id": req_id, "ok": False, "error": str(exc)}


class OracleClient:
    """High-level RPC client supporting both subprocess and nanobind transports."""

    def __init__(
        self,
        transport_or_cli: Union[
            SubprocessNdjsonTransport, NanobindTransport, str, pathlib.Path
        ],
        *,
        timeout: float = 30.0,
    ) -> None:
        if isinstance(
            transport_or_cli, (SubprocessNdjsonTransport, NanobindTransport)
        ):
            self._transport: Union[SubprocessNdjsonTransport, NanobindTransport] = (
                transport_or_cli
            )
            self._owns_transport = False
        else:
            self._transport = SubprocessNdjsonTransport(
                transport_or_cli, timeout=timeout
            )
            self._owns_transport = True

    @classmethod
    def from_nanobind(cls, module: Optional[Any] = None) -> "OracleClient":
        """Create an OracleClient backed by in-process nanobind bindings."""
        client = cls(NanobindTransport(module=module))
        client._owns_transport = True
        return client

    @classmethod
    def from_cli(
        cls,
        cli_path: Union[str, pathlib.Path],
        *,
        timeout: float = 30.0,
        env: Optional[Dict[str, str]] = None,
    ) -> "OracleClient":
        """Create an OracleClient backed by `xgalois_oracle_cli` subprocess."""
        transport = SubprocessNdjsonTransport(cli_path, timeout=timeout, env=env)
        client = cls(transport)
        client._owns_transport = True
        return client

    @property
    def transport(self) -> Union[SubprocessNdjsonTransport, NanobindTransport]:
        return self._transport

    def call(
        self,
        method: str,
        params: Optional[Dict[str, Any]] = None,
        *,
        timeout: Optional[float] = None,
    ) -> Dict[str, Any]:
        resp = self._transport.send_request(method, params, timeout=timeout)
        if not resp.get("ok", False):
            raise OracleRpcError(
                f"Oracle returned error: {resp.get('error', 'unknown error')}",
                method=method,
                request_id=resp.get("id"),
                stderr_tail=self._transport.stderr_tail(),
            )
        result = resp.get("result")
        if not isinstance(result, dict):
            raise OracleRpcError(
                f"Expected dict result from oracle, got {type(result).__name__}",
                method=method,
                request_id=resp.get("id"),
                stderr_tail=self._transport.stderr_tail(),
            )
        return result

    def ping(self) -> Dict[str, Any]:
        return self.call("ping", {})

    def close(self) -> None:
        if self._owns_transport:
            self._transport.close()

    def __enter__(self) -> "OracleClient":
        return self

    def __exit__(self, exc_type: Any, exc: Any, tb: Any) -> None:
        self.close()


__all__ = [
    "OracleRpcError",
    "OracleTransport",
    "SubprocessNdjsonTransport",
    "NanobindTransport",
    "OracleClient",
]
