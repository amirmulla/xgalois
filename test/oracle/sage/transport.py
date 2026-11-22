#!/usr/bin/env python3
"""NDJSON subprocess transport and RPC client for xgalois_oracle_cli.

Pure Python 3 (non-preparsed): works identically under both `sage -python`
and standard `python3`.
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
from typing import Any, Deque, Dict, Optional, Union


class OracleRpcError(RuntimeError):
  """Raised when an RPC call to xgalois_oracle_cli fails or times out."""

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


class OracleClient:
  """High-level RPC client wrapping SubprocessNdjsonTransport."""

  def __init__(
      self,
      transport_or_cli: Union[SubprocessNdjsonTransport, str, pathlib.Path],
      *,
      timeout: float = 30.0,
  ) -> None:
    if isinstance(transport_or_cli, SubprocessNdjsonTransport):
      self._transport = transport_or_cli
      self._owns_transport = False
    else:
      self._transport = SubprocessNdjsonTransport(
          transport_or_cli, timeout=timeout
      )
      self._owns_transport = True

  @property
  def transport(self) -> SubprocessNdjsonTransport:
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
          f"Oracle CLI returned error: {resp.get('error', 'unknown error')}",
          method=method,
          request_id=resp.get("id"),
          stderr_tail=self._transport.stderr_tail(),
      )
    result = resp.get("result")
    if not isinstance(result, dict):
      raise OracleRpcError(
          f"Expected dict result from oracle CLI, got {type(result).__name__}",
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
