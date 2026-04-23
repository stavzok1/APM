"""
Console tee + post-run kernel hints for long pipeline runs.

``dmesg`` is best-effort (Linux only; may require elevated privileges for full buffer).
"""

from __future__ import annotations

import atexit
import os
import signal
import subprocess
import sys
from pathlib import Path
from typing import List, Optional, TextIO


class _TeeTextIO:
    """Write to multiple text streams (typically TTY + log file)."""

    __slots__ = ("_streams",)

    def __init__(self, *streams: TextIO) -> None:
        self._streams = streams

    def write(self, data: str) -> int:
        n = 0
        for s in self._streams:
            n = s.write(data)
        return n

    def flush(self) -> None:
        for s in self._streams:
            s.flush()

    def isatty(self) -> bool:
        return self._streams[0].isatty()

    def fileno(self) -> int:
        return self._streams[0].fileno()


_LOG_STATE: dict = {
    "path": None,
    "restore": None,
    "finalized": False,
}


def _append_dmesg_kill_oom_section(log_path: Path) -> None:
    """Append recent kernel lines mentioning kill / OOM (if any)."""
    banner = "\n" + "=" * 72 + "\n"
    banner += "Kernel ring buffer — lines containing kill / oom (best-effort)\n"
    banner += "=" * 72 + "\n"

    if sys.platform != "linux":
        with open(log_path, "a", encoding="utf-8") as fh:
            fh.write(banner)
            fh.write("(skipped: not Linux — no dmesg)\n")
        return

    try:
        proc = subprocess.run(
            ["dmesg", "-T"],
            capture_output=True,
            text=True,
            timeout=12,
            check=False,
        )
        text = proc.stdout or ""
        err = (proc.stderr or "").strip()
    except FileNotFoundError:
        with open(log_path, "a", encoding="utf-8") as fh:
            fh.write(banner)
            fh.write("`dmesg` not found on PATH.\n")
        return
    except subprocess.TimeoutExpired:
        with open(log_path, "a", encoding="utf-8") as fh:
            fh.write(banner)
            fh.write("`dmesg` timed out.\n")
        return
    except OSError as e:
        with open(log_path, "a", encoding="utf-8") as fh:
            fh.write(banner)
            fh.write(f"dmesg failed: {e}\n")
        return

    keys = ("kill", "killed", "oom", "out of memory")
    hits: List[str] = [
        ln for ln in text.splitlines() if any(k in ln.lower() for k in keys)
    ]

    with open(log_path, "a", encoding="utf-8") as fh:
        fh.write(banner)
        if proc.returncode != 0:
            fh.write(f"dmesg exit code {proc.returncode}\n")
            if err:
                fh.write(f"dmesg stderr: {err[:2000]}\n")
        if not hits:
            fh.write(
                "(no matching lines in current ring buffer — try `sudo dmesg -T | "
                "grep -iE 'kill|oom'` if the process was OOM-killed earlier)\n"
            )
        else:
            tail = hits[-120:]
            fh.write("\n".join(tail))
            fh.write("\n")


def finalize_run_log() -> None:
    """Restore stdio and append dmesg section once."""
    if _LOG_STATE["finalized"]:
        return
    _LOG_STATE["finalized"] = True
    restore = _LOG_STATE["restore"]
    if restore is not None:
        restore()
        _LOG_STATE["restore"] = None
    path = _LOG_STATE["path"]
    if path is not None:
        _append_dmesg_kill_oom_section(Path(path))
        _LOG_STATE["path"] = None


def _signal_finalize(signum: int, frame: Optional[object]) -> None:
    finalize_run_log()
    signal.signal(signum, signal.SIG_DFL)
    os.kill(os.getpid(), signum)


def install_run_tee(log_path: Path) -> None:
    """
    Tee stdout/stderr to *log_path* and register exit/signal hooks.

    ``finalize_run_log()`` runs via atexit and on SIGINT/SIGTERM (not SIGKILL).
    """
    if _LOG_STATE.get("restore") is not None:
        finalize_run_log()
    _LOG_STATE["finalized"] = False

    log_path = Path(log_path)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    fh = open(log_path, "a", encoding="utf-8", buffering=1)

    old_out: TextIO = sys.__stdout__
    old_err: TextIO = sys.__stderr__

    def _restore() -> None:
        try:
            sys.stdout.flush()
            sys.stderr.flush()
        except Exception:
            pass
        sys.stdout = old_out
        sys.stderr = old_err
        try:
            fh.flush()
            fh.close()
        except Exception:
            pass

    sys.stdout = _TeeTextIO(old_out, fh)
    sys.stderr = _TeeTextIO(old_err, fh)

    _LOG_STATE["path"] = str(log_path)
    _LOG_STATE["restore"] = _restore
    atexit.register(finalize_run_log)
    signal.signal(signal.SIGINT, _signal_finalize)
    signal.signal(signal.SIGTERM, _signal_finalize)
