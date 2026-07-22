"""ctypes bridge to the DspTap C ABI, shared by the verification notebooks.

Loads build_capi/libdsptap_capi.{so,dylib,dll} relative to the repo root,
building it first if missing (requires cmake in PATH):

    cmake -B build_capi -S tools/capi
    cmake --build build_capi

The C ABI (tools/capi/) wraps the *same* portable DSP headers the consuming
libraries compile — so the notebooks exercise the real shipping code, not a
Python re-implementation. Exposed primitives: the YIN pitch detector (`Yin`),
the TD-PSOLA shifter (`Psola`), and the peak-locked phase-vocoder shifter
(`Pvoc`, with optional LPC formant preservation).

Copyright 2026 Timothy Place and the DspTap contributors. MIT License.
"""

from __future__ import annotations

import ctypes
import pathlib
import subprocess
import sys

import numpy as np

# The repo root (this file lives in notebooks/).
ROOT = pathlib.Path(__file__).resolve().parent.parent

# Categorical palette for the notebooks (colorblind-safe, fixed assignment
# order — never cycled). Sequential maps use viridis; diverging use RdBu_r.
PALETTE = ["#4269d0", "#efb118", "#ff725c", "#6cc5b0", "#3ca951", "#ff8ab7", "#a463f2"]

_BUILD = ROOT / "build_capi"


def _lib_path() -> pathlib.Path:
    stem = "dsptap_capi"
    names = {"linux": f"lib{stem}.so", "darwin": f"lib{stem}.dylib", "win32": f"{stem}.dll"}
    name = next(v for k, v in names.items() if sys.platform.startswith(k))
    for cand in (_BUILD / name, _BUILD / "Release" / name, _BUILD / "Debug" / name):
        if cand.exists():
            return cand
    return _BUILD / name


def _build_lib() -> None:
    subprocess.run(["cmake", "-B", str(_BUILD), "-S", str(ROOT / "tools" / "capi")],
                   cwd=ROOT, check=True, capture_output=True)
    subprocess.run(["cmake", "--build", str(_BUILD), "--config", "Release", "--parallel"],
                   cwd=ROOT, check=True, capture_output=True)


def load() -> ctypes.CDLL:
    if not _lib_path().exists():
        print("building dsptap_capi ...")
        _build_lib()
    lib = ctypes.CDLL(str(_lib_path()))

    vp = ctypes.c_void_p
    f64p = ctypes.POINTER(ctypes.c_double)
    sigs = {
        "dsptap_yin_create":        ([ctypes.c_int, ctypes.c_int, ctypes.c_int], vp),
        "dsptap_yin_destroy":       ([vp], None),
        "dsptap_yin_set_threshold": ([vp, ctypes.c_double], ctypes.c_int),
        "dsptap_yin_frame_size":    ([vp], ctypes.c_int),
        "dsptap_yin_analyze":       ([vp, f64p, f64p, f64p], ctypes.c_int),
        "dsptap_yin_track":         ([vp, f64p, ctypes.c_int, ctypes.c_int, f64p, ctypes.c_int], ctypes.c_int),
        "dsptap_psola_create":      ([ctypes.c_int], vp),
        "dsptap_psola_destroy":     ([vp], None),
        "dsptap_psola_latency":     ([vp], ctypes.c_int),
        "dsptap_psola_clear":       ([vp], ctypes.c_int),
        "dsptap_psola_process":     ([vp, f64p, f64p, ctypes.c_int, ctypes.c_double, ctypes.c_double],
                                     ctypes.c_int),
        "dsptap_pvoc_create":       ([ctypes.c_int], vp),
        "dsptap_pvoc_destroy":      ([vp], None),
        "dsptap_pvoc_latency":      ([vp], ctypes.c_int),
        "dsptap_pvoc_set_formant":  ([vp, ctypes.c_int], ctypes.c_int),
        "dsptap_pvoc_clear":        ([vp], ctypes.c_int),
        "dsptap_pvoc_process":      ([vp, f64p, f64p, ctypes.c_int, ctypes.c_double], ctypes.c_int),
    }
    for name, (argtypes, restype) in sigs.items():
        fn = getattr(lib, name)
        fn.argtypes = argtypes
        fn.restype = restype
    return lib


_lib = load()


def _f64(x: np.ndarray):
    return np.ascontiguousarray(x, dtype=np.float64)


class Yin:
    """tap::dsp::yin — the double-precision golden profile."""

    def __init__(self, window: int, tau_min: int, tau_max: int, threshold: float | None = None):
        self._h = _lib.dsptap_yin_create(window, tau_min, tau_max)
        if not self._h:
            raise ValueError("bad yin geometry")
        if threshold is not None:
            _lib.dsptap_yin_set_threshold(self._h, threshold)

    def __del__(self):
        if getattr(self, "_h", None):
            _lib.dsptap_yin_destroy(self._h)

    @property
    def frame_size(self) -> int:
        return _lib.dsptap_yin_frame_size(self._h)

    def analyze(self, frame: np.ndarray) -> tuple[float, float]:
        frame = _f64(frame)
        assert frame.size == self.frame_size
        period = ctypes.c_double()
        aper = ctypes.c_double()
        _lib.dsptap_yin_analyze(self._h, frame.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                                ctypes.byref(period), ctypes.byref(aper))
        return period.value, aper.value

    def track(self, x: np.ndarray, hop: int) -> np.ndarray:
        """Periods (samples; 0 = unvoiced) every `hop` samples across x."""
        x = _f64(x)
        out = np.zeros(x.size // hop + 1)
        n = _lib.dsptap_yin_track(self._h, x.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), x.size,
                                  hop, out.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), out.size)
        return out[:max(n, 0)]


class Psola:
    """tap::dsp::psola — TD-PSOLA shifter (caller supplies the period)."""

    def __init__(self, max_period: int):
        self._h = _lib.dsptap_psola_create(max_period)
        if not self._h:
            raise ValueError("bad max_period")

    def __del__(self):
        if getattr(self, "_h", None):
            _lib.dsptap_psola_destroy(self._h)

    @property
    def latency(self) -> int:
        return _lib.dsptap_psola_latency(self._h)

    def process(self, x: np.ndarray, period: float, ratio: float) -> np.ndarray:
        x = _f64(x)
        out = np.zeros_like(x)
        _lib.dsptap_psola_process(self._h, x.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                                  out.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), x.size, period, ratio)
        return out


class Pvoc:
    """tap::dsp::pvoc — peak-locked phase-vocoder shifter."""

    def __init__(self, fft_size: int = 1024, formant: bool = False):
        self._h = _lib.dsptap_pvoc_create(fft_size)
        if not self._h:
            raise ValueError("fft_size must be a power of two >= 64")
        _lib.dsptap_pvoc_set_formant(self._h, int(formant))

    def __del__(self):
        if getattr(self, "_h", None):
            _lib.dsptap_pvoc_destroy(self._h)

    @property
    def latency(self) -> int:
        return _lib.dsptap_pvoc_latency(self._h)

    def process(self, x: np.ndarray, ratio: float) -> np.ndarray:
        x = _f64(x)
        out = np.zeros_like(x)
        _lib.dsptap_pvoc_process(self._h, x.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                                 out.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), x.size, ratio)
        return out
