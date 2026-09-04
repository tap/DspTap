"""ctypes bridge to the DspTap C ABI, shared by the verification notebooks.

Loads build_capi/libdsptap_capi.{so,dylib,dll} relative to the repo root,
building it first if missing (requires cmake in PATH):

    cmake -B build_capi -S tools/capi
    cmake --build build_capi

The C ABI (tools/capi/) wraps the *same* portable DSP headers the consuming
libraries compile — so the notebooks exercise the real shipping code, not a
Python re-implementation. Exposed primitives: the YIN pitch detector (`Yin`),
the TD-PSOLA shifter (`Psola`), and the peak-locked phase-vocoder shifter
(`Pvoc`, with optional LPC formant preservation), the log-mel/PCEN feature
front end (`LogMel`) and the fixed-ratio decimators to 16 kHz (`Decimator`).

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
        "dsptap_log_mel_create":    ([ctypes.c_double, ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int,
                                      ctypes.c_double, ctypes.c_double, ctypes.c_int, ctypes.c_double], vp),
        "dsptap_log_mel_destroy":   ([vp], None),
        "dsptap_log_mel_set_log":   ([vp, ctypes.c_double, ctypes.c_double, ctypes.c_double], ctypes.c_int),
        "dsptap_log_mel_set_pcen":  ([vp, ctypes.c_int, ctypes.c_double, ctypes.c_double, ctypes.c_double,
                                      ctypes.c_double, ctypes.c_double], ctypes.c_int),
        "dsptap_log_mel_reset":     ([vp], ctypes.c_int),
        "dsptap_log_mel_bands":     ([vp], ctypes.c_int),
        "dsptap_log_mel_latency":   ([vp], ctypes.c_int),
        "dsptap_log_mel_contract_version": ([], ctypes.c_int),
        "dsptap_log_mel_process":   ([vp, f64p, ctypes.c_int, f64p, ctypes.c_int], ctypes.c_int),
        "dsptap_decimator_create":  ([ctypes.c_int, ctypes.c_int], vp),
        "dsptap_decimator_destroy": ([vp], None),
        "dsptap_decimator_taps":    ([vp], ctypes.c_int),
        "dsptap_decimator_latency": ([vp], ctypes.c_int),
        "dsptap_decimator_reset":   ([vp], ctypes.c_int),
        "dsptap_decimator_outputs_for": ([vp, ctypes.c_int], ctypes.c_int),
        "dsptap_decimator_process": ([vp, f64p, ctypes.c_int, f64p, ctypes.c_int], ctypes.c_int),
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


class LogMel:
    """tap::dsp::log_mel — the double-precision golden profile of the wake-word front end.

    Defaults are the reference geometry (16 kHz, frame 400, hop 160, FFT 512, 40 HTK mel
    bands 20-7600 Hz, periodic Hann). `pcen=True` switches the output to per-channel energy
    normalization with the paper's defaults; pass a dict to override them.
    """

    def __init__(self, sample_rate: float = 16000.0, frame: int = 400, hop: int = 160, fft_size: int = 512,
                 bands: int = 40, fmin_hz: float = 20.0, fmax_hz: float = 7600.0, sqrt_window: bool = False,
                 preemphasis: float = 0.0, log: tuple[float, float, float] | None = None,
                 pcen: bool | dict | None = None):
        self._h = _lib.dsptap_log_mel_create(sample_rate, frame, hop, fft_size, bands, fmin_hz, fmax_hz,
                                             int(sqrt_window), preemphasis)
        if not self._h:
            raise ValueError("bad log_mel geometry")
        if log is not None:
            if _lib.dsptap_log_mel_set_log(self._h, *log) != 0:
                raise ValueError("bad log constants")
        if pcen:
            p = {"smoother": 0.025, "alpha": 0.98, "delta": 2.0, "power": 0.5, "epsilon": 1e-6}
            if isinstance(pcen, dict):
                p.update(pcen)
            if _lib.dsptap_log_mel_set_pcen(self._h, 1, p["smoother"], p["alpha"], p["delta"], p["power"],
                                            p["epsilon"]) != 0:
                raise ValueError("bad PCEN parameters")
        self.hop = hop

    def __del__(self):
        if getattr(self, "_h", None):
            _lib.dsptap_log_mel_destroy(self._h)

    @property
    def bands(self) -> int:
        return _lib.dsptap_log_mel_bands(self._h)

    @property
    def latency(self) -> int:
        return _lib.dsptap_log_mel_latency(self._h)

    @staticmethod
    def contract_version() -> int:
        return _lib.dsptap_log_mel_contract_version()

    def reset(self) -> None:
        _lib.dsptap_log_mel_reset(self._h)

    def process(self, x: np.ndarray) -> np.ndarray:
        """Features for every completed hop in x, shape (frames, bands); state carries over."""
        x = _f64(x)
        max_frames = x.size // self.hop + 1
        out = np.zeros((max_frames, self.bands))
        n = _lib.dsptap_log_mel_process(self._h, x.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), x.size,
                                        out.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), max_frames)
        return out[:max(n, 0)]


class Decimator:
    """tap::dsp::basic_decimator<float, M> — 32/48/96 kHz to 16 kHz (the float golden model)."""

    def __init__(self, ratio: int, transparent: bool = False):
        self._h = _lib.dsptap_decimator_create(ratio, int(transparent))
        if not self._h:
            raise ValueError("ratio must be 2, 3 or 6")
        self.ratio = ratio

    def __del__(self):
        if getattr(self, "_h", None):
            _lib.dsptap_decimator_destroy(self._h)

    @property
    def taps(self) -> int:
        return _lib.dsptap_decimator_taps(self._h)

    @property
    def latency(self) -> int:
        """Group delay in input samples."""
        return _lib.dsptap_decimator_latency(self._h)

    def reset(self) -> None:
        _lib.dsptap_decimator_reset(self._h)

    def process(self, x: np.ndarray) -> np.ndarray:
        x = _f64(x)
        max_out = _lib.dsptap_decimator_outputs_for(self._h, x.size)
        out = np.zeros(max(max_out, 0))
        n = _lib.dsptap_decimator_process(self._h, x.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), x.size,
                                          out.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), out.size)
        return out[:max(n, 0)]

