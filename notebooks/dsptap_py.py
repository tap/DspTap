"""ctypes bridge to the DspTap C ABI, shared by the verification notebooks.

Loads build_capi/libdsptap_capi.{so,dylib,dll} relative to the repo root,
building it first if missing (requires cmake in PATH):

    cmake -B build_capi -S tools/capi
    cmake --build build_capi

The C ABI (tools/capi/) wraps the *same* portable DSP headers the consuming
libraries compile — so the notebooks exercise the real shipping code, not a
Python re-implementation. Exposed primitives: the real FFT (`RealFFT`, in the
double and float profiles), the YIN pitch detector (`Yin`), the TD-PSOLA
shifter (`Psola`), and the peak-locked phase-vocoder shifter (`Pvoc`, with
optional LPC formant preservation), the log-mel/PCEN feature front end
(`LogMel`) and the fixed-ratio decimators to 16 kHz (`Decimator`).

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
        "dsptap_fft_create":        ([ctypes.c_int, ctypes.c_int], vp),
        "dsptap_fft_destroy":       ([vp], None),
        "dsptap_fft_size":          ([vp], ctypes.c_int),
        "dsptap_fft_num_bins":      ([vp], ctypes.c_int),
        "dsptap_fft_profile":       ([vp], ctypes.c_int),
        "dsptap_fft_sample_bytes":  ([vp], ctypes.c_int),
        "dsptap_fft_backend":       ([], ctypes.c_char_p),
        "dsptap_fft_forward":       ([vp, f64p, f64p], ctypes.c_int),
        "dsptap_fft_inverse":       ([vp, f64p, f64p], ctypes.c_int),
        "dsptap_fft_forward_inplace_raw": ([vp, vp], ctypes.c_int),
        "dsptap_fft_inverse_inplace_raw": ([vp, vp], ctypes.c_int),
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
        try:
            fn = getattr(lib, name)
        except AttributeError as e:
            raise RuntimeError(
                f"{_lib_path()} lacks {name}: build_capi/ is older than tools/capi/dsptap_capi.h. "
                "Delete build_capi/ and import again (it rebuilds on first import).") from e
        fn.argtypes = argtypes
        fn.restype = restype
    return lib


_lib = load()


def _f64(x: np.ndarray):
    return np.ascontiguousarray(x, dtype=np.float64)


class RealFFT:
    """tap::dsp::basic_real_fft<Sample> — `profile` "double" (the golden model) or "float"
    (the embedded profile: Ooura, or the vDSP/CMSIS backend the build selected — see
    `RealFFT.backend()`). "q15" and "q31" are reserved for the fixed-point profiles
    (Stage 3c of docs/audit-fft-and-code-smells.md) and raise until they land.

    `forward`/`inverse` take and return float64 arrays in the header's PACKED layout and
    Ooura's exp(+i) sign convention, whatever the profile (the float profile converts at the
    boundary). `forward_inplace_raw`/`inverse_inplace_raw` run on the profile's native dtype
    with no conversion at all — the path that measures the embedded profile's own arithmetic —
    and the raw inverse is UNSCALED, exactly as in fft.h.

    Packing (fft.h): packed[0] = DC, packed[1] = Nyquist (both real);
    packed[2k] + 1j*packed[2k+1] = bin k for 1 <= k < n/2, with W = exp(+2*pi*i/n).
    `unpack`/`pack` convert to and from a complex array of n/2 + 1 bins in the ENGINEERING
    convention (exp(-2*pi*i/n), what numpy.fft.rfft returns) — i.e. they CONJUGATE, so
    `unpack(forward(x))` matches `numpy.fft.rfft(x)` bin for bin.
    """

    PROFILES = {"double": 0, "float": 1, "q15": 2, "q31": 3}
    _DTYPES = {0: np.float64, 1: np.float32, 2: np.int16, 3: np.int32}

    def __init__(self, size: int, profile: str | int = "double"):
        code = self.PROFILES.get(profile) if isinstance(profile, str) else profile
        if code not in self.PROFILES.values():
            raise ValueError(f"profile must be one of {sorted(self.PROFILES)} "
                             f"(or the codes {sorted(self.PROFILES.values())}), not {profile!r}")
        self._h = _lib.dsptap_fft_create(size, code)
        if not self._h:
            if code in (2, 3):
                raise NotImplementedError("the Q15/Q31 FFT profiles land at Stage 3c")
            raise ValueError(f"size must be a power of two >= 4, not {size!r}")
        self.size = size
        self.profile = next(k for k, v in self.PROFILES.items() if v == code)
        self.dtype = np.dtype(self._DTYPES[code])
        assert self.dtype.itemsize == _lib.dsptap_fft_sample_bytes(self._h)

    def __del__(self):
        if getattr(self, "_h", None):
            _lib.dsptap_fft_destroy(self._h)

    @property
    def num_bins(self) -> int:
        return _lib.dsptap_fft_num_bins(self._h)

    @staticmethod
    def backend() -> str:
        """The float32 engine this build compiled: "ooura", "accelerate" or "cmsis"."""
        return _lib.dsptap_fft_backend().decode()

    def forward(self, x: np.ndarray) -> np.ndarray:
        """Packed spectrum (float64, length n) of n real samples; the float profile rounds x to
        float32 at the boundary."""
        x = _f64(x)
        self._check_len(x)
        out = np.empty_like(x)
        _lib.dsptap_fft_forward(self._h, x.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                                out.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))
        return out

    def inverse(self, packed: np.ndarray) -> np.ndarray:
        """Normalized inverse (scaled by 2/n): inverse(forward(x)) reproduces x."""
        packed = _f64(packed)
        self._check_len(packed)
        out = np.empty_like(packed)
        _lib.dsptap_fft_inverse(self._h, packed.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                                out.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))
        return out

    def forward_inplace_raw(self, data: np.ndarray) -> np.ndarray:
        """In-place forward on the profile's native dtype (no conversion); returns `data`."""
        self._check_raw(data)
        _lib.dsptap_fft_forward_inplace_raw(self._h, data.ctypes.data_as(ctypes.c_void_p))
        return data

    def inverse_inplace_raw(self, data: np.ndarray) -> np.ndarray:
        """In-place UNSCALED inverse on the native dtype (multiply by 2/n yourself); returns
        `data`."""
        self._check_raw(data)
        _lib.dsptap_fft_inverse_inplace_raw(self._h, data.ctypes.data_as(ctypes.c_void_p))
        return data

    def _check_len(self, a: np.ndarray) -> None:
        if a.ndim != 1 or a.size != self.size:
            raise ValueError(f"expected a 1-D array of length {self.size}, got shape {a.shape}")

    def _check_raw(self, data: np.ndarray) -> None:
        if (not isinstance(data, np.ndarray) or data.dtype != self.dtype or data.ndim != 1
                or data.size != self.size or not data.flags.c_contiguous or not data.flags.writeable):
            raise TypeError(f"raw buffers must be 1-D, contiguous, writeable {self.dtype} of length {self.size}")

    @staticmethod
    def unpack(packed: np.ndarray) -> np.ndarray:
        """Packed Ooura spectrum -> complex bins 0..n/2 in the engineering convention
        (conjugated: rfft-compatible)."""
        packed = np.asarray(packed, dtype=np.float64)
        n = packed.size
        bins = np.empty(n // 2 + 1, dtype=np.complex128)
        bins[0] = packed[0]
        bins[-1] = packed[1]
        bins[1:-1] = packed[2::2] - 1j * packed[3::2]
        return bins

    @staticmethod
    def pack(bins: np.ndarray) -> np.ndarray:
        """Engineering-convention complex bins 0..n/2 -> packed Ooura spectrum (float64,
        length n). The imaginary parts of DC and Nyquist are discarded (the packing has no
        slot for them; a real signal's are zero)."""
        bins = np.asarray(bins, dtype=np.complex128)
        n = 2 * (bins.size - 1)
        packed = np.empty(n, dtype=np.float64)
        packed[0] = bins[0].real
        packed[1] = bins[-1].real
        packed[2::2] = bins[1:-1].real
        packed[3::2] = -bins[1:-1].imag
        return packed

    def spectrum(self, x: np.ndarray) -> np.ndarray:
        """Convenience: unpack(forward(x)) — rfft-compatible complex bins."""
        return self.unpack(self.forward(x))


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

