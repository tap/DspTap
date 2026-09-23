"""ctypes bridge to the DspTap C ABI, shared by the verification notebooks.

Loads build_capi/libdsptap_capi.{so,dylib,dll} relative to the repo root,
building it first if missing (requires cmake in PATH):

    cmake -B build_capi -S tools/capi
    cmake --build build_capi

The C ABI (tools/capi/) wraps the *same* portable DSP headers the consuming
libraries compile — so the notebooks exercise the real shipping code, not a
Python re-implementation. Exposed primitives: the real FFT (`RealFFT`, in the
double, float, Q15 and Q31 profiles, the fixed-point ones under both scaling
policies), the YIN pitch detector (`Yin`), the TD-PSOLA
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

    # Every handle is an opaque pointer to a library-private struct (dsptap_capi.h: typedef struct
    # dsptap_fft_s* dsptap_fft, ...): one pointer at the ABI, so c_void_p carries each of them.
    vp = ctypes.c_void_p
    f64p = ctypes.POINTER(ctypes.c_double)
    i32p = ctypes.POINTER(ctypes.c_int)
    sigs = {
        "dsptap_fft_create":        ([ctypes.c_int, ctypes.c_int], vp),
        "dsptap_fft_destroy":       ([vp], None),
        "dsptap_fft_size":          ([vp], ctypes.c_int),
        "dsptap_fft_num_bins":      ([vp], ctypes.c_int),
        "dsptap_fft_profile":       ([vp], ctypes.c_int),
        "dsptap_fft_sample_bytes":  ([vp], ctypes.c_int),
        "dsptap_fft_backend":       ([], ctypes.c_char_p),
        "dsptap_fft_fixed_scaling_exponent": ([vp], ctypes.c_int),
        "dsptap_fft_forward":       ([vp, f64p, f64p], ctypes.c_int),
        "dsptap_fft_inverse":       ([vp, f64p, f64p], ctypes.c_int),
        "dsptap_fft_forward_exp":   ([vp, f64p, f64p, i32p], ctypes.c_int),
        "dsptap_fft_inverse_exp":   ([vp, f64p, f64p, i32p], ctypes.c_int),
        "dsptap_fft_forward_inplace_raw": ([vp, vp], ctypes.c_int),
        "dsptap_fft_inverse_inplace_raw": ([vp, vp], ctypes.c_int),
        "dsptap_fft_forward_inplace_raw_exp": ([vp, vp, i32p], ctypes.c_int),
        "dsptap_fft_inverse_inplace_raw_exp": ([vp, vp, i32p], ctypes.c_int),
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
    """tap::dsp::basic_real_fft<Sample, Scaling> through the C ABI. `profile` is one of

        "double"   the golden model
        "float"    the embedded floating profile (the split-radix engine, or the vDSP/CMSIS
                   backend the build selected — see `RealFFT.backend()`)
        "q15"      std::int16_t, Q0.15 I/O, tap::dsp::scaling::fixed
        "q31"      std::int32_t, Q0.31 I/O, scaling::fixed
        "q15_bfp"  Q15 under scaling::block_floating
        "q31_bfp"  Q31 under scaling::block_floating

    THE EXPONENT (fft.h's fixed-point contract, restated in the bridge's terms). A fixed-point
    transform leaves its buffer scaled by 2^-e and returns e: read the native buffer as
    fractions of full scale (x / 2^15 or x / 2^31) and let G be the "double" profile on the same
    input read the same way — after a forward, G's forward == buffer * 2^e; after an inverse,
    G's UNNORMALIZED inverse == buffer * 2^e (same packing, same exp(+i) sign, up to the kernel's
    rounding noise). Under fixed scaling e is the constant `exponent_bound` (log2 n for Q15,
    log2 n + 1 for Q31); under block floating point 0 <= e <= `exponent_bound`, chosen per block
    from the data's headroom. The fixed-point inverse applies NO 2/n; a raw round trip
    reconstructs x == out * 2^(e_fwd + e_inv + 1 - log2 n). The floating profiles report e == 0
    everywhere an exponent appears.

    `forward`/`inverse` take and return float64 arrays in the header's PACKED layout and
    Ooura's exp(+i) sign convention, whatever the profile, in the DOUBLE profile's units: the
    float profile converts at the boundary; the fixed-point profiles quantize the input to
    their Q format at the C boundary (round half away from zero, saturating at +-full scale)
    and return the native result already multiplied by 2^e, so the array compares directly to
    the double profile's on the same quantized input. The exponent of the last `forward`/
    `inverse` is `last_exponent`; `forward_with_exponent`/`inverse_with_exponent` return
    `(array, e)`. Note that the fixed-point `inverse` is the UNNORMALIZED inverse (no 2/n,
    unlike the floating profiles' `inverse`), and that its input is quantized to [-1, 1), so
    the floating-style `inverse(forward(x))` round trip is a raw-path matter for fixed point.

    `forward_inplace_raw`/`inverse_inplace_raw` run on the profile's native dtype (`dtype`:
    float64, float32, int16 or int32) with no conversion at all — the path that measures the
    embedded profile's own arithmetic — and RETURN e (0 for the floating profiles); the raw
    inverse is UNSCALED, exactly as in fft.h. `to_native` quantizes a float64 signal the way the
    boundary does, for building raw buffers; `from_native` reads one back as fractions times 2^e.

    Packing (fft.h): packed[0] = DC, packed[1] = Nyquist (both real);
    packed[2k] + 1j*packed[2k+1] = bin k for 1 <= k < n/2, with W = exp(+2*pi*i/n).
    `unpack`/`pack` convert to and from a complex array of n/2 + 1 bins in the ENGINEERING
    convention (exp(-2*pi*i/n), what numpy.fft.rfft returns) — i.e. they CONJUGATE, so
    `unpack(forward(x))` matches `numpy.fft.rfft(x)` bin for bin.
    """

    PROFILES = {"double": 0, "float": 1, "q15": 2, "q31": 3, "q15_bfp": 4, "q31_bfp": 5}
    _DTYPES = {0: np.float64, 1: np.float32, 2: np.int16, 3: np.int32, 4: np.int16, 5: np.int32}
    _FRAC_BITS = {0: 0, 1: 0, 2: 15, 3: 31, 4: 15, 5: 31}

    def __init__(self, size: int, profile: str | int = "double"):
        code = self.PROFILES.get(profile) if isinstance(profile, str) else profile
        if code not in self.PROFILES.values():
            raise ValueError(f"profile must be one of {sorted(self.PROFILES)} "
                             f"(or the codes {sorted(self.PROFILES.values())}), not {profile!r}")
        self._h = _lib.dsptap_fft_create(size, code)
        if not self._h:
            # The C ABI returns NULL for a size the profile does not support: the fixed-point
            # profiles are bounded at 65536 (fft/fixed_point.h, k_max_size; the capi reads the
            # constant from the header and refuses above it, so the header's assert is never
            # reached from Python).
            bound = " in [4, 65536] (the fixed-point profiles' bound)" if code >= 2 else " >= 4"
            raise ValueError(f"size must be a power of two{bound}, not {size!r}")
        self.size = size
        self.profile = next(k for k, v in self.PROFILES.items() if v == code)
        self.dtype = np.dtype(self._DTYPES[code])
        self.frac_bits = self._FRAC_BITS[code]
        self.is_fixed_point = code >= 2
        self.is_block_floating = code >= 4
        self.last_exponent = 0
        assert self.dtype.itemsize == _lib.dsptap_fft_sample_bytes(self._h)

    def __del__(self):
        if getattr(self, "_h", None):
            _lib.dsptap_fft_destroy(self._h)

    @property
    def num_bins(self) -> int:
        return _lib.dsptap_fft_num_bins(self._h)

    @property
    def exponent_bound(self) -> int:
        """basic_real_fft<Sample, Scaling>::fixed_scaling_exponent(n): the exponent every
        transform of a fixed-scaling profile returns and the upper bound under block floating
        point; 0 for the floating profiles."""
        return _lib.dsptap_fft_fixed_scaling_exponent(self._h)

    @staticmethod
    def backend() -> str:
        """The float32 engine this build compiled: "split_radix" (the C++20 engine of
        fft/split_radix.h, bit-identical to the Ooura C it replaced; "ooura" until Stage 2c),
        "accelerate" or "cmsis". The double profile is always the split-radix engine; the
        fixed-point profiles are always the portable int32 kernel."""
        return _lib.dsptap_fft_backend().decode()

    def forward_with_exponent(self, x: np.ndarray) -> tuple[np.ndarray, int]:
        """Packed spectrum (float64, length n) of n real samples, and the transform's exponent e.
        The array is already scaled by 2^e (it is the double profile's spectrum of the quantized
        input, up to the kernel's noise); the float profile rounds x to float32 at the boundary,
        the fixed-point profiles quantize it to their Q format."""
        x = _f64(x)
        self._check_len(x)
        out = np.empty_like(x)
        e = ctypes.c_int(0)
        _lib.dsptap_fft_forward_exp(self._h, x.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                                    out.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), ctypes.byref(e))
        self.last_exponent = e.value
        return out, e.value

    def forward(self, x: np.ndarray) -> np.ndarray:
        """forward_with_exponent(x)[0]; the exponent is left in `last_exponent`."""
        return self.forward_with_exponent(x)[0]

    def inverse_with_exponent(self, packed: np.ndarray) -> tuple[np.ndarray, int]:
        """Inverse at the double boundary, and its exponent. Floating profiles: normalized by
        2/n, so inverse(forward(x)) reproduces x. Fixed-point profiles: the packed spectrum is
        quantized to the Q format (fractions of full scale, saturating outside [-1, 1)) and the
        result is the UNNORMALIZED inverse times 2^e — the double profile's raw inverse of the
        same quantized spectrum, with no 2/n."""
        packed = _f64(packed)
        self._check_len(packed)
        out = np.empty_like(packed)
        e = ctypes.c_int(0)
        _lib.dsptap_fft_inverse_exp(self._h, packed.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                                    out.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), ctypes.byref(e))
        self.last_exponent = e.value
        return out, e.value

    def inverse(self, packed: np.ndarray) -> np.ndarray:
        """inverse_with_exponent(packed)[0]; the exponent is left in `last_exponent`."""
        return self.inverse_with_exponent(packed)[0]

    def forward_inplace_raw(self, data: np.ndarray) -> int:
        """In-place forward on the profile's native dtype (no conversion); returns the exponent e
        (0 for the floating profiles). `data` is transformed in place."""
        self._check_raw(data)
        e = ctypes.c_int(0)
        if _lib.dsptap_fft_forward_inplace_raw_exp(self._h, data.ctypes.data_as(ctypes.c_void_p),
                                                   ctypes.byref(e)) != 0:
            raise RuntimeError("dsptap_fft_forward_inplace_raw_exp failed")
        self.last_exponent = e.value
        return e.value

    def inverse_inplace_raw(self, data: np.ndarray) -> int:
        """In-place UNSCALED inverse on the native dtype; returns the exponent e (0 for the
        floating profiles, where you multiply by 2/n yourself; for fixed point the round trip
        is x == out * 2^(e_fwd + e_inv + 1 - log2 n))."""
        self._check_raw(data)
        e = ctypes.c_int(0)
        if _lib.dsptap_fft_inverse_inplace_raw_exp(self._h, data.ctypes.data_as(ctypes.c_void_p),
                                                   ctypes.byref(e)) != 0:
            raise RuntimeError("dsptap_fft_inverse_inplace_raw_exp failed")
        self.last_exponent = e.value
        return e.value

    def to_native(self, x: np.ndarray) -> np.ndarray:
        """A float64 signal in fractions of full scale as a native raw buffer: the plain cast for
        the floating profiles; for the fixed-point profiles the boundary's quantization (round
        half away from zero, saturating at the rails: tap::dsp::detail::round_sat), so
        `forward_inplace_raw(to_native(x))` sees exactly the samples `forward(x)` does."""
        x = _f64(x)
        if not self.is_fixed_point:
            return x.astype(self.dtype)
        scaled = x * float(2 ** self.frac_bits)
        r = np.where(scaled < 0.0, scaled - 0.5, scaled + 0.5)
        info = np.iinfo(self.dtype)
        r = np.clip(np.trunc(r), info.min, info.max)
        return r.astype(self.dtype)

    def from_native(self, data: np.ndarray, exponent: int = 0) -> np.ndarray:
        """A native buffer read as fractions of full scale times 2^exponent (float64), the
        reading the exponent contract is written in; for the floating profiles the widening
        cast times 2^exponent."""
        return np.asarray(data, dtype=np.float64) * float(2.0 ** (exponent - self.frac_bits))

    def _check_len(self, a: np.ndarray) -> None:
        if a.ndim != 1 or a.size != self.size:
            raise ValueError(f"expected a 1-D array of length {self.size}, got shape {a.shape}")

    def _check_raw(self, data: np.ndarray) -> None:
        # The ABI's raw buffers must be aligned for the native type (dsptap_capi.h): an int16
        # view at an odd byte offset is contiguous and writeable but not aligned, so the flag
        # is checked too.
        if (not isinstance(data, np.ndarray) or data.dtype != self.dtype or data.ndim != 1
                or data.size != self.size or not data.flags.c_contiguous or not data.flags.writeable
                or not data.flags.aligned):
            raise TypeError(f"raw buffers must be 1-D, contiguous, aligned, writeable {self.dtype} "
                            f"of length {self.size}")

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

