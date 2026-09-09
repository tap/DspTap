#!/usr/bin/env python3
"""Generates tests/reference/frontend_vectors.h — the independent golden leg
for the two wake-word front-end primitives, log_mel.h and decimate.h.

This is the numpy restatement of both headers' formula-level contracts, and
it is deliberately the ONLY numpy restatement in the family: MuTap's KWS
feature module imports it through the DspTap submodule rather than carrying
a second copy, so there is one source of truth for the feature *values* and
one header for the *formulas* (the ownership rule in MuTap's wake-word plan,
section 5).

log_mel: `Geometry` mirrors tap::dsp::log_mel_geometry field for field, so a
consumer can hand one set of values to this reference and to the C++ engine
(through the C ABI bridge) and compare; `features(x, g)` is what log_mel.h
emits for x at geometry g. Two geometries are emitted over the same
deterministic test signal — six tones plus xorshift noise, exciting every
band, with a +20 dB level step at the midpoint so PCEN's gain tracking is
exercised. `REFERENCE` (16 kHz, frame 400, hop 160, FFT 512, 40 HTK mel bands
20-7600 Hz, periodic Hann, no pre-emphasis) goes to frontend_vectors.h, the
M1 record; `TUNED`, every runtime field off its default, goes to
frontend_vectors_tuned.h together with its geometry as constants, so the
C++ <-> numpy parity is pinned away from the defaults as well as at them.
Plain-log and PCEN paths are emitted for both.

MuTap's tools/ml/kws/kws_features.py imports this module for its parity
self-check and its band-support check only; training features come from the
shipping C++ front end through the bridge (decided 8 September 2026).

Reproducibility: the mel vectors are float64 numpy output, and another numpy
build moves them at the last decimal place — measured 2026-09, numpy 2.5.3
on macOS arm64 against the committed header: <= 1.4e-15 on the log path,
<= 1.2e-14 on the PCEN path, the decimator vectors bit-identical — which the
C++ pins at 1e-13 absorb. Commit a regenerated header only when a contract
changes; a refactor of this script is checked by regenerating before and
after it in one environment and diffing.

decimate: for each ratio (2, 3, 6) and profile, the minimal odd tap count
whose Kaiser design meets the profile's stopband spec with >= 1 dB margin on
a 25 Hz grid is searched here and emitted as a constant the header's table
must match; the expected output is numpy.convolve over the same design on
deterministic xorshift noise, cast to float32 the way the float engine
stores it.

Run from the repo root:  python3 tools/reference/make_frontend_reference.py
(--which reference | tuned | all, default all). Re-run only when a contract
changes; commit the regenerated headers.
"""
from __future__ import annotations

import argparse
import dataclasses
import pathlib

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[2]
OUT = ROOT / "tests" / "reference" / "frontend_vectors.h"
OUT_TUNED = ROOT / "tests" / "reference" / "frontend_vectors_tuned.h"

# ---------------------------------------------------------------------------
# Shared deterministic noise (xorshift32), identical to the C++ tests.
# ---------------------------------------------------------------------------


def xorshift32(count: int, seed: int) -> np.ndarray:
    """Uniform values in [-1, 1) from xorshift32, matching tests' C++ copy."""
    s = seed & 0xFFFFFFFF
    out = np.empty(count, np.float64)
    for i in range(count):
        s ^= (s << 13) & 0xFFFFFFFF
        s ^= s >> 17
        s ^= (s << 5) & 0xFFFFFFFF
        out[i] = (s % 65536 - 32768) / 32768.0
    return out


# ---------------------------------------------------------------------------
# log_mel — formula-level contract (see include/tap/dsp/log_mel.h)
# ---------------------------------------------------------------------------


@dataclasses.dataclass(frozen=True)
class Pcen:
    """tap::dsp::pcen_params, field for field."""

    enabled: bool = False
    smoother: float = 0.025
    alpha: float = 0.98
    delta: float = 2.0
    power: float = 0.5
    epsilon: float = 1e-6


@dataclasses.dataclass(frozen=True)
class Geometry:
    """tap::dsp::log_mel_geometry, field for field; window is "hann" or "sqrt_hann"."""

    sample_rate: float = 16000.0
    frame: int = 400
    hop: int = 160
    fft_size: int = 512
    bands: int = 40
    fmin_hz: float = 20.0
    fmax_hz: float = 7600.0
    window: str = "hann"
    preemphasis: float = 0.0
    log_floor: float = 1e-10
    log_shift: float = 5.0
    log_scale: float = 5.0
    pcen: Pcen = Pcen()

    def valid(self) -> bool:
        """log_mel_geometry::valid(), restated."""
        pow2 = self.fft_size >= 4 and (self.fft_size & (self.fft_size - 1)) == 0
        p = self.pcen
        return (self.sample_rate > 0.0 and self.hop >= 1 and self.frame >= self.hop and pow2
                and self.fft_size >= self.frame and self.bands >= 1 and self.fmin_hz >= 0.0
                and self.fmax_hz > self.fmin_hz and self.fmax_hz <= 0.5 * self.sample_rate
                and self.log_floor > 0.0 and self.log_scale > 0.0 and self.window in ("hann", "sqrt_hann")
                and 0.0 < p.smoother <= 1.0 and 0.0 <= p.alpha <= 1.0 and p.delta > 0.0
                and 0.0 < p.power <= 1.0 and p.epsilon > 0.0)


REFERENCE = Geometry()

# Every runtime field off its default — frame/hop/FFT ratio, band count and
# edges, sqrt-Hann, pre-emphasis, the log affine, every PCEN parameter — so
# the parity pin at this geometry says something the default one does not.
TUNED = Geometry(frame=320, hop=80, fft_size=512, bands=32, fmin_hz=50.0, fmax_hz=7000.0,
                 window="sqrt_hann", preemphasis=0.97, log_floor=1e-8, log_shift=4.0, log_scale=4.0,
                 pcen=Pcen(smoother=0.04, alpha=0.9, delta=1.0, power=0.4, epsilon=1e-5))

TONES = [(150.0, 0.30), (440.0, 0.25), (1000.0, 0.20), (2500.0, 0.15), (4000.0, 0.10), (6500.0, 0.05)]
NOISE_AMP = 0.02
N_SAMPLES = 8000  # 0.5 s -> 50 frames at the reference hop, 100 at the tuned one
STEP_GAIN = 10.0  # +20 dB from the midpoint


def hz_to_mel(f: np.ndarray | float) -> np.ndarray | float:
    return 2595.0 * np.log10(1.0 + np.asarray(f) / 700.0)


def mel_to_hz(m: np.ndarray | float) -> np.ndarray | float:
    return 700.0 * (10.0 ** (np.asarray(m) / 2595.0) - 1.0)


def mel_weights(sr: float, fft: int, bands: int, fmin: float, fmax: float) -> np.ndarray:
    """Unit-peak triangles on the HTK mel scale; shape (bands, fft/2 + 1).

    A band with no FFT bin inside it comes out all-zero here, and log_mel.h
    zero-fills the same band silently: a consumer choosing a geometry asserts
    that every row has a nonzero entry before trusting the features.
    """
    edges = mel_to_hz(np.linspace(hz_to_mel(fmin), hz_to_mel(fmax), bands + 2))
    bins = np.arange(fft // 2 + 1) * sr / fft
    w = np.zeros((bands, bins.size))
    for b in range(bands):
        lo, mid, hi = edges[b], edges[b + 1], edges[b + 2]
        rise = (bins - lo) / (mid - lo)
        fall = (hi - bins) / (hi - mid)
        w[b] = np.maximum(0.0, np.minimum(rise, fall))
    return w


def window(g: Geometry = REFERENCE) -> np.ndarray:
    """Periodic Hann over the frame, or its square root (log_mel.h build_window)."""
    w = 0.5 - 0.5 * np.cos(2.0 * np.pi * np.arange(g.frame) / g.frame)
    return np.sqrt(w) if g.window == "sqrt_hann" else w


def test_signal() -> np.ndarray:
    n = np.arange(N_SAMPLES, dtype=np.float64)
    x = np.zeros(N_SAMPLES)
    for f, a in TONES:
        x += a * np.sin(2.0 * np.pi * f * n / REFERENCE.sample_rate)
    x += NOISE_AMP * xorshift32(N_SAMPLES, 0x2545F491)
    x[N_SAMPLES // 2:] *= STEP_GAIN
    return x


def mel_energies(x: np.ndarray, g: Geometry = REFERENCE) -> np.ndarray:
    """Streaming-aligned mel band powers: frame t ends at sample (t+1)*hop - 1."""
    if not g.valid():
        raise ValueError("invalid log_mel geometry")
    y = x.copy()
    if g.preemphasis != 0.0:
        y[1:] = x[1:] - g.preemphasis * x[:-1]  # x[-1] = 0, as on the stream
    n_frames = y.size // g.hop
    win = window(g)
    w = mel_weights(g.sample_rate, g.fft_size, g.bands, g.fmin_hz, g.fmax_hz)
    padded = np.concatenate([np.zeros(g.frame), y])
    e = np.zeros((n_frames, g.bands))
    for t in range(n_frames):
        end = g.frame + (t + 1) * g.hop
        seg = padded[end - g.frame:end] * win
        spec = np.fft.rfft(seg, n=g.fft_size)  # zero-padded at the end of the frame
        power = spec.real ** 2 + spec.imag ** 2
        e[t] = w @ power
    return e


def log_features(e: np.ndarray, g: Geometry = REFERENCE) -> np.ndarray:
    return (np.log10(e + g.log_floor) + g.log_shift) / g.log_scale


def pcen_features(e: np.ndarray, g: Geometry = REFERENCE) -> np.ndarray:
    p = g.pcen
    out = np.zeros_like(e)
    m = e[0].copy()  # initial smoother state: the first frame's energy
    for t in range(e.shape[0]):
        m = (1.0 - p.smoother) * m + p.smoother * e[t]
        out[t] = (e[t] / (p.epsilon + m) ** p.alpha + p.delta) ** p.power - p.delta ** p.power
    return out


def features(x: np.ndarray, g: Geometry = REFERENCE) -> np.ndarray:
    """What log_mel.h emits for x at g: the PCEN path if g.pcen.enabled, else plain log."""
    e = mel_energies(x, g)
    return pcen_features(e, g) if g.pcen.enabled else log_features(e, g)


# ---------------------------------------------------------------------------
# decimate — formula-level contract (see include/tap/dsp/decimate.h)
# ---------------------------------------------------------------------------

OUT_RATE = 16000.0
PROFILES = {  # name: (stopband_atten_db, passband_hz)
    "economy": (70.0, 7000.0),
    "transparent": (100.0, 7600.0),
}
RATIOS = (2, 3, 6)
N_DEC_INPUT = 3000


def bessel_i0(x: np.ndarray) -> np.ndarray:
    return np.i0(x)


def kaiser_beta(atten_db: float) -> float:
    if atten_db > 50.0:
        return 0.1102 * (atten_db - 8.7)
    if atten_db > 21.0:
        return 0.5842 * (atten_db - 21.0) ** 0.4 + 0.07886 * (atten_db - 21.0)
    return 0.0


def design(taps: int, cutoff_norm: float, beta: float) -> np.ndarray:
    """tap::dsp::design_prototype with num_phases = 1, sum(h) == 1."""
    i = np.arange(taps, dtype=np.float64)
    center = 0.5 * (taps - 1)
    t = i - center
    u = t / center
    w = bessel_i0(beta * np.sqrt(np.maximum(0.0, 1.0 - u * u))) / bessel_i0(beta)
    h = cutoff_norm * np.sinc(cutoff_norm * t) * w
    return h / h.sum()


def response_db(h: np.ndarray, fs: float, f: np.ndarray) -> np.ndarray:
    n = np.arange(h.size)
    z = np.exp(-2j * np.pi * np.outer(f, n) / fs)
    return 20.0 * np.log10(np.abs(z @ h) + 1e-300)


def find_taps(ratio: int, atten_db: float, passband_hz: float) -> int:
    fs_in = OUT_RATE * ratio
    stop_hz = OUT_RATE - passband_hz
    cutoff = 1.0 / ratio  # 2 * (fs_out/2) / fs_in: cutoff at the output Nyquist
    beta = kaiser_beta(atten_db)
    est = int(np.ceil((atten_db - 8.0) / (2.285 * 2.0 * np.pi * (stop_hz - passband_hz) / fs_in)))
    taps = est if est % 2 == 1 else est + 1
    f_stop = np.arange(stop_hz, fs_in / 2.0 + 1.0, 25.0)
    f_pass = np.arange(0.0, passband_hz + 1.0, 25.0)
    while True:
        h = design(taps, cutoff, beta)
        worst_stop = response_db(h, fs_in, f_stop).max()
        worst_pass = np.abs(response_db(h, fs_in, f_pass)).max()
        if worst_stop <= -(atten_db + 1.0) and worst_pass <= 0.1:
            return taps
        taps += 2


def decimate_reference(x: np.ndarray, h: np.ndarray, ratio: int) -> np.ndarray:
    """y[k] = sum_t h[t] x[k*M - t], zero history, one output per M inputs."""
    return np.convolve(x, h, mode="full")[: x.size][::ratio]


# ---------------------------------------------------------------------------
# Emit
# ---------------------------------------------------------------------------


def fmt_array(name: str, values: np.ndarray, ctype: str) -> str:
    flat = np.asarray(values).reshape(-1)
    if ctype == "float":
        items = [f"{v:.9e}f" for v in flat]
    else:
        items = [f"{v:.17e}" for v in flat]
    lines = []
    for i in range(0, len(items), 4):
        lines.append("        " + ", ".join(items[i:i + 4]) + ",")
    body = "\n".join(lines)
    return f"    inline constexpr std::array<{ctype}, {flat.size}> {name} = {{\n{body}\n    }};\n"


HEADER = [
    "// Generated by tools/reference/make_frontend_reference.py — DO NOT EDIT.",
    "// Independent golden reference (numpy float64) for log_mel.h and decimate.h;",
    "// see that script for the contracts it restates and the tolerance arguments.",
    "// SPDX-License-Identifier: MIT",
    "// Copyright 2026 Timothy Place and the DspTap contributors.",
    "// NOLINTBEGIN(readability-identifier-naming)",
    "#pragma once",
    "",
    "#include <array>",
    "#include <cstddef>",
    "",
    "namespace frontend_ref {",
    "",
]
FOOTER = ["} // namespace frontend_ref", "// NOLINTEND(readability-identifier-naming)"]


def write(path: pathlib.Path, parts: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(HEADER + parts + FOOTER) + "\n")
    print(f"wrote {path.relative_to(ROOT)}")


def emit_reference() -> None:
    x = test_signal()
    e = mel_energies(x)
    log_ref = log_features(e)
    pcen_ref = pcen_features(e)

    parts = [
        f"    inline constexpr std::size_t k_mel_frames = {log_ref.shape[0]};",
        f"    inline constexpr std::size_t k_mel_bands  = {REFERENCE.bands};",
        "",
        fmt_array("k_mel_log", log_ref, "double"),
        fmt_array("k_mel_pcen", pcen_ref, "double"),
    ]

    xin = xorshift32(N_DEC_INPUT, 0x9E3779B9).astype(np.float32) * np.float32(0.9)
    parts.append(fmt_array("k_dec_input", xin, "float"))
    for pname, (atten, passband) in PROFILES.items():
        for ratio in RATIOS:
            taps = find_taps(ratio, atten, passband)
            h = design(taps, 1.0 / ratio, kaiser_beta(atten))
            y = decimate_reference(xin.astype(np.float64), h, ratio)
            print(f"{pname:12s} ratio {ratio}: taps {taps}, outputs {y.size}")
            parts.append(f"    inline constexpr std::size_t k_dec_taps_{pname}_{ratio} = {taps};")
            parts.append(fmt_array(f"k_dec_{pname}_{ratio}", y.astype(np.float32), "float"))
    write(OUT, parts)


def emit_tuned() -> None:
    g = TUNED
    x = test_signal()
    e = mel_energies(x, g)
    log_ref = log_features(e, g)
    pcen_ref = pcen_features(e, g)
    sqrt_hann = "true" if g.window == "sqrt_hann" else "false"
    parts = [
        "    // TUNED geometry, every runtime field off its default; the C++ test",
        "    // builds its log_mel_geometry from these so the two sides cannot drift.",
        f"    inline constexpr double      k_tuned_sample_rate   = {g.sample_rate!r};",
        f"    inline constexpr std::size_t k_tuned_frame         = {g.frame};",
        f"    inline constexpr std::size_t k_tuned_hop           = {g.hop};",
        f"    inline constexpr std::size_t k_tuned_fft_size      = {g.fft_size};",
        f"    inline constexpr std::size_t k_tuned_bands         = {g.bands};",
        f"    inline constexpr double      k_tuned_fmin_hz       = {g.fmin_hz!r};",
        f"    inline constexpr double      k_tuned_fmax_hz       = {g.fmax_hz!r};",
        f"    inline constexpr bool        k_tuned_sqrt_hann     = {sqrt_hann};",
        f"    inline constexpr double      k_tuned_preemphasis   = {g.preemphasis!r};",
        f"    inline constexpr double      k_tuned_log_floor     = {g.log_floor!r};",
        f"    inline constexpr double      k_tuned_log_shift     = {g.log_shift!r};",
        f"    inline constexpr double      k_tuned_log_scale     = {g.log_scale!r};",
        f"    inline constexpr double      k_tuned_pcen_smoother = {g.pcen.smoother!r};",
        f"    inline constexpr double      k_tuned_pcen_alpha    = {g.pcen.alpha!r};",
        f"    inline constexpr double      k_tuned_pcen_delta    = {g.pcen.delta!r};",
        f"    inline constexpr double      k_tuned_pcen_power    = {g.pcen.power!r};",
        f"    inline constexpr double      k_tuned_pcen_epsilon  = {g.pcen.epsilon!r};",
        f"    inline constexpr std::size_t k_tuned_frames        = {log_ref.shape[0]};",
        "",
        fmt_array("k_tuned_log", log_ref, "double"),
        fmt_array("k_tuned_pcen", pcen_ref, "double"),
    ]
    write(OUT_TUNED, parts)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--which", choices=("all", "reference", "tuned"), default="all",
                    help="which header(s) to regenerate (default: all)")
    args = ap.parse_args()
    if args.which in ("all", "reference"):
        emit_reference()
    if args.which in ("all", "tuned"):
        emit_tuned()


if __name__ == "__main__":
    main()
