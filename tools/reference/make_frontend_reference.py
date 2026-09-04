#!/usr/bin/env python3
"""Generates tests/reference/frontend_vectors.h — the independent golden leg
for the two wake-word front-end primitives, log_mel.h and decimate.h.

This is the numpy restatement of both headers' formula-level contracts, and
it is deliberately the ONLY numpy restatement in the family: MuTap's KWS
feature module imports it through the DspTap submodule rather than carrying
a second copy, so there is one source of truth for the feature *values* and
one header for the *formulas* (the ownership rule in MuTap's wake-word plan,
section 5).

log_mel: the reference geometry (16 kHz, frame 400, hop 160, FFT 512, 40 HTK
mel bands 20-7600 Hz, periodic Hann, no pre-emphasis) over a deterministic
test signal that excites every band — six tones plus xorshift noise, with a
+20 dB level step at the midpoint so PCEN's gain tracking is exercised. The
expected features are emitted for the plain-log path and the PCEN path.

decimate: for each ratio (2, 3, 6) and profile, the minimal odd tap count
whose Kaiser design meets the profile's stopband spec with >= 1 dB margin on
a 25 Hz grid is searched here and emitted as a constant the header's table
must match; the expected output is numpy.convolve over the same design on
deterministic xorshift noise, cast to float32 the way the float engine
stores it.

Run from the repo root:  python3 tools/reference/make_frontend_reference.py
Re-run only when a contract changes; commit the regenerated header.
"""
from __future__ import annotations

import pathlib

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[2]
OUT = ROOT / "tests" / "reference" / "frontend_vectors.h"

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

SR = 16000.0
FRAME = 400
HOP = 160
FFT = 512
BANDS = 40
FMIN = 20.0
FMAX = 7600.0
LOG_FLOOR = 1e-10
LOG_SHIFT = 5.0
LOG_SCALE = 5.0
PCEN_S = 0.025
PCEN_ALPHA = 0.98
PCEN_DELTA = 2.0
PCEN_R = 0.5
PCEN_EPS = 1e-6

TONES = [(150.0, 0.30), (440.0, 0.25), (1000.0, 0.20), (2500.0, 0.15), (4000.0, 0.10), (6500.0, 0.05)]
NOISE_AMP = 0.02
N_SAMPLES = 8000  # 0.5 s -> 50 frames
STEP_GAIN = 10.0  # +20 dB from the midpoint


def hz_to_mel(f: np.ndarray | float) -> np.ndarray | float:
    return 2595.0 * np.log10(1.0 + np.asarray(f) / 700.0)


def mel_to_hz(m: np.ndarray | float) -> np.ndarray | float:
    return 700.0 * (10.0 ** (np.asarray(m) / 2595.0) - 1.0)


def mel_weights(sr: float, fft: int, bands: int, fmin: float, fmax: float) -> np.ndarray:
    """Unit-peak triangles on the HTK mel scale; shape (bands, fft/2 + 1)."""
    edges = mel_to_hz(np.linspace(hz_to_mel(fmin), hz_to_mel(fmax), bands + 2))
    bins = np.arange(fft // 2 + 1) * sr / fft
    w = np.zeros((bands, bins.size))
    for b in range(bands):
        lo, mid, hi = edges[b], edges[b + 1], edges[b + 2]
        rise = (bins - lo) / (mid - lo)
        fall = (hi - bins) / (hi - mid)
        w[b] = np.maximum(0.0, np.minimum(rise, fall))
    return w


def test_signal() -> np.ndarray:
    n = np.arange(N_SAMPLES, dtype=np.float64)
    x = np.zeros(N_SAMPLES)
    for f, a in TONES:
        x += a * np.sin(2.0 * np.pi * f * n / SR)
    x += NOISE_AMP * xorshift32(N_SAMPLES, 0x2545F491)
    x[N_SAMPLES // 2:] *= STEP_GAIN
    return x


def mel_energies(x: np.ndarray, preemph: float = 0.0) -> np.ndarray:
    """Streaming-aligned mel band powers: frame t ends at sample (t+1)*hop."""
    y = x.copy()
    if preemph != 0.0:
        y[1:] = x[1:] - preemph * x[:-1]
    n_frames = y.size // HOP
    win = 0.5 - 0.5 * np.cos(2.0 * np.pi * np.arange(FRAME) / FRAME)  # periodic Hann
    w = mel_weights(SR, FFT, BANDS, FMIN, FMAX)
    padded = np.concatenate([np.zeros(FRAME), y])
    e = np.zeros((n_frames, BANDS))
    for t in range(n_frames):
        end = FRAME + (t + 1) * HOP
        seg = padded[end - FRAME:end] * win
        spec = np.fft.rfft(seg, n=FFT)  # zero-padded at the end of the frame
        power = spec.real ** 2 + spec.imag ** 2
        e[t] = w @ power
    return e


def log_features(e: np.ndarray) -> np.ndarray:
    return (np.log10(e + LOG_FLOOR) + LOG_SHIFT) / LOG_SCALE


def pcen_features(e: np.ndarray) -> np.ndarray:
    out = np.zeros_like(e)
    m = e[0].copy()  # initial smoother state: the first frame's energy
    for t in range(e.shape[0]):
        m = (1.0 - PCEN_S) * m + PCEN_S * e[t]
        out[t] = (e[t] / (PCEN_EPS + m) ** PCEN_ALPHA + PCEN_DELTA) ** PCEN_R - PCEN_DELTA ** PCEN_R
    return out


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


def main() -> None:
    x = test_signal()
    e = mel_energies(x)
    log_ref = log_features(e)
    pcen_ref = pcen_features(e)

    parts = [
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
        f"    inline constexpr std::size_t k_mel_frames = {log_ref.shape[0]};",
        f"    inline constexpr std::size_t k_mel_bands  = {BANDS};",
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
    parts.append("} // namespace frontend_ref")
    parts.append("// NOLINTEND(readability-identifier-naming)")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join(parts) + "\n")
    print(f"wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
