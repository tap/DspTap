# DspTap

[![CI](https://github.com/tap/DspTap/actions/workflows/ci.yml/badge.svg)](https://github.com/tap/DspTap/actions/workflows/ci.yml)
[![Bench](https://github.com/tap/DspTap/actions/workflows/bench.yml/badge.svg)](https://github.com/tap/DspTap/actions/workflows/bench.yml)
[![Tap House Style](https://github.com/tap/DspTap/actions/workflows/style.yml/badge.svg)](https://github.com/tap/DspTap/actions/workflows/style.yml)

Shared DSP primitives for the **Tap** family of audio libraries. Header-only,
plain portable C++ (C++20, standard library only), no Max/Min or framework
dependency — consumed as a git submodule by the individual libraries.

Today it holds ten primitives (the ninth, the [recursive substrate](#the-recursive-substrate--two-path-half-bands-and-the-guard-band-filter)
under the SampleRateTap family's IIR stages, since 0.3.0; the tenth, the
[Butterworth low-pass stage](#tapdspbutterworth--the-butterworth-low-pass-stage)
over it, since 0.6.0), plus the [FIR substrate](#the-fir-substrate) —
the shared design-math / sample-format / kernel layer under the SampleRateTap
family's two engines, `async` and `bridge` — and the [scalar helpers](#tapdspmathh--periodic-hann-and-decibel-helpers)
(`tap/dsp/math.h`) the consuming libraries call instead of re-typing them:

## `tap::dsp::real_fft` — real FFT with a fixed numeric contract

`include/tap/dsp/fft.h` is one real-FFT contract — the packing, the
`W = exp(+2πi/N)` sign convention, the unnormalized in-place inverse, the
numbers in the tables below — over four sample types, with the engine behind
each floating profile a **template parameter** since Stage 4 of the audit
(`basic_real_fft<Sample, Policy = detail::default_real_fft_policy_t<Sample>>`,
where the second argument is the engine for `float` / `double` — defaulting
to `default_real_fft_engine_t<Sample>`, chosen once by the sample type and
the build — and the scaling policy for Q15 / Q31; any engine can be named
explicitly beside the default in the same binary), and each engine stating
its size range and shareability as contract numbers:

| Engine | Profiles | Selected by | Size range | Shareable | What it is |
|---|---|---|---|---|---|
| srdif (`fft/srdif.h`) | `double`; `float` unless a backend is on | default | 4 … 2^30 | yes | a split-radix decimation-in-frequency kernel of length N/2 (Duhamel–Hollmann 1984; Sorensen–Heideman–Burrus 1986), the bit-reversal permutation and the real post-pass of Cooley–Lewis–Welch 1970 / Sorensen et al. 1987 (the fixed-point engine's structure in floating point), written from the literature; exactly the split-radix operation count, 2N log₂N − 2N − 2 real operations per forward transform; tables from integer arithmetic (no libm), built in the constructor, so without fp-contraction the output bits are one row for every compiler and target CI runs (pinned as output fingerprints in `tests/test_fft_srdif_fingerprint.cpp`). Since tap/DspTap#42, written clean-room; it replaced, at the same contract, the C++20 port of Ooura's `rdft` (`fft/split_radix.h`) that the floating profiles ran from Stage 2b (#31) until #42 (see Provenance) |
| CMSIS-DSP Helium (`fft/backends/cmsis.h`) | `float` | `TAP_DSP_FFT_CMSIS` (default ON exactly when the compiler targets floating-point Helium, `__ARM_FEATURE_MVE & 2`: the Cortex-M55 / M85 class; OFF for every other target, other Cortex-M cores included) | **32 … 4096** (CMSIS-DSP's own init table; the constructor checks the init status since Stage 4, in a debug build — outside this range the library never initializes, and a release build is undefined behaviour with no fault promised: measured under QEMU, N = 16 and 8192 give a wrong spectrum silently, N = 4 a wrong spectrum and a corrupted heap) | no | Arm's radix-4/8 MVE real FFT, re-presented in the same contract to float epsilon; on the M55 the srdif engine (the `m55-ooura` icount key) executes 1.60× (N = 512) / 1.77× (N = 2048) the instructions of the CMSIS build (the `m55` key) over the whole ratchet scenario, transform plus the class's copy loop, 2/N scaling and checksum (`bench/README.md`) |
| Apple vDSP (`fft/backends/accelerate.h`) | `float` | `TAP_DSP_FFT_ACCELERATE` (default ON on Apple) | 4 … 2^20 | no | `vDSP_fft_zrip`, same contract to float epsilon; ~3× faster per transform on Apple Silicon is MuTap's transform-only measurement against the engine the library carried before (tap/MuTap#31), not re-measured in this repo against srdif (the same-binary parity *test* exists since Stage 4; the microbenchmark needs a Mac) |
| int32 radix-4 (`fft/fixed_point.h`) | Q15, Q31 | the sample type (the second template argument is the scaling policy here) | 4 … 65536 | Q31 yes, Q15 no | one kernel over Q1.30 twiddles under two scaling policies, returning an exponent |

The two float32 backends are mutually exclusive, apply to `float` only (double
is always the srdif engine, the golden model), and conjugate imaginary
bins and rescale so every intermediate spectrum matches the default build to
single-precision rounding — so the whole double-precision test battery stays a
valid oracle for the accelerated float paths. `tests/test_fft_backend.cpp`
is typed over the engines the host can build and pins each to the srdif
float engine bin-for-bin at the certified geometries (N = 512, 2048) **in one
binary** — vDSP beside srdif on the macOS leg, CMSIS beside it on the M55
leg; `tests/test_fft_routing.cpp` pins `basic_real_fft` to the srdif engine
byte for byte on every leg (named explicitly) and to the selected default.

**Engine parameter and ABI tag (Stage 4 of the audit, tap/DspTap#35).** The
second template argument is the engine for `float` / `double` and the scaling
policy for Q15 / Q31: `real_fft32` is
`basic_real_fft<float, default_real_fft_engine_t<float>>`; an explicit
`basic_real_fft<float, detail::srdif_rdft<float>>` runs the srdif engine on a
build whose default is vDSP or CMSIS. `basic_real_fft<float |
double, scaling::fixed>`, the pre-Stage-4 spelling, resolved to the default
engine (as a distinct type) for one consumer cycle and cannot be instantiated
since the D4 expiry: a `static_assert` names the one-argument form
(tap/DspTap#40; API break, neither consumer writes it; pinned by the
`fft_compile_fail.*` ctests). The class re-exports the engine's `k_min_size` /
`k_max_size` / `k_is_shareable`, offers `supports_size(n)`, and requires
`supports_size(size)` at construction (`TAP_EXPECTS`: a debug assertion,
STYLE.md §4; `supports_size` is the release-mode query). The build's
selection also opens an inline namespace on `tap::dsp` — `fft_srdif`,
`fft_cmsis` or `fft_vdsp` — that `basic_real_fft`, `pvoc` and `log_mel` are
defined in: lookup is unchanged, but the mangled names of the classes whose
layout follows the selected engine now carry it, so two images built with
different defaults cannot coalesce each other's weak symbols (audit F4;
measured: 50/50 weak symbol names shared between the two builds before,
0/65 after, and in one process an embedder outside the tag ran the other
image's code while the same embedder inside it did not —
`tests/test_fft_abi_tag.cpp`; `docs/fft-design.md`, "Stage 4"). A consumer class that embeds the FFT by
value closes the same exposure by opening the same namespace in its own
(`namespace tap::mu::inline TAP_DSP_FFT_ABI` for MuTap, on its bump).

**Migration note for consumers (the srdif engine, tap/DspTap#42).** What
changed: the floating profiles' engine. `basic_real_fft<double>` and
`basic_real_fft<float>` (where no backend is selected) now hold a
`detail::srdif_rdft<Sample>`; every output bit of those profiles moves, at the
error level of either engine (the srdif engine's rms error against a
quad-precision reference is 4–14 % lower than its predecessor's at every N ≥
128, `fft/srdif.h` and `docs/fft-design.md`), so a consumer's own bit-exact
pins of floating FFT output are re-measured once; its tables take 1.7–1.8× the
memory per object (one allocation of 7.3 kB for a float N = 2048; `fft.h`
states the formula); and the ABI tag is
`fft_srdif` (it was `fft_split_radix`), so an image built against the new tree
does not coalesce with one built against the old. What did not change: the
contract (packing, sign, scale, size range 4 … 2^30, shareability,
allocation-free `noexcept` transforms, tables built in the constructor), the
API, and the Q15 / Q31 profiles, which do not use the floating engine. The
float-I/O-on-double overloads `forward(const float*, float*)` / `inverse(const
float*, float*)` were removed at tap/DspTap#40 (API break; neither consumer
calls them). `tap::dsp` sets no `-ffp-contract` flag and exports none
(Decision D9; the reasoning is in `fft.h`'s class docstring).

**Header-only.** `tap::dsp` is a pure INTERFACE target: no vendored C is
compiled into what ships. A static library (`tap_dsp_fft`, alias
`tap::dsp_fft`) exists only under `TAP_DSP_FFT_CMSIS`, carrying the CMSIS-DSP
objects, and `tap::dsp` links it automatically there. Until Stage 2c of the
audit the same library also carried a vendored C transform on every
platform; a reference copy then served a parity gate alone until Decision D6
deleted it.

Four profiles share the contract; the fixed-point ones (Stage 3b of the
audit, design in [`docs/fft-design.md`](docs/fft-design.md), "The
fixed-point profiles") return an exponent `e` from every transform in place
of a floating scale:

| Profile | Alias | Engine | Forward scale | Noise floor (fixed point: per-bin SNR, full-scale white noise, N = 512) | Target |
|---|---|---|---|---|---|
| `double` | `real_fft` | srdif | X | golden model: 1.60e-16 relative 2-norm error vs the compensated-DFT oracle at N = 256 on x86-64, 1.49 – 1.70e-16 on the Cortex-M legs (pinned at 7.4e-16, `DoubleForwardTracksCompensatedDft`) | desktop, reference |
| `float` | `real_fft32` | srdif (vDSP / CMSIS-Helium backends) | X | 9.73e-8 relative 2-norm error vs double at N = 512 on x86-64 and the soft-float M4, 1.05e-7 on the VFMA legs (pinned at 2.25e-7, `FloatEngineTracksDoubleAtN512`); < 1e-6 at N = 1024 | Cortex-M4F, M33, M55, Hexagon HVX, Apple Silicon |
| `std::int16_t` (Q15), `scaling::fixed` | `real_fft_q15` | int32 radix-4, Q1.30 twiddles | X / N, e = log2 N | 0.29 LSB rms (the output rounding); 66.5 dB at 0 dBFS, 26.6 dB at −40 dBFS | Cortex-M4 (soft-float), M33 |
| `std::int32_t` (Q31), `scaling::fixed` | `real_fft_q31` | same kernel, in place | X / 2N, e = log2 N + 1 | 0.68 LSB rms; 149.4 dB at 0 dBFS, 109.4 dB at −40 dBFS | Cortex-M4, M33, M55 |
| Q15, `scaling::block_floating` | `real_fft_q15_bfp` | same, per-stage headroom scan | X / 2^e, 0 ≤ e ≤ log2 N | 90.6 dB at 0 dBFS, 80.6 dB at −40 dBFS | round-trip consumers |
| Q31, `scaling::block_floating` | `real_fft_q31_bfp` | same | X / 2^e, 0 ≤ e ≤ log2 N + 1 | 160.6 dB at 0 dBFS, 155.7 dB at −40 dBFS | round-trip consumers |

Fixed-point contract in one line: with `G` the double profile on the same
input read as fractions of full scale, `G.forward == data · 2^e` and `G`'s
unnormalized `inverse_inplace == data · 2^e`; the fixed-point `inverse()`
applies no 2/N, and a round trip gives `x == out · 2^(e_fwd + e_inv + 1 − log2 N)`.
`fixed_scaling_exponent(N)` returns the fixed constant; every transform
returns `e` and is `[[nodiscard]]`. The int32 kernel performs no saturating
operation for any input under either policy (Q15 through two guard bits,
Q31 through a one-bit pre-shift under fixed scaling or the headroom rule
under block floating point; pinned by `SaturationFreeWorstCaseDoesNotWrap`);
the Q15 output narrowing clamps at the rail exactly when the true value is
the rail (a full-scale Nyquist alternation, pinned 0.5 LSB). Block floating
point at its full exponent agrees with fixed scaling within a pinned bound
(Q31 4 LSB at the DC index, Q15 1 LSB), bit-identically only when every
stage shifted the fixed amount. The Q15 fixed forward scale is the one
CMSIS-DSP documents for `arm_rfft_q15`; the rest of the convention is stated
in `fft.h` (Decision D3). Floors are the `[ floor ]` rows the battery's
`NoiseFloorTracksWelchModel` prints (output-referred against the double
golden model, N = 512, white noise); the Welch-model derivation and the full
level sweep are in the design note, and the output bit patterns themselves
are pinned per profile (`OutputFingerprintIsPinned`, N = 64 / 512 / 1024 /
2048). The complex kernel is DspTap's own; the real post-pass and pre-pass
around it are the half-length method of Cooley, Lewis and Welch (1970) and
Sorensen et al. (1987), one Q1.30 product per bin pair and exact DC /
Nyquist (design note §1). The deviation maxima the battery pins are
measured at N = 4 / 8 / 16 / 64 / 512 / 1024 / 2048; above that the Q31
block-floating maximum grows with the gap between its exponent and the fixed
constant (31 – 80 LSB at N = 4096 … 65536), so a host-only test pins it, the
F(x) + F(−x) asymmetry and the Q31 fixed maximum per size at 2× the measured
values (design note §3).

```cpp
#include "tap/dsp/fft.h"

tap::dsp::real_fft   fft(1024);   // double, the desktop/golden profile
tap::dsp::real_fft32 fft32(1024); // float,  the embedded / accelerated profile

std::vector<double> x(1024, 0.0);
fft.forward_inplace(x.data());        // packed spectrum, W = exp(+2πi/N)
fft.inverse(x.data(), x.data());      // out-of-place inverse, normalized (2/N)

tap::dsp::real_fft_q15 fft_q15(512);  // Q0.15 I/O, fixed scaling
std::vector<std::int16_t> q(512, 0);
const int e = fft_q15.forward_inplace(q.data()); // e == 9: spectrum == X / 512
```

Key contract points (full detail in the header docstring):

- **Packing** (N/2 + 1 bins): `data[0]` = DC real, `data[1]` = Nyquist real,
  `data[2k]`/`data[2k+1]` = bin *k* real/imag for `1 ≤ k < N/2`.
- **Sign convention** `W = exp(+2πi/N)` — imaginary parts are *conjugated*
  relative to the engineering-convention DFT. Consistent across every operand,
  so spectral products (fast convolution, adaptive-filter regressors) are
  unaffected; conjugate only when importing spectra computed elsewhere.
- **Normalization**: `*_inplace` inverse is unnormalized (multiply by `2/N`);
  the out-of-place `inverse()` applies the `2/N` for you.
- Transforms are `noexcept` and allocation-free after construction — real-time
  safe. Size must be a power of two inside the engine's range (the table
  above; `supports_size(n)`, the mandatory gate for a size that comes from
  configuration — a size outside it is undefined behaviour in a release
  build, with no fault promised), fixed at construction (Q15 allocates an int32
  work buffer of N at construction, Q31 transforms in place).
- One transform at a time per object unless `k_is_shareable` says otherwise
  for the engine in use (srdif and Q31: yes; vDSP, CMSIS and Q15: no).

## `tap::dsp::yin` — YIN pitch detector

`include/tap/dsp/yin.h` implements the time-domain YIN estimator (de Cheveigné
& Kawahara 2002, steps 1–5): squared-difference function, cumulative-mean
normalization, absolute threshold with local-minimum descent, and parabolic
interpolation for sub-sample period precision. Header-only, allocation-free
after construction, `noexcept` on the analysis path.

```cpp
#include "tap/dsp/yin.h"

tap::dsp::yin   det(800, 20, 800);   // window, tau_min, tau_max — double golden model
tap::dsp::yin32 det32(800, 20, 800); // float, the embedded profile

const auto r = det.analyze(frame);   // frame_size() == window + tau_max samples
if (r.voiced()) {
    const double freq = sample_rate / r.period; // fractional-sample period
}
```

Key contract points (full detail in the header docstring):

- **Geometry** fixed at construction: integration `window`, searched lag range
  `[tau_min, tau_max]` (bound from your frequency range as `τ = sr / f`), with
  `window ≥ tau_max`; `analyze()` reads `window + tau_max` samples, oldest first.
- **Result**: fractional period in samples (0 = unvoiced) plus the normalized
  aperiodicity at the chosen lag (global minimum when unvoiced).
- **Threshold**: the paper's absolute threshold on the normalized difference,
  default 0.1, settable at runtime.
- Both precisions run the identical algorithm; the test battery pins sine/
  sawtooth accuracy (sub-cent in double), octave robustness, unvoiced
  rejection, and float/double agreement. The difference-function inner loop is
  the designated Helium-MVE / HVX backend candidate behind this same contract,
  mirroring the FFT's golden-model-plus-backends pattern.

## `tap::dsp::psola` — pitch-synchronous overlap-add shifter

`include/tap/dsp/psola.h` is the real-time TD-PSOLA resynthesis stage:
Hann grains two source periods long, extracted at period-spaced analysis marks
and overlap-added at period/ratio-spaced synthesis marks with sub-sample
(Hermite) placement. Detection-agnostic — the caller supplies the period (from
`tap::dsp::yin` or any tracker); fixed latency of `2 * max_period + 2` samples.

```cpp
tap::dsp::psola shifter(900);            // deepest period it will be given
double y = shifter.process(x, period, ratio); // per sample; ratio 2 = octave up
```

Know what PSOLA is: it resamples the source's **spectral envelope** at the new
harmonic spacing — which is why it preserves formants on voice, and why a pure
tone shifted far from any new harmonic thins toward silence. Feed it
harmonic-rich material; both behaviors are pinned by the tests.

## `tap::dsp::pvoc` — phase-vocoder pitch shifter

`include/tap/dsp/pvoc.h` is an STFT pitch shifter (Hann, 4× overlap, built on
`tap::dsp::real_fft` so the float profile rides the vDSP/CMSIS backends) using
Laroche–Dolson-style **peak-region shifting**: each spectral peak's region is
translated rigidly by an integer bin offset and rotated by one accumulated
residual phase, so phase relationships across the peak stay intact. At
ratio 1 the output reconstructs the input's waveform delayed by exactly one
FFT frame (pinned by the tests). Latency = the FFT size (1024 default). The
size range is the class's own [64, 2^28] intersected with its FFT engine's
(`basic_pvoc<Sample>::k_min_size` / `k_max_size`, `supports_size(n)`): 64 …
2^28 for double on every build and for float on the srdif engine, 64 …
2^20 for float under vDSP, 64 … 4096 for float under CMSIS-DSP.

```cpp
tap::dsp::pvoc shifter(1024);
shifter.set_formant(true);               // optional LPC formant preservation
double y = shifter.process(x, ratio);    // per sample; ratio sampled per hop
```

Optional **formant preservation** (`set_formant`) uses the classic source-filter
method: an LPC spectral envelope (autocorrelation + Levinson–Durbin, order 48)
per analysis frame, with every relocated bin rescaled by
`envelope(target)/envelope(source)` — the excitation moves, the envelope stays.
At ratio 1 the correction is exactly unity, so the identity contract holds
either way.

## `tap::dsp::log_mel` — log-mel / PCEN feature front end

`include/tap/dsp/log_mel.h` is the streaming analysis front end of a keyword
spotter: windowed real FFT (on `tap::dsp::real_fft`, so the float profile
rides the vDSP/CMSIS backends), power spectrum, triangular mel filterbank,
then a floored affine `log10` or per-channel energy normalization (PCEN,
Wang et al. 2017). The header owns the **formula-level** contract as
numbers — HTK mel, unit-peak triangles, periodic Hann or sqrt-Hann, FFT
zero-padded at the end of the frame, frame *t* ending at sample
`(t+1)*hop - 1`, DC excluded, the PCEN recursion and its first-frame
priming — and stamps it with `k_contract_version`. Everything a trainer
tunes (band count, fmin/fmax, log floor/shift/scale, every PCEN parameter,
pre-emphasis) is runtime geometry in `log_mel_geometry`, carried by a trained
model, so retraining never touches the header. Reference geometry: 16 kHz,
400 / 160 / 512, 40 bands 20–7600 Hz. Latency = the frame length. A geometry
from configuration is gated by `basic_log_mel<Sample>::supports_geometry(g)`:
`g.valid()` (engine-independent; any power of two ≥ 4) AND the FFT engine's
`supports_size(g.fft_size)` — under CMSIS-DSP the float profile takes 32 …
4096 only, and a `valid()` geometry outside it (16, 8192) was measured
returning wrong features with no fault in a release build.

```cpp
tap::dsp::log_mel_geometry g;          // the reference geometry
g.pcen.enabled = true;                 // or leave the plain-log path
tap::dsp::log_mel32 fe(g);             // float embedded profile; log_mel is the double golden model
std::vector<float> feats(fe.frames_for(n) * fe.bands());
size_t frames = fe.process(x, n, feats.data(), fe.frames_for(n));   // any chunking, same features
```

Pinned by `tests/test_log_mel.cpp` against the committed numpy restatement
(`tools/reference/make_frontend_reference.py` → `tests/reference/frontend_vectors.h`
at the reference geometry and `frontend_vectors_tuned.h` at a geometry with
every runtime field off its default; the script's `Geometry` mirrors
`log_mel_geometry` field for field and is the only numpy copy of these
formulas in the family — MuTap's KWS feature module imports it through the
submodule for its parity self-check, and trains on the C++ front end through
the C ABI): both paths sample-for-sample at both geometries, chunking
invariance, alignment and latency, the filterbank formulas, PCEN's gain
tracking on a level step and its reset semantics, and float/double agreement
as a measured number (6.7e-7 log, 5.3e-6 PCEN, pinned at 2×). No double
arithmetic on the float path: the RP2350's Cortex-M33 has no FP64.

## `tap::dsp::decimate` — fixed-ratio decimators to 16 kHz

`include/tap/dsp/decimate.h` is the host-rate stage in front of the 16 kHz
front end: `basic_decimator<S, M>` (S a `sample_type`) for M = 2, 3, 6
(32 / 48 / 96 kHz in), in the `bridge` engine's pattern — ratio as a type, Kaiser-windowed sinc from
`kaiser.h` with the cutoff at the output Nyquist and DC gain exactly 1,
`fir_kernels.h`'s `dot_row` over the `sample_traits.h` formats (double golden,
float embedded, Q15 / Q31 with row-sum-preserving quantization). It is deliberately not
`tap::sr::bridge`, whose charter is 44.1 ↔ 48 only; a 44.1 kHz host composes the
`bridge` engine's 44.1 → 48 in front of the by-3 stage. Odd tap counts, integer group delay
`(taps - 1) / 2`, one output as input `k*M` arrives.

| profile     | stopband | passband | taps by 2 / 3 / 6 |
|-------------|----------|----------|-------------------|
| economy     |  70 dB   | 7000 Hz  |  81 / 121 / 239   |
| transparent | 100 dB   | 7600 Hz  | 259 / 389 / 773   |

```cpp
tap::dsp::decimate_by_3 dec;                       // 48 kHz -> 16 kHz, economy, float
std::vector<float> y(dec.outputs_for(n));
dec.process(x, n, y.data());                       // noexcept, allocation-free, chunking-invariant
```

Pinned by `tests/test_decimate.cpp`: the tap counts against the searched
minima, the float output sample-for-sample against the numpy reference,
the passband and stopband numbers measured from the shipped coefficients,
unity DC (exact in Q15), the group delay, and Q15 tracking float within the
format's floor. Q15 stores `M·h` (peak tap 1.0) rather than `h` (peak
1 / M, a sixth of the Q1.14 range at by 6) and divides the M back out in
the single rounding (`finalize_divided`, below), so its quantized tables
attain the economy stopband: −72.3 / −71.7 / −71.0 dB by 2 / 3 / 6, where
`h` itself attained −68.8 / −66.6 / −61.1 dB.

## `tap::dsp::nn` — dense and GRU inference kernels

`include/tap/dsp/nn.h` is the arithmetic of a small recurrent gain network:
`basic_dense<Sample>` (`y = act(W x + b)`, row-major `[out x in]`, linear /
tanh / sigmoid) and `basic_gru<Sample>` (Cho et al. 2014 in PyTorch's
`nn.GRU` convention, gates ordered r, z, n along the `3*hidden` axis). Lifted
from MuTap's learned residual suppressor at its wake-word plan's M3; the
keyword spotter is a different head on the same layers. Weights are stored
as float32 whatever the profile, the storage precision of a trained model's
file, and converted to `Sample` at the point of use; every dot product
accumulates in `Sample` bias-first in ascending input order, so the float
profile contains no double arithmetic. Weight vectors are moved in and
owned (a copied layer is a deep copy); `apply()` / `step()` are noexcept and
allocation-free. `k_contract_version` stamps the formulas.

```cpp
using namespace tap::dsp::nn;
// Weights are taken by value and moved in: std::move hands the loaded vectors over without a copy.
dense32 din(std::move(w_in), std::move(b_in), 64, 56, activation::tanh); // float embedded profile; dense/gru are the double golden model
gru32   cell(std::move(w_ih), std::move(w_hh), std::move(b_ih), std::move(b_hh), 96, 64);
dense32 dout(std::move(w_out), std::move(b_out), 26, 96, activation::sigmoid);
din.apply(features, hidden);  cell.step(hidden);  dout.apply(cell.state(), gains);
```

Pinned by `tests/test_nn.cpp`: the layout on hand-computed numbers, the
activation forms, the gate order by isolating each block through its
biases, the GRU formula against an independent long-double restatement
that sums in the opposite order, reset and copy semantics, noexcept, and
float-tracks-double as a measured number at the suppressor's geometry. The
end-to-end oracle stays in MuTap: its Python parity CI job and its
suppressor's cross-precision pin must be unchanged by the promotion.

## `tap::dsp::subband` — the sub-band round trip

`include/tap/dsp/subband.h` runs a process at a lower rate, where an FFT of
the same size resolves D times finer, and adds back only what the process
*changed* (SampleRateTap's `docs/X4_X8_PLAN.md`, 6.3 and X10):

```text
out[n] = x[n − T] + up(change(down(x)))[n],   T = D·τ_down + D·L_op + τ_up
```

with τ_down the decimator's delay in low-rate frames, L_op the op's latency
in whole low-rate frames and τ_up the interpolator's delay in full-rate
frames. The scaffold holds the full-rate delay line, takes T from its stages
and runs the low-rate op. It is generic over `chain.h`'s `sync_stage`: Down
converts by 1/D and Up by D/1 (inverse integer factors, a `static_assert`),
each a single stage or a `tap::dsp::chain`, and the sample type and channel
count come from the stages. Everything the op leaves alone reaches the
output through the delay line untouched, so the decimator's aliasing in the
bins the op does not touch is never re-injected, no guard band filter is
needed, and bands compose by addition with no crossover filter (Crochiere
and Rabiner 1983; Vaidyanathan 1993; Allen 1977 for an STFT op's overlap-add
modification).

- **Where T comes from** (`exact_delay_stage`): each stage reports an exact
  latency (`latency_output_frames()` or `latency_input_samples()`, as
  `chain.h` reads them) or an exact delay path
  (`delay_path_output_frames()`, the `reports_delay_path` concept: the member
  a near-linear recursive stage reports, so admitting one needs no DspTap
  change), and a chain is walked stage by stage. A stage with neither, such
  as a nonlinear-phase recursive stage, does not compile. A non-whole
  D·τ_down + τ_up is a `static_assert` when both stages' delays are
  compile-time constants, and otherwise is rejected at construction through
  `tap::dsp::raise<std::invalid_argument>`, with `supports(down, up)` as the
  predicate to check first. An odd-length L-th-band FIR delays mL − 1 samples
  at its higher rate, so a chain of integer-factor `nyquist.h` stages is
  always whole.
- **The op** (`subband_op`): a `sample` type and `channels()` as the stages
  declare them, `latency_frames()` (L_op), `process(in, n, change)` writing n
  frames of change (`change[k] = result[k] − in[k − L_op]`, exact zeros for no
  change), `reset()`, and `window_frames()` if the scaffold is to flush. Its
  change must stay in the bins up to p·r_low, the stages' passband (X15: a
  precondition on the op, not enforced; masking belongs to the op).
- **Identity, by construction:** an op that changes nothing gives x[n − T]
  bit for bit, for every bit pattern: where the change is exactly zero the
  delayed input is written as it is (an addition would turn −0 into +0).
- **Composition:** every band is aligned to A, the largest T
  (`round_trip_latency(down, up, op)` gives each band's T before
  construction). A band constructed with `aligned_latency = A` holds its
  extra A − T as ⌊(A − T) / D⌋ frames at its low rate, before Up, and the
  (A − T) mod D rest at the full rate. One band writes
  `subband_output::input_plus_change` and holds the delay line; the others
  write `subband_output::change_only` and hold none, and the caller adds them.
- **Memory (X12):** the delay line holds A·Ch samples: with a 4096-point op
  on SampleRateTap's economy chains, 129 / 258 / 515 KiB of double per channel
  at D = 4 / 8 / 16. A D = 4 band beside a D = 8 band adds its 16480 extra
  frames as 4120 at 12 kHz (32 KiB in double). The rest is a few blocks of
  scratch (the header states each as a formula).
- **Real time:** construction allocates and may reject; `process()`,
  `flush()` and `reset()` are noexcept and allocation-free; any chunking
  writes the same bits, and in-place running (in == out) is allowed.
  `flush()` feeds `window_frames()` zeros, a bound valid from any state, so
  it is zero padding by construction.
- **Formats:** float and double, from the stages. Fixed point is out of scope
  (a `static_assert`): the sum and the delay line would need a Q-format design
  of their own.

```cpp
#include "tap/dsp/subband.h"
// Down, Up: the caller's stages by 1/D and D/1, a stage or a tap::dsp::chain each
// (SampleRateTap's rational::down_by8 / up_by8, say); Op: a subband_op.
tap::dsp::subband_delta<Down, Up, Op> band(std::move(down), std::move(up), std::move(op));
band.process(in, n, out);                 // out[k] = in[k - T] + the change; noexcept, n frames
const std::size_t t = band.latency();     // T, a whole number of full-rate frames
```

| profile | removal floor, worst over the usable band (test stages, 70 dB) | float vs double |
|---|---|---|
| `double` | −65.0 dB (D = 4), −61.3 dB (D = 8) | the golden model |
| `float` | the same within 0.01 dB (0.2 dB at a near-null) | 6.5e-8 max abs (pinned at 1e-7) |
| Q15 / Q31 | out of scope (a `static_assert`) | |

Pinned by `tests/test_subband.cpp`, typed over float and double, on
test-local stages built on `nyquist.h` at 70 dB, as `chain.h` was proven on
`decimate.h` (D = 4 by one 4th-band stage each way, T = 86; D = 8 by ↓4·↓2
and ↑2·↑4, T = 206): the identity bit for bit, −0, subnormals and the
largest finite values included, for one stage, chains, two channels, an op
latency and an aligned latency; a removal at the stages' ripple floor
(−61.3 to −116 dB at the tones removed, the same figures as the designs'
responses predict to 0.01 dB) and failing by at least 46 dB with T − 1 or
T + 1 (2 sin πf), which is what pins T against the signal path, while a tone
above the band passes unchanged; T as the exact sum of the parts' delays
from every form, mixed forms in one chain included; any chunking and in-place
running bit-identical; flush equal to zero padding with both paths drained;
three bands of different T (206, 86 and 61 with `decimate.h`'s by-2) composed
by addition, each change-only band's output its own change delayed by
A − T bit for bit and each band's removal at its floor inside the
composition; an unmasked change re-injecting a transition-band alias at
−21.4 dB (X15's honest limit); a delay-path stage accepted bit-identically
to a latency stage; a stage with neither and a chain holding one refused
(concept negatives); a lagging change path and a non-whole delay rejected
at construction; and float tracking double. `tests/compile_fail/subband_non_integer_delay.cpp`
pins that stages whose delays are compile-time constants and sum to a
fraction of a sample fail on the class's `static_assert`, as the build's
only error. The battery runs whole on every QEMU leg (M4 soft-float, M4F,
M33 and M55, and the M33 and M55 without exceptions or RTTI): nothing is
excluded.

## The recursive substrate — two-path half-bands and the guard band filter

The design math and kernels of the SampleRateTap family's recursive (IIR)
polyphase stages (its `docs/IIR_PLAN.md`: I1, I4, I5, I8, I9, I12, I13,
I14), landed here first: `rational` holds only the stage types. Four
headers over `detail/cx_math.h`, plus `chain.h`'s latency members (below):

- `tap/dsp/iir_design.h` — the designs. H(z) = ½ [A0(z²) + z⁻¹ A1(z²)], a
  half-band of a by-2 stage's higher rate (Valenzuela & Constantinides 1983;
  Vaidyanathan 1993), in two variants: **nonlinear**, both paths cascades of
  first-order allpass sections (the elliptic, power-complementary half-band,
  Gazsi 1985's explicit formulas), and **near-linear**, A0 a pure delay
  z⁻²ᴷ and A1 an allpass approximating it (Renfors & Saramäki 1987's
  structure); the **four-path** stage ¼ Σ z⁻ᵏ Aₖ(z⁴) (harris 2004); and the
  guard band filter's **odd-order Chebyshev II** (Antoniou 2006) and its
  **linear form**, an even-length symmetric (type II) FIR.
- `tap/dsp/iir_tables.h` — the committed designs the compiler cannot
  compute (iterative optimizations), generated by `tools/iir/gen_iir_tables.py`
  in a venv pinned by `tools/iir/requirements.txt` (Python 3.11, numpy
  2.4.6, scipy 1.17.1, mpmath 1.4.1).
- `tap/dsp/iir_kernels.h` — the kernels per `sample_traits` format, and
  `tap/dsp/subnormals.h`, the flush-to-zero precondition.

**Spec points.** A design point is (A dB, passband edge p as a fraction of
the chain's lowest rate) at a design divisor d (the stage's lower rate over
the lowest; IIR_PLAN I13): a stage at d has f_pass = p / (2d) of its higher
rate and its stopband mirrored about its own low Nyquist, and every design
meets −(A + 1 dB) there. The four tabulated points are the family's
profiles' numbers, keyed by the numbers, never by the names
(`k_nonlinear_halfband<70, 3, 8, 2>`, `k_near_linear_halfband<…>`,
`k_fourpath<…>`, `k_guard_band_cheby2<A, num, den, order>`,
`k_guard_band_fir<…>`):

| A, p | (rational's name) | nonlinear order, d = 1 / 2 / 4 / 8 | near-linear N (= K), d = 1 / 2 / 4 / 8 | four-path, d = 2 |
|---|---|---|---|---|
| 70 dB, 1/3 | `super_economy` | 9 / 5 / 5 / 3 | 5 / 2 / 2 / 1 | 6 sections, −89.6 dB |
| 70 dB, 3/8 | `economy` | 9 / 5 / 5 / 3 | 7 / 3 / 2 / 1 | 6 sections, −84.2 dB |
| 70 dB, 19/48 | `balanced` | 9 / 5 / 5 / 3 | 8 / 3 / 2 / 1 | 6 sections, −81.7 dB |
| 120 dB, 5/12 | `transparent` | 17 / 9 / 7 / 5 | 19 / 5 / 3 / 2 | — (I12) |

**Compile-time design (I4).** The nonlinear half-band and the Chebyshev II
are `constexpr` over `tap/dsp/detail/cx_math.h` (sqrt, sin and cos of πt,
atan, exp, expm1, log, log1p, sinh, cosh, asinh, acosh, 10ˣ − 1 in plain
C++20 double, no libm: each within 0.5 to 2.5 ulps of an mpmath leg, pinned
by `tests/test_cx_math.cpp`). Its order comes from the closed form,
10 log₁₀(1 + q⁻ᴺᐟ²/4) ≥ A + 1 with q the nome, and reproduces IIR-0's grid
search at every tabulated point. Three steps are rewritten by algebra to
remove a cancellation (the nome step, IIR-0's amendment; a = (1 − x²)/(1 + x)²;
1 − w²/k as a difference of squares), and every tabulated coefficient is
within **13.2 ulps** of mpmath at 60 digits (IIR-0 promised 15; the prototype
lost up to 94,049). The largest evaluation takes 11,090 of Clang 18's
constant-evaluation steps (the order-17 half-band; the order-11 Chebyshev II
4,749; the exact product in integers roughly doubled both), against MSVC's
documented default of 100,000. The bits are pinned
per design (`tests/test_iir_constexpr.cpp`, IIR-0's bit test moved in), one
row for every compiler CI runs, srdif's rule. A custom spec calls the same
functions at run time. The Chebyshev II's L-inf gain distribution (Jackson
1970) left constexpr as IIR-0 recommended (415k steps at order 11), its
factors committed per spec point and order — and every one is exactly 1.0:
with the sections ordered by increasing pole radius, each at unity DC gain,
every partial cascade of every committed design peaks at DC
(`IirDesign.Cheby2LinfGainsAreUnityAtEveryCommittedPoint`). The linear form's
custom setting designs at run time through `kaiser.h`
(`design_guard_band_linear`, 52 taps where remez needs 40 at 70 dB, 3/8,
order 5's edge). Committed: 16 near-linear half-bands (I13's divisors 1, 2,
4, 8), 3 four-path stages, 20 Chebyshev II gain rows (orders 3 to 11),
19 type II FIRs (remez does not converge at 120 dB order 3, as in IIR-0) and,
since 0.4.0, 15 Q15 rows of the type II FIR (below).

**The linear form in Q15 (0.4.0).** The double type II rows do not survive
Q1.14: quantized as the family quantizes a symmetric row
(`quantize_symmetric_row_preserving_sum`, the first half with its sum
preserved, then mirrored), the 70 dB rows attain −63.5 to −74.8 dB from their
stop edges, and 13 of the 15 miss 70 dB or 0.01 dB of ripple; SampleRateTap's
guard band filter measured −67.1 / −66.8 / −64.6 dB at the three it runs
(`IirDesign.GuardBandFirDoubleRowsAreFormatLimitedInQ15`). So Q15 has rows of
its own, `k_guard_band_fir_q15<A, num, den, order>` (a
`guard_band_fir_q15_design`: the double design to quantize, its stop edge, the
target it was designed for, and the stopband and ripple its Q1.14 taps
attain), at the three 70 dB points and the orders 3 to 11. Each is the
shortest even length at which a Parks–McClellan design (McClellan, Parks &
Rabiner 1973, through `scipy.signal.remez`), targeted at A to A + 30 dB in
0.01 dB steps (one lever: in remez only the band weights' ratio matters),
keeps the double row's spec from the same stop edge once quantized so, on
16384-point grids over the passband and the stopband; at that length, the
design whose Q1.14 taps go deepest. The generator emulates the quantizer bit
for bit (`tools/iir/iir_design.py`), and the battery quantizes the committed
rows with DspTap's own and checks the Q1.14 taps against the emulation's
(`tests/reference/guard_band_q15_taps.h`) before measuring them. Q15 taps (the
double row's), then the Q1.14 taps' stopband and ripple:

| A, p | order 3 | order 5 | order 7 | order 9 | order 11 |
|---|---|---|---|---|---|
| 70 dB, 1/3 | 20 (20): −75.52 dB, 0.00965 dB | 32 (30): −89.23, 0.00898 | 40 (38): −70.40, 0.00984 | 54 (52): −72.82, 0.00939 | 72 (66): −73.12, 0.00982 |
| 70 dB, 3/8 | 26 (26): −76.45, 0.00970 | 42 (40): −90.87, 0.00895 | 50 (50): −72.33, 0.00966 | 72 (64): −74.60, 0.00929 | 96 (86): −75.01, 0.00975 |
| 70 dB, 19/48 | 32 (32): −77.01, 0.00816 | 48 (48): −71.71, 0.00983 | 60 (58): −70.66, 0.00964 | 86 (76): −75.93, 0.00973 | 114 (102): −70.88, 0.00998 |

The search is over rounding outcomes as much as over designs: at the
selected lengths 2 to 2,031 of the 3,001 targets pass (2 for 1/3 order 7 and
for 19/48 order 5), a 0.05 dB grid of targets finds the same lengths, and a
0.25 dB grid finds five of them 2 to 4 taps longer. There is no row at
120 dB: rounding each tap to Q1.14 (an LSB of 2⁻¹⁴, −84.3 dB) adds an
error response whose RMS over frequency is 2⁻¹⁴ √(N/12), −76.5 dB at 72 taps,
and the 120 dB rows quantized so attain −97.4 / −70.9 / −82.1 / −69.2 dB at
orders 5 / 7 / 9 / 11. As SampleRateTap's `rational` states for its FIR stages
(its PLAN.md section 6), Q15 at 120 dB is format-limited, and the 120 dB
promise is float's and Q31's: double, float and Q31 run `k_guard_band_fir`,
whose spec Q1.30 keeps (`IirDesign.GuardBandFirDoubleRowsKeepTheirSpecInQ31`).

**Where the kernels stop** (`fir_kernels.h`'s line): a kernel computes the
outputs of one commutator step for one channel from the inputs its paths
take — `halfband_decimate(table, state, x[2k], x[2k−1])`,
`halfband_interpolate`, `near_linear_decimate` / `_interpolate` (the caller
passes the delay path's x[2(k−K)]), `fourpath_decimate` / `_interpolate`,
`biquad_cascade` (1:1) — and updates the recursive state the caller stores
per channel. One kernel takes a block: `biquad_cascade_block<S, K, Ch>` runs
the floating cascade over n interleaved frames of Ch channels with a
compile-time section count, the coefficients and each channel's state held
in locals for the call. It is the per-sample kernel's expressions on the
same state, so with fp contraction off its output and state are the
per-sample kernel's bit for bit, for any chunking and with the two mixed
(`tests/test_iir_block.cpp`); what it saves is the per-sample kernel's
reload of every coefficient and state word on every sample (0.26 to 0.52 of
the instructions at orders 2 to 8 in the Butterworth stage on the embedded
cores with an FPU, `bench/README.md`). The commutator schedule, the odd-sample latch, the delay path's
history, the stream phase, channel interleaving and flush are the stage
machine's (`rational`'s, IIR-b), as a FIR stage keeps its history and calls
`dot_row`. Coefficient tables in each format (`nonlinear_halfband_table<S>`
and its siblings) are `constexpr`, so a named design's Q15 / Q31 tables round
from the same compile-time doubles.

| format | states | coefficients | per section | output | measured (`tests/test_iir_kernels.cpp`, `test_iir_fixed.cpp`) |
|---|---|---|---|---|---|
| `double` | double | double | — | identity | the golden model: −309 dB from the direct form at worst |
| `float` | float (I14) | float | — | identity | −137.9 dB from double at worst (120 dB, near-linear, white noise; IIR-0 −137.1), −142.6 dB at 70 dB; the guard band cascade −125.4 dB at 120 dB order 11 |
| Q15 | int32 Q2.29 (two guard bits) | Q0.31 (allpass), Q2.30 (second order) | one rounding: floor truncation (one-multiplier allpass), magnitude truncation (every second-order section) | Q2.29 → Q15 once, round half up, saturating | largest state 2.06 of full scale under full-scale square waves and clipped sines; outputs within 0.5 LSB of double; floor −98.1 dB for a −1 dBFS sine, equal to the output's own rounding |
| Q31 | int64 Q2.60 | as Q15 | the same, products exact to 96 bits | Q2.60 → Q31 once, round half up, saturating | largest state 2.06; outputs within 0.5 LSB of double on the same coefficients |

Biquads run in direct form I in fixed point and transposed direct form II in
floating point; the difference nodes (x − y1, which reach 5.22 in L1 at 120
dB) are formed in 64 bits and never stored; sums run in unsigned arithmetic,
so a state pushed past the headroom wraps rather than invoking undefined
behaviour; saturation happens at the output only. Silence after full-scale
drive reaches an exactly zero fixed-point state within 248 (Q15) and 502
(Q31) low-rate frames for the half-bands (bounds 330 and 659, computed from
each section's longest floor-truncated decay), 181 / 351 for the near-linear
ones and 241 / 535 for the guard band cascades. The Q15 and Q31 outputs are
bit-pinned (FNV-1a-64, one row everywhere).

**Flush to zero is the host's job (I9).** No kernel writes or restores the
floating-point control register; each structure kernel reads it as a
debug-build precondition, `TAP_DSP_EXPECTS_SUBNORMALS_FLUSHED(S)`
(`TAP_EXPECTS(subnormals_flushed<S>())` where a bit governs S, nothing
elsewhere, compiled out by `TAP_DSP_ALLOW_SUBNORMALS`, which like
`TAP_EXPECTS` must be one setting per image). `subnormals_flushed<S>()` and
`k_subnormal_control<S>` are public, for a host and its tests:

| target | what the host sets | governs |
|---|---|---|
| x86-64 | MXCSR FTZ (bit 15) and DAZ (bit 6): `_mm_setcsr(_mm_getcsr() \| 0x8040)` on each audio thread | float and double |
| AArch64 (Apple silicon, Arm Linux, Windows on Arm) | FPCR.FZ (bit 24): `msr fpcr` (GCC, Clang), `_WriteStatusReg` (MSVC), on each audio thread | float and double |
| Cortex-M with an FPU (M4F, M33, M55), AArch32 | FPSCR.FZ (bit 24), and FPDSCR (0xE000EF3C) for exception handlers: `__builtin_arm_set_fpscr` or CMSIS's `__set_FPSCR`, at reset | float on the M4F and M33; float and double on the M55 |
| Cortex-M without an FPU (M4), double on the M4F and M33 | nothing: software float | — |
| Hexagon | nothing: USR has no flush-to-zero field (the vendor's V68 Programmer's Reference Manual, 80-N2040-46 Rev. B, Table 2-6; its floating point "defined to support the IEEE floating-point standard", 4.3.4), and QEMU models none | — |

Honest limits, pinned: a host that does not set the bits (or a target with
none) leaves the states at a residue of a few times the smallest subnormal
(2 at 70 dB 3/8, 5 to 7 at 120 dB for the half-bands; 4 to 20 for the Chebyshev
II cascades; `tests/test_iir_subnormals.cpp`), and on x86 silence after
signal runs about ten times slower; under the bits the allpass sections reach
an exactly zero state within `floating_flush_frames()` (I6's count: the
largest pole radius from the format's largest value to its smallest normal,
plus one frame per section), with one exception: the floating near-linear
half-band at 120 dB, 5/12, divisor 1 (N19, nine second-order sections) holds a
zero-input limit cycle at 3.7 to 4.1 times the smallest normal, in float and
double, from 12 to 20 of 39 noise starts (SampleRateTap's IIR-b; `rational`'s
`iir_stage.h` states it and pins its bound). The mechanism is the
**floating Chebyshev II cascade's, which at order 7 or more can hold a
zero-input limit cycle at a few times the smallest normal** (up to 13.1 times
measured, in either direct form: once |a2 y| underflows its product flushes
and the section briefly runs with |a1| > 1), so a floating guard band filter
at those orders never drains to an exactly zero state. Both are honest limits,
kept (SampleRateTap's `docs/IIR_PLAN.md` question 19: no deadband in the
kernels, no design dropped); a consumer states the residue its stage leaves
after its flush window and pins a bound on it, as `rational`'s guard band
filter does at order 7. A non-finite input poisons a
recursive state until the caller zeroes it (I9, decided). And a consumer
built with `-ffast-math` (or `/fp:fast`, or icx's default `-fp-model=fast`)
gets none of these contracts: reassociation and contraction move the float
outputs past the tolerances, `-ffinite-math-only` deletes `std::isnan` and
`std::isfinite`, and an executable linked with it starts with FTZ and DAZ on
(GCC's `crtfastmath.o`; before GCC 13 a shared library built with it did the
same to its host) — the one effect that does what I9 asks of a host.

Pinned by `tests/test_iir_design.cpp` (the orders, ulps, stopbands on the
16384-point grid, power complementarity, group delays against IIR-0's, the
committed rows' promises, the type II FIRs' symmetry and exact zero at
Nyquist; the Q15 rows' lengths, their Q1.14 taps bit for bit against the
generator's emulation, their zero at Nyquist and their attained stopband and
ripple, and the double rows' loss of the spec in Q1.14 and keeping of it in
Q1.30), `tests/test_iir_constexpr.cpp` (the bits, its own target at
`-ffp-contract=off`, where the run-time evaluation equals the compile-time
one too), `tests/test_iir_kernels.cpp` and `test_iir_fixed.cpp` (above),
`tests/test_iir_subnormals.cpp` (its own executable under
`TAP_DSP_ALLOW_SUBNORMALS`) and `tests/test_chain_iir.cpp` (below). Every
case runs on every QEMU leg except nine design sweeps of `test_iir_design.cpp`
(minutes of soft-float emulation on the M33, design-time numbers: 141 s for
the six of 0.3.0, about 1,100 s for the three Q15 measurements of 0.4.0; the
`MAIN_FILTER` names them). The Q15 rows' structural pins, bit-exact Q1.14
taps included, run on every leg. The FIR substrate's pins and the srdif
fingerprints do not move.

## `tap::dsp::butterworth` — the Butterworth low-pass stage

`include/tap/dsp/butterworth.h` (since 0.6.0) designs an even-order
Butterworth low-pass (order 2 … 16: the bilinear transform with pre-warping,
Oppenheim & Schafer sec. 7.2) and runs it as a 1:1 `chain` stage over `Ch`
interleaved channels: `basic_butterworth_lowpass<Sample, Order, Ch>`, float or
double, with `butterworth4<Ch>` / `butterworth4_32<Ch>` and `butterworth8<Ch>`
/ `butterworth8_32<Ch>` for orders 4 and 8. It runs `iir_kernels.h`'s biquad
cascade (the [recursive
substrate](#the-recursive-substrate--two-path-half-bands-and-the-guard-band-filter)
above) per channel, under SampleRateTap's IIR plan's contracts for a recursive
stage in a chain (`docs/IIR_PLAN.md`: I5 latency, I6 flush, I9 subnormals).

```cpp
tap::dsp::basic_butterworth_lowpass<float, 4, 2> lp(0.1); // stereo, -3 dB at 0.1 fs
lp.process(in, frames, out);                             // interleaved, n in, n out
double d = lp.nominal_latency_output_frames();           // group delay at DC: delay a dry path by this
double e = lp.latency_deviation_output_frames();         // how far the passband strays from d
lp.set_cutoff(0.05);                                     // redesign in place: noexcept, allocation-free
```

- **Design:** unity at DC, −3.0103 dB at the cutoff, a double zero at
  Nyquist per section, monotone between; the sections in increasing pole
  radius (`design_butterworth_lowpass<N>(fc)` returns them). In exact
  arithmetic each section is the Audio EQ Cookbook low-pass biquad
  (Bristow-Johnson) at Q_k = 1/α_k, the standard Butterworth Q set
  (0.5412 and 1.3066 at order 4): a cookbook cascade at those Qs printed to
  eight decimals is this design within 2.5e-9 per coefficient (measured,
  `SectionsAreTheCookbookBiquadsAtTheButterworthQs`), so such a cascade can
  be replaced by the stage.
- **Kernel:** `process()` runs `iir_kernels.h`'s `biquad_cascade_block`:
  the cascade a block at a time, two channels together in float at orders 2
  and 4 and one at a time otherwise (the grouping the cores' registers
  hold, measured), and a single frame at order 8 and above through the
  per-sample kernel (the block's copies would cost more than they save
  there). The paths run the same expressions, so the output is bit-identical
  for any chunking and channel count under `fft.h`'s fp-contraction policy
  (contraction off). With contraction on, a compiler may fuse the paths
  differently on a target with a fused multiply-add, and the last bits can
  then depend on the chunking; through 0.6.0, on the per-sample kernel alone,
  they could not. `in` and `out` are the same buffer or disjoint ones.
- **Cutoff range:** [1e-5, ½ − 1e-5] cycles per sample in double and
  [1e-4, ½ − 1e-4] in float (`k_butterworth_min_cutoff<S>`,
  `butterworth_supports_cutoff<S>`); the constructor rejects anything else
  (`validate()`, then `raise<std::invalid_argument>`). Float's bound is
  where its coefficient rounding stops threatening the poles: each
  section's stability margin, 1 + a2 − |a1|, is 4K²/a0 at the low end, and
  rounding a1 and a2 to float takes up to 9e-8 of it, so below fc = 4.7e-5
  a float pole reaches the unit circle (measured, every even order, at
  either end; `FloatRoundingKeepsThePolesInsideDownToTheFloatBound`).
  Double's keeps the flush window in a 32-bit `size_t`: the window grows as
  1/fc toward either end, up to 232,710,082 frames at order 16.
- **`set_cutoff(fc)`:** the constructor's design, in place: the sections,
  the kernel's table, the flush window and the nominal latency are bit for
  bit what a stage constructed at fc holds, and the deviation it reports is
  the same number (the constructor runs the same function; bit for bit under
  `fft.h`'s fp-contraction policy, as every bit-exact claim here). `noexcept`
  and allocation-free (pinned through the heap meter); its precondition,
  `butterworth_supports_cutoff<S>(fc)`, is a `TAP_EXPECTS`, never a throw. It
  leaves the state alone: the recursion continues under the new
  coefficients, and a caller that wants a clean start calls `reset()`. Its
  cost is one design, nearly all of it the flush window: 1.1 to 2.8 µs at
  order 4 and 5.5 to 14 µs at order 8, float or double, at the cutoffs
  0.45/2 … 0.45/16 an oversampling kernel sets (x86-64, a shared 2.1 GHz
  Xeon VM, GCC 13.3 `-O3`, medians of 401 calls over three runs), and 77 µs
  at order 16 in double at its lowest cutoff. The design is double, so a core
  without a double FPU pays libgcc for it: 2.3 M instructions a call at order
  8 on a Cortex-M33 or M4F (2.8 M at most; 0.1 M on a Cortex-M55,
  `bench/README.md`). Through 0.6.0 it also ran the deviation's group-delay
  scan, 80 to 93 % of the call (35.7 M instructions at order 8 on the M33);
  since 0.7.0 the deviation is computed when it is read. Its stack is the
  flush window's matrices, GCC's `-fstack-usage` frames 1,432, 5,272 and
  20,640 bytes at orders 4, 8 and 16 on the Cortex-M33. It is for occasional
  changes — TapTools' kernels redesign their oversampling filters at
  0.45 / os from `set_oversample()` — not for per-sample or per-block
  modulation.
- **Latency (I5):** no exact latency, since a recursive stage has no single
  delay. It reports three numbers in frames (the double design's):
  `group_delay_output_frames(f)`; `nominal_latency_output_frames()`, the DC
  group delay, the delay to give a parallel path (a dry signal, a sibling
  channel) to line it up; and `latency_deviation_output_frames()`, the
  **passband group-delay deviation**, the largest |group delay − nominal| up
  to the cutoff on a 1,024-point grid: the alignment error left after
  compensating the nominal alone (order 4: 2.51 frames on a nominal of 4.02
  at fc = 0.1, 2.39 on 1.31 at 0.25; a Butterworth's group delay peaks at the
  band edge). Since 0.7.0 the deviation is computed on every read, not stored
  (`butterworth_latency_deviation(sections, cutoff)`, the same scan as a free
  function for a caller that keeps the number): 1,025 group-delay evaluations
  in double, 40 to 72 µs at orders 4 and 8 on the x86-64 above and about
  28 M instructions at order 8 on a Cortex-M33, so it is a design-time query,
  not one for the audio thread; it is `const` and has no cache, so threads
  may read it at once. `chain` composes all three with its exact stages.
- **Flush (I6):** `window_frames()` is the count of zero frames after which
  the exact recursion of the coefficients the kernel runs (the design
  rounded to the format) is below the smallest normal value from any state:
  the smallest n with ln(max) + ln ‖Aⁿ‖∞ < ln(min_normal), A the cascade's
  state-transition matrix, plus a frame per section, found by binary lifting
  over A^(2^j) with each power recomputed from A, so the stack holds a few
  matrices whatever the window (780 frames in float and 6,207 in double at
  order 4, fc = 0.1). What the arithmetic leaves after it is the host's:
  usually exact silence; without FTZ a subnormal residue; under FTZ, for some
  designs, a limit cycle just above the smallest normal value, because FTZ
  flushes a2·y on its own while |a1| > 1 lifts the state again (the floating
  biquad cascade's honest limit above). It grows about as 1/fc as the cutoff
  falls: from full scale on x86-64, 124 times `FLT_MIN` at order 8,
  fc = 0.002, 243 times at 0.001 and 2,351 times at float's bound; 99 times
  `DBL_MIN` at order 4, fc = 0.002, and 12,474 times at double's bound
  (in float, 2.8e-35 of full scale, about −690 dB); none for most designs between
  0.01 and 0.3. Pinned below 256 times at the battery's cutoffs, 0.005 to
  0.45, orders 4 and 8, either FTZ setting. So in a chain holding the stage,
  `chain::flush()` equals zero padding up to that residue rather than bit
  for bit (an all-FIR chain keeps the bit-for-bit equality), and what follows
  it is the residue, not exact silence.
- **Subnormals are the host's** (I9, `subnormals.h`): the kernel never
  writes the control register and reads it as a debug-build precondition.
  Non-finite input poisons the state until `reset()`.
- **Formats:** float keeps its states, coefficients and products in float
  (I14). Against double, the largest difference over 48,000 frames of white
  noise at 0.9 of full scale (x86-64): order 4, 2.0e-7 at fc = 0.25, 4.5e-7
  at 0.07, 9.5e-6 at 0.01, 1.1e-4 at 0.001 and 2.4e-3 at float's bound;
  order 8, 2.6e-7, 8.9e-7, 1.2e-5, 8.3e-5 and 3.4e-3. Transposed DF II's
  float error grows as the poles approach z = 1: below about 0.001 the
  float stage is stable but far from the double one, and double is the
  profile for low cutoffs. Q15 / Q31: the kernel has them (direct form I);
  this stage does not offer them yet (a low cutoff needs section scaling).

Pinned by `tests/test_butterworth.cpp`: the magnitude at DC, the cutoff and
Nyquist at orders 4 and 8 over six cutoffs; the sections ordered by pole
radius and equal to the cookbook biquads at the exact Qs within 1e-13 and at
the printed ones within 1e-8; float and double against a double reference
recursion; chunking and channel independence bit for bit; the group delay
against the phase's numerical derivative, the nominal and the passband
deviation; the tail after `window_frames()` below 256 times the smallest
normal value, from full scale and from 10⁻⁶ of the format's largest value,
under the host's FTZ; the window values, and the window at the bounds
within a 32-bit `size_t`; float's bound (every pole inside from it, the
margin gone at 3e-5); the limit cycle FTZ leaves an order-8 design in (the
honest limit); `set_cutoff()` as a fresh construction bit for bit (design,
window, latency members and, after `reset()`, the output) at orders 4 and
8 across an oversampling kernel's cutoffs and the bounds, leaving the state
alone, and asserting outside the range in a Debug build; construction,
`set_cutoff()`, `process()` and `reset()` allocating nothing (the heap meter,
`tap/dsp/analysis/heap_meter.h`); two stereo stages chained as the two in
sequence for any chunking; the nominal latency, deviation and group delay
composed with a decimator; flush equal to zero padding up to the recursive
residue over a chain of the stage and a decimator; the control register never
written; the rejections. `tests/test_iir_block.cpp` pins the block kernel
against the per-sample one, output and state bit for bit: 1 to 8 sections,
1, 2, 3, 6 and 12 channels, empty, single, odd and uneven calls, in place,
the two kernels alternating on one state, a first-order row, a NaN. Without FTZ,
`tests/test_butterworth_subnormals.cpp` pins the tail below the smallest
normal value (in `tap_dsp_iir_subnormals`, built with
`TAP_DSP_ALLOW_SUBNORMALS`). `test_butterworth.cpp` is its own executable,
`tap_dsp_chain_stages`, built without fp contraction (fft.h's policy for
every bit-exact claim; `tests/CMakeLists.txt` says why) and carrying the heap
meter's replacement operators. Every QEMU leg runs both whole.

## The FIR substrate

Five headers carried from **SampleRateTap** (where they design and run the
ASRC's polyphase datapath) and promoted here so **RatioTap**'s fixed-ratio
44.1↔48 converter — now the SampleRateTap family's `bridge` engine, beside the
ASRC as `async` — and any future FIR consumer — shares one implementation,
plus two landed here first for the family's third engine, `rational`
(`nyquist.h`, the L-th-band designer, and `chain.h`, the stage composition
below every engine), plus the FFT's butterfly arithmetic trait
(`fft/fft_arith.h`), which is built over the same sample formats and
documented here until the Stage 3b README rewrite moves it to the FFT
section's profiles table.
The performance-sensitive pieces are regression-gated in the SampleRateTap
family's instruction-count CI (both engines, Cortex-M33/M55, Hexagon, ±3%); treat measured claims in
the header comments as contracts.

### `tap/dsp/kaiser.h` — FIR prototype design

Kaiser-windowed sinc prototype design for L-phase polyphase banks: `bessel_i0`,
`kaiser_beta` (Kaiser's empirical fit), `estimate_taps` (the harris length
estimate), `design_prototype`, and `design_prototype_compensated` — the
zeros-at-k·fs variant with passband droop pre-compensated (closed-form, no
FFT), which turns branch-DC uniformity into exact transmission zeros at every
multiple of the sample rate. Runtime design in double, deliberately not
constexpr (the header's design note does the arithmetic); run it in a
constructor, off the audio path. The note's scope, since 0.3.0: the FIR
designers (here and `nyquist.h`) are runtime double; the recursive designers
of `iir_design.h` are constexpr, because a design there is a handful of
coefficients (at most 16 per half-band, 5 per Chebyshev II section), a few
thousand constant-evaluation steps, not a table of thousands of taps. Also exports `solve_dense`, the small dense
solver the compensated design and the analysis instruments share (noexcept and
allocation-free: it pivots by exchanging rows in the caller's buffers).

### `tap/dsp/nyquist.h` — L-th-band (Nyquist) FIR design

The design math of the SampleRateTap family's `rational` engine (its
`PLAN.md` 2.2), landed here first: `design_nyquist(h, L, beta)` writes the
Kaiser-windowed sinc with cutoff exactly π/L at length N = 2mL − 1
(`nyquist_length`, `nyquist_centre`, `nyquist_nonzero_taps`,
`is_nyquist_length`), the centre tap exactly 1.0 and every L-th tap from it
exactly 0.0 (written as zeros, not libm's 1e−16), every polyphase branch
normalized to DC gain 1 so the whole sums to L — an interpolator's branch 0
is a copy, a decimator uses h / L — and the two halves bit-identical (the
mirror branches' gains can differ by an ulp, so the second half is copied
from the first: a mirrored pair that lands in one quantized row then ties
exactly on every host, and `quantize.h` breaks the tie by index, the same
way everywhere). The response is antisymmetric about the
lower rate's Nyquist, so the transition is symmetric (f_p + f_s = r: the
engine's coverage rule with equality) and a half-band stage computes 2m + 1
MACs per output instead of 4m − 1. `nyquist_response_db`,
`nyquist_worst_stopband_db` and `search_nyquist_m` (the smallest m meeting
a stopband with ≥ 1 dB margin on a grid of `grid_points`, the `bridge`
criterion; the default 1024 points resolve designs up to a few hundred
taps and a long design asks for more — the 8th-band 120 dB candidate reads
−121.7 dB on 1024 points and −119.3 dB on 8192, a pinned test) are the
design-time instruments. Published literature: Mintzer 1982,
Vaidyanathan §4.6.

Pinned by `tests/test_nyquist.cpp`: the exact centre and zeros for L ∈
{2, 3, 4, 6, 8}, per-branch unity, the sum of the L shifted zero-phase
responses equal to 1 at every frequency, the half-band MAC count, and the
searched minima measured for the engine's profile candidates (economy
half-band m = 11, N = 43; economy third-band m = 11, N = 65; transparent
half-band m = 31, N = 123, each with the next-shorter length missing the
spec).

### `tap/dsp/chain.h` — the synchronous-stage concept and `chain<>`

The composition point below the family's engines (rational `PLAN.md`, R6):
the `sync_stage` concept (a `sample` type, the ratio `k_up` / `k_down` and
the channel count `channels()` as compile-time numbers, `process`,
`outputs_for`, `reset`) and `chain<Stages...>`, which requires its stages to
run one channel count (a `static_assert` names the rule), runs them in order
through scratch sized at construction for that count (noexcept,
allocation-free, bit-identical for any chunking),
composes `outputs_for` forward and `frames_needed` (the smallest n with
`outputs_for(n) >= k`; not an equality, an interpolating stage completes L
outputs per input) backward, reports the chain's group delay as an `exact_ratio` at its
output rate from each stage's own (`latency_output_frames()` or
`latency_input_samples()`), and drains the stream with `flush()`: each
stage in order is fed its `window_frames()` of zeros (the zero input after
which its output is silence, `taps - 1` for a decimator) through the stages
after it, `flush_output_frames()` frames in all, bit-identical to zero
padding the chain's input. Since 0.5.0 a chain reports `window_frames()`
itself when every stage does: a bound valid from any state, the last
stage's window and before it each stage's own plus the inputs it needs,
at most, to write the window after it (`ceil(k M / L) + M` for a stage on
its ratio's schedule), the bound `subband.h`'s flush already composed. So
a chain that holds a chain flushes: the inner chain is fed that many zeros
and writes the flattened chain's frames and possibly a few more, every one
of them silence. A chain is written by the caller as a type,
never looked up from a rate pair. The count is never a run-time value: a
stage declares it as a `static constexpr` `channels()` of at least 1, and no
chain constructor takes one (since 0.2.0; 0.1.0's chain took the count at run
time to size its scratch). `basic_decimator` satisfies the concept (`sample`,
`k_up`, `k_down`, `window_frames` and `channels()`, one, added), which is how
the helper is proven on an existing primitive before any new engine depends
on it.

Pinned by `tests/test_chain.cpp`: a chain of the by-2 and by-3 decimators
equals the two run in sequence bit for bit in every sample format, for
chunks of 1 to the whole stream; `outputs_for` exact from every phase;
`frames_needed` the smallest n that reaches k; the stages reset at
construction and read-only through the chain; latency 80/3 output frames (40/2 of the
by-2's outputs through the by-3, plus 60/3), the impulse peak where it
says; flush equal to zero padding bit for bit from every phase in every
format, `flush_output_frames` exact, the chain silent afterwards; a
chain's `window_frames()` composed backward, silencing it from every phase;
a chain holding a chain (first, last and two deep) processing as the
flattened chain and flushing as zero padding in every format, the flattened
chain's flush its prefix and silence after it; reset
bit-exact; the concept rejecting a stage with no channel count, a run-time
one or zero; a two-channel chain keeping its frames interleaved through its
scratch across chunks, and a chain of a stage without a count not a type.
`tests/compile_fail/chain_channels.cpp` pins that a chain of stages running
different counts fails on the chain's `static_assert`, as the build's only
error.

**The nominal latency (since 0.3.0; SampleRateTap's `IIR_PLAN` I5).** A
recursive stage has no exact latency, so a chain holding one reports, in
doubles at its output rate, `nominal_latency_output_frames()` (each stage's
nominal — an IIR stage's DC group delay, a near-linear stage's delay path,
an exact stage's exact latency — converted and summed),
`latency_deviation_output_frames()` (the stages' passband deviations summed,
a bound) and `group_delay_output_frames(f)` (f in cycles per output frame;
each stage's delay at the same frequency in its own frames, an exact stage's
constant). A stage reports the three itself (`reports_nominal_latency`) or an
exact latency; a chain reports them too, so chains nest. The exact
`latency_output_frames()` still compiles only for all-exact chains, the
`sync_stage` concept is unchanged, and `subband.h` reads a near-linear
stage's `delay_path_output_frames()` as before. Pinned by
`tests/test_chain_iir.cpp` on test-local stages over the recursive kernels:
the members against each stage's and against the delay measured from the
chain's impulse response (within 0.01 frames), an all-exact chain's nominal
its exact latency (80/3) with deviation 0, a delay-path chain's nominal its
exact 17/2, and a 1:1 stage in a chain (the guard band filter's shape) in
every format: the ratio, `outputs_for` from every phase, `frames_needed`,
chunking bit for bit, and flush equal to zero padding.
`tests/compile_fail/chain_iir_exact_latency.cpp` pins that a chain with a
recursive stage has no exact latency, the rule's `static_assert` the only
error.

### `tap/dsp/sample_traits.h` — sample formats: double, float, Q15, Q31

The family's sample-format substrate: how each sample type stores
coefficients, accumulates dot products, and rounds/saturates back to samples.
`double` is a sample format because it is the golden model of every
primitive: a traits-based primitive instantiates its reference profile
through the same substrate as its embedded profiles, and the cross-precision
pins measure float/Q15/Q31 against it.

| Type | Coefficients | Accumulation | Output |
|---|---|---|---|
| `double` | double | double | identity (the golden model) |
| `float` | float | double | plain cast |
| `std::int16_t` | Q1.14 | int64, exact | single Q29→Q15 round-half-up, saturating |
| `std::int32_t` | Q1.30 | int64, products pre-shifted to Q45 | single Q45→Q31 round-half-up, saturating |

The Q ladder is spelled as named constants on each fixed-point
specialization (`k_sample_frac_bits`, `k_coeff_frac_bits`,
`k_accum_pre_shift`) with the accumulator format and the single 14-bit
finalize shift derived from them and `static_assert`ed; `k_coeff_scale` is
`2^k_coeff_frac_bits`, and `k_is_fixed_point` is what selects the
fixed-point algorithm in `quantize.h`. Every member is `constexpr`, and the
`sample_type` concept requires a value-initialized accumulator to be the
additive identity.

**Fixed point is a first-class embedded direction, not a legacy path.** The
Q15/Q31 profiles exist for targets where double (sometimes any float) is
unaffordable — SampleRateTap measured its float datapath at ~19× the
instruction count of Q15 on a Cortex-M33 (soft-double accumulation). Expected
deployments include Bluetooth-adjacent conversion (the `bridge` engine) and M33/M55-class
eurorack and pedal targets running TapTools primitives. Per-primitive adoption
is opt-in, and each adoption is its own documented Q-format design: the ladder
of headroom bits, pre-shifts, and the single rounding point is a per-datapath
decision. That is also why these are *traits over raw sample types* rather
than `q15`/`q31` wrapper classes — the arithmetic contract stays visible at
the use site and pinnable by tests, buffers arrive from codecs and C ABIs as
plain `int16_t`/`int32_t`, and the SMLALD kernel's paired loads stay legal.

This header is the format core only. Engine-specific extensions (e.g.
SampleRateTap's inter-phase coefficient blending) derive from these
specializations and refine the `tap::dsp::sample_type` concept.

`finalize_divided<S, D>(acc)` is `finalize()` of `acc / D` under the same
single rounding, for a dot whose coefficients carry a gain of D the format
cannot hold: a Q15 polyphase decimator by M quantizes each of its M branches
at its own unity sum (full Q1.14 precision, where `h / M` in one row would
lose about 20 log10 M dB of stopband) and divides the summed branches by M
here. A power-of-two D widens finalize's shift (exact; `D = 1` *is*
`finalize`); otherwise Q31 takes the exact round-half-up quotient and Q15 a
multiply-back by `round(2^26 / D)`, exact at every multiple of `D·2^14` (all
65 536 DC values pinned for D = 1 … 8) and the exact quotient except within
`D·2^-11` of a rounding boundary. `decimate.h` adopts it for its Q15 tables
(`M·h` in one row), and SampleRateTap's rational engine for its Q15
decimators (each branch at unity).

### `tap/dsp/fft/fft_arith.h` — butterfly arithmetic for the FFT profiles

The sibling trait the fixed-point real FFT is written against (its docstrings
are that kernel's specification — shift-before-butterfly, the magnitude bound,
the rounding count per complex product, the BFP headroom rule, the Q31 input
pre-shift, the twiddle generator; `tests/test_fft_arith.cpp` pins every
number). The profiles table in the FFT section carries the resulting
numbers. One int32 kernel serves both fixed profiles: `fft_arith<std::int32_t>`
carries `mul_coeff` (int32 × Q1.30 → int64, `>> 30` with one round-half-up,
saturating), saturating `add` / `sub`, `shr_round` (round-half-up), and
`headroom_bits` (the block's shared redundant sign bits, via
`std::countl_zero`); `fft_arith<std::int16_t>` is the Q15 I/O width —
`widen` (`<< 14`, two guard bits) and `narrow` (round-half-up, saturating) —
and names the int32 trait as its `work`. Twiddles are
`sample_traits<std::int32_t>::coeff` (Q1.30) for both fixed profiles, and 1.0
is representable. The `float` / `double` specializations are the same names
over plain arithmetic. Everything is `constexpr` and `noexcept`. The kernel
written against it is `fft/fixed_point.h`; its tables (bit reversal, Q1.30
twiddles, the real post-pass coefficients, each one rounding from its double)
are the free functions in `fft/tables.h`, pinned by checksum per certified N.

### `tap/dsp/fir_kernels.h` — dot-product kernels

The FIR hot loops, target-gated the way SampleRateTap's optimization campaign
measured them: `dot_row` (planar; routes Q15 through an eight-lane Helium
VMLALDAVA reduction with a predicated tail on MVE cores — the Cortex-M55 /
M85 class, where GCC 13 leaves the scalar loop at one SMLALBB per tap — and
through a dual-MAC SMLALD loop on DSP-extension Arm cores without Helium; both
bit-exact by construction, since every 16 × 16 product is exact and the int64
sum associative, and both pinned at every tap count from 0 to 40), and the
channel-parallel pair `dot_tile_frame_major` / `dot_rows_frame_major`
(register-blocked 8/4/2/1 tiles over frame-major storage, coefficient
broadcast across channel lanes — bit-exact against the planar path for every
sample type, float included, because lanes are channels, not taps), and
`accumulate_row`, `dot_row`'s accumulation without its finalize, so a caller
can sum several rows under one rounding point (a polyphase decimator that
dots only its nonzero branches; `finalize(accumulate_row(accum{}, …))` is
`dot_row` bit for bit, pinned). The
`TAP_DSP_CHANNEL_PARALLEL` / `TAP_DSP_CP_MIN_CHANNELS` gates encode which
targets prefer which layout.

### `tap/dsp/quantize.h` — row-sum-preserving quantization

`quantize_row_preserving_sum`: quantizes one polyphase branch to a fixed-point
coefficient format while preserving the row's DC sum *exactly*
(largest-remainder distribution of the rounding residual — "the coefficients
of every phase must add to one", R. Bristow-Johnson, music-dsp). Selected by
the trait's `k_is_fixed_point`: plain conversion for double and float. A tap
at the format's rail is never wrapped: each step goes to the largest-remainder
tap that can still move in that direction, so the sum is preserved whenever
one exists. Design-time code.

`quantize_symmetric_row_preserving_sum` (since 0.4.0) is its form for an
even-length symmetric (type II) row: the first half through
`quantize_row_preserving_sum`, then mirrored, in any of the four formats. The
quantized row is symmetric bit for bit, so its zero at Nyquist is exact in
the quantized taps (in integers for Q15 and Q31), and it sums to exactly
twice the half's preserved sum, within one coefficient LSB of the design's
(the half's promise, half an LSB, doubled). Quantizing the whole row instead
can split a mirrored pair where the sum fix-up ties, and lose both.
`tests/test_quantize.cpp` pins the contract on every committed guard band row
and a Kaiser design, and its bit-identity with the half-then-mirror procedure
SampleRateTap's guard band filter ran inline. The guard band filter's Q15
rows (`iir_tables.h`) are designed and measured through it.

### `tap/dsp/analysis/` — measurement instruments

The quality-measurement harness the converter suites share: `sine_analysis.h`
(least-squares single-tone fit, frequency-tracked variant, `snr_db`) and
`multitone_analysis.h` (pink log-spaced `tone_comb`, joint least-squares
multitone fit, `program_weighted_snr_db` — the program-weighted metric with
Fisher-weighted ratio pooling). Instrument floors on exact synthetic signals
are pinned by `tests/test_analysis.cpp`, so a consumer's quality gate never
silently rests on a degraded instrument.

`heap_meter.h` is the footprint instrument: a `heap_meter` scope reports the
bytes allocated through the global `operator new` while it is open — still
live, at peak, in total — and the allocation count. The counting replacements
of every form of the global `operator new` and `operator delete` are emitted
only in the one translation unit that defines
`TAP_DSP_HEAP_METER_REPLACE_GLOBAL_NEW` before including the header, so a
footprint test lives in an executable of its own. The counts are the requested
bytes, not the allocator's, so a footprint formula can be pinned against them
exactly — written per allocation through `heap_meter::allocator_request(n)`,
which models what `std::allocator` asks for (MSVC's on x86 and x64 pads a
request of 4096 bytes or more by `sizeof(void*) + 31`, to align it to 32
bytes). Not under ASan, TSan or MSan, whose runtimes define `operator new`
themselves. Pinned by `tests/test_heap_meter.cpp` (its own target), on every host
and QEMU leg.

## `tap/dsp/math.h` — periodic Hann and decibel helpers

The three scalar formulas the primitives build on, public so consumers call
the same functions instead of re-typing the expressions:
`periodic_hann(i, n)` = `0.5 - 0.5 * std::cos(2.0 * std::numbers::pi * double(i) / double(n))`
(the DFT-even window, denominator `n`, not the symmetric `n - 1` form),
`power_db(r)` = `10.0 * std::log10(r)` and `amplitude_db(r)` =
`20.0 * std::log10(r)`, all in double. The expressions and their association
order are contract points: pvoc and log_mel build their windows from
`periodic_hann` and the analysis instruments report through `power_db`, so a
consumer that swaps its own copy of one of these exact expressions for the
helper moves no bit (`tests/test_math.cpp` pins helper against literal over a
sweep of `i`/`n` and ratios, on every host and QEMU leg). There is no exported
pi: use `std::numbers::pi`, which is the double the helpers use.
`tap::dsp::detail` re-exports the same three functions by using-declaration for
the headers' own call sites.

## `tap/dsp/detail/expects.h` — the two contract idioms

`TAP_EXPECTS(cond)` is the house precondition (STYLE.md §4): `assert` under a
name that says contract — a check in a debug build, nothing in a release
build, never a throw — for what the caller is bound to meet, with the
release-mode predicate beside it (`fft.h`'s power-of-two size has
`supports_size`). `tap::dsp::raise<E>(args...)` is the construction-time
rejection of what the caller could not check without the library's own
arithmetic: `throw E(args...)` where the translation unit has exceptions
(`TAP_DSP_HAS_EXCEPTIONS` is 1), `std::terminate()` where it does not
(`-fno-exceptions`, the default of the Pico SDK, Zephyr and most Cortex-M
firmware trees), so every header compiles either way and a rejected
configuration is an exception on a host and a contract violation on such a
firmware — what the standard library's own allocation failure already is
there. The predicate beside the rejection (a `validate()` that returns the
reason) is what such a firmware checks first. `accelerate.h`'s refused vDSP
table is the one site in this tree; the SampleRateTap engines' constructors,
which validate by throwing today, are the consumers it was written for.
Pinned by `tests/test_expects.cpp` (hosted and on the QEMU legs: `raise`
throws exactly `E` with the caller's text and is `[[noreturn]]`; `TAP_EXPECTS`
does not evaluate its argument in a release build) and by
`expects.PublicHeadersCompileWithoutExceptions`, which compiles every public
header under `-fno-exceptions` on the hosted GCC and Clang legs.

The whole battery runs in that configuration too, with RTTI off as well: a
consumer that forbids both builds this way. `-DTAP_DSP_NO_EXCEPTIONS_RTTI=ON`
compiles this project's own targets with `-fno-exceptions -fno-rtti`
(GoogleTest and the battery, the benches and the C ABI; GCC and Clang only),
and CI runs it on Linux and on the M33 and M55 legs, the srdif fingerprints
included, so the output there is the same bits as everywhere else. In that
build `raise` terminates (pinned by death tests on the host), and the C ABI's
create functions terminate on a refused allocation rather than returning
NULL. No header uses RTTI; the one test that reads it, `TypeNamesCarryTheTag`,
has a compile-time counterpart that holds in every build.

## Notebooks

The notebooks drive the **actual shipping C++** through the C ABI in
`tools/capi/` (ctypes bridge: `notebooks/dsptap_py.py`, which builds
`build_capi/` on first import; the bridge also exposes `LogMel` and
`Decimator` for MuTap's keyword-spotting notebooks). They are committed
executed; re-execute with
`jupyter nbconvert --to notebook --execute --inplace notebooks/<name>.ipynb`
when a primitive's behavior changes.

`notebooks/pitchshift.ipynb` measures the three pitch primitives. It
documents the two findings from the primitives' development: PSOLA's
envelope-resampling nature (why it preserves formants *and* why a pure tone
shifted an octave thins out), and the measured level collapse of naive
phase-vocoder bin remapping vs the shipping peak-locked design — plus the LPC
formant-preservation demo.

`notebooks/fft.ipynb` measures the real FFT through `dsptap_py.RealFFT`
(profiles `"double"` and `"float"`; `unpack`/`pack` convert the header's
packed exp(+i) layout to and from `numpy.fft.rfft`'s convention): the packing
and sign contract against numpy, the float32 profile's per-bin error against
the double golden model vs N — on the profile's own arithmetic via the raw
in-place entry points, not a double round trip — and the round-trip error of
both profiles. With tap/DspTap#30 (Stage 3c) `notebooks/fft.ipynb` also
measures the four fixed-point configurations (`"q15"`, `"q31"`, `"q15_bfp"`,
`"q31_bfp"` in `dsptap_py.RealFFT`, the exponent returned with every
transform): the noise-floor table of the fixed-point design record in
`docs/fft-design.md` ("The fixed-point profiles", §5) re-measured through the
ABI beside the battery's pins, the DC-path bias to N = 65536, the
block-floating exponent on speech-like material, and the battery's output
fingerprints reproduced through the ABI; `sine_analysis.h` /
`multitone_analysis.h` score Q15 / Q31 spans directly from #30 onward.

## Versioning

`tap/dsp/version.h` states this tree's version as `TAP_DSP_VERSION_MAJOR`,
`_MINOR` and `_PATCH`, equal to the CMake `project()` version
(`Version.HeaderMatchesTheProjectVersion` pins the two together), and as one
integer, `TAP_DSP_VERSION`, packed `(major << 16) | (minor << 8) | patch` like
the Tap family's C ABI version probes, so it compares in `#if` and
`static_assert` alike. A tag `vX.Y.Z` marks the commit that carries each
version. DspTap is tagged on its own cadence, never in lockstep with a
consumer: a consumer pins a commit, names the tag in the commit that moves
the pin, and asserts the oldest DspTap it was tested with:

```cpp
#include "tap/dsp/version.h"
static_assert(TAP_DSP_VERSION >= TAP_DSP_VERSION_ENCODE(0, 1, 0),
              "this library needs DspTap 0.1.0 or newer");
```

Before 1.0 a minor version may break a consumer, so a floor states what was
tested, not what will keep working. 0.2.0 is such a break: `chain.h`'s
stages declare their channel count at compile time and the chain's
run-time-count constructor is gone (SampleRateTap's A13). 0.3.0 adds the
sub-band round trip (`subband.h`) and the recursive substrate
(`iir_design.h`, `iir_kernels.h`, `subnormals.h`, `chain.h`'s nominal
latency), and changes nothing an 0.2.0 consumer calls. 0.4.0 adds the linear
guard band filter's Q15 rows (`k_guard_band_fir_q15`) and the symmetric
quantizer they are measured through (`quantize_symmetric_row_preserving_sum`),
and changes nothing an 0.3.0 consumer calls: every committed row before it is
unchanged, bit for bit. 0.5.0 lets a chain report `window_frames()`, so a
chain that holds a chain flushes, and changes nothing an 0.4.0 consumer
calls: every chain that flushed before writes the same frames, and the
sub-band round trip's window is the same number. 0.6.0 adds the Butterworth
low-pass stage (`butterworth.h`: `basic_butterworth_lowpass`, redesigned in
place by its allocation-free `set_cutoff()`, and `design_butterworth_lowpass`)
and changes nothing an 0.5.0 consumer calls: no header an 0.5.0 consumer
includes changes but `version.h`. 0.7.0 makes that stage cheaper to run and to
redesign: `process()` runs the new block kernel
(`iir_kernels.h`'s `biquad_cascade_block`), `set_cutoff()` and the
constructor no longer compute the latency deviation, which
`latency_deviation_output_frames()` now computes on each read (the same
number; `butterworth_latency_deviation` is the scan as a free function), and
the flush window is found with fewer matrix products (the same number at
`fft.h`'s fp-contraction policy). Every signature an 0.6.0 consumer calls is
unchanged. Two things such a consumer can observe: a read of the deviation
now costs a 1,025-point scan in double (a design-time query: keep the
number if you need it on the audio thread), and with fp contraction on, a
stage's last bits may depend on how its input is chunked (with contraction
off they cannot, as before).

## Build

Standalone (builds the srdif engine, plus vDSP on macOS, and runs the
tests):

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
```

CI also runs an emulation-sized selection of the battery
(`tests/bare_metal_main.cpp`) bare-metal under QEMU on four Cortex-M legs:
`cortex-m4-softfp` and `cortex-m4f` (`cmake/arm-cortex-m4-mps2.cmake`, the
FPU flavour selected by `-DTAP_DSP_M4_FPU=ON`), `cortex-m33`
(`cmake/arm-cortex-m33-mps2.cmake`) and `cortex-m55`
(`cmake/arm-cortex-m55-mps3.cmake`, where the CMSIS-DSP Helium FFT backend is
ON — detected, not pinned: the compiler targets MVE-F — and its parity suite
runs against the srdif engine; CI asserts the detected value on every
leg). Every suite compiled into a test
executable runs on the target unless excluded by name in
`tests/CMakeLists.txt` (a negative filter; each exclusion is a budget note).
To run one locally, with `arm-none-eabi-g++` and `qemu-system-arm` on `PATH`:

```sh
cmake -S . -B build-m33 -DCMAKE_BUILD_TYPE=MinSizeRel \
      -DCMAKE_TOOLCHAIN_FILE=cmake/arm-cortex-m33-mps2.cmake
cmake --build build-m33
ctest --test-dir build-m33 --output-on-failure
```

Benchmarks and the per-target instruction-count ratchet build with
`-DTAP_DSP_BUILD_BENCH=ON`; policy and workflow in [`bench/README.md`](bench/README.md).

### As a submodule

```cmake
add_subdirectory(submodules/dsptap)   # or however it is pinned
target_link_libraries(my_dsp PRIVATE tap::dsp)
```

`tap::dsp` is a pure INTERFACE target (the headers; under `TAP_DSP_FFT_CMSIS`
alone it also links the `tap::dsp_fft` static library that carries the CMSIS
objects, and nothing else is ever compiled); it does not build the tests when
added as a subdirectory (`TAP_DSP_BUILD_TESTS` defaults OFF unless
top-level). The per-platform float32 backend defaults follow the target: vDSP
on Apple; CMSIS where the compiler targets floating-point Helium — a compile
check on `__ARM_FEATURE_MVE & 2` under your toolchain's flags, so it does not
matter how the toolchain spells the CPU (`-mcpu=cortex-m55`, `cortex-m85`,
`-march=armv8.1-m.main+mve.fp`, hard or softfp float ABI); the srdif
engine elsewhere, including Cortex-M0+/M4/M7/M33 and bare-metal Cortex-A/R
toolchains. The check sees the C++ compiler with `CMAKE_CXX_FLAGS` only (keep
the CPU flags identical in `CMAKE_C_FLAGS`, which the CMSIS objects use) and
re-runs on every configure; flags it cannot see (per-config or per-target
options) leave it OFF, which costs speed, never the contract. Override with `-DTAP_DSP_FFT_CMSIS=ON|OFF`,
`-DTAP_DSP_FFT_ACCELERATE=OFF` etc.

## Provenance

This code was carried, with textually identical vendored Ooura sources, inside
both **MuTap** (adaptive filtering) and **AmbiTap** (ambisonics / binaural
convolution). The C++ wrappers had begun to diverge — MuTap grew the templated
`basic_real_fft<Sample>` and the CMSIS/vDSP backends; AmbiTap kept an older
double-engine wrapper with no backends — so a bug fix or a new backend in one
would silently miss the other. DspTap is the consolidation: one wrapper, one
contract, one home for the next backend. The unified wrapper is MuTap's
backend-capable `basic_real_fft`, generalized to the `tap::dsp` namespace.

The FIR substrate is the second consolidation wave, moved here from
**SampleRateTap** (its `srt/detail/kaiser.h`, `srt/sample_traits.h` format
core, the dot kernels from `srt/polyphase_filter.h`, the row-sum quantization
from its bank constructor, and the `tests/support/` measurement harness) at
the moment **RatioTap** became the second consumer — the same
extract-on-second-consumer rule that created this repo.

See [`third_party/cmsis-dsp/VENDOR.md`](third_party/cmsis-dsp/VENDOR.md) for
the vendored CMSIS-DSP subset's provenance, and [`NOTICE.md`](NOTICE.md) for
the licenses. The FFT the two libraries carried was Takuya Ooura's `fftsg.c`,
vendored. At Stage 2a (#28) DspTap landed a C++20 port of its `rdft`
(`include/tap/dsp/fft/split_radix.h`), bit-identical to the C, and routed the
floating profiles at it at Stage 2b (#31); at Stage 2c (#32) the vendored C
left the shipping tree, and a reference copy (`fftsg.c` and its float
instantiation `fftsg_float.c`) under `tests/reference/ooura/` served the
bit-identity gate alone until both MuTap and MuTap-Max pinned a tree
containing 2c; Decision D6 (#36) then deleted it and pinned the port's
output as fingerprints. The Q15 / Q31 profiles landed at Stage 3b (#27); their
real post-pass, a transcription of the package's, was re-derived from the
literature under a clean-room procedure in #39 (bit-identical). In #42 the
maintainer replaced the port itself: the floating profiles run the srdif
engine (`include/tap/dsp/fft/srdif.h`), written clean-room from the
literature and DspTap's own fixed-point engine, and every floating output
bit changed once. Since #42 DspTap ships, in the maintainer's judgement
(`NOTICE.md`), no code derived from Ooura's package;
`third_party/ooura/readme.txt` and `LICENSES/LicenseRef-Ooura.txt` stay as
the license record for the older trees consumers pinned (from #28 to #42's
base they carry the port, earlier ones the C itself; `NOTICE.md`).
The procedures, and what each did and did not cover, are in the design note
("Provenance and licensing") and `NOTICE.md`. The
design note is [`docs/fft-design.md`](docs/fft-design.md), filled in as each
stage lands; the plan of record is `docs/audit-fft-and-code-smells.md` (#25).

## License

DspTap's own code is MIT (`LICENSE`). Vendored third-party code keeps its own
license. Since tap/DspTap#42 DspTap ships, in the maintainer's judgement
(`NOTICE.md`), no code derived from Ooura's FFT package. Earlier trees did:
from Stage 2b (#31) until #42 the floating profiles ran a C++ port of its
`fftsg.c` (SPDX `LicenseRef-Ooura AND MIT`), a derivative work whose
redistribution relied on the package's modification grant — its terms, in
`third_party/ooura/readme.txt`, are "You may use, copy, modify this code for
any purpose and without fee. You may distribute this ORIGINAL package." The
readme and `LICENSES/LicenseRef-Ooura.txt` stay permanently as the license
record for those trees. CMSIS-DSP / CMSIS-Core are Apache-2.0 with SPDX
headers retained in every file. The canonical statement, and the maintainer's
reading of what it covers — a judgement call, not legal advice — is
[`NOTICE.md`](NOTICE.md).
