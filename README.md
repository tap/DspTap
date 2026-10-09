# DspTap

Shared DSP primitives for the **Tap** family of audio libraries. Header-only,
plain portable C++ (C++20, standard library only), no Max/Min or framework
dependency — consumed as a git submodule by the individual libraries.

Today it holds seven primitives, plus the [FIR substrate](#the-fir-substrate) —
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
constructor, off the audio path. Also exports `solve_dense`, the small dense
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
the `sync_stage` concept (a `sample` type, the ratio `k_up` / `k_down` as
compile-time numbers, `process`, `outputs_for`, `reset`) and
`chain<Stages...>`, which runs the stages in order through scratch sized at
construction (noexcept, allocation-free, bit-identical for any chunking),
composes `outputs_for` forward and `frames_needed` (the smallest n with
`outputs_for(n) >= k`; not an equality, an interpolating stage completes L
outputs per input) backward, reports the chain's group delay as an `exact_ratio` at its
output rate from each stage's own (`latency_output_frames()` or
`latency_input_samples()`), and drains the stream with `flush()`: each
stage in order is fed its `window_frames()` of zeros (the zero input after
which its output is silence, `taps - 1` for a decimator) through the stages
after it, `flush_output_frames()` frames in all, bit-identical to zero
padding the chain's input. A chain is written by the caller as a type,
never looked up from a rate pair. `basic_decimator` satisfies the concept
(`sample`, `k_up`, `k_down`, `window_frames` added), which is how the helper
is proven on an existing primitive before any new engine depends on it.

Pinned by `tests/test_chain.cpp`: a chain of the by-2 and by-3 decimators
equals the two run in sequence bit for bit in every sample format, for
chunks of 1 to the whole stream; `outputs_for` exact from every phase;
`frames_needed` the smallest n that reaches k; the stages reset at
construction and read-only through the chain; latency 80/3 output frames (40/2 of the
by-2's outputs through the by-3, plus 60/3), the impulse peak where it
says; flush equal to zero padding bit for bit from every phase in every
format, `flush_output_frames` exact, the chain silent afterwards; reset
bit-exact.

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

### `tap/dsp/analysis/` — measurement instruments

The quality-measurement harness the converter suites share: `sine_analysis.h`
(least-squares single-tone fit, frequency-tracked variant, `snr_db`) and
`multitone_analysis.h` (pink log-spaced `tone_comb`, joint least-squares
multitone fit, `program_weighted_snr_db` — the program-weighted metric with
Fisher-weighted ratio pooling). Instrument floors on exact synthetic signals
are pinned by `tests/test_analysis.cpp`, so a consumer's quality gate never
silently rests on a degraded instrument.

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
