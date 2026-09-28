# The real FFT: design note

*September 2026. This is the design note that is too long for a header — the
role `kaiser.h`'s design note plays for the FIR substrate. It is filled in
stage by stage; every entry that still depends on a measurement is marked
`TODO(stage …)` and names the stage whose PR supplies the number. Landed so
far: Stage 2a (the port beside the C, tap/DspTap#28), Stage 3b (the Q15 / Q31
profiles, tap/DspTap#27; its design record, formerly
`docs/fft-fixed-point.md`, is folded in below as "The fixed-point profiles"),
Stage 2b (the routing flip, tap/DspTap#31), Stage 2c (the C leaves the
shipping tree, tap/DspTap#32), Stage 4 (the engine parameter and the ABI tag,
tap/DspTap#35; "Stage 4" below), and the replacement of the floating engine
(tap/DspTap#42: the port of Stage 2a gave way to `fft/srdif.h`, written
clean-room; "The floating engine (srdif)"). Sections written before #42
describe the port where they describe the floating engine at all, and say so;
the contract table and the fp-contraction policy state the tree after #42.
Until every stage has landed, the plan of record is
[`audit-fft-and-code-smells.md`](audit-fft-and-code-smells.md) (parts are
cited as "audit Part N" below), and each PR links its stage.*

## Purpose and scope

`tap::dsp::basic_real_fft<Sample, Policy =
detail::default_real_fft_policy_t<Sample>>` (`include/tap/dsp/fft.h`) is a
four-profile real-FFT contract — `double` (the golden model) and `float` (the
embedded profile) on an engine that is the second template argument since
Stage 4 (the srdif engine, `detail::srdif_rdft`, by default; vDSP or
CMSIS-Helium where the build selects them for `float`), and the Q15 / Q31
fixed-point profiles from `sample_traits.h` on an int32 radix-4 kernel. The
contract is the packing, the sign convention, the scaling, and the numbers in
the table below; the engine behind it changed twice (the vendored Ooura C
until Stage 2b; the C++20 port of the same split-radix algorithm,
bit-identical to it, from Stage 2b, tap/DspTap#31, until tap/DspTap#42; the
srdif engine since, which moved every floating output bit at the error level
of either engine) and the contract did not; the radix-4 int32 kernel carries
fixed point. Consumers hold the transform by value and depend on exactly the
surface audit Part 1 enumerates: `basic_real_fft(size_t)`, `forward_inplace`,
`inverse_inplace`, `forward`, `inverse` (out-of-place, `2/N` applied),
`size()` / `num_bins()`, copyability, the packing, `exp(+i)`, the unnormalized
in-place inverse, and the aliases `real_fft` / `real_fft32`.

Out of scope here: the complex transform (`cdft`), the DCT/DST family, and any
consumer-side spectrum arithmetic (that is the Stage 5 packed-spectrum view).

## Contract summary, per profile

The header owns these numbers and the tests pin them; this table is the
cross-reference, and an entry that names a test is pinned by it. Floating-point
rows are the shipping `fft.h` since tap/DspTap#42 (the srdif engine, which
moved every floating output bit; before it, from the Stage 2b flip,
tap/DspTap#31, the port, whose outputs were bit-identical to the vendored
C's at `5ca3b1c`, set several of the pins below, and each is kept where the
srdif engine measures lower); fixed-point rows are the Stage 3b
kernel as merged in tap/DspTap#27, with the measurements in "The fixed-point
profiles" below.

| Contract point | `double` (`real_fft`) | `float` (`real_fft32`) | Q15 (`real_fft_q15`, `real_fft_q15_bfp`) | Q31 (`real_fft_q31`, `real_fft_q31_bfp`) |
|---|---|---|---|---|
| Engine | `detail::srdif_rdft<double>` (`fft/srdif.h`) by default (`default_real_fft_engine_t<double>`; `detail::split_radix_rdft<double>`, the port, from #31 until #42); since Stage 4 the second template argument, so any `real_fft_engine` can be named | `detail::srdif_rdft<float>` by default (the port from #31 until #42), or `detail::cmsis_real_fft_f32` (`fft/backends/cmsis.h`) / `detail::accelerate_real_fft_f32` (`fft/backends/accelerate.h`) when `TAP_DSP_FFT_CMSIS` / `TAP_DSP_FFT_ACCELERATE` is on; `test_fft_backend.cpp` is typed over the engines the leg can build and pins each to the srdif engine at float epsilon in one binary; `test_fft_routing.cpp` pins the class to the named srdif engine byte for byte on every leg and to the selected default | `detail::fixed_point_rdft<std::int16_t, Scaling>` (`fft/fixed_point.h`); the second argument is the Scaling policy on the fixed-point profiles | `detail::fixed_point_rdft<std::int32_t, Scaling>` |
| Size range (`k_min_size` … `k_max_size`, `supports_size(n)`; Stage 4) | 4 … 2^30 (the srdif engine's uint32 swap offsets reach 2^30 − 2; output fingerprinted 4 … 65536, oracle to 65536, the battery constructs to 2^20; the port's bound, the same number, was its `int` indexing) | srdif 4 … 2^30; vDSP 4 … 2^20; **CMSIS 32 … 4096** (`arm_rfft_fast_init_f32`'s table; status checked since Stage 4) — `EngineRangesAreTheStatedNumbers` as static_asserts in `test_fft_engine.cpp`, `ConstructsAtTheRangeBounds` | 4 … 65536 | 4 … 65536 |
| Packing | `data[0]` = DC, `data[1]` = Nyquist, `data[2k]` / `data[2k+1]` = bin *k* re / im, `1 ≤ k < N/2`; `N/2 + 1` bins | same | same (`ImpulseHasFlatSpectrum`, `DcAndNyquistPacking`, the oracle's closed forms) | same |
| Sign | `W = exp(+2πi/N)`; imaginary parts conjugated vs the engineering DFT | same | same (`SignConventionIsPlusI`, `OnBinSineIsPlusHalfNInTheImaginarySlot`) | same |
| Forward scale | 1 (`A[k] = Σ a[j] W^(jk)`, unnormalized) | 1 | `fixed`: exactly `X / N` in Q0.15, `e = log2 N` (`FixedForwardScaleIsExactlyXOverN`); `block_floating`: `X · 2^-e`, `0 ≤ e ≤ log2 N` returned (`BfpExponentIsWithinRange`) | `fixed`: exactly `X / 2N` in Q0.31, `e = log2 N + 1` — the one-bit input pre-shift (`FixedForwardScaleIsExactlyXOverTwoN`); `block_floating`: `X · 2^-e`, `0 ≤ e ≤ log2 N + 1` |
| Inverse scale | `inverse_inplace` unnormalized (caller applies `2/N`); `inverse()` applies `2/N` | same | the unnormalized inverse over `2^e`, same constant `e` under `fixed` (`FixedInverseCarriesTheSameExponent`); `inverse()` applies no `2/N`; round trip `x == out · 2^(e_fwd + e_inv + 1 − log2 N)` (`RoundTripReconstructsInputPerPolicy`) | same |
| Round-trip identity | `x` reproduced to `1e-12` abs at N = 1024 (`RoundTripReproducesInput`) | `2e-5` abs at N = 1024 (`RoundTripReproducesInput`) | 0.51 reconstructed LSB at N = 1024 (fixed; the output narrowing's half LSB), pinned 1.02 (`RoundTripReproducesInput`) | 3.34 reconstructed LSB at N = 1024 (fixed), pinned 6.7; 41.5 under block floating point (the DC bias), pinned 83 |
| Saturation-free input | n/a (floating point) | n/a | full scale ±1: input placed with 2 guard bits (`<< 14`, not 16); the int32 kernel performs no saturating operation for any input; the output narrowing clamps only at the rail itself (0.5 LSB, the full-scale Nyquist alternation); worst case vs the golden model 0.50 / 0.75 LSB (fixed / BFP), pinned 1.0 / 1.5 (`SaturationFreeWorstCaseDoesNotWrap`) | `fixed`: one input pre-shift (−6 dB); `block_floating`: full scale, the headroom rule; worst case 4.25 / 15.99 LSB at index 0 (fixed / BFP), pinned 8.5 / 32 |
| Noise floor | 1.60e-16 (1.6043e-16) relative 2-norm error of the forward vs the compensated-DFT oracle on full-scale uniform noise at N = 256, measured 2026-09-26 on the srdif engine on x86-64 Linux (g++ 13.3.0 -O3), 1.6972e-16 on the `cortex-m4-softfp` / `cortex-m4f` / `cortex-m33` QEMU legs and 1.4856e-16 on `cortex-m55` (arm-none-eabi-gcc 13.2.1, MinSizeRel; local runs), pinned at 7.4e-16 (`fft_oracle_floor.DoubleForwardTracksCompensatedDft`). The pin was set at 4× the port's value for the libm and fp-contraction spread across hosts, and is kept (the srdif tables are host-independent, contraction is not): the port measured 1.85e-16 (2026-09-23, x86-64 Linux, GCC 13.3.0 and Clang 18.1.3 -O3, glibc 2.39) and, in the CI logs of tap/DspTap#31 (run 35847675386), 1.8928e-16 on `cortex-m4-softfp` / `cortex-m4f` / `cortex-m33` (jobs 107137662986 / 107137663129 / 107137662923) and 1.6643e-16 on `cortex-m55` (107137663056). The Higham correctness envelope the oracle sweeps to N = 65536 is separate (`test_fft_oracle.cpp`; worst measured ratio 0.25 on the srdif engine) | 9.73e-8 (9.7312e-8) relative 2-norm error vs `double` at N = 512 on the srdif engine itself, so the same *engine* is measured on the M55 / macOS backend legs (where `basic_real_fft<float>` is CMSIS / vDSP), measured 2026-09-26 on x86-64 Linux (g++ 13.3.0 -O3, no FMA) and the soft-float `cortex-m4-softfp` leg, 1.0459e-7 on the VFMA legs `cortex-m4f` / `cortex-m33` / `cortex-m55` (local runs); the value moves in the last bits with fp-contraction (D9), no longer with libm. Pinned at 2.25e-7 (`RealFftCrossPrecision.FloatEngineTracksDoubleAtN512`), set at 2× the port's 1.1236e-7 (x86-64, soft-float M4) / 1.1650e-7 (VFMA legs, #31's CI run) and kept; and `< 1e-6` at N = 1024 through `basic_real_fft` (`FloatTracksDouble`). Against a quad-precision reference the srdif engine's float rms is 1.06e-7 / 1.17e-7 / 1.43e-7 at N = 512 / 2048 / 65536 ("The floating engine (srdif)") | fixed: 0.29 LSB rms (the narrow's own rounding) at every N and level; per-bin SNR 69.1 / 66.5 / 60.2 dB at 0 dBFS, N = 256 / 512 / 2048; BFP 87 – 91 dB at 0 dBFS (`NoiseFloorTracksWelchModel`; the full table in "The fixed-point profiles" §5) | fixed: 0.67 – 0.75 LSB rms, level-independent; 152.0 / 149.4 / 142.8 dB at 0 dBFS; BFP 157 – 161 dB |
| Latency | 0 (block transform, no internal delay) | 0 | 0 | 0 |
| Alignment | none required on `Sample*` | none (vDSP's internal split buffers are placed by the wrapper, not the caller) | none | none |
| Shareability across threads (`k_is_shareable`, an engine trait re-exported by the class; Stage 4) | **true**: tables built in the constructor, transforms `const noexcept` | **true** on the srdif engine, **false** on vDSP and CMSIS (per-object scratch); the class's transforms stay non-const on every profile (audit N11) and the trait is the statement; a shareable engine's transforms are `const` and the class `static_assert`s it (`ShareabilityIsTheHeadersNumber` in `test_fft_rt.cpp` and `test_fft_engine.cpp`) | **false**: the in-place int16 API needs the per-object int32 work buffer (a recorded deviation from audit Part 7, which listed both fixed profiles as shareable) | **true**: no mutable state during a transform; the engine's transforms are `const` behind `requires`-constrained overloads since Stage 4 |
| Real-time safety | transforms `noexcept`, allocation-free (`test_fft_rt.cpp`; no exception since the float-I/O-on-double overloads, `[[deprecated]]` at 2b, were removed at the D5 expiry, tap/DspTap#40) | same | same (`TransformsAreNoexcept`, `*AllocatesNothing`, `CopyProducesBitIdenticalOutput`) | same |
| NaN / denormals | NaN propagates to every bin; denormal input is slow on x86 without FTZ, and FTZ differs between the srdif, vDSP and CMSIS builds | same | not applicable / none | not applicable / none |
| Host identity | **one bit pattern for every compiler and target CI runs at `-ffp-contract=off`** since #42: the srdif engine's tables are integer arithmetic (`srdif_trig`, no libm), pinned as FNV-1a-64 output fingerprints at every power of two 4 … 65536 (4 … 4096 on the QEMU legs), **one row** (`tests/test_fft_srdif_fingerprint.cpp`; matched on x86-64 g++ 13.3 / clang++ 18.1 under both glibc dispatches, on the four QEMU legs, and in CI on Windows x64 MSVC and macOS arm64 AppleClang; scope in "The floating engine (srdif)"); with contraction the last bits move (D9, below). Until #42 not across hosts: the port's tables came from libm, whose last bit differs between glibc, newlib, UCRT and Apple (double pins per C library build, "The bit-identity record after D6") | same (until #42 the float row was already one value everywhere) | one bit pattern on every host, pinned per profile, policy, direction and N = 64 / 512 / 1024 / 2048 (`OutputFingerprintIsPinned`, `TwiddleTableChecksumIsPinned`) | same |

## Why split-radix for floating point and a radix-4 int32 kernel for fixed point

**Floating point is split-radix.** From Stage 2a (tap/DspTap#28, routed at
Stage 2b, #31) until tap/DspTap#42 it was Ooura's split-radix algorithm,
transliterated into C++20. The reasons then were contractual, not
aesthetic: every consumer float pin in MuTap and DspTap had been measured on
the vendored C, and the port reproduced the C's table semantics for both
precisions so that `float` stayed *bit-identical* (D10). The alternative —
computing the float twiddles in double "because it is more accurate" — was
measured and was a wash (table max absolute error 1.19e-7 either way;
transform rms relative error at N = 512 was 1.105e-7 with the C's tables and
1.114e-7 with double-computed ones — the audit's Part 6 N2 reviewer-probe
values; the committed re-measurement at Stage 2b was 1.12e-7 on the port),
so it was dropped. The `#define double float` build of the C was replaced by
a `template <std::floating_point Sample>` class whose helpers were private
static members, which removed the 76 global symbols and the lazy table build.

Since tap/DspTap#42 the floating profiles run the srdif engine
(`fft/srdif.h`). The maintainer decided to remove the last code derived from
Ooura's package from what DspTap ships ("Provenance and licensing", below),
which retired D10: every floating output bit moved, once, at the error level
of either engine, and consumers re-measure. The algorithm family stayed
split-radix, for the operation count — the engine attains the split-radix
minimum, 2N log₂N − 2N − 2 real operations per forward transform, which is
the cost on the soft-float Cortex-M4 — and the structure became the
fixed-point engine's in floating point (N/2 complex values, a
decimation-in-frequency kernel, bit reversal, #39's real post-pass), from
the literature (Duhamel and Hollmann 1984; Sorensen, Heideman and Burrus
1986). The twiddle question is answered differently: every table entry comes
from integer arithmetic and every float entry is the correctly rounded
value, so the float-versus-double table choice the port measured no longer
arises, and no libm reaches the output. "The floating engine (srdif)" has
the design and its numbers.

**Fixed point does not reuse Ooura's code, only its contract** (D2; from Stage
3b, #27, until #39 its real post-pass was nonetheless a transcription of the
package's, re-derived from the literature at #39 — "The fixed-point post-pass,
re-derived"). Audit Part 6 (at `5ca3b1c`) established that the port's
split-radix graph would survive only as a data-flow graph with a different
operation order, and that Ooura's first stage (`cftf1st`) derived half its
twiddles at run time as sums of two cosines up to 2.0, which no fixed-point
twiddle format holds. The fixed-point profiles are therefore their own radix-4
decimation-in-frequency complex kernel of length N/2 (radix-2 final stage for
odd `log2 M`) plus the real post-pass of the half-length method (Cooley, Lewis
and Welch 1970; Sorensen et al. 1987), arranged for the int32 arithmetic (§1),
so the packing, the `exp(+i)` convention and the DC / Nyquist slots are
identical. The design, summarized from audit Part 7:

- **One kernel, int32 data, two I/O widths.** Q31 is native; Q15 widens on
  input and narrows on output with the substrate's single round-half-up,
  saturating `finalize`. Fixed 1/N scaling in 16-bit data gives about 25 dB
  per-bin SNR at −40 dBFS for N = 512; in 32-bit data the same scaling leaves
  22 bits at N = 512, a floor below −130 dB. The cost is an int32 work buffer
  of N samples allocated at construction.
- **Q15 input placement uses the spare width as guard bits** (shift left by
  14, not 16), so full-scale ±1 has two guard bits and the worst-case radix-4
  component growth (4√2 per stage under rotation) cannot saturate under fixed
  scaling. Q31 has no spare width: under fixed scaling it takes one input
  pre-shift (−6 dB); under block floating point none is needed.
- **Twiddles are Q1.30 int32** for both profiles, a full table generated in
  double at construction and rounded once (`round_sat`), so 1.0 is
  representable and there is no special case.
- **Multiply is shift-before-multiply with one rounding**: `int32 × Q1.30 →
  int64`, then `>> 30` round-half-up — the single documented rounding point
  per product. Portable C++ first; Armv7E-M / Armv8-M `SMMULR` / `SMMLAR` is
  the designated seam behind the same contract.
- **Two scaling policies as a template parameter** (D3): `scaling::fixed`
  (2 bits per radix-4 stage, 1 per radix-2, one more in the post-pass; output
  exactly `X / N` with a statically known exponent; for magnitude and power
  consumers) and `scaling::block_floating` (a `clz` scan per stage shifts only
  as much as growth requires and returns one exponent per transform; for
  round-trip consumers that would otherwise lose `log2 N` bits twice).
- **The inverse has its own scaling.** The contract's unnormalized inverse
  (the convention of the vendored Ooura C, DspTap ≤ `5ca3b1c`, which every
  engine since keeps) has structural gain N/2 and saturates on the first stage
  if run unscaled in fixed point. Under `fixed` the inverse halves per stage
  and the round-trip factor is a fixed power of two; under `block_floating`
  the round trip returns `x · 2^-(e_fwd + e_inv)`. Whether the output
  convention matches CMSIS-DSP's *documented* q15 / q31 scaling is a
  compatibility decision made in 3b.

The literature this is implemented from, per the house IP policy (published
sources only; no shipping product's behaviour is reverse-engineered):

- P. D. Welch, "A fixed-point fast Fourier transform error analysis," *IEEE
  Transactions on Audio and Electroacoustics*, vol. AU-17, no. 2, pp. 151–157,
  June 1969.
- A. V. Oppenheim and C. J. Weinstein, "Effects of finite register length in
  digital filtering and the fast Fourier transform," *Proceedings of the
  IEEE*, vol. 60, no. 8, pp. 957–976, August 1972.
- A. V. Oppenheim and R. W. Schafer, *Discrete-Time Signal Processing*, 3rd
  ed., Prentice Hall, 2010, Section 9.7, "Effects of Finite Register Length"
  (the FFT round-off analysis and block floating point).
- J. W. Cooley, P. A. W. Lewis and P. D. Welch, "The fast Fourier transform
  algorithm: Programming considerations in the calculation of sine, cosine
  and Laplace transforms," *Journal of Sound and Vibration*, vol. 12, no. 3,
  pp. 315–337, 1970 (the transform of 2M real samples through one M-point
  complex transform: the fixed-point real post-pass and pre-pass).
- H. V. Sorensen, D. L. Jones, M. T. Heideman and C. S. Burrus, "Real-valued
  fast Fourier transform algorithms," *IEEE Transactions on Acoustics,
  Speech, and Signal Processing*, vol. ASSP-35, no. 6, pp. 849–863, June 1987
  (the same half-length method and its inverse).
- For the floating engine: P. Duhamel and H. Hollmann, "Split radix FFT
  algorithm," *Electronics Letters*, vol. 20, no. 1, pp. 14–16, 1984; H. V.
  Sorensen, M. T. Heideman and C. S. Burrus, "On computing the split-radix
  FFT," *IEEE Transactions on Acoustics, Speech, and Signal Processing*,
  vol. ASSP-34, no. 1, pp. 152–156, 1986; the full list is in
  `fft/srdif.h`. (The port the srdif engine replaced at tap/DspTap#42 cited
  Ooura's own references, listed in `third_party/ooura/readme.txt`: Nussbaumer
  1982; Burrus, *Notes on the FFT*.)

## The fixed-point profiles: design record (Stage 3b, tap/DspTap#27)

This section is the Stage 3b design note for `include/tap/dsp/fft/fixed_point.h`
and `include/tap/dsp/fft/tables.h`, folded in from the separate
`docs/fft-fixed-point.md` at Stage 2b (the fold-in carried two corrections
into the contract table above: the Q31 fixed forward scale is `X / 2N`, not
`X / N`, and Q15 is not shareable). Plan of record: audit Part 7 (design),
Part 9 (battery), Part 13 (the `fft_arith` amendments), Decisions D2 and D3.
The arithmetic specification the kernel is written against is
`fft/fft_arith.h`; every number here is one the committed battery
(`tests/test_fft_fixed.cpp`) pins or prints: a `ctest -V` run, or
`tap_dsp_tests --gtest_filter='fft_fixed*'`, reproduces each of them (the
`[ floor ]`, `[ measured ]`, `[ checksum ]` and `[ fingerprint ]` rows, on
every CI host and QEMU leg). The executed notebook and the analysis
instruments that re-measure them on speech-like material are Stage 3c.

Literature (house IP policy: published sources only; no shipping product's
behaviour was measured or reproduced):

- P. D. Welch, "A fixed-point fast Fourier transform error analysis," *IEEE
  Trans. Audio Electroacoust.* AU-17(2), 151–157, 1969.
- A. V. Oppenheim and C. J. Weinstein, "Effects of finite register length in
  digital filtering and the fast Fourier transform," *Proc. IEEE* 60(8),
  957–976, 1972.
- A. V. Oppenheim and R. W. Schafer, *Discrete-Time Signal Processing*, 3rd
  ed., Sec. 9.7 (round-off in the FFT; block floating point).
- J. W. Cooley, P. A. W. Lewis and P. D. Welch, "The fast Fourier transform
  algorithm: Programming considerations in the calculation of sine, cosine
  and Laplace transforms," *J. Sound Vib.* 12(3), 315–337, 1970, and H. V.
  Sorensen, D. L. Jones, M. T. Heideman and C. S. Burrus, "Real-valued fast
  Fourier transform algorithms," *IEEE Trans. Acoust., Speech, Signal
  Process.* ASSP-35(6), 849–863, 1987 — the half-length method the real
  post-pass and pre-pass implement, including the DC / Nyquist pair. The
  fixed-point arrangement of it (§1: one product per bin pair, the
  coefficient table, the roundings) is derived here from their equations;
  the complex kernel is DspTap's own (D2).

### 1. Kernel structure

The N real samples are read in place as M = N/2 complex values
z[j] = x[2j] + i x[2j+1]. The forward transform is

1. a radix-4 decimation-in-frequency complex FFT of length M in the library's
   convention W_M = exp(+2πi/M), with a radix-2 final stage when log2 M is
   odd (N = 4, 16, 64, ...);
2. the bit-reversal permutation of the M complex outputs
   (`make_bit_reversal_table`);
3. a one-bit shift (§2) and the real post-pass, derived below.

The inverse: the one-bit shift (with the input pre-shift, §2) and the real
pre-pass, then the conjugate kernel (same twiddle table, sine term
subtracted), then the permutation.

**Real post-pass: from the equations to the statements.** Library
convention X[k] = Σ x[n] W_N^(nk), W_N = exp(+2πi/N), M = N/2. The kernel
computes Z[k] = Σ_(n<M) z[n] W_M^(nk) of z[n] = x[2n] + i x[2n+1]. Write
E[k] and O[k] for the length-M transforms of the even and odd samples; both
sequences are real, so both transforms are conjugate-symmetric, and
Z = E + iO gives (Cooley, Lewis and Welch 1970; Sorensen et al. 1987)

    E[k] = (Z[k] + conj Z[M−k]) / 2,   O[k] = (Z[k] − conj Z[M−k]) / (2i),
    X[k] = E[k] + W_N^k O[k],   X[M−k] = conj(E[k] − W_N^k O[k])

(the second from E[M−k] = conj E[k], O[M−k] = conj O[k], W_N^(M−k) =
−conj W_N^k). The block entering the pass, after the kernel's shifts s and
the pass's own bit, is a[k] = Z[k] · 2^-(s+1), and the pass must leave
Y[k] = X[k] · 2^-(s+1). With u = a[k] and v = conj a[M−k],

    Y[k] = (u + v)/2 + (−i W_N^k)(u − v)/2 = u − C_k (u − v),   C_k = (1 + i W_N^k) / 2,
    conj Y[M−k] = (u + v)/2 − (−i W_N^k)(u − v)/2 = v + C_k (u − v),

so one complex product G = C_k (u − v) serves the pair (k, M − k):
Y[k] = u − G and Y[M−k] = a[M−k] + conj G. With θ_k = 2πk/N,
C_k = ((1 − sin θ_k)/2, cos θ_k / 2), both parts in [0, 1/2] for
1 ≤ k < M/2, and |C_k| = √((1 − sin θ_k)/2) ≤ 1/√2. The statements, per pair:
`gr = sub(u_r, a[M−k]_r)`, `gi = add(u_i, a[M−k]_i)` (u − v, exact), `rotate`
by C_k (two `mul_coeff` roundings per component, the trait's two-rounding
form), then `a[k] = (sub(u_r, gr), sub(u_i, gi))`, `a[M−k] = (add(a[M−k]_r,
gr), sub(a[M−k]_i, gi))` (exact). Each output component therefore carries
**two roundings** (the product's) on top of the one-bit shift before the
pass; the 1/2 of E and O never becomes a `shr_round`, because writing
(u + v)/2 as u − (u − v)/2 moves it into C_k. The alternatives weighed:
forming E and W_N^k O separately and halving costs a third, biased one-bit
rounding per component (+0.25 LSB mean, 3/4 · 1/12 variance), and
multiplying u and v by the two coefficients (1 − i W_N^k)/2 and
(1 + i W_N^k)/2 costs four `mul_coeff` roundings per component. By the
rounding counts alone (the §5 model's terms; neither alternative was built
and measured) both inject more noise per output component than the two
roundings chosen — 2.75 and 4 units of q²/12 against 2 — the first adds a
coherent +0.25 LSB bias on top, and neither lowers the bound (§3).

- **DC and Nyquist** (k = 0, whose partner is itself): E[0] = Re Z[0],
  O[0] = Im Z[0], W_N^0 = 1, W_N^M = −1, so X[0] = Re Z[0] + Im Z[0] and
  X[M] = Re Z[0] − Im Z[0]: `a[0] = add(a[0], a[1])`, `a[1] = sub(a[0], a[1])`
  (on the entering values), exact, no rounding.
- **Bin N/4** (k = M/2, its own partner): W_N^(N/4) = i, so C = 0 and
  Y[M/2] = a[M/2]; the pass leaves it untouched.

**Real pre-pass (the inverse).** Given X, set E[k] = (X[k] + conj X[M−k])/2,
W_N^k O[k] = (X[k] − conj X[M−k])/2 and Z[k] = E[k] + i O[k]; the
unnormalized length-M inverse kernel then returns Σ_k Z[k] W_M^(−nk) =
M z[n] = (N/2)(x[2n] + i x[2n+1]), which is exactly the golden model's
unnormalized inverse (round-trip gain N/2), so the pre-pass carries no
extra factor. With u = a[k], v = conj a[M−k]: Z[k] = u − conj(C_k)(u − v)
and conj Z[M−k] = v + conj(C_k)(u − v) — the forward's statements with the
conjugate coefficient (`rotate<true>`), the same table. At k = 0,
Z[0] = ((X[0] + X[M])/2, (X[0] − X[M])/2): `shr_round(add(a[0], a[1]), 1)` and
`shr_round(sub(a[0], a[1]), 1)`, one one-bit rounding per component — the
inverse's DC / Nyquist pair has a 1/2 no coefficient can absorb. k = M/2 is
again untouched (conj C = 0).

The derivation is written out in `fixed_point_rdft::real_post_pass`'s
docstring; the arrangement reproduces the output fingerprints and table
checksums the battery pinned before the re-derivation bit for bit (§4, §5).

**Radix-4 DIF butterfly**, for inputs a0..a3 at q, q + L/4, q + L/2,
q + 3L/4 of a sub-transform of length L, t0 = a0 + a2, t1 = a0 − a2,
t2 = a1 + a3, t3 = a1 − a3:

| output | value | written to |
|---|---|---|
| y0 | t0 + t2 | q |
| y1 | (t1 + i t3) · W_L^q | q + L/2 |
| y2 | (t0 − t2) · W_L^2q | q + L/4 |
| y3 | (t1 − i t3) · W_L^3q | q + 3L/4 |

(W_4 = +i forward, −i inverse; conjugate twiddles inverse.) y1 and y2 swap
places relative to the textbook radix-4 (radix-2² placement), so the
mixed-radix digit reversal collapses to a plain bit reversal even with the
trailing radix-2 stage, and one table serves every N. q = 0 skips the
rotation: `mul_coeff` by `k_coeff_one` is the identity, so the skip is
bit-identical, and the last radix-4 stage (L = 4) is multiply-free.

**Arithmetic.** Every operation is `fft_arith<Sample>::work`'s: `mul_coeff`
(int32 × Q1.30 → int64, `>> 30` with one round-half-up; a complex rotation is
two of them per output component, then `add`/`sub`; no fused form),
saturating `add`/`sub`, `shr_round` (round-half-up), `headroom_bits`. The
kernel defines no arithmetic of its own, so the trait's battery
(`tests/test_fft_arith.cpp`) covers every rounding point the transform has.

**Tables** (`fft/tables.h`, all generated in double and rounded once by
`make_coeff`): the kernel twiddles W_M^k for k in [0, M), interleaved
(cos, sin) Q1.30 — the stage with span L reads index k·(M/L) for r = 1, 2, 3,
which stays below 3M/4; the post-pass coefficients C_k = (1 + i W_N^k)/2
for k in [0, N/4), interleaved ((1 − sin θ_k)/2, cos θ_k/2) Q1.30
(`make_real_post_pass_table`: the real part is the integer 2^29 minus the
rounded sin θ_k/2, §4; k = 0 is exactly (1/2, 1/2) and unread, the pass
reads k in [1, N/4)); the bit-reversal table over M. Memory per transform
of size N: 4N bytes of twiddles, 2N bytes of post-pass coefficients, 2N
bytes of permutation, and for Q15 a 4N-byte int32 work buffer.

**Thread rule (a recorded deviation from the plan).** Audit Part 7 lists
`is_shareable` true for both fixed profiles, and the contract table above
carried that as the design until the fold-in. What shipped: one transform at
a time per object. Q15 is not shareable across threads, because the in-place
int16 API needs the per-object int32 work buffer; Q31 touches no object
state during a transform (the caller's buffer and the three const tables
only), but its transforms are declared non-const in this stage to match the
floating engine's shape (the split-radix port's transforms were const, as
the srdif engine's are since #42). Stage
4 states shareability as an engine trait (Q31 true, Q15 false) and decides
transform constness across the two wave-2 engines; nothing is restructured
in 3b. Every transform returns the exponent and is `[[nodiscard]]`: under
block floating point a discarded exponent is a silent scale error.

### 2. Halving count and the exponent

Shift-before-butterfly: each stage's inputs are `shr_round`ed before the
butterfly is computed (the trait's documented condition for the magnitude
bound; shift-after saturates on four rotated full-scale inputs at a 45°
twiddle).

| stage | growth (complex magnitude) | bits shifted under `scaling::fixed` |
|---|---|---|
| radix-4 | 4 | 2 |
| radix-2 | 2 | 1 |
| real post-pass / pre-pass | ≤ 2 per component, ≤ √2 in complex magnitude (see §3) | 1 |
| input pre-shift, first stage only | — | `k_fixed_scaling_input_pre_shift` (Q15: 0, Q31: 1) |

The kernel's stages total log2 M = log2 N − 1 bits; the post-pass adds one;
so the forward exponent is

    e_fixed = log2 N + k_fixed_scaling_input_pre_shift  = fixed_scaling_exponent(N),

and the forward output is exactly X / N (Q15) or X / 2N (Q31) in the
sample's Q format. The inverse pre-pass takes the pre-shift and its own bit
before the pre-pass; the kernel then halves per stage as in the
forward; the inverse exponent is the same constant. With G the double golden
model on the same input read as fractions, `G.forward == data · 2^e` and
`G.inverse_inplace (unnormalized) == data · 2^e`. Since the contract's
unnormalized round trip has gain N/2, a fixed-point round trip reconstructs

    x == out · 2^(e_fwd + e_inv + 1 − log2 N),

for both policies. Adjacent shifts are one `shr_round`: the Q31 pre-shift and
the first stage's two bits are a single 3-bit rounding, and the inverse's
pre-shift and pre-pass bit a single 2-bit one.

**Block floating point.** Before each stage (and, for Q15, before the final
narrow, treated as one more stage of growth 1) the block's `headroom_bits`
h is read and the shift is

    s = clamp(growth + 1 − h, 0, cum_fixed − e),

where cum_fixed is the fixed schedule's cumulative shift through this stage
and e the exponent applied so far. The "+1" leaves one bit of headroom
unconsumed (the trait's rule: a value of −2^20 reports 11 and −2^20 << 11 is
INT32_MIN, where `sub(0, x)` saturates). The upper clamp is what makes
0 ≤ e ≤ fixed_scaling_exponent(N) a contract rather than a typical value: it
binds only when e already equals cum_fixed and h = 0, where the fixed
schedule's bound (§3) applies to the same data with half a bit to spare. An
all-zero block returns e = 0 (the trait has no left shift; a quiet block is
never normalized upward).

**Reaching the constant is not the fixed schedule.** A block can be behind
the fixed schedule and catch up: a Q31 input with one bit of headroom (below
−6.02 dBFS) is shifted 2 at the first stage where fixed shifts 3, and a later
stage where the scan reads h = 0 shifts 3 against fixed's 2. The exponent
ends at the constant with a different rounding history, and the output is
then not bit-identical to the fixed one. A Q15 block always enters with
h ≥ 2 (the widen's guard bits), so every Q15 case that reaches the constant
is such a catch-up; its int32 history differs by the same few LSB32, which is
2^-14 of a Q15 LSB and survives the narrow only where the exact value sits on
a rounding tie. Measured over the battery's sweep (sizes 4 / 8 / 16 / 64 / 512 /
1024 / 2048, levels 0
to −8 dBFS in 0.2 dB steps, constant / tone / square exponentials at M/8 and
M/16 / noise, both directions, plus two Q15 tie cases found by a dense
amplitude sweep): Q31 differs on 140 of the 1162 inputs that reach the
constant, by at most 4 LSB (index 0, the DC path) and 2 LSB elsewhere; Q15 on
4 of 482 (the tie cases), by 1 LSB. Pinned at 2× by
`BfpAtTheFullExponentAgreesWithFixedWithinPin`; bit identity holds, and is
pinned, on the full-scale patterns only (`BfpAtTheFullExponentIsBitIdenticalToFixed`),
where a Q31 block enters with h = 0 and every stage shifts the fixed amount.
The general agreement between the two policies, brought to a common
exponent, is 1.0 LSB (Q15) / 4.0 LSB (Q31), pinned 2 / 8
(`BfpMatchesFixedAfterShift`).

### 3. Headroom: why nothing saturates

Let B be a bound on the complex magnitude of every value entering a stage.
Under shift-before-butterfly each of the four inputs is ≤ B/4 after the
shift, the intermediate sums (a0 ± a2, a1 ± a3) are ≤ B/2, the pre-rotation
sums ≤ B, the rotation is unit-magnitude, and every partial product
|x_r w_r| ≤ B; the outputs are ≤ B. So the complex magnitude never grows
across the kernel, and the bound is set by the input alone.

- **Forward kernel input**: a packed real pair with both components at full
  scale F, |z| ≤ √2 F. Q15: F = 2^29 (Q2.29 after `widen`), B = 2^29.5, 1.5
  bits below the int32 rail. Q31: F = 2^31 would give B = 2^31.5; the 1-bit
  pre-shift makes F = 2^30 and B = 2^30.5, half a bit below the rail.
- **Post-pass (forward).** Entering, after the one-bit shift, every value
  has |a| ≤ A with A = 2^29.5 (Q31 fixed: the kernel's B = 2^30.5 halved;
  Q15 fixed: 2^28.5; block floating point: below). Then, per pair
  (u = a[k], v = conj a[M−k], both ≤ A):
  - u − v: each component |u_r − v_r| ≤ 2A = 2^30.5 < 2^31 (the `sub` /
    `add` forming it cannot saturate), |u − v| ≤ 2A;
  - G = C_k (u − v): each partial product |(u − v)_r · C_r| ≤ 2A · 1/2 =
    2^29.5 (C's parts lie in [0, 1/2]), and |G| ≤ |C_k| · 2A ≤ 2A/√2 = 2^30,
    so G's components are ≤ 2^30 plus their two roundings;
  - Y[k] = u − G and Y[M−k] = a[M−k] + conj G: the pair is
    energy-preserving, |Y[k]|² + |Y[M−k]|² = (|u + v|² + |u − v|²)/2 =
    |u|² + |v|² ≤ 2A², so each output is ≤ √2 A = 2^30 in magnitude (plus
    the product's rounding, one LSB per component) — a growth of at most √2
    in magnitude, 2 per component, which the one bit pays for;
  - DC / Nyquist: |a_r ± a_i| ≤ √2 |a| ≤ 2^30, exact.

  Every intermediate and output is ≤ 2^30.5 + 1, half a bit below the rail
  (Q15: 2^29.5 + 1, 1.5 bits below). The bound on the outputs is attained by
  the golden model itself (a full-scale constant's DC, a full-scale Nyquist
  alternation's X[M], both exactly 2^30 in the Q31 frame); the bound on u − v
  is not reachable from real input: u − v = 2i O[k] · 2^-(s+1) is the odd
  samples' transform alone, |O[k]| ≤ M F, so |u − v| ≤ F = 2^30 (Q31), half a
  bit under the 2A the bullets allow — the forward post-pass is never the
  tight stage, and its worst cases are the constant and alternation
  patterns the battery already drives.
- **Pre-pass (inverse).** The spectrum is arbitrary, so the pre-pass's
  bounds are tight. It takes pre-shift + 1 bits: Q31 input components
  ≤ 2^31 → ≤ 2^29 after two bits (the rails shr_round to ±2^29), |a| ≤ 2^29.5
  = A; Q15 one bit, components ≤ 2^28, A = 2^28.5. The same four bullets
  hold with the conjugate coefficient: u − v ≤ 2^30 per component, partial
  products ≤ 2^29, |G| ≤ 2^30, |Z'[k]| ≤ √2 A = 2^30, and the k = 0 pair
  (a[0] ± a[1]) ≤ 2^30 before its halving, ≤ 2^29 after; the kernel then
  receives |Z'| ≤ 2^30, inside its B. A spectrum whose bin pairs are
  conjugate-antisymmetric at the rails (a[k] = (hi, hi), a[M−k] = (lo, hi):
  v = −u to within an LSB) makes u − v = 2u, |u − v| = 2^30.5 exactly, at
  every k including k = 1 where |C_k| is largest; the battery drives it in
  both sign orders (`antisymmetric pairs` in `SaturationFreeWorstCaseDoesNotWrap`).
- **Block floating point**: the rule shifts by s = growth + 1 − h (h the
  block's headroom before the shift) whenever that is positive, so before
  the shift every component lies in [−2^(31 − h), 2^(31 − h) − 1] and after
  the round-half-up shift in [−2^(30 − growth), 2^(30 − growth)]: closed at
  the top, because rounding can carry 2^(31 − h) − 1 up to the power of two
  itself, which then reads one bit less headroom (2^30 − 1 has h = 1, shifts
  by 1, and becomes 2^29, again h = 1). Complex magnitudes are
  ≤ 2^(30.5 − growth), outputs ≤ 2^30.5; the same half bit. For the
  post-pass and the pre-pass (growth 1) that is components in
  [−2^29, 2^29], |a| ≤ 2^29.5 = A, and the bullets above apply unchanged,
  for both Q formats (the Q15 block is int32 Q2.29 by then, the same
  arithmetic). The one exception is the upper clamp (§2): when the block
  has caught up with the fixed schedule (e = cum_fixed before the pass) and
  reads h = 0, the pass shifts one bit where the rule wanted two; the block
  is then the fixed schedule's data up to accumulated rounding, whose
  magnitude is the fixed bound, B / 2 = 2^29.5 = A again. So the post-pass
  and the pre-pass are saturation-free under both policies, in both
  directions, for both formats.

The worst case the battery drives is therefore the packed pair at full scale
(x[2j] = x[2j+1] = ±full scale, so every complex input has both components
at the rail), including the patterns whose sequence sits on the 45° twiddle
bins, full-scale binary noise, in both directions and under both policies
(`SaturationFreeWorstCaseDoesNotWrap`, sizes 4 / 8 / 16 / 64 / 512 / 1024 /
2048). Measured, as the
test prints it: the int32 kernel performs no saturating operation for any of
them (a clamp or a wrap would show as an output sitting on a rail the golden
model does not reach, and none does); the largest deviation from the golden
model at the returned exponent is 0.500 LSB (Q15 fixed), 0.750 LSB (Q15
block floating), 4.250 LSB (Q31 fixed: the inverse of a full-scale constant,
N = 2048, index 0) and 15.988 LSB (Q31 block floating: the inverse of
full-scale binary noise, N = 1024, index 0), pinned at 1.0 / 1.5 / 8.5 / 32
(unchanged, to the printed digit, by the post-pass re-derivation of
2026-09-26 and by the antisymmetric-pairs pattern it added). Those pins
cover the battery's sweep sizes, N = 4 / 8 / 16 / 64 / 512 / 1024 / 2048.
Above them the same sweep runs one size at a time in two host-only tests
(`SaturationFreeWorstCaseIsPinnedPerLargeSize`,
`RoundingBiasOnNegatedInputIsPinnedPerLargeSize`; compiled out where
`TAP_DSP_TEST_MAX_FFT_N` < 65536, i.e. on the QEMU legs; about 1 s for all
four configurations on the x86-64 host), each size pinned at 2× its own
measured maxima (2026-09-26, x86-64 Linux, GCC 13.3.0 -O3, glibc 2.39):

| N | 4096 | 8192 | 16384 | 32768 | 65536 |
|---|---|---|---|---|---|
| Q31 fixed, max \|out − G/2^e\| (LSB) | 4.704 | 4.250 | 4.250 | 4.741 | 4.724 |
| Q31 fixed, max F(x) + F(−x) (LSB) | 8 | 7 | 8 | 8 | 9 |
| Q31 block floating, max \|out − G/2^e\| (LSB) | 31.004 | 31.992 | 59.003 | 75.997 | 80.009 |
| Q31 block floating, max F(x) + F(−x) (LSB) | 121 | 126 | 137 | 141 | 286 |
| Q31 block floating, pinned deviation bound (LSB) | 62.01 | 63.99 | 118.01 | 152.0 | 160.02 |

The Q31 block-floating maximum grows with the gap between the constant and
the returned exponent (e = 9 … 12 against 13 … 17), always the inverse of
full-scale binary noise at index 0 (the DC-path bias of §5), and follows no
closed form (8192 sits 0.008 LSB under the 32 pin, 16384 nearly doubles it,
and F(x) + F(−x) doubles between 32768 and 65536), which is why it is a
per-size pin and not a law; "15.99 LSB, pinned 32" is a statement about the
sweep sizes only. The mean of F(x) + F(−x) stays at 0.45 – 0.67 LSB (Q31
fixed) and 0.93 – 1.00 LSB (Q31 block floating), inside the N ≤ 2048 bias
pins, and is pinned per size too. Both Q15 profiles stay inside their
N ≤ 2048 pins at every size (0.500 / 0.501 LSB, F(x) + F(−x) ≤ 1 LSB) and
are checked against them. No output sat on a rail the golden model does
not reach at any of those sizes.
The two Q31 maxima sit on the DC path, where the round-half-up biases add
coherently (§5); they are larger than the noise-floor maxima, not "the same
as on noise". The one place a saturating operation does fire is the Q15
output narrowing (`fft_arith<int16_t>::narrow`, round-half-up, saturating),
and it does so exactly when the true value is the rail: the full-scale
Nyquist alternation's X_{N/2}/N is 32767.5 LSB, exact through every shift,
rounds to 32768 and clamps to 32767, the pinned 0.500 LSB rail shortfall.

### 4. Twiddle quantization

Every coefficient is `make_coeff` of a double `cos`/`sin` (round half away
from zero), so |w_q − w| ≤ 0.5 LSB of Q1.30 = 2^-31 on every host. The
kernel twiddles go through libm for the first octant only (k ∈ [0, M/8],
with the k = M/8 sine set equal to its cosine) and reach the other seven
octants by the exact symmetries W^(M/4−k) = i conj W^k, W^(M/2−k) = −conj W^k,
W^(M−k) = conj W^k applied to the Q1.30 integers, so the k ↔ M − k conjugate
symmetry holds bit-exactly by construction rather than by two libm calls
agreeing, and a libm difference can enter through M/8 + 1 evaluations instead
of 2M. The angle 2πk/M is formed in double from the integer k; for the
kernel's M = N/2 ≤ 2^15 its error is < 2^-49 rad, and for the post-pass
table's n ≤ 2^16 < 2^-47 rad, five orders of magnitude below the quantum
either way. 1.0 is representable (2^30), so DC and Nyquist twiddles are exact
and need no special case. Per rotation the
relative error is ≤ √2 · 2^-31; over the ≤ 8 rotating stages of N = 65536
it is ≤ 6 · 10^-9 relative, which on a full-scale Q31 output (2^31) is under
13 LSB and on any Q15 output is invisible. On an on-bin full-scale tone the
measured Q31 error is 0.07–0.19 LSB rms (N = 2048 … 256), i.e. the twiddle
term is below the rounding floor at every N in the battery.

Two things could move a coefficient by one Q1.30 LSB between hosts: a libm
last-bit difference (glibc, newlib, UCRT, Apple), and fp-contraction of a
generator expression (Decision D9: a compiler that fuses `a - b*c` into one
rounding, which is the default on the macOS arm64 leg and on the Cortex-M55
leg, whose FP64 lets g++ contract). Either matters only for a double that
lies within its own error of a rounding boundary, and none lies that close.
Measured in quad precision (`__float128` `sinq` / `cosq` of the double angle
each generator forms, which is the same on every IEEE host; 2026-09-26,
x86-64, GCC 13.3.0 with libquadmath), the libm-evaluated entry nearest a
Q1.30 rounding boundary is 2.79e-5 LSB (2^-15.1) from it in the post-pass
table over every n ≤ 65536 (n = 65536, k = 2915, the sine), and 4.68e-5 LSB
(2^-14.4) in the kernel twiddles over every M ≤ 32768 (M = 32768, k = 767,
the cosine); up to N = 2048 the minima are 2.6e-3 / 2.9e-3 LSB. Counted in
each double's own ulp, the nearest entry is 812 ulp (post-pass) / 393 ulp
(kernel) from its boundary. So every table is the same on any host whose
`sin` and `cos` are within 2^8 ulp, and a fused `a - b*c`, which moves a
result by under one ulp (at most 2^-24 LSB here), could not move an entry
either; glibc's double path gives the rounding of the quad value at all
49,121 evaluated entries. Host identity of the tables is therefore
guaranteed for every size the profiles take (N ≤ 65536), and the pinned
checksums confirm it rather than being what establishes it.

The generators still contain no contractible expression, as D9 asks of all
code here. The post-pass generator (`make_real_post_pass_table`) has no
such expression by construction: its only doubles are the angle (2π/n, a
power-of-two scaling of the double π, times the integer k) and half a libm
`sin` or `cos` (a multiplication by the exact 0.5), and nothing is added to
a product in double. The real part (1 − sin θ)/2 is formed on the integers
as 2^29 − `make_coeff`(sin θ / 2), not as the double 0.5 − 0.5·sin θ (the
`a - b*c` shape D9 names), which also makes the complement C_k real +
C_(N/4−k) imaginary = 2^29 exact by construction. Like the kernel
twiddles it evaluates libm only on the first octant (k ≤ N/8, sine set equal to cosine at N/8) and
reaches k in (N/8, N/4) by swapping the quantized sine and cosine of
N/4 − k. Verified three ways: the reference in `TwiddleTableIsWithinHalfLsb`
(each entry within half an LSB of the double value, and the complement
C_k real + C_(N/4−k) imaginary = 2^29 bit for bit); the re-derived table
reproduces the post-pass checksums pinned on 2026-09-18 at N = 256 / 512 /
2048 on x86-64 (glibc, GCC 13.3.0); and the same checksums on the four QEMU
legs under newlib, including the Cortex-M55 leg whose FP64 lets g++
contract (the post-pass re-derivation's local reproduction of the legs,
2026-09-26). The battery pins each certified N's table checksum
(`table_checksum`, FNV-1a-64 over the bit patterns; `TwiddleTableChecksumIsPinned`)
so a remaining difference is detected, not absorbed, and on top of it the
transforms' own output fingerprints (`OutputFingerprintIsPinned`: FNV-1a-64
over the output block and the exponent, four profile × policy configurations,
both directions, N = 512 and 2048 and, since 2026-09-26, N = 64 and 1024
so the radix-2 final stage is fingerprinted too: thirty-two values), which
also fix the rounding form of the complex multiply, the shift schedule, the
sign convention and the packing as one bit pattern per host. The octant-symmetric
generator reproduces the checksums pinned on 2026-09-18 unchanged. The
kernel's outputs are otherwise integer arithmetic and host-identical.

### 5. The Welch noise model and the measured floors

**Model** (Welch 1969; Oppenheim & Weinstein 1972), for the fixed forward
schedule, output-referred in int32 LSB, per output component, variances
only: an s-bit `shr_round` of an integer injects (1 − 4^-s)/12 on every
component; a rotation injects 2/12 on each output component (two
`mul_coeff` roundings), on 3 of 4 outputs of a radix-4 butterfly (q = 0
skips), on none of the last radix-4 stage (L = 4) or the radix-2 stage;
noise injected before a stage propagates through each later radix-4 stage
with power gain 1/4 (four inputs each shifted 2 bits: 4 · (1/4)²), through a
radix-2 stage with 1/2, and through the post-pass with 1/4: its outputs are
the linear map (u, v) ↦ u(1 − C_k) + v C_k (and the mirror), and
|1 − C_k|² + |C_k|² = (|1 − i W_N^k|² + |1 + i W_N^k|²)/4 = 1, so a noise
(and a white signal) entering it reaches each output component with gain 1
after the shift's 1/4. The post-pass adds its shift's (1 − 1/4)/12 on every
component, and the product G's two `mul_coeff` roundings, 2/12, on both
bins of every pair: on the measured interior bins 1 … M − 1 all but N/4
(C = 0 there), a fraction (M − 2)/(M − 1); DC and Nyquist take no product.
The coefficient's own quantization acts on u − v, whose variance is twice
u's, so its signal-dependent term is twice the kernel rotation's: 8P
relative to the product's rounding term, P the per-component signal
variance at the pass in fractions of full scale: about 10^-4 on full-scale
white noise under fixed scaling, which is the battery's `signal` column at
0 dBFS (−37.92 / −40.60 / −46.93 dBFS per component for Q31 at N = 256 /
512 / 2048, i.e. 1.6e-4 / 8.7e-5 / 2.0e-5; Q15 −31.90 / −34.58 / −40.91,
6.5e-4 / 3.5e-4 / 8.1e-5), so 8P is at most 0.005 of the product's term
and the variance-only sum leaves it out. Summed:

| N | stages | variance-only prediction, LSB32 rms (the battery's `model` column) |
|---|---|---|
| 256 | 4-4-4-2 + post | 0.568 (−191.55 dBFS/component) |
| 512 | 4-4-4-4 + post | 0.584 (−191.31) |
| 2048 | 4-4-4-4-4 + post | 0.584 (−191.31) |

(The battery's `NoiseFloorTracksWelchModel` derives this model independently
from the two headers and prints it beside every measurement; the numbers
here are its printout. Summing the rules above by hand gives the same
0.584 at N = 512: 0.229 from the post-pass (0.0625 from its shift, 0.166
from its product on 254 of 255 interior bins), 0.078 from the last radix-4
stage through the post shift, 0.025 / 0.007 / 0.002 from the earlier ones.
Before the re-derivation of 2026-09-26 the model charged the post-pass's
product to every bin with the kernel rotation's twiddle term: 0.569 /
0.584 / 0.585, −191.53 / −191.31 / −191.30; the measurement did not move.)
The last radix-4 stage and the post-pass dominate; the input pre-shift and
the first stages are attenuated by 4^-(stages) and contribute nothing
measurable, which is why the Q31 floor is the same at every N and the
prediction converges at 0.58 LSB.

What the variance model does not carry is the **round-half-up bias**: an
s-bit `shr_round` has mean error +2^-(s+1) LSB (+0.25 for one bit, +0.125
for two). The 1-bit shift before the post-pass puts +0.25 on every
component entering it, which the pass carries to between 0 and +0.5 per
output component, bin by bin through C_k (+0.5 on DC, 0 on Nyquist:
a_r + a_i and a_r − a_i); the last kernel stage's 4 · 0.125 = 0.5 per component is halved
by the post shift and mixed with bin-dependent phases; earlier stages' biases
arrive rotated and count as variance across bins. Along the unrotated path
(q = 0, r = 0 at every stage) the biases add coherently into the DC bin.

**Measured** (the `[ floor ]` rows of `NoiseFloorTracksWelchModel`: white
noise, forward, output-referred against the double golden model on the same
quantized input, interior bins; rms and ratio at 0 dBFS, the exponent and
SNR per level; the level does not change the floor under fixed scaling, so
the 0 dBFS rms stands for 0 … −60 dBFS within the seed's spread; the on-bin
tone rows print beside them and are pinned by the same test):

| profile, policy | N | e at 0 dBFS | rms error, LSB | measured / model (power) | per-bin SNR at 0 dBFS | at −40 dBFS (e) |
|---|---|---|---|---|---|---|
| Q15 fixed | 256 | 8 | 0.29 | 1.02 | 69.1 dB | 29.7 dB |
| Q15 fixed | 512 | 9 | 0.29 | 1.01 | 66.5 dB | 26.6 dB |
| Q15 fixed | 2048 | 11 | 0.29 | 1.00 | 60.2 dB | 20.3 dB |
| Q31 fixed | 256 | 9 | 0.69 | 1.47 | 152.0 dB | 111.2 dB |
| Q31 fixed | 512 | 10 | 0.68 | 1.35 | 149.4 dB | 109.4 dB |
| Q31 fixed | 2048 | 12 | 0.70 | 1.44 | 142.8 dB | 103.2 dB |
| Q15 block floating | 256 | 5 | 0.29 | 1.01 | 87.2 dB | 77.5 dB (0) |
| Q15 block floating | 512 | 5 | 0.29 | 1.01 | 90.6 dB | 80.6 dB (0) |
| Q15 block floating | 2048 | 7 | 0.29 | 1.03 | 84.2 dB | 86.4 dB (0) |
| Q31 block floating | 256 | 7 | 0.99 | 1.79 | 160.8 dB | 155.7 dB (0) |
| Q31 block floating | 512 | 7 | 1.50 | 2.31 | 160.6 dB | 155.7 dB (1) |
| Q31 block floating | 2048 | 9 | 1.12 | 1.28 | 156.8 dB | 153.8 dB (2) |

Across all levels the Q15 fixed rms is 0.26 – 0.30 LSB and the Q31 fixed rms
0.67 – 0.75 LSB (ratio 1.31 – 1.75); Q31 block floating holds 135 – 136 dB
at −60 dBFS with e = 0, its rms in LSB rising to 2.3 / 3.2 / 6.6 (N = 256 /
512 / 2048) as the output is held at a larger scale. The pinned ratios
(largest measured / model over the sizes, levels and both materials, at 2×):
noise 2.11 / 3.51 / 2.09 / 4.63, tone 0.044 / 0.77 / 1.07 / 1.29 for Q15
fixed / Q31 fixed / Q15 BFP / Q31 BFP.

**Before and after the post-pass re-derivation (2026-09-26).** The
post-pass and pre-pass of the Stage 3b kernel were replaced by the
arrangement §1 derives from Cooley, Lewis and Welch (1970) and Sorensen et
al. (1987). Measured on the same host (x86-64, glibc 2.39, GCC 13.3.0 -O3)
with the same battery, the new code's outputs are the pinned bit patterns
(all sixteen `OutputFingerprintIsPinned` values of N = 512 / 2048 and the
three post-pass table checksums reproduce), so every floor is unchanged to
the printed digit: Q15 fixed 0.26 – 0.30 LSB rms before and after, Q31
fixed 0.67 – 0.75 LSB, Q15 block floating 0.28 – 0.30, Q31 block floating
0.99 / 1.50 / 1.12 LSB at 0 dBFS (N = 256 / 512 / 2048) before and after;
per-bin SNRs as in the table. What moved is the model (above), and with
it the measured / model ratios: Q31 fixed white noise 1.744 → 1.751 (pin
3.49 → 3.51), tone 0.383 → 0.384; Q31 block floating 2.321 → 2.315 (pin
4.65 → 4.63); the Q15 ratios did not move.

Reading it against the model:

- **Q15 fixed** measures the narrow's own rounding (1/√12 = 0.289 LSB16)
  and nothing else: the int32 kernel's 0.58 LSB32 is 2^-14 of a Q15 LSB.
  The per-bin SNR of −40 dBFS white noise at N = 512 is 26.6 dB; Part 7
  quoted "about 25 dB" for exactly this case.
- **Q31 fixed** measures 0.67 – 0.75 LSB rms against the 0.57 – 0.58 LSB
  variance-only prediction, a power ratio of 1.31 – 1.75. The excess is the
  round-half-up bias: the one-bit shift before the post-pass alone puts
  between 0 and +0.5 LSB on every output component (+0.25 on its inputs,
  carried through C_k), the last kernel stage's biases arrive half-shifted with
  bin-dependent phases, and the earlier stages' arrive rotated and count as
  variance across bins. The battery pins the bias directly as the mean of
  F(x) + F(−x) (`RoundingBiasOnNegatedInputIsBounded`): 0.81 LSB (Q31
  fixed), 1.05 LSB (Q31 block floating), pinned 1.63 / 2.1; the sum's
  maximum is 6 / 62 LSB (pinned 12 / 124), at index 0 / 1. Both directions
  are level-independent to within the seed's spread (0.67 – 0.75 at every
  level from 0 to −60 dBFS), which is the Welch statement; the inverse
  direction is pinned through the round trip (3.38 reconstructed LSB, pinned
  6.77) and the saturation sweep (4.25 LSB, pinned 8.5) rather than by a
  floor row of its own.
- **Block floating point** moves the floor with the block: Q15 gains 18 – 24
  dB at full scale and 48 – 66 dB at −40 dBFS over fixed; Q31 gains 9 – 14 dB
  at full scale and holds 135 – 136 dB at −60 dBFS with e = 0, where the
  input's own 2^-31 quantization is the floor. The Q31 BFP rms in LSB grows
  as e falls because the output is held at a larger scale; the SNR is the
  comparable number.
- **The DC bias (honest limit).** Along the unrotated path (q = 0, r = 0 at
  every stage) the biases add coherently into the DC bin: under fixed
  scaling each later 2-bit shift quarters what accumulated, under block
  floating point a stage typically shifts one bit against a magnitude gain
  of two, so the same bias is carried into a finer output LSB instead of
  being quartered. The battery sees it
  as the Q31 maxima all sitting at index 0 or 1: 4.25 LSB (fixed, N = 2048)
  and 15.99 LSB (block floating, N = 1024) in the saturation sweep, 6 / 62
  LSB in F(x) + F(−x), 41.5 reconstructed LSB in the BFP round trip (pinned
  83), with bins 1 – 3 and their mirrors carrying part of it through the
  slowly rotating early twiddles. Above the sweep sizes the host-only
  per-size tests of §3 pin its maxima at N = 4096 … 65536 (Q31 block
  floating: 31.0 – 80.0 LSB deviation, 121 – 286 LSB F(x) + F(−x)); the
  Stage 3c notebook (tap/DspTap#30) measures it from N = 256 to 65536, and
  the reading depends on the unit. Counted in the LSB of the constant exponent,
  the coherent sum at index 0 / 1 of F(x) + F(−x) on the battery's −0 dBFS
  noise is about 4 LSB under both policies and flat in N: fixed scaling 4,
  4, 6, 5, 6, 6, 8, 7, 9 LSB and block floating point 3.5, 3.3, 3.1, 4.0,
  4.2, 4.1, 4.3, 4.4, 4.4 LSB for N = 256 … 65536 (the fixed inverse of a
  full-scale constant levels off at 4.25 – 4.4 LSB the same way). Block
  floating point returns e = 7, 7, 8, 9, 9, 10, 10, 11, 11 there, 2 … 6
  bits under the constant, so in its own finer output LSB the same bias
  reads 14, 26, 25, 32, 67, 66, 137, 141, 279: it grows by the exponent gap
  2^(const − e), not with the stage count (the notebook's numbers,
  re-measured on the header at tap/DspTap#31). Q15 never sees it (below the
  narrow's quantum).
  Convergent (half-even) rounding would remove it; that would be a different
  `shr_round` contract in `fft_arith.h`, pinned by its own numbers, and is
  not this kernel's decision to make.

### 6. The CMSIS-DSP convention (D3, decided)

CMSIS-DSP documents its `arm_rfft_q15` / `arm_rfft_q31` as "internally input
is downscaled by 2 for every stage to avoid saturations", with a forward
output of X / N ("bits to upscale" = log2 N) and an inverse whose output
table lists the same widened format with no upscaling. The decision:

- The **Q15 fixed forward scale is CMSIS's**: X / N in Q0.15.
- The **Q31 fixed forward is X / 2N**, one bit below CMSIS's documented
  scale; the bit is the input pre-shift `fft_arith.h` requires because Q0.31
  has no guard bits, and it is stated as a number rather than hidden.
- The **inverse is not CMSIS's**: it is the floating contract's unnormalized
  inverse (the convention of the vendored C the library started from, kept by
  every engine since) over 2^e, i.e. (1/N) Σ X W^-jk divided by 2 (Q15) or 4
  (Q31) under fixed scaling. The round-trip identity above is the contract.
- **Block floating point has no CMSIS analogue.**

The exponent returned by every transform, not a fixed Q format per N, is the
interface; a consumer porting from CMSIS shifts by
`fixed_scaling_exponent(N) − log2 N` (0 for Q15, 1 for Q31) on the forward
and consults the round-trip identity on the inverse. Only CMSIS's
documentation was read; nothing of its behaviour was measured or reproduced.

### What Stage 3c adds

- `rfft_q15_512`, `rfft_q31_512`, `rfft_q31_2048` ratchet scenarios in
  `bench/`, seeded on `main` (a follow-up after the Stage 2b flip merges, so
  the fixed-point scenarios are seeded once, against the tree that ships).
- The analysis instruments typed over `std::span<const Sample>`, the capi and
  `dsptap_py` exposure of the four profiles, and `notebooks/fft.ipynb`
  executed against them, which re-measures the table in §5 on the shipping
  C++ and adds BFP exponent behaviour on speech-like material.

## The fp-contraction policy

A contract point, not a build detail (D9).

**Since tap/DspTap#42** the policy is unchanged — `tap::dsp` exports no
contraction setting, and bit identity across builds is claimed at
`-ffp-contract=off` only (`fft.h`) — and its evidence is the srdif engine's.
Built without contraction, the engine's output bits are a function of the
source alone (its tables are integer arithmetic, its statements IEEE
operations in a fixed order): one fingerprint row on x86-64 g++ 13.3 and
clang++ 18.1 under both glibc dispatches and on the four QEMU legs
(`tests/test_fft_srdif_fingerprint.cpp`, its own target at
`-ffp-contract=off`; Windows and macOS in CI). With contraction the last
bits move per build and the error statistics do not: the float-vs-double
relative error at N = 512 is 9.7312e-8 on x86-64 without FMA and on the
soft-float M4 leg, 1.0459e-7 on the VFMA legs (M4F, M33, M55), measured
2026-09-26. Everything below, from Stage 2a to tap/DspTap#41, was measured
on the port and the vendored C and is kept as the record; the engine-vs-C
identity it observed at default flags (every statement textually identical
on both sides) has no counterpart now, because there is no second
implementation to compare with.

Measured on the audit's baseline:
`gcc -std=c17` does not contract floating-point expressions into FMAs;
`gcc -std=gnu17` does; `g++` contracts in both `c++20` and `gnu++20`; clang
contracts in every mode, statement-scoped. Neither DspTap nor its consumers set
`-ffp-contract` or `CMAKE_C_EXTENSIONS`, so the C oracle and the port are
contracted differently by default wherever the ISA has FMA (Apple arm64, M55
VFMA, any x86 built with `-march`). Therefore:

1. The Stage 2a parity target compiled both the C and the C++ with
   `-ffp-contract=off`, and *that* was the bit-identity gate, for `double` and
   `float`, forward and inverse, at every power of two from 4 to 65536 plus one
   run at 2^20, same binary and same libm on each host. It ran from Stage 2a
   through Stage 4 and was retired with the reference C at D6; the pinned
   fingerprints that replaced it (`tests/test_fft_split_radix_fingerprint.cpp`
   from D6 until #42, "The bit-identity record after D6" below;
   `tests/test_fft_srdif_fingerprint.cpp` since) were, and are, built at
   `-ffp-contract=off` for the same reason.
2. The default-flags run was informational and pinned at a *measured* bound
   (its target went with the gate at D6; the table stays as the record) —
   a different fusion choice per stage accumulates over `log2 N` stages, and
   "1 ulp" is an assumption, not a number. Measured at Stage 2a on the head
   of tap/DspTap#28 (`8afe6cf`), 2026-09-18, the informational target at
   default flags on both sides, every material, forward and inverse
   (`fft_parity_ooura_default_flags.ReportsMaxUlpVersusOoura`, run with
   `-V` as its own CI step; CI run 35295206578, job ids in the table):

   | Platform (job) | Toolchain | Sizes | max ulp, `double` | max ulp, `float` |
   |---|---|---|---|---|
   | `linux-ooura` (105446230919) | ubuntu-24.04, GNU 13.3.0, x86-64 without `-march` (no FMA) | 4 … 65536, 2^20 | 0 | 0 |
   | `windows-ooura` (105446230835) | windows-2025-vs2026, MSVC 19.51.36256.0 (`/fp:precise`, no contraction) | 4 … 65536, 2^20 | 0 | 0 |
   | `macos-vdsp` (105446230934) | macos-26-arm64, AppleClang 21.0.0.21000101, arm64 (FMA ISA, clang contracts per statement) | 4 … 65536, 2^20 | 0 | 0 |
   | `cortex-m4-softfp` (105446230884) | arm-none-eabi-gcc 13.2.1 (15:13.2.rel1-2), QEMU 8.2.2, soft-float | 4 … 4096 | 0 | 0 |
   | `cortex-m4f` (105446230771) | same toolchain, fpv4-sp-d16 (VFMA) | 4 … 4096 | 0 | 0 |
   | `cortex-m33` (105446230920) | same toolchain, single-precision FPU (VFMA) | 4 … 4096 | 0 | 0 |
   | `cortex-m55` (105446230931) | same toolchain, MVE (VFMA); CMSIS on for `tap::dsp` but the parity binaries do not link it | 4 … 4096 | 0 | 0 |
   | local, this port's development host | g++ 13.3.0 and clang++ 18.1.3, x86-64 without `-march` | 4 … 65536, 2^20 | 0 | 0 |
   | local, at tap/DspTap#31 (2b review) | g++ 13.3.0 `-O3 -march=haswell` (FMA, GCC contracts across statements after inlining) | 4 … 65536, 2^20 | 0 | **> 0 at N = 1024, 4096, 16384, 65536, 2^20** (0 at 2048, 8192, 32768); the informational target's per-value ulp distance reaches 256 / 4 (fwd / inv) at 1024 and 16384 / 65536 at 2^20 on bins near zero; the review's independent harness on `basic_real_fft<float>` vs the C: 56 of 1200 blocks differ at N ≥ 1024, max \|Δ\| / max \|ref\| = 3.6e-7 |
   | local, at tap/DspTap#31 (2b review) | clang++ 18.1.3 `-O3 -march=haswell` (FMA, per-statement contraction) | 4 … 65536, 2^20 | 0 | 0 |
   | local, at tap/DspTap#34 (Stage 6 reviews; #34 touches neither `fft.h` nor `fft/`) | g++ `-O3 -march=x86-64-v3`, default `-ffp-contract=fast` | pvoc's float output through `basic_real_fft<float>` | codegen moved (the engine's fused multiply-add count changed in both precisions) | **float output bits moved with an edit to unrelated code in the same TU** (log_mel.h's constructor), every pvoc function instruction-identical; ~30 unrelated lines appended to the TU made the outputs identical again — the output depends on TU context, not only on flags |

   Zero everywhere, including the four FMA-capable legs (macOS arm64, M4F,
   M33, M55): because every statement is textually identical on the two
   sides, each compiler makes the same fusion choices for both. The bench
   binaries, built Release at default flags, say the same thing: the C and
   the port print identical output checksums on every Ooura key (see the
   instruction-count table). That is a property of these compilers on these
   statements, not a guarantee, which is why the gate stayed at
   `-ffp-contract=off` and the fingerprints are built there. Measured at 2b,
   and not claimed: x86-64 built with `-march` (FMA). g++ 13.3.0 `-O3 -march=haswell` keeps `double` identical
   and moves `float` at N = 1024, 4096, 16384, 65536 and 2^20 (the two rows
   above: a few float ulp relative to the block's peak, 3.6e-7 in the
   review's harness, large per-value ulp counts only on bins near zero);
   clang++ 18.1.3 at the same flags is identical for both precisions. That
   is what the objdump probe predicted: g++ fuses the two sides differently
   after inlining (588 fused instructions in a TU instantiating the port vs
   372 in the two C files, `-O3 -march=haswell`; clang 221 vs 333; both
   0 / 0 with the flag). The gate at `-ffp-contract=off` is 6 / 6 on both
   compilers at `-march=haswell`. The Stage 6 reviews (tap/DspTap#34, last
   row) sharpened this: under GCC's cross-statement contraction on an FMA
   target the engine's output depends on **translation-unit context** —
   what else the TU inlines decides how the engine's inner routines are
    inlined into one another, and with it which products fuse — so an edit to code that
   never touches the FFT can move a consumer's float bits, and the double
   codegen moved in the same experiment even though no double bit was seen
   to move; "double is identical under g++ with FMA" is therefore an
   observation about the builds measured, not a guarantee. Two experiments,
   kept apart: the engine-vs-C gate is bit-identical at `-ffp-contract=off`
   and, at default flags, on every CI platform, MSVC and the four QEMU legs
   included; the fingerprint A/B (pvoc / log_mel through the class, main vs
   branch) has been run on g++ and on the M33 leg and not on MSVC or
   AppleClang. `fft.h`'s D9 paragraph states the same. For #35 itself the
   35b review's A/B of this tree against `main` with `tools/fingerprint`
   (12 lines: pvoc float / double x 4 configs, log_mel float / double x log
   / pcen) is identical at g++ default flags, at `-O3 -DNDEBUG`, at `-O3
   -DNDEBUG -march=x86-64-v3` (FMA) and on the Cortex-M33 leg — the
   inline-namespace move of pvoc and log_mel did not perturb the TU enough
   to trip the #34 effect here, which is evidence, not a guarantee.
3. **Decided at Stage 2b, stated in `fft.h`'s class docstring: `tap::dsp`
   does NOT export `-ffp-contract=off`** (or any contraction setting) to its
   consumers, and the header sets none; the engine is compiled under whatever
   contraction the consumer's compiler applies, alike for every instantiation
   in that build. The evidence for leaving it: the table above (0 ulp at
   default flags on every CI platform, three of them VFMA), the bench
   checksums (C and engine identical on every Ooura key at Release flags), and
   the MuTap fingerprint gate through the 2b bump (14 rows byte-identical,
   float rows included, with no flag anywhere in MuTap's build). The reason
   not to export: an INTERFACE `-ffp-contract=off` reaches every consumer
   translation unit that includes `fft.h` and would pessimize the VFMA / FMA
   targets the float profile exists for (the M55, Apple arm64) across all of
   that code, in exchange for a cross-compiler bit reproducibility that libm's
   last-bit `cos` / `sin` differences already denied across hosts (true of
   the port; since #42 the tables are libm-free, and the reason stands on the
   pessimization alone). What the
   header promises is therefore the observation, not a guarantee: the identity
   held at default flags on these compilers on these statements; g++ on
   x86-64 with `-march` is measured to move the float profile's last bits
   (clang does not) and is not claimed.
   A consumer that needs bit reproducibility across its own compilers sets the
   flag on its own targets.

A related, separate fact, until #42: libm `cos` / `sin` differ in the last
bit between glibc, newlib, UCRT and Apple, so the port's *double* outputs
were not identical across hosts either (the per-libm double rows of "The
bit-identity record after D6"); a platform-independent table was a Stage 3b
candidate (the fixed-point twiddle generator landed there). tap/DspTap#42
closed it for the floating profiles: the srdif engine's tables come from
`srdif_trig`, integer arithmetic, and its fingerprint has one row for every
host.

## Transliteration rules for the port (Stage 2a until tap/DspTap#42; history)

The port and its header were deleted at tap/DspTap#42 (hand-off commit
`17db855`), and the srdif engine that replaced it is not a transliteration
of anything: its own constraints (statement order under D9, the integer
table generator, the leaf blocks) are stated in `fft/srdif.h`. The rules
below applied to the port from Stage 2a (tap/DspTap#28) to #42. They lived
in the engine header (`include/tap/dsp/fft/split_radix.h`) so the next person
would not undo them, and are recorded here because they were what the bit
identity depended on (held by the parity gate from Stage 2a until D6, by the
pinned fingerprints from D6 until #42).

- **Statement fidelity.** Every Ooura statement stayed textually intact. Clang
  contracts FMAs within a statement only; GCC contracts across statements after
  inlining. Refactoring `wk1r * x0r - wk1i * x0i` into a `cmul` helper, a
  lambda, or a policy call would have changed which products fuse and
  silently broken bit identity on one compiler or the other.
- **Explicit `static_cast<double>` on every libm call.** In C, `cos(delta * j)`
  with a float `delta` calls the double `cos`; in C++ `std::cos` of a float
  calls `cosf`, and the tables then differ from the C. The parity test at an N
  where `makewt` takes its `nwh > 4` branch caught this, but only because the
  port kept the C's table semantics.
- **Tables were built once in the constructor**, with the same precision
  semantics the C used for each instantiation (float-typed locals for the
  float instantiation).
- **Index types could widen to `std::size_t`** (every loop bound is a `< m` /
  `> 0` form over non-negative values); `-Wconversion` stayed on to catch a
  missed cast.
- **`m_w.data()` was loaded into a local once per transform** so aliasing
  analysis matched the C, which receives `w` by parameter.
- Scope: `rdft` and the ~2,580 lines it reaches (`makewt`, `makeipt`,
  `makect`, the `bitrv2*` family, `cftfsub` / `cftbsub` and the `cft*` leaves,
  `rftfsub`, `rftbsub`). Not ported: the DCT/DST family, `cdft`, the thread
  scaffolding. Nothing in the reachable set self-recurses (`cftrec4` is a
  `while` plus a `for`).

### What the Stage 2a port did with these rules (tap/DspTap#28, `8afe6cf`)

- **Index types were not widened.** Part 4 permits `std::size_t`; the port
  keeps the C's `int` throughout the helpers so the statements stay textually
  identical, and the single `std::size_t` → `int` narrowing is the
  constructor's. The tables are addressed through raw pointers loaded into
  locals once per transform (`m_ip.data()`, `m_w.data()`), so `a`, `ip` and
  `w` reach the helpers by parameter exactly as in the C. `-Wconversion`,
  `-Wshadow`, `-Wpedantic` and `-Werror` were on for every build.
- **The only textual changes to a kernel statement** are `double` →
  `Sample` in declarations, `const` on the table pointers, the dropped local
  prototypes and `USE_CDFT_THREADS` blocks, and one explicit cast in each of
  `rftfsub` / `rftbsub` (`wkr = static_cast<Sample>(0.5 - c[nc - kk])`, the
  narrowing the C's assignment performs). This was checked mechanically, not
  by eye: a normalizer that strips whitespace, maps `Sample` → `double` and
  undoes only those documented differences finds all 26 mechanically ported
  functions token-identical to `fftsg.c` (`DIFFERENCES: 0`). `makewt` and
  `makect` are the two hand-written functions, per the table-semantics rule.
- **The table-semantics rule is load-bearing, shown by mutation.** With the
  parity alias re-pointed at the port and one `makewt` line changed to
  `std::cos(delta * static_cast<Sample>(j))` (i.e. `cosf` in the float
  instantiation), the three float gates
  (`ForwardIsBitIdenticalToOouraFloat`, `InverseIsBitIdenticalToOouraFloat`,
  `LargeTransformIsBitIdenticalToOouraFloat`) fail and the three double gates
  pass; restored, 7 / 7.
- **The fp-contraction flag is load-bearing, shown by objdump.** The probe
  from the parity file's comment at `-O3 -march=haswell`, counting
  `vfmadd|vfmsub|vfnmadd|vfnmsub`, on the reference C (`fftsg.c` +
  `fftsg_float.c`) and on a TU instantiating the port for both precisions:

  | Compiler | default: C | default: port | `-ffp-contract=off`: C | `-ffp-contract=off`: port |
  |---|---|---|---|---|
  | gcc / g++ 13.3.0 | 372 | 588 | 0 | 0 |
  | clang / clang++ 18.1.3 | 333 | 221 | 0 | 0 |

- **The gate at `8afe6cf`**: memcmp identity, `double` and `float`, forward
  and inverse, N = 4 … 65536 plus 2^20, five materials, on linux (g++ 13.3),
  windows (MSVC 19.51), macos (AppleClang 21, arm64) and, at N ≤ 4096
  against newlib's libm, on the four QEMU legs (arm-none-eabi-gcc 13.2.1);
  locally also on clang++ 18.1. Every leg green on the first CI run of the
  PR.

### Compile-time cost of the flip (audit D8, for the record)

Audit D8 (no explicit-instantiation escape hatch) rested on Part 4's
estimate; from Stage 2b `fft.h` includes `fft/split_radix.h` (2,582 lines)
in every consumer translation unit, so the cost is now measured. At
tap/DspTap#31 on the development container (Intel Xeon @ 2.10 GHz, Ubuntu
24.04, g++ 13.3.0, clang++ 18.1.3; min of 5 compiles, load < 1): a TU that
includes `fft.h` preprocesses to 81,042 lines against 78,558 at `b08f6c6`
(`-std=c++20 -E`); a TU instantiating `real_fft` and `real_fft32` with one
`forward_inplace` each compiles in g++ `-O2` 0.69 → 1.73 s and `-O0` 0.68 →
0.89 s, clang++ `-O2` 0.83 → 1.36 s and `-O0` 0.75 → 0.84 s. The 2b
process review measured the same independently (g++ `-O2` 0.71 → 1.70 s,
`-O0` 0.76 → 0.90 s; clang++ `-O2` 0.82 → 1.35 s, `-O0` unchanged). About
one second per optimized TU that instantiates both profiles, consistent
with the audit's "1.3 s as C++ for the whole file"; D8 stands unless
MuTap's build shows a regression these numbers do not predict.

At tap/DspTap#42 `fft.h` includes `fft/srdif.h` (929 lines) in place of the
port's 2,582. Measured 2026-09-26 on the same kind of container (Intel Xeon
@ 2.10 GHz, 4 vCPUs, g++ 13.3.0, clang++ 18.1.3; min of 5 compiles; load
average 3–4, so read the differences, not the absolutes), a TU that
constructs `real_fft` and `real_fft32` at 512 and calls one
`forward_inplace` each, against the include tree of `7a58ebe` (the base of
#42): preprocessed 52,082 lines against 53,003; g++ `-O2` 1.01 s against
1.18 s, `-O0` 0.57 against 0.56; clang++ `-O2` 0.92 against 0.90, `-O0`
0.52–0.54 against 0.43–0.49 (two sessions, not attributed). A tenth of a
second at most, in either direction: nothing D8 would have to answer.

## Stage 4: the engine parameter and the ABI tag (tap/DspTap#35)

Audit F4 and Decision D4. What landed, what was measured, and what it does
not close.

Since tap/DspTap#42 the floating engine this section calls "the split-radix
engine" is `detail::srdif_rdft` (`fft/srdif.h`), and the ABI tag
`fft_split_radix` is `fft_srdif`. Where a passage states current behaviour
(the spelling table, the ranges, shareability, the test layout) it now names
the srdif engine, with the #35 value beside it; the evidence — symbol
listings, the loader test's output, the ratchet counts — is quoted as
measured at #35 and #41, on the port.

### The engine parameter

`basic_real_fft<Sample, Policy = detail::default_real_fft_policy_t<Sample>>`.
One parameter list, two meanings; every spelling in use compiles with its
meaning (pinned by `tests/test_fft_engine.cpp`), and the one pre-Stage-4
spelling that did not survive the D4 expiry is rejected by name:

| Spelling | Resolves to | Note |
|---|---|---|
| `basic_real_fft<double>`, `real_fft` | `<double, detail::srdif_rdft<double>>` (`detail::split_radix_rdft<double>` from #35 until #42) | always |
| `basic_real_fft<float>`, `real_fft32` | `<float, default_real_fft_engine_t<float>>` | the srdif engine (the port from #35 until #42), or `detail::cmsis_real_fft_f32` under `TAP_DSP_FFT_CMSIS`, or `detail::accelerate_real_fft_f32` under `TAP_DSP_FFT_ACCELERATE`; selected in exactly one place in `fft.h` |
| `basic_real_fft<float, detail::srdif_rdft<float>>` (`detail::split_radix_rdft<float>` from #35 until #42) | that engine | new at #35: an engine named explicitly, beside the accelerated default in the same binary |
| `basic_real_fft<float \| double, scaling::fixed>` | cannot be instantiated | the pre-Stage-4 spelling every profile shared. From #35 until the D4 expiry it resolved to the default engine as a **distinct type** from the one-argument form (same layout and code, different template arguments: +2,949 B x86-64 / +1,995 B M55 of duplicate wrappers when both were instantiated, 35a/F3). **Expired** (tap/DspTap#40): a `static_assert` in the class body whose message names `basic_real_fft<float>` / `basic_real_fft<double>`, and `detail::floating_engine_of` is gone. The template-id can still be named (an alias, a pointer) without error; only instantiating the class fails, so a dead `using` of it survives until something instantiates it. No in-tree writer since the #35 fix pass (the capi's seam uses `detail::default_real_fft_policy_t<Sample>`, so `fft_impl<float>` holds exactly `real_fft32`); none in MuTap or MuTap-Max, which never wrote it |
| `basic_real_fft<std::int16_t \| std::int32_t, Scaling>` | the fixed-point specialization | unchanged; the second argument is the Scaling policy |

Any other second argument on a floating `Sample` must satisfy the concept
`tap::dsp::real_fft_engine<Engine, Sample>`: constructible from the size,
copyable, two `noexcept` in-place transforms, `size()`, and the three contract
constants `k_min_size`, `k_max_size`, `k_is_shareable`. The four engines
satisfy it; `int` and `scaling::fixed` do not. Since the D4 expiry
(tap/DspTap#40) the class rejects `scaling::fixed`, cv-qualified or not, by
name, and that is the only diagnostic: for that spelling `engine` falls back
to the portable engine (the port at #40, srdif since #42) purely as error
recovery, so the rest of the class body stays well-formed while the
unconditional `static_assert` still rejects every instantiation. Measured on
the compile-fail fixture `tests/compile_fail/fft_legacy_spelling.cpp`
(`float`, `double` and `const scaling::fixed`):

| Compiler | Before the recovery (first cut of #40) | Shipped |
|---|---|---|
| clang 18.1.3 | one error, the D4 message | one error, the D4 message |
| g++ 13.3.0 | the D4 message, then 3 cascade errors (no `forward_inplace`, no `k_is_shareable`, a non-constant condition) | one error, the D4 message |

The ctests `fft_compile_fail.LegacySpellingIsRejectedWithTheD4Message.*`
(hosted GCC / Clang legs; `tests/compile_fail/expect_compile_failure.cmake`)
pin it: the fixture must fail to build, its output must contain the D4
message, and it must contain exactly one `error:`. Each condition was shown
to fail when broken: with the named `static_assert` deleted the fixture
compiles (the recovery engine is a working engine, which is why the test
exists); without the recovery g++ reports 5 errors; with a wrong include
path the build fails without the message. The alternative D4 kept
on record — no default, consumers name the engine — was not taken: the
consumers hold `basic_real_fft<Sample>` in forty-odd places and the default
is what makes the build define a build define.

Backends moved to `fft/backends/cmsis.h` and `fft/backends/accelerate.h`,
included by whoever selects them (fft.h under the define, or a translation
unit naming the engine). Every engine is constructed from its size; the
`init()` shape and `make_floating_engine` are gone. The CMSIS object library
is unchanged and stays PIC. The vDSP wrapper's duplicated "Computed at use
rather than cached" paragraph (audit F9; the Stage 6 "two duplicated fft.h
paragraphs" item) is one paragraph now.

### Size ranges as contract numbers

Each engine states `k_min_size` / `k_max_size`; the class re-exports them,
offers `supports_size(n)` (a power of two inside the interval, `constexpr`),
and states `supports_size(size)` as the constructor's precondition through
`TAP_EXPECTS`, checked before the engine is built so the engine's own
assertion is never the first to fire.

| Engine | Range | Where the number comes from |
|---|---|---|
| srdif (`double`, `float`; the port from #35 until #42) | 4 … 2^30 | since #42, the srdif engine's uint32 swap offsets, which reach 2^30 − 2 (`fft/srdif.h`); exercised: output fingerprints 4 … 65536, the oracle 4 … 65536, construction and a round trip at 2^20 (`ConstructsAtTheRangeBounds`). At #35 the same number was the port's int-indexing bound (n an `int`; every index and table offset below 2^31), exercised by bit identity 4 … 65536 and 2^20 (the gate) and the oracle. Above 2^20 the transform is the same statements over larger tables and is not separately measured. 2^30 rather than 2^20 so the precondition does not narrow what a consumer could construct before Stage 4 ("a power of two ≥ 4") |
| vDSP (`float`, Apple) | 4 … 2^20 | vDSP documents no maximum; 2^20 is the bound the float oracle has always stated for this engine and the size the port's parity gate ran to (Stage 2a until D6). Through the class it is swept to 65536 on the macOS leg and pinned bin-for-bin at 512 / 2048; 65536 … 2^20 is inside the range on vDSP's word, not on a measurement here |
| CMSIS-DSP (`float`, Cortex-M55) | **32 … 4096** | `arm_rfft_fast_init_f32` dispatches through a `switch` over exactly {32, …, 4096} and returns `ARM_MATH_ARGUMENT_ERROR` otherwise, leaving the instance uninitialized. Until Stage 4 the wrapper ignored that status and fft.h promised "≥ 4" for every profile: N = 4 HardFaulted in the Stage 2a oracle's first QEMU run on the M55 leg (tap/DspTap#24; audit Part 13) — one outcome of undefined behaviour, not a promise (below). The constructor now checks the status (`TAP_EXPECTS`, debug builds) as well as the range |
| fixed point (Q15, Q31) | 4 … 65536 | as at Stage 3b, now spelled through the same constants |

What construction does out of range: `TAP_EXPECTS` is the house precondition
(STYLE.md §4, `include/tap/dsp/detail/expects.h`): an assertion in a debug
build, nothing in a release build, never a throw (some Tap consumers build
`-fno-exceptions`). Every CI battery is Release / MinSizeRel, so in CI the
check is the constexpr predicate, pinned as `static_assert`s on every leg (N =
16 and 8192 rejected under CMSIS, accepted by the portable engine named in the
same binary (the port at #35, srdif since #42); 32 and 4096 constructed and
round-tripped under CMSIS), plus a Debug-only death test that runs in a local
Debug configure (it passes there). **Release-mode behaviour, decided
(35a/F1):** a violated precondition in a release build is undefined behaviour
exactly as the plain `assert` it replaces was, and no fault is promised. For
the CMSIS engine the instance stays zero-initialized, and the final audit
(2026-09-26) measured what that does under QEMU (mps3-an547) in Release;
re-measured for this record in Release and MinSizeRel at `ddadb74`, one probe
binary per size: N = 16 and 8192 return a wrong spectrum (max |CMSIS −
split-radix| 6.60 at peak 6.52, and 3593 at peak 3593) and **nothing faults**;
N = 4 returns a wrong spectrum and corrupts the heap, so a later unrelated
call fails — a HardFault after the transform in one probe, newlib's `Balloc
succeeded` assertion in the next `printf` in another. (The HardFault at N = 4
that Stage 2a's oracle met is the same undefined behaviour, not a contract.)
Real silicon is not measured. There is deliberately no defined fallback,
because a branch in the transforms would cost the hot path on every call for a
case the precondition excludes and would hand a consumer an object that
silently transforms nothing. `supports_size` is therefore the *mandatory* gate
wherever N comes from configuration. The capi's `dsptap_fft_create` applies it
per profile since the #35 fix pass (before that it hard-coded 4 and the
fixed-point 65536 and never read the float engines' range, so a vDSP or CMSIS
build of the capi would have passed 2^21 or 16 straight through); the CMSIS
wrapper also no longer narrows a size above `k_max_size` into the library's
`uint16_t` argument (65536 would arrive as 0). MuTap's config path is on the
bump checklist below. A repo-wide precondition policy remains its own plan;
`TAP_EXPECTS` is used for fft.h's power-of-two and size-range preconditions
(and the backends' own) and nothing else. `tests/test_fft_oracle.cpp` reads
each profile's sweep range from the class instead of carrying the CMSIS range
under `#if`.

### After the final audit: DspTap's own gates, and the CMSIS default

The final audit (2026-09-26, SYNTHESIS A1–A4) found the rule above applied
to consumers but not to DspTap's own configuration-driven classes, and the
CMSIS backend reaching targets it was not written for. Both are closed here
(tap/DspTap#41):

- **The config-driven classes ask their engine.** `log_mel_geometry::valid()`
  cannot know the engine (the geometry type is engine-independent and sits
  outside the ABI tag) and accepted any power of two ≥ 4; `basic_pvoc` stated
  a fixed [64, 2^28]. On the M55 build, `log_mel32` at `fft_size` 16 and 8192
  — both `valid()` — returns wrong features with no fault (measured under
  QEMU, Release: summed features 13.87 and −24.0 where the double profile at
  the same geometry gives 19.86 and 20.97). Now `basic_log_mel<Sample>::
  supports_geometry(g)` is `g.valid() && fft_type::supports_size(g.fft_size)`
  and is the constructor's `@pre`; `basic_pvoc<Sample>` states `k_min_size =
  max(64, engine min)` and `k_max_size = min(2^28, engine max)` with a
  `constexpr supports_size(n)`, the constructor's `@pre`. Numbers per build:
  log_mel double 4 … 2^30 everywhere, float 32 … 4096 under CMSIS, 4 … 2^20
  under vDSP; pvoc double 64 … 2^28 everywhere, float 64 … 4096 under CMSIS,
  64 … 2^20 under vDSP, 64 … 2^28 on the portable engine (the port at #41,
  srdif since #42). Both classes keep the precondition style they had (a bare
  `assert`; release builds do not check). The capi's `dsptap_log_mel_create`,
  its two setters and `dsptap_pvoc_create` gate on those predicates; the capi
  exposes the double profiles only, whose engine is the portable one on every
  build, so no host result changes (an `int` `fft_size` never exceeds 2^30).
  Pinned by `log_mel_test.SupportsGeometryIsValidityAndTheEngineSizeRange` and
  `pvoc_test.SupportsSizeIsTheEngineRangeInsideTheClassBounds`, typed over
  both profiles on every host; the float instantiations run on the four QEMU
  legs, where the M55 rejects 8192 (and 16, which pvoc's 64 rejects on every
  build).
- **`TAP_DSP_FFT_CMSIS` defaults ON only where the compiler targets
  floating-point Helium.** It was ON for every `CMAKE_SYSTEM_NAME Generic` +
  `arm` target; DspTap's M4 and M33 toolchains escaped only by pinning it OFF,
  and a consumer's own M7 toolchain got the CMSIS engine silently (float range
  4 … 2^30 → 32 … 4096, tag `fft_split_radix` → `fft_cmsis`; the portable tag
  is `fft_srdif` since #42). The root CMakeLists.txt now compiles a check on
  `__ARM_FEATURE_MVE & 2` under the configured flags (a static-library
  `try_compile`, so it needs no linker script on bare metal, whatever the
  toolchain file set `CMAKE_TRY_COMPILE_TARGET_TYPE` to). Configure results,
  before (`ddadb74`) → after, `TAP_DSP_FFT_CMSIS` with no `-D`: M55 (in-tree
  toolchain) ON → ON; M33, M4 soft-float, M4F OFF → OFF (before: by the
  toolchain's pin; after: by detection, the pins deleted); synthetic Cortex-M7
  (`-mcpu=cortex-m7 -mfpu=fpv5-d16 -mfloat-abi=hard`) **ON → OFF**; host
  x86-64 OFF → OFF; MuTap's own M55 toolchain file (`origin/main`, unedited,
  consumer-style `add_subdirectory`) ON → ON. Also measured after:
  `-march=armv8.1-m.main+mve.fp+fp.dp` with `CMAKE_SYSTEM_PROCESSOR
  cortex-m55` ON, `-mcpu=cortex-m85` ON, `cortex-m55 -mfloat-abi=softfp` ON
  (`__ARM_FEATURE_MVE` 3, and `tap_dsp_fft` builds — the review's
  measurement), `cortex-m55+nomve.fp` (integer MVE only) OFF, `cortex-m55
  -mfloat-abi=soft` OFF, an M55 toolchain without
  `CMAKE_TRY_COMPILE_TARGET_TYPE` ON. Bare-metal Cortex-A / R (`Generic` +
  `arm`) targets move ON → OFF like the M7. What the check sees, stated so
  nobody reads more into it: the C++ compiler with `CMAKE_CXX_FLAGS` only —
  not `CMAKE_C_FLAGS` (the CMSIS objects compile as C, so keep the CPU flags
  identical in both, as every toolchain here does), not
  `CMAKE_CXX_FLAGS_<CONFIG>`, not per-target options, and not the caller's
  `CMAKE_REQUIRED_*`, which the check resets (`cmake_push_check_state(RESET)`;
  review repro: a consumer's `CMAKE_REQUIRED_DEFINITIONS
  -D__ARM_FEATURE_MVE=3` on the M33 toolchain selected CMSIS, and
  `CMAKE_REQUIRED_FLAGS -mfloat-abi=soft` on the M55 deselected it; both now
  detect from the toolchain). The cached result is cleared before every check,
  so a re-configure with other flags re-detects (review repro: M55 toolchain,
  then re-configure with M33 `CMAKE_CXX_FLAGS` / `CMAKE_C_FLAGS` — before, the
  M55 result stayed cached and the Helium status line printed; now detection
  is 0, and since the option keeps its cached ON the "ON without MVE-F"
  warning fires). Every blind spot fails safe (OFF, the portable engine) or
  needs deliberately split flags. An explicit `-DTAP_DSP_FFT_CMSIS=ON|OFF`
  still wins: ON on a 32-bit Arm target without MVE-F configures with a
  warning (CMSIS-DSP builds its non-Helium C, same 32 … 4096); ON for a target
  that is not 32-bit Arm — x86, and arm64 / Apple Silicon, which the old
  `arm|ARM` regex let through — is a configure error. The embedded CI job
  asserts the detected value per leg (`TAP_DSP_FFT_CMSIS:BOOL=ON` on the M55,
  `OFF` on the other three), so a broken detection cannot turn the M55 leg
  into a second portable-engine leg unnoticed.
- **Two false statements withdrawn.** The "HardFault" promise (above, and
  in `fft.h` and `fft/backends/cmsis.h`) and "`int` is `std::int32_t` on
  all CI legs" (Shareability, below).

### Shareability as an engine trait

`k_is_shareable` (the plan's `is_shareable`, under the house `k_` prefix): the
portable engine true (the port at #35, srdif since #42), Q31 true, vDSP /
CMSIS / Q15 false. Decision on constness across engines: **shareable implies
const** — a shareable engine's transforms are `const`, and the class
`static_assert`s it, so that half of the trait is what the compiler sees
rather than a comment; the converse is not asserted (an engine with const
transforms over `mutable` scratch that says `k_is_shareable = false` is
accepted, and honest) and is a convention the tests pin: for the six shipped
instantiations const and shareable coincide (`test_fft_rt.cpp`,
`test_fft_engine.cpp`). The Q31 fixed-point transforms became `const` behind
`requires`-constrained overloads (`k_widens` selects the Q15 pair, which goes
through the per-object work buffer and stays non-const). Also noted (35a/F7),
and corrected by the final audit and its review: `int` is `std::int32_t` on
the three hosted CI legs (glibc x86-64, MSVC, Apple arm64) and under clang for
arm-none-eabi (clang 18.1.3 `--target=arm-none-eabi -mcpu=cortex-m55`:
`__INT32_TYPE__` is `int`), so there `basic_real_fft<int>` is the Q31 profile.
Under arm-none-eabi-gcc — all four QEMU legs — `__INT32_TYPE__` is `long int`
(the compiler's own macro, from GCC's newlib-stdint target configuration, not
from newlib's headers), `std::int32_t` is `long`, and `basic_real_fft<int>`
does not compile (the primary template's `static_assert`, measured with
arm-none-eabi-g++ 13.2 for `-mcpu=cortex-m55`, `cortex-m33` and `cortex-m4`;
`basic_real_fft<long>` is the Q31 profile there).
`fft_engine.Int32IsLongOnlyUnderArmNoneEabiGcc` pins which on every leg. Spell
the profile `std::int32_t`. The class's own transforms stay non-const on every
profile (audit N11: the public API's constness must not depend on the selected
engine). The `std::span` overloads Part 9 listed were not added: Stage 4 added
no overload to the transform surface, so there is nothing to equate
(`test_fft_rt.cpp` says so where Part 9 expected the assertion).

### The ABI tag: what it is, the evidence, what it does not change

`fft.h` opens `namespace tap::dsp::inline TAP_DSP_FFT_ABI` around
`basic_real_fft` and the aliases, with `TAP_DSP_FFT_ABI` one of
`fft_srdif` / `fft_cmsis` / `fft_vdsp` from the selection. At #35 the first
was `fft_split_radix` (the plan named them `fft_ooura` / `fft_vdsp` /
`fft_cmsis`, and #35 took `fft_split_radix` because the default had been the
port since 2b); tap/DspTap#42 renamed it with the engine, because the
floating layouts changed, so images built against the two engines do not
coalesce each other's symbols.
`pvoc.h` and `log_mel.h` define `basic_pvoc` and `basic_log_mel` inside the
same inline namespace, because they hold a `basic_real_fft<Sample>` by value.
`k_real_fft_abi_tag` exposes the tag as text.

**Why.** `basic_real_fft<float>`'s object layout follows the selected engine
(`m_engine` *is* the engine), and so does the layout of every class that
holds one by value. Their member functions are weak (COMDAT) definitions in
every image that instantiates them. macOS dyld coalesces weak definitions
across all loaded images — every image's `WEAK_DEF` references to one
mangled name are bound to a single chosen definition, unless the symbol is
not exported (hidden visibility) — and on ELF a `dlopen`ed plugin resolves
its references against the executable's and `RTLD_GLOBAL` objects' exports
before its own. So two images built with different defaults, loaded into one
process, would run one image's code over the other image's layout: F4's
cross-image hazard. Not observed today (MuTap-Max forces vDSP off through
MuTap; the capi is never loaded into Max), real in principle.

**What an inline namespace changes and does not change.** It changes one
thing: the mangled name of every entity declared inside it. Unqualified and
`tap::dsp`-qualified lookup is unchanged (an inline namespace's members are
members of the enclosing namespace for lookup; `using tap::dsp::basic_real_fft;`
in MuTap's `mutap/fft.h` works as before), overload resolution and ADL are
unchanged, the object layouts are unchanged, no code is generated. With the
engine a template argument, `basic_real_fft<float>`'s own symbols already
differ between two builds without the tag (the argument is in the mangled
name); the tag is what separates the embedders'.

**Evidence.** One translation unit (`template class` instantiations of
`basic_real_fft<float>`, `<double>`, `<std::int16_t>` and `basic_pvoc<float>`)
compiled twice for the Cortex-M55 with arm-none-eabi-g++ 13.2.1 (`-O2`,
`-mcpu=cortex-m55 -mfloat-abi=hard`), with and without `-DTAP_DSP_FFT_CMSIS`,
`arm-none-eabi-nm -C` on the objects:

- At `8350f13` (the branch base): **50 of 50** weak (`W` / `V`) symbol names
  are identical in both objects, over two different layouts (each object
  also carries 4 non-weak COMDAT group-section `n` symbols, the
  constructors', not counted), among them

  ```
  W tap::dsp::basic_pvoc<float>::process(float, float)
  W tap::dsp::basic_real_fft<float, tap::dsp::scaling::fixed>::forward_inplace(float*)
  W tap::dsp::basic_real_fft<short, tap::dsp::scaling::fixed>::forward_inplace(short*)
  ```

- At tap/DspTap#35 (on the port; since #42 the default object's names read
  `fft_srdif` and `detail::srdif_rdft<float>` where these read
  `fft_split_radix` and `detail::split_radix_rdft<float>`): **0 of 65** weak
  names are shared (and 0 of all 69 lines, the 4 `n` symbols included). The
  same three read

  ```
  W tap::dsp::fft_split_radix::basic_pvoc<float>::process(float, float)
  W tap::dsp::fft_split_radix::basic_real_fft<float, tap::dsp::detail::split_radix_rdft<float> >::forward_inplace(float*)
  W tap::dsp::fft_split_radix::basic_real_fft<short, tap::dsp::scaling::fixed>::forward_inplace(short*)
  ```
  in the default object and
  ```
  W tap::dsp::fft_cmsis::basic_pvoc<float>::process(float, float)
  W tap::dsp::fft_cmsis::basic_real_fft<float, tap::dsp::detail::cmsis_real_fft_f32>::forward_inplace(float*)
  W tap::dsp::fft_cmsis::basic_real_fft<short, tap::dsp::scaling::fixed>::forward_inplace(short*)
  ```
  in the CMSIS one. The 35a review's recount with a TU that also
  instantiates `basic_log_mel<float>`: 96 of 107 weak names shared at the
  base, 31 of 122 at this branch, every remaining one in `std::` or
  `tap::dsp::detail::` — the engines and tables themselves, whose layout does
  not depend on the define, so their coalescing is harmless.
  `tests/test_fft_engine.cpp` pins the same property from inside the battery
  on every leg: `typeid(real_fft32).name()`, `typeid(pvoc32).name()` and
  `typeid(log_mel32).name()` contain the expected tag and none of the other
  two.

  **And the loader half, measured in a process (`tests/test_fft_abi_tag.cpp`,
  linux and macOS; asked for by the 35a review, whose experiment it
  reproduces).** One translation unit (`tests/abi/abi_tag_image.cpp`) is built
  into two loadable modules: image A with the configured build's defines (tag
  `fft_split_radix` on linux at #35, `fft_srdif` since #42; `fft_vdsp` on
  macOS) and image B with `TAP_DSP_FFT_CMSIS` and a stub engine of a different
  layout on the include path (tag `fft_cmsis`); the host `dlopen`s A then B
  with `RTLD_GLOBAL` and asks each image what its embedders see through
  `extern "C"` probes, the member calls laundered through a volatile pointer
  so the loader, not the optimizer, decides the target. Measured on linux (g++
  13.3): an embedder of `basic_real_fft<float>` defined in plain `tap::dsp`
  and instantiated in image B **ran image A's code** — its `which()` returned
  `fft_split_radix` from inside the `fft_cmsis` image and the bound member saw
  a layout of 80 bytes against the real 184 — which is F4 live (as measured at
  #35); the identical embedder defined inside `inline namespace
  TAP_DSP_FFT_ABI` reported its own tag and its own layout in both images, and
  `pvoc32` / `log_mel32` computed identical, finite checksums in both images
  with both loaded (`TaggedEmbedderIsImmuneToCrossImageCoalescing`; the
  untagged outcome is asserted on the linux/GCC host and printed as `[
  measured ]` elsewhere, `UntaggedEmbedderIsWhatTheTagExistsFor`). Under
  `RTLD_LOCAL` on ELF nothing coalesces (the review measured it; the test runs
  one mode per process because the first `dlopen` decides an image's
  bindings). The test is hosted, non-Windows: the QEMU legs have no loader,
  and a Windows image exports nothing without `__declspec`, so the coalescing
  has no path there.

**The fixed-point profiles are inside the tag, deliberately stated.** A
partial specialization is the same template as its primary and lives in the
same namespace, so `basic_real_fft<std::int16_t, scaling::fixed>` carries
the tag although its layout does not depend on the selection. What that
costs: two images built with different float defaults in one process each
instantiate their own copy of the fixed-point code (a coalescing that would
have been harmless does not happen), and if they ever exchanged a
`real_fft_q15` across their boundary by type the link would fail rather
than silently merge. What it does not cost: any correctness, any
instruction in a transform (the fixed-point icount scenarios are +0.00 % on
every key). The alternative — the fixed-point profiles outside the tag —
requires splitting `basic_real_fft` into two templates or an alias-template
indirection, which breaks `basic_real_fft<std::int16_t, scaling::block_floating>`
as a class-template spelling (partial specializations on it in
`test_fft_rt.cpp`, `::engine` in the capi); not worth it for a hazard the
fixed-point profiles do not have.

**What the tag does not close: MuTap's own embedders, until the bump.**
The classes that hold a `basic_real_fft<Sample>` by value in `tap::mu` (at
MuTap `801204d`, per the 35a review) are `partitioned_fdaf` (fdaf.h:431),
`partitioned_fdkf` (fd_kalman.h:534), `pem_afc` (pem_afc.h:204, plus its
`Core`), `residual_suppressor` (postfilter.h:726) and `nn_suppressor`
(nn_suppressor.h:425); `aec_chain` (postfilter.h:787) is a second-level
embedder that holds its `Canceller` / `Post` as template arguments and so
inherits the tag once they carry it. Their symbols do not carry the tag
until MuTap does the following; until then the hazard F4 describes is
closed for DspTap's own embedders and for `basic_real_fft` itself, and open
one level up in MuTap exactly as before.

**MuTap bump checklist (the pin that picks up #35):**
1. Wrap each of the five headers' class definitions in `namespace
   tap::mu::inline TAP_DSP_FFT_ABI { … }` (a different namespace from
   `tap::dsp::fft_srdif` — `fft_split_radix` before #42 — which is fine: the
   tag only has to appear in the mangled name; the macro is in scope through
   `mutap/fft.h`). postfilter.h holds two of the six classes.
2. **No forward declaration of any tagged class anywhere** — a later
   `template <typename> class partitioned_fdaf;` in plain `tap::mu` declares
   a *different* class. None exists today in MuTap, MuTap-Max or DspTap
   (grepped by the review); it must stay that way. Explicit specialisations
   and `using` through the enclosing namespace remain legal.
3. **Validate N with `basic_real_fft<Sample>::supports_size(n)` wherever a
   block size comes from configuration**, with the CMSIS numbers (32 …
   4096) in the error text: on the M55 a config block size outside that
   range is undefined behaviour in a release build otherwise — measured as
   a silently wrong spectrum (N = 16, 8192) or a corrupted heap (N = 4),
   with no fault promised. For a DspTap class that holds the FFT itself,
   gate on the class's predicate (`basic_log_mel<Sample>::
   supports_geometry`, `basic_pvoc<Sample>::supports_size`), which asks the
   same engine.
4. Expect `partitioned_fdaf<double>` and the other double instantiations to
   carry a tag their layout does not need (the tag is keyed on the float
   default) — the same cost class as the fixed-point paragraph above; no
   correctness effect.
5. Name the Xcode / AppleClang that matters for MuTap-Max on the bump: the
   macOS CI runner here is AppleClang 21.0.0, not 15/16. MuTap-Max's
   externals hold these classes in non-template classes, one dylib each,
   all from one build — no cross-build exposure there.
6. Gate as always: fingerprints byte-identical on every leg.

**Not closed by Stage 4, and not the tag's business:**
- CMSIS-vs-portable-engine parity runs under emulation on the Cortex-M55
  QEMU leg, as it has since that leg landed (main's `test_fft_backend.cpp`
  already compared `basic_real_fft<float>` — CMSIS under the define — to
  the port as reference); what #35 added is that the two engines are
  typed rows in one binary (`fft_backend_parity/cmsis` beside
  `/split_radix`, `/srdif` since #42) and that the named portable-engine
  routing row runs there.
  The remaining gap: QEMU only, never a host, never hardware.
- The vDSP same-binary microbenchmark (above, "Host microbenchmark").

### The ratchet, and one thing it caught

Measured locally before the PR with the CI toolchain (arm-none-eabi-gcc
13.2.1, qemu 8.2.2, the same counting plugin as `bench.yml`; the baseline
binaries rebuilt from `8350f13` reproduce their recorded counts to within
the +35 / +75 instructions `bench/README.md` already records for the DONE
line):

| Key | `rfft_f32_512` | `rfft_f32_2048` | Q15 / Q31 scenarios |
|---|---|---|---|
| m55 (CMSIS), commits 1–4 as first pushed | 58,667,632 (**+12.00 %**) | 61,148,069 (**+11.47 %**) | +0.00 % |
| m55 (CMSIS), final (harness factory) | 52,382,329 (−2) | 54,858,158 (+38) | +5 each |
| m33 | 100,935,605 (−0.01 %) | 115,460,510 (−0.00 %) | Q15 −2,043, Q31 +5 |
| m55-ooura | 89,268,128 (−0.01 %) | 102,783,071 (−0.00 %) | +5 each |
| m4f | 97,138,550 (+6) | 111,278,422 (+6) | +5 each |
| m4-softfp | 1,864,929,145 (+4) | 2,294,359,239 (−1,020) | +5 each |

The +12 % was not in the transform: both float checksums were identical to
the baseline binaries', the CMSIS archive was byte-identical, and the
wrapper's copy-back loops were identical instruction for instruction. It
was the bench's own FNV-1a fold. With the CMSIS engine's constructor
(vector allocation, `memset`, `arm_rfft_fast_init_f32`, the EH cleanup)
inlined into `main` — the baseline had reached it through an out-of-line
`make_floating_engine<float>` — GCC's register allocation of `main` changed
and the fold loop went from 7 to 10 instructions per element (the 64-bit
hash's low word spilled to `[sp, #28]`, the FNV prime `movw r0, #435`
rematerialized every element): +6 per sample per iteration, exactly the
delta. **The defect was the harness**, and the fix is not a re-record:
`bench/icount/icount_main.cpp` constructs the transform under test through
a non-inlined `make_fft()` in both `run()` bodies — the baseline's own shape
— so the count is the transform plus a fixed fold. A first version had put a
`noinline` attribute on the two accelerated engines' constructors instead
(m55 +5 / +45); the hostile review measured the harness fix landing closer
with nothing in the library, and the attribute costing 860 bytes of
MinSizeRel `.text` on the one-construction m55 f32 probe, so shipping code
carries no benchmark-shaped attribute and the split-radix constructor was
never touched. The `.text` probes (MinSizeRel, `size -A`) all stay under
their ceilings: m55 f32 106,725 (−932 vs the recorded 107,657), Q15 27,105
(+8), Q31 26,553 (−40); m33 44,009 / 31,113 / 30,593; m55-ooura 39,281 /
27,105 / 26,553; m4f 44,601 / 31,745 / 31,209; m4-softfp identical to the
record. The harness lesson is in audit Part 13.

### Test layout after Stage 4

- `tests/test_fft_backend.cpp`: typed over `::testing::Types<srdif_f
  [, accelerate_real_fft_f32 | cmsis_real_fft_f32]>` (`split_radix_f` from
  #35 until #42); parity at 512 / 2048, alignment stability and tonal
  accuracy at 512 / 2048 / 4096, per engine, against the srdif engine called
  directly; a `static_assert` that
  the default float engine is one of the rows.
- `tests/test_fft_engine.cpp` (new): the engine-parameter contract as
  compile-time facts on every leg (above).
- `tests/test_fft_routing.cpp`: memcmp class-vs-engine for `double`, for the
  named srdif float engine on every leg (the port until #42), and for the
  selected float default.
- `tests/test_fft_rt.cpp`: `ShareabilityIsTheHeadersNumber` per
  instantiation, with the compiler-checked half.
- `tests/test_fft_parity_ooura.cpp`: untouched at #35, green on every leg —
  the proof that no float or double bit moved through the refactor (deleted
  at D6, tap/DspTap#36, with the reference C; "The bit-identity record after
  D6" below).
- `tests/test_fft_abi_tag.cpp` + `tests/abi/` (hosted, non-Windows): the
  two-image loader test of the tag (above).
- Since #42, beside these: `tests/test_fft_srdif_fingerprint.cpp` and
  `tests/test_fft_srdif.cpp` ("The floating engine (srdif)", "Tests changed
  for the engine").

## Size and instruction counts, per target

Seeded by the Stage 1b ratchet from the vendored C, then re-measured at 2b
(port), 3b (fixed point), 4 (engine parameter) and tap/DspTap#42 (the srdif
engine). Numbers here carry the SHA, toolchain and QEMU versions; the gate
itself lives in `bench/baselines.json` and the CI job, never in this file.
The bench key `m55-ooura` (the Cortex-M55 with `TAP_DSP_FFT_CMSIS=OFF`, i.e.
the portable float engine on the M55) keeps the name it was seeded under,
because the baselines are keyed on it; it has measured the port since 2b and
the srdif engine since #42. Every float-engine figure below the #42 entries
is the port's (or the vendored C's, where the table says so).

### `.text` per profile at N = 512

**The srdif engine (tap/DspTap#42; 2026-09-26, local; confirmed by CI on `main`).**
MinSizeRel float probe (`tap_dsp_size_probe_rfft_f32_512`, `size -A`
`.text`), arm-none-eabi-gcc 13.2.1 (15:13.2.rel1-2), 2026-09-27 (the tree
after the fix pass for review A of #42), the port's Stage 2c figure in
parentheses: `m4-softfp` 31,081 (53,289), `m4f` 28,921 (44,601), `m33`
28,305 (44,001), `m55-ooura` 28,161 (39,281). The float profile no
longer links libm's double `sin` / `cos` (the tables are integer
arithmetic), which is most of the 11–22 kB. Ceilings re-recorded to
measured + 3 %, rounded up to 64 bytes: 32,064 / 29,824 / 29,184 / 29,056 (from
54,912 / 45,952 / 45,376 / 40,512). The Q15 / Q31 probes and the `m55`
(CMSIS) probe are unchanged. The push-to-`main` bench run on the #42
squash, `72977aa` (run 36348787347), measured the same four figures to the
byte (`bench/README.md`, "Sizes").

**Post-pass re-derivation (2026-09-26; local, not a CI run).** The
fixed-point post-pass / pre-pass re-derived from the literature (§1),
measured locally with the CI toolchain (arm-none-eabi-gcc 13.2.1
(15:13.2.rel1-2), QEMU 8.2.2 (1:8.2.2+ds-0ubuntu1.18)), MinSizeRel, `size -A`
`.text`, the float probe unchanged. Against the Stage 2c record (and the
Stage 4 local figures where they differ):

| Target | Q15 | Q31 | ceilings q15 / q31 |
|---|---|---|---|
| `m4-softfp` | 31,473 (unchanged) | 31,001 → 30,937 | 32,448 / 31,936 |
| `m4f` | 31,681 (2c; 31,745 at Stage 4) → 31,681 | 31,209 → 31,145 | 32,640 / 32,192 |
| `m33` | 31,105 (2c; 31,113 at Stage 4) → 31,113 | 30,633 (2c; 30,593 at Stage 4) → 30,585 | 32,064 / 31,552 |
| `m55` | 27,097 (2c; 27,105 at Stage 4) → 27,097 | 26,593 (2c; 26,553 at Stage 4) → 26,553 | 27,968 / 27,392 |
| `m55-ooura` | 27,097 (2c; 27,105 at Stage 4) → 27,097 | 26,593 (2c; 26,553 at Stage 4) → 26,553 | 27,968 / 27,392 |

Every probe is under its ceiling; the float probes read 53,289 / 44,601 /
44,009 / 106,725 / 39,281, the Stage 4 figures. No ceiling is re-recorded.

**Stage 2c (this PR): the fixed-point columns, and the ceilings.** Measured
on the recording run of the 2c branch (`21b3488`, bench run 35861317022,
2026-09-23, arm-none-eabi-gcc 13.2.1 (15:13.2.rel1-2), QEMU 8.2.2 (1:8.2.2+ds-0ubuntu1.18), ubuntu-24.04), MinSizeRel,
`size -A` `.text`. The float probe is byte-identical to the Stage 2b figure
on every key (the shipping float path did not change; the `_c` probe is gone
with the C). The fixed-point probes are `basic_real_fft<std::int16_t>` and
`<std::int32_t>` under `scaling::fixed`: the same int32 kernel at both I/O
widths, so Q15 is ~500 bytes over Q31 (the narrow/widen glue), 27-32 KB per
profile against the float engine's 39-53 KB, and identical on the two M55
keys (no backend applies). The ceilings are now set in `bench.yml`
(`text_ceiling_f32/q15/q31`), policy measured + 3 % rounded up to 64 bytes
(`bench/README.md`, "Sizes"): the Part 4 assertion for the float
instantiation on the M55 leg is `40,512` bytes on
`m55-ooura` (the engine) and `110,912` on `m55` (CMSIS-DSP).

| Target | float (engine) | Q15 | Q31 | ceilings f32 / q15 / q31 |
|---|---|---|---|---|
| `m4-softfp` | 53,289 (split-radix) | 31,473 | 31,001 | 54,912 / 32,448 / 31,936 |
| `m4f` | 44,601 (split-radix) | 31,681 | 31,209 | 45,952 / 32,640 / 32,192 |
| `m33` | 44,001 (split-radix) | 31,105 | 30,633 | 45,376 / 32,064 / 31,552 |
| `m55` | 107,657 (CMSIS-DSP Helium) | 27,097 | 26,593 | 110,912 / 27,968 / 27,392 |
| `m55-ooura` | 39,281 (split-radix) | 27,097 | 26,593 | 40,512 / 27,968 / 27,392 |


**Stage 2b (the flip; tap/DspTap#31).** Measured on the seeding run of the
2b branch (`b078a80`, bench run 35844483811, 2026-09-23, arm-none-eabi-gcc
13.2.1 (15:13.2.rel1-2), QEMU 8.2.2 (1:8.2.2+ds-0ubuntu1.18), ubuntu-24.04).
The shipping probe is `basic_real_fft<float>` (the split-radix engine, or
CMSIS on `m55`); the `_c` probe is the vendored C through the bench adapter
(`reference_c_bench_adapter`), informational until Stage 2c. The `.text`
ceilings in `bench.yml` are still 0 (not gated): the seeding commit (#26) did
not set them and this PR does not either; a ceiling policy (slack, what a
re-record looks like) is a follow-up.

| Target | float, shipping (engine) | float, `_c` (vendored C via adapter) | shipping / C | SHA / job |
|---|---|---|---|---|
| `m4-softfp` | 53,289 (split-radix) | 53,225 | 1.0012 | `b078a80`, job 107127217360 |
| `m4f` | 44,601 (split-radix) | 44,793 | 0.9957 | job 107127217630 |
| `m33` | 44,001 (split-radix) | 44,185 | 0.9958 | job 107127217777 |
| `m55` (CMSIS on) | 107,657 (CMSIS-DSP Helium) | 39,449 | 2.7290 (vs the C) | job 107127217635 |
| `m55-ooura` | 39,281 (split-radix) | 39,449 | 0.9957 | job 107127217645 |

The shipping probe's `.text` is the Stage 2a port probe's to the byte on
every Ooura key (53,289 / 44,601 / 44,009 → 44,001 / 39,281; the M33 figure
moved by 8 bytes with the wrapper). The `_c` probe is a different binary from
#28's C probe (the adapter carries its own copy loop and constructor), so
its column is not the 51,729 / 43,153 / 42,601 / 38,505 recorded below; the
port-vs-C size comparison of record stays the Stage 2a one. The CMSIS probe
lost 24 bytes (107,681 → 107,657): the wrapper no longer carries the two
unused Ooura workspace vectors under a backend define.

**Stage 2a (for the record; the numbers the 2b comparison was judged against).**

Bytes in the `.text` row of `arm-none-eabi-size -A` on the MinSizeRel size
probe (`bench/size_probe.cpp`: startup + one transform + what it pulls in;
no stdio). The port column is informational at Stage 2a (not a ceiling until
2b routes the port). Measured on the head of tap/DspTap#28 (`8afe6cf`), bench
run 35295262684, 2026-09-18, arm-none-eabi-gcc 13.2.1 (15:13.2.rel1-2),
QEMU 8.2.2 (1:8.2.2+ds-0ubuntu1.18), ubuntu-24.04.

| Target | float, C | float, port | port / C | Q15 | Q31 | double (host only) | SHA / toolchain |
|---|---|---|---|---|---|---|---|
| `m4-softfp` | 51,729 | 53,289 | 1.0302 | | | n/a | `8afe6cf`, gcc 13.2.1, run 35295262684 job 105446394297 |
| `m4f` | 43,153 | 44,601 | 1.0336 | | | n/a | same, job 105446394232 |
| `m33` | 42,601 | 44,009 | 1.0331 | | | n/a | same, job 105446394276 |
| `m55` (CMSIS on) | 107,681 (CMSIS-DSP Helium, not Ooura) | 39,281 | 0.3648 (vs CMSIS) | | | n/a | same, job 105446393983 |
| `m55-ooura` | 38,505 | 39,281 | 1.0202 | | | n/a | same, job 105446394192 |

The port costs 0.8 – 1.6 KB more `.text` than the C on the Ooura keys
(+2.0 % to +3.4 %); the two `m55` rows carry the same port probe, so the
39,281 is one number measured twice.

Pre-ratchet reference points from the audit (thumbv8.1m, hard float,
rdft-reachable text, `--gc-sections`): float 15.7 KB at `-Os`, 19.9 KB at
`-O2`; double (soft-float) 39.9 KB at `-O2`.

### Instructions per scenario

**The srdif engine (tap/DspTap#42; 2026-09-27, local, the tree after the fix
pass for review A of #42; confirmed to the instruction by the push-to-`main`
bench run on the #42 squash, `72977aa`, run 36348787347).**
`scripts/icount.py` as `bench.yml` runs it; the port's baselines → the
srdif engine's counts:

| Scenario | `m4-softfp` | `m4f` | `m33` | `m55` (CMSIS, unchanged) | `m55-ooura` |
|---|---|---|---|---|---|
| `rfft_f32_512` | 1,864,929,141 → 1,814,303,702 (−2.71 %) | 97,138,544 → 92,575,609 (−4.70 %) | 100,945,841 → 92,979,212 (−7.89 %) | 52,382,331 → 52,382,329 (−0.00 %) | 89,276,321 → 83,933,381 (−5.98 %) |
| `rfft_f32_2048` | 2,294,360,259 → 2,243,212,021 (−2.23 %) | 111,278,416 → 106,282,752 (−4.49 %) | 115,465,626 → 106,773,786 (−7.53 %) | 54,858,120 → 54,858,158 (+0.00 %) | 102,784,096 → 97,602,563 (−5.04 %) |

The float baselines on the four portable-engine keys are re-recorded to
these counts (six of the eight moved beyond the −3 % band, and every
checksum changed with the output bits); the fixed-point scenarios measured
+0.49 … +1.15 % against theirs with no change to that code and are not
re-recorded. Against CMSIS-DSP on the M55 the portable engine now executes
1.60× (N = 512) / 1.78× (N = 2048) the instructions (the port 1.70× /
1.87×). How the engine got there is in "The floating engine (srdif)".

**Post-pass re-derivation (2026-09-26; local, not a CI run).** Same
toolchain and QEMU as the sizes above, the counting plugin built as
`bench.yml` builds it, `scripts/icount.py` in compare mode against
`bench/baselines.json` (the Stage 2c recording of the previous post-pass):

| Scenario | `m4-softfp` | `m4f` | `m33` | `m55` | `m55-ooura` |
|---|---|---|---|---|---|
| `rfft_q15_512` | 746,311,273 → 752,864,296 (+0.88 %) | 753,215,637 → 760,026,298 (+0.90 %) | 752,189,619 → 758,994,136 (+0.90 %) | 684,157,511 → 688,580,080 (+0.65 %) | 684,157,511 → 688,580,080 (+0.65 %) |
| `rfft_q31_512` | 709,161,574 → 714,608,662 (+0.77 %) | 715,019,447 → 720,988,333 (+0.83 %) | 714,626,257 → 719,567,047 (+0.69 %) | 657,595,969 → 665,145,823 (+1.15 %) | 657,595,969 → 665,145,823 (+1.15 %) |
| `rfft_q31_2048` | 869,187,566 → 873,980,755 (+0.55 %) | 876,486,909 → 881,802,089 (+0.61 %) | 875,801,377 → 880,074,125 (+0.49 %) | 806,142,635 → 813,689,587 (+0.94 %) | 806,142,635 → 813,689,587 (+0.94 %) |

The float scenarios read −0.01 … +0.00 % (the Stage 4 figures). Every
fixed-point scenario is inside the ±3 % band, so nothing is re-recorded
(D11: a re-record is for a move beyond the band). On `m33` the increase is
+1,661 instructions per 512-point transform for Q15 (6,804,517 over the
scenario's 4,096 transforms, 2,048 per direction) and +1,206 for Q31
(4,940,790 over 4,096); `rfft_q31_2048` +4,173 (4,272,748 over 1,024). The
old and new loops execute the same operations per bin pair (eight
saturating `add` / `sub`, four `mul_coeff`, as the identical outputs
require), so the operations are not the difference; the code GCC generates
around them is. Attributed on `m33` from per-instruction execution counts
(a per-PC counting QEMU plugin on the icount binaries built from `main` at
`0eb09fa` and from this tree with the bench's `-O3` flags; local,
2026-09-26, totals within 0.001 % of the ratchet's):

- the post-pass and pre-pass loop bodies (the instructions `main()`
  executes 260,096 = 2,048 × 127 times, the two loops together) grow from
  351 to 375 instructions for Q31 and from 353 to 378 for Q15, +12 per bin
  pair per transform, +6,242,304 / +6,502,400 over the scenario;
- the disassembly diff of those bodies: the int64 clamps are the same
  compare-and-subtract pairs, but the old loops branched on most of them
  to an out-of-line clamp block that never runs (three instructions
  executed), where the new loops if-convert eight more into an IT block of
  two predicated moves (five executed): IT-predicated clamps 16 → 24,
  branching clamps 27 → 20 (Q31; Q15 25 → 20), plus the register copies the
  predicated form needs (plain `mov` 17 → 28);
- the rest nets against it: one other loop in `main()` (where the bench
  harness and the inlined transform sit together) is one instruction
  shorter for Q31 (executed 1,048,576 times) and one longer for Q15
  (524,288 times), and outside `main()` the table build at construction is
  236,562 (Q31) / 236,528 (Q15) instructions cheaper (libm and the
  soft-double helpers run on the first octant only); the kernel's radix-4
  stages execute the same count to the instruction.

`rotate` inlines at `-O3`; its out-of-line call at `-Os` (the test legs'
MinSizeRel) is not in what the ratchet counts. A
variant that loads the
four inputs into locals first measured worse on `m33` (+1.24 % /
+1.32 % / +0.94 %) and was not kept.

**Stage 2c (this PR): the fixed-point scenarios seeded.** Same run as the
sizes above (35861317022 at `21b3488`), compare mode on every key: the float
scenarios at +0.00 % against the 2b baselines on all five keys (the C left the
shipping tree — it remained under `tests/reference/ooura/` as the parity oracle
until D6 deleted it — and the shipping path did not move), the three
fixed-point scenarios recorded as new baselines (`bench/README.md`, one row
per key and scenario). Executed guest instructions for the whole scenario
binary, as for the float rows: 2^20 samples per direction through
`forward()` + `inverse()` with the exponent folded after each block, construction and the
checksum.

| Scenario | `m4-softfp` | `m4f` | `m33` | `m55` | `m55-ooura` |
|---|---|---|---|---|---|
| `rfft_q15_512` | 746,311,273 | 753,215,637 | 752,189,619 | 684,157,511 | 684,157,511 |
| `rfft_q31_512` | 709,161,574 | 715,019,447 | 714,626,257 | 657,595,969 | 657,595,969 |
| `rfft_q31_2048` | 869,187,566 | 876,486,909 | 875,801,377 | 806,142,635 | 806,142,635 |
| `rfft_q31_512` / `rfft_f32_512` | 0.38 | 7.36 | 7.08 | 12.55 | 7.37 |

Read against the float engine: the int32 kernel costs 7.1-7.4x the split-radix
float count on the FPU cores (M4F, M33, M55 with the engine) and 12.6x the
CMSIS-DSP Helium count on `m55`, but only 0.38x the soft-float M4's float
count — the profile the fixed-point kernel exists for. Q15 costs 4-5 % more
than Q31 at the same N (the widen/narrow passes); the 2048-point Q31
transform costs 1.23x the 512-point one for the same sample throughput
(log2 N stages: 11 vs 9, plus the post-pass). These are the numbers every
later `SMMULR`, Helium or table-layout change to `fft/fixed_point.h` is
ratcheted against (D11). Host reference (x86-64, GCC 13.3 `-O3`, callgrind):
434 M / 407 M / 495 M instructions for `rfft_q15_512` / `rfft_q31_512` /
`rfft_q31_2048`, of which the two radix-4 stage kernels are 71-80 %.


**Stage 2b (the flip; tap/DspTap#31): the baselines re-recorded to the
port.** Same run as the sizes above (35844483811 at `b078a80`, seed mode on
the four split-radix keys, compare mode on `m55`). "C" is the seeded
baseline (`df482d1`, the vendored C through `basic_real_fft`); "port" is the
new baseline, `basic_real_fft` routed at the engine; the `_c` sibling of the
same run reproduced the seeded C counts to within +28 … +34 instructions
(the adapter's prologue) with output checksums identical to the port's on
every Ooura key. `bench/README.md` carries the per-row record.

| Scenario | `m4-softfp` | `m4f` | `m33` | `m55` (CMSIS) | `m55-ooura` |
|---|---|---|---|---|---|
| `rfft_f32_512`, C (seeded `df482d1`) | 1,868,441,244 | 98,090,666 | 102,248,169 | 52,382,331 | 94,561,954 |
| `rfft_f32_512`, port (re-recorded, run 35844483811) | 1,864,929,141 | 97,138,544 | 100,945,841 | 52,382,331 (not re-recorded; measured 52,382,366, +0.00 %) | 89,276,321 |
| ratio port / C | 0.9981 | 0.9903 | 0.9873 | 1.0000 | 0.9441 |
| `rfft_f32_2048`, C (seeded `df482d1`) | 2,296,984,479 | 111,859,257 | 116,385,409 | 54,858,120 | 107,806,480 |
| `rfft_f32_2048`, port (re-recorded, run 35844483811) | 2,294,360,259 | 111,278,416 | 115,465,626 | 54,858,120 (not re-recorded; measured 54,858,195, +0.00 %) | 102,784,096 |
| ratio port / C | 0.9989 | 0.9948 | 0.9921 | 1.0000 | 0.9534 |
| checksums, port vs `_c` | identical | identical | identical | differ (CMSIS ≠ Ooura, expected) | identical |

Against the ±3 % band: the port executes 0.1 – 1.3 % fewer instructions than
the C on the M4 and M33 keys (inside the band) and 4.7 – 5.6 % fewer on
`m55-ooura` (IMPROVED beyond the band, hence the re-record per D11 rather
than an absorb). The port's counts through `basic_real_fft` differ from the
Stage 2a `_port` adapter's by 0.00 – 0.26 % (`m55-ooura` 512: 89,047,005 →
89,276,321; the wrapper's own prologue and the DONE line), which is why the
bare keys were re-recorded from a run of
the routed class rather than copied from the 2a table. `m55` (CMSIS) is
untouched by the flip: +35 / +75 instructions out of 52 / 55 million, the
wrapper's two dropped vectors and three more characters in the printed
engine name, +0.00 % against its baseline, which stays.

**Stage 2a (for the record).**

Executed guest instructions for the whole scenario binary (2^20 samples per
direction through `forward()` + `inverse()`, the out-of-place surface with
its copy loop and 2/N scaling, plus construction and the checksum fold;
`bench/README.md`). The C rows are the seeded baselines (`bench/baselines.json`,
seeded at `df482d1`, run 35281280300) and were reproduced at +0.00 % on the
head of tap/DspTap#28. The `_port` rows are the Stage 2a port measured in
the same run (35295262684, 2026-09-18, arm-none-eabi-gcc 13.2.1, QEMU 8.2.2),
informational: nothing is routed at the port until 2b. "ratio" is
port / C on that key; "checksums" says whether the two binaries' FNV-1a
output fingerprints agree, i.e. whether the port is bit-identical to the C
at the bench's default flags (Release, i.e. CMake's GNU default `-O3 -DNDEBUG` —
no toolchain file overrides it — with VFMA on the M4F / M33 / M55).

| Scenario | `m4-softfp` | `m4f` | `m33` | `m55` (C = CMSIS) | `m55-ooura` | Source |
|---|---|---|---|---|---|---|
| `rfft_f32_512` (C) | 1,868,441,244 | 98,090,666 | 102,248,169 | 52,382,331 | 94,561,954 | seeded `df482d1`; +0.00 % at `8afe6cf` |
| `rfft_f32_512_port` | 1,864,904,581 | 97,126,274 | 100,833,265 | 89,407,453 | 89,047,005 | `8afe6cf`, run 35295262684 |
| ratio port / C | 0.9981 | 0.9902 | 0.9862 | 1.7068 (vs CMSIS) | 0.9417 | |
| checksums | identical | identical | identical | differ (CMSIS ≠ Ooura, expected) | identical | |
| `rfft_f32_2048` (C) | 2,296,984,479 | 111,859,257 | 116,385,409 | 54,858,120 | 107,806,480 | seeded `df482d1`; +0.00 % at `8afe6cf` |
| `rfft_f32_2048_port` | 2,294,343,851 | 111,270,714 | 115,396,527 | 102,870,129 | 102,644,849 | `8afe6cf`, run 35295262684 |
| ratio port / C | 0.9989 | 0.9947 | 0.9915 | 1.8752 (vs CMSIS) | 0.9521 | |
| checksums | identical | identical | identical | differ (CMSIS ≠ Ooura, expected) | identical | |
| `rfft_f64_512` (host-class only) | n/a | n/a | n/a | n/a | n/a | |
| `rfft_q15_512`, `rfft_q31_512`, `rfft_q31_2048` | seeded at Stage 2c: the table at the top of this section | | | | | |

Read against the ±3 % ratchet 2b will apply: the port executes 0.1 % to
5.8 % *fewer* instructions than the C on every Ooura key, well inside the
band on the low side (the two-sided gate would flag an improvement beyond
3 % on `m55-ooura`, 0.9417 / 0.9521, so 2b re-records rather than absorbs
it; `m33` at 0.9862 / 0.9915 and the M4 keys are inside). Against CMSIS-DSP
Helium on the deployed `m55` profile the port is 1.7 – 1.9× the count (and
the C itself 1.81× / 1.97×: the 2b seeding run 35844483811 counts the `m55`
key's `_c` pair at 94,561,988 / 52,382,366 and 107,806,511 / 54,858,195,
C / CMSIS, which is the figure the CMSIS wrapper's comment and the README
quote in place of the consumers' transform-only "~3×"); `m55` stays on
CMSIS and the port is its fallback, as today. The port's
count on the `m55` key (89,407,453) differs from the same port binary's
count on `m55-ooura` (89,047,005) by 0.4 %: the two builds differ only in
what `tap::dsp` links, so this is startup and layout, not the transform.

### Host microbenchmark (`bench/bench_fft.cpp`, informational)

The srdif engine's host timings, against the port's, are in "The floating
engine (srdif)", "Host timings (informational)" (tap/DspTap#42, a same-machine
A/B on x86-64: at `-O3` float 0.5–5.9 % slower and double −3.4 … +3.6 %; at
`-O3 -march=native` float 9–18 % faster and double −7.5 … +3.5 %). What
follows is the
port against the vendored C.

Measured at Stage 2b (2026-09-23) on the development container: Intel Xeon
@ 2.10 GHz (4 vCPUs, a shared cloud VM), Ubuntu 24.04, g++ 13.3.0, Release
(`-O3 -DNDEBUG`, no `-march`), `tap_dsp_bench_fft` min of 25 reps of 2^20
samples per direction; `basic_real_fft` (the routed port) and `reference_c`
(the vendored C through the bench adapter) built in two configure trees and
run alternately, twice each; 1-minute load average 2.7 – 6.7 (other builds
were finishing on the machine), so only the port/C ratio is worth reading,
and the numbers are one machine's:

| Scenario | port, forward / inverse ns | C, forward / inverse ns | port / C, round trip |
|---|---|---|---|
| `rfft_f32_512` | 1482 – 1483 / 1523 – 1526 | 1510 – 1512 / 1598 – 1605 | 0.966 – 0.968 |
| `rfft_f32_2048` | 7134 – 7296 / 7219 – 7301 | 7158 – 7160 / 7502 – 7513 | 0.984 – 0.990 |
| `rfft_f64_512` | 1438 – 1449 / 1628 – 1631 | 1434 – 1435 / 1634 – 1656 | 0.992 – 1.004 |

The port is at parity with the C on x86-64 without `-march` (the inverse
is a few percent faster, the double forward a few tenths of a percent
slower: the same statements, one compiler, differently inlined), which is
what the instruction counts on the Cortex-M keys also say. Both binaries
print the same output checksums (`f32=-4753.16699`, `f64=-1903.0533955268113`).
`TODO(after stage 4, needs a Mac)`: portable engine vs vDSP same-binary
*microbenchmark* on Apple Silicon. Stage 4 made the same-binary comparison
possible (`basic_real_fft<float, detail::srdif_rdft<float>>` since #42,
`detail::split_radix_rdft<float>` before, beside the vDSP default in one
binary, which is how `test_fft_backend.cpp` now runs
parity on the macOS leg), but the number itself is a wall-clock measurement
on a Mac that this repo's linux development environment cannot make. Of the
consumers' "~3× faster / ~3× fewer instructions vs autovectorized Ooura"
claims, the M55 one is replaced in `fft.h` and the README by the bench's own
whole-scenario figure (C / CMSIS 1.81× / 1.97× on the `m55` key, above); the
Apple one is stated there as what it is, MuTap's transform-only measurement
on the vendored C (tap/MuTap#31), and is not re-measured in this repo until
someone runs `tap_dsp_bench_fft` twice on a Mac (`-DTAP_DSP_FFT_ACCELERATE`
ON and OFF) and records the machine.

## The floating engine (srdif)

`include/tap/dsp/fft/srdif.h`, `detail::srdif_rdft<Sample>`: the engine
`basic_real_fft<double>` runs on every build and `basic_real_fft<float>` runs
wherever no accelerated backend is selected. It replaced, at the same
contract (tap/DspTap#42), the port of Ooura's `rdft` (`fft/split_radix.h`)
that the floating profiles ran from Stage 2b (#31). It was written under a
clean-room procedure: from the published literature the header cites and
from DspTap's own fixed-point engine (`fft/fixed_point.h`, `fft/tables.h`),
without reference to the engine it replaced, to the package that engine
came from, or to any other FFT library ("Provenance and licensing", "The
floating engine, replaced clean-room", records the procedure and what it
did not cover). The contract table and the fp-contraction policy above state
the tree after #42; the Stage 4 record and the size and instruction-count
tables keep their dated measurements of the port, with the #42 values added
beside them.

### Structure, and why

The maintainer's decision was to generalize the fixed-point engine's
structure to floating point, and the engine does exactly that: the N real
samples as M = N/2 complex values, a decimation-in-frequency complex kernel
of length M, the bit-reversal permutation, and the real post-pass /
inverse pre-pass of Cooley, Lewis and Welch (1970) and Sorensen, Jones,
Heideman and Burrus (1987) in the arrangement tap/DspTap#39 derived for the
fixed-point engine (`fixed_point_rdft::real_post_pass`), with the scaling
removed: per bin pair (k, M − k), `G = C_k (u − v)`, `X[k] = u − G`,
`X[M − k] = conj(v + G)` with `C_k = (1 + i W_N^k) / 2`; `X[0], X[M] =
Re Z[0] ± Im Z[0]`; bin N/4 untouched; the inverse the same statements with
`conj C_k`, and its DC / Nyquist pair halved.

The complex kernel is split-radix rather than the fixed-point engine's
radix-4, for the operation count. On the soft-float Cortex-M4 every
floating operation is a library call of 40–60 instructions and the count is
the cost; split-radix (Duhamel and Hollmann 1984; the decimation-in-frequency
program and its count in Sorensen, Heideman and Burrus 1986) needs
4M log₂M − 6M + 8 real operations for the kernel, and the engine attains that
count exactly (measured through an instrumented sample type): with the
post-pass's 12 per bin pair, 2N log₂N − 2N − 2 per forward transform
(8,190 at N = 512, 40,958 at 2048) and two more per inverse.

On the FPU cores the count is loads and stores as much as arithmetic, and
the kernel is arranged around that:

- **The fused two-level pass.** The split-radix step on a block of length l
  computes, for j < l/4, the butterfly on x[j], x[j + l/4], x[j + l/2],
  x[j + 3l/4]; its two sums are the inputs of the half-length block. The
  half-length block's own butterfly at j (quarter l/8) reads x[j],
  x[j + l/8], x[j + l/4], x[j + 3l/8] — exactly the half-block outputs of the
  level-l butterflies at j and j + l/8. So one group of eight values
  x[j + r·l/8], r = 0 … 7, runs two levels between one load and one store:
  three butterflies per 16 loads and 16 stores instead of 48 of each. Each
  block of l ≥ 32 is processed that way and recurses on the five blocks the
  two levels leave (l/4, l/8, l/8, l/4, l/4). This is the index map of the
  split-radix decomposition (Duhamel and Vetterli 1990 give the general
  form), read two levels at a time.
- **Special butterflies** where the twiddle is trivial: j = 0 (additions
  only) and j = l/8 (the eighth turn, two multiplications per product
  instead of four), the two cases Sorensen et al. count; inside a group the
  level-l/2 butterfly at j = l/16 is its own eighth-turn case and the level-l
  twiddles there are W₁₆¹, W₁₆³, W₁₆³ and W₁₆⁹, spelled as literals.
- **Twiddle pairing.** W^(l/4 − j) = i·conj(W^j) and W^(3(l/4 − j)) =
  −i·conj(W^(3j)): the group l/8 − j uses the six twiddles of the group j
  with cos and sin exchanged and signs flipped (exact), so the table holds
  j < l/16 only and the second half of each pass walks it backwards.
- **Compile-time blocks and register leaves.** Blocks of 32 and 64 run
  with every offset a constant (one pointer, immediate addressing, fully
  unrolled groups) from their own small tables; blocks of 16 and fewer are
  leaves held in registers through every remaining level. The run-time
  recursion handles l ≥ 128 only, and its groups run two per loop step so
  the eight row pointers serve both. (As #42 shipped it. The Hexagon tuning
  below replaced the register leaves with one split-radix level at a time in
  memory down to pairs, and runs the double profile's run-time groups one
  per step.)
- **The permutation** keeps no index table: a bit-reversed counter walks
  beside the natural index (the reverse-carry increment of Gold and Rader,
  *Digital Processing of Signals*, 1969: adding one flips an index's
  trailing ones and the zero above them, so the reversed counter flips the
  mirror image of those bits), and splitting the index into its two end bit
  pairs and its middle (one of the arrangements Karp's 1996 survey
  describes) lets one counter step over the middle serve sixteen indices:
  with M = 2^b and i = A·M/4 + 4y + B (A, B < 4), rev(i) = rev₂(B)·M/4 +
  4·rev′(y) + rev₂(A), so the six (A, B) with A < rev₂(B) swap whatever y
  is, the four with A = rev₂(B) swap when y < rev′(y), and the other six are
  their partners. Below M = 16 the counter walks every index. Until the
  fix pass for review A of #42 this was a precomputed list of swap pairs
  (M uint32 words of heap); the counter costs fewer instructions on every
  key (below) and, alone, is 1.7–2.4× faster on x86-64 from N = 65536 up
  (the list, 2 bytes per sample, streamed through the cache beside the
  data).
- **No SLP vectorization under GCC.** A `#pragma GCC optimize
  ("no-tree-slp-vectorize")` around the class: GCC's basic-block
  vectorizer packs the (re, im) halves into vector registers and pays a
  lane shuffle for every multiplication by i and every complex product.
  That was the host double slowdown review A of #42 measured (double forward
  at N = 512: 1,416 ns with it, 1,228 ns without, same session, `-O3`; 1,438
  and 1,126 ns at `-march=native`) and 1.9–2.3 % of the M55 float scenario's
  instructions; the M4, M4F and M33 have no vector unit and compile the same
  code either way. Clang's SLP vectorizer helps on this code and is left on.
  Nothing reassociates, so the output bits are the same with and without
  (the fingerprints match).
- **Out-of-line transforms.** `forward_inplace` / `inverse_inplace` are
  `noinline` (a performance attribute only): inlined into a caller's loop,
  the permutation and the post-pass lost their registers to the caller's
  live values (+1.4 % on the M4F scenario, measured).

### The tables: integer trigonometry, no libm

Every table value comes from `srdif_trig`: cos, sin, 1 − cos and 1 − sin of
2πk/2^L for the first octant, in unsigned 64-bit fixed point — θ as the
exact 128-bit product of k and ⌊2π·2⁶¹⌋, Horner's rule over the Maclaurin
coefficients ⌊2⁶⁴/n!⌋ (cos to 20!, sin to 19!), 1 − cos as θ²·(1 − cos)/θ² so
small angles keep their relative precision — rounded once to the profile
(ties to even on the 64-bit mantissa), and taken to the other octants by
exact symmetry. Against libquadmath, exhaustive over every k of every
2^L, L = 3 … 30 (every angle any N ≤ 2^30 builds; review A of #42 ran
it, 2.7·10⁸ angles × four quantities): the mantissa is within 13.44 units
of 2⁻⁶⁴ relative (1 − cos at L = 30; 12.34 at L = 19, 13.07 at L = 24);
**every float entry is the correctly rounded value**, with no misrounding
at any L ≤ 30 (the nearest exact value to a float tie is 547 such units
away, first at L = 25; 21,970 up to L = 24); every double entry is within
0.5 + 13.44/2048 = 0.5 + 2^−7.25 ulp (0.5052 ulp measured), and 0.04–0.17 %
of them, by function, are the neighbour of the correctly rounded double
(0.08 % of the 2^20 post table). Only the post table is evaluated: the
kernel and block twiddles are read out of its imaginary parts (cos θₖ/2 for
every k < N/4, doubled exactly), which gives the bits a table evaluated at
their own resolution would hold, since `srdif_trig` computes the same thing
at (k, L) and (2k, L + 1).

The point is host identity. A table from `std::cos` / `std::sin` is only as
portable as the last bit of each libm, and the engine this one replaced
needed its output pins per C library and per glibc CPU dispatch. Here the
table is a function of (k, L), the arithmetic is IEEE operations in a fixed
order, and with fp-contraction off the output bits are a function of the
source alone: `tests/test_fft_srdif_fingerprint.cpp` pins FNV-1a-64 of the
forward and inverse output at every power of two 4 … 65536, both profiles,
and **one row** matched every compiler and target CI runs: x86-64 Linux
(g++ 13.3 and clang++ 18.1, glibc 2.39 under either libm dispatch), all
four QEMU legs (soft-float M4, M4F, M33, M55, newlib; N ≤ 4096 there), and,
in CI, Windows x64 (MSVC, no contraction by default) and macOS arm64
(AppleClang at `-ffp-contract=off`, the srdif engine named explicitly);
review A of #42 read those logs, and compared the table bytes themselves
across g++ `-O0` / `-O3` / `-march=native`, clang++ and the four legs. CI
runs the linux check twice, the second time under glibc's SSE2 dispatch, as
the evidence that no libm reaches the output. A mismatch is a finding, not
a row to add. The claim's scope is those builds: the test asserts
`FLT_EVAL_METHOD == 0` (an x87 build evaluates in extended precision), the
no-contraction flags cover GNU, Clang, AppleClang and IntelLLVM (the last
also at `-fp-model=precise`, since icx defaults to `-fp-model=fast`; no CI
leg builds with icx, so that flag is untested), MSVC
on ARM64 has no leg and is unverified, and a consumer's `-ffast-math` or
flush-to-zero mode is outside it. A side effect on the Cortex-M legs: the float
profile links no libm `sin` / `cos` at all, and the MinSizeRel float probe's
`.text` fell by 11–22 kB (below).

Construction evaluates N/8 angles (the post-pass table's), each a few
hundred instructions of 32 × 32 → 64 multiplies on a Cortex-M4: 52 k
instructions for the whole constructor at N = 512 and 204 k at 2048 on the
M4F (103 k and 396 k before the fix pass, which evaluated the kernel
table's octant separately and built the swap list). The size precondition
is checked (`assert`) before any table is sized from N. At N = 2^30 (float,
x86-64, the engine alone) construction takes 12.8 s at a 3.67 GB peak,
against 39.3 s and 5.77 GB before (the swap list was built bit by bit,
O(M log M)).

### Memory

One heap allocation per object, sized once and built in place (no
temporary at construction): N/2 Samples of post-pass coefficients, the
compile-time blocks' twiddles (12 Samples at N = 64, 48 from N = 128) and
the run-time pass's (12·(N/32 − 1) Samples from N = 256), i.e.

  heap = sizeof(Sample) × (N/2 for N ≤ 32; 44 at N = 64; 112 at N = 128;
  7N/8 + 36 from N = 256) bytes

(`srdif_rdft::heap_bytes`; fft.h states it beside the size range, and
`fft_rt_srdif.*HeapIsOneAllocationOfTheStatedSize` pins the count and the
bytes). `sizeof` is 32 B for the engine and 40 B for `basic_real_fft`
on LP64 (72 B before #42, 112 B before the fix pass). Measured (x86-64,
libstdc++):

| N | float heap | double heap | engine it replaced (float / double) | ratio (float / double) |
|---:|---:|---:|---:|---:|
| 512 | 1,936 B | 3,872 B | 1,100 / 2,124 B | 1.76 / 1.82 |
| 2048 | 7,312 B | 14,624 B | 4,236 / 8,332 B | 1.73 / 1.76 |
| 65536 | 229,520 B | 459,040 B | 131,808 / 262,880 B | 1.74 / 1.75 |

1.73–1.82 times the predecessor, the price of the pre-rounded post-pass
coefficients and the fused pass's twelve-value entries (whose values are
also in the post table; reading them from there would save 3N/8 Samples
and cost the pass its contiguous entries). Before the fix pass for review A
of #42 the object also held the swap list (M uint32 words) and made nine
allocations, five of them temporaries interleaved with the four kept
tables: 2,896 / 4,832 B at N = 512, 2.2 (double) to 2.8 (float) times the
predecessor.

### Accuracy against a quad-precision reference

The maintainer's targets sheet's method, with its harness (a program over
the public `basic_real_fft` API only, kept outside the tree with the sheet):
`__float128` reference (O(N²) DFT for N ≤ 8192, an independent radix-2
FFT above), full-scale input rounded once to the profile; *rms* =
‖y − ref‖₂ / ‖ref‖₂ over the N packed words, *max* = max|y − ref| / max|ref|;
white noise pooled over eight trials (seeds 0x5EED0001 … 8), max over the
worst trial; forward, unnormalized inverse of the exact spectrum rounded to
the profile, and the round trip through `inverse()`. x86-64, g++ 13.3.0
`-std=gnu++20 -O3 -DNDEBUG`, no `-march` (no FMA), 2026-09-26. The last
two columns are this engine over its predecessor on the same trials.

White noise, float:

| N | ref | fwd rms | fwd max | inv rms | inv max | rt rms | rt max | rms ratio fwd/inv/rt | max ratio fwd/inv/rt |
|---:|:---:|---:|---:|---:|---:|---:|---:|---|---|
| 4 | DFT | 4.16e-08 | 7.71e-08 | 3.11e-08 | 6.11e-08 | 4.79e-08 | 7.47e-08 | 1.00/0.92/0.93 | 1.00/1.00/1.00 |
| 8 | DFT | 5.47e-08 | 7.63e-08 | 4.57e-08 | 7.32e-08 | 7.59e-08 | 1.39e-07 | 1.00/0.94/1.01 | 1.00/1.00/1.00 |
| 16 | DFT | 6.82e-08 | 1.09e-07 | 6.98e-08 | 1.10e-07 | 8.81e-08 | 1.57e-07 | 1.00/0.97/1.02 | 1.00/1.00/1.00 |
| 32 | DFT | 7.05e-08 | 1.10e-07 | 7.61e-08 | 1.62e-07 | 1.04e-07 | 1.96e-07 | 0.97/0.97/1.00 | 1.00/1.14/0.95 |
| 64 | DFT | 8.75e-08 | 1.36e-07 | 8.00e-08 | 1.36e-07 | 1.20e-07 | 3.01e-07 | 0.94/0.90/1.04 | 0.90/0.68/1.47 |
| 128 | DFT | 9.08e-08 | 1.44e-07 | 8.87e-08 | 2.01e-07 | 1.20e-07 | 2.99e-07 | 0.93/0.95/0.90 | 1.06/0.83/1.00 |
| 256 | DFT | 9.71e-08 | 1.44e-07 | 9.71e-08 | 2.17e-07 | 1.34e-07 | 2.98e-07 | 0.91/0.89/0.86 | 0.78/0.93/0.70 |
| 512 | DFT | 1.06e-07 | 1.74e-07 | 1.02e-07 | 2.04e-07 | 1.44e-07 | 3.58e-07 | 0.92/0.89/0.90 | 1.17/0.81/1.00 |
| 1024 | DFT | 1.12e-07 | 1.52e-07 | 1.10e-07 | 2.59e-07 | 1.54e-07 | 3.58e-07 | 0.93/0.94/0.94 | 1.03/0.86/1.00 |
| 2048 | DFT | 1.17e-07 | 1.87e-07 | 1.16e-07 | 2.75e-07 | 1.62e-07 | 4.17e-07 | 0.91/0.91/0.91 | 1.06/0.80/0.87 |
| 4096 | DFT | 1.23e-07 | 1.58e-07 | 1.22e-07 | 2.96e-07 | 1.71e-07 | 4.17e-07 | 0.90/0.90/0.88 | 0.62/0.79/0.82 |
| 8192 | DFT | 1.29e-07 | 2.02e-07 | 1.28e-07 | 3.29e-07 | 1.80e-07 | 4.77e-07 | 0.94/0.93/0.95 | 1.17/0.90/1.00 |
| 16384 | FFT | 1.33e-07 | 1.77e-07 | 1.33e-07 | 3.55e-07 | 1.86e-07 | 5.36e-07 | 0.91/0.91/0.90 | 0.85/0.92/0.90 |
| 32768 | FFT | 1.38e-07 | 1.78e-07 | 1.38e-07 | 3.87e-07 | 1.94e-07 | 5.36e-07 | 0.92/0.92/0.93 | 0.87/0.91/0.82 |
| 65536 | FFT | 1.43e-07 | 1.81e-07 | 1.42e-07 | 4.28e-07 | 2.00e-07 | 6.56e-07 | 0.92/0.92/0.92 | 0.90/0.88/1.10 |
| 2^20 | FFT | 1.61e-07 | 1.99e-07 | 1.60e-07 | 5.09e-07 | 2.25e-07 | 7.15e-07 | 0.93/0.93/0.94 | 0.93/0.95/0.92 |

White noise, double:

| N | ref | fwd rms | fwd max | inv rms | inv max | rt rms | rt max | rms ratio fwd/inv/rt | max ratio fwd/inv/rt |
|---:|:---:|---:|---:|---:|---:|---:|---:|---|---|
| 4 | DFT | 5.54e-17 | 8.41e-17 | 5.32e-17 | 9.73e-17 | 8.50e-17 | 1.30e-16 | 1.00/1.00/1.00 | 1.00/1.00/1.00 |
| 8 | DFT | 1.02e-16 | 1.49e-16 | 1.01e-16 | 1.65e-16 | 9.99e-17 | 1.52e-16 | 0.98/0.95/0.71 | 1.00/1.00/0.67 |
| 16 | DFT | 1.30e-16 | 2.09e-16 | 1.28e-16 | 2.74e-16 | 1.59e-16 | 3.04e-16 | 0.94/0.97/0.96 | 0.81/0.83/0.90 |
| 32 | DFT | 1.36e-16 | 1.89e-16 | 1.20e-16 | 2.19e-16 | 1.84e-16 | 3.34e-16 | 1.00/0.89/0.99 | 1.00/0.76/0.86 |
| 64 | DFT | 1.52e-16 | 2.28e-16 | 1.39e-16 | 3.72e-16 | 2.03e-16 | 3.41e-16 | 0.93/0.88/0.83 | 1.00/1.17/0.61 |
| 128 | DFT | 1.68e-16 | 2.65e-16 | 1.68e-16 | 3.53e-16 | 2.19e-16 | 5.62e-16 | 0.93/0.96/0.92 | 1.00/0.99/1.26 |
| 256 | DFT | 1.81e-16 | 2.47e-16 | 1.82e-16 | 4.50e-16 | 2.41e-16 | 5.64e-16 | 0.92/0.93/0.89 | 0.94/1.06/1.01 |
| 512 | DFT | 1.89e-16 | 2.26e-16 | 1.89e-16 | 4.46e-16 | 2.61e-16 | 5.70e-16 | 0.91/0.92/0.90 | 0.95/0.78/0.85 |
| 1024 | DFT | 2.06e-16 | 3.16e-16 | 2.04e-16 | 5.08e-16 | 2.81e-16 | 6.67e-16 | 0.93/0.92/0.93 | 1.00/1.09/0.92 |
| 2048 | DFT | 2.18e-16 | 2.84e-16 | 2.17e-16 | 5.46e-16 | 3.00e-16 | 7.77e-16 | 0.91/0.92/0.90 | 0.83/0.95/0.78 |
| 4096 | DFT | 2.28e-16 | 3.04e-16 | 2.25e-16 | 5.49e-16 | 3.17e-16 | 7.78e-16 | 0.92/0.91/0.93 | 0.97/0.83/0.88 |
| 8192 | DFT | 2.39e-16 | 3.17e-16 | 2.35e-16 | 6.63e-16 | 3.30e-16 | 8.88e-16 | 0.92/0.92/0.93 | 0.91/0.86/1.00 |
| 16384 | FFT | 2.47e-16 | 3.12e-16 | 2.47e-16 | 6.60e-16 | 3.43e-16 | 9.99e-16 | 0.91/0.91/0.92 | 0.86/0.88/1.00 |
| 32768 | FFT | 2.57e-16 | 3.04e-16 | 2.54e-16 | 7.83e-16 | 3.56e-16 | 1.11e-15 | 0.91/0.90/0.91 | 0.88/0.98/0.91 |
| 65536 | FFT | 2.66e-16 | 3.04e-16 | 2.65e-16 | 7.72e-16 | 3.69e-16 | 9.99e-16 | 0.90/0.91/0.91 | 0.77/0.91/0.90 |
| 2^20 | FFT | 2.99e-16 | 3.59e-16 | 2.98e-16 | 9.00e-16 | 4.17e-16 | 1.33e-15 | 0.91/0.91/0.92 | 0.97/0.80/0.92 |

The rms targets (no worse than the predecessor within 5 % measurement
noise) are met at every N in both profiles and every direction: 4–14 %
lower from N = 128 up, and at most 4 % higher below that (float round trip
at N = 64). The max targets (1.25× the predecessor) are met in every noise
cell but two: the float round trip at N = 64 (3.01e-7 vs 2.04e-7, 1.47×)
and the double round trip at N = 128 (5.62e-16 vs 4.46e-16, 1.26×). Both
are the maximum over one draw of eight trials, and review A of #42
settled them with a paired A/B (its measurement, not this note's: the
targets harness run on both trees, the same inputs to both engines, a
`__float128` radix-2 reference; 400 trials at N = 64, 128, 256, 512 and
4096, 60 at 65536 (240 for the off-bin tone); noise, random-frequency on-
and off-bin tones, random-position impulses). Every (profile, N, signal,
direction) cell has a new/old mean-error ratio ≤ 1.00 for rms, and ≤ 1.04
for max (0.98 after 240 more trials at 65536); the srdif engine is worse on
0–35 % of paired trials in almost every cell, 52 % at most (the float
impulse round trip at N = 64, equal means). The two misses are draws: the
float round-trip max at N = 64 has mean ratio 0.93 (p95 0.89), the double
one at N = 128 0.95 (p95 0.86).

Deterministic signals (on-bin tone, off-bin tone, impulse at x[1],
constant; one realization each): the constant is bit-exact against the
reference (all errors 0) in both profiles at every N, as it was. Of the
other three, 516 cells (N = 4 … 65536 × two profiles × three directions ×
rms / max) with a nonzero predecessor value: the srdif engine is lower in
310 of them and above the targets' bars (1.05 × rms, 1.25 × max) in 54:
21 for the on-bin tone (the double forward at N = 128 … 512 is the largest
rms excess, 1.85e-16 vs 1.10–1.18e-16, where the predecessor happened to
cancel; at N ≥ 1024 the srdif engine is the lower one, 1.84–2.39e-16 vs
1.96–2.52e-16), 13 for the off-bin tone, and 20 for the impulse (the
round-trip maximum on the zero-valued samples, 4e-9 vs 1.1e-9 in float and
7.8e-18 vs 2.1e-18 in double at N = 65536, while the impulse round-trip rms
is 6–8 % lower from N = 16384 up). A single realization's error is a draw from the same
rounding statistics the noise rows average, which is why the targets'
headline is the pooled noise. Review A's paired sweeps say the same of the
deterministic cells: at N = 128 the on-bin tone at k0 = 13 (the sheet's
⌊N/10⌋ + 1, the normalized frequency of its N = 256 and 512 cells) is the
srdif engine's single worst of the 62 bins, 1.68× rms and 2.78× max, with a
mean over all bins of 0.93×; at N = 4096, over 2046 bins, the float forward
rms is 1.5× worse at no bin and 1.5× better at 56 (double: 2 and 34); the
impulse round-trip maximum's growth is specific to x[1] (over random
positions at 65536 the mean ratio is 0.65). Under `-march=native` (FMA
contraction) the noise and tone rms ratios are 0.83–0.96, and contraction
moves the srdif engine's own rms by −1 to −7 %.

### Instruction counts and `.text`

`scripts/icount.py` as `bench.yml` runs it (fresh Release build per key,
arm-none-eabi-gcc 13.2.1, qemu-system-arm 8.2.2, the pinned plugin
header), 2026-09-27, the tree after the fix pass for review A of #42. The
gate was "every float scenario on every key that runs this engine at or
below its baseline":

| key | `rfft_f32_512`: baseline → srdif | Δ | `rfft_f32_2048`: baseline → srdif | Δ |
|---|---:|---:|---:|---:|
| `m4-softfp` | 1,864,929,141 → 1,814,303,702 | −2.71 % | 2,294,360,259 → 2,243,212,021 | −2.23 % |
| `m4f` | 97,138,544 → 92,575,609 | −4.70 % | 111,278,416 → 106,282,752 | −4.49 % |
| `m33` | 100,945,841 → 92,979,212 | −7.89 % | 115,465,626 → 106,773,786 | −7.53 % |
| `m55-ooura` | 89,276,321 → 83,933,381 | −5.98 % | 102,784,096 → 97,602,563 | −5.04 % |
| `m55` (CMSIS, unchanged) | 52,382,331 → 52,382,329 | −0.00 % | 54,858,120 → 54,858,158 | +0.00 % |

The fix pass moved only the float scenarios of the four portable-engine
keys, all down, from the first srdif record (512 / 2048): `m4-softfp`
1,815,063,718 / 2,244,317,934 (−0.04 / −0.05 %), `m4f` 94,003,130 /
108,030,219 (−1.52 / −1.62 %), `m33` 94,382,323 / 108,470,414 (−1.49 /
−1.56 %), `m55-ooura` 86,188,216 / 100,091,803 (−2.62 / −2.49 %). Of that,
construction (half the angles, no swap list) is −0.05 / −0.18 % on the M4F
(51 k / 192 k instructions), SLP off is −1.9 … −2.3 % on the M55 and
nothing elsewhere, and the permutation without an index table the rest,
about −1.5 % on the M4F and M33 (a first version that split off one end
bit per side instead of two measured +2.0 … +2.7 % there, and one that
addressed the sixteen indices through eight row pointers −0.8 … −1.0 %;
both were replaced).

The float baselines are re-recorded to these counts (`bench/README.md`);
the fixed-point scenarios measure +0.49 … +1.15 % against theirs with no
change to that code, inside the band, and are not. Against CMSIS-DSP on the
M55 the srdif engine now executes 1.60× (N = 512) / 1.78× (N = 2048) the
instructions (the predecessor 1.70× / 1.87×).

How the engine got there (development measurements on the `m4f` key, N =
512 / 2048, each a change to the one before; the counts include the
harness's ~14.6 k instructions per forward + inverse pair):

| step | `m4f` 512 | `m4f` 2048 |
|---|---:|---:|
| first version: split-radix recursion to length 2, the post-pass, a swap-list permutation | +29.5 % | +27.4 % |
| register leaves for blocks ≤ 16 (≈ 300 calls per transform at N = 512 gone) | +13.7 % | +13.8 % |
| permutation through FP registers (a `memcpy` of the complex pair through stack temporaries had cost 29 instructions per swap) | +7.0 % | +7.8 % |
| transforms out of line | +5.4 % | +6.4 % |
| post-pass reads each pair once (the compiler re-loaded the partner after the stores, for aliasing) | +4.4 % | +5.5 % |
| the fused two-level pass | +1.6 % | +1.7 % |
| post-pass two bin pairs per step | +0.5 % | +0.8 % |
| compile-time blocks of 32 and 64 | −2.9 % | −2.1 % |
| run-time fused pass two groups per step | −2.9 % | −2.7 % |

Tried and not kept: 64-bit integer moves in the permutation (GCC splits
them into the same four accesses; no change). The soft-float key was below
its baseline from the leaves step on, because the operation count is the
split-radix minimum from the first version.

`.text` of the MinSizeRel float probe (`tap_dsp_size_probe_rfft_f32_512`,
`size -A`): `m4-softfp` 31,081 B (53,289 before), `m4f` 28,921 (44,601),
`m33` 28,305 (44,001), `m55-ooura` 28,161 (39,281); the libm double
`sin` / `cos` the predecessor's tables linked are gone (the fix pass took
another 1.1 kB, the swap-list builder and three table vectors). The four
float ceilings are re-recorded to measured + 3 % (32,064 / 29,824 / 29,184 /
29,056); the Q15 / Q31 probes and the `m55` CMSIS probe are unchanged.

### Host timings (informational)

A same-machine A/B of the targets sheet's `hostbench.cpp` (`bench/bench_fft.cpp`'s
workload over five sizes: out-of-place `forward()`, scaled `inverse()`,
batches of 2^20 samples), built `g++ 13.3 -std=c++20 -O3 -DNDEBUG` with and
without `-march=native` (Xeon 2.1 GHz with AVX-512 and FMA), pinned to one
core, 11 runs, the median of the per-run medians in ns per transform. The
predecessor's medians are review A of #42's (it ran both trees; this tree
cannot build the predecessor); this engine's were measured the same way on
the same machine the next day, alternating each run with review A's own
binary of the previous srdif tree, whose medians came out within 2 % of
review A's for it (so the two sessions compare). Predecessor → srdif:

`-O3`:

| N | float forward | float inverse | double forward | double inverse |
|---:|---:|---:|---:|---:|
| 256 | 597 → 602 (+0.8 %) | 628 → 650 (+3.5 %) | 541 → 549 (+1.5 %) | 616 → 606 (−1.6 %) |
| 512 | 1,276 → 1,349 (+5.7 %) | 1,337 → 1,416 (+5.9 %) | 1,227 → 1,228 (+0.1 %) | 1,395 → 1,348 (−3.4 %) |
| 1024 | 2,900 → 2,927 (+0.9 %) | 2,988 → 3,064 (+2.6 %) | 2,671 → 2,734 (+2.4 %) | 2,932 → 2,881 (−1.7 %) |
| 2048 | 6,160 → 6,458 (+4.9 %) | 6,390 → 6,653 (+4.1 %) | 5,980 → 6,071 (+1.5 %) | 6,540 → 6,355 (−2.8 %) |
| 4096 | 13,746 → 13,819 (+0.5 %) | 14,035 → 14,307 (+1.9 %) | 12,925 → 13,386 (+3.6 %) | 14,161 → 14,302 (+1.0 %) |

`-O3 -march=native`:

| N | float forward | float inverse | double forward | double inverse |
|---:|---:|---:|---:|---:|
| 256 | 532 → 461 (−13.4 %) | 570 → 466 (−18.1 %) | 497 → 500 (+0.4 %) | 538 → 510 (−5.3 %) |
| 512 | 1,161 → 1,039 (−10.6 %) | 1,233 → 1,053 (−14.6 %) | 1,128 → 1,126 (−0.2 %) | 1,198 → 1,138 (−5.0 %) |
| 1024 | 2,464 → 2,242 (−9.0 %) | 2,651 → 2,263 (−14.7 %) | 2,413 → 2,456 (+1.8 %) | 2,630 → 2,434 (−7.5 %) |
| 2048 | 5,808 → 5,042 (−13.2 %) | 6,111 → 5,061 (−17.2 %) | 5,265 → 5,411 (+2.8 %) | 5,755 → 5,407 (−6.0 %) |
| 4096 | 11,983 → 10,762 (−10.2 %) | 12,631 → 10,796 (−14.5 %) | 11,527 → 11,928 (+3.5 %) | 12,391 → 12,119 (−2.2 %) |

The double forward, the golden model's slowest cell, is +0.1 … +3.6 % at
`-O3` and −0.2 … +3.5 % at `-march=native`. Review A measured the previous
srdif tree at +9 … +18 % and +26 … +28 % there; the difference is GCC's SLP
vectorizer ("Structure, and why"), whose shuffles cost both directions
alike (same session, N = 512 double: forward 1,416 → 1,228 ns, inverse 1,514
→ 1,348 at `-O3`; 1,438 → 1,126 and 1,386 → 1,138 at `-march=native`). Why
it looked like a forward problem: the predecessor's double forward was
12 % cheaper than its own inverse (1,227 vs 1,395 ns at N = 512), while
this engine's two directions run the same kernel, conjugated, and cost the
same to within 10 %; the same slowdown therefore showed as twice the loss
against the forward. Float at `-O3` stays 0.5–5.9 % slower (16 SSE
registers against the M4F's 32 single-precision ones, for which the
register leaves are sized). Not a gate; recorded so no host speed claim is
made.

### Tests changed for the engine

- `tests/test_fft_srdif_fingerprint.cpp` + `tests/support/srdif_fingerprints.h`
  (new): the output fingerprints above, their own target at
  `-ffp-contract=off` (IntelLLVM also at `-fp-model=precise`), asserting
  `FLT_EVAL_METHOD == 0`; `ci.yml` runs them verbosely and requires the one
  row ("every compiler and target CI runs") under both glibc dispatches.
- `tests/test_fft_rt.cpp`: the allocation guard also sums the bytes
  requested, and `fft_rt_srdif.{Float,Double}HeapIsOneAllocationOfTheStatedSize`
  pin the heap formula at N = 4 … 4096: one allocation of exactly
  `heap_bytes(N)` per constructed engine.
- `tests/test_fft_srdif.cpp` (new): the integer trigonometry against libm to
  the last bits (1 / 2 / 4 / 6 ulp pins on 0 / 1 / 2 / 3 measured on glibc;
  1 − cos 4 on newlib), exact octant endpoints, the leaf literals equal to
  the generator's values, round-half-to-even in `to_sample`.
- `tests/test_fft_routing.cpp`: the sizes now cross the engine's code paths
  (register leaves only at N = 4 … 32, the compile-time blocks at 64 and 128,
  the run-time pass from 256 up, and the post-pass's single-pair and paired
  forms).
- `tests/test_fft_oracle.cpp`: μ is one rounding (≤ 1.02 u) for these tables;
  the Higham constant (4, i.e. 8u) is kept; worst measured ratio 0.155
  float / 0.25 double (0.17 / 0.25 before).
- Renames only (the engine type, the ABI tag `fft_srdif`, the backend string
  `srdif`) in `test_fft_engine.cpp`, `test_fft_backend.cpp`
  (`ForwardMatchesOoura` is `ForwardMatchesTheReferenceEngine`),
  `test_fft_rt.cpp`, `test_fft.cpp`, `test_log_mel.cpp`, the ABI stub and
  image, and the measured values in comments. **No tolerance and no pin
  was loosened**; every pin that the predecessor's error set is kept where
  the srdif engine measures lower (the round-trip pin of
  `ConstructsAtTheRangeBounds`, 2.86e-6 against 3.28e-7 at 2^20;
  `FloatEngineTracksDoubleAtN512`, 2.25e-7 against 9.73e-8 / 1.05e-7;
  `DoubleForwardTracksCompensatedDft`, 7.4e-16 against 1.49–1.70e-16).


### Hexagon (clang) tuning

The engine met "must not regress" on every key the ratchet measured at #42,
all four Cortex-M cores under GCC 13.2. DspTap's main consumer also ships on
Qualcomm Hexagon (clang 19, v68 with HVX-128), where the srdif engine of
`72977aa` executed more instructions than the engine it replaced; the bench
gained a `hexagon` key for it (`cmake/hexagon-linux-musl.cmake`, the
CodeLinaro clang 19.1.5 toolchain, `-mv68 -mhvx -mhvx-length=128b`, static
musl, a plugin-enabled `qemu-hexagon` 8.2.2), and the engine was tuned until
every floating Hexagon scenario is below the predecessor's count, without
costing any Cortex-M key. Done under the same clean-room rules as the
engine: from this tree, its own disassembly and counts, and the literature
the header cites.

**What the Hexagon key counts.** `qemu-hexagon` translates a packet (up to
four instructions issued together) as one guest instruction, so the plugin
counts **packets**, not instructions: a per-address profile (a scratch
plugin, one counter per translated instruction) holds counts at
packet-start addresses only, and they sum to the ratchet's total. The
Hexagon count therefore rewards what the scheduler can pack side by side,
not only what it executes. In the disassembly the float and double
arithmetic takes two of a packet's four slots at most (with the negations
and shifts that share them), the loads and stores the other two; a double
multiply is six instructions (`dfmpyfix` twice, `dfmpyll`, `dfmpylh` twice,
`dfmpyhh`), a double add one.

**Where the count went at `72977aa`** (`rfft_f32_512`, 57.66 M packets):
the class's out-of-place copies, which clang compiles to calls of musl's
`memcpy`, 17.0 M (29 %; 33.7 M of the double scenario's 110.7 M), and the
harness's fold 6.1 M, both the same for either engine (the fixed-point
scenarios, which do not run srdif, measure the same on both trees); the
engine the rest. clang's inliner had kept the butterflies' callers out of
line: the fused pass called `group<…>` per group, each register leaf went
through a stack array of values (`leaf_blocks<L>(cv*)`) and a call, and the
index-sequence lambdas of `block<32 | 64>` and `leaf<L>` were calls. GCC
inlines differently and the Cortex-M keys had not shown it.

**Measured** with `scripts/icount.py` as `bench.yml` runs it (a fresh
Release build per key), on one machine, 2026-09-28. The Hexagon counts carry
a per-process offset: user-mode qemu passes the guest its environment and the
binary's path, and the count moves with them (the review of the tuning
measured +731 packets for a 44-character longer path, and −1,584 under
`env -i`). The offset is constant across the scenarios of one setup. Here it
read +1,250 against the setup that measured the predecessor (the fixed-point
scenarios, which do not run srdif, read the same up to that offset), so the
predecessor column below is that setup's count + 1,250. The review rebuilt
`7a58ebe` at its own path and read the predecessor's float scenarios 28
packets below these cells, so read them as ±28. No percentage moves.
Against CI the offsets were +1,723 (this setup) and +1,596 (the review's).
The Cortex-M float counts reproduce `bench/baselines.json` exactly. The
fixed-point rows read the same +0.49 … +1.15 % in-band drift as on `main`:
they were seeded at #32 and have not been re-recorded since.

| key | scenario | `72977aa` | tuned | Δ | predecessor | tuned vs predecessor |
|---|---|---:|---:|---:|---:|---:|
| `hexagon` | `rfft_f32_512` | 57,660,946 | 45,936,146 | −20.33 % | 49,734,654 | −7.64 % |
| `hexagon` | `rfft_f32_2048` | 65,686,018 | 51,610,626 | −21.43 % | 54,836,437 | −5.88 % |
| `hexagon` | `rfft_f64_512` | 110,685,320 | 99,941,512 | −9.71 % | 103,241,260 | −3.20 % |
| `m4-softfp` | `rfft_f32_512` | 1,814,303,702 | 1,813,126,102 | −0.06 % | 1,864,929,141 | −2.78 % |
| `m4-softfp` | `rfft_f32_2048` | 2,243,212,021 | 2,241,935,605 | −0.06 % | 2,294,360,259 | −2.28 % |
| `m4f` | `rfft_f32_512` | 92,575,609 | 91,545,465 | −1.11 % | 97,138,544 | −5.76 % |
| `m4f` | `rfft_f32_2048` | 106,282,752 | 105,151,232 | −1.06 % | 111,278,416 | −5.51 % |
| `m33` | `rfft_f32_512` | 92,979,212 | 92,624,908 | −0.38 % | 100,945,841 | −8.24 % |
| `m33` | `rfft_f32_2048` | 106,773,786 | 106,248,986 | −0.49 % | 115,465,626 | −7.98 % |
| `m55-ooura` | `rfft_f32_512` | 83,933,381 | 83,628,229 | −0.36 % | 89,276,321 | −6.33 % |
| `m55-ooura` | `rfft_f32_2048` | 97,602,563 | 97,109,507 | −0.51 % | 102,784,096 | −5.52 % |

Unchanged, to the instruction: the `m55` key (CMSIS-DSP, 52,382,329 /
54,858,158) and every fixed-point scenario on every key, Hexagon included
(159,978,609 / 157,304,960 / 184,892,304 for Q15 512 / Q31 512 / Q31 2048).
The Cortex-M float baselines are re-recorded to these counts
(`bench/README.md`); the Hexagon key is seeded from CI. Against CMSIS-DSP
on the M55 the srdif engine now executes 1.60× (N = 512) / 1.77×
(N = 2048) the instructions (1.60× / 1.78× before).

**What paid, step by step** (Hexagon, each step on top of the one before;
`rfft_f32_512` / `rfft_f32_2048` / `rfft_f64_512`, millions of packets):

| step | f32 512 | f32 2048 | f64 512 |
|---|---:|---:|---:|
| `72977aa` | 57.66 | 65.69 | 110.69 |
| clang only: butterflies, groups, leaves, post-pass pairs and swaps forced inline (`TAP_DSP_SRDIF_INLINE`) | 49.67 | 55.24 | 105.83 |
| post-pass: two bin pairs with every load before any store (`post_two`) | 48.26 | 53.80 | 104.54 |
| the first level of a block of 16 (and of 8, 4 in double) in memory, register leaves of 8 (float) / 2 (double) | 48.15 | 53.70 | 103.02 |
| forced inline also the index-sequence lambdas and the compile-time blocks | 46.10 | 51.72 | 100.44 |
| double: the run-time fused pass one group per loop step | 46.10 | 51.72 | 100.06 |
| `mul_sub`: `a b − c d` spelled so clang fuses it as a multiply-subtract | 45.88 | 51.58 | 100.09 |
| no register leaves at all: levels in memory down to pairs, both profiles | 45.94 | 51.61 | 100.08 |
| `mul_sub`'s split spelling for float only | 45.94 | 51.61 | 99.94 |

Each change measured alone against the final tree (the change undone,
everything else kept; Hexagon f32 512 / f32 2048 / f64 512, then the
Cortex-M keys):

- **Forced inlining (clang only).** Without it: 54.13 / 62.12 / 107.12 M
  (+17.8 / +20.4 / +7.2 %). Forced under GCC as well: `m4f` +9.1 / +6.8 %,
  `m33` +5.7 / +4.5 %, `m55-ooura` +6.3 / +4.9 % (spills in the larger
  bodies), `m4-softfp` −0.22 / −0.14 %; so GCC keeps its own choices. Not
  forced under `-Os` / `-Oz` (`__OPTIMIZE_SIZE__`): forced, the Hexagon
  MinSizeRel float probe's `.text` grows from 238,628 to 263,204 bytes.
- **`post_two`.** Pair by pair the compiler cannot move the second pair's
  loads above the first pair's stores (a store through `pk` may alias a
  load through `pj` as far as it knows), so the two pairs' arithmetic
  cannot share packets. Without it: 47.47 / 53.17 / 100.97 M (+3.3 / +3.0 /
  +1.0 %). The Cortex-M keys pay 0.00 – 0.02 % for it.
- **Levels in memory down to pairs.** Register leaves of at most 16 / 8 / 4
  values restored in the final tree: Hexagon 46.10 / 45.88 / 45.89 (f32 512),
  51.79 / 51.58 / 51.55 (f32 2048), 101.30 / 100.63 / 100.12 (f64), against
  45.94 / 51.61 / 99.94 for pairs; `m4f` 92.54 / 92.00 / 91.85 and 106.30 /
  105.83 / 105.58 against 91.55 / 105.15; `m33` and `m55-ooura` 0.3 – 0.5 %
  below leaves of 8, `m4-softfp` flat. A leaf of 16 holds 32 values: the
  M4F's whole single-precision register file and all 32 of Hexagon's general
  registers, of which a double takes two. Pairs cost Hexagon float 0.1 %
  against its best (leaves of 4 or 8) and are the best everywhere else, in
  one code path.
- **One group per step in double.** Two groups per step: 100.32 M (+0.38 %).
  One per step in float too: 46.16 / 52.12 M (+0.48 / +0.99 %), so float
  keeps two.
- **`mul_sub`.** clang contracts within a statement and, on `a * b − c * d`,
  fuses the left product and negates the right, `fma(a, b, −(c d))`: a
  separate negation (`togglebit`) in the two slots the float arithmetic
  needs. `ab − c * d` with `ab` a statement of its own fuses as
  `fma(−c, d, ab)`, one `Rx −= sfmpy(Rs, Rt)`. Without it: float 46.24 /
  51.83 M (+0.65 / +0.43 %). Double keeps the plain expression (Hexagon has
  no double fused multiply-add; the split spelling only moved the schedule,
  +0.14 %). Without contraction both spellings are the same two products
  and one subtraction; GCC, which contracts after SSA across statements,
  compiles both the same (identical Cortex-M counts).

Tried and not kept:

- All eight loads of a group before its stores (the level-l butterfly at
  j + l/8 no longer waits for the stores of the one at j): +1.8 / +3.2 %
  float, +0.1 % double on Hexagon at the time; the extra live values cost
  more than the freedom.
- The permutation's six unconditional exchanges as two batches of three,
  every load before any store (measured at the `mul_sub` step): Hexagon
  float −0.39 / −0.37 %, but double +0.34 % (−0.05 % in batches of two),
  `m55-ooura` +1.69 / +1.47 % and `m4-softfp` +0.10 / +0.08 %.
- Register-leaf sizes of 16, 8 and 4 (above).

**Output bits.** Unchanged wherever the fingerprints are defined: every
change reorders loads, stores, calls and loop steps, or (`mul_sub`) splits a
statement without changing its operations, so without contraction every
output is the same operations in the same order. `tap_dsp_srdif_fingerprint`
(`-ffp-contract=off`) passes unchanged on x86-64 g++ 13.3 and clang++ 18.1,
on the four QEMU legs, and on Hexagon (clang 19.1.5 under `qemu-hexagon`,
every N from 4 to 65536, both profiles; at `72977aa` as well): **Hexagon
matches the one row**, though CI does not run it there. Under default flags,
which contract, the bench checksums move where an FMA exists: the Cortex-M
FPU keys' `rfft_f32_512` / `rfft_f32_2048` lines (`m4f`, `m33` and
`m55-ooura` alike) from `0xd26ebf9b45534325` / `0xd06aa150d5bd9325` to
`0x322a64b478d14325` / `0xb7eeda0471060d25` (GCC fuses across the former
register leaves differently once their levels run through memory; `post_two`
and `mul_sub` alone leave the GCC lines as they were), and Hexagon float
from `0x8407c06ac8ac8325` / `0xe3dffd98aa4cf725` to `0x46adbb5e3ebf2325` /
`0x7ecf6080a0f94125` (`mul_sub`). Without an FMA nothing moves: x86-64
(`0x5caf505901a97f25`, g++ and clang++), `m4-softfp`, Hexagon double. The
error statistics are the contraction policy's business ("The fp-contraction
policy"); no accuracy re-run was needed, since no bit moved at
`-ffp-contract=off`.

**Contract.** Unchanged: packing, sign, scale, sizes 4 … 2^30,
shareability, `noexcept` and allocation-free transforms, the tables and the
heap formula (nothing about the tables changed). The Hexagon test battery
(`tap_dsp_tests` and the side targets built with the Hexagon toolchain file,
run by ctest through `qemu-hexagon`, sweeps capped at 4096) passes 423 of
426; the three failures are `fft_abi_tag`'s `dlopen` tests, which a static
musl image cannot run ("Dynamic loading not supported").

**`.text`** of the MinSizeRel float probe: `m4-softfp` 31,465 B (31,081
before), `m4f` 28,985 (28,921), `m33` 28,385 (28,305), `m55-ooura` 28,289
(28,161), all under their ceilings, which are not re-recorded; Hexagon
238,628 (239,044); the Q15 / Q31 probes and the `m55` CMSIS probe
unchanged.

**Host timings** (informational; `bench/bench_fft.cpp`, `-O3 -DNDEBUG`,
g++ 13.3 and clang++ 18.1, Xeon 2.8 GHz VM (4 vCPUs, load about 1), pinned to one core, 11 runs
alternating the two trees, the median of the per-run minima, ns per
transform; `72977aa` → tuned):

| compiler | scenario | forward | inverse |
|---|---|---:|---:|
| g++ | `rfft_f32_512` | 1,714 → 1,682 (−1.9 %) | 1,788 → 1,780 (−0.5 %) |
| g++ | `rfft_f32_2048` | 8,259 → 8,022 (−2.9 %) | 8,605 → 8,398 (−2.4 %) |
| g++ | `rfft_f64_512` | 1,852 → 1,759 (−5.0 %) | 1,965 → 1,879 (−4.4 %) |
| clang++ | `rfft_f32_512` | 1,252 → 1,283 (+2.5 %) | 1,396 → 1,288 (−7.8 %) |
| clang++ | `rfft_f32_2048` | 6,136 → 6,279 (+2.3 %) | 6,753 → 6,339 (−6.1 %) |
| clang++ | `rfft_f64_512` | 1,516 → 1,404 (−7.4 %) | 1,556 → 1,449 (−6.8 %) |

Nothing is slower by more than 5 % (the clang++ float forward, +2.3 …
+2.5 %, is the only cell that rose); faster by 5 % or more: the g++ double
forward and every clang++ double and float-inverse cell. A different
machine from #42's table above, so the two tables do not compare cell by
cell.

**A finding for the class, not the engine.** On Hexagon a third of every
floating scenario is `basic_real_fft`'s out-of-place copy (`forward()` and
`inverse()` copy the input before transforming in place), which clang turns
into calls of musl's `memcpy`: 17.0 M of `rfft_f32_512`'s count, 16.8 M of
`rfft_f32_2048`'s, 33.7 M of `rfft_f64_512`'s, the same for any engine.
Consumers that transform in place do not pay it; the ratchet does. Left
alone here: it is outside the engine, and the fixed-point class copies the
same way.

## Provenance and licensing

### Where the code came from

The vendored `fftsg.c` was textually identical between MuTap and AmbiTap
before DspTap consolidated the wrappers; the FFT's Tap lineage is the
AmbiTap/MuTap one. The C++20 port landed at Stage 2a (tap/DspTap#28) as
`include/tap/dsp/fft/split_radix.h`, beside the C and routed nowhere, with
the banner below in place; Stage 2b (tap/DspTap#31) routed `basic_real_fft`
at it; Stage 2c moved the C out of the shipping tree. Since 2c the library is
header-only (`tap::dsp` a pure INTERFACE target; `tap_dsp_fft` exists only
under `TAP_DSP_FFT_CMSIS`, carrying the CMSIS objects). From 2c until D6 a
reference copy of `fftsg.c` — with `fftsg_float.c`, because the float half
of the parity gate compared against `rdft_f` and would otherwise have been
port-vs-port — lived under `tests/reference/ooura/`, compiled by the gate
alone and declared to it by `tests/reference/ooura_rdft.h`. D6's condition
(both MuTap and MuTap-Max pin a tree containing 2c: MuTap `801204d` pins
DspTap `8350f13`, MuTap-Max `544e756` pins that MuTap) was met after Stage 4,
and D6 deleted the reference copy, its declaration header and the gate
(`test_fft_parity_ooura.cpp`, both its targets). The bit identity the gate
held was pinned from D6 until #42 by output fingerprints ("The bit-identity
record after D6" below), beside the routing proof (`test_fft_routing.cpp`)
and the independent oracle (`test_fft_oracle.cpp`). Provenance was visible,
while the port shipped, through: the port header's attribution banner,
`NOTICE.md`, `readme.txt` at `third_party/ooura/readme.txt` (D6), the
fingerprint pins with the record of their equality to the C, and one
glossary line in MuTap's `docs/itu-compliance.md` mapping "measured on
Ooura" to "the vendored C at DspTap ≤ `5ca3b1c` and the bit-identical port
from the Stage 2b SHA onward". Consumer *code* comments dropped "Ooura
packing" in favour of the numeric definition in the Stage 5 view. The engine
was named for what it was (`detail::split_radix_rdft`, D7), not for its
author; attribution was carried by the banner and the notices, not the
identifier. `detail::srdif_rdft` follows the same rule (split-radix
decimation in frequency).

**Since tap/DspTap#42, in the maintainer's judgement (`NOTICE.md`), DspTap
ships no code derived from the package.** #42 deleted the port and its
fingerprint pins; the floating profiles run the srdif engine, written
clean-room ("The floating engine, replaced clean-room" below), and the
fixed-point engine's post-pass was re-derived at #39 (next subsection). What
remains of the package in the repo is `readme.txt` at `third_party/ooura/`,
with `LICENSES/LicenseRef-Ooura.txt`, kept as the license record for the
trees consumers pinned before #42: from Stage 2a (#28) to #42's base they
carry the port, from #27 until #39 the fixed-point transcription, before #32
the C itself, and from #32 to before D6 (#36) the C as a test-only reference
copy under `tests/reference/ooura/` (`NOTICE.md`, "Which trees carry what").
The records of the port — this note's transliteration rules, bit-identity
record, licensing statement and banner below — are kept as history.

### The fixed-point post-pass, re-derived (tap/DspTap#39)

From Stage 3b (#27) until #39 the fixed-point engine's real post-pass, the
inverse's pre-pass, their coefficient table and the DC/Nyquist handling
were Ooura's `rftfsub` / `rftbsub` / `makect` formulas transcribed into
`fft_arith`'s operations, and `fixed_point.h` said so while its banner read
`SPDX: MIT` and `NOTICE.md` named only the port as a derivative (final
audit, 2026-09-26). The maintainer chose to re-derive the pass from the
literature rather than extend the `LicenseRef-Ooura` marking. The procedure,
as #39's history and review record it, including where it fell short:

1. **Hand-off.** The PR's first commit removes the transcribed code and
   every passage in the code, the tests and the fixed-point design record
   that stated its formulas, leaving 20 `<<CLEAN-ROOM>>` markers (the tree
   does not build there). It was written by the party that had read the
   package. **What it left, and should not have** (hostile review B of #39):
   the pre-existing post-pass table checksums (`TwiddleTableChecksumIsPinned`,
   N = 256 / 512 / 2048) and the 16 pre-existing output fingerprints
   (`OutputFingerprintIsPinned`, N = 512 / 2048) stayed pinned in
   `tests/test_fft_fixed.cpp` — hash values, not formulas, but together a
   bit-exact black-box check of the removed table and arithmetic — and
   structural remarks that are not formulas: the Welch model's post-pass
   stage `{1, 1, 1.0}` (a one-bit shift, unit power gain, a rotation on every
   output) in the test and the notebook, "the pre-pass takes pre + 1 bits",
   "+0.25 from the final one-bit shift on every component", and §5's
   "0.229 from the post-pass".
2. **Brief.** The implementer — a separate agent — received the half-length
   method as Cooley, Lewis & Welch (1970) and Sorensen et al. (1987) state
   it, rewritten in DspTap's convention (E/O decomposition, the (k, M − k)
   pair, DC/Nyquist, bin N/4, the inverse), the numeric contract (packing,
   exponents, the one-bit growth budget, saturation freedom, the floors to
   hold), and a list it could not access: `fft/split_radix.h`;
   `third_party/ooura/**` and any `fftsg*.c` / `fft4g*.c`, in the tree or in
   history; the hand-off diff and the history of the four affected files;
   this note's port and provenance sections, the audit doc and `NOTICE.md`;
   and every other FFT library's source, the vendored CMSIS-DSP included.
   The brief did not contain the (1 ± i W)/2 fold, the ½ − ½ sin
   coefficient or the transcription's statement order (review B); it did
   say that the pair (k, M − k) follows from one complex product W_N^k O[k].
3. **Access statement** (the implementer's reports). None of the listed
   material was viewed. Tree-wide greps for the markers and for pin values
   searched the listed files' contents and surfaced no content from them,
   only the audit doc's file name, which it filtered out; one grep of
   `README.md` printed one line of its Provenance section, about when Stage
   3b landed; heading lists were read to find the forbidden sections'
   boundaries. The same machine held further copies of the transcribed
   header outside the list (the main clone, the audit worktree, a
   reviewer's worktree, and about 87 files in the session's shared
   scratchpad); asked afterwards, the implementer reported opening none of
   them — the one scratchpad file it read was the conventions file the
   brief named — while its worktree shared the main clone's git directory.
   On the pins: it read `tests/test_fft_fixed.cpp` in full before writing
   code, so it saw the hex values; it wrote the arrangement and the table
   generator before building anything; the first run matched every checksum
   and all 16 fingerprints, and nothing was changed after any comparison;
   the rejected arrangements were weighed on paper and never coded. Before
   choosing, it noticed that §5's surviving "0.229" equals 2.75/12, the
   count of the arrangement it then chose. It also noted that, as a model,
   it may have met widely published real-FFT code in training.
4. **Result.** The implementer's arrangement: with u = Z[k] and
   v = conj Z[M − k] after the pass's one-bit shift, G = C_k (u − v) with
   C_k = (1 + i W_N^k)/2, then Z[k] ← u − G and Z[M − k] ← Z[M − k] + conj G
   (the inverse the same statements with conj C_k). Rounding cost per
   output component, simulated by review B in the trait's integer
   arithmetic (N = 512, 50 800 components, units of q²/12): 2.08 for this
   form; 2.81 with a +0.246 LSB mean bias for E and W·O formed separately
   with E halved; 4.10 for A·u + B·v. It is not the only form at that cost:
   v + D (u − v) with D = (1 − i W_N^k)/2 = 1 − C_k also takes two roundings,
   stays saturation-free, and gives the same integers as the chosen form
   except at exact half-LSB ties (0 differing components in the
   simulation), because round-half-up satisfies round(x − y) = x − round(y)
   away from ties. That is why bit identity with the transcription was to
   be expected of any single-product, two-rounding arrangement over a
   once-rounded table, and it is what the result showed: every output
   fingerprint and table checksum pinned before #39 reproduced unchanged,
   on the host and on all four newlib legs (and in CI on Windows and
   macOS).

**What the record supports, and what it does not.** The text — the code,
the comments, the table generator and the derivation in the docstrings —
was written without access to the package or to the removed code, from the
published method. It is not different code: side by side with `rftfsub` /
`rftbsub` (review B), the dataflow, the statement order, the signs, the
coefficient values and the table layout are the same; what differs is the
text (indexing by bin, pointer aliases, the kernel's `rotate<Inverse>`, DC
before the loop). The table generator's first-octant evaluation and its
½ − ½ sin split, which differ from the removed transcription, coincide with
`makect`'s structure; here they come from `make_twiddle_table`'s octant fold
and D9's no-contraction rule. The fold (u + v)/2 = u − (u − v)/2 into one
product is a derivation step beyond the cited equations, and the same step
Ooura's code takes. The bit identity was checked against the old pins the
hand-off left in the tree, not obtained blind.

The maintainer's judgement (`NOTICE.md`): the fixed-point engine is DspTap's
own and MIT. Not legal advice.

### The floating engine, replaced clean-room (tap/DspTap#42)

From Stage 2b (tap/DspTap#31, `bbfa48d`) until tap/DspTap#42 every floating
transform a consumer ran was the port of Ooura's `rdft`
(`include/tap/dsp/fft/split_radix.h`), a derivative work shipped under
`LicenseRef-Ooura AND MIT` whose redistribution relied on the package's
modification grant ("The licensing statement", below). On 2026-09-26 the
maintainer decided to remove the last code derived from the package from
what DspTap ships rather than rely on that grant for a derivative it does
not expressly license for distribution — #39 had done the same for the
fixed-point post-pass — and, with it, not to contact the author ("Draft
email to the author", below). #42 replaced the port with
`detail::srdif_rdft` (`fft/srdif.h`) at the same contract. The procedure, as
#42's description and its hand-off commits state it:

1. **Hand-off** (`17db855` and `05c81f1`, written by the party that had read
   the port). Removed: `include/tap/dsp/fft/split_radix.h` (2,606 lines);
   its bit-exact output pins, `tests/support/split_radix_fingerprints.h`
   and `tests/test_fft_split_radix_fingerprint.cpp` with its CMake target —
   unlike #39's hand-off, which left the per-N pins in the tree, no per-N
   oracle of the removed engine remained, but one whole-scenario checksum
   did: `bench/README.md`'s `rfft_f32_512` line (`0x662dd085b5b88325`, the
   x86-64 value), which the implementer read and rewrote (`e4b4ede`); the
   srdif engine's checksum differs from it; and every comment naming the
   port's internal routines or table recurrences: `test_fft_oracle.cpp`'s
   account of how the port formed its twiddles and its citation of
   `fftsg.c`'s statement of the inverse, `test_fft_routing.cpp`'s list of the
   sizes the port's dispatch branched on, and the routine names in `fft.h`'s
   and `tools/fingerprint/fingerprint.cpp`'s D9 paragraphs (`17db855`) and
   in this note's fp-contraction record (`05c81f1`). `<<CLEAN-ROOM>>`
   markers took their place; the tree did not build at `17db855`.
2. **Brief.** The implementer, a separate agent, was forbidden: every copy
   of the removed engine or of the package on the machine (the main clone,
   the MuTap and MuTap-Max trees and their submodules, MuTap's fingerprint
   files, other worktrees); git history before the hand-off, and the
   hand-off diffs themselves; this note's port and provenance sections ("Why
   split-radix …", "Transliteration rules for the port", "Provenance and
   licensing"), `NOTICE.md` and the audit doc; and every other FFT library's
   source, the vendored CMSIS-DSP included. It worked from DspTap's own
   `fft/fixed_point.h` and `fft/tables.h`, the contract text, the
   literature the brief named (the split-radix papers, for a half-length
   DIF kernel, bit reversal and post-pass: the family and pipeline were the
   brief's choice, written by the party that had read the port, and are
   also the port's; independence is claimed for the arrangement and the
   text, not the family) and what it cites itself (the list in
   `fft/srdif.h`), and a targets sheet of numbers (the port's accuracy
   against a quad-precision reference, its instruction counts, `.text`,
   heap and host timings; no code, but its heap formula implies the port's
   two-table layout, and its deterministic cells and the notebook's error
   tables are three-digit fingerprints of the port's output).
3. **Access statement** (the implementer's report, in #42's description).
   No forbidden item was opened. The near misses it reported: listing this
   note's section headings (to find the forbidden sections' boundaries);
   `git log --oneline` subject lines, which include the history before the
   hand-off; one tree-wide `git grep -c` that printed only a count; and a
   re-install of the shared pre-commit hook.
4. **What the hand-off did not remove**, stated so the record is not read
   as more than it is. In files the brief allowed, the implementer could
   read:
   - history prose naming the removed engine and its origin: `fft.h`'s file
     and class comments (the port as "the C++20 transliteration of Ooura's
     rdft", bit-identical to the vendored C; its size range as "the
     int-indexing bound of Ooura's arithmetic"; the vendored C's lazy table
     build), `README.md` (the engine table's row, the migration note's
     "instead of Ooura's `ip`/`w` workspace", Provenance and License),
     `CLAUDE.md`, the root `CMakeLists.txt`, `bench.yml`, `ci.yml`,
     `bench/README.md` (with the port's whole-scenario checksum, item 1),
     and `notebooks/fft.ipynb`'s executed output (prose, and per-N error
     tables of the port);
   - structural facts at that level: that the port was a split-radix `rdft`
     (the header's name, the type `detail::split_radix_rdft`); that it built
     its twiddle tables in the constructor from libm `cos` / `sin` (the
     contract table's host-identity row, this note's fp-contraction record,
     `CLAUDE.md`); that its statements were textually the C's and that GCC
     inlined its inner routines into one another (the fp-contraction
     record, after `05c81f1` removed the routines' names); the test name
     `ForwardMatchesOoura` in `test_fft_backend.cpp`; and the port's
     measured numbers in this note's size and instruction-count tables
     (numbers the targets sheet also gave); and the fp-contraction record's
     list of the sizes at which FMA moved the port's output (> 0 at N = 1024,
     4096, 16384, 65536; 0 at 2048, 8192, 32768), a weak hint of a
     size-parity dispatch;
   - the fixed-point engine's real post-pass
     (`fixed_point_rdft::real_post_pass`), which the brief directed the
     implementer to generalize, and which #39's record ("What the record
     supports, and what it does not", above) finds the same as the package's
     `rftfsub` / `rftbsub` in dataflow, statement order, signs, coefficient
     values and table layout. The srdif engine's real post-pass and inverse
     pre-pass are that arrangement in floating point with the scaling
     removed, so they are — by construction, not by independent derivation
     — the published half-length method in the arrangement the package's
     code also takes. The kernel, the permutation and the tables do not come
     by that path: the fixed-point kernel is radix-4, and the srdif engine's
     split-radix kernel, permutation (a precomputed swap list until the fix
     pass for review A, a reverse-carry counter since) and integer table
     generator were written from the literature its header cites.

   Outside the allowed files:
   - in files the brief forbade but the worktree still held (this note's
     port and provenance sections, the audit doc, `NOTICE.md`): the port's
     routine names, `cftf1st`'s run-time twiddle derivation, a verbatim
     `cftf161` statement and the size dispatch map; for these the control
     was the brief and the access statement, not removal;
   - review A's threads, which the fix pass (`5571445`) worked from: their
     author had built and measured the port from `main`; they quote no code
     of it, and what they suggest (a heap formula as a contract number, an
     O(√M) seed table from Karp's survey in place of the swap list, the
     kernel twiddles read from the post table's imaginary parts) was taken
     as the formula in `fft.h`, a reverse-carry counter with no table at
     all, and that reading.
5. **Structural comparison.** The hostile provenance review of #42 compared
   the srdif engine against Ooura's upstream `fftsg.c`, which the reviewer,
   unlike the implementer, was allowed to read (and the port, which follows
   it routine for routine). The review read `srdif.h` at `51914d0`; its
   srdif-side claims were re-checked here against the tree after the fix
   pass for review A (`5571445`), which changed the permutation and where
   the kernel twiddles come from (the `srdif.h` line numbers below are this
   commit's). The two share what the brief and the literature fix:
   split-radix in decimation in frequency (`fftsg.c`'s own header says
   "decimation: frequency, radix: split-radix"), a half-length complex
   kernel followed by the bit-reversal permutation and the
   Cooley–Lewis–Welch post-pass, the trivial-twiddle butterflies (j = 0, the
   eighth turn, the W₁₆ constants) and the symmetry W^(l/4 − j) = i·conj W^j
   that halves a twiddle table. They differ in every implementation choice.
   - *Passes.* Each of Ooura's passes (`cftf1st`, `cftmdl1`) is one
     split-radix level plus an untwiddled radix-2 step on the half, whose
     twiddles it defers into a second block type (`cftmdl2`, `cftf162`,
     `cftf082`), walked iteratively from the end of the array (`cftrec4` /
     `cfttree`), with separate inverse routines and a conjugating bit
     reversal. srdif's pass (`group`, `fused_pass`, `srdif.h:762-842`) is
     two full split-radix levels on groups of eight values, every block
     untwiddled, recursing first half first (`kernel`, `:914-946`: l/4,
     l/8, l/8, l/4, l/4), one butterfly for both directions (`bf`,
     `:701-754`).
   - *Leaves.* Ooura's are hand-written 16- and 8-point routines in two
     variants each (64 points split 16/16/16/16, `cftfx41`); srdif's are one
     template recursion held in registers (`leaf_blocks`, `:873-901`; 64
     points split 16/8/8/16/16, `block`, `:847-868`).
   - *Kernel twiddles.* Ooura's come from libm in per-level sub-tables of
     (cos θ, sin θ, cos 3θ, −sin 3θ) with secant factors (`makewt`), and
     `cftf1st` interpolates half of them at run time. srdif's hold twelve
     values per entry (W^j, W^3j, W^(j+M/8), W^(3j+3M/8), W^2j, W^6j) at the
     full resolution M, one table strided per level, with separate 32- and
     64-point tables for the compile-time blocks, nothing derived at run
     time; since the fix pass they are read at construction out of the post
     table's imaginary parts (`fill_kernel_table`, `:581-607`: cos θ/2
     doubled exactly, sin θ as the cosine of the complement, quarter turns
     by exchange and negation), whose values come from integer Taylor
     series.
   - *Post-pass table.* Ooura's (`makect`) is ½cos over the quarter turn,
     the first octant folded, with ½ − ½sin formed at use; srdif's
     (`fill_post_table`, `:546-567`) interleaves a once-rounded (1 − sin)/2
     with cos/2. The second column is the same function as Ooura's `c` (a
     value any table for this pass holds), formed from integer series rather
     than libm; reading the kernel twiddles out of it is srdif's own
     (`makewt` builds Ooura's kernel table apart). Both keep kernel and
     post-pass values in one array: Ooura's `w` (N/2 words, beside the `ip`
     seed table), srdif's single allocation (7N/8 + 36 Samples) since the
     fix pass.
   - *Permutation.* Ooura's (`bitrv2`, with fixed `bitrv216` / `bitrv208`
     for the 16- and 8-point cases) enumerates pairs of reversed middle
     indices from an O(√N) seed table (`makeipt`) in a nested j < k loop,
     with 16 (or 8) unrolled swaps per pair and a separate diagonal.
     srdif's (`permute`, `:621-660`) keeps no table: a reverse-carry
     counter (Gold and Rader) steps over the middle bits and, with the two
     end bit pairs split off, one step serves sixteen indices — six swaps
     unconditional, four when y < rev′(y); below M = 16 the counter walks
     every index. Splitting end bits off the middle so that one step serves
     several swaps is common to the two and is one of the arrangements
     Karp's survey describes; the enumeration differs (seed-table pairs
     against a counter), and srdif has no fixed-size variants.

   No identifier or comment is shared. The one statement-for-statement
   passage is the real post-pass pair (`post_pair`, `:982-1005`, against
   `rftfsub` / `rftbsub`: the same eight statements in the same order, one
   product's operands commuted; srdif's loop takes k = 1 alone and then two
   pairs per step, Ooura's one per step), which item 4 records as inherited
   through #39; beyond it, only idioms any implementation writes alike (the
   complex swap, `swap2`, and the eighth-turn product). The review's
   side-by-side table, with line numbers at `51914d0`, is on #42; its rows
   12 (bit reversal) and 16 (kernel twiddle table) describe the tree before
   the fix pass.
6. **Result.** Every floating output bit changed, at the error level of
   either engine: rms error against a quad-precision reference 4–14 % below
   the port's from N = 128 up, at most 4 % above it below that; two
   white-noise max-error cells and 54 of 516 single-realization cells above
   the targets' bars, which review A's paired A/B shows to be sampling
   draws: no (profile, N, signal, direction) cell is worse in the mean
   ("The floating engine (srdif)", "Accuracy against a quad-precision
   reference"). Their acceptance #42 puts to the maintainer. One bit-exact
   oracle of the port remained, a whole-scenario checksum (item 1); the
   srdif engine does not match it, and nothing was tuned toward the port's
   output.

**What the record supports, and what it does not.** The text of
`fft/srdif.h` — code, table generator, derivations — was written without
access to the port, the package or any other FFT library's source, from the
published literature and DspTap's own fixed-point engine. The real
post-pass is the stated exception: it inherits #39's arrangement, which #39
records as arithmetically the package's. The access statement is the
implementer's own report. The maintainer's judgement (`NOTICE.md`): since
#42 DspTap ships no code derived from the package. Not legal advice.

### The Hexagon tuning pass, clean-room (2026-09-28)

After #42 merged, MuTap's Hexagon ratchet (clang 19, v68 + HVX) measured the
srdif engine 12–13 % above the port on its chain workloads. No DspTap key
measured Hexagon, so "must not regress" had never been checked there. The
maintainer chose to tune srdif before the consumer took it. The procedure
followed #42's:

1. **What the orchestrator did, which had read the port.** It added the
   `hexagon` bench key (toolchain file, `icount.py` target, `bench.yml` job).
   It measured the port's Hexagon counts once, from a detached checkout of
   `7a58ebe` built with the same harness, then deleted that checkout. It also
   deleted every port build product it had made while diagnosing on Hexagon:
   a public-API microbenchmark's binary and its disassembly. What crossed to
   the implementer, all numbers, verbatim from the brief:
   - the port's scenario counts: 49,733,404 / 54,835,187 / 103,240,010 for
     `rfft_f32_512` / `rfft_f32_2048` / `rfft_f64_512`;
   - "Per transform pair (forward + inverse) at N = 512 on Hexagon float, a
     public-API microbenchmark read about 14.4 k instructions for the
     removed engine against 18.3 k for srdif; double about 31.4 k against
     35.0 k";
   - "at N = 64 srdif was already 1–2 % below", a second measurement of the
     port, at another size, that locates the gap.

   The brief also stated that the ratchet counts instructions, not packets.
   That was wrong, and the implementer showed it (item 4).
2. **The brief.** The implementer was barred from:
   - every copy of the port and the package, and git history before `17db855`;
   - the main DspTap clone, the MuTap and MuTap-Max trees (their submodules
     carry the port), the other worktrees, and the orchestrator's scratch
     except the brief, the conventions file, and #42's accuracy harness and
     targets sheet (`clean-engine/TARGETS.md`, `clean-engine/targets/**`,
     which use only the public API and were allowed at #42 too);
   - `NOTICE.md` and the audit doc;
   - every section of this note except the contract, fixed-point,
     fp-contraction, Stage 4, size/count and srdif sections;
   - every tap/DspTap and tap/MuTap pull request and review, and every other
     FFT library's source.

   The brief listed hypotheses to measure, all about the compiler, not about
   the port: clang ignoring srdif's GCC pragma, inlining boundaries, index
   arithmetic, the complex-multiply form, and the cost of the permutation.
3. **Access statement (the implementer's report).** It opened nothing on the
   list. The one exception: broad `grep`s over this file printed single lines
   from sections it could not read. These were:
   - every heading;
   - three history lines;
   - line 971 of the transliteration-rules history, which names the port's
     routines (`makect`, the `bitrv2*` family, `cftfsub` / `cftbsub`, the
     `cft*` leaves);
   - a line of the bit-identity record naming `bitrv2` / `bitrv2conj` and the
     port's 512-point leaf;
   - lines of the provenance and library-comparison sections describing
     srdif's own routines and a 16-point leaf's operation count.

   These were names and one-line prose, not code. It reported using none of
   them. Its history operations were a `git show --stat` of `28befbd` and a
   `git archive` of it for the baseline builds.
4. **Result.** "Hexagon (clang) tuning" in "The floating engine (srdif)"
   records the changes, what was tried and the counts:
   - forced inlining under clang;
   - post-pass pairs loaded before either is stored;
   - blocks of 16 and fewer run one level at a time in memory;
   - one fused-pass group per loop step for double;
   - one multiply-subtract form for float.

   Every Hexagon floating scenario ends 3–8 % below the port, and every
   Cortex-M float key 0.06–1.1 % below #42's baselines. No output bit moved
   at `-ffp-contract=off`.

   The implementer found that the `hexagon` key counts packets (up to four
   instructions each), not instructions, so it rewards packing as well as
   instruction count.

### Comparison with other FFT libraries (2026-09-27)

Review B compared the srdif engine with Ooura's package. A separate
provenance audit on #42 (at `e2e04f6`) asked whether DspTap's FFT engines
copy or closely follow any other well-known FFT implementation. It covered
`fft/srdif.h`, the Q15 / Q31 engine (`fft/fixed_point.h`, `fft/tables.h`,
`fft/fft_arith.h`), `fft.h`, `fft/spectrum.h` and the two backend
wrappers. Its full report, with every score and side-by-side, is on #42.
It is a statement of fact, not legal advice.

**Sources.** All come from official locations:

- FFTW 3.3.10, the release tarball (sha256 `56c93254…e26467`; its md5
  equals the published `.md5sum`). The scan covered the generated
  `dft/scalar/codelets`, `rdft/scalar/r2cf`, `r2cb` and `r2r` codelets
  (230 files) and the hand-written `kernel/`, `dft/`, `rdft/` and `api/`
  (179 files).
- KissFFT, GitHub `mborgerding/kissfft` at `e5e3fac`.
- pocketfft, GitHub `mreineck/pocketfft`: the C version on `master` at
  `81d171a` and the C++ header on `cpp` at `c90e55b`.
- PFFFT, the author's Bitbucket repository at `0aec032`, with the FFTPACK
  translation it bundles.
- CMSIS-DSP, the vendored sources plus upstream's cfft / rfft / radix /
  bit-reversal sources at the pinned `918014f`. The vendored files are
  byte-identical to upstream.

Numerical Recipes was not downloaded (its licence forbids it). `four1`,
`realft` and `twofft` were compared as the book describes them, from the
auditor's knowledge.

**Method.**

- **Signature greps in both directions:** each library's characteristic
  identifiers, macros, constants and comment phrases.
- **A MOSS-style token scan:** comments stripped; k-gram coverage and
  winnowed containment.
  - *normalized* mode, every identifier and literal replaced, at
    k = 15 / 20 / 25;
  - *raw* mode, identifiers kept, at k = 8 / 12.
- **Calibration of the scan** (coverage at normalized k = 20 / raw k = 12):
  - Ooura's `fftsg.c` against the port that transliterated it scores
    96 % / 95 %;
  - PFFFT against the FFTPACK it rewrote into SIMD macros scores
    9 % / 8 %;
  - independent pairs (KissFFT against pocketfft, pocketfft against
    `fftsg.c`) score 0–3 % / 0 %;
  - DspTap's non-FFT headers against the FFT files score up to 11 % / 2 %.
- **A structural side-by-side** of each engine against its closest
  comparators: decomposition, loops, twiddles, bit reversal, leaves and the
  real pass.

The scan did not flag the post-pass pair that #42 item 5 records as
statement-for-statement Ooura's. Its operands are spelled differently. A
low score therefore rules out copied text, not paraphrase, which is what
the structural reading is for.

**Results.** Against every comparator, the raw scan covers at most 1.0 %
of any DspTap FFT file at k = 12, and at most 0.7 % of `srdif.h` and
`fixed_point.h` at k = 8. That is below the same-author baseline. Three
normalized hits exceed the baseline, and each was read:

- `fixed_point.h` against PFFFT (22 %) and CMSIS (13 %): runs of
  `x = f(a, b);`. The trait calls `work::add(…)` normalize like
  `VADD(…)` and `vaddq(…)`, but the operands differ.
- `tables.h` against pocketfft C (23 %): the octant symmetry fill (below).
- `srdif.h` against FFTW codelets (5 %): generic loads and complex
  products.

No comparator-characteristic identifier occurs in the engines. None of the
following appears:

- FFTW's `KP…` constants, `FMA` / `FNMS` macros and `E` / `R` types;
- KissFFT's `C_MUL` and `super_twiddles`;
- FFTPACK's `ido`, `cc`, `ch` and `wa`;
- CMSIS's `CoefA` / `CoefB`;
- Numerical Recipes' `wpr`, `wpi`, `wtemp` and `h1r` … `h2i`.

**Verdict per comparator.**

- **FFTW: not derived.**
  - FFTW is a planner over genfft-generated straight-line codelets, with
    halfcomplex or r2c layouts and real Cooley–Tukey steps (`hc2hc`,
    `hc2c`).
  - It has no half-length complex transform with a post-pass, and no bit
    reversal.
  - Shared are only the multiplication by i through a temporary (srdif
    `fill_kernel_table`, `kernel/trig.c` `real_cexp`) and the split-radix
    operation count of a 16-point leaf (144 + 24, the algorithm's).
- **KissFFT: not derived.**
  - KissFFT is a mixed-radix decimation-in-time recursion.
  - Its real pass is `kiss_fftr`'s sum and difference with a halving,
    not one folded coefficient.
  - It has different twiddles and no packed DC / Nyquist.
- **pocketfft: not derived.**
  - pocketfft uses FFTPACK-style Stockham passes and a native real radix,
    with no post-pass.
  - Its twiddles come from a floating minimax polynomial for cos − 1 and a
    √n table.
  - The one textual parallel is forced by the mathematics. `tables.h`
    (`make_twiddle_table`, `:120-124`) and `calc_first_half` both write
    `table[2k] = table[2j + 1]; table[2k + 1] = table[2j];`: the identity
    W^(q−k) = i·conj W^k in an interleaved table. The loop control differs.
- **PFFFT: not derived.** Its scalar path is FFTPACK's real radix; its
  ordered output layout (DC at [0], Nyquist at [1]) is a convention, not
  code.
- **CMSIS-DSP: not derived.** It is the same pipeline family
  (half-length CFFT, bit reversal, split pass; for Q31, radix-4 DIF with
  per-stage scaling). The kernels differ:
  - f32: radix-8 with table bit reversal;
  - q31: truncating shifts, a twiddle-outer loop and a separate last stage;
  - the split pass: per bin, with A and B tables and four products per
    component.

  Two details are forced by the mathematics, not taken from CMSIS:
  - the radix-2² output placement, which a radix-4 DIF followed by plain
    bit reversal requires;
  - the ½(1 − sin θ, ±cos θ) coefficient, the ½-folded twiddle of either
    sign convention.

  `backends/cmsis.h` and `backends/accelerate.h` wrap the vendor APIs and
  copy nothing.
- **Numerical Recipes (knowledge only): not derived, as far as such a
  comparison reaches.** NR's `four1` is a radix-2 decimation in time with
  the trigonometric recurrence; `realft` uses the h1 / h2 sum-difference
  form. Shared are:
  - the reverse-carry idea behind srdif's small-M permutation, differently
    expressed (Gold and Rader, cited);
  - the half-length method;
  - the output contract (exp(+i), DC / Nyquist in the first pair, N/2
    round-trip gain), which DspTap kept from the pre-#42 contract.

**What the comparison supports.** No passage in DspTap's FFT engines is
statement-for-statement close to FFTW, KissFFT, pocketfft, PFFFT, CMSIS-DSP
or, as far as can be judged without its text, Numerical Recipes. The shared
elements are textbook-forced formulas and idioms and contract conventions.
The one statement-for-statement parallel in DspTap's FFT remains the real
post-pass pair against Ooura's `rftfsub` / `rftbsub` (item 5 above).

**Not covered:** older releases; libraries without public source (vDSP,
IPP); smaller libraries (FFTS, muFFT, meow_fft, djbfft); and the text of
Numerical Recipes.

### Consumer follow-ups

Ride the pin bumps; none is DspTap's to make. Line numbers below were read
from MuTap `origin/main` `0f6a7f0` and MuTap-Max `origin/main` `8850144` on
2026-09-23 (`git show origin/main:<path>`); re-verify against the tree each
bump starts from.

- **After tap/DspTap#42 (the srdif engine): MuTap re-measures and rewrites
  its notices; MuTap-Max follows through MuTap.** Read from MuTap
  `origin/main` `eb773e9` and MuTap-Max `origin/main` `4cfcca3` on
  2026-09-26.
  - Output bits. Every floating output moves (double on every build, float
    wherever the portable engine runs: MuTap's linux, Windows, macOS — MuTap
    forces vDSP off — M33, Hexagon and M55 split-radix legs; the M55 CMSIS
    leg's float path is unchanged, its double path is not). The fingerprint
    gate (tap/MuTap#64) re-records all nine legs from one CI run with its
    artifact procedure (`scripts/fingerprints.sh record`, the
    `fingerprints-<leg>` artifacts); MuTap's own instruction-count baselines
    on its portable-engine keys move with them. The float numbers
    `docs/itu-compliance.md` certifies "on Ooura" (its glossary: the
    vendored C at DspTap ≤ `5ca3b1c` and the bit-identical port from the
    Stage 2b SHA) are re-measured on the srdif engine, and the glossary gains
    the bump's DspTap SHA as the point from which they are srdif numbers.
  - ABI tag and backend string. `.github/workflows/ci.yml` lines 276, 327
    and 406 assert `backend=split_radix abi=fft_split_radix` in the
    `m55-split-radix`, `m33` and `hexagon` fingerprint logs; after the bump
    the harness prints `backend=srdif abi=fft_srdif` (the capi's backend
    string and `k_real_fft_abi_tag`), so the three lines change with it.
    `include/mutap/fft.h`'s range `static_assert`s hold unchanged (4 … 2^30
    for the portable engine); its comments naming `tap/dsp/fft/split_radix.h`
    and `fft_split_radix` (lines 11, 51, 59, 71, 131), and the same in
    `fdaf.h:40` and `fd_kalman.h:20`, are re-pointed at `fft/srdif.h` and
    `fft_srdif`.
  - Notices and prose. `THIRD_PARTY_NOTICES.md` lines 13–72 name the port as
    what compiles into every consumer; after the bump nothing derived from
    the package does, and the section becomes `NOTICE.md`'s history (the
    port until DspTap #42, the fixed-point transcription #27 … #38, the
    readme and `LICENSES/LicenseRef-Ooura.txt` kept as the license record
    for older pins). Likewise `README.md` lines 40–42 and 393–394,
    `CMakeLists.txt` lines 9–12, `HANDOFF.md:45`, `ci.yml` comments at lines
    34 and 240, and `docs/optimization.md:16` (as history). MuTap-Max:
    `THIRD_PARTY_NOTICES.md` lines 16–39 and 169, `README.md` lines 159–174,
    `CMakeLists.txt:77` and both externals' `CMakeLists.txt:23` describe the
    port as compiled in; rewritten at its re-pin past MuTap's bump. No code
    change there: MuTap-Max names no engine (grepped).
  - Memory. The srdif engine holds 1.73–1.82× the port's heap per object,
    in one allocation (N = 2048: 7,312 B float, 14,624 B double; the
    formula is in `fft.h`); a consumer with many FFT objects re-checks its
    budget. `sizeof(basic_real_fft)` is 40 B on LP64 (72 B before).
- **After tap/DspTap#40 (D4 expiry, D5): nothing to change.** Neither
  removed API is used by MuTap or MuTap-Max (grepped at MuTap `origin/main`,
  MuTap `edf160e` and MuTap-Max `origin/main`; MuTap builds against the tree
  at both SHAs with `-DMUTAP_WERROR=ON`). A future consumer that stages float
  data through the double FFT keeps its own `double` buffer and calls the
  same-type transforms.
- **After #39: MuTap (and MuTap-Max through it) pins past the re-derivation.**
  No code or number changes (no consumer uses the Q15/Q31 FFT, and the
  output is bit-identical), but MuTap's `THIRD_PARTY_NOTICES.md` says what
  remains of the package is the derived port, which is incomplete for any
  DspTap pin from #27 to before #39 (that tree's `fixed_point.h` is the
  transcription, under `SPDX: MIT`) and exact again once the pin is past
  #39. The next routine bump closes it; until then the notice is short by
  one file.
- **MuTap `.github/workflows/ci.yml`**, job `branchless-parity` (lines
  360–392; its `run:` block starts at 369) — compiles
  `submodules/dsptap/third_party/ooura/fftsg*.c` by path, so it breaks by
  construction at the 2c bump (the files moved) and is rewritten header-only.
  What changes, all inside the `run:` block:
  - lines 373–375, the comment "Compiles DspTap's C by path; Stage 2c (P18)
    rewrites this job when the C goes. Ooura is C: compile with gcc so g++
    does not name-mangle rdft/rdft_f." — goes (this is that rewrite);
  - lines 376–377, `gcc -O2 -c submodules/dsptap/third_party/ooura/fftsg.c -o
    /tmp/fftsg.o` and `gcc -O2 -c
    submodules/dsptap/third_party/ooura/fftsg_float.c -o /tmp/fftsg_float.o`
    — go;
  - line 380, the `g++` line: `/tmp/fftsg.o /tmp/fftsg_float.o` — go, so the
    compile line reads `g++ -std=c++20 -O2 -DMUTAP_SUPPRESSOR_BRANCHLESS=$bl
    -Iinclude -Isubmodules/dsptap/include tests/fingerprint_harness.cpp -o
    /tmp/parity_$bl` (the harness is `tests/fingerprint_harness.cpp`; it
    never references `rdft`/`rdft_f`, and hostile review B of tap/DspTap#32
    verified that it compiles, links and runs header-only against this tree
    for `bl=0/1` with the 14 `FINGERPRINT` lines identical).

  Nothing replaces the dropped lines: the job's subject is the suppressor's
  two forms, which never needed the C after 2b; DspTap's own parity gate was
  the C-vs-port check (the pinned fingerprints since D6). Not recommended:
  compiling the reference files from
  `submodules/dsptap/tests/reference/ooura/`, which would re-couple MuTap CI
  to a test fixture D6 deletes.
- **MuTap `THIRD_PARTY_NOTICES.md`** — rewritten with `NOTICE.md`'s
  statement, because the notice now lives in a header compiled into every
  external and no `submodules/dsptap/third_party/ooura/fftsg*.c` exists to
  cite: line 10 ("the Ooura FFT, which compiles into every consumer via
  `MuTap::fft`") becomes "via the header `tap::dsp` provides"; the section
  "Ooura FFT — `third_party/ooura/fftsg.c`" (lines 16–30) names what ships
  (the C++20 port `submodules/dsptap/include/tap/dsp/fft/split_radix.h`, a
  derivative work under `LicenseRef-Ooura AND MIT` whose redistribution
  relies on the modification grant), quotes the notice verbatim, points at
  `submodules/dsptap/third_party/ooura/readme.txt`, drops the "permissive
  grant" sentence (lines 28–30) and names the reference C under
  `submodules/dsptap/tests/reference/ooura/` as not shipped.
- **MuTap `README.md`** — lines 8–10, "Header-only C++20 apart from one tiny
  static target (`MuTap::fft`, the vendored Ooura `fftsg.c`)", is false twice
  over after 2c (MuTap's CMake has no `MuTap::fft` target and nothing
  compiled ships): "Header-only C++20." Line 375 of the tree listing,
  `third_party/ooura/   vendored Ooura FFT (see THIRD_PARTY_NOTICES.md)`,
  is re-pointed at what the submodule carries (the port under
  `submodules/dsptap/include/tap/dsp/fft/`, the readme at
  `submodules/dsptap/third_party/ooura/readme.txt`).
- **MuTap `docs/optimization.md`** — line 15 cites "DspTap's
  `third_party/ooura/fftsg_float.c`" for the autovectorization measurement;
  that was the vendored C at that pin, so it is cited as history ("the
  vendored Ooura C, since replaced by the bit-identical port
  `fft/split_radix.h`"), not by a path that no longer exists.
- **MuTap `docs/itu-compliance.md`** — the glossary line above, at the
  Stage 2b bump.
- **MuTap CMake** — no change: MuTap links `tap::dsp`, now INTERFACE-only,
  and never named `tap_dsp_fft`/`tap::dsp_fft` (grepped; MuTap-Max neither),
  so no alias is kept.
- **MuTap-Max** — audit Part 4 keeps MuTap-Max's per-external comments
  "until a bump proves they can go", and Part 8 item 7 asks for that proof;
  this PR is it: with `TAP_DSP_FFT_CMSIS` off (every MuTap-Max build)
  `tap::dsp` is a pure INTERFACE target with no C source under it (consumer
  link line `c++ -O3 -DNDEBUG main.cpp.o -o consumer`, no `.a` under the
  subdirectory). At the MuTap-Max re-pin, on `origin/main` `8850144`:
  - comments premised on a compiled Ooura target, to drop or rewrite as
    "header-only; `tap::dsp` is INTERFACE": `CMakeLists.txt` lines 74–75
    ("Header-only apart from the vendored Ooura FFT (MuTap::fft), which the
    `mutap` interface target carries") and 83–84 ("`mutap` (interface
    target: headers + cxx_std_20 + MuTap::fft)", "`MuTap::fft` (the vendored
    Ooura rdft the FDAF core calls)"), and
    `source/projects/mutap.aec_tilde/CMakeLists.txt` /
    `mutap.afc_tilde/CMakeLists.txt` lines 22–24 ("tap::dsp (the vendored
    Ooura rdft that tap::mu::basic_real_fft calls)");
  - two Windows workarounds whose stated reason is gone, to re-justify or
    remove in the bump's Windows build: `CMakeLists.txt` lines 21–25,
    `CMAKE_MSVC_RUNTIME_LIBRARY` forced to the static CRT "because the MuTap
    core's Ooura FFT target uses CMP0091-NEW, so it would otherwise default
    to the dynamic /MD and clash (LNK4098)" — no policy-NEW compiled target
    exists under `tap::dsp` any more; if min-api or another target still
    needs it the comment says which, otherwise the line goes; and lines
    29–34, `/FIalgorithm` guarded to `COMPILE_LANGUAGE:CXX` because
    "force-including a C++ STL header into the C sources (Ooura FFT) trips
    STL1003" — no C TU remains under `tap::dsp`; keep the guard only if
    another C source is in the build, and reword the comment either way;
  - notices: nothing to write unless MuTap-Max redistributes `NOTICE.md`
    separately (it picks the notice up transitively through MuTap).

- **At the first MuTap bump past D6** (read from MuTap `origin/main` `6a0d260`,
  2026-09-23): `THIRD_PARTY_NOTICES.md` lines 54–60, the "**Not shipped:**"
  paragraph, describe the reference copy under
  `submodules/dsptap/tests/reference/ooura/` as present and D6 as future;
  after the bump no copy exists. The paragraph becomes one sentence: DspTap
  carried a test-only reference copy of `fftsg.c` from its Stage 2c to its
  Decision D6 and carries none since; `readme.txt` is the record. MuTap's
  `docs/itu-compliance.md` glossary ("measured on Ooura") stays true and
  does not move. MuTap-Max (`544e756`) names no path under `tests/reference/`.

- **After the contract-safety PR (#41): nothing is required of MuTap or
  MuTap-Max; verified, not assumed** (MuTap `origin/main` `0385ea9`,
  MuTap-Max `4cfcca3`, read 2026-09-26). MuTap's own
  `cmake/arm-cortex-m55-mps3.cmake`, unedited, configures DspTap with the
  MVE-F check true and `TAP_DSP_FFT_CMSIS` ON, so its default M55 leg (and
  its `backend=cmsis abi=fft_cmsis` assertion), its split-radix M55 leg
  (explicit `-DTAP_DSP_FFT_CMSIS=OFF`, still honoured) and its `m55` icount
  key are unaffected; `ci.yml`'s comment "ON by default on the bare-metal
  M55 profile" stays true. `cmake/arm-cortex-m33-mps2.cmake:40`'s
  `set(TAP_DSP_FFT_CMSIS OFF CACHE …)` is now redundant — harmless; drop it
  (and reword its reason) at a routine bump if wanted. MuTap uses
  `log_mel` only through the capi / bridge (double profile, split-radix on
  every build: no reachable `fft_size` is newly rejected) and `pvoc` not at
  all in C++ (grepped); MuTap-Max builds for hosts only (CMSIS never ON),
  names neither class, and its `CMakeLists.txt:90` comment ("OFF for every
  host build") stays true. No mangled name, layout or output bit changes.

### The bit-identity record after D6

History since tap/DspTap#42, which deleted the port together with the pins
this subsection records (`17db855`). The srdif engine's own pins are
`tests/test_fft_srdif_fingerprint.cpp` / `tests/support/srdif_fingerprints.h`:
one row for every compiler and target CI runs, because its tables come from
integer arithmetic
("The floating engine (srdif)"). The record is kept because consumers pinned
the trees it describes (DspTap `0db95b6`, #36, through the base of #42,
`7a58ebe`), and the pins and the recipe at the end are how anyone re-checks
such a tree; the files it names are read there with `git show 7a58ebe:<path>`.

The parity gate (`tests/test_fft_parity_ooura.cpp`, with an informational
default-flags twin) ran from Stage 2a (tap/DspTap#28) through Stage 4
(tap/DspTap#35) on the three hosts and the four QEMU legs, memcmp identity
between the engine and the reference C at `-ffp-contract=off`. D6 retired it
with the C. What held D10 from D6 until #42 was
`tests/test_fft_split_radix_fingerprint.cpp` (own target, `-ffp-contract=off`,
MSVC at its default `/fp:precise`): FNV-1a-64 folds over the IEEE bit patterns
of `detail::split_radix_rdft`'s outputs, forward and inverse, for `double` and
`float`, at every power of two from 4 to 65536 (the QEMU legs stop at 4096,
`TAP_DSP_TEST_MAX_FFT_N`), over four libm-free materials; procedure and pins
in `tests/support/split_radix_fingerprints.h`. Every power of two, not a
sample of them: Ooura's dispatch branched on N (`cftf040` / `cftb040` at 8,
`cftf161` and `bitrv216` at 32, `cftfx41` at 64 and 128, the odd-log4 halves
of `bitrv2` / `bitrv2conj`, `cftleaf`'s 512-point leaf at odd powers ≥ 2048),
and the first version of this test, at five even powers, left about 480 lines
of the engine unexecuted (gcov) and passed a real rounding change in `cftf161`
(`y10r = wn4r * (x0r - x0i)` rewritten as `wn4r * x0r - wn4r * x0i`) that the
full set catches in ten cells (review of tap/DspTap#36).

**The invariant: the float row.** `split_radix.h` had no precision-specific
branch (no `if constexpr`, no `is_same`): every kernel statement was the same
template code for both precisions. The `float` instantiation rounded each
double libm result to float, and on every configuration measured that
absorbed the libm differences, so `float` was one value everywhere — glibc
under both CPU dispatches, newlib with and without a double-precision FPU,
the UCRT, Apple's libm, g++ and clang. A float pin that moved was therefore
always a change to the engine's arithmetic. `double` also agrees everywhere
up to N = 64 — not because those twiddles are exact (only N = 4 reads none;
N = 16 already reads cos π/4, ½ cos π/8 and ½ sin π/8) but because every libm
measured rounds them alike; it first differs at N = 128.

**Why the double rows were per C library build.** The port built its
tables from libm's `cos` / `sin` / `atan` exactly as `fftsg.c` does (D10), so
a last-bit libm difference moves the `double` outputs — and moved the C's
identically. A row identifies a libm *build including its run-time
dispatch*: x86-64 glibc selects FMA/AVX2 or SSE2 implementations of
`sin` / `cos` / `atan` by CPU, and the two differ from N = 8192 up, so glibc
is a pair of rows, not one value (a pre-Haswell machine, or a VM that masks
AVX2/FMA, takes the SSE2 row). A run passed when every cell equalled one row,
and printed which; the linux CI job ran the test a second time under
`GLIBC_TUNABLES=glibc.cpu.hwcaps=-AVX2,-FMA` and checks that each run matched
its row, so both are exercised on every run. The rows are keyed on the
macros that identify the libm (`__GLIBC__` and `__x86_64__`; `__NEWLIB__`,
`__arm__` and `__ARM_FP & 0x8`; `_MSC_VER` and `_M_X64`; `__APPLE__` and
`__aarch64__`), not on its version: glibc's row was measured on 2.39 and
holds until a glibc changes those functions, at which point the failure
message says so. A platform none of them names fails the double test after
printing its values, to be recorded, not skipped. A failure message
distinguishes the two cases: a double cell that matches no row while the
float pin at the same N holds points at the platform's libm version or
dispatch (check against the C, below, and add a row); a float pin that
moves means the kernel changed.

| N | float (one row, every configuration) | double N ≤ 64 (every configuration) |
|---|---|---|
| 4 | `ef9b60ebf0e8f587` / `7d9b9e4480f66924` | `f014008109d04382` / `a4e32575aed1c801` |
| 8 | `96742dab1b577d2c` / `02f73648a60330ca` | `85c71bbb0b4690b6` / `477a3c910df5bca9` |
| 16 | `a222e5c5b80bba3b` / `3a182210d01399a7` | `ea32e166ad289557` / `4a89b6e54a2b629e` |
| 32 | `06e774c9f71d5af5` / `68592652a0ecee43` | `b8aec28ef6d673eb` / `e8caa0863bb9b639` |
| 64 | `31eb423adbb65986` / `9d67d026c4f48ecc` | `5644ff7c7b0e6419` / `4d1bee09c0ebc556` |
| 128 | `c2f53e32cca6945f` / `d2da49eb3feccae5` | — (per library, below) |
| 256 | `1dc8ccaa8f8db08d` / `50306aa83541329d` | — (per library, below) |
| 512 | `3272fd6659a91c9a` / `77384618c779efa6` | — (per library, below) |
| 1024 | `fd0ffea8357b8ddd` / `526e4272bdaab623` | — (per library, below) |
| 2048 | `2e8a6f5fb458520d` / `23fc4df37b433aec` | — (per library, below) |
| 4096 | `6723317886e7c1f2` / `7a318e19506b1c24` | — (per library, below) |
| 8192 | `a5edafd7ec823192` / `530ea64869bf67ae` | — (per library, below) |
| 16384 | `e562a21922c4c589` / `ed246b24b4c2f3df` | — (per library, below) |
| 32768 | `5ce3d674258b6035` / `3ad5668f242936a3` | — (per library, below) |
| 65536 | `141744d674ac1c64` / `8c9d260bedee7129` | — (per library, below) |

| double, N | glibc, FMA/AVX2 dispatch | glibc, SSE2 dispatch | newlib, soft double (M4, M4F, M33) | newlib, DP FPU (M55) | MSVC x64, UCRT | arm64 macOS |
|---|---|---|---|---|---|---|
| 128 | `6a7199cab6054194` / `e27bedfa6fec4919` | `6a7199cab6054194` / `e27bedfa6fec4919` | `f94fa452b6b6631b` / `73c4adb6801f242a` | `f94fa452b6b6631b` / `73c4adb6801f242a` | `f94fa452b6b6631b` / `73c4adb6801f242a` | `f94fa452b6b6631b` / `73c4adb6801f242a` |
| 256 | `4e70f1158efcfe3d` / `b5d51404f905fea6` | `4e70f1158efcfe3d` / `b5d51404f905fea6` | `e59fe4d3e6326eb3` / `43318e1983699fd3` | `5c6c45c5898a558a` / `43318e1983699fd3` | `cabb356393c2dac7` / `5c0eeb0579868a88` | `1a682072ef54d745` / `a9e60433894b7013` |
| 512 | `21e35dd0330a809a` / `d5c2231d159d4698` | `21e35dd0330a809a` / `d5c2231d159d4698` | `5c2479007e97b9c6` / `00306b85c97aee64` | `8b67b1235ed466a9` / `b5065c9484181d81` | `5563d637de5e961d` / `10e1a0bd60bce082` | `d0a19e500f844c5c` / `76d38e82464248af` |
| 1024 | `f24f19a7fe3274a8` / `54dce063fb4161eb` | `f24f19a7fe3274a8` / `54dce063fb4161eb` | `878a01ec9aa4e1bd` / `068c7471dbde8796` | `f5992e25957ac958` / `6c2e7e635f21adf7` | `dd44bf873268ba83` / `abc8af10c100bee6` | `5cb137231c0a46e9` / `ff6e0741acfebeb3` |
| 2048 | `f3a7d6278d939926` / `1c95174b8f085c8e` | `f3a7d6278d939926` / `1c95174b8f085c8e` | `4f5075e73328a14a` / `731f6158e19aadd7` | `c1c49092ec9a4212` / `6ee7a9907f767089` | `e8a241c435bb2601` / `ca0d22bb81356f5a` | `379d6b9eef0fbe66` / `f80a8d3a5cf07d68` |
| 4096 | `753d7cfca461e82b` / `8c9e917a44067ff2` | `753d7cfca461e82b` / `8c9e917a44067ff2` | `58c61a3b29d4b55e` / `2d0482b12bc4e246` | `c42e9fb613c1c2c9` / `81d554f525c82f66` | `4e0db6e49731d286` / `bcbcc885236c3585` | `7909ca3d02d78252` / `4a9a5c439158323d` |
| 8192 | `6165e9ec3f1269a7` / `32d82388dc2a21fc` | `cc1adb16b8385bf9` / `655eced7f8b256f4` | not run | not run | `25213b8387ce17df` / `184a6995c63e5507` | `dab35c43e9136ed6` / `503f3b059e7ac263` |
| 16384 | `b54f372489309a93` / `3efc33f7c7e1e0ad` | `22e9355d131725ac` / `e8b892b022107e49` | not run | not run | `2e24dc611386f95a` / `b6f8aab7dc17962d` | `bb8b79e8470c05cf` / `0496245fe2a1a0e4` |
| 32768 | `59c3d563577f7644` / `7d40ace7f287957e` | `0b12f1160600a08b` / `33cb520f47ba6d5b` | not run | not run | `128e42969938f8d4` / `1fdf8d272bf2b6a1` | `588bd3585488aef2` / `032c156afe537e02` |
| 65536 | `b7909fad30a275af` / `e72fb50d719ef181` | `23f50b8c94aecb5b` / `aa919e1c04c85b21` | not run | not run | `b02d408a57e2ebde` / `288c58dad9580f4b` | `bae4e82f8d79d549` / `bfa3a82a87bbe16d` |

**How the pins were shown to be the C's.** (1) Locally, on the tree of
`6f6f77f` (the last `main` that carried the C) with a scratch TU that ran
this procedure on the reference C's `rdft` / `rdft_f` beside the engine,
built with the gate's own reference library at `-ffp-contract=off`, and then
again on the retiring PR's own commits with the C present: C equal to engine
in every cell on x86-64 g++ 13.3.0 and clang++ 18.1.3 (identical), under
both glibc dispatches, and on the M4 soft-float, M4F, M33 and M55 QEMU legs
(arm-none-eabi-gcc 13.2.1, QEMU 8.2.2). (2) In CI, the retiring PR kept the
C through its pinning commits and added
`fft_parity_ooura.ReferenceCHasThePinnedFingerprints` to the gate — the same
procedure on the C at all fifteen sizes, required equal to the engine's, to
the float row and to one of the platform's double rows — and ran it
verbosely on every leg. The Windows and macOS rows were taken from that
run (CI run 35920111160; both sides printed, equal in every cell), and the
commit that pinned them (`df187b0`) asserted the C against every pin on all
seven legs before the C was deleted. (3) Upstream: `fftsg.c` from a fresh
download of Ooura's `fft.tgz`, unmodified, compiled with a float
instantiation made the way `fftsg_float.c` made it, reproduces both glibc
rows (default and under the tunable above) and the float row at all fifteen
sizes.

**Re-verifying a pre-#42 tree against upstream, if it is ever needed** (a
consumer pinned to such a tree meets a new platform, or a doubt). Not for
the srdif engine, which has no upstream:

1. Download `fft.tgz` from Ooura's published page (currently
   `https://www.kurims.kyoto-u.ac.jp/~ooura/fft.tgz`; the 2001 address in the
   notice, `http://momonga.t.u-tokyo.ac.jp/~ooura/fft.html`, is historical).
   The archive checked at D6: `fft.tgz` SHA-256
   `52bb637c70b971958ec79c9c8752b1df5ff0218a4db4510e60826e0cb79b5296` (gzip
   header dated 2006-12-28), `fft/fftsg.c` SHA-256
   `21f8ea961f13284b0272798af99c5ad6081cd30f339c5d4bcf82c05b19fbdc2a`. The
   copy DspTap carried (`tests/reference/ooura/fftsg.c` at `6f6f77f`, git
   blob `5cc0a00196b483438b591525e0574f1582e0c07b`, SHA-256
   `063df78915a0fb66d788d8c3866b97b0fe9c4a54f9d2b1d86a137406ed761c26`)
   differs from it only by Tap's banner at the top and stripped trailing
   whitespace (checked with a whitespace-insensitive diff at D6).
2. Compile `fftsg.c` as C at `-ffp-contract=off`, and a float instantiation
   of it: `#define double float` plus the `_f` renames of every external
   function around `#include "fftsg.c"`, as `fftsg_float.c` did (`git show
   6f6f77f:tests/reference/ooura/fftsg_float.c`, blob
   `8640dea39184851d1cb23e0cf298227805c55d33`).
3. Run `tap::dsp::test::fingerprint_of` (the support header, `git show
   7a58ebe:tests/support/split_radix_fingerprints.h`) over a thin
   adapter that calls `rdft` / `rdft_f` with `ip[0] = 0` on the first call
   and the workspace geometry `readme.txt` prescribes (the adapter the
   parity gate had: `git show 6f6f77f:tests/test_fft_parity_ooura.cpp`),
   compiled at `-ffp-contract=off`, and compare with the pins for that
   platform (`matching_double_row` picks the row; on x86-64 glibc run it
   under both dispatches). Alternatively check out `6f6f77f` whole and run
   its gate beside the engine under test.

### The licensing statement

`NOTICE.md` is the canonical statement — it is what ships — and this section
gives the reasoning behind it. When the reasoning changes, `NOTICE.md`
changes in the same PR.

Ooura's terms are stated in `third_party/ooura/readme.txt`, the only upstream
license text (the banner at the top of the formerly vendored `fftsg.c` was
Tap's copy of it; the upstream file carries no notice of its own). They are:

> You may use, copy, modify this code for any purpose and without fee. You
> may distribute this ORIGINAL package.

What that grants: use, copying and modification of the code, for any purpose
and without fee; and distribution of the *original* package. Distribution of a
modified derivative is not expressly granted.

Until Stage 2c DspTap shipped one source file of the package, `fftsg.c`, plus
its `readme.txt`; from 2c until D6 it shipped `readme.txt` and carried
`fftsg.c` outside the shipping tree, under `tests/reference/ooura/`, for the
parity gate; since D6 it carries `readme.txt` alone (with
`LICENSES/LicenseRef-Ooura.txt`, the notice in the REUSE layout). The
reference `fftsg.c` was textually identical to the 2006-12-28 `fft.tgz` except
for a provenance banner Tap added at the top (the upstream file carries no
notice of its own; the notice lives in `readme.txt`) and stripped trailing
whitespace (checked at D6 against a fresh download; hashes in "The
bit-identity record after D6"). That was a partial copy of the original
package with the notice attached, not a modified transform; the maintainer's
reading was that it was within the intent of the distribution grant, and the
draft email below was to put the question to the author (never sent; below).
With the copy gone, only the derivative question below remained, until #42
removed that too.

What shipped from Stage 2b (tap/DspTap#31) until tap/DspTap#42, the C++20
port of `rdft`, was a **derivative work, not the ORIGINAL package**, and its
redistribution relied on the **modification grant** ("modify this code for
any purpose"). Precedent exists but is not relied on:
WebRTC/Chromium ship a modified `fft4g.c` under
`common_audio/third_party/ooura/`, and their `LICENSE` quotes a broader
notice — "You may use, copy, modify and distribute this code for any purpose
(include commercial use) and without fee. Please refer to this package when
you modify this code." — that is not in the `fft.tgz` readme. DspTap has not
traced the origin of that text and does not rely on it; whether it is the
author's, and applies to `fftsg.c`, was the most useful question the draft
email asked. None of those projects relicense the derived code, and neither
did DspTap: the port header carried Ooura's notice verbatim as the governing
terms for the derived portion, with SPDX `LicenseRef-Ooura AND MIT` (MIT only
for the wrapper and DspTap's additions), and a line stating that it is a
derivative work with the modifications copyright and date.
`third_party/ooura/readme.txt` stays at that path permanently as the license
record for the trees that carried derived code (every DspTap tree from Stage
2a, #28, to the base of #42, `7a58ebe`; and from #27 until #39 the
fixed-point transcription too); `fftsg.c` moved to the test reference tree at
Stage 2c (with `fftsg_float.c`, for the float half of the gate) and was
deleted at D6. "permissive" and "public" — the words the earlier notice
used — overstated the grant and are not used.

**Since tap/DspTap#42** the question the grant had to answer does not arise
for what ships: in the maintainer's judgement (`NOTICE.md`), DspTap ships no
code derived from the package (the port was deleted at #42, "The floating
engine, replaced clean-room"; the fixed-point transcription was re-derived
at #39). The analysis above stays because consumers pinned trees that carry
the port, and for those trees it is still the statement. The maintainer
decided on 2026-09-26 to replace the derived code rather than to ask the
author; the draft email below was never sent.

This is a maintainer judgement call, not legal advice.

### The port header's banner

Landed verbatim in `include/tap/dsp/fft/split_radix.h` at Stage 2a
(tap/DspTap#28), plus one line naming the modifications and their date, and
deleted with the header at tap/DspTap#42 (`17db855`); recorded here as what
the trees from #28 to `7a58ebe` carry. The four house lines came first
(STYLE.md §3: `@file`, `@brief`, SPDX, copyright), with both copyright
holders as lines 4–5, then a `//` prose block carrying the notice verbatim
and the derivative statement. The prose block was a **documented exception**
to the three-line house banner; the wave-2 port PR tagged it as such under
"Notes for the reviewer" so its tidy reviewer expected it. The banner names
no path that changes: `readme.txt`'s path is fixed by D6.

```cpp
/// @file split_radix.h
/// @brief Split-radix real FFT engine: a C++20 transliteration of Ooura's rdft.
// SPDX-License-Identifier: LicenseRef-Ooura AND MIT
// Copyright(C) 1996-2001 Takuya OOURA
// Copyright 2026 Timothy Place and the DspTap contributors (modifications).
// Derived from fftsg.c in Takuya Ooura's General Purpose FFT Package; this is
// a derivative work, not the ORIGINAL package.
//
// The original package's notice, which governs the derived portion (the
// transform itself), verbatim from its readme.txt:
//
//     Copyright(C) 1996-2001 Takuya OOURA
//         email: ooura@mmm.t.u-tokyo.ac.jp
//         download: http://momonga.t.u-tokyo.ac.jp/~ooura/fft.html
//         You may use, copy, modify this code for any purpose and
//         without fee. You may distribute this ORIGINAL package.
//
// The wrapper and DspTap's additions are MIT (see LICENSE). The full
// statement is NOTICE.md; the original package's readme.txt is kept in-tree
// at third_party/ooura/readme.txt, and LICENSES/LicenseRef-Ooura.txt carries
// the notice text the SPDX reference denotes.
```

`LicenseRef-Ooura` is a project-local SPDX license reference. It denotes the
`Copyright:` block of `third_party/ooura/readme.txt` (lines 140–145), which
`LICENSES/LicenseRef-Ooura.txt` reproduces verbatim in the REUSE layout so the
reference resolves to a file; that is also why the readme stays in-tree.
Since #42 no file DspTap ships carries the identifier; the two files stay as
the license record for the trees that do (`NOTICE.md`).

### Draft email to the author (not sent; superseded by the replacement, 2026-09-26)

**Resolution (2026-09-26): not sent.** The maintainer decided to remove the
derived code instead of asking the author for a statement on derivative
distribution: #39 re-derived the fixed-point post-pass and tap/DspTap#42
replaced the port ("The floating engine, replaced clean-room"), so no
derivative remains in what DspTap ships and the question has no subject.
Nothing was ever sent, by the maintainer or on the maintainer's behalf. The
draft stays below as the record of what would have been asked; `NOTICE.md`
records the resolution.

Per audit Part 5, the maintainer decides whether to ask the author for an
explicit statement on derivative distribution; the decision and any outcome
are recorded in `NOTICE.md`. **Status at Stage 2c (2026-09-23): not yet
attempted** — nothing has been sent, by the maintainer or on the maintainer's
behalf, and the Stage 2c PR flags the open decision to the maintainer rather
than taking it. This is the draft. Addresses: `ooura@kurims.kyoto-u.ac.jp` first
(published on the author's current homepage at kurims.kyoto-u.ac.jp, last
updated 2006-12-28, the site that serves the current `fft.tgz`, and used in
the `fft2d` readme), then `ooura@mmm.t.u-tokyo.ac.jp` (the address in the
`fft.tgz` notice) as secondary.

> **Subject:** Permission question: distributing a modified derivative of the
> General Purpose FFT Package (fftsg.c)
>
> Dear Dr. Ooura,
>
> Thank you for the General Purpose FFT Package. I maintain DspTap, a small
> open-source (MIT) library of DSP primitives for audio, at
> https://github.com/tap/DspTap. It carried your `fftsg.c` for some time,
> textually unchanged apart from a provenance comment we added at the top,
> with your `readme.txt` alongside.
>
> We would now like to carry a C++ translation of the `rdft` path of
> `fftsg.c` — the same algorithm and the same statements, rewritten as a C++
> template so one source serves float and double — and to distribute it as
> part of the library. Your notice grants use, copying and modification for
> any purpose, and distribution of the original package. Because a
> translation is a modified derivative rather than the original package, I
> want to ask directly rather than assume: are you willing to permit
> distribution of such a derivative, with your copyright notice and terms
> retained verbatim in the file and your readme kept alongside it?
>
> A related question: the WebRTC project distributes a modified `fft4g.c`
> under a notice that reads "You may use, copy, modify and distribute this
> code for any purpose (include commercial use) and without fee. Please refer
> to this package when you modify this code." That wording is broader than
> the one in the `fft.tgz` readme. Is it yours, and does it apply to
> `fftsg.c` as well?
>
> The derived code would remain under your terms (it would not be relicensed),
> and would be clearly marked as a modified derivative with a link back to
> your package. A one-line reply either way would be very welcome, and I will
> record it in the project's notices.
>
> With thanks and best regards,
>
> Timothy Place
> https://github.com/tap/DspTap

**Outcome:** not sent; superseded by the replacement of the derived code
(maintainer decision, 2026-09-26; `NOTICE.md`).

## Changelog of contract-affecting SHAs

Every SHA that moves a documented contract point, per profile, with the
consumer pins that were re-measured. The port (Stage 2b) was designed to add
*no* numeric row here for `double` or `float`, and it added none: its row
records the engine change and the deprecation, with every consumer pin
byte-identical. tap/DspTap#42 is the first row that moves floating output
bits.

| SHA | Stage | Profile(s) | What moved | Consumer pins re-measured |
|---|---|---|---|---|
| tap/DspTap#27 (`b08f6c6` on `main`) | 3b | Q15, Q31 | the profiles exist: `basic_real_fft<std::int16_t \| std::int32_t, Scaling>` returning an exponent, `scaling::fixed` / `scaling::block_floating`, the numbers in the table above | none (no consumer on fixed point) |
| tap/DspTap#31 (`bbfa48d` on `main`) | 2b | `double`, `float` | **no output bit**: `basic_real_fft` routes to `detail::split_radix_rdft` instead of the vendored C (bit-identical, both precisions); tables built in the constructor (no first-call cost); `forward(const float*, float*)` / `inverse(const float*, float*)` on the double engine `[[deprecated]]` (D5); fp-contraction policy stated (D9: no export) | MuTap fingerprint harness: 14 rows byte-identical, float rows included; `test_float32`, `test_g168`, `test_nn_suppressor` unchanged; DspTap icount baselines re-recorded to the port's counts (`bench/README.md`) |
| tap/DspTap#32 (`8350f13` on `main`) | 2c | none numerically | **no output bit and no contract point**: the vendored C leaves the shipping tree (`fftsg.c` and `fftsg_float.c` to `tests/reference/ooura/`, D6); `tap::dsp` is a pure INTERFACE target and `tap_dsp_fft` exists only under `TAP_DSP_FFT_CMSIS`; the `extern "C"` `rdft`/`cdft`/`rdft_f`/`cdft_f` declarations leave `fft.h` (a consumer that took them from there no longer links — none did: MuTap and MuTap-Max were grepped); the capi's `dsptap_fft_backend()` returns `"split_radix"` where it returned `"ooura"`; the Q15/Q31 ratchet scenarios are seeded and the `.text` ceilings set | MuTap fingerprint harness (scratch build of MuTap `0f6a7f0` with this tree as the submodule, `-DMUTAP_WERROR=ON`): 14 rows byte-identical to the current pin `b08f6c6`; DspTap icount at +0.00 % on every float key |
| tap/DspTap#35 | 4 | none numerically | **no output bit** (the Ooura gate and the class-vs-engine memcmp are green on every leg; fixed-point icount checksums identical): `basic_real_fft`'s second template argument is the engine for the floating profiles (default `default_real_fft_engine_t<Sample>`); `k_min_size` / `k_max_size` / `supports_size` / `k_is_shareable` added; construction requires `supports_size(size)` (`TAP_EXPECTS`, debug) — the one narrowing is the CMSIS engine's 32 … 4096, which was already the library's behaviour (undefined outside it); `basic_real_fft`, `basic_pvoc`, `basic_log_mel` and the aliases move into `inline namespace fft_split_radix | fft_cmsis | fft_vdsp` (mangled names change: every consumer image is rebuilt on its bump, and two images with different defaults no longer share symbols); backends move to `fft/backends/`; the Q31 engine's transforms are `const`. Fingerprint A/B of this tree vs `main` (`tools/fingerprint`, 12 lines, review 35b): identical at g++ default flags, at `-O3 -DNDEBUG`, at `-O3 -DNDEBUG -march=x86-64-v3`, and on the Cortex-M33 leg | none yet (MuTap's bump: fingerprints must be byte-identical; its own embedders adopt the tag in `tap::mu` — checklist above) |
| tap/DspTap#39 | post-program | Q15, Q31 | **no output bit**: the real post-pass / pre-pass, its table and the DC/Nyquist handling re-derived from the literature under a clean-room procedure ("The fixed-point post-pass, re-derived", including what the hand-off left in the tree), bit-identical to the transcription it replaces (every pinned fingerprint and checksum unchanged, as any single-product two-rounding arrangement predicts); fingerprints added at N = 64 and 1024 (the radix-2 stage); the inverse's antisymmetric-pair worst case added to the saturation sweep; the Q31 deviation and F(x) + F(−x) maxima pinned per size at N = 4096 … 65536 by two host-only tests (Q31 block floating 31–80 / 121–286 LSB measured, pinned at 2×); Welch-model ratio pins re-derived; Q15/Q31 icount +0.49 … +1.15 %, inside the band | none numerically (no consumer on fixed point); provenance: MuTap's `THIRD_PARTY_NOTICES.md` ("what remains of the package is the derived port") is incomplete for any pin in #27 … #38 and exact again once it pins past #39 (Consumer follow-ups) |
| tap/DspTap#41 | final audit A1–A4 | `float` (the CMSIS build); `log_mel`, `pvoc` | **no output bit**. Build: `TAP_DSP_FFT_CMSIS` defaults ON only where the compiler targets MVE-F (`__ARM_FEATURE_MVE & 2`) — for a consumer's non-Helium bare-metal Arm toolchain (Cortex-M7, M33, …, and Cortex-A / R) that previously defaulted ON, the float engine moves from CMSIS-DSP to split-radix: size range 32 … 4096 → 4 … 2^30, ABI tag `fft_cmsis` → `fft_split_radix`; unchanged for every in-tree leg and for MuTap's M55 (measured). Contract points: `basic_log_mel<Sample>::supports_geometry` (valid() AND the engine's size predicate) is the constructor's `@pre`; `basic_pvoc<Sample>`'s range is derived, `[max(64, engine min), min(2^28, engine max)]` (float 64 … 4096 under CMSIS, 64 … 2^20 under vDSP), with `k_min_size` / `k_max_size` / `supports_size`; the capi's `log_mel` / `pvoc` create and setters gate on them. The narrowing is only where release builds were already undefined behaviour. Documentation: no fault is promised out of range (measured wrong output, or heap corruption at N = 4); `int` is Q31 on the hosts and under clang, not under arm-none-eabi-gcc | none needed (MuTap `0385ea9` / MuTap-Max `4cfcca3` read: nothing changes for them; "Consumer follow-ups") |
| tap/DspTap#36 | D6 | none numerically | **no output bit and no contract point**: the reference C (`tests/reference/ooura/`), its declaration header and the parity gate with its informational twin are deleted; `tests/test_fft_split_radix_fingerprint.cpp` pins the engine's output bits at every power of two from 4 to 65536 (one float row; double per C library build, glibc as an FMA/SSE2 dispatch pair), measured equal to the C's on every leg ("The bit-identity record after D6"); nothing under `include/` changes but comments; the test-only cache variable `TAP_DSP_PARITY_MAX_N` is renamed `TAP_DSP_TEST_MAX_FFT_N` | none needed (no shipping code changed; icount ratchet expected +0.00 % on every key) |
| tap/DspTap#40 | D4 expiry, D5 | `double`, `float` | **API break, no output bit**: `forward(const float*, float*)` / `inverse(const float*, float*)` on the double engine are deleted (D5; they allocated per call and were not `noexcept`), and `basic_real_fft<float \| double, scaling::fixed>` is a `static_assert` naming the one-argument form (D4 expiry; `detail::floating_engine_of` deleted). Every remaining spelling instantiates the same types as before, so no instruction a consumer can still write changes other than assertion line numbers (`TAP_EXPECTS` `__LINE__` immediates at `-O0` without `NDEBUG`: fft.h's line numbers moved; objects byte-identical to `ddadb74` at `-O2`/`-O3`/`-Os` on g++, clang, M33, M55 and M4, review of #40); the class now fails with the D4 message alone (error-recovery `engine`), pinned by the `fft_compile_fail.*` ctests | **consumers verified unaffected**: no use of either API in MuTap `origin/main`, MuTap `edf160e` (MuTap-Max's pin) or MuTap-Max `origin/main` (`git grep`); MuTap at both SHAs builds against this tree with `-DMUTAP_WERROR=ON`; no pin to re-measure |
| tap/DspTap#42 | the floating engine replaced | `double`, `float` (the portable engine: every build for `double`; `float` wherever no backend is selected — linux, Windows, the M4 / M4F / M33 legs, the M55 with `TAP_DSP_FFT_CMSIS=OFF`, any build with the backends off); `log_mel`, `pvoc` and every consumer class through them | **output bits change** for `double` and `float` on the portable engine: `detail::split_radix_rdft` (the port of Ooura's `rdft`, `fft/split_radix.h`, deleted) → `detail::srdif_rdft` (`fft/srdif.h`, written clean-room), at the error level of either engine (rms vs a quad-precision reference 4–14 % lower from N = 128; the misses are in "The floating engine (srdif)"). vDSP (macOS) and CMSIS (M55) float paths unchanged. ABI tag `fft_split_radix` → `fft_srdif` (mangled names of `basic_real_fft`, `basic_pvoc`, `basic_log_mel` and every tagged embedder change; images built against the two trees do not coalesce); the capi's `dsptap_fft_backend()` returns `"srdif"` for `"split_radix"`. Unchanged: packing, sign, scale, size range 4 … 2^30, shareability, `noexcept` / allocation-free transforms, the API. Output bits are now one fingerprint row for every compiler and target CI runs at `-ffp-contract=off` (no libm); heap per object 1.73–1.82× in one allocation (`fft.h` states the formula; `sizeof(basic_real_fft)` 72 → 40 B on LP64); float icount −2.2 … −7.9 % on the four portable-engine keys (re-recorded) and MinSizeRel float `.text` −11 … −22 kB (ceilings re-recorded); host x86-64, same-machine A/B (informational): `-O3` float +0.5 … +5.9 %, double −3.4 … +3.6 %; `-march=native` float −9 … −18 %, double −7.5 … +3.5 %. Accuracy: no cell worse in the mean (review A's paired A/B). Provenance: in the maintainer's judgement DspTap ships no code derived from Ooura's package (`NOTICE.md`) | **every consumer re-measures** its floating pins once: MuTap's fingerprint gate (tap/MuTap#64) re-records all nine legs from one CI run with its artifact procedure, its ABI / backend assertions (`backend=split_radix abi=fft_split_radix` → `backend=srdif abi=fft_srdif`) and certified float numbers follow, and its notices drop the port ("Consumer follow-ups"); MuTap-Max via MuTap (notices only) |
