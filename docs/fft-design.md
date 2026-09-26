# The real FFT: design note

*September 2026. This is the design note that is too long for a header — the
role `kaiser.h`'s design note plays for the FIR substrate. It is filled in
stage by stage; every entry that still depends on a measurement is marked
`TODO(stage …)` and names the stage whose PR supplies the number. Landed so
far: Stage 2a (the port beside the C, tap/DspTap#28), Stage 3b (the Q15 / Q31
profiles, tap/DspTap#27; its design record, formerly `docs/fft-fixed-point.md`,
is folded in below as "The fixed-point profiles"), Stage 2b (the routing
flip, tap/DspTap#31), Stage 2c (the C leaves the shipping tree, tap/DspTap#32)
and Stage 4 (the engine parameter and the ABI tag, tap/DspTap#35; "Stage 4"
below). Until every stage has landed, the plan of record is
[`audit-fft-and-code-smells.md`](audit-fft-and-code-smells.md) (parts are
cited as "audit Part N" below), and each PR links its stage.*

## Purpose and scope

`tap::dsp::basic_real_fft<Sample, Policy = detail::default_real_fft_policy_t<Sample>>`
(`include/tap/dsp/fft.h`) is a four-profile real-FFT contract — `double` (the
golden model) and `float` (the embedded profile) on an engine that is the
second template argument since Stage 4
(the C++20 split-radix engine by default; vDSP or CMSIS-Helium where the
build selects them for `float`), and the Q15 / Q31 fixed-point profiles from
`sample_traits.h` on an int32 radix-4 kernel. The contract is the packing, the
sign convention, the scaling, and the numbers in the table below; the engine
behind it changed (the vendored C until Stage 2b, the port of the same
split-radix algorithm since, bit-identical to it; the radix-4 int32 kernel for
fixed point) and the contract did not. Consumers hold
the transform by value and depend on exactly the surface audit Part 1
enumerates: `basic_real_fft(size_t)`, `forward_inplace`, `inverse_inplace`,
`forward`, `inverse` (out-of-place, `2/N` applied), `size()` / `num_bins()`,
copyability, the packing, `exp(+i)`, the unnormalized in-place inverse, and the
aliases `real_fft` / `real_fft32`.

Out of scope here: the complex transform (`cdft`), the DCT/DST family, and any
consumer-side spectrum arithmetic (that is the Stage 5 packed-spectrum view).

## Contract summary, per profile

The header owns these numbers and the tests pin them; this table is the
cross-reference, and an entry that names a test is pinned by it. Floating-point
rows are the shipping `fft.h` since the Stage 2b flip (tap/DspTap#31: the
split-radix engine, whose outputs are bit-identical to the vendored C's at
`5ca3b1c`, so no floating number moved); fixed-point rows are the Stage 3b
kernel as merged in tap/DspTap#27, with the measurements in "The fixed-point
profiles" below.

| Contract point | `double` (`real_fft`) | `float` (`real_fft32`) | Q15 (`real_fft_q15`, `real_fft_q15_bfp`) | Q31 (`real_fft_q31`, `real_fft_q31_bfp`) |
|---|---|---|---|---|
| Engine | `detail::split_radix_rdft<double>` (`fft/split_radix.h`) by default (`default_real_fft_engine_t<double>`); since Stage 4 the second template argument, so any `real_fft_engine` can be named | `detail::split_radix_rdft<float>` by default, or `detail::cmsis_real_fft_f32` (`fft/backends/cmsis.h`) / `detail::accelerate_real_fft_f32` (`fft/backends/accelerate.h`) when `TAP_DSP_FFT_CMSIS` / `TAP_DSP_FFT_ACCELERATE` is on; `test_fft_backend.cpp` is typed over the engines the leg can build and pins each to the split-radix engine at float epsilon in one binary; `test_fft_routing.cpp` pins the class to the named split-radix engine byte for byte on every leg and to the selected default | `detail::fixed_point_rdft<std::int16_t, Scaling>` (`fft/fixed_point.h`); the second argument is the Scaling policy on the fixed-point profiles | `detail::fixed_point_rdft<std::int32_t, Scaling>` |
| Size range (`k_min_size` … `k_max_size`, `supports_size(n)`; Stage 4) | 4 … 2^30 (the int-indexing bound; bit identity gated to 2^20, oracle to 65536) | split-radix 4 … 2^30; vDSP 4 … 2^20; **CMSIS 32 … 4096** (`arm_rfft_fast_init_f32`'s table; status checked since Stage 4) — `EngineRangesAreTheStatedNumbers` as static_asserts in `test_fft_engine.cpp`, `ConstructsAtTheRangeBounds` | 4 … 65536 | 4 … 65536 |
| Packing | `data[0]` = DC, `data[1]` = Nyquist, `data[2k]` / `data[2k+1]` = bin *k* re / im, `1 ≤ k < N/2`; `N/2 + 1` bins | same | same (`ImpulseHasFlatSpectrum`, `DcAndNyquistPacking`, the oracle's closed forms) | same |
| Sign | `W = exp(+2πi/N)`; imaginary parts conjugated vs the engineering DFT | same | same (`SignConventionIsPlusI`, `OnBinSineIsPlusHalfNInTheImaginarySlot`) | same |
| Forward scale | 1 (`A[k] = Σ a[j] W^(jk)`, unnormalized) | 1 | `fixed`: exactly `X / N` in Q0.15, `e = log2 N` (`FixedForwardScaleIsExactlyXOverN`); `block_floating`: `X · 2^-e`, `0 ≤ e ≤ log2 N` returned (`BfpExponentIsWithinRange`) | `fixed`: exactly `X / 2N` in Q0.31, `e = log2 N + 1` — the one-bit input pre-shift (`FixedForwardScaleIsExactlyXOverTwoN`); `block_floating`: `X · 2^-e`, `0 ≤ e ≤ log2 N + 1` |
| Inverse scale | `inverse_inplace` unnormalized (caller applies `2/N`); `inverse()` applies `2/N` | same | the unnormalized inverse over `2^e`, same constant `e` under `fixed` (`FixedInverseCarriesTheSameExponent`); `inverse()` applies no `2/N`; round trip `x == out · 2^(e_fwd + e_inv + 1 − log2 N)` (`RoundTripReconstructsInputPerPolicy`) | same |
| Round-trip identity | `x` reproduced to `1e-12` abs at N = 1024 (`RoundTripReproducesInput`) | `2e-5` abs at N = 1024 (`RoundTripReproducesInput`) | 0.51 reconstructed LSB at N = 1024 (fixed; the output narrowing's half LSB), pinned 1.02 (`RoundTripReproducesInput`) | 3.34 reconstructed LSB at N = 1024 (fixed), pinned 6.7; 41.5 under block floating point (the DC bias), pinned 83 |
| Saturation-free input | n/a (floating point) | n/a | full scale ±1: input placed with 2 guard bits (`<< 14`, not 16); the int32 kernel performs no saturating operation for any input; the output narrowing clamps only at the rail itself (0.5 LSB, the full-scale Nyquist alternation); worst case vs the golden model 0.50 / 0.75 LSB (fixed / BFP), pinned 1.0 / 1.5 (`SaturationFreeWorstCaseDoesNotWrap`) | `fixed`: one input pre-shift (−6 dB); `block_floating`: full scale, the headroom rule; worst case 4.25 / 15.99 LSB at index 0 (fixed / BFP), pinned 8.5 / 32 |
| Noise floor | 1.85e-16 relative 2-norm error of the forward vs the compensated-DFT oracle on full-scale uniform noise at N = 256, measured 2026-09-23 on x86-64 Linux (GCC 13.3.0 and Clang 18.1.3 -O3, glibc 2.39; both print the same value), pinned at 4× for the libm and fp-contraction spread across hosts (`fft_oracle_floor.DoubleForwardTracksCompensatedDft`); the spread the pin exists for is in the CI logs of tap/DspTap#31 (run 35847675386): 1.8928e-16 on `cortex-m4-softfp` / `cortex-m4f` / `cortex-m33` (jobs 107137662986 / 107137663129 / 107137662923) and 1.6643e-16 on `cortex-m55` (107137663056), newlib's libm and each leg's contraction; the hosted legs run ctest non-verbose and print no value. The Higham correctness envelope the oracle sweeps to N = 65536 is separate (`test_fft_oracle.cpp`) | 1.12e-7 relative 2-norm error vs `double` at N = 512 on the engine itself, so the same *engine* is measured on the M55 / macOS backend legs (where `basic_real_fft<float>` is CMSIS / vDSP), measured as above, pinned at 2× (`RealFftCrossPrecision.FloatEngineTracksDoubleAtN512`); the value moves in the last bits with libm and fp-contraction (D9), not the engine: 1.1236e-7 on x86-64 and the soft-float `cortex-m4-softfp` leg, 1.1650e-7 on the VFMA legs `cortex-m4f` / `cortex-m33` / `cortex-m55` (same run and jobs), all inside the pin — the audit's Part 6 N2 probe value was 1.105e-7 on other material; and `< 1e-6` at N = 1024 through `basic_real_fft` (`FloatTracksDouble`) | fixed: 0.29 LSB rms (the narrow's own rounding) at every N and level; per-bin SNR 69.1 / 66.5 / 60.2 dB at 0 dBFS, N = 256 / 512 / 2048; BFP 87 – 91 dB at 0 dBFS (`NoiseFloorTracksWelchModel`; the full table in "The fixed-point profiles" §5) | fixed: 0.67 – 0.75 LSB rms, level-independent; 152.0 / 149.4 / 142.8 dB at 0 dBFS; BFP 157 – 161 dB |
| Latency | 0 (block transform, no internal delay) | 0 | 0 | 0 |
| Alignment | none required on `Sample*` | none (vDSP's internal split buffers are placed by the wrapper, not the caller) | none | none |
| Shareability across threads (`k_is_shareable`, an engine trait re-exported by the class; Stage 4) | **true**: tables built in the constructor, transforms `const noexcept` | **true** on the split-radix engine, **false** on vDSP and CMSIS (per-object scratch); the class's transforms stay non-const on every profile (audit N11) and the trait is the statement; a shareable engine's transforms are `const` and the class `static_assert`s it (`ShareabilityIsTheHeadersNumber` in `test_fft_rt.cpp` and `test_fft_engine.cpp`) | **false**: the in-place int16 API needs the per-object int32 work buffer (a recorded deviation from audit Part 7, which listed both fixed profiles as shareable) | **true**: no mutable state during a transform; the engine's transforms are `const` behind `requires`-constrained overloads since Stage 4 |
| Real-time safety | transforms `noexcept`, allocation-free (`test_fft_rt.cpp`; the float-I/O-on-double overloads are the exception, `[[deprecated]]` since 2b and removed after one consumer cycle, D5) | same | same (`TransformsAreNoexcept`, `*AllocatesNothing`, `CopyProducesBitIdenticalOutput`) | same |
| NaN / denormals | NaN propagates to every bin; denormal input is slow on x86 without FTZ, and FTZ differs between the split-radix, vDSP and CMSIS builds | same | not applicable / none | not applicable / none |
| Host identity | not across hosts: libm's `cos` / `sin` last bit differs between glibc, newlib, UCRT and Apple (below) | same | one bit pattern on every host, pinned per profile, policy, direction and N = 64 / 512 / 1024 / 2048 (`OutputFingerprintIsPinned`, `TwiddleTableChecksumIsPinned`) | same |

## Why split-radix for floating point and a radix-4 int32 kernel for fixed point

**Floating point stays on Ooura's split-radix algorithm**, transliterated into
C++20 (Stage 2a). The reasons are contractual, not aesthetic: every consumer
float pin in MuTap and DspTap was measured on the vendored C, and the port
reproduces the C's table semantics for both precisions so that `float` stays
*bit-identical* (D10). The alternative — computing the float twiddles in double
"because it is more accurate" — was measured and is a wash (table max absolute
error 1.19e-7 either way; transform rms relative error at N = 512 is 1.105e-7
with the C's tables and 1.114e-7 with double-computed ones — the audit's Part 6
N2 reviewer-probe values; the committed re-measurement at Stage 2b is 1.12e-7
on the engine as it ships, pinned in the table above), so it was dropped. The `#define double float` build of the C is replaced by a
`template <std::floating_point Sample>` class whose helpers are private static
members, which removes the 76 global symbols and the lazy table build.

**Fixed point does not reuse Ooura's code, only its contract** (D2). Audit
Part 6 established that the split-radix graph survives only as a data-flow
graph with a different operation order, and that Ooura's first stage
(`cftf1st`) derives half its twiddles at run time as sums of two cosines up to
2.0, which no fixed-point twiddle format holds. The fixed-point profiles are
therefore their own radix-4 decimation-in-frequency complex kernel of length
N/2 (radix-2 final stage for odd `log2 M`) plus the real post-pass of the
half-length method (Cooley, Lewis and Welch 1970; Sorensen et al. 1987),
arranged for the int32 arithmetic (§1), so the packing, the `exp(+i)`
convention and the DC / Nyquist slots are identical. The design, summarized from audit Part 7:

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
- **The inverse has its own scaling.** Ooura's unnormalized inverse has
  structural gain N/2 and saturates on the first stage if run unscaled in fixed
  point. Under `fixed` the inverse halves per stage and the round-trip factor
  is a fixed power of two; under `block_floating` the round trip returns
  `x · 2^-(e_fwd + e_inv)`. Whether the output convention matches CMSIS-DSP's
  *documented* q15 / q31 scaling is a compatibility decision made in 3b.

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
- Ooura's own references for the split-radix algorithm are listed in
  `third_party/ooura/readme.txt` (Nussbaumer 1982; Burrus, *Notes on the FFT*).

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
floating engine's shape (the split-radix port's transforms are const). Stage
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
`G.inverse_inplace (unnormalized) == data · 2^e`. Since Ooura's unnormalized
round trip has gain N/2, a fixed-point round trip reconstructs

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
a rounding tie. Measured over the battery's sweep (sizes 4 … 2048, levels 0
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
- **Block floating point**: the rule keeps h ≥ growth + 1 before each stage,
  so input components are < 2^(30 − growth), complex magnitudes
  < 2^(30.5 − growth), outputs ≤ 2^30.5; the same half bit. For the post-pass
  and the pre-pass (growth 1) that is h ≥ 2 after the shift: components in
  [−2^29, 2^29), |a| < 2^29.5 = A, and the bullets above apply unchanged,
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
(`SaturationFreeWorstCaseDoesNotWrap`, sizes 4 … 2048). Measured, as the
test prints it: the int32 kernel performs no saturating operation for any of
them (a clamp or a wrap would show as an output sitting on a rail the golden
model does not reach, and none does); the largest deviation from the golden
model at the returned exponent is 0.500 LSB (Q15 fixed), 0.750 LSB (Q15
block floating), 4.250 LSB (Q31 fixed: the inverse of a full-scale constant,
N = 2048, index 0) and 15.988 LSB (Q31 block floating: the inverse of
full-scale binary noise, N = 1024, index 0), pinned at 1.0 / 1.5 / 8.5 / 32
(unchanged, to the printed digit, by the post-pass re-derivation of
2026-09-26 and by the antisymmetric-pairs pattern it added). The pins cover
the battery's N = 4 … 2048. The same sweep run one size at a time beyond it
(2026-09-26, x86-64 GCC 13.3.0, not committed as a test) keeps Q31 fixed
inside its pin — 4.70 / 4.25 / 4.25 / 4.74 / 4.72 LSB at N = 4096 / 8192 /
16384 / 32768 / 65536 — while the Q31 block-floating maximum grows with the
gap between the constant and the returned exponent: 31.0 / 32.0 / 59.0 /
76.0 / 80.0 LSB, always the inverse of full-scale binary noise at index 0
(the DC-path bias of §5), so "15.99 LSB, pinned 32" is a statement about
N ≤ 2048 only; F(x) + F(−x) likewise reaches 121 / 126 / 137 / 141 / 286
LSB there against the 124 pin. No output sat on a rail the golden model
does not reach at any of those sizes.
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

Two things can move a coefficient by one Q1.30 LSB between hosts: a libm
last-bit difference (glibc, newlib, UCRT, Apple), and fp-contraction of a
generator expression (Decision D9: a compiler that fuses `a - b*c` into one
rounding, which is the default on the macOS arm64 leg and on the Cortex-M55
leg, whose FP64 lets g++ contract). Either can carry a double that lies
within 2^-31 of a rounding boundary onto the other side. The generators
therefore contain no contractible expression. The post-pass generator
(`make_real_post_pass_table`) meets it by construction: its only doubles are
the angle (2π/n, a power-of-two scaling of the double π, times the integer
k) and half a libm `sin` or `cos` (a multiplication by the exact 0.5), and
nothing is ever added to a product in double; the real part
(1 − sin θ)/2 is formed on the integers as 2^29 − `make_coeff`(sin θ / 2),
never as the double 0.5 − 0.5·sin θ, which is exactly the `a - b*c` shape
an FMA-contracting compiler fuses. Like the kernel twiddles it evaluates
libm only on the first octant (k ≤ N/8, sine set equal to cosine at N/8) and
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
variance at the pass in fractions of full scale (about 10^-3 on full-scale
white noise under fixed scaling). Summed:

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
  slowly rotating early twiddles. Beyond N = 2048 no committed test measures
  it; the Stage 3c notebook (tap/DspTap#30) does, to N = 65536, and the
  reading depends on the unit. Counted in the LSB of the constant exponent,
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
- The **inverse is not CMSIS's**: it is Ooura's unnormalized inverse over
  2^e, i.e. (1/N) Σ X W^-jk divided by 2 (Q15) or 4 (Q31) under fixed
  scaling. The round-trip identity above is the contract.
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

A contract point, not a build detail (D9). Measured on the audit's baseline:
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
   fingerprints that replaced it (`tests/test_fft_split_radix_fingerprint.cpp`,
   "The bit-identity record after D6" below) are built at
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
   | local, at tap/DspTap#34 (Stage 6 reviews; #34 touches neither `fft.h` nor `fft/`) | g++ `-O3 -march=x86-64-v3`, default `-ffp-contract=fast` | pvoc's float output through `basic_real_fft<float>` | codegen moved (float `cftrec4` 194 → 168 `vfmadd`, double 224 → 254) | **float output bits moved with an edit to unrelated code in the same TU** (log_mel.h's constructor), every pvoc function instruction-identical; ~30 unrelated lines appended to the TU made the outputs identical again — the output depends on TU context, not only on flags |

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
   what else the TU inlines decides how `cftmdl2` is inlined into `cftrec4`
   / `cftleaf`, and with it which products fuse — so an edit to code that
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
   last-bit `cos` / `sin` differences already deny across hosts. What the
   header promises is therefore the observation, not a guarantee: the identity
   held at default flags on these compilers on these statements; g++ on
   x86-64 with `-march` is measured to move the float profile's last bits
   (clang does not) and is not claimed.
   A consumer that needs bit reproducibility across its own compilers sets the
   flag on its own targets.

A related, separate fact: libm `cos` / `sin` differ in the last bit between
glibc, newlib, UCRT and Apple, so *double* outputs are not identical across
hosts today either; a platform-independent table is a Stage 3b candidate (the
fixed-point twiddle generator lands there) and otherwise out of scope.

## Transliteration rules for the port

These live in the engine header (`include/tap/dsp/fft/split_radix.h`) so the
next person does not undo them; they are recorded here because they are what
the bit identity depends on (held by the parity gate from Stage 2a until D6,
by the pinned fingerprints since).

- **Statement fidelity.** Every Ooura statement stays textually intact. Clang
  contracts FMAs within a statement only; GCC contracts across statements after
  inlining. Refactoring `wk1r * x0r - wk1i * x0i` into a `cmul` helper, a
  lambda, or a policy call changes which products fuse and silently breaks bit
  identity on one compiler or the other.
- **Explicit `static_cast<double>` on every libm call.** In C, `cos(delta * j)`
  with a float `delta` calls the double `cos`; in C++ `std::cos` of a float
  calls `cosf`, and the tables then differ from the C. The parity test at an N
  where `makewt` takes its `nwh > 4` branch catches this, but only because the
  port keeps the C's table semantics.
- **Tables are built once in the constructor**, with the same precision
  semantics the C uses for each instantiation (float-typed locals for the
  float instantiation).
- **Index types may widen to `std::size_t`** (every loop bound is a `< m` /
  `> 0` form over non-negative values); `-Wconversion` stays on to catch a
  missed cast.
- **Load `m_w.data()` into a local once per transform** so aliasing analysis
  matches the C, which receives `w` by parameter.
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

## Stage 4: the engine parameter and the ABI tag (tap/DspTap#35)

Audit F4 and Decision D4. What landed, what was measured, and what it does
not close.

### The engine parameter

`basic_real_fft<Sample, Policy = detail::default_real_fft_policy_t<Sample>>`.
One parameter list, two meanings, so every spelling in use keeps compiling
with its meaning (pinned by `tests/test_fft_engine.cpp`):

| Spelling | Resolves to | Note |
|---|---|---|
| `basic_real_fft<double>`, `real_fft` | `<double, detail::split_radix_rdft<double>>` | always |
| `basic_real_fft<float>`, `real_fft32` | `<float, default_real_fft_engine_t<float>>` | the split-radix engine, or `detail::cmsis_real_fft_f32` under `TAP_DSP_FFT_CMSIS`, or `detail::accelerate_real_fft_f32` under `TAP_DSP_FFT_ACCELERATE`; selected in exactly one place in `fft.h` |
| `basic_real_fft<float, detail::split_radix_rdft<float>>` | that engine | new: an engine named explicitly, beside the accelerated default in the same binary |
| `basic_real_fft<float \| double, scaling::fixed>` | the selected engine | the pre-Stage-4 spelling every profile shared; accepted, resolves to the default engine, and is a **distinct type** from the one-argument form (same layout and code, different template arguments: +2,949 B x86-64 / +1,995 B M55 of duplicate wrappers when both are instantiated, 35a/F3). No in-tree writer since the #35 fix pass (the capi's seam uses `detail::default_real_fft_policy_t<Sample>`, so `fft_impl<float>` holds exactly `real_fft32`); **expires after one consumer cycle** (D4), then a `static_assert` |
| `basic_real_fft<std::int16_t \| std::int32_t, Scaling>` | the fixed-point specialization | unchanged; the second argument is the Scaling policy |

Any other second argument on a floating `Sample` must satisfy the concept
`tap::dsp::real_fft_engine<Engine, Sample>`: constructible from the size,
copyable, two `noexcept` in-place transforms, `size()`, and the three
contract constants `k_min_size`, `k_max_size`, `k_is_shareable`. The four
engines satisfy it; `int` and `scaling::fixed` do not (the class resolves
the legacy spelling before the concept is checked). The alternative D4 kept
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
| split-radix (`double`, `float`) | 4 … 2^30 | the int-indexing bound of Ooura's arithmetic (n is an `int`; every index and table offset stays below 2^31). Exercised: bit identity 4 … 65536 and 2^20 (the gate), the oracle 4 … 65536. Above 2^20 the transform is the same statements over larger tables and is not separately measured. 2^30 rather than 2^20 so the precondition does not narrow what a consumer could construct before Stage 4 ("a power of two ≥ 4") |
| vDSP (`float`, Apple) | 4 … 2^20 | vDSP documents no maximum; 2^20 is the bound the float oracle has always stated for this engine and the size the split-radix gate runs to. Through the class it is swept to 65536 on the macOS leg and pinned bin-for-bin at 512 / 2048; 65536 … 2^20 is inside the range on vDSP's word, not on a measurement here |
| CMSIS-DSP (`float`, Cortex-M55) | **32 … 4096** | `arm_rfft_fast_init_f32` dispatches through a `switch` over exactly {32, …, 4096} and returns `ARM_MATH_ARGUMENT_ERROR` otherwise, leaving the instance uninitialized. Until Stage 4 the wrapper ignored that status and fft.h promised "≥ 4" for every profile: N = 4 hard-faulted on the M55 leg (audit Part 13, found by the oracle). The constructor now checks the status (`TAP_EXPECTS`) as well as the range |
| fixed point (Q15, Q31) | 4 … 65536 | as at Stage 3b, now spelled through the same constants |

What construction does out of range: `TAP_EXPECTS` is the house
precondition (STYLE.md §4, `include/tap/dsp/detail/expects.h`): an
assertion in a debug build, nothing in a release build, never a throw
(some Tap consumers build `-fno-exceptions`). Every CI battery is Release /
MinSizeRel, so in CI the check is the constexpr predicate, pinned as
`static_assert`s on every leg (N = 16 and 8192 rejected under CMSIS,
accepted by the split-radix engine named in the same binary; 32 and 4096
constructed and round-tripped under CMSIS), plus a Debug-only death test
that runs in a local Debug configure (it passes there). **Release-mode
behaviour, decided (35a/F1):** a violated precondition in a release build
is undefined behaviour exactly as the plain `assert` it replaces was — for
the CMSIS engine the instance stays zero-initialized and the first transform
HardFaults; there is deliberately no defined fallback, because a branch in
the transforms would cost the hot path on every call for a case the
precondition excludes and would hand a consumer an object that silently
transforms nothing. `supports_size` is therefore the *mandatory* gate
wherever N comes from configuration. The capi's `dsptap_fft_create` applies
it per profile since the #35 fix pass (before that it hard-coded 4 and the
fixed-point 65536 and never read the float engines' range, so a vDSP or
CMSIS build of the capi would have passed 2^21 or 16 straight through);
the CMSIS wrapper also no longer narrows a size above `k_max_size` into the
library's `uint16_t` argument (65536 would arrive as 0). MuTap's config path
is on the bump checklist below. A repo-wide precondition policy remains its
own plan; `TAP_EXPECTS` is used for fft.h's power-of-two and size-range
preconditions (and the backends' own) and nothing else.
`tests/test_fft_oracle.cpp` reads each profile's sweep range from the class
instead of carrying the CMSIS range under `#if`.

### Shareability as an engine trait

`k_is_shareable` (the plan's `is_shareable`, under the house `k_` prefix):
split-radix true, Q31 true, vDSP / CMSIS / Q15 false. Decision on constness
across engines: **shareable implies const** — a shareable engine's
transforms are `const`, and the class `static_assert`s it, so that half of
the trait is what the compiler sees rather than a comment; the converse is
not asserted (an engine with const transforms over `mutable` scratch that
says `k_is_shareable = false` is accepted, and honest) and is a convention
the tests pin: for the six shipped instantiations const and shareable
coincide (`test_fft_rt.cpp`, `test_fft_engine.cpp`). The Q31 fixed-point
transforms became `const` behind `requires`-constrained overloads
(`k_widens` selects the Q15 pair, which goes through the per-object work
buffer and stays non-const). Also noted (35a/F7): `int` is `std::int32_t` on
every CI target, so `basic_real_fft<int>` compiles and is the Q31 profile. The class's own transforms stay non-const on every profile
(audit N11: the public API's constness must not depend on the selected
engine). The `std::span` overloads Part 9 listed were not added: Stage 4
added no overload to the transform surface, so there is nothing to equate
(`test_fft_rt.cpp` says so where Part 9 expected the assertion).

### The ABI tag: what it is, the evidence, what it does not change

`fft.h` opens `namespace tap::dsp::inline TAP_DSP_FFT_ABI` around
`basic_real_fft` and the aliases, with `TAP_DSP_FFT_ABI` one of
`fft_split_radix` / `fft_cmsis` / `fft_vdsp` from the selection (the plan
named them `fft_ooura` / `fft_vdsp` / `fft_cmsis`; `fft_ooura` predates 2c,
when the default became the split-radix engine, and the name is true today).
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

- At tap/DspTap#35: **0 of 65** weak names are shared (and 0 of all 69
  lines, the 4 `n` symbols included). The same three read

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
  reproduces).** One translation unit (`tests/abi/abi_tag_image.cpp`) is
  built into two loadable modules: image A with the configured build's
  defines (tag `fft_split_radix` on linux, `fft_vdsp` on macOS) and image B
  with `TAP_DSP_FFT_CMSIS` and a stub engine of a different layout on the
  include path (tag `fft_cmsis`); the host `dlopen`s A then B with
  `RTLD_GLOBAL` and asks each image what its embedders see through
  `extern "C"` probes, the member calls laundered through a volatile
  pointer so the loader, not the optimizer, decides the target. Measured on
  linux (g++ 13.3): an embedder of `basic_real_fft<float>` defined in plain
  `tap::dsp` and instantiated in image B **ran image A's code** — its
  `which()` returned `fft_split_radix` from inside the `fft_cmsis` image and
  the bound member saw a layout of 80 bytes against the real 184 — which is
  F4 live; the identical embedder defined inside `inline namespace
  TAP_DSP_FFT_ABI` reported its own tag and its own layout in both images,
  and `pvoc32` / `log_mel32` computed identical, finite checksums in both
  images with both loaded (`TaggedEmbedderIsImmuneToCrossImageCoalescing`;
  the untagged outcome is asserted on the linux/GCC host and printed as
  `[ measured ]` elsewhere, `UntaggedEmbedderIsWhatTheTagExistsFor`). Under
  `RTLD_LOCAL` on ELF nothing coalesces (the review measured it; the test
  runs one mode per process because the first `dlopen` decides an image's
  bindings). The test is hosted, non-Windows: the QEMU legs have no loader,
  and a Windows image exports nothing without `__declspec`, so the
  coalescing has no path there.

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
1. Wrap each of the five headers' class definitions in
   `namespace tap::mu::inline TAP_DSP_FFT_ABI { … }` (a different namespace
   from `tap::dsp::fft_split_radix`, which is fine: the tag only has to
   appear in the mangled name; the macro is in scope through
   `mutap/fft.h`). postfilter.h holds two of the six classes.
2. **No forward declaration of any tagged class anywhere** — a later
   `template <typename> class partitioned_fdaf;` in plain `tap::mu` declares
   a *different* class. None exists today in MuTap, MuTap-Max or DspTap
   (grepped by the review); it must stay that way. Explicit specialisations
   and `using` through the enclosing namespace remain legal.
3. **Validate N with `basic_real_fft<Sample>::supports_size(n)` wherever a
   block size comes from configuration**, with the CMSIS numbers (32 …
   4096) in the error text: on the M55 a config block size outside that
   range is a release-mode HardFault at the first `process()` otherwise.
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
- CMSIS-vs-split-radix parity runs under emulation on the Cortex-M55 QEMU
  leg, as it has since that leg landed (main's `test_fft_backend.cpp`
  already compared `basic_real_fft<float>` — CMSIS under the define — to
  the split-radix reference); what #35 adds is that the two engines are
  typed rows in one binary (`fft_backend_parity/cmsis` beside
  `/split_radix`) and that the named split-radix routing row runs there.
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

- `tests/test_fft_backend.cpp`: typed over `::testing::Types<split_radix_f
  [, accelerate_real_fft_f32 | cmsis_real_fft_f32]>`; parity at 512 / 2048,
  alignment stability and tonal accuracy at 512 / 2048 / 4096, per engine,
  against the split-radix engine called directly; a `static_assert` that
  the default float engine is one of the rows.
- `tests/test_fft_engine.cpp` (new): the engine-parameter contract as
  compile-time facts on every leg (above).
- `tests/test_fft_routing.cpp`: memcmp class-vs-engine for `double`, for
  the named split-radix float engine on every leg, and for the selected
  float default.
- `tests/test_fft_rt.cpp`: `ShareabilityIsTheHeadersNumber` per
  instantiation, with the compiler-checked half.
- `tests/test_fft_parity_ooura.cpp`: untouched, green on every leg — the
  proof that no float or double bit moved through the refactor (deleted at
  D6 with the reference C; "The bit-identity record after D6" below).
- `tests/test_fft_abi_tag.cpp` + `tests/abi/` (hosted, non-Windows): the
  two-image loader test of the tag (above).

## Size and instruction counts, per target

Seeded by the Stage 1b ratchet from the vendored C, then re-measured at 2b
(port), 3b (fixed point) and 4 (engine parameter). Numbers here carry the
SHA, toolchain and QEMU versions; the gate itself lives in
`bench/baselines.json` and the CI job, never in this file.

### `.text` per profile at N = 512

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
(D11: a re-record is for a move beyond the band). The +0.5 … +1.2 % is the
post-pass loop's code generation — its roundings sit where the previous
one's did, as the identical outputs show —
about 1,200 instructions per 512-point transform on `m33`, some 10 per bin
pair, most of it the int64 saturation sequences of the eight `add` / `sub`
per pair and the register pressure around them. A variant that loads the
four inputs into locals first measured worse on `m33` (+1.24 % /
+1.32 % / +0.94 %) and was not kept.

**Stage 2c (this PR): the fixed-point scenarios seeded.** Same run as the
sizes above (35861317022 at `21b3488`), compare mode on every key: the float
scenarios at +0.00 % against the 2b baselines on all five keys (the C left the
shipping tree — it remained under `tests/reference/ooura/` as the parity oracle
until D6 deleted it — and the shipping path did not move), the three
fixed-point scenarios recorded as new baselines (`bench/README.md`, one row
per key and scenario). Executed guest instructions for the whole scenario
binary, as for the float rows: 2^20 samples per direction through `forward()`
+ `inverse()` with the exponent folded after each block, construction and the
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
`TODO(after stage 4, needs a Mac)`: split-radix vs vDSP same-binary
*microbenchmark* on Apple Silicon. Stage 4 made the same-binary comparison
possible (`basic_real_fft<float, detail::split_radix_rdft<float>>` beside
the vDSP default in one binary, which is how `test_fft_backend.cpp` now runs
parity on the macOS leg), but the number itself is a wall-clock measurement
on a Mac that this repo's linux development environment cannot make. Of the
consumers' "~3× faster / ~3× fewer instructions vs autovectorized Ooura"
claims, the M55 one is replaced in `fft.h` and the README by the bench's own
whole-scenario figure (C / CMSIS 1.81× / 1.97× on the `m55` key, above); the
Apple one is stated there as what it is, MuTap's transform-only measurement
on the vendored C (tap/MuTap#31), and is not re-measured in this repo until
someone runs `tap_dsp_bench_fft` twice on a Mac (`-DTAP_DSP_FFT_ACCELERATE`
ON and OFF) and records the machine.

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
(`test_fft_parity_ooura.cpp`, both its targets). What remains of Ooura's
package in the repo is `readme.txt` at `third_party/ooura/`, permanently, as
the license record for the derived engine. The bit identity the gate held is
pinned since by output fingerprints ("The bit-identity record after D6"
below), beside the routing proof (`test_fft_routing.cpp`) and the
independent oracle (`test_fft_oracle.cpp`). Provenance stays visible
through: the port header's attribution banner, `NOTICE.md`, `readme.txt` at
`third_party/ooura/readme.txt` (D6), the fingerprint pins with the record of
their equality to the C, and one glossary line in MuTap's
`docs/itu-compliance.md` mapping "measured on Ooura" to "the vendored C at
DspTap ≤ `5ca3b1c` and the bit-identical port from the Stage 2b SHA
onward". Consumer *code* comments drop "Ooura packing"
in favour of the numeric definition in the Stage 5 view. The engine is named
for what it is (`detail::split_radix_rdft`, D7), not for its author;
attribution is carried by the banner and the notices, not the identifier.

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

### Consumer follow-ups

Ride the pin bumps; none is DspTap's to make. Line numbers below were read
from MuTap `origin/main` `0f6a7f0` and MuTap-Max `origin/main` `8850144` on
2026-09-23 (`git show origin/main:<path>`); re-verify against the tree each
bump starts from.

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

### The bit-identity record after D6

The parity gate (`tests/test_fft_parity_ooura.cpp`, with an informational
default-flags twin) ran from Stage 2a (tap/DspTap#28) through Stage 4
(tap/DspTap#35) on the three hosts and the four QEMU legs, memcmp identity
between the engine and the reference C at `-ffp-contract=off`. D6 retired it
with the C. What holds D10 since is `tests/test_fft_split_radix_fingerprint.cpp`
(own target, `-ffp-contract=off`, MSVC at its default `/fp:precise`): FNV-1a-64
folds over the IEEE bit patterns of `detail::split_radix_rdft`'s outputs,
forward and inverse, for `double` and `float`, at every power of two from 4 to
65536 (the QEMU legs stop at 4096, `TAP_DSP_TEST_MAX_FFT_N`), over four
libm-free materials; procedure and pins in
`tests/support/split_radix_fingerprints.h`. Every power of two, not a sample
of them: Ooura's dispatch branches on N (`cftf040` / `cftb040` at 8,
`cftf161` and `bitrv216` at 32, `cftfx41` at 64 and 128, the odd-log4 halves
of `bitrv2` / `bitrv2conj`, `cftleaf`'s 512-point leaf at odd powers ≥ 2048),
and the first version of this test, at five even powers, left about 480
lines of the engine unexecuted (gcov) and passed a real rounding change in
`cftf161` (`y10r = wn4r * (x0r - x0i)` rewritten as
`wn4r * x0r - wn4r * x0i`) that the full set catches in ten cells (review of
tap/DspTap#36).

**The invariant: the float row.** `split_radix.h` has no precision-specific
branch (no `if constexpr`, no `is_same`): every kernel statement is the same
template code for both precisions. The `float` instantiation rounds each
double libm result to float, and on every configuration measured that
absorbed the libm differences, so `float` is one value everywhere — glibc
under both CPU dispatches, newlib with and without a double-precision FPU,
the UCRT, Apple's libm, g++ and clang. A float pin that moves is therefore
always a change to the engine's arithmetic. `double` also agrees everywhere
up to N = 64 — not because those twiddles are exact (only N = 4 reads none;
N = 16 already reads cos π/4, ½ cos π/8 and ½ sin π/8) but because every libm
measured rounds them alike; it first differs at N = 128.

**Why the double rows are per C library build.** The engine builds its
tables from libm's `cos` / `sin` / `atan` exactly as `fftsg.c` does (D10), so
a last-bit libm difference moves the `double` outputs — and moved the C's
identically. A row identifies a libm *build including its run-time
dispatch*: x86-64 glibc selects FMA/AVX2 or SSE2 implementations of
`sin` / `cos` / `atan` by CPU, and the two differ from N = 8192 up, so glibc
is a pair of rows, not one value (a pre-Haswell machine, or a VM that masks
AVX2/FMA, takes the SSE2 row). A run passes when every cell equals one row,
and prints which; the linux CI job runs the test a second time under
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

**Re-verifying against upstream, if it is ever needed** (a pin moves, a new
platform, a doubt):

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
3. Run `tap::dsp::test::fingerprint_of` (the support header) over a thin
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
parity gate; since D6 it carries `readme.txt` alone. The reference `fftsg.c`
was textually identical to the 2006-12-28 `fft.tgz` except for a provenance
banner Tap added at the top (the upstream file carries no notice of its own;
the notice lives in `readme.txt`) and stripped trailing whitespace (checked
at D6 against a fresh download; hashes in "The bit-identity record after
D6"). That was a partial copy of the original package with the notice
attached, not a modified transform; the maintainer's reading was that it was
within the intent of the distribution grant, and the draft email below puts
the question to the author. With the copy gone, only the derivative question
below remains.

What ships, the C++20 port of `rdft`, is a **derivative work, not the
ORIGINAL package**, and its redistribution relies on the **modification grant**
("modify this code for any purpose"). Precedent exists but is not relied on:
WebRTC/Chromium ship a modified `fft4g.c` under
`common_audio/third_party/ooura/`, and their `LICENSE` quotes a broader
notice — "You may use, copy, modify and distribute this code for any purpose
(include commercial use) and without fee. Please refer to this package when
you modify this code." — that is not in the `fft.tgz` readme. DspTap has not
traced the origin of that text and does not rely on it; whether it is the
author's, and applies to `fftsg.c`, is the most useful question the draft
email asks. None of those projects relicense the derived code, and neither
does DspTap: the port header carries Ooura's notice verbatim as the governing
terms for the derived portion, with SPDX `LicenseRef-Ooura AND MIT` (MIT only
for the wrapper and DspTap's additions), and a line stating that it is a
derivative work with the modifications copyright and date.
`third_party/ooura/readme.txt` stays at that path permanently as the license
record for the derived code; `fftsg.c` moved to the test reference tree at
Stage 2c (with `fftsg_float.c`, for the float half of the gate) and was
deleted at D6. "permissive" and "public" — the words the earlier notice
used — overstated the grant and are not used.

This is a maintainer judgement call, not legal advice.

### The port header's banner

Landed verbatim in `include/tap/dsp/fft/split_radix.h` at Stage 2a
(tap/DspTap#28), plus one line naming the modifications and their date. The
four house lines come first (STYLE.md §3: `@file`, `@brief`, SPDX,
copyright), with both copyright holders as lines 4–5, then a `//` prose block
carrying the notice verbatim and the derivative statement. The prose block is
a **documented exception** to the three-line house banner; the wave-2 port PR
tags it as such under "Notes for the reviewer" so its tidy reviewer expects
it. The banner names no path that changes: `readme.txt`'s path is fixed by D6.

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

### Draft email to the author — for the maintainer to send

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

**Outcome:** recorded in `NOTICE.md` (sent on / not sent because; reply / no
reply by date), not here. As of Stage 2c: not yet attempted.

## Changelog of contract-affecting SHAs

Every SHA that moves a documented contract point, per profile, with the
consumer pins that were re-measured. The port (Stage 2b) was designed to add
*no* numeric row here for `double` or `float`, and it added none: its row
records the engine change and the deprecation, with every consumer pin
byte-identical.

| SHA | Stage | Profile(s) | What moved | Consumer pins re-measured |
|---|---|---|---|---|
| tap/DspTap#27 (`b08f6c6` on `main`) | 3b | Q15, Q31 | the profiles exist: `basic_real_fft<std::int16_t \| std::int32_t, Scaling>` returning an exponent, `scaling::fixed` / `scaling::block_floating`, the numbers in the table above | none (no consumer on fixed point) |
| tap/DspTap#31 (`bbfa48d` on `main`) | 2b | `double`, `float` | **no output bit**: `basic_real_fft` routes to `detail::split_radix_rdft` instead of the vendored C (bit-identical, both precisions); tables built in the constructor (no first-call cost); `forward(const float*, float*)` / `inverse(const float*, float*)` on the double engine `[[deprecated]]` (D5); fp-contraction policy stated (D9: no export) | MuTap fingerprint harness: 14 rows byte-identical, float rows included; `test_float32`, `test_g168`, `test_nn_suppressor` unchanged; DspTap icount baselines re-recorded to the port's counts (`bench/README.md`) |
| tap/DspTap#32 (`8350f13` on `main`) | 2c | none numerically | **no output bit and no contract point**: the vendored C leaves the shipping tree (`fftsg.c` and `fftsg_float.c` to `tests/reference/ooura/`, D6); `tap::dsp` is a pure INTERFACE target and `tap_dsp_fft` exists only under `TAP_DSP_FFT_CMSIS`; the `extern "C"` `rdft`/`cdft`/`rdft_f`/`cdft_f` declarations leave `fft.h` (a consumer that took them from there no longer links — none did: MuTap and MuTap-Max were grepped); the capi's `dsptap_fft_backend()` returns `"split_radix"` where it returned `"ooura"`; the Q15/Q31 ratchet scenarios are seeded and the `.text` ceilings set | MuTap fingerprint harness (scratch build of MuTap `0f6a7f0` with this tree as the submodule, `-DMUTAP_WERROR=ON`): 14 rows byte-identical to the current pin `b08f6c6`; DspTap icount at +0.00 % on every float key |
| tap/DspTap#35 | 4 | none numerically | **no output bit** (the Ooura gate and the class-vs-engine memcmp are green on every leg; fixed-point icount checksums identical): `basic_real_fft`'s second template argument is the engine for the floating profiles (default `default_real_fft_engine_t<Sample>`); `k_min_size` / `k_max_size` / `supports_size` / `k_is_shareable` added; construction requires `supports_size(size)` (`TAP_EXPECTS`, debug) — the one narrowing is the CMSIS engine's 32 … 4096, which was already the library's behaviour (undefined outside it); `basic_real_fft`, `basic_pvoc`, `basic_log_mel` and the aliases move into `inline namespace fft_split_radix | fft_cmsis | fft_vdsp` (mangled names change: every consumer image is rebuilt on its bump, and two images with different defaults no longer share symbols); backends move to `fft/backends/`; the Q31 engine's transforms are `const`. Fingerprint A/B of this tree vs `main` (`tools/fingerprint`, 12 lines, review 35b): identical at g++ default flags, at `-O3 -DNDEBUG`, at `-O3 -DNDEBUG -march=x86-64-v3`, and on the Cortex-M33 leg | none yet (MuTap's bump: fingerprints must be byte-identical; its own embedders adopt the tag in `tap::mu` — checklist above) |
| tap/DspTap#39 | post-program | Q15, Q31 | **no output bit**: the real post-pass / pre-pass, its table and the DC/Nyquist handling re-derived from the literature under a clean-room procedure ("The fixed-point post-pass, re-derived", including what the hand-off left in the tree), bit-identical to the transcription it replaces (every pinned fingerprint and checksum unchanged, as any single-product two-rounding arrangement predicts); fingerprints added at N = 64 and 1024 (the radix-2 stage); the inverse's antisymmetric-pair worst case added to the saturation sweep; the Q31 block-floating maximum scoped to N ≤ 2048 (31–80 LSB measured at 4096 … 65536, stated, not pinned); Welch-model ratio pins re-derived; Q15/Q31 icount +0.49 … +1.15 %, inside the band | none numerically (no consumer on fixed point); provenance: MuTap's `THIRD_PARTY_NOTICES.md` ("what remains of the package is the derived port") is incomplete for any pin in #27 … #38 and exact again once it pins past #39 (Consumer follow-ups) |
| tap/DspTap#36 | D6 | none numerically | **no output bit and no contract point**: the reference C (`tests/reference/ooura/`), its declaration header and the parity gate with its informational twin are deleted; `tests/test_fft_split_radix_fingerprint.cpp` pins the engine's output bits at every power of two from 4 to 65536 (one float row; double per C library build, glibc as an FMA/SSE2 dispatch pair), measured equal to the C's on every leg ("The bit-identity record after D6"); nothing under `include/` changes but comments; the test-only cache variable `TAP_DSP_PARITY_MAX_N` is renamed `TAP_DSP_TEST_MAX_FFT_N` | none needed (no shipping code changed; icount ratchet expected +0.00 % on every key) |
