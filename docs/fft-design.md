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

`tap::dsp::basic_real_fft<Sample, Engine>` (`include/tap/dsp/fft.h`) is a
four-profile real-FFT contract — `double` (the golden model) and `float` (the
embedded profile) on an engine that is a template parameter since Stage 4
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
| Host identity | not across hosts: libm's `cos` / `sin` last bit differs between glibc, newlib, UCRT and Apple (below) | same | one bit pattern on every host, pinned per profile, policy, direction and N = 512 / 2048 (`OutputFingerprintIsPinned`, `TwiddleTableChecksumIsPinned`) | same |

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
N/2 (radix-2 final stage for odd `log2 N`) plus Ooura's real post-pass
formulas, so the packing, the `exp(+i)` convention and the DC / Nyquist slots
are identical. The design, summarized from audit Part 7:

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
- The real post-pass formulas and the DC/Nyquist glue are transcribed from
  Ooura's `rdft` / `rftfsub` / `rftbsub` (`third_party/ooura/fftsg.c`); the
  complex kernel is not Ooura's (D2).

### 1. Kernel structure

The N real samples are read in place as M = N/2 complex values
z[j] = x[2j] + i x[2j+1]. The forward transform is

1. a radix-4 decimation-in-frequency complex FFT of length M in the library's
   convention W_M = exp(+2πi/M), with a radix-2 final stage when log2 M is
   odd (N = 4, 16, 64, ...);
2. the bit-reversal permutation of the M complex outputs
   (`make_bit_reversal_table`);
3. Ooura's real post-pass (`rftfsub`) over bins 1 … N/4 − 1 paired with
   their mirrors N/2 − k, bin N/4 untouched, and the DC/Nyquist glue
   a[0], a[1] = a[0] + a[1], a[0] − a[1].

The inverse mirrors `rdft`'s `isgn < 0` path: glue
a[0], a[1] = (a[0] + a[1])/2, (a[0] − a[1])/2, `rftbsub`, then the conjugate
kernel (same twiddle table, sine term subtracted), then the permutation.

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
which stays below 3M/4; the post-pass pairs (0.5 − 0.5 sin 2πk/N,
0.5 cos 2πk/N) for k in [0, N/4) as Ooura's `makect` defines them (|w| ≤ 1/√2);
the bit-reversal table over M. Memory per transform of size N: 4N bytes of
twiddles, 2N bytes of post-pass coefficients, 2N bytes of permutation, and
for Q15 a 4N-byte int32 work buffer.

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
| real post-pass / pre-pass | ≤ 2 (see §3) | 1 |
| input pre-shift, first stage only | — | `k_fixed_scaling_input_pre_shift` (Q15: 0, Q31: 1) |

The kernel's stages total log2 M = log2 N − 1 bits; the post-pass adds one;
so the forward exponent is

    e_fixed = log2 N + k_fixed_scaling_input_pre_shift  = fixed_scaling_exponent(N),

and the forward output is exactly X / N (Q15) or X / 2N (Q31) in the
sample's Q format. The inverse pre-pass takes the pre-shift and its own bit
before the glue and `rftbsub`; the kernel then halves per stage as in the
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
- **Post-pass**: after its 1-bit shift the values are ≤ B/2; with
  w = 0.5(1 − sin θ) + 0.5 i cos θ, |1 − w|² + |w|² = 1 and
  |1 − w| + |w| = √(0.5(1+sin θ)) + √(0.5(1−sin θ)) ≤ √2, so
  X = (1 − w) Z[k] + w conj Z[M−k] is ≤ √2 · B/2 ≤ B; the intermediate
  x = Z[k] − conj Z[M−k] is ≤ B and |y| = |w||x| ≤ B/√2. The glue's sums
  a[0] ± a[1] are ≤ 2 · B/2 = B per component.
- **Inverse**: the pre-pass takes pre-shift + 1 bits. Q31: input components
  ≤ 2^31 → ≤ 2^29 after two bits, |X| ≤ 2^29.5, |Z'| ≤ √2 · 2^29.5 = 2^30 =
  the kernel's B. Q15: one bit, |Z'| ≤ 2^29.
- **Block floating point**: the rule keeps h ≥ growth + 1 before each stage,
  so input components are < 2^(30 − growth), complex magnitudes
  < 2^(30.5 − growth), outputs ≤ 2^30.5; the same half bit.

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
full-scale binary noise, N = 1024, index 0), pinned at 1.0 / 1.5 / 8.5 / 32.
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
therefore contain no contractible expression: the post-pass coefficient is
written `0.5 * (1.0 - sin)`, the same double as `0.5 - 0.5 * sin` under
round-to-nearest (halving the singly rounded `1 - sin` is exact) but a
product of a difference, which no `-ffp-contract` mode fuses; verified by
inspecting the `-O2 -mfma` object code of both generators under g++ 13.3 and
clang 18.1 (no fused multiply-add remains, where the original expression
produced one under each). The battery pins each certified N's table checksum
(`table_checksum`, FNV-1a-64 over the bit patterns; `TwiddleTableChecksumIsPinned`)
so a remaining difference is detected, not absorbed, and on top of it the
transforms' own output fingerprints (`OutputFingerprintIsPinned`: FNV-1a-64
over the output block and the exponent, four profile × policy configurations,
both directions, N = 512 and 2048, sixteen values), which also fix the
rounding form of the complex multiply, the shift schedule, the sign
convention and the packing as one bit pattern per host. The octant-symmetric
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
radix-2 stage with 1/2, and through the post-pass with 1/4
((|1 − w|² + |w|²)/4 after its 1-bit shift); the post-pass adds its own
1-bit shift (3/48) and two rotations (2/12). Summed:

| N | stages | variance-only prediction, LSB32 rms (the battery's `model` column) |
|---|---|---|
| 256 | 4-4-4-2 + post | 0.569 (−191.53 dBFS/component) |
| 512 | 4-4-4-4 + post | 0.584 (−191.31) |
| 2048 | 4-4-4-4-4 + post | 0.585 (−191.30) |

(The battery's `NoiseFloorTracksWelchModel` derives this model independently
from the two headers and prints it beside every measurement; the numbers
here are its printout. Summing the rules above by hand gives the same
0.584 at N = 512: 0.229 from the post-pass, 0.078 from the last radix-4
stage through the post shift, 0.025 / 0.007 / 0.002 from the earlier ones.)
The last radix-4 stage and the post-pass dominate; the input pre-shift and
the first stages are attenuated by 4^-(stages) and contribute nothing
measurable, which is why the Q31 floor is the same at every N and the
prediction converges at 0.58 LSB.

What the variance model does not carry is the **round-half-up bias**: an
s-bit `shr_round` has mean error +2^-(s+1) LSB (+0.25 for one bit, +0.125
for two). The final 1-bit shift's +0.25 lands on every output component
unrotated; the last kernel stage's 4 · 0.125 = 0.5 per component is halved
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
| Q31 fixed | 256 | 9 | 0.69 | 1.46 | 152.0 dB | 111.2 dB |
| Q31 fixed | 512 | 10 | 0.68 | 1.35 | 149.4 dB | 109.4 dB |
| Q31 fixed | 2048 | 12 | 0.70 | 1.44 | 142.8 dB | 103.2 dB |
| Q15 block floating | 256 | 5 | 0.29 | 1.01 | 87.2 dB | 77.5 dB (0) |
| Q15 block floating | 512 | 5 | 0.29 | 1.01 | 90.6 dB | 80.6 dB (0) |
| Q15 block floating | 2048 | 7 | 0.29 | 1.03 | 84.2 dB | 86.4 dB (0) |
| Q31 block floating | 256 | 7 | 0.99 | 1.79 | 160.8 dB | 155.7 dB (0) |
| Q31 block floating | 512 | 7 | 1.50 | 2.32 | 160.6 dB | 155.7 dB (1) |
| Q31 block floating | 2048 | 9 | 1.12 | 1.28 | 156.8 dB | 153.8 dB (2) |

Across all levels the Q15 fixed rms is 0.26 – 0.30 LSB and the Q31 fixed rms
0.67 – 0.75 LSB (ratio 1.31 – 1.74); Q31 block floating holds 135 – 136 dB
at −60 dBFS with e = 0, its rms in LSB rising to 2.3 / 3.2 / 6.6 (N = 256 /
512 / 2048) as the output is held at a larger scale. The pinned ratios
(largest measured / model over the sizes, levels and both materials, at 2×):
noise 2.11 / 3.49 / 2.09 / 4.65, tone 0.044 / 0.77 / 1.07 / 1.29 for Q15
fixed / Q31 fixed / Q15 BFP / Q31 BFP.

Reading it against the model:

- **Q15 fixed** measures the narrow's own rounding (1/√12 = 0.289 LSB16)
  and nothing else: the int32 kernel's 0.58 LSB32 is 2^-14 of a Q15 LSB.
  The per-bin SNR of −40 dBFS white noise at N = 512 is 26.6 dB; Part 7
  quoted "about 25 dB" for exactly this case.
- **Q31 fixed** measures 0.67 – 0.75 LSB rms against the 0.57 – 0.59 LSB
  variance-only prediction, a power ratio of 1.31 – 1.74. The excess is the
  round-half-up bias: the final one-bit shift alone puts +0.25 LSB on every
  component, the last kernel stage's biases arrive half-shifted with
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

1. The Stage 2a parity target compiles both the C and the C++ with
   `-ffp-contract=off`, and *that* is the bit-identity gate, for `double` and
   `float`, forward and inverse, at every power of two from 4 to 65536 plus one
   run at 2^20, same binary and same libm on each host.
2. The default-flags run is informational and pinned at a *measured* bound —
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
   statements, not a guarantee, which is why the gate stays at
   `-ffp-contract=off`. Measured at 2b, and not claimed: x86-64 built with
   `-march` (FMA). g++ 13.3.0 `-O3 -march=haswell` keeps `double` identical
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
   observation about the builds measured, not a guarantee. Everything is
   still bit-identical at `-ffp-contract=off`, on clang, on MSVC and on the
   four QEMU legs. `fft.h`'s D9 paragraph states the same.
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
the parity gate depends on.

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
| `basic_real_fft<float \| double, scaling::fixed>` | the selected engine | the pre-Stage-4 spelling every profile shared (`tools/capi/dsptap_capi.cpp`'s generic seam writes it); accepted, resolves to the default engine, and is a **distinct type** from the one-argument form (same layout and code, different template arguments) |
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
that runs in a local Debug configure (it passes there). A violated
precondition in a release build is undefined behaviour exactly as the plain
`assert` it replaces was; the release-mode tool is `supports_size`, which
the capi's pattern (check, then construct) already follows. A repo-wide
precondition policy remains its own plan; `TAP_EXPECTS` is used for fft.h's
power-of-two and size-range preconditions (and the backends' own) and
nothing else. `tests/test_fft_oracle.cpp` reads each profile's sweep range
from the class instead of carrying the CMSIS range under `#if`.

### Shareability as an engine trait

`k_is_shareable` (the plan's `is_shareable`, under the house `k_` prefix):
split-radix true, Q31 true, vDSP / CMSIS / Q15 false. Decision on constness
across engines: an engine's transforms are `const` exactly when it is
shareable, and the class `static_assert`s it, so the trait is what the
compiler sees rather than a comment — the Q31 fixed-point transforms became
`const` behind `requires`-constrained overloads (`k_widens` selects the
Q15 pair, which goes through the per-object work buffer and stays
non-const). The class's own transforms stay non-const on every profile
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

- At `8350f13` (the branch base): **54 of 54** weak symbols have identical
  names in both objects, over two different layouts, among them

  ```
  W tap::dsp::basic_pvoc<float>::process(float, float)
  W tap::dsp::basic_real_fft<float, tap::dsp::scaling::fixed>::forward_inplace(float*)
  W tap::dsp::basic_real_fft<short, tap::dsp::scaling::fixed>::forward_inplace(short*)
  ```

- At tap/DspTap#35: **0 of 69** do. The same three read

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
  in the CMSIS one. `tests/test_fft_engine.cpp` pins the same property from
  inside the battery on every leg: `typeid(real_fft32).name()`,
  `typeid(pvoc32).name()` and `typeid(log_mel32).name()` contain the
  expected tag and none of the other two.

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

**What the tag does not close, recorded.**
- MuTap's `fdaf<float>`, `fd_kalman`, `pem_afc`, `postfilter`, `nn_suppressor`
  embed `basic_real_fft<Sample>` by value in `tap::mu`; their symbols do not
  carry the tag until MuTap opens the same inline namespace in its own
  namespace (`namespace tap::mu::inline TAP_DSP_FFT_ABI`), which is a
  one-line change per header on the bump that pins this tree. Until then
  the hazard F4 describes is closed for DspTap's own embedders and for
  `basic_real_fft` itself, and open one level up in MuTap exactly as before.
- CMSIS-vs-split-radix parity runs on the Cortex-M55 QEMU leg only (the plan
  said "compile-only on M55"; it in fact runs there, under emulation, with
  both engines in one binary — `fft_backend_parity/cmsis` beside
  `fft_backend_parity/split_radix`), nowhere on a host and nowhere on
  hardware. Recorded, not closed.
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
  proof that no float or double bit moved through the refactor.

## Size and instruction counts, per target

Seeded by the Stage 1b ratchet from the vendored C, then re-measured at 2b
(port), 3b (fixed point) and 4 (engine parameter). Numbers here carry the
SHA, toolchain and QEMU versions; the gate itself lives in
`bench/baselines.json` and the CI job, never in this file.

### `.text` per profile at N = 512

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

**Stage 2c (this PR): the fixed-point scenarios seeded.** Same run as the
sizes above (35861317022 at `21b3488`), compare mode on every key: the float
scenarios at +0.00 % against the 2b baselines on all five keys (the C left the
shipping tree — it remains under `tests/reference/ooura/` as the parity oracle
until D6's condition is met — and the shipping path did not move), the three
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
under `TAP_DSP_FFT_CMSIS`, carrying the CMSIS objects), and what remains of
Ooura's package in the repo is `readme.txt` at `third_party/ooura/` plus the
reference copy of `fftsg.c` — with `fftsg_float.c`, because the float half
of the parity gate compares against `rdft_f` and would otherwise be
port-vs-port — under `tests/reference/ooura/`, compiled by the gate alone
(`tests/CMakeLists.txt`) and declared to it by `tests/reference/ooura_rdft.h`.
See `NOTICE.md` for exactly how the reference file differs from upstream. D6
retires the reference copy after both MuTap and MuTap-Max pin a tree
containing 2c; with it goes the gate (`test_fft_parity_ooura.cpp`), leaving
the routing proof (`test_fft_routing.cpp`) and the independent oracle
(`test_fft_oracle.cpp`). Provenance stays visible after the port through four
things: the port header's attribution banner, `NOTICE.md`, the parity test
against `tests/reference/ooura/fftsg.c` (while it lasts;
`readme.txt` stays at `third_party/ooura/readme.txt` permanently, D6), and
one glossary line in MuTap's `docs/itu-compliance.md` mapping "measured on
Ooura" to "the vendored C at DspTap ≤ `5ca3b1c` and the bit-identical port
from the Stage 2b SHA onward". Consumer *code* comments drop "Ooura packing"
in favour of the numeric definition in the Stage 5 view. The engine is named
for what it is (`detail::split_radix_rdft`, D7), not for its author;
attribution is carried by the banner and the notices, not the identifier.

### Consumer follow-ups

Ride the pin bumps; none is DspTap's to make. Line numbers below were read
from MuTap `origin/main` `0f6a7f0` and MuTap-Max `origin/main` `8850144` on
2026-09-23 (`git show origin/main:<path>`); re-verify against the tree each
bump starts from.

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
  two forms, which never needed the C after 2b; DspTap's own parity gate is
  the C-vs-port check. Not recommended: compiling the reference files from
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

### The licensing statement

`NOTICE.md` is the canonical statement — it is what ships — and this section
gives the reasoning behind it. When the reasoning changes, `NOTICE.md`
changes in the same PR.

Ooura's terms are stated in `third_party/ooura/readme.txt`, the only upstream
license text (the banner at the top of the vendored `fftsg.c` is Tap's copy of
it; the upstream file carries no notice of its own). They are:

> You may use, copy, modify this code for any purpose and without fee. You
> may distribute this ORIGINAL package.

What that grants: use, copying and modification of the code, for any purpose
and without fee; and distribution of the *original* package. Distribution of a
modified derivative is not expressly granted.

Until Stage 2c DspTap shipped one source file of the package, `fftsg.c`, plus
its `readme.txt`; since 2c it ships `readme.txt` and carries `fftsg.c` outside
the shipping tree, under `tests/reference/ooura/`, for the parity gate. The
reference `fftsg.c` is textually identical to the 2006-12-28 `fft.tgz` except
for a provenance banner Tap added at the top (the upstream file carries no
notice of its own; the notice lives in `readme.txt`; the banner's pointer to
the readme was updated when the file moved) and stripped trailing
whitespace. That is a partial copy of the original package with the notice
attached, not a modified transform; the maintainer's reading is that it is
within the intent of the distribution grant, and the draft email below puts
the question to the author.

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
Stage 2c (with `fftsg_float.c`, for the float half of the gate) and is retired
afterwards (D6). "permissive" and "public" — the words the earlier notice
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
> https://github.com/tap/DspTap. It has carried your `fftsg.c` for some time,
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
| tap/DspTap#35 | 4 | none numerically | **no output bit** (the Ooura gate and the class-vs-engine memcmp are green on every leg; fixed-point icount checksums identical): `basic_real_fft`'s second template argument is the engine for the floating profiles (default `default_real_fft_engine_t<Sample>`); `k_min_size` / `k_max_size` / `supports_size` / `k_is_shareable` added; construction requires `supports_size(size)` (`TAP_EXPECTS`, debug) — the one narrowing is the CMSIS engine's 32 … 4096, which was already the library's behaviour (undefined outside it); `basic_real_fft`, `basic_pvoc`, `basic_log_mel` and the aliases move into `inline namespace fft_split_radix | fft_cmsis | fft_vdsp` (mangled names change: every consumer image is rebuilt on its bump, and two images with different defaults no longer share symbols); backends move to `fft/backends/`; the Q31 engine's transforms are `const` | none yet (MuTap's bump: fingerprints must be byte-identical; its own embedders should adopt the tag in `tap::mu`) |
