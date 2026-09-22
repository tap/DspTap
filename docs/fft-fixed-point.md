# The fixed-point real FFT: design note (Stage 3b)

Status: the Stage 3b design note for `include/tap/dsp/fft/fixed_point.h` and
`include/tap/dsp/fft/tables.h`. It is folded into `docs/fft-design.md` (the
sections marked `TODO(stage 3b)` there) at the wave-3 merge; until then this
file is the record. Plan of record: `docs/audit-fft-and-code-smells.md`,
Part 7 (design), Part 9 (battery), Part 13 (the `fft_arith` amendments),
Decisions D2 and D3. The arithmetic specification the kernel is written
against is `fft/fft_arith.h`; every number here is one the committed
battery (`tests/test_fft_fixed.cpp`) pins or prints: a `ctest -V` run, or
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

## 1. Kernel structure

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
`is_shareable` true for both fixed profiles and `docs/fft-design.md`'s
shareability row carries that as the design. What shipped: one transform at
a time per object. Q15 is not shareable across threads, because the in-place
int16 API needs the per-object int32 work buffer; Q31 touches no object
state during a transform (the caller's buffer and the three const tables
only), but its transforms are declared non-const in this stage to match the
floating engine's shape (the split-radix port's transforms are const). Stage
4 states shareability as an engine trait (Q31 true, Q15 false) and decides
transform constness across the two wave-2 engines; nothing is restructured
in 3b. Every transform returns the exponent and is `[[nodiscard]]`: under
block floating point a discarded exponent is a silent scale error.

## 2. Halving count and the exponent

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

## 3. Headroom: why nothing saturates

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

## 4. Twiddle quantization

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

## 5. The Welch noise model and the measured floors

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
  of two and the bias roughly doubles per stage instead. The battery sees it
  as the Q31 maxima all sitting at index 0 or 1: 4.25 LSB (fixed, N = 2048)
  and 15.99 LSB (block floating, N = 1024) in the saturation sweep, 6 / 62
  LSB in F(x) + F(−x), 41.5 reconstructed LSB in the BFP round trip (pinned
  83), with bins 1 – 3 and their mirrors carrying part of it through the
  slowly rotating early twiddles. Beyond N = 2048 no committed test measures
  it (the Stage 3c notebook does, to N = 65536, where it keeps growing with
  the stage count). Q15 never sees it (below the narrow's quantum).
  Convergent (half-even) rounding would remove it; that would be a different
  `shr_round` contract in `fft_arith.h`, pinned by its own numbers, and is
  not this kernel's decision to make.

## 6. The CMSIS-DSP convention (D3, decided)

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

## 7. What Stage 3c and the wave-3 merge add

- `rfft_q15_512`, `rfft_q31_512`, `rfft_q31_2048` ratchet scenarios in
  `bench/` (owned by the port agent this wave), seeded on `main`.
- This note folded into `docs/fft-design.md` §"Derivations" and the size /
  instruction-count tables from the QEMU legs. Corrections the fold-in must
  carry into that file's contract table: the Q31 forward-scale row says
  "exactly `X / N` in Q31"; the shipped contract is `X / 2N` (the one-bit
  input pre-shift, `fixed_scaling_exponent(N) = log2 N + 1`). The
  shareability row's "(design: true)" becomes Q31 true, Q15 false (§1). The
  `TODO(stage 3b)` cells take the pinned numbers from §2, §3 and §5.
- The analysis instruments typed over `std::span<const Sample>`, the capi and
  `dsptap_py` exposure of the four profiles, and `notebooks/fft.ipynb`
  executed against them, which re-measures the table in §5 on the shipping
  C++ and adds BFP exponent behaviour on speech-like material.
