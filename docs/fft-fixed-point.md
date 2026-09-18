# The fixed-point real FFT: design note (Stage 3b)

Status: the Stage 3b design note for `include/tap/dsp/fft/fixed_point.h` and
`include/tap/dsp/fft/tables.h`. It is folded into `docs/fft-design.md` (the
sections marked `TODO(stage 3b)` there) at the wave-3 merge; until then this
file is the record. Plan of record: `docs/audit-fft-and-code-smells.md`,
Part 7 (design), Part 9 (battery), Part 13 (the `fft_arith` amendments),
Decisions D2 and D3. The arithmetic specification the kernel is written
against is `fft/fft_arith.h`; every number here was produced by the shipping
C++ against `basic_real_fft<double>` (the scratch programs are described in
the PR, not committed; the executed notebook and the analysis instruments
that will re-measure them are Stage 3c).

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
binds only when e already equals cum_fixed and h = 0, i.e. when the block is
at full scale and is, up to rounding noise, the fixed schedule's own data,
whose bound (§3) has half a bit to spare. An all-zero block returns e = 0
(the trait has no left shift; a quiet block is never normalized upward).

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
bins, in both directions and under both policies. Measured: no saturating
operation in any of them; the outputs match the golden model to 0.5 LSB
(Q15) and 2.3 / 5.6 LSB (Q31 forward / inverse), the same as on noise.

## 4. Twiddle quantization

Every coefficient is `make_coeff` of a double `cos`/`sin` (round half away
from zero), so |w_q − w| ≤ 0.5 LSB of Q1.30 = 2^-31 on every host, and the
angle 2πk/M formed in double is exact to < 2^-48 rad for M ≤ 2^15 (three
orders of magnitude below the quantum). 1.0 is representable (2^30), so DC
and Nyquist twiddles are exact and need no special case. Per rotation the
relative error is ≤ √2 · 2^-31; over the ≤ 8 rotating stages of N = 65536
it is ≤ 6 · 10^-9 relative, which on a full-scale Q31 output (2^31) is under
13 LSB and on any Q15 output is invisible. On an on-bin full-scale tone the
measured Q31 error is 0.07–0.19 LSB rms (N = 2048 … 256), i.e. the twiddle
term is below the rounding floor at every N in the battery.

Host libm last-bit differences can move a double that lies within 2^-31 of
a rounding boundary to the other side; the battery pins each certified N's
table checksum (`table_checksum`, FNV-1a-64 over the bit patterns) so such
a difference is detected, not absorbed. The kernel's outputs are otherwise
integer arithmetic and host-identical.

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

| N | stages | variance-only prediction, LSB32 rms |
|---|---|---|
| 256 | 4-4-4-2 + post | 0.527 |
| 512 | 4-4-4-4 + post | 0.545 |
| 2048 | 4-4-4-4-4 + post | 0.545 |

The last radix-4 stage and the post-pass dominate; the input pre-shift and
the first stages are attenuated by 4^-(stages) and contribute nothing
measurable, which is why the Q31 floor is the same at every N and the
prediction converges at 0.55 LSB.

What the variance model does not carry is the **round-half-up bias**: an
s-bit `shr_round` has mean error +2^-(s+1) LSB (+0.25 for one bit, +0.125
for two). The final 1-bit shift's +0.25 lands on every output component
unrotated; the last kernel stage's 4 · 0.125 = 0.5 per component is halved
by the post shift and mixed with bin-dependent phases; earlier stages' biases
arrive rotated and count as variance across bins. Along the unrotated path
(q = 0, r = 0 at every stage) the biases add coherently into the DC bin.

**Measured** (white noise, forward, output-referred; the level does not
change the floor under fixed scaling, so one row per N covers 0 … −60 dBFS;
the full per-level table with the on-bin tone and the inverse is in the PR
body and will be re-measured by the Stage 3c notebook):

| profile, policy | N | rms error, LSB | of which mean bias | max, LSB | per-bin SNR at 0 dBFS | at −40 dBFS |
|---|---|---|---|---|---|---|
| Q15 fixed | 256 | 0.29 | — | 0.50 | 68.7 dB | 28.6 dB |
| Q15 fixed | 512 | 0.30 | — | 0.50 | 65.8 dB | 25.9 dB |
| Q15 fixed | 2048 | 0.29 | — | 0.50 | 60.1 dB | 20.3 dB |
| Q31 fixed | 256 | 0.72 | +0.26 | 2.5 | 151.0 dB | 110.6 dB |
| Q31 fixed | 512 | 0.69 | +0.22 | 1.9 | 148.9 dB | 108.9 dB |
| Q31 fixed | 2048 | 0.71 | +0.22 | 2.8 | 142.7 dB | 102.8 dB |
| Q15 block floating | 256 | 0.29 (e = 5 … 0) | — | 0.50 | 86.9 dB | 77.3 dB |
| Q15 block floating | 512 | 0.28 (e = 5 … 0) | — | 0.50 | 90.4 dB | 80.5 dB |
| Q15 block floating | 2048 | 0.29 (e = 7 … 0) | — | 0.50 | 84.1 dB | 86.5 dB |
| Q31 block floating | 256 | 0.95 … 2.2 (e = 7 … 0) | — | 4.9 … 9.7 | 160.7 dB | 155.6 dB |
| Q31 block floating | 512 | 1.6 … 3.3 (e = 7 … 0) | — | 14 … 11 | 159.6 dB | 155.3 dB |
| Q31 block floating | 2048 | 1.1 … 6.6 (e = 9 … 0) | — | 12 … 28 | 156.8 dB | 153.8 dB |

Reading it against the model:

- **Q15 fixed** measures the narrow's own rounding (1/√12 = 0.289 LSB16)
  and nothing else: the int32 kernel's 0.55 LSB32 is 2^-14 of a Q15 LSB.
  The per-bin SNR of −40 dBFS white noise at N = 512 is 25.9 dB; Part 7
  quoted "about 25 dB" for exactly this case.
- **Q31 fixed** measures 0.69–0.73 LSB rms against the 0.53–0.55 LSB
  variance-only prediction. Removing the mean bias (+0.22–0.28 LSB, the
  final shift's +0.25) leaves 0.65–0.68 LSB of spread about the mean,
  about 20% above the variance model, from the earlier stages' biases
  arriving with rotated phases. The inverse measures 0.76–0.83 LSB (its
  pre-pass rotations act on unshifted data, so their 2/12 propagates through
  the whole kernel instead of being the last term). Both directions are
  level-independent to within the seed's spread (0.65–0.76 at every level
  from 0 to −60 dBFS), which is the Welch statement.
- **Block floating point** moves the floor with the block: Q15 gains 18–24 dB
  at full scale and 48–66 dB at −40 dBFS over fixed; Q31 gains 8–14 dB at
  full scale and holds 135–136 dB at −60 dBFS with e = 0, where the input's
  own 2^-31 quantization is the floor. The Q31 BFP rms in LSB grows as e
  falls because the output is held at a larger scale; the SNR is the
  comparable number.
- **The DC bias (honest limit).** Under fixed scaling the DC bin's mean
  error is +1.0 … +2.1 LSB (Q31, N = 256 … 65536): each stage's 4 · 0.125
  sums coherently but each later 2-bit shift quarters it, so it grows by
  about one LSB every few stages. Under block floating point a stage
  typically shifts one bit against a magnitude gain of two, so the DC bias
  doubles per stage: ~15 LSB at N = 2048 and ~140 LSB at N = 65536 on
  full-scale noise (−143 dBFS; 3% of the block's noise power), with bins
  1–3 and their mirrors carrying about half of that through the
  slowly-rotating early twiddles. Q15 never sees it (below the narrow's
  quantum). Convergent (half-even) rounding would remove it; that would be a
  different `shr_round` contract in `fft_arith.h`, pinned by its own numbers,
  and is not this kernel's decision to make.

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
  instruction-count tables from the QEMU legs.
- The analysis instruments typed over `std::span<const Sample>`, the capi and
  `dsptap_py` exposure of the four profiles, and `notebooks/fft.ipynb`
  executed against them, which re-measures the table in §5 on the shipping
  C++ and adds BFP exponent behaviour on speech-like material.
