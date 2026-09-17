# The real FFT: design note

*Skeleton, September 2026. This is the design note that is too long for a
header — the role `kaiser.h`'s design note plays for the FIR substrate. It is
filled in stage by stage; every entry that depends on a measurement is marked
`TODO(stage …)` and names the stage whose PR supplies the number. Until every
stage has landed, the plan of record is
[`audit-fft-and-code-smells.md`](audit-fft-and-code-smells.md) (parts are
cited as "audit Part N" below), and each PR links its stage.*

## Purpose and scope

`tap::dsp::basic_real_fft<Sample>` (`include/tap/dsp/fft.h`) is one real-FFT
contract over four sample profiles — `double` (the golden model), `float` (the
embedded profile, optionally accelerated by vDSP or CMSIS-Helium), and the
Q15 / Q31 fixed-point profiles from `sample_traits.h`. The contract is the
packing, the sign convention, the scaling, and the numbers in the table below;
the engine behind it changes (vendored C today, a C++20 port of the same
split-radix algorithm next, a radix-4 int32 kernel for fixed point) and the
contract does not. Consumers hold the transform by value and depend on
exactly the surface audit Part 1 enumerates: `basic_real_fft(size_t)`,
`forward_inplace`, `inverse_inplace`, `forward`, `inverse` (out-of-place,
`2/N` applied), `size()` / `num_bins()`, copyability, the packing, `exp(+i)`,
the unnormalized in-place inverse, and the aliases `real_fft` / `real_fft32`.

Out of scope here: the complex transform (`cdft`), the DCT/DST family, and any
consumer-side spectrum arithmetic (that is the Stage 5 packed-spectrum view).

## Contract summary, per profile

The header owns these numbers and the tests pin them; this table is the
cross-reference. Floating-point rows are the shipping `fft.h` at `5ca3b1c`;
fixed-point rows are the Part 7 design and fill in at Stage 3b.

| Contract point | `double` (`real_fft`) | `float` (`real_fft32`) | Q15 (`real_fft_q15`) | Q31 (`real_fft_q31`) |
|---|---|---|---|---|
| Packing | `data[0]` = DC, `data[1]` = Nyquist, `data[2k]` / `data[2k+1]` = bin *k* re / im, `1 ≤ k < N/2`; `N/2 + 1` bins | same | same (`TODO(stage 3b)` pin) | same (`TODO(stage 3b)` pin) |
| Sign | `W = exp(+2πi/N)`; imaginary parts conjugated vs the engineering DFT | same | same (`TODO(stage 3b)` pin) | same (`TODO(stage 3b)` pin) |
| Forward scale | 1 (`A[k] = Σ a[j] W^(jk)`, unnormalized) | 1 | `fixed`: exactly `X / N` in Q15; `block_floating`: `X · 2^-e`, `e` returned — `TODO(stage 3b)` | `fixed`: exactly `X / N` in Q31; `block_floating`: as Q15 — `TODO(stage 3b)` |
| Inverse scale | `inverse_inplace` unnormalized (caller applies `2/N`); `inverse()` applies `2/N` | same | own per-stage halving; round-trip factor a fixed power of two — `TODO(stage 3b)` | same — `TODO(stage 3b)` |
| Round-trip identity | `x` reproduced to `1e-12` abs at N = 1024 (`RoundTripReproducesInput`) | `2e-5` abs at N = 1024 (`RoundTripReproducesInput`) | `TODO(stage 3b)` per scaling policy | `TODO(stage 3b)` per scaling policy |
| Saturation-free input | n/a (floating point) | n/a | full scale ±1: input placed with 2 guard bits (`<< 14`, not 16) — `TODO(stage 3b)` measured worst case | `fixed`: one input pre-shift (−6 dB); `block_floating`: full scale — `TODO(stage 3b)` |
| Noise floor | `TODO(stage 2a)` against the compensated-DFT oracle | rms relative error vs double `< 1e-6` at N = 1024 (`FloatTracksDouble`); measured 1.105e-7 at N = 512 (audit Part 4) | `TODO(stage 3b)` per-bin floor vs level at N = 256 / 512 / 2048, vs the Welch model | `TODO(stage 3b)` |
| Latency | 0 (block transform, no internal delay) | 0 | 0 | 0 |
| Alignment | none required on `Sample*` | none (vDSP's internal split buffers are placed by the wrapper, not the caller) | none | none |
| Shareability across threads | not today: Ooura tables are built lazily on the first transform (audit F6); `is_shareable` engine trait true for the port — `TODO(stage 2a, 4)` | not today: the vDSP / CMSIS engines carry scratch; false for those engines — `TODO(stage 4)` | `TODO(stage 3b)` (design: true) | `TODO(stage 3b)` (design: true) |
| Real-time safety | transforms `noexcept`, allocation-free (the float-I/O-on-double overloads are the exception and are slated for removal, D5) | same | `TODO(stage 3b)` | `TODO(stage 3b)` |
| NaN / denormals | NaN propagates to every bin; denormal input is slow on x86 without FTZ, and FTZ differs between the Ooura, vDSP and CMSIS builds | same | not applicable / none | not applicable / none |

## Why split-radix for floating point and a radix-4 int32 kernel for fixed point

**Floating point stays on Ooura's split-radix algorithm**, transliterated
into C++20 (Stage 2a). The reasons are contractual, not aesthetic: every
consumer float pin in MuTap and DspTap was measured on the vendored C, and the
port reproduces the C's table semantics for both precisions so that `float`
stays *bit-identical* (D10). The alternative — computing the float twiddles in
double "because it is more accurate" — was measured and is a wash (table max
absolute error 1.19e-7 either way; transform rms relative error at N = 512 is
1.105e-7 with the C's tables and 1.114e-7 with double-computed ones), so it
was dropped. The `#define double float` build of the C is replaced by a
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
- A. V. Oppenheim and R. W. Schafer, *Discrete-Time Signal Processing*, the
  finite-precision FFT and block-floating-point sections.
- Ooura's own references for the split-radix algorithm are listed in
  `third_party/ooura/readme.txt` (Nussbaumer 1982; Burrus, *Notes on the FFT*).

### Derivations (`TODO(stage 3b)`)

- Halving count per policy and the resulting forward / inverse exponents.
- Headroom: the 4√2 per-stage bound, the packed-pair pattern at a 45° twiddle,
  and why two guard bits suffice for Q15 under fixed scaling.
- Twiddle quantization: `max |w_q − w| ≤ 0.5 LSB` of Q1.30 and its
  contribution to the floor.
- The Welch noise model per profile and policy against the measured floors at
  N = 256 / 512 / 2048 for white noise and an on-bin tone, from the executed
  `notebooks/fft.ipynb` (Stage 3c).

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
   "1 ulp" is an assumption, not a number.
   `TODO(stage 2a)`: the measured max-ulp deviation at default flags, per host
   (linux / windows / macos) and per QEMU leg, with compiler versions.
3. `fft.h` states whether `tap::dsp` exports `-ffp-contract=off` as an
   INTERFACE compile option (bit reproducibility across compilers) or leaves it
   to the consumer (VFMA speed on the M55). Today it is silently "whatever the
   consumer does". `TODO(stage 1)`: the decision and its wording.

A related, separate fact: libm `cos` / `sin` differ in the last bit between
glibc, newlib, UCRT and Apple, so *double* outputs are not identical across
hosts today either; a platform-independent table is a Stage 3a candidate and
otherwise out of scope.

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

## Size and instruction counts, per target

Seeded by the Stage 1b ratchet from the vendored C, then re-measured at 2b
(port), 3b (fixed point) and 4 (engine parameter). Numbers here carry the
SHA, toolchain and QEMU versions; the gate itself lives in
`bench/baselines.json` and the CI job, never in this file.

### `.text` per profile at N = 512 (`TODO(stage 1b, 2b, 3b)`)

| Target | float, C | float, port | Q15 | Q31 | double (host-class only) | SHA / toolchain |
|---|---|---|---|---|---|---|
| `m4-softfp` | | | | | n/a | |
| `m4f` | | | | | n/a | |
| `m33` | | | | | n/a | |
| `m55` (CMSIS on) | | | | | n/a | |
| `m55-ooura` | | | | | n/a | |

Pre-ratchet reference points from the audit (thumbv8.1m, hard float,
rdft-reachable text, `--gc-sections`): float 15.7 KB at `-Os`, 19.9 KB at
`-O2`; double (soft-float) 39.9 KB at `-O2`.

### Instructions per forward + inverse, per scenario (`TODO(stage 1b, 2b, 3b)`)

| Scenario | `m4-softfp` | `m4f` | `m33` | `m55` | `m55-ooura` | Baseline SHA |
|---|---|---|---|---|---|---|
| `rfft_f32_512` | | | | | | |
| `rfft_f32_2048` | | | | | | |
| `rfft_f64_512` (host-class only) | n/a | n/a | n/a | n/a | n/a | |
| `rfft_q15_512` | | | | | | |
| `rfft_q31_512` | | | | | | |
| `rfft_q31_2048` | | | | | | |

### Host microbenchmark (`bench/bench_fft.cpp`, informational)

`TODO(stage 2b)`: port vs C, min-of-N wall clock, with machine and date.
`TODO(stage 4)`: Ooura-port vs vDSP same-binary on Apple Silicon. The existing
"~3× faster / ~3× fewer instructions vs autovectorized Ooura" claims in
`fft.h` and the README are the consumers' measurements on the vendored C and
are re-measured here, not carried forward on trust.

## Provenance and licensing

### Where the code came from

The vendored `fftsg.c` was carried byte-for-byte inside MuTap and AmbiTap
before DspTap consolidated the wrappers; the FFT's Tap lineage is the
AmbiTap/MuTap one. Provenance stays visible after the port through four
things: the port header's attribution banner, `NOTICE.md`, the parity test
against `tests/reference/ooura/` (where `fftsg.c` and `readme.txt` move at
Stage 2c), and one glossary line in MuTap's `docs/itu-compliance.md` mapping
"measured on Ooura" to "the vendored C at DspTap ≤ `5ca3b1c` and the
bit-identical port from the Stage 2b SHA onward". Consumer *code* comments drop
"Ooura packing" in favour of the numeric definition in the Stage 5 view. The
engine is named for what it is (`detail::split_radix_rdft`, D7), not for its
author; attribution is carried by the banner and the notices, not the
identifier.

### The licensing statement

`NOTICE.md` carries the same statement; when one changes, the other changes in the same PR.

Ooura's terms, stated in full in `third_party/ooura/readme.txt` and in the
`fftsg.c` banner, are:

> You may use, copy, modify this code for any purpose and without fee. You
> may distribute this ORIGINAL package.

What that grants: use, copying and modification of the code, for any purpose
and without fee; and distribution of the *original* package. Distribution of a
modified derivative is not expressly granted. Today DspTap distributes the
original `fftsg.c` unmodified (`fftsg_float.c` is a 62-line DspTap wrapper
that `#include`s it under a type remap), which the distribution grant covers.

Going forward, the C++20 port of `rdft` is a **derivative work, not the
ORIGINAL package**, and its redistribution relies on the **modification grant**
("modify this code for any purpose") together with the widespread practice of
shipping modified Ooura code with the notice retained (WebRTC / Chromium under
`common_audio/third_party/ooura/`, with a local-modifications README) —
practice plus the author's informally reported lack of objection, not license
text, and none of those projects relicense the derived code. Neither does
DspTap: the port header carries Ooura's notice verbatim as the governing terms
for the derived portion, with SPDX `LicenseRef-Ooura AND MIT` (MIT only for the
wrapper and DspTap's additions), and a line stating that it is a derivative
work with the modifications copyright and date. `third_party/ooura/readme.txt`
stays in-tree permanently as the license record for the derived code, after
`fftsg.c` itself has moved to the test reference tree and been retired (D6).
"permissive" and "public" — the words the earlier notice used — overstated
the grant and are not used.

This is a maintainer judgement call, not legal advice.

### The port header's banner, ready to paste

```cpp
/// @file split_radix.h
/// @brief Split-radix real FFT engine: a C++20 transliteration of Ooura's rdft.
// SPDX-License-Identifier: LicenseRef-Ooura AND MIT
// Copyright(C) 1996-2001 Takuya OOURA
//     email: ooura@mmm.t.u-tokyo.ac.jp
//     download: http://momonga.t.u-tokyo.ac.jp/~ooura/fft.html
//     You may use, copy, modify this code for any purpose and
//     without fee. You may distribute this ORIGINAL package.
// Modifications Copyright 2026 Timothy Place and the DspTap contributors;
// this is a derivative work, not the ORIGINAL package. The derived portion
// (the transform itself) is governed by Ooura's notice above; the wrapper and
// DspTap's additions are MIT (see LICENSE). Full statement in NOTICE.md and
// docs/fft-design.md; the original package's readme.txt is kept in-tree at
// third_party/ooura/readme.txt.
```

`LicenseRef-Ooura` is a project-local SPDX license reference whose text is
`third_party/ooura/readme.txt`'s copyright section; that is why the readme
stays in-tree.

### Draft email to the address in Ooura's notice — for the maintainer to send

Per audit Part 5, the maintainer attempts the address in the notice before
Stage 2c merges and records the outcome here either way. Nothing has been
sent; this is the draft. Whether to send it is the maintainer's decision. The
address in the notice is `ooura@mmm.t.u-tokyo.ac.jp` and dates from
1996–2001, so a bounce is a likely outcome and is itself the record.

> **Subject:** Permission question: distributing a modified derivative of the
> General Purpose FFT Package (fftsg.c)
>
> Dear Professor Ooura,
>
> Thank you for the General Purpose FFT Package. I maintain DspTap, a small
> open-source (MIT) library of DSP primitives for audio, at
> https://github.com/tap/DspTap. It has shipped your `fftsg.c` unmodified for
> some time, with your readme and notice intact.
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
> The derived code would remain under your terms (it would not be relicensed),
> and would be clearly marked as a modified derivative with a link back to
> your package. A one-line reply either way would be very welcome, and I will
> record it in the project's notices.
>
> With thanks and best regards,
>
> Timothy Place
> https://github.com/tap/DspTap

**Outcome:** `TODO(stage 2c)` — sent on / not sent because; reply / bounce /
no reply by date.

## Changelog of contract-affecting SHAs

Every SHA that moves a documented contract point, per profile, with the
consumer pins that were re-measured. Empty until the first such change; the
port (Stage 2b) is designed to add *no* row here for `double` or `float`.

| SHA | Stage | Profile(s) | What moved | Consumer pins re-measured |
|---|---|---|---|---|
| | | | | |
