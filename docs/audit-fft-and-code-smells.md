# Audit: the FFT layer and its relatives across DspTap

*September 2026. Baseline at the time of the audit: `5ca3b1c`, builds warning-free with
`-DTAP_DSP_WERROR=ON`, 160/160 tests pass on Linux/Ooura. Revision 2: Parts 3-5 rewritten
after the adversarial review recorded in Part 6.*

The trigger was one smell: the float Ooura build is produced by `#define double float` plus
forty symbol-renaming `#define`s in `third_party/ooura/fftsg_float.c`. This document records
what else in the repo belongs to the same family, ranks it, and lays out a staged plan whose end
state is one FFT implementation that supports `double`, `float`, and the fixed-point profiles
from `sample_traits.h`. The *shape* of that FFT is left as explicit decision points at the end.

---

## Part 1 — The FFT layer

### F1. `#define double float` and the forty renames (`third_party/ooura/fftsg_float.c`) — the complaint
Textual retargeting of a 3325-line C file. It works only because `<math.h>` is included first,
depends on every file-scope symbol being enumerated by hand, and cannot be extended to a third
type. It is the root of F2, F3 and F7.

### F2. Ooura leaks 38 externally-visible functions into the global C namespace, twice
`fftsg.c` declares none of its helpers `static` (`makewt`, `bitrv2`, `cftf161`, ...); with the
`_f` copies that is 76 global symbols per binary. The rename table in F1 exists *because* of
this. Actual exposure, stated precisely: Max externals are MODULE bundles bound per image
(macOS two-level namespace; Windows exports nothing without `__declspec`), so the realistic
collision is a Linux shared object such as the capi, whose C objects are not visibility-hidden.
Real, but not urgent; the port removes it as a side effect.

### F3. A "header-only" library that requires a compiled static library
README and CLAUDE.md say header-only; in fact `tap::dsp` links `tap_dsp_fft`, and every
consumer CMake (MuTap, MuTap-Max per-external) carries comments explaining the extra target.
Of the 3325 lines compiled twice, the wrappers call one function (`rdft`); `ddct/ddst/dfct/dfst`
and the pthread/Win32 threading scaffolding are compiled and exported for nothing.

### F4. Backend selection by preprocessor inside the class body (`fft.h`)
`TAP_DSP_FFT_FLOAT_BACKEND` is derived from two build defines and then consulted in four
places inside `basic_real_fft`: constructor, both transforms, and the member list. Each site is
`#if ... if constexpr (std::is_same_v<Sample, float>) { ...; return; } #endif` followed by the
Ooura path. Consequences:
- The **object layout of `basic_real_fft<float>` depends on a build define** (`m_engine` exists
  or not; `m_ip/m_w` are allocated or not). The defines are INTERFACE compile definitions, so
  every TU inside one CMake build agrees; the hazard is **cross-image**: two separately built
  images that both export the weak template symbols and land in one process (on macOS dyld
  coalesces weak definitions across images unless visibility is hidden). Not observed today
  (MuTap-Max forces vDSP off through MuTap, the capi is never loaded into Max, AmbiTap keeps
  its own wrapper), real in principle, and the same hazard exists one level up for every class
  that embeds the FFT by value (`basic_pvoc<float>`, MuTap's `fdaf<float>`, ...).
- `fft_engine_noop` and `float_engine_t` exist only to make the `#if` compile for `double`.
- There is no way to instantiate Ooura and vDSP in the *same* binary, so backend parity is a
  CI-matrix property rather than a test.

### F5. Type genericity by enumeration
`static_assert(is_same_v<Sample,float> || is_same_v<Sample,double>)` at `fft.h:334`, and the
same gate in `yin.h:66`, `psola.h:60`, `pvoc.h:67`, `log_mel.h:156`. `decimate.h` already
shows the house alternative (`template <sample_type S>`). Any fixed-point FFT hits these walls
in pvoc and log_mel first.

### F6. Construction does not construct: Ooura tables are built lazily on the first transform
`m_ip[0] = 0; // triggers Ooura table init on first call`. `rdft` then runs `makewt`/`makect`
(a `sin`/`cos` per twiddle) inside the first `forward_inplace`, i.e. inside the `noexcept`,
"allocation-free, real-time safe" path. No allocation, but a first-call cost spike on the audio
thread, and a contract the header does not state.

### F7. The float-I/O-on-double overloads allocate per call
`forward(const float*, float*)` and `inverse(...)` on the double engine build a `std::vector`
on every call, in a class whose docstring leads with real-time safety. No consumer in DspTap,
MuTap or MuTap-Max calls them (the stated consumer is AmbiTap). They are a cast the caller can
do; they should not be members of the RT primitive.

### F8. The spectrum packing is re-derived at every consumer
Bin arithmetic (`data[0]`=DC, `data[1]`=Nyquist, `data[2k]/[2k+1]`) is hand-written in
`pvoc.h` (four sites, including the sign-convention flip), `log_mel.h`, `test_pvoc.cpp`, and in
MuTap: 25 sites in `fdaf.h`, 28 in `fd_kalman.h`, 15 in `nn_suppressor.h`, 9 in `postfilter.h`.
Each also recomputes `N/2+1` instead of calling `num_bins()`. There is no shared packed-bin
view. This is the single seam a redesigned FFT most needs, and it should land *before* the FFT
changes.

### F9. Documentation debt inside `fft.h`
The vDSP alignment story is told three times; the paragraph beginning "Computed at use rather
than cached" appears verbatim twice (`fft.h:232`, `:271`). The class comment, the member
comment and the test file each carry the same measured table. One home, cross-referenced.

### F10. Stale names and dangling references in the build files
- `NOTICE.md` refers to `include/dsptap/` and `DSPTAP_FFT_CMSIS`; the real names are
  `include/tap/dsp/` and `TAP_DSP_FFT_CMSIS`.
- `cmake/arm-cortex-m55-mps3.cmake` links `../platform/mps3_an547.ld` and
  `armv8m_startup.c`, neither of which exists in this repo, and sets `MUTAP_BARE_METAL`, a
  MuTap variable. The CI leg builds only the static library, so the link line is never
  exercised and the file cannot actually produce a test binary here.
- MuTap's bare-metal gtest filter still names `real_fft_test/*` and `fft_backend_parity`,
  tests that moved here (MuTap-side follow-up).

### What the consumers actually depend on (the surface the redesign must keep)
Across pvoc, log_mel, the tests, MuTap (fdaf, fd_kalman, pem_afc, postfilter, nn_suppressor)
and MuTap's test support: the constructor `basic_real_fft(size_t)`, `forward_inplace`,
`inverse_inplace`, `forward`, `inverse` (out-of-place, `2/N` applied), `size()`/`num_bins()`,
copyability (held by value in every consumer), the packing, the `exp(+i)` convention (pvoc is
the only consumer that would break if it flipped), the unnormalized in-place inverse (pvoc folds
`2/N` into its COLA gain), and the aliases `real_fft`/`real_fft32`. Nothing consumes the
float-I/O-on-double overloads or the raw `cdft` declarations.

---

## Part 2 — Same-family smells elsewhere in the repo

Ranked by severity. Every item was read in the source; line numbers are at `5ca3b1c`.

### High
- **`long m_n` sample counters** (`psola.h:189`, `pvoc.h:365`). `long` is 32-bit on Cortex-M
  and on Windows; at 48 kHz it overflows (signed UB) after ~12.4 hours, the subsequent
  `m_n % size` goes negative, and the `static_cast<size_t>` makes it an out-of-bounds write.
  The embedded profile is exactly where `long` is 32 bits. Fixed-width type, or wrap modulo
  the ring explicitly.
- **`sample_traits` has no `double` specialization, and the test pins that as intentional**
  (`test_sample_traits.cpp:108`, `static_assert(!sample_type<double>)`). CLAUDE.md says double
  is the golden model; the one primitive built on the traits (`decimate`) therefore has a
  *float* golden model and no double instantiation. A traits-based FFT cannot instantiate its
  golden profile until this contradiction is resolved one way or the other.
- **`fir_kernels.h` is the `fft.h` pattern again**: a feature macro (`TAP_DSP_Q15_SMLALD`)
  wrapping `if constexpr (is_same_v<S, int16_t>)` inside the body of a
  `template <sample_type S>`. The SMLALD path is never compiled on a host, and contains an
  unbraced loop body the style gate cannot see (`:98`, `:147`). `TAP_DSP_CHANNEL_PARALLEL` and
  `TAP_DSP_CP_MIN_CHANNELS` are defined and exported but read by nothing in DspTap or MuTap;
  README says they select a layout.
- **`kaiser.h:135` `solve_dense` is `noexcept` and allocates a `std::vector`**: `bad_alloc`
  becomes `std::terminate`. Its sibling `design_prototype_compensated` gets this right.
- **pvoc's float profile is a double profile with a float FFT.** `m_mag`, `m_true_bin`,
  `m_env` are `std::vector<double>`; `atan2`/`sqrt`/`round` and the region rotate run in
  double per bin. Only the Levinson recursion is documented as double. On an FP32-only M55
  this soft-floats every bin around a Helium transform; README's "the float profile rides the
  vDSP/CMSIS backends" oversells it.

### Medium
- **`sample_traits` is a FIR-dot-product trait, not an arithmetic customization point.** One
  operation (`mac` + `finalize`), no sample<->double conversion, no add/sub/mul-by-coefficient
  returning a sample, no shift/rounding helper, and the Q ladder is bare literals (`16384.0`,
  `(1 << 13)) >> 14`, `>> 16`) rather than named fraction-bit constants. Tests bypass it with
  hand-rolled `sample_from<>` casts. The `sample_type` concept omits `k_coeff_scale`, which
  `quantize.h` requires.
- **STYLE.md prescribes `TAP_EXPECTS`/`TAP_ENSURES`; the repo has zero uses.** Every
  precondition is a bare `assert` (`fft.h:340`, `yin.h:78`, `pvoc.h:84`, `psola.h:72`,
  `nn.h:92`, `decimate.h:127`); in release, bad geometry is silent UB (nn: out-of-bounds read
  on a wrong weight size). The capi re-implements the yin/psola/pvoc checks by hand because the
  headers offer nothing callable.
- **Duplication with no shared home**: `k_pi` re-literalled in `psola.h`, `pvoc.h`,
  `test_pvoc.cpp`, `test_fft_backend.cpp` while other headers use `std::numbers::pi`; periodic
  Hann hand-rolled in `pvoc.h`, `log_mel.h`, `test_pvoc.cpp`; dB conversion inline in four
  files; three Gaussian-elimination solvers (`kaiser.h`, `sine_analysis.h` 3x3, and the
  duplicated tone tracker between `sine_analysis.h` and `multitone_analysis.h`); `run_sine`,
  `measure_hz`, `cents` copy-pasted verbatim between `test_psola.cpp` and `test_pvoc.cpp`; two
  different generators both named `xorshift32` (`test_log_mel.cpp`, `test_nn.cpp`).
- **`decimate_profile::taps<M>()`** is an `if constexpr` chain whose `else` swallows any `M`;
  the charter `static_assert` lives in a different struct.
- **`quantize.h:49`** picks its algorithm by `is_floating_point_v` rather than a trait
  property, and its `+1/-1` correction bypasses the saturating `make_coeff`.
- **`nn.h`** hardcodes `std::vector<float>` weights whatever `Sample` is, and converts per
  element inside the inner loop.
- **Integer vocabulary**: `int` geometry in yin/pvoc/fft, `long` in psola/pvoc, `std::size_t`
  in log_mel/decimate/nn; unqualified `size_t` in half the headers. Every public accessor is a
  `static_cast<size_t>(m_int)` round trip.
- **`sine_analysis.h`/`multitone_analysis.h`** take `std::span<const float>` only, so the
  instruments cannot score the double golden model or a Q15 output without a copy;
  `program_weighted_snr_db` carries a commented-out unused parameter in its public signature;
  `solve_dense` lives in `kaiser.h` and the analysis headers include the filter-design header
  to get it.
- **capi**: `void*` handles (no type safety at the C boundary), exceptions can cross
  `extern "C"` (`new`, `make_unique`, `resize` in `process`), and every entry point is
  double-only, so the notebooks cannot measure the float or fixed-point profiles at all.

### Low (recorded so they are not rediscovered)
- `kaiser.h`: `\param` docs with pre-rename camelCase names (`attenDb`, `numPhases`, ...);
  `for (int pass = 0; pass < 1; ++pass)`; `std::max`/`std::swap` without `<algorithm>`/
  `<utility>`; unnamed tuning literals.
- `psola.h`: integer geometry stored as `Sample` and recovered by cast; a `std::cos` per
  output sample per grain and four double-modulo ring indexes per sample in the hot loop.
- `pvoc.h`: `m_lpc_time = m_frame` (vector copy-assignment in the noexcept frame path, safe
  only because sizes match); `push_back` in the per-frame path relying on an unrelated
  `reserve`; default geometry argument on a primitive whose siblings require it; magic
  `1e-4`, `3`, `±2` in the peak picker.
- `log_mel.h`: `band_weight(b, k)` with no bounds check on `b`; geometry held as `double` and
  re-cast to `Sample` every frame; `memmove` of `frame - hop` per hop.
- `sample_traits.h`: `round_sat` on NaN is UB; `clamp_sat` hardcodes `int64`; nothing is
  `constexpr`; unused `<cmath>`.
- Tests: `DotRowMatchesReferenceAccumulation` is tautological (the reference is the same
  `mac` loop); no kernel-level cross-precision test for Q15/Q31; per-type tolerances expressed
  three different ways across files.
- README: `basic_decimator<Sample, M>` vs header's `S`; the `nn` example copies where the
  header advertises move-in.

---

## Part 3 — Plan (revision 2)

The stated goal is one real-FFT implementation over `double`, `float` and, when a consumer
needs it, the fixed-point profiles. Revision 1 mixed that goal with unrelated hygiene, made the
float numbers move in the same PR as the port, promised M55 gates the repo cannot run, and
carried a fixed-point design that does not survive arithmetic (Part 6). Revision 2 is the
minimum path to "no `#define double float`, no global symbols, header-only, one template",
with everything else either sequenced behind it or moved to Appendix A.

**Rollout rule for every stage.** DspTap PR, squash to `main`; MuTap PR that bumps the pin *and*
carries that stage's MuTap-side change in one PR; MuTap-Max pin. Each stage names which MuTap
tests must be **unchanged** (the fingerprint harness, bit-identical) and which are
**re-measured**. DspTap has no tags; "released" means both consumers pin a tree containing it.

### Stage 0 — The counter bug, alone
`long m_n` in `psola.h:189` and `pvoc.h:365` (32-bit on Cortex-M and Windows; signed overflow
after ~12.4 h at 48 kHz, then an out-of-bounds write). Two files, a fixed-width counter or an
explicit modulo wrap, and a test that seeds the counter near 2^31 and crosses it. Bump MuTap
immediately. Nothing else rides in this PR.

### Stage 1 — Consumer and CI prerequisites (MuTap PR, DspTap PR)
- MuTap: `ci.yml:176` passes `-DMUTAP_FFT_CMSIS=OFF`; the option is `TAP_DSP_FFT_CMSIS`, so the
  "Ooura fallback on M55" leg has been rebuilding CMSIS. Fix the flag. Prune the bare-metal
  gtest filter of the `real_fft_test/*` and `fft_backend_parity` names that moved here.
  Generalize `tests/branchless_parity_check.cpp` into a fingerprint harness over the fdaf,
  fd_kalman, pem_afc, postfilter and nn_suppressor chain outputs, runnable at any pin; this is
  the bit-identity gate every later MuTap bump uses.
- DspTap: make the M55 leg honest. The toolchain file links `../platform/*` files that live in
  MuTap, and CI builds only the static library with tests off. Decision: strip the link flags,
  declare the toolchain compile-only, and add a compile-only object target that instantiates
  the engine for `float` and `double` (Stage 2a onward) so the M55 leg compiles the port. No
  sentence in this plan calls an M55 *runtime* property a DspTap gate; the on-target gates are
  MuTap's M33 and M55 QEMU legs, which run the float suppressor suite on the Ooura float path
  once the flag above is fixed. Porting the platform files here is deferred to Appendix A.
- DspTap: state the fp-contraction policy (Part 4) as a contract point in `fft.h` before the
  port lands, so the parity target and the consumers are measured against a written rule.

### Stage 2 — The port, in three PRs
**2a. Add the engine beside the C, route nothing.** `detail::split_radix_rdft<Sample>` per
Part 4, transliterated statement-for-statement, tables computed *exactly as the C computes
them for each precision* (float-typed locals for the float instantiation, libm calls with
explicit `static_cast<double>` arguments). Gate: a parity test TU that compiles both C files
and the C++ with `-ffp-contract=off` and requires **bit identity for double and for float**,
forward and inverse, at every power of two from 4 to 65536 plus one run at 2^20, on the
linux, windows and macos legs (same binary, same libm; cross-platform identity is not claimed
and does not hold today either). A second run at default flags is informational and pinned at
a *measured* bound. An independent oracle: closed-form vectors (impulse, DC, Nyquist, on-bin
tone) plus a compensated-summation DFT at small N, because `long double` is `double` on MSVC
and Apple arm64. A host microbenchmark, port vs C at N=512/2048 float/double, as a CI
artifact. M55 leg compiles the instantiations. MuTap: nothing moves; no bump needed.

**2b. Flip routing.** One line in `fft.h`; `test_fft_backend.cpp`'s `ooura_ref` re-pointed at the
ported float engine (on Linux/Windows the backend test becomes port-vs-port, which is fine
only because 2a's parity test exists); README rewritten in the same PR. MuTap bump: fingerprint
harness bit-identical (double rows *and*, because 2a kept float table semantics, the float
rows), icount ratchet at 0% delta on m33 and hexagon (the ratchet's ±3% is not the gate here;
0% is, because nothing numeric changed), `test_float32`, `test_g168`, `test_nn_suppressor`
unchanged. Rollback is a one-line revert.

**2c. Remove the C.** Delete `fftsg_float.c`; move `fftsg.c` and `readme.txt` to
`tests/reference/ooura/`; `tap_dsp_fft` exists only when `TAP_DSP_FFT_CMSIS` is on. MuTap:
rewrite the CI job that compiles `submodules/dsptap/third_party/ooura/fftsg*.c` by path
(`ci.yml:319-324`) header-only; rewrite `THIRD_PARTY_NOTICES.md`. Licensing text per Part 5.

### Stage 3 — Engine as an explicit parameter, with an ABI tag
`basic_real_fft<Sample, Engine = default_real_fft_engine_t<Sample>>`. The default alias lives
in one place, selected by the one remaining build define. To close F4 where it actually lives
(the consumers that embed the FFT by value), the selection also opens an inline namespace ABI
tag on `tap::dsp` (`inline namespace fft_ooura {}` / `fft_vdsp {}` / `fft_cmsis {}`), so two
images built with different defaults cannot coalesce each other's symbols. Backends move to
`fft/backends/accelerate.h` and `fft/backends/cmsis.h`, included by whoever selects them; the
CMSIS object library stays PIC. Consumer-facing transforms stay **non-const**; shareability is
an engine trait (`Engine::is_shareable`), true for the ported engine, false for the two
scratch-carrying backends. Typed tests over the engines available on the host: Ooura vs vDSP
same-binary on macOS; CMSIS compile-only on M55 (it cannot run on a host, and nobody runs
CMSIS-vs-Ooura parity anywhere today; that gap is recorded, not closed, by this stage).

### Stage 4 — The packed-spectrum view, DspTap only
`fft/spectrum.h`: a non-owning view whose docstring carries the numeric definition
(`bin[k] = a[2k] + i a[2k+1]`, DC at `a[0]`, Nyquist at `a[1]`, `W = exp(+2πi/N)`, inverse
unnormalized). Primary accessors are **native**: `dc()`, `nyquist()`, `re(k)`, `im(k)`,
`power(k)`, `num_bins()`. A convention-flipping accessor, if kept, is named unmistakably
(`bin_engineering(k)`) and is not the default. Migrate `pvoc.h`, `log_mel.h`, `test_pvoc.cpp`,
gated by the existing pinned tests. MuTap adopts the view per header when each is next
touched, each such PR gated by the fingerprint harness and 0% icount, because those 77 sites
do native-convention complex products by hand in ratcheted hot loops and a wholesale rewrite
is not the mechanical migration revision 1 called it. This stage is **not** a prerequisite for
Stage 2: the packing is a kept contract.

### Stage 5 — Hygiene (after the port, not before)
`detail/math.h` (pi via `std::numbers`, periodic Hann, dB helpers), `tests/support/` for the
copy-pasted helpers, `solve_dense` noexcept, kaiser includes and doc drift, NOTICE names,
the two duplicated `fft.h` paragraphs, the unbraced SMLALD loops. Gate the Hann/pi
consolidation with the fingerprint harness: association order must be preserved, and
log_mel's numpy pin only bites above 1e-6. Do **not** delete `TAP_DSP_CHANNEL_PARALLEL` /
`TAP_DSP_CP_MIN_CHANNELS` until SampleRateTap and RatioTap (not on disk) have been grepped;
README says they are consumed there. `TAP_EXPECTS` lands here for `fft.h`'s power-of-two
precondition only; a repo-wide precondition policy is its own plan.

### Gates, in one table

| Stage | DspTap gate | MuTap gate | Rollback |
|---|---|---|---|
| 0 | new overflow test | pin bump, suite green | revert 2 files |
| 1 | M55 leg compiles instantiations | fallback leg builds Ooura float; filter pruned | n/a |
| 2a | bit identity double+float, 3 hosts, `-ffp-contract=off`; oracle; bench artifact | none (no bump) | delete header |
| 2b | existing battery on the port; backend test re-pointed | fingerprint identical; icount 0%; float pins unchanged | one-line revert |
| 2c | build without the C | CI job rewritten; notices | re-pin |
| 3 | typed engine tests; macOS same-binary parity | fingerprint identical | re-pin |
| 4 | pinned pvoc/log_mel tests | per-header fingerprint + icount 0% | per header |
| 5 | fingerprint on Hann/pi change | pin bump | per item |

---

## Part 4 — The C++20 port in detail (revision 2)

### Scope
`rdft` reaches about 2,580 of `fftsg.c`'s 3,325 lines: `makewt`, `makeipt`, `makect`,
`bitrv2`, `bitrv2conj`, `bitrv216`, `bitrv216neg`, `bitrv208`, `bitrv208neg`, `cftfsub`,
`cftbsub`, `cftf1st`, `cftb1st`, `cftrec4`, `cfttree`, `cftleaf`, `cftmdl1`, `cftmdl2`,
`cftfx41`, `cftf161`, `cftf162`, `cftf081`, `cftf082`, `cftf040`, `cftb040`, `cftx020`,
`rftfsub`, `rftbsub`. Not ported: the DCT/DST family, `cdft`, and the thread scaffolding.
`cftrec4` is a `while` plus a `for` over `cfttree`/`cftleaf`, not a recursion (revision 1 was
wrong); nothing in the reachable set self-recurses.

### Shape
- `template <std::floating_point Sample> class split_radix_rdft` in
  `include/tap/dsp/fft/split_radix.h`, namespace `tap::dsp::detail`. Every helper is a private
  static member function; the 76 global symbols disappear. The concept replaces the
  `static_assert` enumeration.
- Members: size, bit-reversal table, trig table (`std::vector<Sample>`). Index types may be
  widened to `std::size_t` (checked: every loop bound is a `< m` / `> 0` form over non-negative
  values); keep `-Wconversion` on to catch a missed cast. Load `m_w.data()` into a local once
  per transform so aliasing analysis matches the C, which receives `w` by parameter.
- **Every Ooura statement stays textually intact.** Clang contracts FMAs within a statement
  only; GCC across statements after inlining. Refactoring `wk1r * x0r - wk1i * x0i` into a
  `cmul` helper, a lambda, or a policy call changes which products fuse and silently breaks bit
  identity on one compiler or the other. This rule is stated in the header.
- **Every libm call takes an explicit `static_cast<double>`.** In C, `cos(delta * j)` with a
  float `delta` calls the double `cos`; in C++ `std::cos` of a float calls `cosf`, and the tables
  then differ from the C. The parity test at N where `makewt` takes the `nwh > 4` branch catches
  this, but only because 2a keeps the C's table semantics.
- Engine surface: constructor from size, `forward_inplace(Sample*) const`,
  `inverse_inplace(Sample*) const`, `size()`. Const is sound for the ported engine: the only
  mutable state in the C is the lazy `nw`/`nc` re-init in `rdft` and the thread code, both gone.
  `basic_real_fft` keeps the consumer-facing surface (Part 1), adds `std::span` overloads,
  `[[nodiscard]]` on size queries, `std::size_t` with one narrowing point, `TAP_EXPECTS` on the
  power-of-two precondition, and drops the raw `cdft` declarations.
- Tables are built in the constructor, once, so F6 closes. They are built **with the same
  precision semantics the C uses for each instantiation** (see next section).

### The float twiddles: measured, and the decision
Revision 1 claimed the float build's twiddles were degraded because `#define double float`
retargets `makewt`'s locals, and proposed computing them in double "because it is more
accurate", accepting that float would no longer be bit-identical. Measured (Part 6, item N2):
table max absolute error 1.19e-7 either way; transform rms relative error vs double at N=512
is 1.105e-7 with the C's tables and 1.114e-7 with double-computed tables (worse), N=16 worse,
N=65536 6% better. It is a wash. Decision: **the port reproduces the C's table semantics for
both precisions and float stays bit-identical.** That is what lets Stage 2b bump MuTap with
every float pin unchanged. Double-computed tables are dropped from the plan. A separate note
worth keeping: libm `cos`/`sin` differ in the last bit between glibc, newlib, UCRT and Apple,
which is why `test_log_mel.cpp` carries a per-platform double tolerance; a
platform-independent table (compensated evaluation or a checked-in generator) is the only way
to make *double* outputs identical across hosts, and is Appendix A material.

### fp-contraction policy (a contract point, not a build detail)
Measured: `gcc -std=c17` does not contract; `gcc -std=gnu17` does; `g++` contracts in both
`c++20` and `gnu++20`; clang contracts in every mode, statement-scoped. Neither DspTap nor
MuTap nor MuTap-Max sets `-ffp-contract` or `CMAKE_C_EXTENSIONS`, so the C oracle and the
port are contracted differently by default wherever the ISA has FMA (Apple arm64, M55 VFMA,
any x86 built with `-march`). Therefore: (a) the parity target compiles both sides with
`-ffp-contract=off` and that is the bit-identity gate; (b) the default-flags run is
informational with a measured bound, because a different fusion choice per stage accumulates
over log2 N stages and "1 ulp" is an assumption; (c) `fft.h` states whether `tap::dsp` exports
`-ffp-contract=off` as an INTERFACE option (bit reproducibility across compilers) or leaves it
to the consumer (VFMA speed on M55). Today it is silently "whatever the consumer does", and
MuTap's "bit-identical" double rows are protected only if this holds in MuTap's build.

### Header-only, and why the escape hatch is gone
Header-only, unconditionally for the ported engine. `tap::dsp` is a true INTERFACE target
*unless* `TAP_DSP_FFT_CMSIS` is on, in which case a PIC object library carries the CMSIS C;
there is no `install()` rule in the repo, so the `$<INSTALL_INTERFACE>` lines are dead and
"header-only" is a statement to `add_subdirectory` consumers only. Keep MuTap-Max's
per-external comments until a bump proves they can go.

The `extern template` escape hatch from revision 1 is dropped: [temp.explicit] exempts inline
functions from suppression, and members defined in the class body are inline, so the hatch
does nothing unless every heavy member is defined out-of-class without `inline` (verified on
g++ and clang++). It is also unnecessary: `fftsg.c` compiles in 0.6-0.8 s as C and 1.3 s as
C++ for the whole file; MuTap has about 20 TUs that reach `fft.h`, two instantiations each,
well under a minute of CPU. Code size on M55 (thumbv8.1m, hard float, rdft-reachable text):
float 15.7 KB at `-Os`, 19.9 KB at `-O2`; double (soft-float) 39.9 KB at `-O2`. Header-only is
neutral for size because the toolchain already uses `--gc-sections`. Add a `.text` assertion
for the float instantiation to the M55 leg; MuTap's pico2w job is the model.

### Performance
The backends' quoted gains are "vs autovectorized Ooura", i.e. the C compiled as C. A port
with member access and vector storage may vectorize differently (SLP over the straight-line
leaves is what makes Ooura fast). Before the C is deleted: the host microbenchmark from 2a,
MuTap's icount ratchet on m33 and hexagon with the port swapped in, and the M55 size
assertion. Thresholds are decided before 2b merges.

### Other contract points to write down
No alignment requirement on `Sample*` (Ooura needs none; a future MVE backend must not add
one quietly). NaN propagates to every bin (no data-dependent branches). Denormal inputs are
slow on x86 without FTZ, and FTZ differs between Ooura, vDSP and CMSIS builds, which can show
as cross-backend differences on near-silent input.

---

## Part 5 — Naming and provenance (revision 2)

**The license is narrower than revision 1 said.** `readme.txt:141-145` and the `fftsg.c` banner:
"You may use, copy, modify this code for any purpose and without fee. You may distribute this
ORIGINAL package." Modification is granted; distribution is granted for the *original*
package; distribution of a modified derivative is not expressly granted. Industry practice
(WebRTC/Chromium ship heavily modified Ooura under `common_audio/third_party/ooura/` with the
notice retained and a local-modifications README; Ooura is widely reported as having no
objection) is practice plus an informal statement, not license text, and none of those
projects relicense the derived code. Therefore:
- The port header's banner carries Ooura's notice verbatim as the governing terms for the
  derived portion, with SPDX `LicenseRef-Ooura AND MIT` (MIT only for the wrapper and the
  additions), and a line stating that this is a derivative work, not the original package,
  with the modifications copyright and date.
- `NOTICE.md` states that redistribution of the derivative relies on the modification grant,
  quotes the full notice, replaces "permissive/public" (an overstatement already), and keeps
  `readme.txt` in-tree after `fftsg.c` leaves the main tree.
- Before Stage 2c merges, attempt the address in the notice for an explicit statement on
  derivative distribution and record the outcome either way. This is a maintainer judgement
  call, not legal advice.
- MuTap's `THIRD_PARTY_NOTICES.md` is rewritten in the 2c bump, because the notice will then
  live in a header compiled into every external.

**Naming, settled with one correction.** The engine is `detail::split_radix_rdft` (D7).
`fftsg` is genuinely split-radix (`readme.txt:14`). Use one house token for everything
consumer-facing (`real_fft`, `real_fft32`, `fft/backends/accelerate_real_fft32.h`), not the
three spellings revision 1 had. Provenance stays visible: the attribution banner, `NOTICE.md`,
the parity test against `tests/reference/ooura/`, and one glossary line in MuTap's
`docs/itu-compliance.md`, whose certified float numbers are described as "measured on Ooura":
that line maps "Ooura" to "the vendored C at DspTap ≤ 5ca3b1c and the bit-identical port from
the Stage 2b SHA onward" rather than scrubbing the word. Consumer *code* comments drop "Ooura
packing" in favour of the Stage 4 view's numeric definition.

---

## Part 6 — Adversarial review record

Two independent hostile reviews (numerics/port, and process/ecosystem) plus a third pass by
the author, against the revision 1 text. Every item below was verified against source or by a
probe before being accepted; items the reviews raised that did not survive verification are
not listed. **P** = process review, **N** = numerics review, **A** = author.

### Blockers (all resolved in revision 2)
- **P1/N1 — Float numerics changed inside the port PR, with no float oracle left.** Revision 1
  simultaneously required float bit identity (Stage 3 gate) and declared float would differ
  "by design" (Part 4), then deleted `fftsg_float.c`. Downstream float pins that would have
  moved: MuTap `test_float32.cpp:92-191`, `test_g168.cpp`, the macOS repeat-until-fail rows,
  `test_nn_suppressor.cpp:270-305`, the Python parity float profile, the M33 leg, the icount
  baselines, and `docs/itu-compliance.md:655-668`. Resolution: 2a reproduces the C's table
  semantics; float is bit-identical; both C files stay as oracles until 2c.
- **N2 — The "more accurate twiddles" justification was measured and is a wash** (numbers in
  Part 4). Resolution: dropped.
- **N3 — "Same `-ffp-contract` setting" is necessary, not sufficient, and CMake does not give
  it.** GCC C++ contracts in ISO mode; clang contracts statement-scoped; a helper-function
  refactor changes fusion on clang. Resolution: `-ffp-contract=off` parity target, statement
  fidelity rule, policy stated in `fft.h`.
- **N4/P4 — "Bit-identical on all four CI legs" and every M55 gate were unrunnable.** The M55
  leg is compile-only with tests off; the toolchain links files that live in MuTap; libm
  differs per host so cross-host identity does not hold today. Resolution: same-binary
  identity on three hosts; M55 compile-only with size assertion; on-target gates are MuTap's.
- **N5/P17 — `test_fft_backend.cpp` uses `rdft_f` as its oracle and dies with
  `fftsg_float.c`.** Resolution: re-pointed in 2b.
- **P5 — MuTap's Ooura-on-M55 CI leg is a no-op** (`-DMUTAP_FFT_CMSIS=OFF` against an option
  named `TAP_DSP_FFT_CMSIS`), so "MuTap's compliance batteries are the second gate" was hollow
  for the profile the port changes. Resolution: Stage 1.
- **P2 — Stage 1's MuTap half was not mechanical.** The 77 sites compute native-convention
  products by hand in icount-ratcheted loops; an engineering-convention view is a
  sign-sensitive rewrite of five certified headers. Resolution: native accessors, DspTap-only
  migration, MuTap per header when touched, and the view is no longer a prerequisite.
- **P3 — Licensing overstated.** Resolution: Part 5.
- **N6/N7/N8 — The fixed-point design did not survive arithmetic.** Radix-4 output scaling
  overflows the twiddle multiply in int32 for Q15 (17-bit data × Q1.14) and in int64 for Q31;
  a 16-point leaf grows 4× twice before write-back; per-radix-4 /4 still saturates on a
  full-scale packed-pair pattern at a 45° twiddle; `rftfsub` and the DC step have gain up to 3
  and 2; Ooura's inverse has structural gain N/2 so "inverse unscaled" saturates on the first
  stage; and fixed 1/N scaling gives ~25 dB per-bin SNR at −40 dBFS and ~5 dB at −60 dBFS for
  log_mel's N=512, i.e. quiet speech becomes rounding noise. Also `cftf1st` derives half its
  twiddles at run time as `csc1 * (wd1r + w[k])`, a sum of two cosines up to 2.0 that does not
  fit Q1.14, so "the same table" is false for the first stage. Resolution: fixed point moves to
  Appendix A with these constraints as its entry conditions.

### Should-fix (all adopted)
- **P6/N15/A — Stage 2's coupling to D1 was artificial**; the engine is `std::floating_point`
  until fixed point. D1 also had the coefficient type wrong (`coeff` must be `double` or the
  double FFT cannot be bit-identical). Moved to Appendix A, corrected.
- **P7/A — The counter bug shipped inside a hygiene PR.** Now Stage 0, alone.
- **P8 — Deleting `TAP_DSP_CHANNEL_PARALLEL` without grepping SampleRateTap/RatioTap.** Held.
- **P9/N12/A — D4 did not close F4 for the embedding consumers.** ABI tag added; F4 reworded to
  cross-image.
- **P10 — F2 severity overstated for MODULE consumers.** Reworded.
- **P11 — "True INTERFACE target" was conditional and contradicted by the escape hatch; no
  `install()` exists.** Reworded; hatch dropped.
- **P12 — Fixed point is scope creep against every consumer's written plan** (MuTap
  `wake-word-plan.md:978-981`: float32 on every target including RP2350, Q15 front end only if
  a measured M33 count misses, "DspTap has no fixed-point FFT"). Appendix A.
- **P14 — Stage 3 was not rollback-able.** Split into 2a/2b/2c.
- **P15/N17/A — No performance, size or compile-time gate.** Added (Part 4).
- **P18 — A MuTap CI job compiles the vendored C by path.** Stage 2c.
- **N9 — `extern template` does nothing for in-class-defined members** (verified on both
  compilers). Dropped, with the measured compile cost that makes it unnecessary.
- **N10 — Header-only is size-neutral** because `--gc-sections` is already on. Stated, with
  measured sizes.
- **N11 — Const transforms would make the public API's constness depend on the selected
  engine.** Public methods stay non-const; `Engine::is_shareable`.
- **N13 — Runtime twiddle derivation in `cftf1st`.** Recorded under fixed point.
- **N14 — `std::cos(float)` resolves to `cosf`.** Explicit-cast rule in Part 4.
- **N16 — Downstream float pins enumerated** (DspTap `test_fft.cpp:37,209`,
  `test_fft_backend.cpp:90,249`, `test_log_mel.cpp:97,135,389` at 2× margin, `test_pvoc.cpp:81,183,209`;
  MuTap `test_nn_suppressor.cpp:122,231,305`, `test_float32.cpp:325`). Moot once float is
  bit-identical, kept as the list to re-measure if any future numeric change is made.
- **P23 — `TAP_EXPECTS` everywhere plus per-primitive `valid()` was unrelated to the FFT.**
  Shrunk to `fft.h`.
- **P24 — "One typed suite runs all three engines in one binary" was false** (CMSIS cannot
  run on hosts). Reworded.
- **P25/N25 — Shareable plans are an API change not in the plan.** Removed from the port;
  `is_shareable` trait is the only residue.
- **P20 — Three spellings for one thing in the naming.** One house token.
- **P21/P22 — D5 gated on a repo not on disk; D6 "one release" undefined.** `[[deprecated]]`
  one cycle; "released" defined by consumer pins.
- **N18-N22 — Index widening is safe; `long double` is not a better oracle everywhere; run
  parity once at 2^20; NaN/denormal/alignment stated as contract points.** All in Part 4.
- **A — Part 4 called `cftrec4` recursive and counted 2,400 lines.** Loop; 2,580.

### Not adopted, with reasons
- **P20's suggestion to keep Ooura's name in the detail identifier** (`detail::ooura_split_radix`).
  D7 was settled with the maintainer; provenance is carried by the banner, NOTICE, the parity
  test path and the itu-compliance glossary line instead.
- **P13's suggestion to make the spectrum view entirely optional.** Kept as Stage 4, DspTap
  only, because pvoc and log_mel duplicate the layout today and the view's docstring is where
  the numeric definition lives once "Ooura contract" leaves the consumer headers.

---

## Appendix A — Deferred: fixed point, and the substrate questions it drags in

Deferred until a consumer measures the need (MuTap's wake-word plan says float32 on every
target, RP2350 included). Entry conditions, so the work is not re-derived:
- **D1 (was Stage 2).** Is `double` a sample format? If yes, `sample_traits<double>` has
  `coeff = double`, `accum = double`; `test_sample_traits.cpp:108` flips; `test_decimate.cpp:192`
  becomes Q15-vs-double. Its FIR `mac`/`finalize` shape (int64 accumulate, one rounding) is
  not what an FFT needs; the FFT wants a second trait (`add`, `sub`, `mul_coeff`,
  `shift_round`, named fraction-bit constants), so "extend sample_traits" means "add a sibling".
- **D2. One structure or two.** The ported split-radix can be reused as a *data-flow graph*
  only: Q15 needs a shift before each radix-4 sub-stage or int64 products inside the 16-point
  leaves; Q31 needs shift-before-multiply or 32×32→high-half multiplies. Either way the
  operation order is not Ooura's. `cftf1st`'s runtime twiddle sums (up to 2.0) need a
  materialized first-stage table or int32 twiddle arithmetic. Block floating point is feasible
  (sub-blocks finish at different times, so per-block exponents plus a final normalization
  pass), just bookkeeping.
- **D3. Scales, both directions, as numbers.** Forward per-sub-stage input shifts plus one /2
  at the real post-pass; the inverse needs its own per-stage scaling because Ooura's inverse
  has structural gain N/2; the round-trip identity, the caller's required pre-shift and the
  saturation-free input level (one headroom bit, or an input contract of ≤ −3 dBFS) are all
  header numbers with a worst-case-pattern test. Per-bin noise floor pinned: fixed 1/N in Q15
  is unfit for log_mel at N=512 (≈25 dB per-bin SNR at −40 dBFS); the profile is Q31 data with
  Q15 at the I/O boundary, or BFP.
- **CMSIS q15/q31 as an engine** requires vendoring more of CMSIS (the subset here is f32 only)
  and adopting its size-dependent output scaling; it can only be tested once the house
  contract above exists.
- **On-target runtime gates in DspTap** require porting MuTap's three platform files and a
  one-shot gtest harness (`TAP_DSP_BARE_METAL`, modeled on MuTap `tests/CMakeLists.txt:7-19`).
- **Platform-independent double tables** (Part 4) if cross-host bit identity of the double
  golden model is ever wanted.

---

## Decisions (revision 2)

**D1, D2, D3** — moved to Appendix A; not on the port's critical path.

**D4. Engine as a template parameter with a build-selected default, plus an inline-namespace
ABI tag derived from the selection.** The default alone leaves the layout-by-define hazard in
every class that embeds the FFT by value; the tag closes it. Alternative kept on record: no
default, consumers name the engine.

**D5. Float-I/O-on-double overloads: `[[deprecated]]` for one consumer cycle, then deleted.**
Not gated on AmbiTap, which is not on disk and keeps its own wrapper per README.

**D6. `fftsg.c` and `readme.txt` move to `tests/reference/ooura/` at 2c and are deleted after
both MuTap and MuTap-Max pin a tree containing 2c.** `readme.txt` stays in-tree regardless
(Part 5).

**D7. Settled: `detail::split_radix_rdft`**, one house token (`real_fft`) for everything
consumer-facing, provenance per Part 5.

**D8. No explicit-instantiation escape hatch.** Verified ineffective for in-class definitions
and unnecessary at measured compile cost. Revisit only if MuTap's test build shows a
regression the numbers in Part 4 do not predict.

**D9 (new). fp-contraction is a written contract point.** Parity gate at `-ffp-contract=off`;
whether `tap::dsp` exports that flag to consumers is decided in Stage 1 and stated in `fft.h`.

**D10 (new). Float stays bit-identical to the vendored C.** The port reproduces the C's table
semantics for both precisions; no numeric change ships with the port.
