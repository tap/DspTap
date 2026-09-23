# Audit: the FFT layer and its relatives across DspTap

*September 2026. Baseline at the time of the audit: `5ca3b1c`, builds warning-free with
`-DTAP_DSP_WERROR=ON`, 160/160 tests pass on Linux/Ooura. Revision 2: Parts 3-5 rewritten
after the adversarial review recorded in Part 6. Revision 3: fixed point reinstated as
Stage 3 (Part 7); documentation, testing and embedded-CI plans added (Parts 8-10). Revision 3.1:
amendments from the wave-1 hostile reviews (Part 13).*

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

## Part 3 — Plan (revision 3)

The stated goal is one real-FFT contract over `double`, `float`, Q15 and Q31, with the
fixed-point profiles a committed deliverable: fixed-point projects not yet started need FFT
support, so they are on the main line, not in an appendix. Revision 1 mixed the port with
unrelated hygiene, made the float numbers move in the same PR as the port, promised M55 gates
the repo cannot run, and carried a fixed-point design that does not survive arithmetic
(Part 6). Revision 3 keeps revision 2's port sequence, puts fixed point back as Stage 3 with
the design in Part 7, and adds the documentation plan (Part 8), the test plan (Part 9) and the
embedded CI legs that make the fixed-point gates real (Part 10).

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
- DspTap: make the embedded legs real. The toolchain file links `../platform/*` files that
  live in MuTap, and CI builds only the static library with tests off. Decision (Part 10):
  bring the platform files in, add the bare-metal one-shot gtest harness, and run the
  emulation-sized selection under QEMU on Cortex-M4 (soft-float and M4F), Cortex-M33 and
  Cortex-M55. Landed in wave 1 (DspTap #17): all four legs green at first attempt, 14-98 s of
  emulation each; the CMSIS Helium parity suites ran for the first time anywhere on the M55
  leg. Fixed point's gates (Stage 3) run on those legs. The test selection on target must be a
  NEGATIVE filter (exclude the named slow double suites) so that new suites run by default;
  a positive filter silently drops every later suite (Part 13).
- DspTap: state the fp-contraction policy (Part 4) as a contract point in `fft.h` before the
  port lands, so the parity target and the consumers are measured against a written rule.

### Stage 1b — Benchmarks and the performance ratchet, baselined on the C
Revision 2 had no performance stage: a host microbenchmark "as an artifact", an M55 size
assertion, and an instruction-count idea marked optional. That is not a gate. This stage adds
the measurement infrastructure **before** the port so that Stage 2b is ratcheted against the
vendored C rather than against nothing. Design in Part 11.

- `bench/icount/`: one deterministic binary per scenario, compile-time selected (no argv on bare
  metal), xorshift input, a checksum on every output so nothing is dead-code-eliminated;
  scenarios `rfft_f32_512`, `rfft_f32_2048`, `rfft_f64_512` (host-class targets only), each
  a forward+inverse loop sized so the transforms dominate plan construction. The C and the port
  build as two binaries per scenario from Stage 2a on (`TAP_DSP_BENCH_ENGINE`), so the ratchet
  for 2b is a direct comparison, not a memory of a number.
- `bench/baselines.json` keyed by target (`m4-softfp`, `m4f`, `m33`, `m55`, `m55-ooura`), seeded
  from the vendored C on the Part 10 legs and committed with the SHA and toolchain versions.
- `scripts/icount.py` and the QEMU TCG counting plugin, copied from MuTap (same license, same
  plugin API pin), with `mps2-an386` added as a target.
- `bench/bench_fft.cpp`: host wall-clock microbenchmark, min-of-N, port vs C while the C exists
  and vs the recorded number afterwards. Informational on GitHub runners (noisy), the tool for
  the Apple vDSP and desktop claims locally.
- Gate: on every push and pull request, each target's scenarios within ±3% of baseline, the
  family's tolerance. `--update` is a deliberate commit that records the measured delta and the
  reason in `bench/README.md`; a regression that is "expected" is written down, never absorbed.
- The `.text` ceilings from Part 10 item 4 live in the same job.

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
and Apple arm64. The Stage 1b bench builds port-vs-C pairs on every leg and reports them; nothing is
ratcheted yet because nothing is routed. MuTap: nothing moves; no bump needed.

**2b. Flip routing.** One line in `fft.h`; `test_fft_backend.cpp`'s `ooura_ref` re-pointed at the
ported float engine (on Linux/Windows the backend test becomes port-vs-port, which is fine
only because 2a's parity test exists); README rewritten in the same PR. MuTap bump: fingerprint
harness bit-identical (double rows *and*, because 2a kept float table semantics, the float
rows), icount ratchet at 0% delta on m33 and hexagon (the ratchet's ±3% is not the gate here;
0% is, because nothing numeric changed), `test_float32`, `test_g168`, `test_nn_suppressor`
unchanged. DspTap's own ratchet (Stage 1b): the port within ±3% of the C's instruction count
on every QEMU leg and within the `.text` ceilings; the measured deltas go into
`bench/README.md` with the flip. Rollback is a one-line revert. *Amended after the bump: the
"0% delta" prediction did not hold — the counts fell on m33 and hexagon and were re-recorded
per D11; see Part 13, "Stage 2b's MuTap gate".*

**2c. Remove the C.** Delete `fftsg_float.c`; move `fftsg.c` and `readme.txt` to
`tests/reference/ooura/`; `tap_dsp_fft` exists only when `TAP_DSP_FFT_CMSIS` is on. MuTap:
rewrite the CI job that compiles `submodules/dsptap/third_party/ooura/fftsg*.c` by path
(`ci.yml:319-324`) header-only; rewrite `THIRD_PARTY_NOTICES.md`. Licensing text per Part 5.

### Stage 3 — Fixed point: Q15 and Q31 profiles of the same contract
Design in Part 7; this is the sequence. Three PRs.

**3a. Substrate.** `sample_traits<double>` (`coeff = double`, `accum = double`; D1) so the
concept covers all four sample types, plus a sibling trait for butterfly arithmetic
(`fft_arith<Sample>`: `mul_coeff` with the single documented rounding, saturating `add`/`sub`,
`shr_round`, `headroom_bits`), because the FIR trait's `mac`/`finalize` shape is the wrong one
for an FFT. Named fraction-bit constants replace the bare literals. `test_sample_traits.cpp:108`
flips; `test_decimate.cpp:192` becomes Q15-vs-double; the numpy pin does not move.

**3b. The fixed-point kernel** (`fft/fixed_point.h`): one int32 complex kernel with Q1.30
twiddles, two scaling policies (fixed, block floating point), Ooura's real post-pass formulas
so the packing and the `exp(+i)` convention are identical, Q15 and Q31 I/O widths.
`basic_real_fft<std::int16_t>` / `<std::int32_t>` route to it; aliases `real_fft_q15`,
`real_fft_q31`. Gates: the Part 9 fixed-point battery on the hosts, and the same battery on the
M4 (soft-float), M4F, M33 and M55 QEMU legs; new ratchet scenarios `rfft_q15_512`,
`rfft_q31_512`, `rfft_q31_2048` seeded in the same PR, so the `SMMULR` and Helium seams that
come later have a number to beat.

**3c. Consumers and instruments.** The analysis instruments take `std::span<const Sample>`
so they can score Q15/Q31 output; the capi and `dsptap_py` expose all four profiles; the
executed notebook measures the noise floors (Part 8). MuTap: nothing moves (no consumer is on
fixed point yet); the fixed-point projects start from the header's numbers.

### Stage 4 — Engine as an explicit parameter, with an ABI tag
`basic_real_fft<Sample, Engine = default_real_fft_engine_t<Sample>>`. The default alias lives
in one place, selected by the one remaining build define. To close F4 where it actually lives
(the consumers that embed the FFT by value), the selection also opens an inline namespace ABI
tag on `tap::dsp` (`inline namespace fft_ooura {}` / `fft_vdsp {}` / `fft_cmsis {}`), so two
images built with different defaults cannot coalesce each other's symbols. Backends move to
`fft/backends/accelerate.h` and `fft/backends/cmsis.h`, included by whoever selects them; the
CMSIS object library stays PIC. Consumer-facing transforms stay **non-const**; shareability is
an engine trait (`Engine::is_shareable`), true for the ported engine, false for the two
scratch-carrying backends. Each engine states its supported size range and construction
checks it (the CMSIS engine supports 32 … 4096 only; today `fft.h` ignores the init status,
Part 13). Typed tests over the engines available on the host: Ooura vs vDSP
same-binary on macOS; CMSIS compile-only on M55 (it cannot run on a host, and nobody runs
CMSIS-vs-Ooura parity anywhere today; that gap is recorded, not closed, by this stage). The
`m55` and `m55-ooura` baseline keys are what make a backend regression visible.

**Landed** (wave 4, tap/DspTap#35; design record in `docs/fft-design.md`, "Stage 4").
Deviations from the text above, each recorded there:
- **Tag names.** `fft_split_radix` / `fft_cmsis` / `fft_vdsp`, not `fft_ooura` / … : the plan's
  name predates 2c, when the default became the split-radix engine. Opened as
  `namespace tap::dsp::inline TAP_DSP_FFT_ABI` in `fft.h`, `pvoc.h` and `log_mel.h`.
- **Second parameter, two meanings.** `basic_real_fft<Sample, Policy>`: the engine for the
  floating profiles (default `default_real_fft_engine_t<Sample>`, one selection point in
  `fft.h`), the Scaling policy for the fixed-point ones, so no spelling in MuTap or the capi
  changes; `basic_real_fft<float | double, scaling::fixed>` (the pre-Stage-4 spelling; no
  in-tree code writes it since the #35 fix pass) resolves to the selected engine and is a
  distinct type from the one-argument form, tolerated for one consumer cycle (D4).
- **The fixed-point profiles carry the tag** (a partial specialization lives beside its
  primary) although their layout does not depend on the selection; cost stated in the design
  note (a harmless coalescing prevented; a link error instead of a silent merge).
- **`is_shareable` is `k_is_shareable`** (house `k_` prefix), and shareable implies const: a
  shareable engine's transforms are `const` (Q31's became `const` behind `requires`-constrained
  overloads), `static_assert`ed by the class, the converse a convention the tests pin; the
  class's transforms stay non-const (N11).
- **CMSIS parity is not compile-only on the M55: it runs there**, under QEMU, as it has since
  that leg landed (main's backend test already compared the CMSIS default to the split-radix
  reference); what #35 adds is the two engines as typed rows in one binary
  (`fft_backend_parity/cmsis` beside `/split_radix`) and the named split-radix routing row on
  that leg. Nowhere on a host and nowhere on hardware; recorded, not closed.
- **Construction "checks" the range as the house precondition** (`TAP_EXPECTS`: debug
  assertion, nothing in release, STYLE.md §4), with `supports_size(n)` as the constexpr,
  release-mode predicate pinned by `static_assert`s on every leg; the CMSIS init status is
  checked the same way and the `uint16_t` narrowing is guarded. Not a release-mode check:
  out of range in release remains a precondition violation (decided, 35a/F1; no defined
  fallback in the transforms), and `supports_size` is the mandatory gate wherever N comes
  from configuration — the capi applies it per profile since the fix pass, MuTap's config
  path is on the bump checklist in `docs/fft-design.md` (with the real embedder list:
  `partitioned_fdaf`, `partitioned_fdkf`, `pem_afc` + `Core`, `residual_suppressor`,
  `nn_suppressor`; `aec_chain` inherits the tag through its template arguments; no forward
  declarations of the tagged classes anywhere).
- **No `std::span` overloads** (Part 9 listed them as a Stage 4 assertion): none were added,
  so none are asserted.
- **The Apple same-binary *microbenchmark*** still needs a Mac; the same-binary parity *test*
  exists (the macOS leg's typed backend suite).
- **`TAP_EXPECTS` and the two duplicated `fft.h` paragraphs** (Stage 6 items that live inside
  `fft.h`) landed here; nothing else of Stage 6.

### Stage 5 — The packed-spectrum view, DspTap only
`fft/spectrum.h`: a non-owning view whose docstring carries the numeric definition
(`bin[k] = a[2k] + i a[2k+1]`, DC at `a[0]`, Nyquist at `a[1]`, `W = exp(+2πi/N)`, inverse
unnormalized). Primary accessors are **native**: `dc()`, `nyquist()`, `re(k)`, `im(k)`,
`power(k)`, `num_bins()`. A convention-flipping accessor, if kept, is named unmistakably
(`bin_engineering(k)`) and is not the default. Migrate `pvoc.h`, `log_mel.h`, `test_pvoc.cpp`,
gated by the existing pinned tests. MuTap adopts the view per header when each is next
touched, each such PR gated by the fingerprint harness and 0% icount, because those 77 sites
do native-convention complex products by hand in ratcheted hot loops and a wholesale rewrite
is not the mechanical migration revision 1 called it. This stage is **not** a prerequisite for
Stage 2 or 3: the packing is a kept contract, and the fixed-point profiles present it too.

### Stage 6 — Hygiene (after the port, not before)
`detail/math.h` (pi via `std::numbers`, periodic Hann, dB helpers), `tests/support/` for the
copy-pasted helpers, `solve_dense` noexcept, kaiser includes and doc drift, NOTICE names,
the two duplicated `fft.h` paragraphs, the unbraced SMLALD loops. Gate the Hann/pi
consolidation with the fingerprint harness: association order must be preserved, and
log_mel's numpy pin only bites above 1e-6. Do **not** delete `TAP_DSP_CHANNEL_PARALLEL` /
`TAP_DSP_CP_MIN_CHANNELS` until SampleRateTap and RatioTap (not on disk) have been grepped;
README says they are consumed there. `TAP_EXPECTS` lands here for `fft.h`'s power-of-two
precondition only; a repo-wide precondition policy is its own plan. (Landed at Stage 4
instead, #35, together with the size-range precondition and the two duplicated `fft.h`
paragraphs: `include/tap/dsp/detail/expects.h`, a debug assertion per STYLE.md §4.) The capi's pre-existing
defects from Part 2 (`void*` handles, exceptions crossing `extern "C"` in the old entry points,
allocation in the decimator's process path) are owned here too.

**Landed (wave 4, tap/DspTap#34), with these deviations.**
- *Delegated to Stage 4* (same wave, which owns `fft.h` and `fft/`): the two duplicated
  `fft.h` paragraphs (F9) and `TAP_EXPECTS` for the power-of-two precondition. Nothing
  under `fft.h` / `fft/` was touched here; `fft/tables.h` already takes pi from
  `std::numbers`. Also handed to Stage 4: fft.h's D9 paragraph must state the TU-context
  dependence of contraction after inlining under g++ `-march` with FMA (below and Part 13):
  the engine's output then depends on the including TU's inlining decisions, for float, and
  the double engine's codegen moved too, so "double does not move" is an observation on the
  inputs measured, not a guarantee.
- *`detail/math.h`* holds `k_pi`, `periodic_hann` and `power_db` / `amplitude_db`; pvoc,
  psola, log_mel and the two analysis instruments call it. The gate is the fingerprint tool
  (Part 13), which gained psola lines for this: 16/16 identical with linux g++ 13 at its
  default x86-64 flags (-O3, and -O0), g++ -O2 `-march=x86-64-v3`, g++ -O3
  `-march=x86-64-v3 -ffp-contract=off`, clang++ 18 -O3 with and without `-march=x86-64-v3`
  and with `-ffp-contract=fast`, and the QEMU builds M33 and M55 (MinSizeRel, the CI build
  type) and M33 and M4F (Release, VFMA); MSVC and AppleClang were not A/B'd. The one A/B that
  moves is g++ -O3 `-march=x86-64-v3` at its default `-ffp-contract=fast`, where the
  log_mel.h change moves pvoc's *float* lines through the tool's shared single-TU FFT
  instantiation (contraction after inlining; pvoc.h's own change alone is identical there).
  The move is TU-context dependent: appending ~30 lines of raw-FFT hashing to the tool's TU
  makes main and head identical again, and disassembly locates the difference in
  `split_radix_rdft::cftrec4` / `cftleaf` inlining, not in pvoc. That is the configuration
  D9 leaves unclaimed, in a form D9 does not yet describe (handed to Stage 4 above); the
  tool's header now says how to read it.
- *`tests/support/`* keeps two synthesizers, not one: the bin-exact `tone` (FFT oracles) and
  the Hz-at-a-rate `sine` (audio batteries) associate differently, and merging them would move
  the pinned numbers. Likewise `mt19937_signal` stays beside `random_signal` for the three
  batteries measured on it. Every migrated sequence was shown byte-identical to the copy it
  replaced.
- *`solve_dense`* became allocation-free (pivot rows exchanged in place) rather than losing
  `noexcept`; bit-identical to the permutation form.
- *`frontend_vectors.h`* is not reformatted (`git diff -w` would not be empty) and gets no
  directory-local `.clang-format` (it would disable formatting for `ooura_rdft.h` in the same
  directory); the generator brackets both generated headers in `// clang-format off/on`.
- *The tidy gap for `bench/`* (and `tools/capi`, `tools/fingerprint`) is closed inside
  `tests/CMakeLists.txt` with excluded-from-ALL object targets, so no taphouse follow-up is
  needed.
- *The capi contract change* (Part 2's capi defects). Handles are typed (`typedef struct
  dsptap_fft_s* dsptap_fft;` ...) instead of `void*`: binary-compatible (one pointer each; the
  ctypes bridge's `c_void_p` is unchanged), but a source break for C++ code that stores a
  handle as `void*` and for C or C++ function pointers typed on `void*` handles. New NULL
  returns where the headers' preconditions were previously unchecked: `dsptap_pvoc_create`
  above 2^28, `dsptap_psola_create` at 2^26 or more, `dsptap_yin_create` when
  window + tau_max exceeds INT_MAX (yin.h's `@pre` gained that bound; the header still only
  asserts its other preconditions, TAP_EXPECTS being Stage 4's). Every entry point is
  `noexcept`; `dsptap_decimator_process` no longer allocates (staging fixed at create);
  `DSPTAP_API` is dllexport only under `DSPTAP_BUILDING` and dllimport for consumers.
- *Not in Stage 6's list, and still present*: Part 2's three Gaussian-elimination solvers
  (`kaiser.h`'s `solve_dense`, the 3x3 in `sine_analysis.h`, and the tone tracker duplicated
  between `sine_analysis.h` and `multitone_analysis.h`).
- *`frontend_vectors_tuned.h`* was not regenerated: the generator's numpy output drifts in
  the last digit of 1,422 lines in this environment (the effect its docstring describes), so
  the committed file is main's bytes plus the six marker and comment lines, inserted by hand;
  `frontend_vectors.h` is the regenerated output, byte for byte.
- `TAP_DSP_CHANNEL_PARALLEL` / `TAP_DSP_CP_MIN_CHANNELS` are kept: SampleRateTap and RatioTap
  are still not on disk to grep.

### Gates, in one table

| Stage | DspTap gate | MuTap gate | Rollback |
|---|---|---|---|
| 0 | new overflow test | pin bump, suite green | revert 2 files |
| 1 | M4/M33/M55 legs run the emulation-sized selection | fallback leg builds Ooura float; filter pruned | n/a |
| 1b | baselines seeded from the C on all QEMU legs; ratchet job green | none | delete `bench/` |
| 2a | bit identity double+float, 3 hosts + 4 QEMU legs, `-ffp-contract=off`; oracle; bench artifact | none (no bump) | delete header |
| 2b | existing battery on the port; backend test re-pointed; ratchet ±3% vs C; `.text` ceilings | fingerprint identical; icount 0%; float pins unchanged | one-line revert |
| 2c | build without the C | CI job rewritten; notices | re-pin |
| 3a | traits battery incl. double; decimate re-pinned Q15-vs-double | pin bump, suite green | re-pin |
| 3b | Part 9 fixed-point battery on hosts and all four QEMU legs; fixed-point scenarios seeded; `.text` ceilings | none | delete header |
| 3c | instruments typed; capi/notebook executed | none | per item |
| 4 | typed engine tests; macOS same-binary parity — **landed** (#35): typed over the engines the leg builds, CMSIS beside split-radix on the M55, vDSP beside it on macOS; icount +0.00 % on every key without a re-record | fingerprint identical (on the bump; MuTap's own embedders adopt the tag in `tap::mu`) | re-pin |
| 5 | pinned pvoc/log_mel tests | per-header fingerprint + icount 0% | per header |
| 6 | fingerprint on Hann/pi change | pin bump | per item |

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
to make *double* outputs identical across hosts; it is a Stage 3a candidate because the
fixed-point twiddles are generated there anyway, and is otherwise out of scope.

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
leaves is what makes Ooura fast). Stage 1b exists so this is measured, not argued: the
instruction-count ratchet on every QEMU leg with the C as the seeded baseline, the `.text`
ceilings, and the host microbenchmark. Load `m_w.data()` into a local once per transform so
aliasing analysis matches the C, which receives `w` by parameter. Thresholds are the family's
±3%, decided here, not after the fact.

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
  fit Q1.14, so "the same table" is false for the first stage. Resolution: the fixed-point
  design is rebuilt in Part 7 with these constraints as inputs (one int32 kernel, Q1.30
  twiddles, shift-before-multiply, BFP for round-trip consumers, its own inverse scaling).

### Should-fix (all adopted)
- **P6/N15/A — Stage 2's coupling to D1 was artificial**; the engine is `std::floating_point`
  until fixed point. D1 also had the coefficient type wrong (`coeff` must be `double` or the
  double FFT cannot be bit-identical). Moved to Stage 3a, corrected.
- **P7/A — The counter bug shipped inside a hygiene PR.** Now Stage 0, alone.
- **P8 — Deleting `TAP_DSP_CHANNEL_PARALLEL` without grepping SampleRateTap/RatioTap.** Held.
- **P9/N12/A — D4 did not close F4 for the embedding consumers.** ABI tag added; F4 reworded to
  cross-image.
- **P10 — F2 severity overstated for MODULE consumers.** Reworded.
- **P11 — "True INTERFACE target" was conditional and contradicted by the escape hatch; no
  `install()` exists.** Reworded; hatch dropped.
- **P12 — Fixed point is scope creep against every consumer's written plan** (MuTap
  `wake-word-plan.md:978-981`). Overruled by the maintainer: fixed-point projects not yet
  started need FFT support. Kept on the main line as Stage 3; MuTap's plan is unaffected
  because MuTap does not consume it.
- **P14 — Stage 3 was not rollback-able.** Split into 2a/2b/2c.
- **P15/N17/A — No performance, size or compile-time gate.** Added (Part 4).
- **P18 — A MuTap CI job compiles the vendored C by path.** Stage 2c.
- **N9 — `extern template` does nothing for in-class-defined members** (verified on both
  compilers). Dropped, with the measured compile cost that makes it unnecessary.
- **N10 — Header-only is size-neutral** because `--gc-sections` is already on. Stated, with
  measured sizes.
- **N11 — Const transforms would make the public API's constness depend on the selected
  engine.** Public methods stay non-const; `Engine::is_shareable`.
- **N13 — Runtime twiddle derivation in `cftf1st`.** The fixed-point kernel does not reuse
  Ooura's first stage; it has a full Q1.30 table (Part 7).
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

## Part 7 — Fixed-point design (Stage 3)

Inputs: the arithmetic constraints from Part 6 (N6-N8, N13), the substrate's rule that fixed
point is traits over raw sample types with the Q ladder visible at the use site, and the
literature: Welch 1969, "A fixed-point fast Fourier transform error analysis" (IEEE Trans.
Audio Electroacoust.); Oppenheim & Weinstein 1972, "Effects of finite register length in
digital filtering and the fast Fourier transform" (Proc. IEEE); block floating point per
Oppenheim & Schafer, *Discrete-Time Signal Processing*, the finite-precision FFT section. No
shipping product's behaviour is reverse-engineered; CMSIS-DSP's *documented* output-scaling
convention is noted only where compatibility with it is a stated goal.

### One kernel, two I/O widths
- **Internal data is int32 for both profiles.** Q31 for the Q31 profile; the Q15 profile widens
  on input and narrows (round-half-up, saturating, the substrate's `finalize` rule) on output.
  This is what makes the Q15 profile usable: fixed 1/N scaling in 16-bit data gives about
  25 dB per-bin SNR at −40 dBFS for N=512 (Part 6, N8); in 32-bit data the same scaling leaves
  22 bits at N=512, a floor below −130 dB. The cost is an int32 work buffer of N allocated at
  construction; the in-place API is preserved at the caller's int16 buffer.
- **Q15 input placement uses the spare width as guard bits.** Input is shifted left by 14, not
  16, so full-scale ±1 has two guard bits and the worst-case radix-4 growth (component bound
  4√2 per stage under rotation) cannot saturate under fixed scaling without any input-level
  contract. The Q31 profile has no spare width: under fixed scaling it takes one input
  pre-shift (−6 dB, one bit of 31) and the header says so; under block floating point no
  pre-shift is needed because scaling happens only when the data has grown.
- **Twiddles are Q1.30 int32** (`sample_traits<std::int32_t>::coeff`) for both profiles, a full
  table for the complex kernel and the post-pass, generated in double at construction and
  rounded once (`round_sat`). 1.0 is representable, so no special case. The Q1.14 idea from
  revision 1 is dropped: it was the wrong width for an int32 kernel and could not hold Ooura's
  run-time first-stage twiddle sums.
- **Multiply is shift-before-multiply with one rounding.** `int32 × Q1.30 → int64`, then
  `>> 30` with round-half-up: the single documented rounding point per product. Portable C++
  first. The Armv7E-M/Armv8-M `SMMULR`/`SMMLAR` (32×32 high-half with rounding) round at
  bit 32, not bit 30, so they can NOT reproduce this operation bit-exactly (verified: ~12% of
  random products differ; Part 13). Unlike `SMLALD` for the FIR kernels, an Arm fast path here
  is a separately pinned profile or a second named operation with the hardware's rounding
  point that the kernel selects explicitly.
- **Structure is its own radix-4 DIF complex kernel of length N/2 with a radix-2 final stage
  for odd log2, plus Ooura's real post-pass formulas.** Part 6 established that Ooura's
  split-radix graph can be reused only as a data-flow graph with a different operation order,
  and that its first stage derives twiddles at run time; a dedicated kernel is simpler to
  scale, simpler to prove against Welch's model, and shares the *contract* (packing, `exp(+i)`,
  DC/Nyquist slots) rather than the code. Bit-reversal table and twiddle generator live in
  `fft/tables.h`, shared with nothing else by design.

### Two scaling policies, as a template parameter
- **`scaling::fixed`.** Shift-before-butterfly by 2 bits per radix-4 stage (1 per radix-2) and
  one more in the real post-pass, so the forward output is exactly `X / N` in the sample's Q
  format, deterministic, with a statically known exponent. For magnitude and power consumers
  (log_mel-class front ends, detectors).
- **`scaling::block_floating`.** Before each stage, a `clz` scan of the block finds the
  available headroom; the stage shifts only as much as growth requires and accumulates the
  exponent, which the transform returns. Sub-blocks of a DIF kernel finish together per stage,
  so the bookkeeping is one exponent per transform. For round-trip consumers (overlap-add
  resynthesis, adaptive filters), where fixed scaling's forward-and-inverse halving would
  discard log2(N) bits twice.
- **The inverse has its own scaling.** Ooura's unnormalized inverse has structural gain N/2 and
  saturates on the first stage if run unscaled in fixed point (Part 6, N7). Under `fixed` the
  inverse halves per stage the same way and the round-trip factor is a fixed power of two,
  stated in the header and pinned by a test; under `block_floating` the round trip returns
  `x · 2^-(e_fwd + e_inv)` with both exponents returned. Whether to match CMSIS's documented
  q15/q31 output scaling is a compatibility decision made in 3b, not an accident.

### Numbers the header states and the tests pin (Part 9)
Per profile and policy: forward scale; inverse scale; round-trip identity; saturation-free
input level; twiddle quantization (max |w_q − w| ≤ 0.5 LSB of Q1.30); per-bin noise floor vs
input level at N = 256/512/2048 for white noise and an on-bin tone, against Welch's model and
against the double golden model; DC of a constant exact; worst-case growth pattern does not
saturate; NaN not applicable; no denormals; latency 0; `noexcept`, allocation-free; copyable;
`is_shareable` true.

### Target notes
Cortex-M4 (Armv7E-M) and Cortex-M33 (Armv8-M Mainline) both define `__ARM_FEATURE_DSP`, so
`SMMULR`/`SMLALD`-class instructions exist on both; the M4 without FPU is the profile's reason
to exist. The M55 has Helium (MVE) with `vqdmulh`-class Q31 lanes; that is a later backend
behind the same contract, as CMSIS f32 is today. Measured `.text` and instruction counts per
leg are CI artifacts from 3b onward (Part 10).

---

## Part 8 — Documentation plan

The house rule is that the header owns the contract as numbers and the tests pin them. The
FFT work adds three documentation surfaces and touches three more.

1. **Header docstrings are the contract.** `fft.h` carries, per profile: packing, sign
   convention, forward/inverse scale, round-trip identity, saturation-free input level,
   noise-floor numbers, latency, alignment (none), NaN/denormal behaviour, fp-contraction
   policy, thread-shareability, real-time safety, and the measured performance table with the
   target and date. Each contract bullet ends with the name of the test that pins it, the
   pattern README already uses. `@param` per STYLE.md, not `\param`. The engine header
   (`fft/split_radix.h`) carries the transliteration rules from Part 4 (statement fidelity,
   explicit `static_cast<double>` on libm calls) so the next person does not undo them.
2. **`docs/fft-design.md`** is the design note that is too long for a header, in the role
   `kaiser.h`'s design note plays for the FIR substrate: why split-radix for floating point and
   a radix-4 kernel for fixed point, the halving-count and headroom derivations, the
   twiddle-quantization argument, the Welch noise model against the measured floors, the
   contraction policy and its measurements, the per-target size and instruction-count tables,
   the licensing and provenance statement, and a changelog of contract-affecting SHAs. This
   audit document is the plan of record until every stage has landed; each PR links its stage.
3. **README FFT section rewrite** at Stage 2b, extended at 3b: the profiles table
   (double / float / Q15 / Q31 with scale, floor, target), the backend table, the two-line
   migration note for consumers (what words change, what does not), and the "header-only
   unless CMSIS is on" statement. The intro's "seven primitives" line stays; profiles are not
   primitives.
4. **Provenance and licensing** per Part 5: the port header's banner, `NOTICE.md`, MuTap's
   `THIRD_PARTY_NOTICES.md`, and one glossary line in MuTap's `docs/itu-compliance.md`.
5. **CLAUDE.md** updates at 2c and 3b: the header-only statement, the four-profile ladder in
   the design-discipline section, the embedded CI legs in the build section, and one line in
   the adding-a-primitive checklist: "state each profile's numbers, and name the on-target leg
   that runs its battery".
6. **Verification layer.** `tools/capi` and `notebooks/dsptap_py.py` gain FFT entry points for
   all four profiles (today the notebooks cannot measure the FFT at all), and
   `notebooks/fft.ipynb` is committed *executed*: float-vs-double, the per-bin noise floor vs
   input level for Q15 and Q31 under both scaling policies (the Welch curve, measured on the
   shipping C++), BFP exponent behaviour on speech-like material, and the backend comparison on
   the host that runs it.
7. **Consumer-side docs** ride the pin bumps: MuTap-Max's per-external CMake comments about the
   compiled FFT library go when a bump proves they can; MuTap's `docs/optimization.md` M55
   section is corrected with the CI flag fix (Stage 1).

---

## Part 9 — Unit testing plan

Structure first, then what each file pins. Everything is GoogleTest, typed where the house
pattern is typed, and named for the promise it pins.

### Files
- `tests/test_fft.cpp` (existing, widened): the contract battery, `TYPED_TEST_SUITE` over
  `float, double, std::int16_t, std::int32_t`, with a per-type traits helper supplying the
  tolerance, the forward scale and the input conversion so the same `RoundTripReproducesInput`,
  `ImpulseHasFlatSpectrum`, `DcAndNyquistPacking`, `SignConventionIsPlusI`,
  `ParsevalEnergyConservation` run on every profile. Cross-precision: `FloatTracksDouble`
  (existing), `Q15TracksDouble`, `Q31TracksDouble`, each a measured number pinned at 2×
  margin, the log_mel pattern.
- `tests/test_fft_parity_ooura.cpp` (Stage 2a): its **own CMake target** so it can carry
  `-ffp-contract=off` for both the C and the C++; bit identity for `double` and `float`,
  forward and inverse, N = 4 … 65536 plus one run at 2^20, on broadband, on-bin tone, impulse,
  DC, and a full-scale alternating pattern. A second target at default flags reports the
  measured max-ulp deviation as an informational number, never as a pass/fail with an
  assumed bound. Runs on the three hosts and, size-limited to N ≤ 4096, on the four QEMU legs
  (same-binary identity against newlib's libm is exactly the property that matters there).
- `tests/test_fft_oracle.cpp`: the independent oracle. Closed-form vectors (impulse, DC,
  Nyquist, on-bin cosine and sine, two-tone) with exact answers, and a compensated-summation
  DFT (double-double, not `long double`, which is `double` on MSVC and Apple arm64) for
  N ≤ 256, driven against all four profiles with the profile's documented scale applied.
- `tests/test_fft_fixed.cpp` (Stage 3b): saturation-free worst case (the packed-pair pattern
  at a 45° twiddle from Part 6, plus full-scale ±1 and `INT16_MIN` inputs) does not wrap;
  forward scale exactly `X/N` under `fixed` on an integer-valued impulse and DC; round-trip
  identity per policy as the header states it; BFP exponent reconstructs the fixed-scaling
  result bit-for-bit after shifting; per-bin noise floor vs level at −0/−20/−40/−60 dBFS for
  N = 256/512/2048 pinned against the Welch-model number and the double golden model;
  twiddle quantization ≤ 0.5 LSB; rounding symmetry on negated input (the half-up bias is
  bounded and stated); Q15 and Q31 agree with each other to the Q15 floor.
- `tests/test_fft_backend.cpp` (existing, re-pointed at 2b, typed over engines at Stage 4):
  parity, alignment stability and tonal accuracy as today; `ooura_ref` becomes the ported
  float engine; on macOS the Ooura-vs-vDSP comparison is same-binary.
- `tests/test_fft_rt.cpp`: `static_assert(noexcept(...))` on every transform (the nn pattern);
  an operator-new counting guard proving no allocation inside the transforms for every
  profile and engine; copy and copy-assignment leave the source and copy bit-identical;
  `is_shareable` is what the header says; `std::span` overloads equal the pointer overloads.
- `tests/test_spectrum.cpp` (Stage 5): the view's numeric definition against forward outputs,
  including the DC/Nyquist slots and the sign of `im(k)` for a sine.
- `tests/support/` (Stage 6, but started at 2a for the new files): one `random_signal`, one
  `xorshift32`, one tone synthesizer, one dB helper; the three existing copies migrate later.
- `tests/bare_metal_main.cpp` + `TAP_DSP_BARE_METAL` one-shot mode (Part 10): a positive filter
  of the emulation-sized selection, judged on the gtest summary text.

### Rules
- Every contract bullet in a header names its test; every test's name is the promise.
- Fixed seeds, no wall-clock, no filesystem (bare metal has none).
- Tolerances are measured numbers with a stated margin, never round numbers.
- A cross-precision test compares to the **double** golden model, not to a sibling profile,
  except the one Q15-vs-Q31 agreement test whose purpose is the I/O-width narrowing.
- The benchmark target is informational and lives under `bench/`, not `tests/`; the
  on-target `.text` and instruction-count assertions live in CI, not in gtest.
- MuTap's fingerprint harness (Stage 1) is the downstream bit-identity gate for every bump.

---

## Part 10 — Embedded CI: Cortex-M4, M33 and M55 under QEMU

### What exists to copy
MuTap already runs its float suite on QEMU: `cmake/arm-cortex-m33-mps2.cmake`
(`-mcpu=cortex-m33 -mthumb -mfloat-abi=hard`, `--specs=rdimon.specs -nostartfiles`,
`qemu-system-arm -M mps2-an505 -semihosting -kernel`), `platform/armv8m_startup.c` (161 lines:
vector table, `MSPLIM`, CPACR enable for CP10/CP11, bss clear, semihosting init, static
constructors, `exit(main())`), `platform/mps2_an505.ld` (CODE at 0x10000000, DATA at
0x38000000, 4 MB each), `platform/mps3_an547.ld` for the M55, and a one-shot gtest mode in
`tests/CMakeLists.txt` (`gtest_disable_pthreads`, `GTEST_HAS_POSIX_RE=0`,
`GTEST_HAS_STREAM_REDIRECTION=0`, `GTEST_HAS_FILE_SYSTEM=0`, `PASS_REGULAR_EXPRESSION` on a
completion marker because semihosting does not propagate exit codes reliably). Its M33
selection runs in about 3 minutes of emulation. DspTap's M55 toolchain is a copy of MuTap's
with the platform files missing.

### What adding the legs entails
1. **`platform/` in DspTap**: the startup file with the `MSPLIM` write guarded by
   `__ARM_ARCH_8M_MAIN__` (it does not exist on Armv7E-M; the CPACR write is harmless on a core
   without an FPU), `mps2_an505.ld`, `mps3_an547.ld`, and a new **`mps2_an386.ld`** for QEMU's
   Cortex-M4 board model (`-M mps2-an386`; code SRAM at 0x00000000 and data SRAM at 0x20000000,
   sizes to be taken from QEMU's `hw/arm/mps2.c` when the script is written, not from memory).
2. **Toolchain files**: `arm-cortex-m4-mps2.cmake` in two flavours selected by a cache
   variable, `-mcpu=cortex-m4 -mfloat-abi=soft` (no FPU: the fixed-point profile's reason to
   exist; everything float is soft) and `-mcpu=cortex-m4 -mfpu=fpv4-sp-d16 -mfloat-abi=hard`
   (M4F); `arm-cortex-m33-mps2.cmake` as MuTap's; the existing M55 file with its link line
   made valid. All set `TAP_DSP_BARE_METAL ON` (the variable is renamed from `MUTAP_`).
3. **Test harness**: the one-shot mode and `tests/bare_metal_main.cpp` with a positive filter:
   the fixed-point battery in full; the float contract battery; the parity test at N ≤ 4096;
   decimate, log_mel and nn float suites; the double suites limited to the small-N contract
   tests (soft-float double is a correctness check here, not a profile). Budget: under
   5 minutes of emulation per leg, measured and written into the job comment as MuTap does.
4. **CI matrix**: four jobs, `cortex-m4-softfp`, `cortex-m4f`, `cortex-m33`, `cortex-m55`
   (the last replacing today's compile-only job; CMSIS stays on there and its parity test
   finally runs somewhere). Each installs `gcc-arm-none-eabi` and `qemu-system-arm` from apt
   (Ubuntu 24.04 ships QEMU 8.2, which has `mps2-an386`, `mps2-an505` and `mps3-an547`),
   configures `MinSizeRel`, builds, runs `ctest`, then runs `arm-none-eabi-size` on the test
   binary and on a size-probe object that instantiates each profile alone, asserting per-leg
   `.text` ceilings for the float and Q31 real FFT at N = 512.
5. **Instruction counts** are Stage 1b, not optional: MuTap's `scripts/icount.py` and TCG plugin
   pattern with `bench/baselines.json`, scenarios per Part 11, ratcheted at ±3% per leg. This
   is the gate that catches a port that vectorizes worse than the C, and later a fixed-point
   seam that does not pay for itself.
6. **Side effects worth knowing**: both M4 and M33 define `__ARM_FEATURE_DSP`, so
   `fir_kernels.h`'s `SMLALD` path is compiled and *run* for the first time in this repo (the
   unbraced-loop style defect in it becomes visible to a compiler); M4 soft-float is the first
   leg where `double` costs enough that any accidental double in a float or fixed path shows
   up in the emulation budget; `-Wconversion` on `int16x2_t` casts needs checking under the
   M4's GCC.
7. **Cost**: roughly 10–15 minutes wall per leg, dominated by the toolchain install and the
   gtest build; emulation is minutes. Run on `push` and `pull_request` like the rest of CI;
   the icount jobs, if adopted, only on `pull_request`.

---

## Part 11 — Benchmarks and the performance ratchet (Stage 1b)

### Why a ratchet and not a benchmark
Wall-clock numbers on shared CI runners are noise; instruction counts under QEMU's TCG plugin
are deterministic to the instruction, which is why MuTap, SampleRateTap and RatioTap all
gate on them. DspTap is the primitive layer under those ratchets, so a regression here shows
up downstream as a mysterious chain-level delta someone "updates" past. The primitive needs
its own gate, per profile, per target, seeded before anything changes.

### Scenarios
One binary per scenario, selected by compile definitions (bare metal has no argv), all built
from `bench/icount/icount_main.cpp` the way MuTap's are:

| scenario | profile | N | targets | from |
|---|---|---|---|---|
| `rfft_f32_512`, `rfft_f32_2048` | float | 512, 2048 | all legs | Stage 1b |
| `rfft_f64_512` | double | 512 | host-class only (soft-float double is not a profile) | Stage 1b |
| `rfft_q15_512` | Q15 I/O | 512 | all legs | Stage 3b |
| `rfft_q31_512`, `rfft_q31_2048` | Q31 | 512, 2048 | all legs | Stage 3b |

Each scenario constructs one plan, then runs a forward+inverse loop over an xorshift-generated
signal for enough iterations that construction is under 1% of the count, folding every
output into a checksum that is printed (defeats dead-code elimination, pins determinism). No
`<random>`, no `std::vector` growth in the loop, no double anywhere in the float and fixed
scenarios (the M4 soft-float leg would otherwise measure libgcc).

`TAP_DSP_BENCH_ENGINE` selects the engine under test: `reference_c` (the vendored Ooura, while
it exists), `split_radix` (the port), `fixed_point`, `cmsis`, `accelerate`. From Stage 2a until
2c each float scenario builds twice, C and port, and the job prints both and their ratio; that
ratio is the Stage 2b gate.

### Baselines and targets
`bench/baselines.json` keyed by target: `m4-softfp`, `m4f`, `m33`, `m55` (CMSIS on, the
deployed profile), `m55-ooura` (the fallback), and, if a Hexagon leg is ever added here,
`hexagon`. Seeded by running `scripts/icount.py --update` in the target's CI environment and
committing, with the SHA, GCC and QEMU versions recorded in `bench/README.md` next to the
number. The plugin header is pinned to the QEMU version Ubuntu ships, digest-verified, as
MuTap does.

### Policy
- Tolerance ±3%, the family's number. Tighter is not better: link-order and libgcc changes
  move counts by fractions of a percent legitimately.
- The ratchet runs on push and pull request on every QEMU leg. A red ratchet is a failing
  check, not a warning.
- `--update` is a commit on its own, with the measured before/after per target and the reason
  in `bench/README.md`. Expected regressions (a correctness fix that costs instructions) are
  written down; there is no silent absorb.
- Size runs in the same job: `arm-none-eabi-size` on a probe object per profile, ceilings per
  leg in the workflow file, recorded the same way.
- The host microbenchmark (`bench/bench_fft.cpp`, min-of-N wall clock, port vs C or vs the
  recorded number) is informational in CI and the local tool for desktop and Apple vDSP
  claims. Its numbers go into `docs/fft-design.md` with the machine and date, never into a
  gate.

### What the ratchet catches, by stage
2b: a port that vectorizes worse than the C. 3b: nothing yet (it seeds), but every later
`SMMULR`, Helium or table-layout change to the fixed-point kernel has a number to beat. 4: a
backend regression on the deployed M55 profile, or the fallback quietly getting slower. 6: the
hygiene pass cannot slow the hot path unnoticed.

### Sharing
The plugin source and `icount.py` would then exist in four repos. Copy now (they are small and
MIT); a later taphouse-style consolidation is the right home, out of scope here.

---

## Part 12 — Parallelization: what a swarm can and cannot do here

### What serializes
- **The port is one hand.** Bit identity depends on statement-for-statement fidelity across
  ~2,580 lines; two agents transliterating different function groups produce two conventions
  in one file, and the parity test reports failure without location. One agent ports; a
  second writes the parity target, the oracle and the fp-contract flag handling *first*, so the
  porter works test-first against a red suite.
- **The critical path is linear**: Stage 1 embedded legs → 2a → 2b → 2c → 4. Each step's gate is
  a CI run; the development container has neither the Arm toolchain nor QEMU, so agents write
  those stages and only CI proves them. Every such PR needs a watcher agent driving it to green.
- **Contention files**: `fft.h`, root and tests `CMakeLists.txt`, `.github/workflows/ci.yml`,
  `README.md`. One named owner per shared file per wave; the bench workflow lives in its own
  `bench.yml` so Stage 1b does not collide with the CI-legs owner.
- **Consumer bumps are serial by the rollout rule** (DspTap squash → MuTap pin+code → MuTap-Max
  pin), and the fingerprint harness runs once per bump.

### Waves
| Wave | Parallel work items (one agent each) | Depends on |
|---|---|---|
| 1 | Stage 0 counter bug · Stage 1 MuTap (flag typo, filter, fingerprint harness) · Stage 1 DspTap (platform files, M4/M33 toolchains, one-shot harness, four-job matrix) · Stage 1b scaffold (scenarios vs the C, plugin, script, `bench.yml`) · Stage 2a tests written first against the C · Stage 3a substrate · Stage 5 spectrum view · licensing/NOTICE + `docs/fft-design.md` skeleton · capi + notebook FFT exposure for the current profiles | nothing |
| 2 | Stage 2a port (one agent, test-first against wave 1's suite) · Stage 3b fixed-point kernel (one agent) + its Part 9 battery (partner agent, test-first; may start in wave 1 on a branch stacked on 3a) · Stage 1b baselines seeded once the QEMU legs are green | wave 1 merged |
| 3 | Stage 2b flip + ratchet · Stage 2c removal · Stage 3b merge · Stage 3c instruments/capi/notebook for fixed point | wave 2 |
| 4 | Stage 4 engine parameter + ABI tag · Stage 6 hygiene · measured numbers into docs · MuTap/MuTap-Max bumps (serial) | wave 3 |

Roughly fifteen serial PRs become four waves; the port and the fixed-point kernel are the two
long poles and neither shortens with more agents.

### Rules for a swarm on this repo
- Every agent runs the hosted build and `ctest` with `-DTAP_DSP_WERROR=ON` before pushing;
  embedded and bench legs are proven by CI, and the agent that opened the PR watches it.
- Two hostile reviewers per PR before merge (numerics, and process/downstream), as was done
  for this plan; findings are verified against source before being acted on.
- No agent edits a contention file it does not own in the current wave; needed changes go to
  the owner as a request.
- Merge order inside a wave follows the table's left-to-right order to minimize rebase churn.

---

## Part 13 — Wave-1 hostile-review amendments

Eighteen reviews (two per PR) on the nine wave-1 PRs. Findings that change this plan, as
opposed to the PRs, are recorded here; the per-PR findings live on the PRs.

- **On-target test selection must be a negative filter.** #17's positive filter would have
  silently dropped every suite added by #24, #22 and #19; a guard on the selected count cannot
  see omissions. Stage 1 text amended; the fix lands in #17.
- **Passing test output is invisible in both repos' CI** (`ctest --output-on-failure`), so
  measured numbers written into PR bodies were unverifiable from logs and the informational
  ulp table from #24 recorded nothing anywhere. Pattern for both repos: run the informational
  and fingerprint targets with `-V` as their own step and exclude them from the battery step.
  GitHub's default `run` shell is `bash -e {0}` without `pipefail`; steps that pipe through
  `tee` must set it.
- **The bench checksum must be an integer fold over bit patterns.** A floating running sum
  absorbs 1-ulp output differences and cannot fingerprint the port against the C (Stage 2b's
  whole purpose). Seeding must run only on pushes to `main`, never from a pull request, and an
  emptied key must fail against the base branch's baselines.
- **`SMMULR` is not a bit-exact seam for `fft_arith::mul_coeff`** (rounding at bit 32 vs 30).
  Part 7 amended.
- **The arithmetic trait's docstring is the 3b kernel's spec and must state**: shift-before-
  butterfly (2 bits per radix-4 stage, 1 per radix-2, 1 before the real post-pass) and the
  resulting bound; the complex-multiply rounding count (one rounding per real product, two
  per complex product); that a BFP kernel never shifts by the full headroom; the Q31
  fixed-scaling input pre-shift; that `headroom_bits` of an all-zero block is 31 and of `{-1}`
  is 31. Lands in #22.
- **The packed-spectrum view must not foreclose int16/int32**: `power()` needs a promoted
  type and the conjugating accessor a floating-point constraint. Lands in #19.
- **Fingerprint tooling.** DspTap needs a committed same-host A/B fingerprint tool for pvoc
  and log_mel (golden hashes would fail across the three hosts' libm by design); MuTap's
  harness must hash the estimate output stream too. The Stage 6 gate presupposes the DspTap
  tool.
- **Licensing facts corrected** (#18): upstream `fftsg.c` carries no notice; the banner on the
  vendored copy is Tap's and the file is whitespace-normalized, so it is not "unmodified";
  WebRTC's Ooura LICENSE quotes a broader notice ("use, copy, modify and distribute … include
  commercial use") that appears in neither of Ooura's current tarballs; Ooura's current
  published address is at Kyoto University, not the 2001 Tokyo one. NOTICE, README, the design
  note and the draft email are corrected in #18; a `LICENSES/LicenseRef-Ooura.txt` carries the
  notice text so the SPDX reference resolves.
- **Branch protection on `main` is off**; none of the CI legs gate anything. Recommend
  requiring the seven CI jobs plus drift and clang-tidy once wave 1 lands, and the five
  `icount <key>` jobs after the seeding commit (never the artifact-merge job).
- **The capi's audit defects had no owning stage.** Stage 6 amended.
- **The plan branch itself needs a PR** so the design note's link resolves; opened with these
  amendments (#25).
- **CMSIS backend size range (found by #24's oracle on the M55 leg).** `fft.h`'s CMSIS wrapper
  ignores `arm_rfft_fast_init_f32`'s return status; CMSIS-DSP supports N = 32 … 4096 only, so
  the documented contract "power of two, ≥ 4" is undefined behaviour under `TAP_DSP_FFT_CMSIS`
  outside that range (the existing battery started at 64 and never saw it; N = 4 hard-faults).
  Owned by Stage 4: the engine reports its supported range, construction checks it, and the
  backend's docstring states it as a contract number. Until then the oracle carries the range.
  Done at #35 (`fft/backends/cmsis.h`: `k_min_size = 32`, `k_max_size = 4096`, the init
  status checked; the oracle reads the class's range).
- **Contraction and the single-TU fingerprint (wave 4, Stage 6).** `tools/fingerprint`
  compiles pvoc, log_mel and psola into one TU, so they share one float FFT instantiation.
  Under g++'s GNU-mode default `-ffp-contract=fast` on FMA hardware, contraction happens after
  inlining, and a change to one header can move another primitive's float lines with no change
  to that primitive's arithmetic (measured: the log_mel Hann call-site change moved pvoc float
  at `-march=x86-64-v3`; the move depends on the TU's contents and sits in the FFT engine's
  `cftrec4` / `cftleaf` inlining, float and double codegen alike). An A/B at such flags reads
  as "codegen changed"; arithmetic identity is proven at `-ffp-contract=off` or at flags that
  do not contract across statements, and the QEMU builds (GCC, VFMA) were identical as well.
  The same holds for any consumer TU that combines primitives; MuTap's harness today includes
  only `fft.h` and `nn.h`, which Stage 6 leaves untouched.
- **Wave-4 amendments (Stage 4, tap/DspTap#35).** (1) The ABI tag is named for what ships:
  `fft_split_radix` / `fft_cmsis` / `fft_vdsp` (`fft_ooura` predated 2c). (2) The engine
  argument alone already separates `basic_real_fft<float>`'s own weak symbols between two
  differently-built images; the tag is what separates the *embedders'* (`basic_pvoc`,
  `basic_log_mel`, MuTap's `partitioned_fdaf` …), which is the F4 hazard as stated; measured
  with `arm-none-eabi-nm` on one TU built for the M55 with and without `TAP_DSP_FFT_CMSIS`: 50 of
  50 weak symbol names shared before, 0 of 65 after; and in one process (`tests/test_fft_abi_tag.cpp`,
  two images from one TU, `RTLD_GLOBAL`) an untagged embedder in the second image ran the
  first image's code (layout seen 80 vs real 184) while the tagged one was immune. MuTap's
  embedders live in `tap::mu` and must open the same inline namespace on their bump (checklist
  in the design note). (3) The fixed-point profiles carry the tag because a
  partial specialization lives beside its primary; their layout does not depend on the
  selection, so the tag costs a duplicate instantiation per differently-built image and buys
  nothing they needed. (4) **The ratchet counts the harness's own fold, and the fold's register
  allocation is a function of everything inlined into `main`.** Moving the CMSIS engine's
  construction from an out-of-line `make_floating_engine<float>` to an inlined constructor
  moved the `m55` float keys by +12.00 % / +11.47 % with the transform byte-identical (checksums
  equal, CMSIS archive identical, wrapper loops identical instruction for instruction): GCC
  spilled the fold's 64-bit hash and rematerialized the FNV prime, +6 instructions per sample.
  **The harness was the defect**: the fold's codegen depended on what was inlined before it,
  and the recorded baselines had reached the CMSIS constructor through an out-of-line
  `make_floating_engine<float>`. Resolved at #35 without a re-record by making
  `bench/icount/icount_main.cpp` construct the transform under test through a non-inlined
  `make_fft()` in both `run()` bodies — the baseline's shape — so the count is the transform
  plus a fixed fold: m55 −2 / +38 instructions against the float baselines, +5 (the call) on
  every fixed-point scenario, every key inside the band. A first version had instead put a
  `noinline` attribute on the two accelerated engines' constructors; the hostile review
  measured the harness fix landing closer with nothing in the library, and the attribute
  costing 860 bytes of MinSizeRel `.text` on the m55 f32 probe; shipping code carries no
  benchmark-shaped attribute.
  (5) "CMSIS compile-only on M55" (Stage 4 text, P24) understated it: the CMSIS rows run under
  QEMU on that leg beside the split-radix rows; still nowhere on a host or on hardware.
- **Stage 2b's MuTap gate "icount ratchet at 0% delta on m33 and hexagon (0% is the gate,
  because nothing numeric changed)" was a prediction that did not hold** (wave 3; recorded at
  the Stage 2c fix pass, tap/DspTap#32). Measured on the MuTap bump to DspTap `ae0c027`
  (tap/MuTap#54; PR run 35866422618 against the previous `main` run 35804748507): every
  scenario's output checksum unchanged and the fingerprints identical on all legs, but the
  instruction counts moved downward — m33 chain −1.17/−1.15%, fdkf −2.46/−2.03%, shadow
  −3.05/−3.05% (crossing the two-sided ±3% band), suppressor −1.13/−1.14%; hexagon −0.004 to
  −0.86%; m55 (CMSIS) ≤ 0.003%. Cause: the header-only port compiles to fewer instructions
  than the vendored C on GCC 13.2 Cortex-M and on Hexagon clang, the same direction DspTap's
  own Stage 2b re-record measured (tap/DspTap#31). "Nothing numeric changed" was true of the
  outputs and said nothing about the instruction stream. Disposition: re-recorded per D11 with
  a written commit in the MuTap bump. The gate for a bump is "fingerprints byte-identical on
  every leg"; the ratchet is re-recorded when the engine changes. The Stage 2b paragraph in
  Part 3 keeps its text and carries a one-line pointer here.

---

## Decisions (revision 3)

**D1. `double` becomes a sample format. Settled.** `sample_traits<double>` with
`coeff = double`, `accum = double`, landed in Stage 3a; the concept covers four types; decimate's
Q15 oracle is re-pinned against double. Not on the port's (Stage 2) critical path.

**D2. Fixed point is its own radix-4 kernel over int32 data, sharing the contract, not
Ooura's code. Settled.** Part 7 has the reasons (operation order, run-time first-stage twiddles,
headroom). Both profiles run the same kernel; Q15 is an I/O width.

**D3. Two scaling policies as a template parameter, each with both directions' scales as
header numbers. Settled.** `fixed` for analysis consumers, `block_floating` for round-trip
consumers; the inverse has its own per-stage scaling; CMSIS compatibility of the output
convention is decided in 3b.

**D4. Engine as a template parameter with a build-selected default, plus an inline-namespace
ABI tag derived from the selection. Landed (#35).** The default alone leaves the layout-by-define
hazard in every class that embeds the FFT by value; the tag closes it — for DspTap's own
embedders at #35, for MuTap's on its bump. Tag names `fft_split_radix` / `fft_cmsis` /
`fft_vdsp`; the second template argument keeps meaning "scaling policy" on the fixed-point
profiles, which therefore carry the tag too (recorded in Part 13). **Expiry, D5-style (35a/F3):**
`basic_real_fft<float | double, scaling::fixed>`, the pre-Stage-4 spelling, resolves to the
selected engine as a second type with identical code (+2,949 B x86-64 / +1,995 B M55 of duplicate
wrappers when both are instantiated); no in-tree writer remains after the #35 fix pass (the
capi's seam uses `detail::default_real_fft_policy_t`), and the resolution is tolerated for one
consumer cycle — until MuTap and MuTap-Max pin a tree containing #35 — then
`detail::floating_engine_of<Sample, scaling::fixed>` goes and the spelling becomes a
`static_assert`. Alternative kept on record: no default, consumers name the engine.

**D5. Float-I/O-on-double overloads: `[[deprecated]]` for one consumer cycle, then deleted.**
Not gated on AmbiTap, which is not on disk and keeps its own wrapper per README.

**D6. `fftsg.c` moves to `tests/reference/ooura/` at 2c and is deleted after both MuTap and
MuTap-Max pin a tree containing 2c. `readme.txt` stays at `third_party/ooura/readme.txt`
permanently** (settled in the wave-1 review; it is the license record for the derived code, and
one fixed path is what the port header's banner can cite).

**D7. Settled: `detail::split_radix_rdft`**, one house token (`real_fft`) for everything
consumer-facing, provenance per Part 5.

**D8. No explicit-instantiation escape hatch.** Verified ineffective for in-class definitions
and unnecessary at measured compile cost. Revisit only if MuTap's test build shows a
regression the numbers in Part 4 do not predict.

**D9 (new). fp-contraction is a written contract point.** Parity gate at `-ffp-contract=off`;
whether `tap::dsp` exports that flag to consumers is decided in Stage 1 and stated in `fft.h`.

**D10 (new). Float stays bit-identical to the vendored C.** The port reproduces the C's table
semantics for both precisions; no numeric change ships with the port.

**D11 (new). Performance is a ratcheted gate from Stage 1b onward.** Instruction counts per
scenario per QEMU leg at ±3%, `.text` ceilings per leg, baselines seeded from the vendored C
before the routing flip; wall-clock is informational only. A regression is a failing check
and an accepted regression is a written commit.
