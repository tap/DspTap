# Audit: the FFT layer and its relatives across DspTap

*September 2026. Baseline at the time of the audit: `5ca3b1c`, builds warning-free with
`-DTAP_DSP_WERROR=ON`, 160/160 tests pass on Linux/Ooura.*

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
`_f` copies that is 76 global symbols per binary. Any other Ooura user in the same process
(Max hosts and third-party externals commonly vendor this exact file) is a duplicate-symbol or
silent-interposition hazard. The rename table in F1 exists *because* of this.

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
  or not; `m_ip/m_w` are allocated or not). Two translation units built with different defines
  and linked together are an ODR violation with a real crash mode. MuTap-Max builds through
  MuTap (which forces vDSP off) while the capi and any other consumer take DspTap's default
  (vDSP on for Apple): the hazard is one CMake cache away.
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

## Part 3 — Plan

Ordered so that each stage is independently landable, leaves the tree green, and reduces the
blast radius of the next one. Numeric contracts do not move until Stage 3, and Stage 3 is
gated on bit-exactness.

### Stage 0 — Hygiene that touches no numerics (one PR)
Fix `long m_n` (both files), `solve_dense` noexcept, NOTICE.md names, the dangling toolchain
references and the `MUTAP_` variable, delete the dead `TAP_DSP_CHANNEL_PARALLEL` macros (or
make something read them), brace the SMLALD tail loops, de-duplicate the `fft.h` paragraph,
add `<algorithm>`/`<utility>` to kaiser. Add `include/tap/dsp/detail/math.h` with `k_pi`
delegating to `std::numbers::pi`, `hann_periodic`, `db_from_power`/`db_from_amplitude`, and
migrate the four copies. Add a `tests/support/` for the copy-pasted test helpers.

### Stage 1 — One packed-spectrum view (one PR here, one in MuTap)
`include/tap/dsp/spectrum.h`: a non-owning view over the Ooura-packed buffer with `dc()`,
`nyquist()`, `re(k)`, `im(k)`, `bin(k)` returning a small complex value in the *engineering*
convention (the conjugation lives in exactly one place), `num_bins()`, and `power(k)`.
Migrate `pvoc.h`, `log_mel.h`, `test_pvoc.cpp`; then MuTap's five headers. After this the
packing contract has one definition and the FFT can change its internals without the
consumers knowing.

### Stage 2 — Preconditions and the `double` question in the substrate (one PR)
Introduce `TAP_EXPECTS` (STYLE.md already names it) and replace the bare asserts; give
every primitive a callable `valid()` or a checked factory so the capi stops re-deriving
rules. Decide and land the `sample_traits<double>` question (Decision D1 below); until it is
settled, a traits-based FFT has no golden profile.

### Stage 3 — Port Ooura's `rdft` path to a C++ template (the core PR)
*(Detailed design in Part 4; naming in Part 5.)*

`include/tap/dsp/detail/ooura_rdft.h`: `template <class Sample> class ooura_rdft`, a
mechanical transliteration of *only* the functions `rdft` reaches (`makewt`, `makeipt`,
`makect`, `bitrv2*`, `cftfsub/cftbsub`, `cftf1st/cftb1st`, `cftrec4`, `cftleaf`, `cftmdl1/2`,
`cftfx41`, `cftf161/162/081/082/040`, `cftb040`, `cftx020`, `rftfsub/rftbsub`), same operation
order, tables built in the constructor (fixes F6), all helpers private members. Ooura's terms
permit modification; the header keeps his copyright banner and NOTICE.md records the port.

Gate: a parity test that instantiates the template and the vendored C side by side and
requires **bit-identical** output for `double` and `float` at every power of two from 4 to
65536, on broadband, on-bin tone, impulse and DC material, forward and inverse. Two notes for
that test: (a) build both with the same `-ffp-contract` setting or the comparison is not
meaningful; (b) if a platform cannot be made bit-exact, the fallback is a documented
1-ulp bound, not a loosened tolerance. The existing `test_fft.cpp`, `test_fft_backend.cpp`,
and MuTap's compliance batteries are the second gate.

Then `basic_real_fft<double/float>` routes to the template, `fftsg_float.c` is deleted, and
`fftsg.c` moves to `tests/reference/` as the oracle (deleted after one release). `tap_dsp_fft`
exists only when the CMSIS backend is on; the library becomes header-only as documented.
This closes F1, F2, F3, F6, F7 (the overloads move to a free function in a non-RT header).

### Stage 4 — Backends as a customization point, not a preprocessor branch
`basic_real_fft<Sample, Engine = default_real_fft_engine_t<Sample>>`, where `Engine` is any
type with `init(n)`, `forward_inplace`, `inverse_inplace`. Ooura is `ooura_rdft<Sample>`;
vDSP and CMSIS become `include/tap/dsp/backends/{accelerate,cmsis}_rfft32.h`, included
only by the consumer or by the one build-selected alias line. The build define shrinks to
choosing the *default* engine and nothing else touches it; the class layout no longer depends
on it (closes F4); Ooura and vDSP can be instantiated in one binary, so backend parity
becomes an ordinary test rather than a CI matrix.

### Stage 5 — Fixed point
Extend `sample_traits` with the vocabulary an FFT needs (from/to double, `add`/`sub`,
`mul_coeff` with the documented rounding, named fraction-bit constants, a `halve()` or
shift-with-rounding). Then instantiate the same ported Ooura structure over an arithmetic
policy for `int16_t` and `int32_t`: Q15/Q31 data, twiddles in the traits' existing
coefficient formats (Q1.14 / Q1.30, which is exactly what a twiddle table wants), int32/int64
intermediates, fixed scaling per butterfly (radix-4 leaves `/4`, radix-8 `/8`, radix-2 `/2`) so
the output is deterministic and the headroom is a number in the header. Gate: the double
battery as oracle within the format's floor, plus the pinned-behaviour tests the fixed-point
substrate already uses. CMSIS `arm_rfft_q15/q31` can become an optional Stage 4 engine for
that profile later, once the house contract exists to test it against.

---

## Part 4 — The C++20 port in detail

### Scope
`rdft` reaches about 2,400 of `fftsg.c`'s 3,325 lines: `makewt`, `makeipt`, `makect`,
`bitrv2`, `bitrv2conj`, `bitrv216`, `bitrv216neg`, `bitrv208`, `bitrv208neg`, `cftfsub`,
`cftbsub`, `cftf1st`, `cftb1st`, `cftrec4`, `cftleaf`, `cftmdl1`, `cftmdl2`, `cftfx41`,
`cftf161`, `cftf162`, `cftf081`, `cftf082`, `cftf040`, `cftb040`, `cftx020`, `rftfsub`,
`rftbsub`. The DCT/DST family (`ddct`, `ddst`, `dfct`, `dfst`, `dctsub`, `dstsub`), the
complex entry point `cdft`, and the pthread/Win32 scaffolding (`cftrec4_th` and friends) are
not ported. `cftrec4` is recursive to depth log4(N); it is bounded and allocation-free and
stays recursive.

### Shape
- One class template, `template <std::floating_point Sample> class split_radix_rdft`, in
  `include/tap/dsp/fft/split_radix.h`, namespace `tap::dsp::detail`. (Name discussed in
  Part 5.) Every helper is a private static member function, so the 76 global symbols
  disappear and the rename table has no reason to exist. The concept replaces the
  `static_assert` enumeration.
- Members: `std::size_t m_n`, the bit-reversal table `m_ip` and the trig table `m_w`
  (`std::vector<Sample>`). Index types may be widened to `std::size_t` freely; they do not
  affect the floating-point bit pattern. Arithmetic order inside every butterfly is preserved
  exactly, because bit-exactness against the C is the gate.
- Public surface of the engine: constructor from size, `forward_inplace(Sample*) const`,
  `inverse_inplace(Sample*) const`, `size()`. `basic_real_fft` keeps the consumer-facing
  surface (see Part 1) and adds `std::span<Sample>` overloads, `[[nodiscard]]` on the size
  queries, `std::size_t` throughout with a single narrowing point, and `TAP_EXPECTS` on the
  power-of-two precondition. The raw `cdft`/`cdft_f` declarations leave the public header; a
  complex transform, if a consumer ever needs one, is its own class.

### Tables are built in the constructor, in double, once
This closes F6 (lazy first-call initialization on the RT path). It also fixes a defect in the
current float build that the audit surfaced only by reading the macro: `#define double float`
retargets the *locals* inside `makewt` too, so the float twiddles are computed from a
float-rounded `delta = atan(1)/nwh` and their argument error grows with the table index. The
port computes every twiddle in double and rounds once into `Sample`.

Consequence for the gate: the **double** port must be bit-identical to the vendored C; the
**float** port will differ from `fftsg_float.c` by design, because it is more accurate. The
float gate is therefore the existing float-tracks-double test plus a small ulp bound against
the double port, not bit identity against the old float build.

### Transforms become `const`
With immutable tables, `forward_inplace` and `inverse_inplace` are `const` member functions.
That states in the type what the docstring promises in prose, and it lets one plan be shared
across threads. The vDSP and CMSIS engines carry scratch buffers and stay non-const, which is
itself a useful signal about which backends are shareable.

### Two mechanical passes, one gate
1. Port with plain `+ - *` and prove bit-exactness at every power of two from 4 to 65536, on
   broadband, on-bin tone, impulse and DC material, forward and inverse. Build the C and the
   C++ with the same `-ffp-contract` setting or the comparison is meaningless; if a platform
   cannot be made bit-exact, the documented fallback is a 1-ulp bound, never a loosened
   tolerance.
2. Only then introduce the arithmetic policy that fixed point needs (Stage 5), under the same
   gate. Designing the policy before the port exists means designing it blind; the gate makes
   the second pass safe.

### Header-only, with one escape hatch
Header-only. The templates belong in headers because that is what the rest of the library
is, and because a consumer that instantiates only `float` then compiles only the float code,
which the current static library cannot do. `tap::dsp` becomes a true INTERFACE target and
the compiled-library explanations disappear from MuTap's and MuTap-Max's CMake. The CMSIS
backend still compiles C, so a small object library exists only when that option is on.

The cost is compile time: the 2,400 lines are heavily unrolled radix-8/16 leaves, and MuTap
pulls `fft.h` into every external and test translation unit through its umbrella header. The
escape hatch is the standard one: `extern template class split_radix_rdft<double>;` (and
`<float>`) behind an opt-in CMake option that adds one `.cpp` with the explicit
instantiations. Default off; measure first; turn on only if the numbers say so. No C++20
modules: the Max toolchain and the submodule consumers are not ready for them.

### Cleanup the port makes possible
- **Shareable plans.** MuTap's chains hold several FFTs of one size, each with its own tables.
  With immutable tables the plan can be a value that copies cheaply (explicit sharing, no
  global cache), which matters on the M55.
- **Independent oracle.** Every FFT test today compares Ooura to Ooura. A naive long-double
  DFT at small sizes is an oracle that is not the implementation under test.
- **Typed tests over engines.** With the engine as a template parameter (Stage 4), one typed
  suite runs the split-radix engine, vDSP and CMSIS in the same binary; backend parity stops
  being a CI-matrix property.
- **File layout.** Public `fft.h`; then `fft/split_radix.h`, `fft/spectrum.h` (the bin view
  from Stage 1), `fft/backends/accelerate.h`, `fft/backends/cmsis.h`. A backend header is
  included by whoever selects it, not by everyone.
- **Fixed-point twiddles fall out of the substrate.** Q1.14 and Q1.30 are already
  `sample_traits<int16_t/int32_t>::coeff`, so the fixed-point twiddle table needs no new
  format design.
- **Do not fold the port into the same PR as the bin view or the preconditions work.** The
  bit-exact gate is convincing only when the diff around it is boring.

---

## Part 5 — Naming: is this still "Ooura"?

What lands is a C++ port of one algorithm from Takuya Ooura's package (the split-radix real
DFT, "Fast Version III"), not a fork of the package. Two different things want names.

**Attribution stays, and stays prominent.** Ooura's terms permit modification and require
the copyright notice. The header keeps his banner, `NOTICE.md` records the port and what was
dropped, and the docs say "a C++ port of Takuya Ooura's `fftsg.c` `rdft`" in the first
paragraph. The house IP policy (implement from published literature, cite it) also wants the
provenance visible, so nothing should read as if the algorithm were new.

**Identifiers should say what the code is, not whose it was.** Recommendation:
- Engine class `split_radix_rdft<Sample>` in `include/tap/dsp/fft/split_radix.h`. Descriptive,
  greppable, and it cannot be confused with the many `ooura_fft.h` / `namespace ooura`
  headers in the wild (WebRTC, Chromium, assorted audio SDKs). Linkage collision is not the
  concern (everything is inside `tap::dsp` and the port removes the global symbols); reader
  confusion and search results are.
- The public API is unchanged: `basic_real_fft`, `real_fft`, `real_fft32`, and later
  `real_fft_q15` / `real_fft_q31`.
- The words "Ooura packing" and "Ooura contract" leave the consumer headers. They become
  "the DspTap packed spectrum" and "the DspTap real-FFT contract (inherited from Ooura's
  `rdft`)", and the Stage 1 spectrum view is where the layout is defined. Consumers stop
  naming the vendor of a layout they never see.
- "Ooura" remains in exactly three places: the attribution banner, `NOTICE.md`, and the
  bit-exact parity test, which compares against `tests/reference/ooura/fftsg.c` until that
  file is deleted (D6).
- Not recommended: a coined brand (`tapfft`, `taprfft`) that suggests a novel algorithm, or
  keeping `ooura_rdft` as the type name, which invites the "is this the same as the fork I
  already have" question the rename is meant to remove.

---

## Decisions to make before Stage 3 (the "shape" discussion)

**D1. Is `double` a sample format?** Today the substrate says no, on purpose, and CLAUDE.md
says double is every primitive's golden model. Options: (a) add `sample_traits<double>`
(float coefficients-as-double, double accumulation) and make it the golden profile for the
FIR substrate too, which re-baselines decimate's oracle; (b) keep the substrate float-golden
and let the FFT carry its own `fft_traits<Sample>` with a double specialization. (a) is one
rule for the whole repo; (b) is less churn. Recommendation: (a).

**D2. One algorithm or two for fixed point?** (a) The ported Ooura split-radix over an
arithmetic policy, one bit-reversal, one twiddle table, one test battery; fixed scaling only
(the depth-first recursion makes block floating point awkward). (b) A separate simpler
radix-2² DIF kernel for the integer profiles that could do block floating point. (a) keeps
the golden-model-plus-profiles story literally true; (b) buys BFP dynamic range at the cost of
two implementations. Recommendation: (a) first; (b) only if a consumer measures the need.

**D3. What is the fixed-point forward's scale?** Ooura's forward is an unnormalized sum and
overflows Q15 immediately. The fixed-point contract has to be either "forward scaled by 1/N,
inverse unscaled" or expose the scale as a number per profile. This is a documented contract
point, not an implementation detail, and it should be chosen so that a future CMSIS q15
engine can re-present it.

**D4. Engine as a template parameter or build-only?** Stage 4 proposes the template
parameter with a build-selected default. The cost is one more template argument visible in
error messages; the benefit is same-binary backend parity tests and no layout-by-define.
Recommendation: template parameter.

**D5. Keep or drop the float-I/O-on-double overloads?** No consumer in the three repos in
scope uses them. Recommendation: move to a free function `transform_as_double(...)` in a
non-RT header and delete from the class, after checking AmbiTap.

**D6. How long does `fftsg.c` stay?** Recommendation: as a test-only oracle for one release
after Stage 3 lands in the consumers, then deleted; the bit-exact parity test is what makes
the deletion safe.

**D7. The engine's name. Settled: `split_radix_rdft`.** Part 5 has the reasoning; Ooura stays
in the attribution banner, NOTICE.md and the parity test only. Rejected: `ooura_rdft`
(collides with the world's forks in every search) and a coined brand (hides provenance the
IP policy wants visible). Consumer-facing wording moves from "Ooura packing / Ooura
contract" to "the DspTap packed spectrum / real-FFT contract" as part of Stage 1.

**D8. Compile-time escape hatch.** Header-only by default; an opt-in explicit-instantiation
`.cpp` behind a CMake option only if measured compile time in MuTap's test build says so.
Recommendation: land header-only, measure, decide.
