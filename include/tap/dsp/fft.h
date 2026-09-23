/// @file fft.h
/// @brief Real FFT with a fixed numeric contract: double, float, Q15 and Q31 profiles over an explicit engine.
// SPDX-License-Identifier: MIT
// Copyright 2025-2026 Timothy Place and the DspTap contributors.
//
// Extracted from the Tap family DSP libraries (MuTap's adaptive-filtering FFT
// and AmbiTap's binaural convolution FFT), which each carried a byte-identical
// copy of the vendored Ooura transform under a diverging wrapper. This is the
// consolidated wrapper: one numeric contract (Ooura's, carried since Stage 2b
// by the bit-identical C++20 port), one place to add a faster engine. See
// README.md for the provenance and migration notes.
//
// Four profiles share the contract (packing, exp(+i), unnormalized inverse):
// double (the golden model) and float (the embedded floating profile) run an
// ENGINE that is a template parameter of the class (Stage 4 of
// docs/audit-fft-and-code-smells.md, Decision D4): the split-radix engine in
// fft/split_radix.h, the C++20 transliteration of Ooura's rdft that Stage 2a
// landed bit-identical to the vendored C, is the default for double always
// and for float unless the build selected an accelerated engine
// (fft/backends/cmsis.h on the Cortex-M55, fft/backends/accelerate.h on
// Apple). std::int16_t (Q15) and std::int32_t (Q31) run the int32
// fixed-point kernel in fft/fixed_point.h, whose transforms return an
// exponent in place of a floating scale (Stage 3b).

#pragma once

#include <algorithm>
#include <bit>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <type_traits>
#include <utility>
#include <vector>

#include "tap/dsp/detail/expects.h"
#include "tap/dsp/fft/fixed_point.h"
#include "tap/dsp/fft/split_radix.h"

// Header-only. No vendored C is compiled into what ships: the split-radix
// engine replaced Ooura's rdft at Stage 2b, and Stage 2c moved the C out of
// the shipping tree (docs/audit-fft-and-code-smells.md, Part 3); Decision D6
// then deleted the reference copy the bit-identity gate compiled, and the
// engine's identity to the C is pinned since as output fingerprints
// (tests/test_fft_split_radix_fingerprint.cpp). tap::dsp is a pure INTERFACE
// target unless TAP_DSP_FFT_CMSIS is on, in which case it links the CMSIS-DSP
// objects (root CMakeLists.txt).

// ----------------------------------------------------------------------------
// THE ONE PLACE the build selects the float default engine and, with it, the
// ABI tag (Stage 4). AT MOST ONE of the two defines may be set (they are
// mutually exclusive), and each applies ONLY to float — double always defaults
// to the split-radix engine, the golden model, so the double-precision
// reference battery is unaffected:
//   TAP_DSP_FFT_CMSIS       CMSIS-DSP Helium on the bare-metal Cortex-M55
//                           (fft/backends/cmsis.h; N = 32 … 4096 only)
//   TAP_DSP_FFT_ACCELERATE  Apple's vDSP (Accelerate) on macOS
//                           (fft/backends/accelerate.h)
// Each accelerated engine re-presents the split-radix engine's EXACT numeric
// contract (Ooura's: same packed layout, exp(+i) sign convention,
// unnormalized inverse), so every intermediate spectrum matches the default
// build to float epsilon and the whole float32 test battery stays a valid
// oracle (tests/test_fft_backend.cpp, typed over the engines the host can
// build). Against the vendored C, to which the split-radix engine is
// bit-identical: on the M55 icount key the C executes 1.81x (N = 512) /
// 1.97x (N = 2048) the instructions of the CMSIS build over the whole
// ratchet scenario (bench/README.md, run 35844483811); the "~3x faster on
// Apple Silicon" figure is MuTap's transform-only measurement on the C
// (tap/MuTap#31), not re-measured here (docs/fft-design.md).
//
// TAP_DSP_FFT_ABI names the inline namespace this header opens on tap::dsp
// (below) from the selection: fft_split_radix, fft_cmsis or fft_vdsp. The
// selection is visible in the mangled name of every class defined inside it,
// which is the point — see "The ABI tag" in the class docstring.
// ----------------------------------------------------------------------------
#if defined(TAP_DSP_FFT_CMSIS) && defined(TAP_DSP_FFT_ACCELERATE)
#error "TAP_DSP_FFT_CMSIS and TAP_DSP_FFT_ACCELERATE are mutually exclusive"
#endif
#if defined(TAP_DSP_FFT_CMSIS)
#include "tap/dsp/fft/backends/cmsis.h"
#define TAP_DSP_FFT_ABI fft_cmsis
#define TAP_DSP_FFT_ABI_NAME "fft_cmsis"
#elif defined(TAP_DSP_FFT_ACCELERATE)
#include "tap/dsp/fft/backends/accelerate.h"
#define TAP_DSP_FFT_ABI fft_vdsp
#define TAP_DSP_FFT_ABI_NAME "fft_vdsp"
#else
#define TAP_DSP_FFT_ABI fft_split_radix
#define TAP_DSP_FFT_ABI_NAME "fft_split_radix"
#endif

namespace tap::dsp {

    /// What basic_real_fft requires of an engine (Stage 4): constructed from
    /// the transform size (which allocates and builds everything), copyable,
    /// two noexcept in-place transforms over Sample* (returning void for the
    /// floating engines, the exponent for the fixed-point one), size(), and
    /// three contract numbers — the supported size range, as the bounds of a
    /// power-of-two interval, and whether one object may be transformed
    /// through by two threads at once. Satisfied by detail::split_radix_rdft,
    /// detail::accelerate_real_fft_f32, detail::cmsis_real_fft_f32 and
    /// detail::fixed_point_rdft.
    template <typename Engine, typename Sample>
    concept real_fft_engine = std::copyable<Engine> && std::constructible_from<Engine, std::size_t>
                              && requires(Engine& e, const Engine& ce, Sample* a) {
                                     { ce.size() } noexcept -> std::same_as<std::size_t>;
                                     { e.forward_inplace(a) } noexcept;
                                     { e.inverse_inplace(a) } noexcept;
                                     requires std::is_same_v<decltype(Engine::k_min_size), const std::size_t>;
                                     requires std::is_same_v<decltype(Engine::k_max_size), const std::size_t>;
                                     requires std::is_same_v<decltype(Engine::k_is_shareable), const bool>;
                                 };

    namespace detail {

        /// The engine a floating profile runs when the caller names none:
        /// the split-radix engine, except float under one of the two build
        /// defines above. Specialized in exactly one place (here); every
        /// other header reads it through default_real_fft_engine_t.
        template <std::floating_point Sample>
        struct default_real_fft_engine {
            using type = split_radix_rdft<Sample>;
        };
#if defined(TAP_DSP_FFT_CMSIS)
        template <>
        struct default_real_fft_engine<float> {
            using type = cmsis_real_fft_f32;
        };
#elif defined(TAP_DSP_FFT_ACCELERATE)
        template <>
        struct default_real_fft_engine<float> {
            using type = accelerate_real_fft_f32;
        };
#endif

    } // namespace detail

    /// The engine basic_real_fft<Sample> (one argument) runs for a floating
    /// Sample in this build: detail::split_radix_rdft<Sample>, or for float
    /// the accelerated engine the build selected (the block above).
    template <std::floating_point Sample>
    using default_real_fft_engine_t = typename detail::default_real_fft_engine<Sample>::type;

    namespace detail {

        /// basic_real_fft's second template argument when the caller writes
        /// none: the selected engine for the floating profiles, scaling::fixed
        /// for the fixed-point ones (one parameter list, two meanings; the
        /// class docstring, "Two parameter lists in one").
        template <typename Sample>
        struct default_real_fft_policy {
            using type = scaling::fixed;
        };
        template <std::floating_point Sample>
        struct default_real_fft_policy<Sample> {
            using type = default_real_fft_engine_t<Sample>;
        };
        template <typename Sample>
        using default_real_fft_policy_t = typename default_real_fft_policy<Sample>::type;

        /// The engine a floating profile's second argument names: itself, or
        /// — when it is scaling::fixed, the spelling every profile shared
        /// before Stage 4 — the selected engine.
        template <typename Sample, typename Policy>
        struct floating_engine_of {
            using type = Policy;
        };
        template <typename Sample>
        struct floating_engine_of<Sample, scaling::fixed> {
            using type = default_real_fft_engine_t<Sample>;
        };

    } // namespace detail

} // namespace tap::dsp

// The ABI tag: an inline namespace on tap::dsp named for the build's float
// default engine (TAP_DSP_FFT_ABI above). Everything whose object layout
// follows that selection lives inside it — basic_real_fft here, and the
// classes that embed a default-engine basic_real_fft by value (pvoc.h,
// log_mel.h) — so its mangled names differ between two images built with
// different defaults and cannot be coalesced into each other (audit F4;
// the class docstring below has the full statement).
namespace tap::dsp::inline TAP_DSP_FFT_ABI {

    /// The tag this translation unit was built under, as text ("fft_split_radix",
    /// "fft_cmsis" or "fft_vdsp"); the tests pin it against the build define.
    inline constexpr const char* k_real_fft_abi_tag = TAP_DSP_FFT_ABI_NAME;

    /// Real FFT with a fixed numeric contract, parameterized over the sample
    /// type and, for the floating profiles, the engine: float or double (this
    /// primary template) run the engine named by the second argument, which
    /// defaults to default_real_fft_engine_t<Sample>; std::int16_t and
    /// std::int32_t (the Q15 and Q31 fixed-point profiles, the specialization
    /// below) take the Scaling policy in the second position instead and
    /// return an exponent from every transform.
    ///
    /// Two parameter lists in one (Stage 4, Decision D4). The second
    /// parameter means "engine" for a floating Sample and "scaling policy"
    /// for a fixed-point one, so that every spelling in use keeps compiling
    /// and keeps its meaning: basic_real_fft<double> and <float> (the
    /// selected engine), basic_real_fft<float, detail::split_radix_rdft<float>>
    /// (an engine named explicitly — on a macOS or M55 build, the split-radix
    /// engine beside the accelerated default in the same binary, which is
    /// what makes engine parity a test rather than a CI-matrix property),
    /// basic_real_fft<std::int16_t> / <std::int32_t, scaling::block_floating>
    /// (the fixed-point profiles, unchanged; note that `int` is std::int32_t
    /// on every CI target, so basic_real_fft<int> compiles and IS the Q31
    /// profile — the docstrings say std::int32_t and mean it), and
    /// basic_real_fft<float, scaling::fixed> / <double, scaling::fixed>, the spelling all four
    /// profiles shared before Stage 4 (the capi's generic seam writes it):
    /// for a floating Sample, scaling::fixed in the second position names the
    /// selected engine. That legacy spelling is a distinct type from
    /// basic_real_fft<float> (same layout, same code, different template
    /// arguments — measured +2,949 bytes of .text on x86-64 g++ -O2 and
    /// +1,995 on the M55 when both are instantiated, the class wrappers
    /// duplicated over one shared engine); nothing in DspTap writes it since
    /// the #35 fix pass (the capi's generic seam uses
    /// detail::default_real_fft_policy_t<Sample>) and nothing in MuTap or
    /// MuTap-Max ever did. EXPIRY, D5-style: the resolution
    /// (detail::floating_engine_of<Sample, scaling::fixed>) is tolerated for
    /// one consumer cycle — until MuTap and MuTap-Max have pinned a tree
    /// containing #35 — and then becomes a static_assert naming the
    /// one-argument form. Any other second argument on a floating Sample
    /// must satisfy real_fft_engine<Engine, Sample> (static_assert).
    ///
    /// Engine contract numbers, read from the engine and re-exported here:
    ///   - k_min_size / k_max_size, the power-of-two size range, and
    ///     supports_size(n), the constexpr predicate over it. Construction
    ///     requires supports_size(size) (TAP_EXPECTS: a debug assertion and
    ///     nothing in a release build, the house precondition style,
    ///     STYLE.md §4 and detail/expects.h). Stated so nobody reads more
    ///     into it: in a release build a size outside the range is a
    ///     precondition violation, i.e. undefined behaviour (for the CMSIS
    ///     engine a HardFault on the first transform; for the others whatever
    ///     the arithmetic does past its tables), and there is deliberately no
    ///     defined fallback in the transforms. supports_size is therefore the
    ///     MANDATORY gate wherever N comes from configuration rather than a
    ///     constant: the capi's dsptap_fft_create applies it per profile, and
    ///     a consumer's config path must (docs/fft-design.md, "MuTap bump
    ///     checklist"). `SupportsSizeIsThePowerOfTwoInterval` pins the
    ///     predicate on every leg. Per engine:
    ///       split-radix (double, float)   4 … 2^30  (the int-indexing bound
    ///                                     of Ooura's arithmetic; bit identity
    ///                                     to the C was gated to 2^20 until
    ///                                     D6 and is fingerprinted to 65536,
    ///                                     the oracle sweeps to 65536)
    ///       vDSP (float, Apple)           4 … 2^20  (fft/backends/accelerate.h)
    ///       CMSIS-DSP (float, Cortex-M55) 32 … 4096 (CMSIS's own table; below
    ///                                     32 or above 4096 the library's init
    ///                                     fails, which the constructor now
    ///                                     checks instead of ignoring — audit
    ///                                     Part 13; N = 4 hard-faulted)
    ///       fixed point (Q15, Q31)        4 … 65536
    ///     `EngineRangesAreTheStatedNumbers`, `SupportsSizeIsThePowerOfTwoInterval`,
    ///     `ConstructsAtTheRangeBounds` (tests/test_fft_engine.cpp).
    ///   - k_is_shareable: whether two threads may run transforms through ONE
    ///     object concurrently. True for the split-radix engine (its tables
    ///     are built in the constructor and its transforms are const) and for
    ///     Q31 (no mutable state during a transform); false for the two
    ///     accelerated engines (per-object scratch) and for Q15 (the int32
    ///     work buffer behind the in-place int16 API). Shareable IMPLIES
    ///     const: a shareable engine's transforms are const, and the class
    ///     static_asserts that, so that half of the trait is
    ///     compiler-checked. The converse is not asserted — an engine may
    ///     declare const transforms over mutable scratch and say
    ///     k_is_shareable = false, which is honest — and is a convention the
    ///     tests pin: for the six shipped instantiations const and shareable
    ///     coincide (`ShareabilityIsTheHeadersNumber`, tests/test_fft_rt.cpp
    ///     and tests/test_fft_engine.cpp). This class's own transforms stay
    ///     non-const on every profile (audit N11: the public API's constness
    ///     must not depend on the selected engine) — the trait, not the
    ///     signature, is the statement.
    ///
    /// The ABI tag (Stage 4, audit F4). basic_real_fft<float>'s layout follows
    /// the selected engine (its m_engine IS the engine), and so does the
    /// layout of every class that holds a basic_real_fft<float> by value
    /// (basic_pvoc<float>, basic_log_mel<float>, MuTap's fdaf<float>, …).
    /// Their member functions are weak (COMDAT) symbols in every image that
    /// instantiates them, and a dynamic loader may coalesce weak definitions
    /// across images — macOS dyld binds every image's WEAK_DEF references to
    /// one chosen definition; on ELF a plugin's references resolve against
    /// the executable's exported definitions first — so two images built
    /// with different defaults, loaded into one process, would run one
    /// image's code over the other image's layout. With the engine a
    /// template argument, basic_real_fft<float>'s own symbols already differ
    /// between the two builds (the argument is in the mangled name). The
    /// embedders' do not, which is what the tag is for: fft.h, pvoc.h and
    /// log_mel.h define those classes inside `inline namespace
    /// TAP_DSP_FFT_ABI` (fft_split_radix / fft_cmsis / fft_vdsp), so their
    /// mangled names carry the selection and the loader sees two unrelated
    /// symbols. Measured with arm-none-eabi-nm on one translation unit built
    /// twice for the M55 (docs/fft-design.md, "Stage 4"). What the tag does
    /// not change: unqualified and tap::dsp-qualified names resolve exactly
    /// as before (an inline namespace's members are members of the enclosing
    /// namespace for lookup, `using tap::dsp::basic_real_fft;` included);
    /// and the fixed-point profiles carry the tag too, because a partial
    /// specialization lives beside its primary template — their layout does
    /// not depend on the selection, so for them the tag prevents a coalescing
    /// that would have been harmless, at the cost of one instantiation per
    /// differently-built image and a link error rather than a silent merge if
    /// two such images ever exchange a real_fft_q15 across their boundary.
    /// A consumer whose own class embeds a basic_real_fft<float> by value has
    /// the same exposure and closes it the same way, in its own namespace
    /// (MuTap: a follow-up on its bump, recorded in the design note).
    ///
    /// Routing (Stage 2b, generalized at Stage 4):
    ///   - double  -> detail::split_radix_rdft<double> (fft/split_radix.h) by
    ///                default. The golden model.
    ///   - float   -> detail::split_radix_rdft<float> by default, unless the
    ///                build defines TAP_DSP_FFT_CMSIS (CMSIS-DSP Helium,
    ///                bare-metal Cortex-M55) or TAP_DSP_FFT_ACCELERATE (Apple
    ///                vDSP), in which case the default is that engine, which
    ///                re-presents the same contract to float epsilon
    ///                (tests/test_fft_backend.cpp, typed over the engines the
    ///                host can build — same-binary on macOS).
    ///   - Q15/Q31 -> detail::fixed_point_rdft (the specialization below).
    /// The split-radix engine is the C++20 transliteration of Ooura's rdft
    /// and is BIT-IDENTICAL to the C it replaced for both precisions
    /// (Decision D10: gated against the reference C from Stage 2a until D6
    /// deleted it, pinned since by tests/test_fft_split_radix_fingerprint.cpp),
    /// so the flip changed no output bit of any consumer built at default
    /// fp-contraction (measured on every CI platform) or with clang and FMA;
    /// the one measured exception is a g++ x86-64 build with -march (FMA),
    /// which moves the float profile's last bits at N >= 1024 (the
    /// fp-contraction policy below, D9). tests/test_fft_routing.cpp pins
    /// that this class's output is byte-identical to the engine's, for the
    /// split-radix engine on every build (named explicitly) and for whatever
    /// the default resolves to.
    ///
    /// FFT size must be a power of 2 inside the engine's range (above), fixed
    /// at construction. Workspace (bit-reversal and trig tables) is allocated
    /// AND BUILT in the constructor (the vendored C built its tables lazily
    /// on the first transform; the engine does not, audit item F6), so the
    /// first transform costs what every later one costs; the transforms
    /// themselves are noexcept and allocation-free, so they are safe on a
    /// real-time audio thread (tests/test_fft_rt.cpp). No alignment
    /// requirement on the data pointer. NaN propagates to every bin (no
    /// data-dependent branches). Latency 0. Copyable, and a copy is
    /// bit-identical to its source; the object is held by value in every
    /// consumer.
    ///
    /// Packing after a forward transform of N real samples (N/2 + 1 bins):
    ///   - bin[0].real    = data[0]   (DC;      imag is zero, not stored)
    ///   - bin[N/2].real  = data[1]   (Nyquist; imag is zero, not stored)
    ///   - bin[k].real    = data[2k], bin[k].imag = data[2k+1]  for 1 <= k < N/2
    ///
    /// Sign convention: the forward transform is A[k] = sum_j a[j]*W^(jk) with
    /// W = exp(+2*pi*i/N) — the imaginary parts are CONJUGATED relative to the
    /// engineering-convention DFT (exp(-2*pi*i/N)). Spectral products (e.g.
    /// fast convolution, FDAF regressor accumulation) are unaffected as long
    /// as every operand uses this class; conjugate if importing spectra
    /// computed elsewhere.
    ///
    /// The raw inverse is unnormalized: inverse_inplace() must be followed by
    /// a 2/N scaling for a round trip, which inverse() applies for you.
    ///
    /// Noise floor (docs/fft-design.md, contract table; measured 2026-09-23
    /// on x86-64 Linux, GCC 13.3.0 and Clang 18.1.3 -O3, glibc 2.39, the
    /// engine as routed here): double against the compensated-summation DFT
    /// oracle, relative 2-norm error of the forward on full-scale uniform
    /// noise at N = 256, 1.85e-16, pinned at 4x (7.4e-16) for the libm and
    /// fp-contraction spread across hosts
    /// (`fft_oracle_floor.DoubleForwardTracksCompensatedDft`); float against
    /// double, same metric at N = 512 on the split-radix engine itself,
    /// 1.12e-7 (the audit's Part 6 N2 probe value of 1.105e-7, re-measured on
    /// what ships; the two runs differ in material, not in engine), pinned at
    /// 2x (2.25e-7) (`RealFftCrossPrecision.FloatEngineTracksDoubleAtN512`),
    /// and < 1e-6 at N = 1024 through basic_real_fft (`FloatTracksDouble`).
    ///
    /// fp-contraction policy (Decision D9; measurements in docs/fft-design.md,
    /// "The fp-contraction policy"). tap::dsp does NOT export an
    /// -ffp-contract setting to its consumers and this header sets none: the
    /// engine's arithmetic is compiled under whatever contraction the
    /// consumer's compiler applies, alike for every instantiation in that
    /// build. The bit-identity gate against the reference C ran at
    /// -ffp-contract=off on both sides (Stage 2a until D6; the fingerprints
    /// that replaced it are built the same way), and at default flags the
    /// identity was MEASURED to hold as well on every CI platform (0 ulp on
    /// linux, windows and macOS arm64 and on the Cortex-M4, M4F, M33 and M55
    /// legs, three of them VFMA; identical bench checksums between the C and
    /// the engine on every bench key while the C was built there; MuTap's
    /// fingerprint rows byte-identical through the flip, float rows
    /// included, with no flag set anywhere), because the
    /// engine's statements are textually the C's and each compiler fuses both
    /// sides alike. That is an observation, not a guarantee: g++ on x86-64
    /// built with -march (FMA) is measured to fuse the two sides differently
    /// (float moves by a few ulp at N >= 1024; clang with the same flags is
    /// identical) and is not claimed. That is the first experiment, the
    /// engine-vs-reference-C gate (retired at D6): bit-identical at
    /// -ffp-contract=off on both sides and, at default flags, on every CI
    /// platform including MSVC and the four QEMU legs.
    ///
    /// The second experiment, stronger, measured on tap/DspTap#34 (Stage 6,
    /// which touches neither fft.h nor fft/), is a fingerprint A/B of pvoc
    /// and log_mel through this class, main vs branch, run on g++ and on the
    /// M33 leg and NOT on MSVC or AppleClang: under g++ -O3
    /// -march=x86-64-v3 at the default -ffp-contract=fast, the output of
    /// basic_real_fft depends on the TRANSLATION-UNIT CONTEXT, not only on
    /// the flags — an edit to unrelated code in the same TU (log_mel.h's
    /// constructor) changed GCC's inlining of cftmdl2 into the split-radix
    /// engine's cftrec4 / cftleaf (float cftrec4 went from 194 to 168
    /// vfmadd instructions, double from 224 to 254) and moved pvoc's float
    /// output bits while every pvoc function was instruction-identical;
    /// appending ~30 unrelated lines to the TU made the outputs identical
    /// again. The double codegen moved in that experiment too, so "double
    /// does not move under g++ with FMA" (the 2b measurement) is an
    /// observation, not a guarantee.
    ///
    /// Why no export: an INTERFACE -ffp-contract=off would reach every
    /// consumer translation unit that includes this header and would
    /// pessimize the VFMA / FMA targets the float profile exists for (the
    /// M55, Apple arm64) for the whole of that code, in exchange for a
    /// cross-compiler bit reproducibility that libm's last-bit cos/sin
    /// differences between glibc, newlib, UCRT and Apple already deny across
    /// hosts. A consumer that needs bit reproducibility across its own
    /// compilers sets -ffp-contract=off on its own targets; the contract this
    /// header makes is the one above.
    ///
    /// Thread rule: one transform at a time per object unless k_is_shareable
    /// (above) says otherwise for the engine in use.
    template <typename Sample, typename Policy = detail::default_real_fft_policy_t<Sample>>
    class basic_real_fft {
        static_assert(std::is_same_v<Sample, float> || std::is_same_v<Sample, double>,
                      "basic_real_fft supports float and double (split radix, or the engine named) and std::int16_t "
                      "/ std::int32_t (Q15 / Q31)");

      public:
        /// The engine this instantiation runs (see the class docstring).
        using engine = typename detail::floating_engine_of<Sample, Policy>::type;
        static_assert(real_fft_engine<engine, Sample>,
                      "basic_real_fft<float | double, Engine>: Engine must satisfy tap::dsp::real_fft_engine (or be "
                      "scaling::fixed, the pre-Stage-4 spelling of the selected engine)");
        static_assert(std::is_void_v<decltype(std::declval<engine&>().forward_inplace(std::declval<Sample*>()))>,
                      "a floating engine's transforms return void (the fixed-point engine returns the exponent)");

        /// The engine's power-of-two size range, as contract numbers.
        static constexpr std::size_t k_min_size = engine::k_min_size;
        static constexpr std::size_t k_max_size = engine::k_max_size;
        /// Whether two threads may run transforms through one object at once.
        static constexpr bool k_is_shareable = engine::k_is_shareable;
        static_assert(
            !k_is_shareable || requires(const engine& e, Sample* a) {
                e.forward_inplace(a);
                e.inverse_inplace(a);
            }, "a shareable engine's transforms are const");

        /// True exactly for the sizes the constructor accepts: a power of two
        /// in [k_min_size, k_max_size]. The release-mode counterpart of the
        /// constructor's precondition.
        [[nodiscard]] static constexpr bool supports_size(std::size_t n) noexcept {
            return n >= k_min_size && n <= k_max_size && std::has_single_bit(n);
        }

        /// @pre supports_size(size) (TAP_EXPECTS).
        explicit basic_real_fft(std::size_t size)
            : m_size(static_cast<int>(size))
            , m_engine(checked(size)) {}

        [[nodiscard]] std::size_t size() const noexcept { return static_cast<std::size_t>(m_size); }
        [[nodiscard]] std::size_t num_bins() const noexcept { return static_cast<std::size_t>(m_size / 2 + 1); }

        /// In-place forward FFT: time-domain Sample[size] -> packed spectrum Sample[size].
        void forward_inplace(Sample* data) noexcept { m_engine.forward_inplace(data); }

        /// In-place inverse FFT: packed spectrum Sample[size] -> time-domain Sample[size],
        /// UNSCALED — multiply by 2/size for a normalized round trip.
        void inverse_inplace(Sample* data) noexcept { m_engine.inverse_inplace(data); }

        /// Out-of-place forward FFT. Output may alias input.
        void forward(const Sample* input, Sample* output) noexcept {
            copy(input, output);
            forward_inplace(output);
        }

        /// Out-of-place inverse FFT, scaled by 2/size so forward() -> inverse()
        /// reproduces the input. Output may alias input.
        void inverse(const Sample* input, Sample* output) noexcept {
            copy(input, output);
            inverse_inplace(output);
            const Sample scale = Sample(2) / static_cast<Sample>(m_size);
            for (int i = 0; i < m_size; ++i) {
                output[i] *= scale;
            }
        }

        /// Float-I/O convenience on the DOUBLE engine: run the double-precision
        /// transform over float buffers (copy in, transform, copy out), so a
        /// caller holding float data can use the double FFT without maintaining
        /// its own double staging buffer — e.g. AmbiTap's binaural HRTF analysis,
        /// which wants double-precision spectra from float impulse responses.
        /// Only the double instantiation offers these; basic_real_fft<float>
        /// already takes float in the same-type forward()/inverse() above. These
        /// allocate a staging buffer (a setup-time path), unlike the noexcept,
        /// allocation-free in-place transforms.
        ///
        /// DEPRECATED (Decision D5, audit item F7): removed after one consumer
        /// cycle, i.e. once MuTap and MuTap-Max have pinned a tree containing
        /// this deprecation. No consumer on disk calls them (grepped at Stage
        /// 2b: DspTap's tools/capi, MuTap, MuTap-Max); the replacement is the
        /// caller's own double staging buffer and the same-type transforms.
        [[deprecated("basic_real_fft<double>::forward(const float*, float*) allocates per call and is removed "
                     "after one consumer cycle (Decision D5); stage through a double buffer instead")]]
        void forward(const float* input, float* output)
            requires std::is_same_v<Sample, double>
        {
            std::vector<double> buf(static_cast<size_t>(m_size));
            for (int i = 0; i < m_size; ++i) {
                buf[static_cast<size_t>(i)] = static_cast<double>(input[i]);
            }
            forward_inplace(buf.data());
            for (int i = 0; i < m_size; ++i) {
                output[i] = static_cast<float>(buf[static_cast<size_t>(i)]);
            }
        }

        /// Float-I/O inverse on the double engine, scaled by 2/size like inverse().
        /// DEPRECATED with forward(const float*, float*) above (Decision D5).
        [[deprecated("basic_real_fft<double>::inverse(const float*, float*) allocates per call and is removed "
                     "after one consumer cycle (Decision D5); stage through a double buffer instead")]]
        void inverse(const float* input, float* output)
            requires std::is_same_v<Sample, double>
        {
            std::vector<double> buf(static_cast<size_t>(m_size));
            for (int i = 0; i < m_size; ++i) {
                buf[static_cast<size_t>(i)] = static_cast<double>(input[i]);
            }
            inverse_inplace(buf.data());
            const double scale = 2.0 / static_cast<double>(m_size);
            for (int i = 0; i < m_size; ++i) {
                output[i] = static_cast<float>(buf[static_cast<size_t>(i)] * scale);
            }
        }

      private:
        /// The constructor's precondition, checked before the engine is built
        /// (so the engine's own assertion is never the first to fire).
        static std::size_t checked(std::size_t n) noexcept {
            TAP_EXPECTS(supports_size(n));
            return n;
        }

        void copy(const Sample* input, Sample* output) noexcept {
            if (input != output) {
                for (int i = 0; i < m_size; ++i) {
                    output[i] = input[i];
                }
            }
        }

        int    m_size;
        engine m_engine;
    };

    /// The fixed-point profiles: Q15 (std::int16_t) and Q31 (std::int32_t)
    /// I/O over one int32 radix-4 kernel (fft/fixed_point.h), same packing and
    /// exp(+i) convention as the floating profiles, with the scale carried by
    /// an exponent instead of a floating mantissa. The second template
    /// argument is the Scaling policy here (the primary template's is the
    /// engine; "Two parameter lists in one" above); the engine is always
    /// detail::fixed_point_rdft<Sample, Scaling>, which has no build-selected
    /// alternative — the tag on this specialization is the primary's, not a
    /// statement about these profiles. Contract, as numbers (each pinned by
    /// the named test of the Stage 3b battery, tests/test_fft_fixed.cpp and
    /// the widened tests/test_fft.cpp):
    ///
    ///  - Exponent. Every transform returns e. Read the buffer as fractions of
    ///    full scale (Q0.15 / Q0.31) and let G be basic_real_fft<double> on the
    ///    same input read the same way: after forward_inplace, G's forward
    ///    result == data * 2^e (same packing); after inverse_inplace, G's
    ///    UNNORMALIZED inverse result == data * 2^e. The fixed-point inverse()
    ///    applies NO 2/N (unlike the floating profiles); a round trip
    ///    reconstructs x == out * 2^(e_fwd + e_inv + 1 - log2 N) under both
    ///    policies. `RoundTripReproducesInput`, `RoundTripReconstructsInputPerPolicy`,
    ///    `FixedForwardScaleIsExactlyXOverN`, `FixedForwardScaleIsExactlyXOverTwoN`.
    ///  - scaling::fixed (default): e == fixed_scaling_exponent(N)
    ///    == log2 N + fft_arith<Sample>::k_fixed_scaling_input_pre_shift in
    ///    both directions: log2 N for Q15 (output exactly X / N), log2 N + 1
    ///    for Q31 (X / 2N; the one-bit input pre-shift is the price of no
    ///    guard bits, fft_arith.h). `FixedExponentIsTheStatedConstant`,
    ///    `FixedInverseCarriesTheSameExponent`.
    ///  - scaling::block_floating: 0 <= e <= fixed_scaling_exponent(N),
    ///    data-dependent; the kernel never shifts more than the stage's
    ///    growth requires beyond the headroom present, and never consumes
    ///    the last bit. Brought to a common exponent, the block-floating
    ///    output agrees with the fixed one within a pinned bound (Q15 1.0
    ///    LSB, Q31 4.0 LSB measured; pinned 2 / 8). At the full exponent it
    ///    is bit-identical to the fixed output only when every stage shifted
    ///    the fixed amount (the full-scale patterns, pinned); in general the
    ///    two schedules reach the same exponent through different rounding
    ///    histories (an input with headroom at entry shifts less at the
    ///    first stage and catches up later) and differ by up to 4 LSB (Q31,
    ///    at the DC index; 2 LSB elsewhere) or 1 LSB (Q15, where the int32
    ///    history's few LSB32 survive the narrow only on a rounding tie),
    ///    measured and pinned at 2x. `BfpExponentIsWithinRange`,
    ///    `BfpMatchesFixedAfterShift`, `BfpAtTheFullExponentIsBitIdenticalToFixed`,
    ///    `BfpAtTheFullExponentAgreesWithFixedWithinPin`, `SilenceIsSilence`.
    ///  - Saturation. The int32 kernel performs no saturating operation for
    ///    any input under either policy: Q15 through the two guard bits of
    ///    the widened Q2.29 data, Q31 through the pre-shift (fixed) or the
    ///    headroom rule (block floating); the worst case is the packed pair
    ///    at full scale under a 45-degree twiddle. The Q15 output narrowing
    ///    (round-half-up, saturating) clamps at the rail exactly when the
    ///    true value is the rail: a full-scale Nyquist alternation lands on
    ///    32767.5 LSB and comes back as 32767, the pinned 0.5 LSB rail
    ///    shortfall. Largest deviation from the golden model over the
    ///    adversarial full-scale sweep, both directions: Q15 0.50 / 0.75 LSB
    ///    (fixed / block floating), Q31 4.25 / 15.99 LSB, at index 0 for
    ///    Q31; pinned at 1.0 / 1.5 / 8.5 / 32. `SaturationFreeWorstCaseDoesNotWrap`.
    ///  - Host identity. The fixed-point output is integer arithmetic over a
    ///    checksum-pinned table: for a fixed input it is one bit pattern on
    ///    every host, pinned per profile, policy, direction and N = 512 /
    ///    2048. `OutputFingerprintIsPinned`, `TwiddleTableChecksumIsPinned`.
    ///  - Noise floor (output-referred, against the double golden model on
    ///    the same quantized input; N = 256 / 512 / 2048, 0 to -60 dBFS; the
    ///    numbers are the `[ floor ]` rows `NoiseFloorTracksWelchModel`
    ///    prints, forward, white noise): Q15 fixed 0.26 - 0.30 LSB rms at
    ///    every N and level (the narrow's own rounding; per-bin SNR 69.1 /
    ///    66.5 / 60.2 dB at 0 dBFS); Q31 fixed 0.67 - 0.75 LSB rms (152.0 /
    ///    149.4 / 142.8 dB), level-independent, i.e. SNR falls 20 dB per
    ///    20 dB of level; block floating point keeps 84 - 91 dB (Q15) and
    ///    157 - 161 dB (Q31) at 0 dBFS and does not lose the low-level
    ///    signal (78 - 86 dB Q15, 154 - 156 dB Q31 at -40 dBFS). Welch's
    ///    variance model predicts 0.57 - 0.59 LSB32 for the kernel; the
    ///    measured/model ratio of 1.31 - 1.74 (Q31 fixed) is the
    ///    round-half-up bias, largest at index 0 under block floating point
    ///    (fft/fixed_point.h, "Honest limit"). `NoiseFloorTracksWelchModel`,
    ///    `RoundingBiasOnNegatedInputIsBounded`, `Q15TracksDouble`, `Q31TracksDouble`,
    ///    `Q15AndQ31AgreeToTheQ15Floor`.
    ///  - CMSIS-DSP compatibility (Decision D3): CMSIS documents its q15 /
    ///    q31 real FFTs as "downscaled by 2 for every stage", i.e. a forward
    ///    output of X / N with log2 N bits to upscale. The Q15 fixed forward
    ///    here is that same X / N; the Q31 fixed forward is X / 2N, one bit
    ///    below, because of the pre-shift; the inverse here is Ooura's
    ///    unnormalized inverse over 2^e, i.e. (1/N) sum X W^-jk divided by 2
    ///    (Q15) or 4 (Q31), which is not CMSIS's inverse table; block
    ///    floating point has no CMSIS analogue. The exponent, not a fixed
    ///    Q format per N, is the contract here; nothing of CMSIS's
    ///    behaviour was measured or reproduced, only its documentation read.
    ///  - Size k_min_size = 4 <= N <= k_max_size = 65536, a power of two,
    ///    fixed at construction; supports_size(N) says so and the constructor
    ///    requires it (TAP_EXPECTS). Q15 allocates an int32 work buffer of N
    ///    at construction (in-place API preserved at the caller's int16
    ///    buffer); Q31 transforms in place. Transforms are noexcept and
    ///    allocation-free, the object is copyable, there is no alignment
    ///    requirement. `TransformsAreNoexcept`, `ForwardInplaceAllocatesNothing`
    ///    (and the inverse and out-of-place forms),
    ///    `CopyProducesBitIdenticalOutput`, `OutOfPlaceIsCopyThenInPlace`.
    ///  - Thread rule, as the engine trait k_is_shareable (Stage 4): Q15
    ///    false — the in-place int16 API needs the per-object int32 work
    ///    buffer; Q31 true — its transforms touch the caller's buffer and the
    ///    tables only, and are const in the engine so the compiler checks it
    ///    (a recorded deviation from audit Part 7, which listed both fixed
    ///    profiles as shareable). This class's transforms are non-const on
    ///    every profile (N11). `ShareabilityIsTheHeadersNumber`.
    ///  - Every transform returns the exponent and is [[nodiscard]]: under
    ///    block floating point a discarded e is a silent scale error.
    ///  - Latency 0; no NaN or denormal behaviour to state (integer data).
    ///
    /// Per-profile noise floors and the Welch-model derivation are in
    /// docs/fft-design.md ("The fixed-point profiles") and the README
    /// profiles table.
    template <typename Sample, typename Scaling>
        requires fft_arith<Sample>::k_is_fixed_point
    class basic_real_fft<Sample, Scaling> {
      public:
        using engine = detail::fixed_point_rdft<Sample, Scaling>;
        static_assert(real_fft_engine<engine, Sample>, "the fixed-point engine satisfies the engine concept");

        /// The engine's power-of-two size range, as contract numbers.
        static constexpr std::size_t k_min_size = engine::k_min_size;
        static constexpr std::size_t k_max_size = engine::k_max_size;
        /// Whether two threads may run transforms through one object at once
        /// (Q31 true, Q15 false; the docstring's thread rule).
        static constexpr bool k_is_shareable = engine::k_is_shareable;
        static_assert(
            !k_is_shareable || requires(const engine& e, Sample* a) {
                e.forward_inplace(a);
                e.inverse_inplace(a);
            }, "a shareable engine's transforms are const");

        /// The constant exponent of scaling::fixed (and the upper bound of
        /// scaling::block_floating) for size n; the tests use it in place of
        /// magic numbers.
        [[nodiscard]] static constexpr int fixed_scaling_exponent(std::size_t n) noexcept {
            return engine::fixed_scaling_exponent(n);
        }

        /// True exactly for the sizes the constructor accepts: a power of two
        /// in [k_min_size, k_max_size].
        [[nodiscard]] static constexpr bool supports_size(std::size_t n) noexcept {
            return n >= k_min_size && n <= k_max_size && std::has_single_bit(n);
        }

        /// @pre supports_size(size) (TAP_EXPECTS).
        explicit basic_real_fft(std::size_t size)
            : m_engine(checked(size)) {}

        [[nodiscard]] std::size_t size() const noexcept { return m_engine.size(); }
        [[nodiscard]] std::size_t num_bins() const noexcept { return m_engine.size() / 2 + 1; }

        /// In-place forward FFT: Sample[size] -> packed spectrum Sample[size],
        /// scaled by 2^-e. @return e, the output's scale: under block floating
        /// point a discarded exponent is a silent scale error, hence nodiscard.
        [[nodiscard]] int forward_inplace(Sample* data) noexcept { return m_engine.forward_inplace(data); }

        /// In-place inverse FFT: packed spectrum -> Sample[size], the
        /// UNNORMALIZED inverse scaled by 2^-e. @return e.
        [[nodiscard]] int inverse_inplace(Sample* data) noexcept { return m_engine.inverse_inplace(data); }

        /// Out-of-place forward FFT: copy, then forward_inplace. Output may
        /// alias input. @return e.
        [[nodiscard]] int forward(const Sample* input, Sample* output) noexcept {
            copy(input, output);
            return forward_inplace(output);
        }

        /// Out-of-place inverse FFT: copy, then inverse_inplace. NO 2/N is
        /// applied (the exponent carries the scale, unlike the floating
        /// profiles' inverse()). Output may alias input. @return e.
        [[nodiscard]] int inverse(const Sample* input, Sample* output) noexcept {
            copy(input, output);
            return inverse_inplace(output);
        }

      private:
        static std::size_t checked(std::size_t n) noexcept {
            TAP_EXPECTS(supports_size(n));
            return n;
        }

        void copy(const Sample* input, Sample* output) noexcept {
            if (input != output) {
                std::copy_n(input, m_engine.size(), output);
            }
        }

        engine m_engine;
    };

    /// Double-precision real FFT — the desktop/golden-model profile.
    using real_fft = basic_real_fft<double>;

    /// Single-precision real FFT — the embedded real-time profile (Cortex-M55,
    /// Hexagon HVX), where hardware floating point is single-precision only;
    /// on the engine the build selected (default_real_fft_engine_t<float>).
    using real_fft32 = basic_real_fft<float>;

    /// Q15 real FFT (Q0.15 I/O), fixed scaling: e == log2 N, output exactly X / N.
    using real_fft_q15 = basic_real_fft<std::int16_t>;
    /// Q31 real FFT (Q0.31 I/O), fixed scaling: e == log2 N + 1, output exactly X / 2N.
    using real_fft_q31 = basic_real_fft<std::int32_t>;
    /// Q15 real FFT, block floating point: 0 <= e <= log2 N, data-dependent.
    using real_fft_q15_bfp = basic_real_fft<std::int16_t, scaling::block_floating>;
    /// Q31 real FFT, block floating point: 0 <= e <= log2 N + 1, data-dependent.
    using real_fft_q31_bfp = basic_real_fft<std::int32_t, scaling::block_floating>;

} // namespace tap::dsp::inline TAP_DSP_FFT_ABI
