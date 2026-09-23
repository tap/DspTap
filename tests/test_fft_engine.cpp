// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// THE ENGINE-PARAMETER CONTRACT (Stage 4 of docs/audit-fft-and-code-smells.md,
// Decision D4): what fft.h promises about basic_real_fft's second template
// argument — the build-selected default, the ABI tag derived from it, every
// spelling that must keep compiling, the per-engine size range with its
// supports_size predicate, and the k_is_shareable trait. Compile-time facts
// are pinned as static_asserts (so a violation fails the build of the
// battery, on every leg); constructions at the range bounds run for real.
//
// Where this runs, and what each leg proves:
//   - every host and the M4 / M4F / M33 legs: the default float engine is the
//     split-radix engine, the tag is fft_split_radix, the split-radix range
//     is 4 … 2^30 and N = 16 and 8192 are accepted;
//   - the M55 leg (TAP_DSP_FFT_CMSIS): the default is the CMSIS engine, the
//     tag is fft_cmsis, its range is 32 … 4096, N = 16 and 8192 are REJECTED
//     by supports_size while the split-radix engine named explicitly in the
//     same binary accepts them, and construction at 32 and 4096 works;
//   - the macOS leg (TAP_DSP_FFT_ACCELERATE): the default is vDSP, the tag is
//     fft_vdsp, its range 4 … 2^20.
// The precondition itself (TAP_EXPECTS) is a debug assertion: every CI
// battery is Release / MinSizeRel, so the death test at the end compiles in
// only for a local Debug configure with death tests available, and the
// release-mode statement is supports_size, which is what the rest of this
// file pins. A stub engine with a narrowed range (narrow_engine) proves on
// the host that the class reads the ENGINE's numbers rather than a
// hard-coded 4, which is the mechanism the M55 leg relies on.

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <type_traits>
#include <vector>

#include <gtest/gtest.h>

#include "support/signals.h"
#include "tap/dsp/fft.h"
#include "tap/dsp/fft/split_radix.h"
#include "tap/dsp/log_mel.h"
#include "tap/dsp/pvoc.h"

#if defined(__cpp_rtti)
#include <typeinfo>
#endif

namespace {

    using split_radix_f = tap::dsp::detail::split_radix_rdft<float>;
    using split_radix_d = tap::dsp::detail::split_radix_rdft<double>;
    using fixed_q15     = tap::dsp::detail::fixed_point_rdft<std::int16_t, tap::dsp::scaling::fixed>;
    using fixed_q31_bfp = tap::dsp::detail::fixed_point_rdft<std::int32_t, tap::dsp::scaling::block_floating>;

    // ------------------------------------------------------------------------
    // The default selection and the tag follow the one build define.
    // ------------------------------------------------------------------------
#if defined(TAP_DSP_FFT_CMSIS)
    using expected_float_default           = tap::dsp::detail::cmsis_real_fft_f32;
    constexpr const char* k_expected_tag   = "fft_cmsis";
    constexpr std::size_t k_expected_min   = 32;
    constexpr std::size_t k_expected_max   = 4096;
    constexpr bool        k_expected_share = false;
#elif defined(TAP_DSP_FFT_ACCELERATE)
    using expected_float_default           = tap::dsp::detail::accelerate_real_fft_f32;
    constexpr const char* k_expected_tag   = "fft_vdsp";
    constexpr std::size_t k_expected_min   = 4;
    constexpr std::size_t k_expected_max   = std::size_t{1} << 20;
    constexpr bool        k_expected_share = false;
#else
    using expected_float_default           = split_radix_f;
    constexpr const char* k_expected_tag   = "fft_split_radix";
    constexpr std::size_t k_expected_min   = 4;
    constexpr std::size_t k_expected_max   = std::size_t{1} << 30;
    constexpr bool        k_expected_share = true;
#endif

    static_assert(std::is_same_v<tap::dsp::default_real_fft_engine_t<float>, expected_float_default>,
                  "the float default is the engine the build define selects");
    static_assert(std::is_same_v<tap::dsp::default_real_fft_engine_t<double>, split_radix_d>,
                  "double always defaults to the split-radix engine");
    static_assert(std::is_same_v<tap::dsp::real_fft32::engine, expected_float_default>);
    static_assert(std::is_same_v<tap::dsp::real_fft::engine, split_radix_d>);
    static_assert(std::is_same_v<tap::dsp::real_fft32, tap::dsp::basic_real_fft<float, expected_float_default>>,
                  "the one-argument spelling IS the explicit spelling of the default engine");

    TEST(fft_engine, AbiTagNamesTheSelectedDefault) {
        EXPECT_STREQ(tap::dsp::k_real_fft_abi_tag, k_expected_tag);
    }

#if defined(__cpp_rtti)
    // The tag is in the TYPE, not only in a string: typeid's name (mangled
    // on GCC / Clang, decorated on MSVC) spells the inline namespace for the
    // class and for the two classes that embed it by value. This is the
    // property the loader sees (nm evidence in docs/fft-design.md, "Stage 4").
    TEST(fft_engine, TypeNamesCarryTheTag) {
        EXPECT_NE(std::strstr(typeid(tap::dsp::real_fft32).name(), k_expected_tag), nullptr)
            << typeid(tap::dsp::real_fft32).name();
        EXPECT_NE(std::strstr(typeid(tap::dsp::pvoc32).name(), k_expected_tag), nullptr)
            << typeid(tap::dsp::pvoc32).name();
        EXPECT_NE(std::strstr(typeid(tap::dsp::log_mel32).name(), k_expected_tag), nullptr)
            << typeid(tap::dsp::log_mel32).name();
        // And the tag is the selection, not a constant: the string for a
        // different default never appears.
        for (const char* other : {"fft_split_radix", "fft_cmsis", "fft_vdsp"}) {
            if (std::strcmp(other, k_expected_tag) != 0) {
                EXPECT_EQ(std::strstr(typeid(tap::dsp::real_fft32).name(), other), nullptr) << other;
            }
        }
    }
#endif

    // ------------------------------------------------------------------------
    // Every spelling in use keeps compiling, with its meaning.
    // ------------------------------------------------------------------------
    using legacy_float  = tap::dsp::basic_real_fft<float, tap::dsp::scaling::fixed>;
    using legacy_double = tap::dsp::basic_real_fft<double, tap::dsp::scaling::fixed>;
    static_assert(std::is_same_v<legacy_float::engine, tap::dsp::real_fft32::engine>,
                  "scaling::fixed on a floating Sample names the selected engine (the pre-Stage-4 spelling)");
    static_assert(std::is_same_v<legacy_double::engine, split_radix_d>);
    static_assert(!std::is_same_v<legacy_float, tap::dsp::real_fft32>,
                  "... and is a distinct type from the one-argument form, as documented");
    static_assert(
        std::is_same_v<tap::dsp::real_fft_q15, tap::dsp::basic_real_fft<std::int16_t, tap::dsp::scaling::fixed>>);
    static_assert(std::is_same_v<tap::dsp::real_fft_q15::engine, fixed_q15>);
    static_assert(std::is_same_v<tap::dsp::real_fft_q31_bfp::engine, fixed_q31_bfp>);
    static_assert(std::is_same_v<tap::dsp::basic_real_fft<float, split_radix_f>::engine, split_radix_f>,
                  "an engine named explicitly is the engine");

    // The concept admits the four engines and nothing that is not one.
    static_assert(tap::dsp::real_fft_engine<split_radix_f, float>);
    static_assert(tap::dsp::real_fft_engine<split_radix_d, double>);
    static_assert(tap::dsp::real_fft_engine<expected_float_default, float>);
    static_assert(tap::dsp::real_fft_engine<fixed_q15, std::int16_t>);
    static_assert(tap::dsp::real_fft_engine<fixed_q31_bfp, std::int32_t>);
    static_assert(!tap::dsp::real_fft_engine<int, float>);
    static_assert(!tap::dsp::real_fft_engine<tap::dsp::scaling::fixed, float>,
                  "the legacy spelling is not an engine; the class resolves it before the concept is checked");
    static_assert(!tap::dsp::real_fft_engine<split_radix_d, float>, "an engine is typed over its sample");

    TEST(fft_engine, LegacySpellingComputesWhatTheDefaultComputes) {
        constexpr std::size_t n = 512;
        const auto            x = tap::dsp::test::random_signal<float>(n, 0x9E3779B9u);
        tap::dsp::real_fft32  a(n);
        legacy_float          b(n);
        std::vector<float>    via_a = x;
        std::vector<float>    via_b = x;
        a.forward_inplace(via_a.data());
        b.forward_inplace(via_b.data());
        EXPECT_EQ(std::memcmp(via_a.data(), via_b.data(), n * sizeof(float)), 0);
    }

    // ------------------------------------------------------------------------
    // Size ranges: the stated numbers, per engine, and the class re-exports
    // the engine's rather than its own.
    // ------------------------------------------------------------------------
    static_assert(split_radix_f::k_min_size == 4 && split_radix_f::k_max_size == (std::size_t{1} << 30));
    static_assert(split_radix_d::k_min_size == 4 && split_radix_d::k_max_size == (std::size_t{1} << 30));
    static_assert(tap::dsp::real_fft::k_min_size == split_radix_d::k_min_size
                  && tap::dsp::real_fft::k_max_size == split_radix_d::k_max_size);
    static_assert(tap::dsp::real_fft32::k_min_size == k_expected_min
                      && tap::dsp::real_fft32::k_max_size == k_expected_max,
                  "the float default's range is the selected engine's");
    static_assert(tap::dsp::real_fft_q15::k_min_size == 4 && tap::dsp::real_fft_q15::k_max_size == 65536);
    static_assert(tap::dsp::real_fft_q31::k_min_size == 4 && tap::dsp::real_fft_q31::k_max_size == 65536);
    static_assert(tap::dsp::real_fft_q15_bfp::k_max_size == 65536 && tap::dsp::real_fft_q31_bfp::k_max_size == 65536);

    // The M55 requirement, as compile-time facts on every leg: 16 and 8192
    // are rejected by the CMSIS default and accepted by the split-radix
    // engine named explicitly in the same binary.
    using split_radix_fft32 = tap::dsp::basic_real_fft<float, split_radix_f>;
    static_assert(split_radix_fft32::supports_size(16) && split_radix_fft32::supports_size(8192));
#if defined(TAP_DSP_FFT_CMSIS)
    static_assert(!tap::dsp::real_fft32::supports_size(16) && !tap::dsp::real_fft32::supports_size(8192),
                  "CMSIS-DSP's arm_rfft_fast_init_f32 handles 32 … 4096 only");
    static_assert(tap::dsp::real_fft32::supports_size(32) && tap::dsp::real_fft32::supports_size(4096));
#else
    static_assert(tap::dsp::real_fft32::supports_size(16) && tap::dsp::real_fft32::supports_size(8192));
#endif
    static_assert(!tap::dsp::real_fft32::supports_size(0) && !tap::dsp::real_fft32::supports_size(2)
                  && !tap::dsp::real_fft32::supports_size(12) && !tap::dsp::real_fft32::supports_size(4095));

    /// An engine that differs from the split-radix engine only in a narrower
    /// stated range: the class's predicate (and, in a debug build, its
    /// precondition) must follow the engine's numbers. This is the mechanism
    /// the CMSIS leg relies on, exercised where CMSIS cannot be built.
    struct narrow_engine : split_radix_f {
        using split_radix_f::split_radix_f;
        static constexpr std::size_t k_min_size = 32;
        static constexpr std::size_t k_max_size = 4096;
    };
    using narrow_fft32 = tap::dsp::basic_real_fft<float, narrow_engine>;
    static_assert(tap::dsp::real_fft_engine<narrow_engine, float>);
    static_assert(narrow_fft32::k_min_size == 32 && narrow_fft32::k_max_size == 4096);
    static_assert(!narrow_fft32::supports_size(16) && !narrow_fft32::supports_size(8192));
    static_assert(narrow_fft32::supports_size(32) && narrow_fft32::supports_size(4096));

    template <typename Fft>
    void expect_supports_size_is_the_interval() {
        for (std::size_t n = 0; n <= (std::size_t{1} << 17); ++n) {
            const bool expected = std::has_single_bit(n) && n >= Fft::k_min_size && n <= Fft::k_max_size;
            ASSERT_EQ(Fft::supports_size(n), expected) << "n=" << n;
        }
    }

    TEST(fft_engine, SupportsSizeIsThePowerOfTwoInterval) {
        expect_supports_size_is_the_interval<tap::dsp::real_fft>();
        expect_supports_size_is_the_interval<tap::dsp::real_fft32>();
        expect_supports_size_is_the_interval<split_radix_fft32>();
        expect_supports_size_is_the_interval<narrow_fft32>();
        expect_supports_size_is_the_interval<tap::dsp::real_fft_q15>();
        expect_supports_size_is_the_interval<tap::dsp::real_fft_q31_bfp>();
    }

    /// Acceptance at a size, not a floor: construct, forward, scaled inverse,
    /// and require the input back within a loose bound (the per-N floors are
    /// pinned in test_fft.cpp and test_fft_oracle.cpp; a wrong engine or a
    /// broken table at a new size is off by orders of magnitude).
    template <typename Fft>
    void expect_constructs_and_round_trips(std::size_t n) {
        ASSERT_TRUE(Fft::supports_size(n)) << "n=" << n;
        Fft                fft(n);
        const auto         x = tap::dsp::test::random_signal<float>(n, 0x2545F491u, 0.5);
        std::vector<float> spectrum(n);
        std::vector<float> back(n);
        fft.forward(x.data(), spectrum.data());
        fft.inverse(spectrum.data(), back.data());
        float worst = 0.0f;
        for (std::size_t i = 0; i < n; ++i) {
            worst = std::fmax(worst, std::fabs(back[i] - x[i]));
        }
        EXPECT_LT(worst, 1e-3f) << "n=" << n << " (acceptance bound; the floors are pinned elsewhere)";
    }

    // The bounds are constructed for real on every leg. The float default's
    // upper bound is capped at 4096 (what the emulated legs' data regions and
    // the CMSIS range have in common; the host sweeps reach 65536 in the
    // oracle), the explicit split-radix engine runs the M55 requirement's
    // 16 and 8192 everywhere.
    TEST(fft_engine, ConstructsAtTheRangeBounds) {
        expect_constructs_and_round_trips<tap::dsp::real_fft32>(tap::dsp::real_fft32::k_min_size);
        expect_constructs_and_round_trips<tap::dsp::real_fft32>(
            std::min<std::size_t>(tap::dsp::real_fft32::k_max_size, 4096));
        expect_constructs_and_round_trips<split_radix_fft32>(16);
        expect_constructs_and_round_trips<split_radix_fft32>(8192);
        expect_constructs_and_round_trips<narrow_fft32>(32);
        expect_constructs_and_round_trips<narrow_fft32>(4096);
    }

    // ------------------------------------------------------------------------
    // Shareability: the header's numbers, and the compiler-checked half.
    // ------------------------------------------------------------------------
    static_assert(tap::dsp::real_fft::k_is_shareable, "double: tables built in the constructor, const transforms");
    static_assert(tap::dsp::real_fft32::k_is_shareable == k_expected_share,
                  "float: true on the split-radix engine, false on the two scratch-carrying backends");
    static_assert(std::is_same_v<decltype(tap::dsp::real_fft32::k_is_shareable), const bool>);
    static_assert(!tap::dsp::real_fft_q15::k_is_shareable && !tap::dsp::real_fft_q15_bfp::k_is_shareable,
                  "Q15: the int32 work buffer behind the in-place int16 API");
    static_assert(tap::dsp::real_fft_q31::k_is_shareable && tap::dsp::real_fft_q31_bfp::k_is_shareable,
                  "Q31: no mutable state during a transform");

    /// A shareable engine's transforms are callable on a const object; a
    /// non-shareable one's are not (so the trait cannot quietly drift from
    /// the signatures).
    template <typename Engine, typename Sample>
    constexpr bool k_transforms_are_const = requires(const Engine& e, Sample* a) {
        e.forward_inplace(a);
        e.inverse_inplace(a);
    };
    static_assert(k_transforms_are_const<split_radix_f, float> && k_transforms_are_const<split_radix_d, double>);
    static_assert(k_transforms_are_const<fixed_q31_bfp, std::int32_t>);
    static_assert(!k_transforms_are_const<fixed_q15, std::int16_t>);
    static_assert(k_transforms_are_const<expected_float_default, float> == k_expected_share);

    TEST(fft_engine, ShareabilityIsTheHeadersNumber) {
        EXPECT_TRUE(tap::dsp::real_fft::k_is_shareable);
        EXPECT_EQ(tap::dsp::real_fft32::k_is_shareable, k_expected_share);
        EXPECT_FALSE(tap::dsp::real_fft_q15::k_is_shareable);
        EXPECT_TRUE(tap::dsp::real_fft_q31::k_is_shareable);
    }

    // ------------------------------------------------------------------------
    // The precondition itself: a debug assertion (TAP_EXPECTS). Compiled in
    // only where it can fire and be caught — a Debug configure on a host with
    // death tests — since every CI battery is Release / MinSizeRel and the
    // QEMU legs have no death tests.
    // ------------------------------------------------------------------------
#if !defined(NDEBUG) && defined(GTEST_HAS_DEATH_TEST) && GTEST_HAS_DEATH_TEST
    TEST(fft_engine_debug, ConstructionOutsideTheRangeAssertsInADebugBuild) {
        EXPECT_DEATH({ narrow_fft32 fft(16); }, "");
        EXPECT_DEATH({ narrow_fft32 fft(8192); }, "");
        EXPECT_DEATH({ tap::dsp::real_fft fft(12); }, ""); // not a power of two
        EXPECT_DEATH({ tap::dsp::real_fft_q15 fft(2); }, "");
    }
#endif

} // namespace
