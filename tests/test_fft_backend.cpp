// SPDX-License-Identifier: MIT
// Copyright 2025-2026 Timothy Place and the DspTap contributors.
//
// Parity oracle for the float32 FFT engines, typed over the engines this
// host can build (Stage 4 of docs/audit-fft-and-code-smells.md: the engine
// is a template parameter of basic_real_fft, so an accelerated engine and
// the split-radix engine are instantiated side by side in ONE binary and
// compared there, rather than parity being a property of the CI matrix):
//
//   - detail::split_radix_rdft<float>       everywhere (the port-vs-port row,
//                                           trivially exact, a live guard that
//                                           the class adds nothing);
//   - detail::accelerate_real_fft_f32       on the macOS leg (TAP_DSP_FFT_ACCELERATE),
//                                           the same-binary comparison the plan
//                                           asked for;
//   - detail::cmsis_real_fft_f32            on the Cortex-M55 QEMU leg
//                                           (TAP_DSP_FFT_CMSIS), where the
//                                           CMSIS-vs-split-radix rows have run
//                                           under emulation since that leg
//                                           landed (main's version compared
//                                           basic_real_fft<float> — CMSIS under
//                                           the define — to the reference);
//                                           since Stage 4 they are typed rows
//                                           beside the split-radix ones in one
//                                           binary. Never on a host and never on
//                                           hardware: recorded, not closed.
//
// Each row pins basic_real_fft<float, Engine> to the reference float engine —
// detail::split_radix_rdft<float> called directly, the C++20 transliteration
// of Ooura's rdft_f that is bit-identical to the C it replaced
// (tests/test_fft_parity_ooura.cpp, against the reference copy under
// tests/reference/ooura/) — bin-for-bin at the two certified geometries
// (512-pt canceller, 2048-pt suppressor analysis), and holds it to the
// alignment-stability and tonal-accuracy gates below at 512 / 2048 / 4096.
// A static_assert pins that the build's DEFAULT float engine is one of the
// rows, so real_fft32 is always covered.
//
// Before Stage 2b the reference here was the raw rdft_f of the then-vendored
// fftsg_float.c; 2b re-pointed it at the ported engine (a change of oracle
// only, because the Stage 2a gate holds the engine bit-identical to the C);
// since 2c the C is not in the shipping tree (Decision D6) and this file
// needs nothing from it.
//
// The reconciliation under test: both accelerated engines use the
// engineering convention exp(-i2*pi/N), while the contract is Ooura's
// exp(+i2*pi/N) with an unnormalized inverse. Each wrapper conjugates the
// imaginary bins and rescales so every intermediate spectrum matches the
// reference; if that ever drifts, the bin-for-bin comparison below fails
// loudly. The rest of the FFT contract (packing, the +i sign convention,
// Parseval, float-tracks-double) is covered by test_fft.cpp on the default
// engine.

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <numbers>
#include <string>
#include <type_traits>
#include <vector>

#include <gtest/gtest.h>

#include "support/signals.h"
#include "tap/dsp/fft.h"
#include "tap/dsp/fft/split_radix.h"

namespace {

    using split_radix_f = tap::dsp::detail::split_radix_rdft<float>;

    /// The engines this translation unit can instantiate: the split-radix
    /// engine always, plus the accelerated engine the build selected (its
    /// header is included by fft.h under the define, and its library is on
    /// the link line only then).
    using host_engines = ::testing::Types<split_radix_f
#if defined(TAP_DSP_FFT_ACCELERATE)
                                          ,
                                          tap::dsp::detail::accelerate_real_fft_f32
#elif defined(TAP_DSP_FFT_CMSIS)
                                          ,
                                          tap::dsp::detail::cmsis_real_fft_f32
#endif
                                          >;

    /// The default float engine is always one of the rows above.
    template <typename Engine, typename List>
    struct is_listed;
    template <typename Engine, typename... Engines>
    struct is_listed<Engine, ::testing::Types<Engines...>>
        : std::bool_constant<(std::is_same_v<Engine, Engines> || ...)> {};
    static_assert(is_listed<tap::dsp::default_real_fft_engine_t<float>, host_engines>::value,
                  "the build's default float engine must be one of the typed rows");

    struct engine_names {
        template <typename Engine>
        static std::string GetName(int) { // NOLINT(readability-identifier-naming): gtest's required spelling
            if constexpr (std::is_same_v<Engine, split_radix_f>) {
                return "split_radix";
            }
#if defined(TAP_DSP_FFT_ACCELERATE)
            else if constexpr (std::is_same_v<Engine, tap::dsp::detail::accelerate_real_fft_f32>) {
                return "vdsp";
            }
#elif defined(TAP_DSP_FFT_CMSIS)
            else if constexpr (std::is_same_v<Engine, tap::dsp::detail::cmsis_real_fft_f32>) {
                return "cmsis";
            }
#endif
            else {
                return "unknown";
            }
        }
    };

    // The float reference, independent of the class and of any engine
    // parameter: the split-radix engine called directly. Bit-identical to
    // the rdft_f this struct called before the Stage 2b flip.
    struct ooura_ref {
        split_radix_f m_engine;
        explicit ooura_ref(int n)
            : m_engine(static_cast<size_t>(n)) {}
        void forward(float* a) { m_engine.forward_inplace(a); }
    };

    std::vector<float> broadband(int n, unsigned seed) {
        std::vector<float>         x(static_cast<size_t>(n));
        tap::dsp::test::xorshift32 rng(seed);
        for (auto& v : x) {
            v = rng.next_unit_f() * 0.1f;
        }
        return x;
    }

    constexpr int k_certified_geometries[] = {512, 2048};
    constexpr int k_stability_geometries[] = {512, 2048, 4096};

    template <typename Engine>
    class fft_backend_parity : public ::testing::Test {};
    TYPED_TEST_SUITE(fft_backend_parity, host_engines, engine_names);

    // Forward transform of the class over this engine must match the
    // split-radix reference to single-precision rounding, in the packed
    // layout and sign convention the whole DSP chain assumes. The test keeps
    // its name: the reference is Ooura's transform, bit for bit.
    TYPED_TEST(fft_backend_parity, ForwardMatchesOoura) {
        for (const int n : k_certified_geometries) {
            const auto         x   = broadband(n, 0x9E3779B9u);
            std::vector<float> ref = x;
            ooura_ref          oref(n);
            oref.forward(ref.data());

            tap::dsp::basic_real_fft<float, TypeParam> fft(static_cast<size_t>(n));
            std::vector<float>                         got = x;
            fft.forward_inplace(got.data());

            float peak = 0.0f;
            for (float v : ref) {
                peak = std::fmax(peak, std::fabs(v));
            }
            float max_err = 0.0f;
            for (int i = 0; i < n; ++i) {
                max_err = std::fmax(max_err, std::fabs(got[static_cast<size_t>(i)] - ref[static_cast<size_t>(i)]));
            }
            // De-risk harness measured 1.3e-7..2.0e-7 relative; 5e-6 is a wide,
            // stable band that still catches a convention or scaling regression.
            EXPECT_LT(max_err, 5e-6f * peak) << "N=" << n << " peak=" << peak;
        }
    }

    // Round trip through the class reproduces the input: locks the inverse
    // scaling (the N/2 factor in the CMSIS path, the x0.25 in vDSP's)
    // independently of the forward.
    TYPED_TEST(fft_backend_parity, RoundTripReproducesInput) {
        for (const int n : k_certified_geometries) {
            tap::dsp::basic_real_fft<float, TypeParam> fft(static_cast<size_t>(n));
            const auto                                 x = broadband(n, 0x1234567u);

            std::vector<float> spectrum(static_cast<size_t>(n));
            std::vector<float> back(static_cast<size_t>(n));
            fft.forward(x.data(), spectrum.data());
            fft.inverse(spectrum.data(), back.data());

            for (int i = 0; i < n; ++i) {
                EXPECT_NEAR(back[static_cast<size_t>(i)], x[static_cast<size_t>(i)], 2e-5f)
                    << "N=" << n << " sample " << i;
            }
        }
    }

    template <typename Engine>
    class fft_alignment_stability : public ::testing::Test {};
    TYPED_TEST_SUITE(fft_alignment_stability, host_engines, engine_names);

    // THE GATE fft_backend_parity CANNOT BE. Parity builds one engine, in one
    // process, at whatever address the allocator happened to pick, and compares
    // it to Ooura once — so an engine that returns DIFFERENT bits depending on
    // where its scratch buffers land passes it every time, by comparing an
    // arbitrary draw.
    //
    // That is not hypothetical: Apple's vDSP dispatches on 64-byte alignment
    // and the two paths disagree (tap/MuTap#31 — at N=2048, 140/60 across 200
    // processes on an M1, tracked exactly to whether the split-complex halves
    // were 64-byte aligned). Downstream, MuTap's residual suppressor turned
    // that into a ~97 dB swing in a certified compliance row, and it read as
    // CI flake for two weeks.
    //
    // So: build many engines with the heap deliberately shifted between
    // constructions, run identical input through each, and require every
    // output to be BIT-identical — not close, identical. An engine is a
    // function of its input or it is not usable as one.
    TYPED_TEST(fft_alignment_stability, OutputDoesNotDependOnBufferAddress) {
        for (const int n : k_stability_geometries) {
            const auto x = broadband(n, 0x9E3779B9u);

            std::vector<float> reference;
            // Odd, growing spacer allocations walk the engine's buffers through
            // every alignment class the allocator can produce. Kept alive so each
            // engine really does land somewhere new.
            std::vector<std::vector<char>> spacers;

            for (int trial = 0; trial < 32; ++trial) {
                spacers.emplace_back(static_cast<size_t>(4 * trial + 1), char{});

                tap::dsp::basic_real_fft<float, TypeParam> fft(static_cast<size_t>(n));
                std::vector<float>                         got = x;
                fft.forward_inplace(got.data());

                if (trial == 0) {
                    reference = got;
                    continue;
                }
                for (size_t i = 0; i < got.size(); ++i) {
                    // memcmp rather than ==, so the comparison is over the exact
                    // bits (== would call +0.0 and -0.0 equal, and any NaN
                    // unequal to itself).
                    ASSERT_EQ(std::memcmp(&got[i], &reference[i], sizeof(float)), 0)
                        << "N=" << n << " trial " << trial << ", bin " << i
                        << ": the forward transform is not a function of its input — the engine is dispatching on "
                           "something other than the data (buffer alignment is the known case; see k_align_bytes in "
                           "fft/backends/accelerate.h)";
                }
            }
        }
    }

    // A copy must keep the property: basic_real_fft is held by value in the
    // chains that use it, so a copy whose storage landed at a different
    // alignment would reintroduce the bug at the copy site.
    TYPED_TEST(fft_alignment_stability, CopiesAgreeWithTheirSource) {
        for (const int n : k_stability_geometries) {
            const auto x = broadband(n, 0x9E3779B9u);

            tap::dsp::basic_real_fft<float, TypeParam> original(static_cast<size_t>(n));
            std::vector<float>                         from_original = x;
            original.forward_inplace(from_original.data());

            std::vector<std::vector<char>> spacers;
            for (int trial = 0; trial < 8; ++trial) {
                spacers.emplace_back(static_cast<size_t>(8 * trial + 3), char{});

                tap::dsp::basic_real_fft<float, TypeParam> copy =
                    original; // NOLINT(performance-unnecessary-copy-initialization)
                std::vector<float> got = x;
                copy.forward_inplace(got.data());

                for (size_t i = 0; i < got.size(); ++i) {
                    ASSERT_EQ(std::memcmp(&got[i], &from_original[i], sizeof(float)), 0)
                        << "N=" << n << " copy " << trial << " disagrees with its source at bin " << i;
                }
            }
        }
    }

    template <typename Engine>
    class fft_tonal_accuracy : public ::testing::Test {};
    TYPED_TEST_SUITE(fft_tonal_accuracy, host_engines, engine_names);

    // THE GATE THAT MAKES accelerate.h's k_skew_bytes SELF-CHECKING.
    //
    // Transform a pure ON-BIN tone — one large bin, N-1 numerically empty ones
    // — and require the engine to track a double-precision reference through
    // the empty ones. That is the material where the two vDSP kernels diverge:
    // 64-byte-aligned buffers put the MEDIAN bin 65% away from truth at N=2048,
    // while the skewed placement the wrapper uses, and Ooura, both track it to
    // float epsilon (~1e-07).
    //
    // Asserting on the MEDIAN is the point. Max error is useless here — divide
    // any float32 noise by a numerically empty bin and it saturates, for every
    // engine — and peak-normalized absolute error is worse than useless,
    // because dividing by the one loud bin reports a spectrum that is wrong
    // everywhere except its peak as excellent. Neither can see the defect this
    // test exists to catch, which is exactly why fft_backend_parity (broadband
    // material, peak-normalized bound) passed throughout.
    //
    // If a future SDK dispatches differently and the skew stops selecting the
    // accurate kernel, this fails loudly instead of degrading in silence.
    TYPED_TEST(fft_tonal_accuracy, EmptyBinsTrackADoubleReference) {
        for (const int n : k_stability_geometries) {
            const double pi = std::numbers::pi;

            // Exactly on bin n/16, so no leakage lifts the empty bins above the
            // noise floor and hides the effect.
            std::vector<float> x(static_cast<size_t>(n));
            for (int i = 0; i < n; ++i) {
                x[static_cast<size_t>(i)] =
                    static_cast<float>(0.5 * std::sin(2.0 * pi * (n / 16.0) * i / static_cast<double>(n)));
            }

            // float -> double is exact, so both engines see identical input values.
            std::vector<double>              ref(x.begin(), x.end());
            tap::dsp::basic_real_fft<double> dfft(static_cast<size_t>(n));
            dfft.forward_inplace(ref.data());

            std::vector<float>                         got = x;
            tap::dsp::basic_real_fft<float, TypeParam> ffft(static_cast<size_t>(n));
            ffft.forward_inplace(got.data());

            std::vector<double> rel;
            rel.reserve(ref.size());
            for (size_t i = 0; i < ref.size(); ++i) {
                const double d = std::fabs(static_cast<double>(got[i]) - ref[i]);
                rel.push_back(std::fabs(ref[i]) > 0.0 ? d / std::fabs(ref[i]) : 0.0);
            }
            std::sort(rel.begin(), rel.end());
            const double median = rel[rel.size() / 2];

            // Measured ~1.1e-07..1.9e-07 on every good engine (Ooura, and vDSP at
            // the skewed placement); 0.65 and worse on the 64-byte-aligned vDSP
            // kernel. 1e-5 sits ~50x above the good case and orders below the bad
            // one, so it discriminates without being brittle.
            EXPECT_LT(median, 1e-5)
                << "N=" << n << ": the median bin is " << median
                << " away from the double-precision reference, i.e. most of this spectrum is wrong. "
                << "On Apple that is the signature of the 64-byte-aligned vDSP kernel — check "
                << "k_skew_bytes in fft/backends/accelerate.h and whether the SDK's dispatch rule "
                   "has moved.";
        }
    }

} // namespace
