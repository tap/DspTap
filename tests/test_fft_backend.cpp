// SPDX-License-Identifier: MIT
// Copyright 2025-2026 Timothy Place and the DspTap contributors.
//
// Parity oracle for the optional float32 FFT backends (TAP_DSP_FFT_CMSIS on the
// M55, TAP_DSP_FFT_ACCELERATE on Apple; see include/tap/dsp/fft.h and README.md).
// When one is active, basic_real_fft<float> routes through that backend; this
// test pins it to Ooura's rdft_f — the golden model — bin-for-bin at the two
// certified geometries (512-pt canceller, 2048-pt suppressor analysis).
//
// The reconciliation under test: both backends use the engineering convention
// exp(-i2*pi/N), while our contract is Ooura's exp(+i2*pi/N) with an
// unnormalized inverse. Each wrapper conjugates the imaginary bins and rescales
// so every intermediate spectrum matches Ooura; if that ever drifts, the
// bin-for-bin comparison below fails loudly.
//
// When no backend is defined (the default desktop/Hexagon/Linux build),
// basic_real_fft<float> already IS Ooura, so these become Ooura-vs-Ooura
// identity checks — trivially exact, and a live guard that the reference path
// keeps matching the wrapper. The rest of the FFT contract (packing, the +i
// sign convention, Parseval, float-tracks-double) is covered by test_fft.cpp,
// which exercises the same backend.

#include <cmath>
#include <cstddef>
#include <cstring>
#include <vector>

#include <gtest/gtest.h>

#include "tap/dsp/fft.h"

namespace {

    // Direct Ooura float reference, independent of the C++ wrapper's backend:
    // rdft_f is the extern "C" symbol regardless of TAP_DSP_FFT_CMSIS, so it is
    // always the golden model to compare the wrapper against.
    struct ooura_ref {
        int                m_n;
        std::vector<int>   m_ip;
        std::vector<float> m_w;
        explicit ooura_ref(int n)
            : m_n(n)
            , m_ip(2 + static_cast<size_t>(std::sqrt(static_cast<double>(n) / 2.0)) + 1)
            , m_w(static_cast<size_t>(n) / 2) {
            m_ip[0] = 0;
        }
        void forward(float* a) { rdft_f(m_n, 1, a, m_ip.data(), m_w.data()); }
    };

    std::vector<float> broadband(int n, unsigned seed) {
        std::vector<float> x(static_cast<size_t>(n));
        unsigned           s = seed;
        for (auto& v : x) {
            s ^= s << 13;
            s ^= s >> 17;
            s ^= s << 5;
            v = (static_cast<float>(s) / 2147483648.0f - 1.0f) * 0.1f;
        }
        return x;
    }

    class fft_backend_parity : public ::testing::TestWithParam<int> {};

    // Forward transform of the wrapper (CMSIS when TAP_DSP_FFT_CMSIS) must match
    // the Ooura reference to single-precision rounding, in the packed layout
    // and sign convention the whole DSP chain assumes.
    TEST_P(fft_backend_parity, ForwardMatchesOoura) {
        const int n = GetParam();

        const auto         x   = broadband(n, 0x9E3779B9u);
        std::vector<float> ref = x;
        ooura_ref          oref(n);
        oref.forward(ref.data());

        tap::dsp::basic_real_fft<float> fft(static_cast<size_t>(n));
        std::vector<float>              got = x;
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

    // Round trip through the wrapper reproduces the input: locks the inverse
    // scaling (the N/2 factor in the CMSIS path) independently of the forward.
    TEST_P(fft_backend_parity, RoundTripReproducesInput) {
        const int n = GetParam();

        tap::dsp::basic_real_fft<float> fft(static_cast<size_t>(n));
        const auto                      x = broadband(n, 0x1234567u);

        std::vector<float> spectrum(static_cast<size_t>(n));
        std::vector<float> back(static_cast<size_t>(n));
        fft.forward(x.data(), spectrum.data());
        fft.inverse(spectrum.data(), back.data());

        for (int i = 0; i < n; ++i) {
            EXPECT_NEAR(back[static_cast<size_t>(i)], x[static_cast<size_t>(i)], 2e-5f) << "sample " << i;
        }
    }

    INSTANTIATE_TEST_SUITE_P(CertifiedGeometries, fft_backend_parity, ::testing::Values(512, 2048));

    class fft_alignment_stability : public ::testing::TestWithParam<int> {};

    // THE GATE fft_backend_parity CANNOT BE. Parity builds one engine, in one
    // process, at whatever address the allocator happened to pick, and compares
    // it to Ooura once — so a backend that returns DIFFERENT bits depending on
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
    // output to be BIT-identical — not close, identical. A backend is a
    // function of its input or it is not usable as one.
    TEST_P(fft_alignment_stability, OutputDoesNotDependOnBufferAddress) {
        const int  n = GetParam();
        const auto x = broadband(n, 0x9E3779B9u);

        std::vector<float> reference;
        // Odd, growing spacer allocations walk the engine's buffers through
        // every alignment class the allocator can produce. Kept alive so each
        // engine really does land somewhere new.
        std::vector<std::vector<char>> spacers;

        for (int trial = 0; trial < 32; ++trial) {
            spacers.emplace_back(static_cast<size_t>(4 * trial + 1), char{});

            tap::dsp::basic_real_fft<float> fft(static_cast<size_t>(n));
            std::vector<float>              got = x;
            fft.forward_inplace(got.data());

            if (trial == 0) {
                reference = got;
                continue;
            }
            for (size_t i = 0; i < got.size(); ++i) {
                // memcmp rather than ==, so the comparison is over the exact
                // bits (== would call +0.0 and -0.0 equal, and any NaN
                // unequal to itself), and rather than std::bit_cast, which
                // needs C++20 and is not available on every leg's test build.
                ASSERT_EQ(std::memcmp(&got[i], &reference[i], sizeof(float)), 0)
                    << "trial " << trial << ", bin " << i << ": the forward transform is not a function of its input — "
                    << "the backend is dispatching on something other than the data (buffer alignment is the known "
                    << "case; see fft.h k_align_bytes)";
            }
        }
    }

    // A copy must keep the property: basic_real_fft is held by value in the
    // chains that use it, so a copy whose storage landed at a different
    // alignment would reintroduce the bug at the copy site.
    TEST_P(fft_alignment_stability, CopiesAgreeWithTheirSource) {
        const int  n = GetParam();
        const auto x = broadband(n, 0x9E3779B9u);

        tap::dsp::basic_real_fft<float> original(static_cast<size_t>(n));
        std::vector<float>              from_original = x;
        original.forward_inplace(from_original.data());

        std::vector<std::vector<char>> spacers;
        for (int trial = 0; trial < 8; ++trial) {
            spacers.emplace_back(static_cast<size_t>(8 * trial + 3), char{});

            tap::dsp::basic_real_fft<float> copy = original; // NOLINT(performance-unnecessary-copy-initialization)
            std::vector<float>              got  = x;
            copy.forward_inplace(got.data());

            for (size_t i = 0; i < got.size(); ++i) {
                ASSERT_EQ(std::memcmp(&got[i], &from_original[i], sizeof(float)), 0)
                    << "copy " << trial << " disagrees with its source at bin " << i;
            }
        }
    }

    INSTANTIATE_TEST_SUITE_P(CertifiedGeometries, fft_alignment_stability, ::testing::Values(512, 2048, 4096));

    class fft_tonal_accuracy : public ::testing::TestWithParam<int> {};

    // THE GATE THAT MAKES fft.h's k_skew_bytes SELF-CHECKING.
    //
    // Transform a pure ON-BIN tone — one large bin, N-1 numerically empty ones
    // — and require the backend to track a double-precision reference through
    // the empty ones. That is the material where the two vDSP kernels diverge:
    // 64-byte-aligned buffers put the MEDIAN bin 65% away from truth at N=2048,
    // while the skewed placement fft.h uses, and Ooura, both track it to float
    // epsilon (~1e-07).
    //
    // Asserting on the MEDIAN is the point. Max error is useless here — divide
    // any float32 noise by a numerically empty bin and it saturates, for every
    // backend — and peak-normalized absolute error is worse than useless,
    // because dividing by the one loud bin reports a spectrum that is wrong
    // everywhere except its peak as excellent. Neither can see the defect this
    // test exists to catch, which is exactly why fft_backend_parity (broadband
    // material, peak-normalized bound) passed throughout.
    //
    // If a future SDK dispatches differently and the skew stops selecting the
    // accurate kernel, this fails loudly instead of degrading in silence.
    TEST_P(fft_tonal_accuracy, EmptyBinsTrackADoubleReference) {
        const int    n  = GetParam();
        const double pi = 3.14159265358979323846;

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

        std::vector<float>              got = x;
        tap::dsp::basic_real_fft<float> ffft(static_cast<size_t>(n));
        ffft.forward_inplace(got.data());

        std::vector<double> rel;
        rel.reserve(ref.size());
        for (size_t i = 0; i < ref.size(); ++i) {
            const double d = std::fabs(static_cast<double>(got[i]) - ref[i]);
            rel.push_back(std::fabs(ref[i]) > 0.0 ? d / std::fabs(ref[i]) : 0.0);
        }
        std::sort(rel.begin(), rel.end());
        const double median = rel[rel.size() / 2];

        // Measured ~1.1e-07..1.9e-07 on every good backend (Ooura, and vDSP at
        // the skewed placement); 0.65 and worse on the 64-byte-aligned vDSP
        // kernel. 1e-5 sits ~50x above the good case and orders below the bad
        // one, so it discriminates without being brittle.
        EXPECT_LT(median, 1e-5) << "N=" << n << ": the median bin is " << median
                                << " away from the double-precision reference, i.e. most of this spectrum is wrong. "
                                << "On Apple that is the signature of the 64-byte-aligned vDSP kernel — check "
                                << "k_skew_bytes in fft.h and whether the SDK's dispatch rule has moved.";
    }

    INSTANTIATE_TEST_SUITE_P(CertifiedGeometries, fft_tonal_accuracy, ::testing::Values(512, 2048, 4096));

} // namespace
