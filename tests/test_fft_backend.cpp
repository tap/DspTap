// SPDX-License-Identifier: MIT
// Copyright 2025-2026 Timothy Place and the DspTap contributors.
//
// Parity oracle for the optional float32 FFT backends (DSPTAP_FFT_CMSIS on the
// M55, DSPTAP_FFT_ACCELERATE on Apple; see include/dsptap/fft.h and README.md).
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
#include <vector>

#include <gtest/gtest.h>

#include "dsptap/fft.h"

namespace {

    // Direct Ooura float reference, independent of the C++ wrapper's backend:
    // rdft_f is the extern "C" symbol regardless of DSPTAP_FFT_CMSIS, so it is
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

    // Forward transform of the wrapper (CMSIS when DSPTAP_FFT_CMSIS) must match
    // the Ooura reference to single-precision rounding, in the packed layout
    // and sign convention the whole DSP chain assumes.
    TEST_P(fft_backend_parity, ForwardMatchesOoura) {
        const int n = GetParam();

        const auto         x   = broadband(n, 0x9E3779B9u);
        std::vector<float> ref = x;
        ooura_ref          oref(n);
        oref.forward(ref.data());

        dsptap::basic_real_fft<float> fft(static_cast<size_t>(n));
        std::vector<float>            got = x;
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

        dsptap::basic_real_fft<float> fft(static_cast<size_t>(n));
        const auto                    x = broadband(n, 0x1234567u);

        std::vector<float> spectrum(static_cast<size_t>(n));
        std::vector<float> back(static_cast<size_t>(n));
        fft.forward(x.data(), spectrum.data());
        fft.inverse(spectrum.data(), back.data());

        for (int i = 0; i < n; ++i) {
            EXPECT_NEAR(back[static_cast<size_t>(i)], x[static_cast<size_t>(i)], 2e-5f) << "sample " << i;
        }
    }

    INSTANTIATE_TEST_SUITE_P(CertifiedGeometries, fft_backend_parity, ::testing::Values(512, 2048));

} // namespace
