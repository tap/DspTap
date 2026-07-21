// SPDX-License-Identifier: MIT
// Copyright 2025-2026 Timothy Place and the DspTap contributors.
//
// Locks down the tap::dsp::basic_real_fft contract: round-trip identity, the
// packed-spectrum layout (DC / Nyquist in slots 0 and 1), the documented
// W = exp(+2*pi*i/N) sign convention, Parseval energy conservation, and
// float/double cross-precision agreement. FDAF correctness later depends on
// every one of these staying exactly as documented in fft.h.

#include <cmath>
#include <numbers>
#include <random>
#include <vector>

#include <gtest/gtest.h>

#include "tap/dsp/fft.h"

namespace {

    template <typename Sample>
    std::vector<Sample> random_signal(size_t n, unsigned seed) {
        std::mt19937                           gen(seed);
        std::uniform_real_distribution<double> dist(-1.0, 1.0);
        std::vector<Sample>                    x(n);
        for (auto& v : x) {
            v = static_cast<Sample>(dist(gen));
        }
        return x;
    }

    template <typename Sample>
    constexpr double k_tolerance = 0.0;
    template <>
    constexpr double k_tolerance<double> = 1e-12;
    template <>
    constexpr double k_tolerance<float> = 2e-5;

    template <typename Sample>
    class real_fft_test : public ::testing::Test {};

    using sample_types = ::testing::Types<float, double>;
    TYPED_TEST_SUITE(real_fft_test, sample_types);

    TYPED_TEST(real_fft_test, SizeAndBinCount) {
        tap::dsp::basic_real_fft<TypeParam> fft(1024);
        EXPECT_EQ(fft.size(), 1024u);
        EXPECT_EQ(fft.num_bins(), 513u);
    }

    TYPED_TEST(real_fft_test, RoundTripReproducesInput) {
        constexpr size_t n = 1024;

        tap::dsp::basic_real_fft<TypeParam> fft(n);
        const auto                          x = random_signal<TypeParam>(n, 42);

        std::vector<TypeParam> spectrum(n);
        std::vector<TypeParam> back(n);
        fft.forward(x.data(), spectrum.data());
        fft.inverse(spectrum.data(), back.data());

        for (size_t i = 0; i < n; ++i) {
            EXPECT_NEAR(back[i], x[i], k_tolerance<TypeParam>) << "sample " << i;
        }
    }

    TYPED_TEST(real_fft_test, RoundTripInPlaceAndAliased) {
        constexpr size_t n = 256;

        tap::dsp::basic_real_fft<TypeParam> fft(n);
        const auto                          x = random_signal<TypeParam>(n, 7);

        // forward()/inverse() explicitly permit output aliasing input.
        auto buf = x;
        fft.forward(buf.data(), buf.data());
        fft.inverse(buf.data(), buf.data());

        for (size_t i = 0; i < n; ++i) {
            EXPECT_NEAR(buf[i], x[i], k_tolerance<TypeParam>) << "sample " << i;
        }
    }

    TYPED_TEST(real_fft_test, ImpulseHasFlatSpectrum) {
        constexpr size_t n = 64;

        tap::dsp::basic_real_fft<TypeParam> fft(n);
        std::vector<TypeParam>              x(n, TypeParam(0));
        x[0] = TypeParam(1);

        fft.forward_inplace(x.data());

        EXPECT_NEAR(x[0], 1.0, k_tolerance<TypeParam>); // DC
        EXPECT_NEAR(x[1], 1.0, k_tolerance<TypeParam>); // Nyquist
        for (size_t k = 1; k < n / 2; ++k) {
            EXPECT_NEAR(x[2 * k], 1.0, k_tolerance<TypeParam>) << "bin " << k << " real";
            EXPECT_NEAR(x[2 * k + 1], 0.0, k_tolerance<TypeParam>) << "bin " << k << " imag";
        }
    }

    TYPED_TEST(real_fft_test, DcAndNyquistPacking) {
        constexpr size_t n = 64;

        tap::dsp::basic_real_fft<TypeParam> fft(n);

        // Constant input: all energy in DC = data[0].
        std::vector<TypeParam> dc(n, TypeParam(1));
        fft.forward_inplace(dc.data());
        EXPECT_NEAR(dc[0], static_cast<double>(n), k_tolerance<TypeParam> * n);
        EXPECT_NEAR(dc[1], 0.0, k_tolerance<TypeParam> * n);

        // Alternating +1/-1: all energy in Nyquist = data[1].
        std::vector<TypeParam> nyq(n);
        for (size_t i = 0; i < n; ++i) {
            nyq[i] = (i % 2 == 0) ? TypeParam(1) : TypeParam(-1);
        }
        fft.forward_inplace(nyq.data());
        EXPECT_NEAR(nyq[0], 0.0, k_tolerance<TypeParam> * n);
        EXPECT_NEAR(nyq[1], static_cast<double>(n), k_tolerance<TypeParam> * n);
    }

    // The documented sign convention (W = exp(+2*pi*i/N)): a pure cosine at
    // bin k lands entirely in the bin's real part (+N/2), a pure sine lands
    // in the bin's IMAG part with a PLUS sign (+N/2) — the conjugate of the
    // engineering-convention DFT, where it would be -N/2.
    TYPED_TEST(real_fft_test, SignConventionIsPlusI) {
        constexpr size_t n = 128;
        constexpr size_t k = 5;

        tap::dsp::basic_real_fft<TypeParam> fft(n);
        const double w = 2.0 * std::numbers::pi * static_cast<double>(k) / static_cast<double>(n);

        std::vector<TypeParam> cosine(n);
        std::vector<TypeParam> sine(n);
        for (size_t j = 0; j < n; ++j) {
            cosine[j] = static_cast<TypeParam>(std::cos(w * static_cast<double>(j)));
            sine[j]   = static_cast<TypeParam>(std::sin(w * static_cast<double>(j)));
        }
        fft.forward_inplace(cosine.data());
        fft.forward_inplace(sine.data());

        const double half = static_cast<double>(n) / 2.0;
        const double tol  = k_tolerance<TypeParam> * static_cast<double>(n);
        EXPECT_NEAR(cosine[2 * k], half, tol);
        EXPECT_NEAR(cosine[2 * k + 1], 0.0, tol);
        EXPECT_NEAR(sine[2 * k], 0.0, tol);
        EXPECT_NEAR(sine[2 * k + 1], half, tol); // +N/2, not -N/2

        // And nothing leaks into any other bin.
        for (size_t bin = 1; bin < n / 2; ++bin) {
            if (bin == k) {
                continue;
            }
            EXPECT_NEAR(cosine[2 * bin], 0.0, tol) << "cos leak, bin " << bin;
            EXPECT_NEAR(sine[2 * bin], 0.0, tol) << "sin leak, bin " << bin;
        }
    }

    TYPED_TEST(real_fft_test, ParsevalEnergyConservation) {
        constexpr size_t n = 512;

        tap::dsp::basic_real_fft<TypeParam> fft(n);
        const auto                          x = random_signal<TypeParam>(n, 1234);

        double time_energy = 0.0;
        for (size_t i = 0; i < n; ++i) {
            time_energy += static_cast<double>(x[i]) * static_cast<double>(x[i]);
        }

        auto spectrum = x;
        fft.forward_inplace(spectrum.data());
        double freq_energy = static_cast<double>(spectrum[0]) * static_cast<double>(spectrum[0])
                             + static_cast<double>(spectrum[1]) * static_cast<double>(spectrum[1]);
        for (size_t k = 1; k < n / 2; ++k) {
            const double re = static_cast<double>(spectrum[2 * k]);
            const double im = static_cast<double>(spectrum[2 * k + 1]);
            freq_energy += 2.0 * (re * re + im * im);
        }
        freq_energy /= static_cast<double>(n);

        EXPECT_NEAR(freq_energy, time_energy, k_tolerance<TypeParam> * time_energy * static_cast<double>(n));
    }

    // The two precisions are the same algorithm on the same data: float must
    // track double to single-precision rounding, which is what makes the
    // float64 desktop oracle meaningful for the float32 embedded targets.
    TEST(RealFftCrossPrecision, FloatTracksDouble) {
        constexpr size_t n = 1024;

        const auto         xd = random_signal<double>(n, 99);
        std::vector<float> xf(n);
        for (size_t i = 0; i < n; ++i) {
            xf[i] = static_cast<float>(xd[i]);
        }

        tap::dsp::real_fft   fft64(n);
        tap::dsp::real_fft32 fft32(n);
        auto                 sd = xd;
        auto                 sf = xf;
        fft64.forward_inplace(sd.data());
        fft32.forward_inplace(sf.data());

        double err = 0.0;
        double ref = 0.0;
        for (size_t i = 0; i < n; ++i) {
            const double d = static_cast<double>(sf[i]) - sd[i];
            err += d * d;
            ref += sd[i] * sd[i];
        }
        EXPECT_LT(std::sqrt(err / ref), 1e-6);
    }

    // The double-engine float-I/O convenience overloads (used by AmbiTap's HRTF
    // analysis): a float-buffer round trip through the double FFT reproduces the
    // input to float precision.
    TEST(RealFftFloatIO, RoundTripReproducesInput) {
        constexpr size_t n = 256;

        tap::dsp::real_fft fft(n); // double engine
        const auto         xd = random_signal<double>(n, 555);
        std::vector<float> xf(n);
        for (size_t i = 0; i < n; ++i) {
            xf[i] = static_cast<float>(xd[i]);
        }

        std::vector<float> spectrum(n);
        std::vector<float> back(n);
        fft.forward(xf.data(), spectrum.data());
        fft.inverse(spectrum.data(), back.data());

        for (size_t i = 0; i < n; ++i) {
            EXPECT_NEAR(back[i], xf[i], 2e-5f) << "sample " << i;
        }
    }

    // The float-I/O forward is exactly the double forward with a float cast at
    // the end — i.e. it really runs on the double engine, not a float one.
    TEST(RealFftFloatIO, RunsOnTheDoubleEngine) {
        constexpr size_t n = 256;

        tap::dsp::real_fft fft(n);
        const auto         xd = random_signal<double>(n, 777);
        std::vector<float> xf(n);
        for (size_t i = 0; i < n; ++i) {
            xf[i] = static_cast<float>(xd[i]);
        }

        std::vector<double> spec_d(n);
        for (size_t i = 0; i < n; ++i) {
            spec_d[i] = static_cast<double>(xf[i]);
        }
        fft.forward_inplace(spec_d.data());

        std::vector<float> spec_f(n);
        fft.forward(xf.data(), spec_f.data());

        for (size_t i = 0; i < n; ++i) {
            EXPECT_FLOAT_EQ(spec_f[i], static_cast<float>(spec_d[i])) << "bin " << i;
        }
    }

} // namespace
