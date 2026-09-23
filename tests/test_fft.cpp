// SPDX-License-Identifier: MIT
// Copyright 2025-2026 Timothy Place and the DspTap contributors.
//
// Locks down the tap::dsp::basic_real_fft contract: round-trip identity, the
// packed-spectrum layout (DC / Nyquist in slots 0 and 1), the documented
// W = exp(+2*pi*i/N) sign convention, Parseval energy conservation, and
// cross-precision agreement against the double golden model. FDAF
// correctness later depends on every one of these staying exactly as
// documented in fft.h.
//
// The typed suite runs on all four profiles (Part 9 of
// docs/audit-fft-and-code-smells.md): float and double, and the Q15 / Q31
// fixed-point profiles (Stage 3b) through the per-profile trait `profile`
// below, which supplies the input conversion, the forward scale read from
// the returned exponent, and the tolerance. The float and double checks are
// numerically what they were before the widening — same signals, same
// tolerances, same comparisons; the trait is the identity for them.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <numbers>
#include <random>
#include <type_traits>
#include <vector>

#include <gtest/gtest.h>

#include "support/signals.h"
#include "tap/dsp/fft.h"
#include "tap/dsp/fft/split_radix.h"

namespace {

    using tap::dsp::test::sample_scale;

    template <typename Sample>
    std::vector<Sample> random_signal(size_t n, unsigned seed) {
        std::mt19937                           gen(seed);
        std::uniform_real_distribution<double> dist(-1.0, 1.0);
        std::vector<Sample>                    x(n);
        for (auto& v : x) {
            v = sample_scale<Sample>::from_double(dist(gen));
        }
        return x;
    }

    // Absolute tolerance in the sample's own units (fractions of full scale
    // for the fixed profiles). Float and double: the golden battery's
    // numbers, unchanged. Q15 / Q31: the largest |got - expected| over the
    // closed-form tests below (ImpulseHasFlatSpectrum, DcAndNyquistPacking,
    // SignConventionIsPlusI) measured 2026-09-22 on x86-64 Linux, GCC
    // 13.3.0 -O3, on the Stage 3b kernel as merged in tap/DspTap#27: 0.5 LSB
    // for Q15 (the full-scale Nyquist alternation's exact 32767.5 LSB,
    // clamped at the rail: the battery's pinned rail case), 1.0 LSB for Q31
    // (the on-bin sine's rounding); pinned at 2 LSB each. See profile<> for
    // the exponent handling.
    template <typename Sample>
    constexpr double k_tolerance = 0.0;
    template <>
    constexpr double k_tolerance<double> = 1e-12;
    template <>
    constexpr double k_tolerance<float> = 2e-5;
    template <>
    constexpr double k_tolerance<std::int16_t> = 2.0 / 32768.0; // 2 LSB of Q0.15
    template <>
    constexpr double k_tolerance<std::int32_t> = 2.0 / 2147483648.0; // 2 LSB of Q0.31

    // ------------------------------------------------------------------------
    // Per-profile trait. A forward transform's result is read as
    // to_double(data[i]) * 2^e, with e the exponent the fixed-point profiles
    // return (0 for float and double, whose transforms return void); the
    // inverse used here is the UNNORMALISED one in every profile, so a round
    // trip is x == back * 2^(e_fwd + e_inv + 1 - log2 n) in all four (for the
    // floating profiles that is the 2/n their inverse() applies).
    // ------------------------------------------------------------------------
    template <typename Sample>
    struct profile {
        using fft = tap::dsp::basic_real_fft<Sample>;

        static constexpr bool k_fixed_point = std::is_integral_v<Sample>;

        static double to_double(Sample v) { return sample_scale<Sample>::to_double(v); }
        static Sample from_double(double v) { return sample_scale<Sample>::from_double(v); }

        static int forward(fft& f, Sample* data) {
            if constexpr (k_fixed_point) {
                return f.forward_inplace(data);
            }
            else {
                f.forward_inplace(data);
                return 0;
            }
        }
        static int inverse_unnormalised(fft& f, Sample* data) {
            if constexpr (k_fixed_point) {
                return f.inverse_inplace(data);
            }
            else {
                f.inverse_inplace(data);
                return 0;
            }
        }
        /// Tolerance on a value the golden battery compared at k_tolerance * n:
        /// a spectrum value of magnitude ~n in the floating profiles, which the
        /// fixed profiles hold scaled by 2^-e, so their tolerance does not grow.
        static double tolerance_at_n(size_t n) {
            if constexpr (k_fixed_point) {
                return k_tolerance<Sample>;
            }
            else {
                return k_tolerance<Sample> * static_cast<double>(n);
            }
        }
    };

    template <typename Sample>
    class real_fft_test : public ::testing::Test {};

    using sample_types = ::testing::Types<float, double, std::int16_t, std::int32_t>;
    TYPED_TEST_SUITE(real_fft_test, sample_types);

    /// EXPECT_NEAR that also tracks the largest |got - expected| seen, so a
    /// test can print the measured number its pin was taken from (in the
    /// sample's units; for the fixed profiles the pin is stated in LSBs).
    class near_tracker {
      public:
        void check(double got, double expected, double tol, const char* what, size_t i) {
            m_max = std::max(m_max, std::fabs(got - expected));
            EXPECT_NEAR(got, expected, tol) << what << " " << i;
        }
        double max() const { return m_max; }

        template <typename Sample>
        void report(const char* test) const {
            if constexpr (profile<Sample>::k_fixed_point) {
                std::printf("[ measured ] %s %s: max |got - expected| = %.4f LSB (tolerance %.4f LSB)\n",
                            sizeof(Sample) == 2 ? "Q15" : "Q31", test, m_max / sample_scale<Sample>::k_lsb,
                            k_tolerance<Sample> / sample_scale<Sample>::k_lsb);
            }
        }

      private:
        double m_max = 0.0;
    };

    constexpr int log2_of(size_t n) {
        int l = 0;
        while ((size_t{1} << l) < n) {
            ++l;
        }
        return l;
    }

    /// forward() then the unnormalised inverse, reconstructed to the input's
    /// scale: back[i] * 2^(e_fwd + e_inv + 1 - log2 n), as doubles.
    template <typename Sample>
    std::vector<double> round_trip(size_t n, const std::vector<Sample>& x, bool in_place_aliased) {
        using p = profile<Sample>;
        typename p::fft     fft(n);
        std::vector<Sample> buf   = x;
        int                 e_fwd = 0;
        int                 e_inv = 0;
        if (in_place_aliased) {
            // forward()/inverse() explicitly permit output aliasing input.
            if constexpr (p::k_fixed_point) {
                e_fwd = fft.forward(buf.data(), buf.data());
                e_inv = fft.inverse(buf.data(), buf.data());
            }
            else {
                fft.forward(buf.data(), buf.data());
                fft.inverse(buf.data(), buf.data()); // the 2/n-normalised inverse, as before
            }
        }
        else {
            std::vector<Sample> spectrum(n);
            if constexpr (p::k_fixed_point) {
                e_fwd = fft.forward(x.data(), spectrum.data());
                e_inv = fft.inverse(spectrum.data(), buf.data());
            }
            else {
                fft.forward(x.data(), spectrum.data());
                std::vector<Sample> back(n);
                fft.inverse(spectrum.data(), back.data()); // the 2/n-normalised inverse, as before
                buf = back;
            }
        }
        std::vector<double> out(n);
        const double        scale = p::k_fixed_point ? std::ldexp(1.0, e_fwd + e_inv + 1 - log2_of(n)) : 1.0;
        for (size_t i = 0; i < n; ++i) {
            out[i] = p::to_double(buf[i]) * scale;
        }
        return out;
    }

    /// The round-trip tolerance: the profile's k_tolerance for float and
    /// double (as before); for the fixed profiles the reconstruction unit is
    /// the output LSB times the same power of two (under fixed scaling the
    /// trip discards log2 n bits twice: Q15 at n = 1024 reconstructs in
    /// steps of 2^-4, the honest number fft.h states), times the pinned
    /// number of those units. Measured 2026-09-18 (x86-64 Linux, GCC 13.3.0
    /// -O3, the Stage 3b kernel as merged in tap/DspTap#27) over
    /// RoundTripReproducesInput (n = 1024) and
    /// RoundTripInPlaceAndAliased (n = 256): Q15 0.5054 / 0.5098 (the
    /// output narrowing's half LSB), Q31 3.342 / 2.199 (the kernel's rounding
    /// noise over two transforms); pinned at 2x the larger.
    template <typename Sample>
    constexpr double k_round_trip_units = 0.0;
    template <>
    constexpr double k_round_trip_units<std::int16_t> = 1.02;
    template <>
    constexpr double k_round_trip_units<std::int32_t> = 6.7;

    template <typename Sample>
    double round_trip_tolerance(size_t n) {
        if constexpr (profile<Sample>::k_fixed_point) {
            const int e = profile<Sample>::fft::fixed_scaling_exponent(n);
            return k_round_trip_units<Sample> * sample_scale<Sample>::k_lsb * std::ldexp(1.0, 2 * e + 1 - log2_of(n));
        }
        else {
            return k_tolerance<Sample>;
        }
    }

    TYPED_TEST(real_fft_test, SizeAndBinCount) {
        tap::dsp::basic_real_fft<TypeParam> fft(1024);
        EXPECT_EQ(fft.size(), 1024u);
        EXPECT_EQ(fft.num_bins(), 513u);
    }

    /// The measured round-trip error in reconstructed LSBs (fixed profiles).
    template <typename Sample>
    void report_round_trip(const char* test, size_t n, double max_error) {
        if constexpr (profile<Sample>::k_fixed_point) {
            const int    e    = profile<Sample>::fft::fixed_scaling_exponent(n);
            const double unit = sample_scale<Sample>::k_lsb * std::ldexp(1.0, 2 * e + 1 - log2_of(n));
            // %lu, not %zu: newlib's printf on the QEMU legs has no %zu.
            std::printf("[ measured ] %s %s n=%lu: max error %.4f reconstructed LSB (unit %.3g; pin %.3f)\n",
                        sizeof(Sample) == 2 ? "Q15" : "Q31", test, static_cast<unsigned long>(n), max_error / unit,
                        unit, k_round_trip_units<Sample>);
        }
    }

    TYPED_TEST(real_fft_test, RoundTripReproducesInput) {
        constexpr size_t n = 1024;

        const auto x    = random_signal<TypeParam>(n, 42);
        const auto back = round_trip<TypeParam>(n, x, false);
        const auto tol  = round_trip_tolerance<TypeParam>(n);

        near_tracker t;
        for (size_t i = 0; i < n; ++i) {
            t.check(back[i], profile<TypeParam>::to_double(x[i]), tol, "sample", i);
        }
        report_round_trip<TypeParam>("RoundTripReproducesInput", n, t.max());
        EXPECT_GT(tol, 0.0) << "unmeasured pin";
    }

    TYPED_TEST(real_fft_test, RoundTripInPlaceAndAliased) {
        constexpr size_t n = 256;

        const auto x    = random_signal<TypeParam>(n, 7);
        const auto back = round_trip<TypeParam>(n, x, true);
        const auto tol  = round_trip_tolerance<TypeParam>(n);

        near_tracker t;
        for (size_t i = 0; i < n; ++i) {
            t.check(back[i], profile<TypeParam>::to_double(x[i]), tol, "sample", i);
        }
        report_round_trip<TypeParam>("RoundTripInPlaceAndAliased", n, t.max());
        EXPECT_GT(tol, 0.0) << "unmeasured pin";
    }

    TYPED_TEST(real_fft_test, ImpulseHasFlatSpectrum) {
        using p            = profile<TypeParam>;
        constexpr size_t n = 64;

        typename p::fft        fft(n);
        std::vector<TypeParam> x(n, TypeParam(0));
        x[0] = p::from_double(1.0);

        const int    e     = p::forward(fft, x.data());
        const double scale = std::ldexp(1.0, -e);
        const double tol   = k_tolerance<TypeParam>;

        near_tracker t;
        t.check(p::to_double(x[0]), 1.0 * scale, tol, "dc", 0);
        t.check(p::to_double(x[1]), 1.0 * scale, tol, "nyquist", 1);
        for (size_t k = 1; k < n / 2; ++k) {
            t.check(p::to_double(x[2 * k]), 1.0 * scale, tol, "bin real", k);
            t.check(p::to_double(x[2 * k + 1]), 0.0, tol, "bin imag", k);
        }
        t.report<TypeParam>("ImpulseHasFlatSpectrum");
        EXPECT_GT(tol, 0.0) << "unmeasured pin";
    }

    TYPED_TEST(real_fft_test, DcAndNyquistPacking) {
        using p            = profile<TypeParam>;
        constexpr size_t n = 64;

        typename p::fft fft(n);
        const double    tol = p::tolerance_at_n(n);
        near_tracker    t;

        // Constant input: all energy in DC = data[0]. The expectation is built
        // from the constant as quantised (1 - 2^-15 for Q15, 1 - 2^-31 for
        // Q31), so the fixed profiles' tolerance is spent on the kernel, not
        // on the test's own input rounding.
        const TypeParam        one = p::from_double(1.0);
        std::vector<TypeParam> dc(n, one);
        const int              e_dc = p::forward(fft, dc.data());
        t.check(p::to_double(dc[0]), static_cast<double>(n) * p::to_double(one) * std::ldexp(1.0, -e_dc), tol,
                "dc slot", 0);
        t.check(p::to_double(dc[1]), 0.0, tol, "dc slot", 1);

        // Alternating +1/-1: all energy in Nyquist = data[1]; X_{N/2} is
        // (N/2) (x_even - x_odd) on the quantised values.
        const TypeParam        minus_one = p::from_double(-1.0);
        std::vector<TypeParam> nyq(n);
        for (size_t i = 0; i < n; ++i) {
            nyq[i] = (i % 2 == 0) ? one : minus_one;
        }
        const int e_nyq = p::forward(fft, nyq.data());
        t.check(p::to_double(nyq[0]), 0.0, tol, "nyquist slot", 0);
        t.check(p::to_double(nyq[1]),
                static_cast<double>(n / 2) * (p::to_double(one) - p::to_double(minus_one)) * std::ldexp(1.0, -e_nyq),
                tol, "nyquist slot", 1);
        t.report<TypeParam>("DcAndNyquistPacking");
        EXPECT_GT(tol, 0.0) << "unmeasured pin";
    }

    // The documented sign convention (W = exp(+2*pi*i/N)): a pure cosine at
    // bin k lands entirely in the bin's real part (+N/2), a pure sine lands
    // in the bin's IMAG part with a PLUS sign (+N/2) — the conjugate of the
    // engineering-convention DFT, where it would be -N/2.
    TYPED_TEST(real_fft_test, SignConventionIsPlusI) {
        using p            = profile<TypeParam>;
        constexpr size_t n = 128;
        constexpr size_t k = 5;

        typename p::fft fft(n);
        const double    w = 2.0 * std::numbers::pi * static_cast<double>(k) / static_cast<double>(n);

        std::vector<TypeParam> cosine(n);
        std::vector<TypeParam> sine(n);
        for (size_t j = 0; j < n; ++j) {
            cosine[j] = p::from_double(std::cos(w * static_cast<double>(j)));
            sine[j]   = p::from_double(std::sin(w * static_cast<double>(j)));
        }
        const int e_cos = p::forward(fft, cosine.data());
        const int e_sin = p::forward(fft, sine.data());

        const double half_cos = static_cast<double>(n) / 2.0 * std::ldexp(1.0, -e_cos);
        const double half_sin = static_cast<double>(n) / 2.0 * std::ldexp(1.0, -e_sin);
        const double tol      = p::tolerance_at_n(n);
        near_tracker t;
        t.check(p::to_double(cosine[2 * k]), half_cos, tol, "cos re", k);
        t.check(p::to_double(cosine[2 * k + 1]), 0.0, tol, "cos im", k);
        t.check(p::to_double(sine[2 * k]), 0.0, tol, "sin re", k);
        t.check(p::to_double(sine[2 * k + 1]), half_sin, tol, "sin im (+N/2, not -N/2)", k);

        // And nothing leaks into any other bin.
        for (size_t bin = 1; bin < n / 2; ++bin) {
            if (bin == k) {
                continue;
            }
            t.check(p::to_double(cosine[2 * bin]), 0.0, tol, "cos leak, bin", bin);
            t.check(p::to_double(sine[2 * bin]), 0.0, tol, "sin leak, bin", bin);
        }
        t.report<TypeParam>("SignConventionIsPlusI");
        EXPECT_GT(tol, 0.0) << "unmeasured pin";
    }

    // Parseval's relative tolerance for the fixed profiles: the rounding
    // noise adds power, so the relative error is of the order of the per-bin
    // noise-to-signal ratio. Measured 2026-09-18 (x86-64 Linux, GCC 13.3.0
    // -O3, the Stage 3b kernel as merged in tap/DspTap#27) at n = 512 on
    // full-scale uniform noise: Q15
    // 7.15e-6, Q31 1.83e-9; pinned at 2x.
    template <typename Sample>
    constexpr double k_parseval_relative = 0.0;
    template <>
    constexpr double k_parseval_relative<std::int16_t> = 1.43e-5;
    template <>
    constexpr double k_parseval_relative<std::int32_t> = 3.7e-9;

    TYPED_TEST(real_fft_test, ParsevalEnergyConservation) {
        using p            = profile<TypeParam>;
        constexpr size_t n = 512;

        typename p::fft fft(n);
        const auto      x = random_signal<TypeParam>(n, 1234);

        double time_energy = 0.0;
        for (size_t i = 0; i < n; ++i) {
            time_energy += p::to_double(x[i]) * p::to_double(x[i]);
        }

        auto         spectrum    = x;
        const int    e           = p::forward(fft, spectrum.data());
        const double scale       = std::ldexp(1.0, e); // spectrum values are X * 2^-e
        double       freq_energy = p::to_double(spectrum[0]) * scale * (p::to_double(spectrum[0]) * scale)
                             + p::to_double(spectrum[1]) * scale * (p::to_double(spectrum[1]) * scale);
        for (size_t k = 1; k < n / 2; ++k) {
            const double re = p::to_double(spectrum[2 * k]) * scale;
            const double im = p::to_double(spectrum[2 * k + 1]) * scale;
            freq_energy += 2.0 * (re * re + im * im);
        }
        freq_energy /= static_cast<double>(n);

        if constexpr (p::k_fixed_point) {
            std::printf("[ measured ] %s ParsevalEnergyConservation: relative error %.3e (pin %.3e)\n",
                        sizeof(TypeParam) == 2 ? "Q15" : "Q31", std::fabs(freq_energy - time_energy) / time_energy,
                        k_parseval_relative<TypeParam>);
            EXPECT_GT(k_parseval_relative<TypeParam>, 0.0) << "unmeasured pin";
            EXPECT_NEAR(freq_energy, time_energy, k_parseval_relative<TypeParam> * time_energy);
        }
        else {
            EXPECT_NEAR(freq_energy, time_energy, k_tolerance<TypeParam> * time_energy * static_cast<double>(n));
        }
    }

    // ------------------------------------------------------------------------
    // Cross-precision: every profile against the double golden model.
    // ------------------------------------------------------------------------

    /// The relative 2-norm error of a profile's forward transform against the
    /// double golden model on the same signal: sqrt(sum |s - g * 2^-e|^2 /
    /// sum |g * 2^-e|^2) over the packed spectrum, with g the double result
    /// on the profile's own quantised input.
    template <typename Sample>
    double relative_error_vs_double(size_t n, unsigned seed) {
        using p               = profile<Sample>;
        const auto          x = random_signal<Sample>(n, seed);
        std::vector<double> xd(n);
        for (size_t i = 0; i < n; ++i) {
            xd[i] = p::to_double(x[i]);
        }

        tap::dsp::real_fft fft64(n);
        typename p::fft    fftp(n);
        auto               sd = xd;
        auto               sp = x;
        fft64.forward_inplace(sd.data());
        const int    e     = p::forward(fftp, sp.data());
        const double scale = std::ldexp(1.0, -e);

        double err = 0.0;
        double ref = 0.0;
        for (size_t i = 0; i < n; ++i) {
            const double g = sd[i] * scale;
            const double d = p::to_double(sp[i]) - g;
            err += d * d;
            ref += g * g;
        }
        return std::sqrt(err / ref);
    }

    /// The relative 2-norm error of a float forward transform against the
    /// double golden model on the same signal (float -> double is exact, so
    /// both engines see identical input values).
    template <typename FloatEngine>
    double float_engine_error_vs_double(size_t n, unsigned seed) {
        const auto         xd = random_signal<double>(n, seed);
        std::vector<float> xf(n);
        for (size_t i = 0; i < n; ++i) {
            xf[i] = static_cast<float>(xd[i]);
        }
        std::vector<double> sd(n);
        for (size_t i = 0; i < n; ++i) {
            sd[i] = static_cast<double>(xf[i]);
        }
        auto sf = xf;

        tap::dsp::real_fft fft64(n);
        FloatEngine        fft32(n);
        fft64.forward_inplace(sd.data());
        fft32.forward_inplace(sf.data());

        double err = 0.0;
        double ref = 0.0;
        for (size_t i = 0; i < n; ++i) {
            const double d = static_cast<double>(sf[i]) - sd[i];
            err += d * d;
            ref += sd[i] * sd[i];
        }
        return std::sqrt(err / ref);
    }

    // The two precisions are the same algorithm on the same data: float must
    // track double to single-precision rounding, which is what makes the
    // float64 desktop oracle meaningful for the float32 embedded targets.
    // Through basic_real_fft, so on a backend build (CMSIS, vDSP) it is the
    // backend that is held to the bound.
    TEST(RealFftCrossPrecision, FloatTracksDouble) {
        EXPECT_LT(float_engine_error_vs_double<tap::dsp::real_fft32>(1024, 99), 1e-6);
    }

    // The float profile's noise floor as a MEASURED number (docs/fft-design.md,
    // contract table, "Noise floor"): the same metric at N = 512 on the
    // split-radix engine itself, which is basic_real_fft<float> wherever no
    // backend define is active and is bit-identical to the vendored C's float
    // build on every leg, so the number is the shipping float profile's on
    // every host and does not change on the M55 / macOS backend legs.
    // Measured 2026-09-23 (x86-64 Linux, GCC 13.3.0 and Clang 18.1.3 -O3,
    // glibc 2.39, the engine as routed at Stage 2b; both compilers print the
    // same value): 1.1236e-7, against the audit's Part 6 N2 probe value of
    // 1.105e-7 (different material; the probe was never a committed test).
    // Pinned at 2x: the float rounding sequence is fixed by the statements,
    // and the VFMA legs and libm differences move the last bits, not the rms.
    constexpr double k_float_engine_tracks_double_512 = 2.25e-7;

    TEST(RealFftCrossPrecision, FloatEngineTracksDoubleAtN512) {
        const double err = float_engine_error_vs_double<tap::dsp::detail::split_radix_rdft<float>>(512, 99);
        std::printf("[ measured ] FloatEngineTracksDoubleAtN512: relative 2-norm error %.4e (pin %.3e)\n", err,
                    k_float_engine_tracks_double_512);
        EXPECT_GT(k_float_engine_tracks_double_512, 0.0) << "unmeasured pin";
        EXPECT_LT(err, k_float_engine_tracks_double_512) << "measured " << err;
    }

    // The fixed-point profiles against the same golden model: the relative
    // 2-norm error of the fixed-scaling forward at n = 1024 on full-scale
    // uniform noise, a measured number pinned at 2x (the log_mel pattern).
    // Q15's error is the output narrowing (2^-15 / sqrt(12) per value on a
    // spectrum whose per-component RMS is sqrt(1/3 / 2n) = 0.0128: predicted
    // 6.9e-4); Q31's is the int32 kernel's rounding noise, orders of
    // magnitude lower. Measured 2026-09-18 (x86-64 Linux, GCC 13.3.0 -O3,
    // the Stage 3b kernel as merged in tap/DspTap#27): Q15 6.949e-4, Q31
    // 5.458e-8.
    constexpr double k_q15_tracks_double = 1.4e-3;
    constexpr double k_q31_tracks_double = 1.1e-7;

    TEST(RealFftCrossPrecision, Q15TracksDouble) {
        const double err = relative_error_vs_double<std::int16_t>(1024, 99);
        std::printf("[ measured ] Q15TracksDouble: relative 2-norm error %.3e (pin %.3e)\n", err, k_q15_tracks_double);
        EXPECT_GT(k_q15_tracks_double, 0.0) << "unmeasured pin";
        EXPECT_LT(err, k_q15_tracks_double) << "measured " << err;
    }

    TEST(RealFftCrossPrecision, Q31TracksDouble) {
        const double err = relative_error_vs_double<std::int32_t>(1024, 99);
        std::printf("[ measured ] Q31TracksDouble: relative 2-norm error %.3e (pin %.3e)\n", err, k_q31_tracks_double);
        EXPECT_GT(k_q31_tracks_double, 0.0) << "unmeasured pin";
        EXPECT_LT(err, k_q31_tracks_double) << "measured " << err;
    }

    // The double-engine float-I/O convenience overloads (used by AmbiTap's HRTF
    // analysis): a float-buffer round trip through the double FFT reproduces the
    // input to float precision. The overloads are [[deprecated]] since Stage
    // 2b (Decision D5: removed after one consumer cycle); these two tests keep
    // pinning them until they go, so the warning is silenced here and nowhere
    // else.
#if defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#elif defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable : 4996)
#endif
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
#if defined(__GNUC__)
#pragma GCC diagnostic pop
#elif defined(_MSC_VER)
#pragma warning(pop)
#endif

} // namespace
