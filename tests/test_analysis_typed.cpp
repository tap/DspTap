// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// The typed analysis instruments (Stage 3c of docs/audit-fft-and-code-smells.md):
// sine_analysis.h and multitone_analysis.h over any contiguous range of float /
// double / Q15 / Q31 samples. Its own translation unit and CMake target, like
// the Ooura parity gate of Stage 2a, so that it can carry -ffp-contract=off:
// the pins below are same-binary bit-identity comparisons, and GCC's default
// cross-statement contraction (-ffp-contract=fast) fuses multiply-adds by
// inlining and unrolling context, so two textually identical functions can
// differ by an ulp on a host with FMA hardware (measured: g++ 13 -O3 -mfma on
// the sub-span and tracked fits; clang, which contracts per expression, and
// MSVC, which does not contract, agree bit for bit either way). With
// contraction off both copies are the same sequence of IEEE operations on
// every host, which is what "byte for byte" means. The smoke battery of the
// instruments themselves stays in test_analysis.cpp.

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <numbers>
#include <span>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include "tap/dsp/analysis/multitone_analysis.h"
#include "tap/dsp/analysis/sine_analysis.h"
#include "tap/dsp/sample_traits.h"

namespace {

    namespace an = tap::dsp::analysis;

    // ------------------------------------------------------------------------
    // Stage 3c (docs/audit-fft-and-code-smells.md): the instruments are typed
    // over any contiguous range of float / double / Q15 / Q31 samples. Two
    // promises are pinned below. (1) The floating path is byte for byte what
    // it was before the templates: `legacy` holds the pre-Stage-3c entry
    // points verbatim (the std::span<const float> fit_sine / fit_sine_tracked
    // of sine_analysis.h and the program_weighted_snr_db prologue of
    // multitone_analysis.h, whose joint-fit internals are unchanged and are
    // called as they are), and the tests compare the new instantiations
    // against them bit for bit in the same binary. A same-binary A/B is the
    // valid form of this pin: a golden hash would fail across the CI hosts'
    // libms by design (Part 13), while two copies of the same expressions in
    // one TU under one set of flags compute the same bits. (2) The integer
    // instantiations read Q0.15 / Q0.31 fractions: bit-identical to the
    // double instrument on the same fractions, and within the quantisation
    // noise of the double instrument on the unquantised signal (measured,
    // pinned with margin). Two things keep (1) exact: the floating read in
    // the header is the bare cast the legacy code evaluated (the first
    // version multiplied by an exact 1.0 step, and the macOS arm64 leg moved
    // residual_rms by one ulp through a fused multiply-subtract), and this
    // TU is compiled with -ffp-contract=off (file comment).
    // ------------------------------------------------------------------------
    namespace legacy {

        an::sine_fit fit_sine(std::span<const float> x, double freq_norm) {
            const double w       = 2.0 * std::numbers::pi * freq_norm;
            double       m[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
            double       rhs[3]  = {0, 0, 0};
            for (std::size_t i = 0; i < x.size(); ++i) {
                const double s        = std::sin(w * static_cast<double>(i));
                const double c        = std::cos(w * static_cast<double>(i));
                const double basis[3] = {s, c, 1.0};
                for (int r = 0; r < 3; ++r) {
                    for (int q = 0; q < 3; ++q) {
                        m[r][q] += basis[r] * basis[q];
                    }
                    rhs[r] += basis[r] * static_cast<double>(x[i]);
                }
            }
            int order[3] = {0, 1, 2};
            for (int col = 0; col < 3; ++col) {
                int piv = col;
                for (int r = col + 1; r < 3; ++r) {
                    if (std::abs(m[order[r]][col]) > std::abs(m[order[piv]][col])) {
                        piv = r;
                    }
                }
                std::swap(order[col], order[piv]);
                const int p = order[col];
                for (int r = col + 1; r < 3; ++r) {
                    const int    rr = order[r];
                    const double f  = m[rr][col] / m[p][col];
                    for (int q = col; q < 3; ++q) {
                        m[rr][q] -= f * m[p][q];
                    }
                    rhs[rr] -= f * rhs[p];
                }
            }
            double sol[3];
            for (int col = 2; col >= 0; --col) {
                const int p = order[col];
                double    v = rhs[p];
                for (int q = col + 1; q < 3; ++q) {
                    v -= m[p][q] * sol[q];
                }
                sol[col] = v / m[p][col];
            }
            an::sine_fit fit;
            fit.amplitude = std::hypot(sol[0], sol[1]);
            fit.phase     = std::atan2(sol[1], sol[0]);
            fit.dc        = sol[2];
            double sq     = 0.0;
            for (std::size_t i = 0; i < x.size(); ++i) {
                const double s = std::sin(w * static_cast<double>(i));
                const double c = std::cos(w * static_cast<double>(i));
                const double r = static_cast<double>(x[i]) - (sol[0] * s + sol[1] * c + sol[2]);
                sq += r * r;
            }
            fit.residual_rms = std::sqrt(sq / static_cast<double>(x.size()));
            fit.freq_norm    = freq_norm;
            return fit;
        }

        an::sine_fit fit_sine_tracked(std::span<const float> x, double freq_norm_guess) {
            double            f    = freq_norm_guess;
            const std::size_t half = x.size() / 2;
            for (int iter = 0; iter < 4; ++iter) {
                const an::sine_fit a         = fit_sine(x.first(half), f);
                const an::sine_fit b         = fit_sine(x.subspan(half), f);
                const double       two_pi    = 2.0 * std::numbers::pi;
                const double       predicted = a.phase + two_pi * f * static_cast<double>(half);
                const double       dphi      = std::remainder(b.phase - predicted, two_pi);
                f += dphi / (two_pi * static_cast<double>(half));
            }
            return fit_sine(x, f);
        }

        double program_weighted_snr_db(std::span<const float> tail, const an::tone_comb& comb, double /*fsIn*/,
                                       double fs_out) {
            std::vector<double> work(tail.begin(), tail.end());
            const std::size_t   k = comb.freq_hz.size();
            std::vector<double> nus(k);
            for (std::size_t t = 0; t < k; ++t) {
                nus[t] = comb.freq_hz[t] / fs_out;
            }
            std::vector<an::tone_fit> fits(k);
            an::joint_fit_residual_power(work, nus, fits);
            std::vector<double> lone(work.size());
            double              resid_power = 0.0;
            for (int round = 0; round < 2; ++round) {
                std::vector<double> resid(work);
                for (std::size_t i = 0; i < work.size(); ++i) {
                    double model = 0.0;
                    for (std::size_t t = 0; t < k; ++t) {
                        const double w = 2.0 * std::numbers::pi * nus[t] * static_cast<double>(i);
                        model += fits[t].a * std::sin(w) + fits[t].b * std::cos(w);
                    }
                    resid[i] -= model;
                }
                double rho_num = 0.0, rho_den = 0.0;
                for (std::size_t t = 0; t < k; ++t) {
                    const double w = 2.0 * std::numbers::pi * nus[t];
                    for (std::size_t i = 0; i < lone.size(); ++i) {
                        lone[i] = resid[i] + fits[t].a * std::sin(w * static_cast<double>(i))
                                  + fits[t].b * std::cos(w * static_cast<double>(i));
                    }
                    const double refined = an::track_tone_freq(lone, nus[t]);
                    const double wt      = comb.amplitude[t] * nus[t];
                    rho_num += wt * wt * (refined / nus[t]);
                    rho_den += wt * wt;
                }
                const double rho = rho_num / rho_den;
                for (std::size_t t = 0; t < k; ++t) {
                    nus[t] *= rho;
                }
                resid_power = an::joint_fit_residual_power(work, nus, fits);
            }
            double signal = 0.0;
            for (const auto& f : fits) {
                signal += f.power();
            }
            return 10.0 * std::log10(signal / resid_power);
        }

    } // namespace legacy

    /// Bitwise equality of two doubles (NaN-safe, -0.0 != +0.0): the A/B pin.
    bool same_bits(double a, double b) {
        std::uint64_t ua = 0;
        std::uint64_t ub = 0;
        std::memcpy(&ua, &a, sizeof a);
        std::memcpy(&ub, &b, sizeof b);
        return ua == ub;
    }

    ::testing::AssertionResult fits_are_bit_identical(const an::sine_fit& a, const an::sine_fit& b) {
        const std::array<std::pair<const char*, std::pair<double, double>>, 5> fields{{
            {"amplitude", {a.amplitude, b.amplitude}},
            {"phase", {a.phase, b.phase}},
            {"dc", {a.dc, b.dc}},
            {"residual_rms", {a.residual_rms, b.residual_rms}},
            {"freq_norm", {a.freq_norm, b.freq_norm}},
        }};
        for (const auto& [name, v] : fields) {
            if (!same_bits(v.first, v.second)) {
                return ::testing::AssertionFailure()
                       << name << " differs: " << v.first << " vs " << v.second << " (bit patterns differ)";
            }
        }
        return ::testing::AssertionSuccess();
    }

    /// The double-domain signal quantised to a fixed-point sample type by the
    /// substrate's round_sat (half away from zero, saturating), the rule the
    /// capi boundary and tests/support/signals.h use.
    template <typename I>
    std::vector<I> quantise(const std::vector<double>& x) {
        constexpr double scale = static_cast<double>(std::int64_t{1} << tap::dsp::sample_traits<I>::k_sample_frac_bits);
        std::vector<I>   q(x.size());
        for (std::size_t i = 0; i < x.size(); ++i) {
            q[i] = tap::dsp::detail::round_sat<I>(x[i] * scale);
        }
        return q;
    }

    template <typename I>
    std::vector<double> fractions(const std::vector<I>& q) {
        constexpr double scale = static_cast<double>(std::int64_t{1} << tap::dsp::sample_traits<I>::k_sample_frac_bits);
        std::vector<double> d(q.size());
        for (std::size_t i = 0; i < q.size(); ++i) {
            d[i] = static_cast<double>(q[i]) / scale;
        }
        return d;
    }

    // ------------------------------------------------------------------------
    // Stage 3c: the typed instruments.
    // ------------------------------------------------------------------------

    // The float and double instantiations compute exactly what the
    // pre-template std::span<const float> instrument computed: the same
    // expressions over the same values (a floating sample times the exact
    // 1.0 step is the sample), so the bits agree, for a vector, a span and a
    // span subrange, at exponent 0. Two signals: an exact tone (the existing
    // smoke case) and a tone with noise and DC, so every solver branch and
    // the residual path see non-trivial values.
    TEST(SineAnalysisTyped, FloatingSpansAreBitIdenticalToThePreTemplateInstrument) {
        const double       nu = 997.0 / 48000.0;
        std::vector<float> clean(16384);
        std::vector<float> dirty(16384);
        std::uint32_t      state = 0x2545F491u;
        for (std::size_t i = 0; i < clean.size(); ++i) {
            const double t = 2.0 * std::numbers::pi * nu * static_cast<double>(i) + 0.3;
            state ^= state << 13;
            state ^= state >> 17;
            state ^= state << 5;
            const double u = static_cast<double>(state) / 2147483648.0 - 1.0;
            clean[i]       = static_cast<float>(0.5 * std::sin(t));
            dirty[i]       = static_cast<float>(0.5 * std::sin(t) + 0.01 * u + 0.02);
        }
        for (const auto* x : {&clean, &dirty}) {
            const std::span<const float> xs(*x);
            EXPECT_TRUE(fits_are_bit_identical(an::fit_sine(*x, nu), legacy::fit_sine(xs, nu)));
            EXPECT_TRUE(fits_are_bit_identical(an::fit_sine(xs, nu), legacy::fit_sine(xs, nu)));
            EXPECT_TRUE(fits_are_bit_identical(an::fit_sine(xs, nu, 0), legacy::fit_sine(xs, nu)));
            EXPECT_TRUE(fits_are_bit_identical(an::fit_sine(xs.subspan(1000, 5000), nu),
                                               legacy::fit_sine(xs.subspan(1000, 5000), nu)));
            EXPECT_TRUE(fits_are_bit_identical(an::fit_sine_tracked(*x, nu * (1.0 + 5e-6)),
                                               legacy::fit_sine_tracked(xs, nu * (1.0 + 5e-6))));
            // The double instantiation on the widened floats: the same values, the same bits.
            const std::vector<double> xd(x->begin(), x->end());
            EXPECT_TRUE(fits_are_bit_identical(an::fit_sine(xd, nu), legacy::fit_sine(xs, nu)));
            EXPECT_TRUE(fits_are_bit_identical(an::fit_sine_tracked(xd, nu * (1.0 + 5e-6)),
                                               legacy::fit_sine_tracked(xs, nu * (1.0 + 5e-6))));
        }
    }

    // A Q15 / Q31 span is read as Q0.15 / Q0.31 fractions: the integer
    // instantiation is bit-identical to the double instrument on the same
    // fractions (the read is exact), and agrees with the double instrument on
    // the UNQUANTISED signal to within the quantisation noise: the residual
    // is the quantiser's own LSB / sqrt(12) (0.289 LSB rms) and the amplitude
    // and DC move by a few thousandths of an LSB (the residual / sqrt(N/2)
    // it implies, 3.2e-3 LSB). Measured 2026-09-23 on x86-64 Linux (Ubuntu
    // 24.04, glibc 2.39, GCC 13.3.0 -O3) for a 0.5 full-scale tone at 997 Hz
    // / 48 kHz over 16384 samples (the exact-signal double fit scores 266 dB):
    //   Q15: residual 0.2886 LSB rms, SNR 92.07 dB (theory 1.76 + 6.02 * 15
    //        - 6.02 = 92.06), |amplitude - exact| 1.69e-3 LSB, |dc| 2.60e-3 LSB
    //   Q31: residual 0.2873 LSB rms, SNR 188.44 dB (theory 188.39),
    //        |amplitude - exact| 2.04e-4 LSB, |dc| 2.58e-3 LSB
    // Pinned at 2x on the residual band and the amplitude / DC deviations
    // and at -0.5 dB on the SNR; every value is printed for a -V run to record.
    template <typename I>
    void quantised_sine_agrees(double residual_lo_lsb, double residual_hi_lsb, double snr_min_db,
                               double amplitude_max_lsb, double dc_max_lsb) {
        constexpr double lsb =
            1.0 / static_cast<double>(std::int64_t{1} << tap::dsp::sample_traits<I>::k_sample_frac_bits);
        const double        nu = 997.0 / 48000.0;
        std::vector<double> exact(16384);
        for (std::size_t i = 0; i < exact.size(); ++i) {
            exact[i] = 0.5 * std::sin(2.0 * std::numbers::pi * nu * static_cast<double>(i) + 0.3);
        }
        const std::vector<I> q = quantise<I>(exact);

        const an::sine_fit fit_q = an::fit_sine(q, nu);
        const an::sine_fit fit_f = an::fit_sine(fractions(q), nu);
        const an::sine_fit fit_x = an::fit_sine(exact, nu);
        EXPECT_TRUE(fits_are_bit_identical(fit_q, fit_f)) << "the integer read is not the exact fraction";
        EXPECT_TRUE(fits_are_bit_identical(an::fit_sine(std::span<const I>(q), nu), fit_f));
        EXPECT_TRUE(fits_are_bit_identical(an::fit_sine_tracked(q, nu * (1.0 + 5e-6)),
                                           an::fit_sine_tracked(fractions(q), nu * (1.0 + 5e-6))));

        const double residual_lsb = fit_q.residual_rms / lsb;
        const double amp_dev_lsb  = std::abs(fit_q.amplitude - fit_x.amplitude) / lsb;
        const double dc_dev_lsb   = std::abs(fit_q.dc - fit_x.dc) / lsb;
        std::printf("[ measured ] %s quantised sine: residual %.4f LSB rms (1/sqrt(12) = %.4f), snr %.2f dB, "
                    "|amplitude - exact| %.3e LSB, |dc - exact| %.3e LSB, exact-signal fit snr %.2f dB\n",
                    sizeof(I) == 2 ? "Q15" : "Q31", residual_lsb, 1.0 / std::sqrt(12.0), an::snr_db(fit_q), amp_dev_lsb,
                    dc_dev_lsb, an::snr_db(fit_x));
        EXPECT_GT(residual_lsb, residual_lo_lsb);
        EXPECT_LT(residual_lsb, residual_hi_lsb);
        EXPECT_GT(an::snr_db(fit_q), snr_min_db);
        EXPECT_LT(amp_dev_lsb, amplitude_max_lsb);
        EXPECT_LT(dc_dev_lsb, dc_max_lsb);

        // The scale exponent: the same block read at 2^e is the same fit at
        // 2^e times the level (snr unchanged), for e of either sign.
        for (const int e : {-3, 2}) {
            const an::sine_fit fit_e = an::fit_sine(q, nu, e);
            const double       k     = std::ldexp(1.0, e);
            EXPECT_TRUE(same_bits(fit_e.amplitude, fit_q.amplitude * k));
            EXPECT_TRUE(same_bits(fit_e.residual_rms, fit_q.residual_rms * k));
            EXPECT_TRUE(same_bits(fit_e.dc, fit_q.dc * k));
            EXPECT_TRUE(same_bits(fit_e.phase, fit_q.phase));
        }
    }

    TEST(SineAnalysisTyped, Q15SpanAgreesWithTheDoubleInstrumentToTheQuantisationNoise) {
        quantised_sine_agrees<std::int16_t>(0.2886 / 2.0, 0.2886 * 2.0, 92.07 - 0.5, 1.69e-3 * 2.0, 2.60e-3 * 2.0);
    }

    TEST(SineAnalysisTyped, Q31SpanAgreesWithTheDoubleInstrumentToTheQuantisationNoise) {
        quantised_sine_agrees<std::int32_t>(0.2873 / 2.0, 0.2873 * 2.0, 188.44 - 0.5, 2.04e-4 * 2.0, 2.58e-3 * 2.0);
    }

    // The multitone instrument: the float instantiation is bit-identical to
    // the pre-template one, the integer instantiations equal the double
    // instrument on the same fractions bit for bit, and a Q15 / Q31 comb tail
    // scores at the format's quantisation floor against the exact comb.
    // Measured 2026-09-23 (x86-64 Linux, glibc 2.39, GCC 13.3.0 -O3), 12 pink
    // tones, 32768 samples at 48 kHz, peak sum 0.5: Q15 83.86 dB, Q31
    // 180.15 dB (the float tail 151.88 dB); pinned at -1 dB.
    TEST(MultitoneAnalysisTyped, TypedTailsMatchTheDoubleInstrumentAndReachTheirQuantisationFloor) {
        const an::tone_comb comb = an::tone_comb::pink(12, 40.0, 20000.0, 0.5);
        std::vector<double> exact(32768);
        for (std::size_t i = 0; i < exact.size(); ++i) {
            exact[i] = comb.sample_at(i, 48000.0);
        }
        std::vector<float> tail32(exact.size());
        for (std::size_t i = 0; i < exact.size(); ++i) {
            tail32[i] = static_cast<float>(exact[i]);
        }
        const double snr_new    = an::program_weighted_snr_db(tail32, comb, 48000.0, 48000.0);
        const double snr_legacy = legacy::program_weighted_snr_db(tail32, comb, 48000.0, 48000.0);
        EXPECT_TRUE(same_bits(snr_new, snr_legacy)) << snr_new << " vs " << snr_legacy;
        EXPECT_TRUE(
            same_bits(an::program_weighted_snr_db(std::span<const float>(tail32), comb, 48000.0, 48000.0), snr_legacy));

        const auto   q15     = quantise<std::int16_t>(exact);
        const auto   q31     = quantise<std::int32_t>(exact);
        const double snr_q15 = an::program_weighted_snr_db(q15, comb, 48000.0, 48000.0);
        const double snr_q31 = an::program_weighted_snr_db(q31, comb, 48000.0, 48000.0);
        EXPECT_TRUE(same_bits(snr_q15, an::program_weighted_snr_db(fractions(q15), comb, 48000.0, 48000.0)));
        EXPECT_TRUE(same_bits(snr_q31, an::program_weighted_snr_db(fractions(q31), comb, 48000.0, 48000.0)));
        std::printf(
            "[ measured ] program-weighted snr on the quantised comb: float %.2f dB, Q15 %.2f dB, Q31 %.2f dB\n",
            snr_legacy, snr_q15, snr_q31);
        EXPECT_GT(snr_q15, 83.86 - 1.0);
        EXPECT_GT(snr_q31, 180.15 - 1.0);
    }

} // namespace
