// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// Locks down the tap::dsp::basic_decimator contract: the pinned tap counts,
// sample-for-sample agreement of the float profile with the committed numpy
// reference for every ratio and profile (and of the double golden model to
// the reference's own float32 rounding), chunking invariance, exact output
// counting, unity DC gain (exact in Q15 through row-sum-preserving
// quantization), the stated passband/stopband numbers measured from the
// shipped coefficients, the integer group delay, the Q15 / Q31 profiles
// tracking the double golden model within each format's own floor, and a
// byte-wise bit pin of every fixed-point coefficient table. The float/double
// contract tests are typed over both profiles.

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstring>
#include <numbers>
#include <span>
#include <type_traits>
#include <vector>

#include <gtest/gtest.h>

#include "reference/frontend_vectors.h"
#include "tap/dsp/decimate.h"

namespace {

    using tap::dsp::basic_decimator;
    using tap::dsp::decimate_profile;

    template <typename S>
    class decimate_test : public ::testing::Test {};
    using floating_types = ::testing::Types<float, double>;
    TYPED_TEST_SUITE(decimate_test, floating_types, );

    // Float coefficient storage against float64 reference coefficients: the
    // same -90 dB floor RatioTap pins at 3e-5 on a 0.9 peak.
    constexpr float k_reference_tolerance = 3e-5f;

    template <std::size_t M>
    void expect_matches_reference(const decimate_profile& p, const float* ref, std::size_t n_ref) {
        basic_decimator<float, M> dec(p);
        const auto&               x = frontend_ref::k_dec_input;
        ASSERT_EQ(dec.outputs_for(x.size()), n_ref);
        std::vector<float> y(n_ref);
        ASSERT_EQ(dec.process(x.data(), x.size(), y.data()), n_ref);
        float worst = 0.0f;
        for (std::size_t i = 0; i < n_ref; ++i) {
            worst = std::max(worst, std::abs(y[i] - ref[i]));
        }
        EXPECT_LT(worst, k_reference_tolerance) << "ratio " << M << " taps " << dec.taps();
    }

    TEST(Decimate, TapCountsMatchTheSearchedMinima) {
        EXPECT_EQ(decimate_profile::economy().taps<2>(), frontend_ref::k_dec_taps_economy_2);
        EXPECT_EQ(decimate_profile::economy().taps<3>(), frontend_ref::k_dec_taps_economy_3);
        EXPECT_EQ(decimate_profile::economy().taps<6>(), frontend_ref::k_dec_taps_economy_6);
        EXPECT_EQ(decimate_profile::transparent().taps<2>(), frontend_ref::k_dec_taps_transparent_2);
        EXPECT_EQ(decimate_profile::transparent().taps<3>(), frontend_ref::k_dec_taps_transparent_3);
        EXPECT_EQ(decimate_profile::transparent().taps<6>(), frontend_ref::k_dec_taps_transparent_6);
    }

    TEST(Decimate, MatchesNumpyReference) {
        using namespace frontend_ref;
        expect_matches_reference<2>(decimate_profile::economy(), k_dec_economy_2.data(), k_dec_economy_2.size());
        expect_matches_reference<3>(decimate_profile::economy(), k_dec_economy_3.data(), k_dec_economy_3.size());
        expect_matches_reference<6>(decimate_profile::economy(), k_dec_economy_6.data(), k_dec_economy_6.size());
        expect_matches_reference<2>(decimate_profile::transparent(), k_dec_transparent_2.data(),
                                    k_dec_transparent_2.size());
        expect_matches_reference<3>(decimate_profile::transparent(), k_dec_transparent_3.data(),
                                    k_dec_transparent_3.size());
        expect_matches_reference<6>(decimate_profile::transparent(), k_dec_transparent_6.data(),
                                    k_dec_transparent_6.size());
    }

    // The double golden model against the numpy float64 reference, which is
    // committed rounded to float32: the disagreement is the reference's own
    // rounding. Measured 2026-09 at 5.75e-8 over every ratio and profile
    // (half a float ulp at the signal's level); pinned at 2x. This bounds
    // the double path to one float ulp of numpy, no tighter: a double
    // reference (the generator already emits k_mel_log as double) would pin
    // it near 1e-12 — a recorded follow-up for the reference-script owner,
    // since the numpy pin was frozen for this change.
    TEST(Decimate, DoubleMatchesNumpyReferenceToItsFloatRounding) {
        using namespace frontend_ref;
        const auto check = [](auto ratio, const decimate_profile& p, const float* ref, std::size_t n_ref) {
            constexpr std::size_t      M = decltype(ratio)::value;
            basic_decimator<double, M> dec(p);
            const std::vector<double>  x(k_dec_input.begin(), k_dec_input.end());
            std::vector<double>        y(n_ref);
            ASSERT_EQ(dec.process(x.data(), x.size(), y.data()), n_ref);
            double worst = 0.0;
            for (std::size_t i = 0; i < n_ref; ++i) {
                worst = std::max(worst, std::abs(y[i] - static_cast<double>(ref[i])));
            }
            EXPECT_LT(worst, 1.15e-7) << "ratio " << M << " taps " << dec.taps();
        };
        using two   = std::integral_constant<std::size_t, 2>;
        using three = std::integral_constant<std::size_t, 3>;
        using six   = std::integral_constant<std::size_t, 6>;
        check(two{}, decimate_profile::economy(), k_dec_economy_2.data(), k_dec_economy_2.size());
        check(three{}, decimate_profile::economy(), k_dec_economy_3.data(), k_dec_economy_3.size());
        check(six{}, decimate_profile::economy(), k_dec_economy_6.data(), k_dec_economy_6.size());
        check(two{}, decimate_profile::transparent(), k_dec_transparent_2.data(), k_dec_transparent_2.size());
        check(three{}, decimate_profile::transparent(), k_dec_transparent_3.data(), k_dec_transparent_3.size());
        check(six{}, decimate_profile::transparent(), k_dec_transparent_6.data(), k_dec_transparent_6.size());
    }

    TYPED_TEST(decimate_test, ChunkingIsBitIdentical) {
        using sample                 = TypeParam;
        const auto&                x = frontend_ref::k_dec_input;
        const std::vector<sample>  xs(x.begin(), x.end());
        basic_decimator<sample, 3> whole;
        std::vector<sample>        ref(whole.outputs_for(xs.size()));
        whole.process(xs.data(), xs.size(), ref.data());
        for (const std::size_t chunk :
             {std::size_t{1}, std::size_t{2}, std::size_t{7}, std::size_t{255}, std::size_t{1000}}) {
            basic_decimator<sample, 3> dec;
            std::vector<sample>        y;
            for (std::size_t i = 0; i < xs.size(); i += chunk) {
                const std::size_t   n = std::min(chunk, xs.size() - i);
                std::vector<sample> part(dec.outputs_for(n));
                ASSERT_EQ(dec.process(xs.data() + i, n, part.data()), part.size());
                y.insert(y.end(), part.begin(), part.end());
            }
            ASSERT_EQ(y.size(), ref.size()) << "chunk " << chunk;
            for (std::size_t i = 0; i < ref.size(); ++i) {
                ASSERT_EQ(y[i], ref[i]) << "chunk " << chunk << " output " << i;
            }
        }
    }

    TEST(Decimate, OutputsForCountsExactly) {
        basic_decimator<float, 3> dec;
        EXPECT_EQ(dec.outputs_for(0), 0U);
        EXPECT_EQ(dec.outputs_for(1), 1U); // x[0] emits
        EXPECT_EQ(dec.outputs_for(3), 1U);
        EXPECT_EQ(dec.outputs_for(4), 2U);
        float              in[4] = {0, 0, 0, 0};
        std::vector<float> out(2);
        EXPECT_EQ(dec.process(in, 4, out.data()), 2U); // consumed x[0..3]; phase now 1
        EXPECT_EQ(dec.outputs_for(1), 0U);             // x[4]
        EXPECT_EQ(dec.outputs_for(2), 0U);             // x[5]
        EXPECT_EQ(dec.outputs_for(3), 1U);             // x[6] emits
        EXPECT_EQ(dec.outputs_for(6), 2U);             // x[9] too
    }

    TYPED_TEST(decimate_test, DcGainIsUnity) {
        using sample = TypeParam;
        basic_decimator<sample, 6> dec;
        std::vector<sample>        x(dec.taps() * 2, sample{0.5});
        std::vector<sample>        y(dec.outputs_for(x.size()));
        dec.process(x.data(), x.size(), y.data());
        EXPECT_NEAR(y.back(), 0.5, 1e-6);
    }

    TEST(Decimate, DcGainIsExactInQ15) {
        basic_decimator<std::int16_t, 3> q15;
        std::vector<std::int16_t>        xq(q15.taps() * 2, std::int16_t{16384});
        std::vector<std::int16_t>        yq(q15.outputs_for(xq.size()));
        q15.process(xq.data(), xq.size(), yq.data());
        EXPECT_EQ(yq.back(), std::int16_t{16384}) << "row-sum-preserving quantization keeps DC exact";
    }

    template <std::size_t M>
    double response_db(const basic_decimator<float, M>& dec, double f_hz) {
        const double         fs = basic_decimator<float, M>::traits::k_input_rate_hz;
        std::complex<double> acc{0.0, 0.0};
        const auto           h = dec.coefficients();
        for (std::size_t n = 0; n < h.size(); ++n) {
            const double w = -2.0 * std::numbers::pi * f_hz * static_cast<double>(n) / fs;
            acc += static_cast<double>(h[n]) * std::complex<double>{std::cos(w), std::sin(w)};
        }
        return 20.0 * std::log10(std::abs(acc) + 1e-300);
    }

    template <std::size_t M>
    void expect_profile_met(const decimate_profile& p) {
        basic_decimator<float, M> dec(p);
        const double              fs   = basic_decimator<float, M>::traits::k_input_rate_hz;
        const double              stop = 16000.0 - p.passband_hz;
        for (double f = 0.0; f <= p.passband_hz; f += 25.0) {
            EXPECT_NEAR(response_db(dec, f), 0.0, 0.1) << "passband at " << f << " Hz, ratio " << M;
        }
        for (double f = stop; f <= fs / 2.0; f += 25.0) {
            EXPECT_LT(response_db(dec, f), -p.stopband_atten_db) << "stopband at " << f << " Hz, ratio " << M;
        }
    }

    TEST(Decimate, EconomyMeetsItsNumbers) {
        expect_profile_met<2>(decimate_profile::economy());
        expect_profile_met<3>(decimate_profile::economy());
        expect_profile_met<6>(decimate_profile::economy());
    }

    TEST(Decimate, TransparentMeetsItsNumbers) {
        expect_profile_met<2>(decimate_profile::transparent());
        expect_profile_met<3>(decimate_profile::transparent());
        expect_profile_met<6>(decimate_profile::transparent());
    }

    TYPED_TEST(decimate_test, GroupDelayIsHalfTheTaps) {
        using sample = TypeParam;
        basic_decimator<sample, 3> dec;
        EXPECT_EQ(dec.latency_input_samples(), (dec.taps() - 1) / 2);
        std::vector<sample> x(dec.taps() + 3, sample{0});
        x[0] = sample{1};
        std::vector<sample> y(dec.outputs_for(x.size()));
        dec.process(x.data(), x.size(), y.data());
        const auto peak = static_cast<std::size_t>(std::max_element(y.begin(), y.end()) - y.begin());
        // y[k] = h[3k]; the centre tap (index 60) lands on output 20 exactly.
        EXPECT_EQ(peak * 3, dec.latency_input_samples());
        EXPECT_EQ(y[peak], dec.coefficients()[dec.latency_input_samples()]);
    }

    // The fixed-point profiles against the double golden model (the house
    // rule: cross-precision compares to double, never to a sibling profile).
    // The input is scaled to half full scale: at 0.9 the output overshoots
    // 1.0 on this noise and the fixed formats saturate, which is their
    // contract, not tracking error. Returns the worst deviation in full-scale
    // units.
    template <typename S>
    double worst_deviation_from_double(double full_scale) {
        const auto&                x = frontend_ref::k_dec_input;
        basic_decimator<double, 3> fd;
        basic_decimator<S, 3>      fq;
        std::vector<double>        xd(x.size());
        std::vector<S>             xq(x.size());
        for (std::size_t i = 0; i < x.size(); ++i) {
            xd[i] = 0.5 * static_cast<double>(x[i]);
            xq[i] = static_cast<S>(std::llround(xd[i] * full_scale));
        }
        std::vector<double> yd(fd.outputs_for(xd.size()));
        std::vector<S>      yq(fq.outputs_for(xq.size()));
        fd.process(xd.data(), xd.size(), yd.data());
        fq.process(xq.data(), xq.size(), yq.data());
        double worst = 0.0;
        for (std::size_t i = 0; i < yd.size(); ++i) {
            worst = std::max(worst, std::abs(static_cast<double>(yq[i]) / full_scale - yd[i]));
        }
        return worst;
    }

    TEST(Decimate, Q15TracksDoubleWithinTheFormatFloor) {
        // Measured 2026-09 at 1.57e-4 (Q1.14 coefficient rounding over 121
        // taps, plus one output LSB); pinned at 2x. The same number as the
        // former Q15-vs-float pin: float sits 2.5e-8 from double here.
        EXPECT_LT(worst_deviation_from_double<std::int16_t>(32768.0), 3.2e-4) << "Q15 vs double, in full-scale units";
    }

    TEST(Decimate, Q31TracksDoubleWithinTheFormatFloor) {
        // Measured 2026-09 at 2.45e-9 (Q1.30 coefficient rounding over 121
        // taps, the 16-bit product pre-shift, one output LSB); pinned at 2x.
        EXPECT_LT(worst_deviation_from_double<std::int32_t>(2147483648.0), 4.9e-9)
            << "Q31 vs double, in full-scale units";
    }

    // FNV-1a-64 over the coefficient bytes in memory order (little-endian on
    // every CI host and target; the pin is byte-wise so a one-LSB move in
    // any tap changes it).
    template <typename C>
    std::uint64_t fnv1a64(std::span<const C> c) {
        std::uint64_t h = 0xcbf29ce484222325ULL;
        for (const C v : c) {
            unsigned char bytes[sizeof(C)];
            std::memcpy(bytes, &v, sizeof bytes);
            for (const unsigned char b : bytes) {
                h ^= b;
                h *= 0x100000001b3ULL;
            }
        }
        return h;
    }

    template <typename S, std::size_t M>
    void expect_table_pinned(const decimate_profile& p, std::size_t taps, std::int64_t sum, std::uint64_t fnv) {
        basic_decimator<S, M> dec(p);
        const auto            c = dec.coefficients();
        ASSERT_EQ(c.size(), taps);
        std::int64_t s = 0;
        for (const auto v : c) {
            s += v;
        }
        EXPECT_EQ(s, sum) << "DC gain exactly 1 in the format's unity";
        EXPECT_EQ(fnv1a64(c), fnv) << "format " << sizeof(S) * 8 << " ratio " << M << " taps " << taps;
    }

    // The committed bit pin for the fixed-point tables: every Q15 and Q31
    // coefficient table, both profiles, all three ratios — row sum (the
    // "DC gain exactly 1" contract, 2^14 / 2^30) and a byte-wise FNV-1a-64.
    // Measured 2026-09 on the Stage 3a substrate and identical on the
    // substrate before it (bare Q-ladder literals): this is the evidence
    // that the named constants changed no bit, kept in-tree so a later
    // kernel or backend change that moves a coefficient LSB inside the
    // tolerance pins is caught here.
    TEST(Decimate, FixedPointTablesAreBitPinned) {
        const auto eco = decimate_profile::economy();
        const auto tra = decimate_profile::transparent();
        expect_table_pinned<std::int16_t, 2>(eco, 81, 16384, 0x9c824c3cc6602f9dULL);
        expect_table_pinned<std::int16_t, 2>(tra, 259, 16384, 0x57fcb157f4ee1c25ULL);
        expect_table_pinned<std::int16_t, 3>(eco, 121, 16384, 0x6b81f5469b8d019cULL);
        expect_table_pinned<std::int16_t, 3>(tra, 389, 16384, 0xbb59dcc474a57d61ULL);
        expect_table_pinned<std::int16_t, 6>(eco, 239, 16384, 0x8810d24884cf9bc5ULL);
        expect_table_pinned<std::int16_t, 6>(tra, 773, 16384, 0xa69d8469e463b919ULL);
        expect_table_pinned<std::int32_t, 2>(eco, 81, 1073741824, 0x8ba26dffc4b03946ULL);
        expect_table_pinned<std::int32_t, 2>(tra, 259, 1073741824, 0x980e87e92485f881ULL);
        expect_table_pinned<std::int32_t, 3>(eco, 121, 1073741824, 0xa8406e283f4cef8eULL);
        expect_table_pinned<std::int32_t, 3>(tra, 389, 1073741824, 0x3550a72e43996bd7ULL);
        expect_table_pinned<std::int32_t, 6>(eco, 239, 1073741824, 0xe20c46ad95d31338ULL);
        expect_table_pinned<std::int32_t, 6>(tra, 773, 1073741824, 0x4cd54a0678ad8bcaULL);
    }

    TEST(Decimate, Q15SaturatesInsteadOfWrapping) {
        // A full-scale step: the linear-phase lowpass pre-rings below zero by up
        // to the Gibbs undershoot (about 9 % of the step; measured 6.8 % here),
        // then settles at full scale. Saturation, never a wrap: no output goes
        // anywhere near the wrapped value of an overflow.
        basic_decimator<std::int16_t, 2> dec;
        std::vector<std::int16_t>        x(dec.taps() * 2, std::int16_t{32767});
        std::vector<std::int16_t>        y(dec.outputs_for(x.size()));
        dec.process(x.data(), x.size(), y.data());
        for (const auto v : y) {
            EXPECT_GT(v, std::int16_t{-3300}) << "Gibbs ringing only, no wrap-around";
        }
        EXPECT_EQ(y.back(), std::int16_t{32767});
    }

} // namespace
