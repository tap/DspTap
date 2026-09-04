// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// Locks down the tap::dsp::basic_decimator contract: the pinned tap counts,
// sample-for-sample agreement of the float golden model with the committed
// numpy reference for every ratio and profile, chunking invariance, exact
// output counting, unity DC gain (exact in Q15 through row-sum-preserving
// quantization), the stated passband/stopband numbers measured from the
// shipped coefficients, the integer group delay, and the Q15 profile
// tracking float within the format's own floor.

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>
#include <numbers>
#include <vector>

#include <gtest/gtest.h>

#include "reference/frontend_vectors.h"
#include "tap/dsp/decimate.h"

namespace {

    using tap::dsp::basic_decimator;
    using tap::dsp::decimate_profile;

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

    TEST(Decimate, ChunkingIsBitIdentical) {
        const auto&               x = frontend_ref::k_dec_input;
        basic_decimator<float, 3> whole;
        std::vector<float>        ref(whole.outputs_for(x.size()));
        whole.process(x.data(), x.size(), ref.data());
        for (const std::size_t chunk :
             {std::size_t{1}, std::size_t{2}, std::size_t{7}, std::size_t{255}, std::size_t{1000}}) {
            basic_decimator<float, 3> dec;
            std::vector<float>        y;
            for (std::size_t i = 0; i < x.size(); i += chunk) {
                const std::size_t  n = std::min(chunk, x.size() - i);
                std::vector<float> part(dec.outputs_for(n));
                ASSERT_EQ(dec.process(x.data() + i, n, part.data()), part.size());
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

    TEST(Decimate, DcGainIsUnity) {
        basic_decimator<float, 6> dec;
        std::vector<float>        x(dec.taps() * 2, 0.5f);
        std::vector<float>        y(dec.outputs_for(x.size()));
        dec.process(x.data(), x.size(), y.data());
        EXPECT_NEAR(y.back(), 0.5f, 1e-6f);

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

    TEST(Decimate, GroupDelayIsHalfTheTaps) {
        basic_decimator<float, 3> dec;
        EXPECT_EQ(dec.latency_input_samples(), (dec.taps() - 1) / 2);
        std::vector<float> x(dec.taps() + 3, 0.0f);
        x[0] = 1.0f;
        std::vector<float> y(dec.outputs_for(x.size()));
        dec.process(x.data(), x.size(), y.data());
        const auto peak = static_cast<std::size_t>(std::max_element(y.begin(), y.end()) - y.begin());
        // y[k] = h[3k]; the centre tap (index 60) lands on output 20 exactly.
        EXPECT_EQ(peak * 3, dec.latency_input_samples());
        EXPECT_FLOAT_EQ(y[peak], dec.coefficients()[dec.latency_input_samples()]);
    }

    TEST(Decimate, Q15TracksFloatWithinTheFormatFloor) {
        // The input is scaled to half full scale: at 0.9 the float output
        // overshoots 1.0 on this noise and Q15 saturates, which is the format's
        // contract, not tracking error. Measured 2026-09 at 1.6e-4 (Q1.14
        // coefficient rounding over 121 taps, plus one output LSB); pinned at 2x.
        const auto&                      x = frontend_ref::k_dec_input;
        basic_decimator<float, 3>        ff;
        basic_decimator<std::int16_t, 3> fq;
        std::vector<float>               xh(x.size());
        std::vector<std::int16_t>        xq(x.size());
        for (std::size_t i = 0; i < x.size(); ++i) {
            xh[i] = 0.5f * x[i];
            xq[i] = static_cast<std::int16_t>(std::lround(xh[i] * 32768.0f));
        }
        std::vector<float>        yf(ff.outputs_for(xh.size()));
        std::vector<std::int16_t> yq(fq.outputs_for(xq.size()));
        ff.process(xh.data(), xh.size(), yf.data());
        fq.process(xq.data(), xq.size(), yq.data());
        double worst = 0.0;
        for (std::size_t i = 0; i < yf.size(); ++i) {
            worst = std::max(worst, std::abs(static_cast<double>(yq[i]) / 32768.0 - static_cast<double>(yf[i])));
        }
        EXPECT_LT(worst, 3.2e-4) << "Q15 vs float, in full-scale units";
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
