// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// Locks down the tap::dsp::basic_log_mel contract: sample-for-sample
// agreement with the committed numpy reference on both the plain-log and
// PCEN paths, chunking invariance of the streaming API, frame alignment and
// latency, the filterbank's stated formulas (HTK mel, unit-peak triangles,
// DC excluded), PCEN's gain tracking and reset semantics, geometry
// validation, and float/double cross-precision agreement as a measured
// number. MuTap's keyword spotter is trained against exactly these numbers.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <numbers>
#include <vector>

#include <gtest/gtest.h>

#include "reference/frontend_vectors.h"
#include "tap/dsp/log_mel.h"

namespace {

    using tap::dsp::basic_log_mel;
    using tap::dsp::log_mel_geometry;

    // The reference generator's signal (tools/reference/make_frontend_reference.py):
    // six tones plus xorshift32 noise, +20 dB from the midpoint. Built in double
    // and cast, so both profiles see the same samples to within their rounding.
    constexpr double k_sr        = 16000.0;
    constexpr size_t k_n         = 8000;
    constexpr double k_noise_amp = 0.02;
    constexpr double k_step_gain = 10.0;

    struct tone {
        double freq;
        double amp;
    };
    constexpr tone k_tones[] = {{150.0, 0.30},  {440.0, 0.25},  {1000.0, 0.20},
                                {2500.0, 0.15}, {4000.0, 0.10}, {6500.0, 0.05}};

    class xorshift32 {
      public:
        explicit xorshift32(std::uint32_t seed)
            : m_s(seed) {}
        double next() noexcept {
            m_s ^= m_s << 13;
            m_s ^= m_s >> 17;
            m_s ^= m_s << 5;
            return (static_cast<double>(m_s % 65536U) - 32768.0) / 32768.0;
        }

      private:
        std::uint32_t m_s;
    };

    std::vector<double> reference_signal() {
        std::vector<double> x(k_n, 0.0);
        for (const tone& t : k_tones) {
            for (size_t n = 0; n < k_n; ++n) {
                x[n] += t.amp * std::sin(2.0 * std::numbers::pi * t.freq * static_cast<double>(n) / k_sr);
            }
        }
        xorshift32 rng(0x2545F491U);
        for (size_t n = 0; n < k_n; ++n) {
            x[n] += k_noise_amp * rng.next();
        }
        for (size_t n = k_n / 2; n < k_n; ++n) {
            x[n] *= k_step_gain;
        }
        return x;
    }

    template <typename Sample>
    std::vector<Sample> as(const std::vector<double>& x) {
        std::vector<Sample> y(x.size());
        std::transform(x.begin(), x.end(), y.begin(), [](double v) { return static_cast<Sample>(v); });
        return y;
    }

    template <typename Sample>
    std::vector<Sample> run(basic_log_mel<Sample>& fe, const std::vector<Sample>& x) {
        std::vector<Sample> out(fe.frames_for(x.size()) * fe.bands());
        const size_t        n = fe.process(x.data(), x.size(), out.data(), fe.frames_for(x.size()));
        out.resize(n * fe.bands());
        return out;
    }

    // Pinned cross-language tolerances (feature units; the plain-log features
    // span [-1, 1]). Measured 2026-09 on this signal, C++ vs the numpy float64
    // reference: double 1.5e-14 (log) / 3.4e-14 (PCEN) — FFT and log rounding
    // only; float 6.7e-7 (log) / 5.3e-6 (PCEN). Pinned at ~2x the worse path.
    template <typename Sample>
    constexpr double reference_tolerance() {
        return std::is_same_v<Sample, double> ? 1e-13 : 1.2e-5;
    }

    template <typename Sample>
    class log_mel_test : public ::testing::Test {};
    using profiles = ::testing::Types<float, double>;
    TYPED_TEST_SUITE(log_mel_test, profiles);

    TYPED_TEST(log_mel_test, MatchesNumpyReferenceLogPath) {
        basic_log_mel<TypeParam> fe{log_mel_geometry{}};
        const auto               out = run(fe, as<TypeParam>(reference_signal()));
        ASSERT_EQ(out.size(), frontend_ref::k_mel_log.size());
        double worst = 0.0;
        for (size_t i = 0; i < out.size(); ++i) {
            worst = std::max(worst, std::abs(static_cast<double>(out[i]) - frontend_ref::k_mel_log[i]));
        }
        EXPECT_LT(worst, reference_tolerance<TypeParam>()) << "max |C++ - numpy| on the log path";
    }

    TYPED_TEST(log_mel_test, MatchesNumpyReferencePcenPath) {
        log_mel_geometry g;
        g.pcen.enabled = true;
        basic_log_mel<TypeParam> fe{g};
        const auto               out = run(fe, as<TypeParam>(reference_signal()));
        ASSERT_EQ(out.size(), frontend_ref::k_mel_pcen.size());
        double worst = 0.0;
        for (size_t i = 0; i < out.size(); ++i) {
            worst = std::max(worst, std::abs(static_cast<double>(out[i]) - frontend_ref::k_mel_pcen[i]));
        }
        EXPECT_LT(worst, reference_tolerance<TypeParam>()) << "max |C++ - numpy| on the PCEN path";
    }

    TYPED_TEST(log_mel_test, ChunkingNeverChangesAFeature) {
        const auto               x = as<TypeParam>(reference_signal());
        basic_log_mel<TypeParam> whole{log_mel_geometry{}};
        const auto               ref = run(whole, x);
        for (const size_t chunk : {size_t{1}, size_t{7}, size_t{160}, size_t{333}, size_t{4096}}) {
            basic_log_mel<TypeParam> fe{log_mel_geometry{}};
            std::vector<TypeParam>   out(ref.size(), TypeParam(0));
            size_t                   frames = 0;
            for (size_t i = 0; i < x.size(); i += chunk) {
                const size_t n = std::min(chunk, x.size() - i);
                frames += fe.process(x.data() + i, n, out.data() + frames * fe.bands(), fe.frames_for(n));
            }
            ASSERT_EQ(frames * fe.bands(), ref.size()) << "chunk " << chunk;
            for (size_t i = 0; i < ref.size(); ++i) {
                ASSERT_EQ(out[i], ref[i]) << "chunk " << chunk << " feature " << i;
            }
        }
    }

    TYPED_TEST(log_mel_test, ProcessHopEqualsProcess) {
        const auto               x = as<TypeParam>(reference_signal());
        basic_log_mel<TypeParam> a{log_mel_geometry{}};
        basic_log_mel<TypeParam> b{log_mel_geometry{}};
        const auto               ref = run(a, x);
        std::vector<TypeParam>   frame(b.bands());
        for (size_t t = 0; t < ref.size() / b.bands(); ++t) {
            b.process_hop(x.data() + t * b.hop(), frame.data());
            for (size_t k = 0; k < b.bands(); ++k) {
                ASSERT_EQ(frame[k], ref[t * b.bands() + k]);
            }
        }
    }

    TYPED_TEST(log_mel_test, FrameAlignmentAndLatency) {
        // An impulse at sample 0 lies inside frames 0 and 1 only: frame t spans
        // [(t+1)*hop - frame, (t+1)*hop), i.e. [-240, 160) and [-80, 320).
        log_mel_geometry         g;
        basic_log_mel<TypeParam> fe{g};
        EXPECT_EQ(fe.latency_samples(), g.frame);
        std::vector<TypeParam> x(g.hop * 4, TypeParam(0));
        x[0]                       = TypeParam(1);
        const auto   out           = run(fe, x);
        const double floor_feature = (std::log10(g.log_floor) + g.log_shift) / g.log_scale; // -1
        auto         frame_energy  = [&](size_t t) {
            double e = 0.0;
            for (size_t b = 0; b < g.bands; ++b) {
                e += static_cast<double>(out[t * g.bands + b]) - floor_feature;
            }
            return e;
        };
        EXPECT_GT(frame_energy(0), 1.0);
        EXPECT_GT(frame_energy(1), 1.0);
        EXPECT_NEAR(frame_energy(2), 0.0, 1e-6);
        EXPECT_NEAR(frame_energy(3), 0.0, 1e-6);
    }

    TYPED_TEST(log_mel_test, FilterbankIsHtkMelWithUnitPeakTrianglesAndNoDc) {
        log_mel_geometry         g;
        basic_log_mel<TypeParam> fe{g};
        EXPECT_NEAR(fe.band_edge_hz(0), g.fmin_hz, 1e-9);
        EXPECT_NEAR(fe.band_edge_hz(g.bands + 1), g.fmax_hz, 1e-9);
        // Edges are equally spaced in HTK mel.
        const double lo = 2595.0 * std::log10(1.0 + g.fmin_hz / 700.0);
        const double hi = 2595.0 * std::log10(1.0 + g.fmax_hz / 700.0);
        for (size_t i = 0; i <= g.bands + 1; ++i) {
            const double mel = lo + (hi - lo) * static_cast<double>(i) / static_cast<double>(g.bands + 1);
            EXPECT_NEAR(fe.band_edge_hz(i), 700.0 * (std::pow(10.0, mel / 2595.0) - 1.0), 1e-9);
        }
        for (size_t b = 0; b < g.bands; ++b) {
            EXPECT_EQ(fe.band_weight(b, 0), TypeParam(0)) << "DC weight in band " << b;
            TypeParam peak = TypeParam(0);
            for (size_t k = 0; k < g.bins(); ++k) {
                peak = std::max(peak, fe.band_weight(b, k));
            }
            EXPECT_LE(peak, TypeParam(1));
            EXPECT_GT(peak, TypeParam(0.25)) << "band " << b << " has no bin near its centre";
        }
        // A bin sitting exactly on a band centre carries weight 1: band 35's centre
        // at the default geometry is not on the grid, so test the formula directly.
        const double f_mid   = fe.band_edge_hz(21);
        const double f_lo    = fe.band_edge_hz(20);
        const auto   k_below = static_cast<size_t>(std::floor(f_mid * static_cast<double>(g.fft_size) / g.sample_rate));
        const double f_k     = static_cast<double>(k_below) * g.sample_rate / static_cast<double>(g.fft_size);
        EXPECT_NEAR(static_cast<double>(fe.band_weight(20, k_below)), (f_k - f_lo) / (f_mid - f_lo), 1e-6);
    }

    TYPED_TEST(log_mel_test, PcenTracksALevelStep) {
        // A stationary tone at level A for 200 frames then 10 A for 200 frames.
        // With the smoother settled (M = E) the recursion's steady state is the
        // closed form out(E) = (E / (eps + E)^alpha + delta)^r - delta^r, which
        // keeps a deliberate residual level dependence of E^0.02 at alpha 0.98. Both
        // halves are pinned to that closed form after five time constants, and
        // the step's first frame must stand well above it.
        log_mel_geometry g;
        g.pcen.enabled = true;
        basic_log_mel<TypeParam> fe{g};
        basic_log_mel<TypeParam> energies{log_mel_geometry{}}; // the log path recovers E
        const size_t             frames = 400;
        std::vector<TypeParam>   x(frames * g.hop);
        for (size_t n = 0; n < x.size(); ++n) {
            const double a = n < x.size() / 2 ? 0.05 : 0.5;
            x[n] =
                static_cast<TypeParam>(a * std::sin(2.0 * std::numbers::pi * 1000.0 * static_cast<double>(n) / k_sr));
        }
        const auto out = run(fe, x);
        const auto lg  = run(energies, x);
        // The band holding the tone: the loudest band of a settled quiet frame.
        size_t b = 0;
        for (size_t k = 1; k < g.bands; ++k) {
            if (lg[100 * g.bands + k] > lg[100 * g.bands + b]) {
                b = k;
            }
        }
        auto energy_at = [&](size_t t) {
            return std::pow(10.0, static_cast<double>(lg[t * g.bands + b]) * g.log_scale - g.log_shift) - g.log_floor;
        };
        auto steady = [&](double e) {
            return std::pow(e / std::pow(g.pcen.epsilon + e, g.pcen.alpha) + g.pcen.delta, g.pcen.power)
                   - std::pow(g.pcen.delta, g.pcen.power);
        };
        const double quiet = out[(frames / 2 - 1) * g.bands + b];
        const double loud  = out[(frames - 1) * g.bands + b];
        const double jump  = out[(frames / 2) * g.bands + b];
        EXPECT_NEAR(quiet, steady(energy_at(frames / 2 - 1)), 0.01 * quiet);
        EXPECT_NEAR(loud, steady(energy_at(frames - 1)), 0.01 * loud);
        EXPECT_GT(jump, 3.0 * quiet) << "the +20 dB step is visible on its first frame";
    }

    TYPED_TEST(log_mel_test, ResetMakesTheNextFrameAFirstFrame) {
        log_mel_geometry g;
        g.pcen.enabled             = true;
        g.preemphasis              = 0.97;
        const auto               x = as<TypeParam>(reference_signal());
        basic_log_mel<TypeParam> fresh{g};
        basic_log_mel<TypeParam> reused{g};
        const auto               ref = run(fresh, x);
        (void)run(reused, x);
        reused.reset();
        const auto again = run(reused, x);
        ASSERT_EQ(again.size(), ref.size());
        for (size_t i = 0; i < ref.size(); ++i) {
            ASSERT_EQ(again[i], ref[i]) << "feature " << i;
        }
    }

    TYPED_TEST(log_mel_test, PreemphasisIsAppliedOnTheStream) {
        // With c = 1 a constant input pre-emphasizes to a single impulse at
        // sample 0 (x[-1] = 0), so only the frames containing sample 0 have energy.
        log_mel_geometry g;
        g.preemphasis = 1.0;
        basic_log_mel<TypeParam> fe{g};
        std::vector<TypeParam>   x(g.hop * 4, TypeParam(0.5));
        const auto               out           = run(fe, x);
        const double             floor_feature = (std::log10(g.log_floor) + g.log_shift) / g.log_scale;
        double                   tail          = 0.0;
        for (size_t b = 0; b < g.bands; ++b) {
            tail += static_cast<double>(out[3 * g.bands + b]) - floor_feature;
        }
        EXPECT_NEAR(tail, 0.0, 1e-6);
    }

    TYPED_TEST(log_mel_test, GeometryValidation) {
        EXPECT_TRUE(log_mel_geometry{}.valid());
        log_mel_geometry g;
        g.fft_size = 400; // not a power of two
        EXPECT_FALSE(g.valid());
        g          = log_mel_geometry{};
        g.fft_size = 256; // smaller than the frame
        EXPECT_FALSE(g.valid());
        g     = log_mel_geometry{};
        g.hop = 401; // larger than the frame
        EXPECT_FALSE(g.valid());
        g         = log_mel_geometry{};
        g.fmax_hz = 9000.0; // above Nyquist
        EXPECT_FALSE(g.valid());
        g           = log_mel_geometry{};
        g.log_floor = 0.0;
        EXPECT_FALSE(g.valid());
        g               = log_mel_geometry{};
        g.pcen.smoother = 0.0;
        EXPECT_FALSE(g.valid());
        EXPECT_EQ(log_mel_geometry::k_contract_version, 1U);
    }

    TEST(LogMelCrossPrecision, FloatTracksDouble) {
        // Measured 2026-09 on the reference signal: 6.7e-7 (log), 5.3e-6 (PCEN).
        // Pinned at 2x so a regression in the float path is caught, not absorbed.
        const auto xd = reference_signal();
        const auto xf = as<float>(xd);
        for (const bool pcen : {false, true}) {
            log_mel_geometry g;
            g.pcen.enabled = pcen;
            basic_log_mel<double> fd{g};
            basic_log_mel<float>  ff{g};
            const auto            od    = run(fd, xd);
            const auto            of    = run(ff, xf);
            double                worst = 0.0;
            for (size_t i = 0; i < od.size(); ++i) {
                worst = std::max(worst, std::abs(od[i] - static_cast<double>(of[i])));
            }
            EXPECT_LT(worst, pcen ? 1.1e-5 : 1.4e-6) << (pcen ? "PCEN" : "log") << " path";
        }
    }

} // namespace
