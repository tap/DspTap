// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// Pins chain.h's latency members for chains holding recursive stages
// (SampleRateTap's docs/IIR_PLAN.md, I5) and a 1:1 stage in a chain (I7's
// note), on test-local stages built over tap/dsp/iir_kernels.h, as chain.h
// was proven on decimate.h before an engine depended on it:
//   - nominal_latency_output_frames(), latency_deviation_output_frames() and
//     group_delay_output_frames(f) compose over the stages, each stage's
//     converted to the chain's output frames and its frequency to the stage's
//     own; an exact stage contributes its exact latency at every frequency;
//     the chain's group delay matches the delay measured from its impulse
//     response;
//   - an all-exact chain's nominal is its exact latency, its deviation 0;
//   - a near-linear stage's delay path is its nominal, and a chain of such
//     stages stays an exact_delay_stage for subband.h (whose battery runs
//     unchanged); a chain with a nonlinear stage is not one, and its exact
//     latency_output_frames() does not compile
//     (tests/compile_fail/chain_iir_exact_latency.cpp);
//   - a 1:1 stage (the guard band filter's shape) composes: the ratio,
//     outputs_for from every phase, frames_needed, bit-identical chunking and
//     flush equal to zero padding, in every sample format.
// The float and double cases run under the host's flush-to-zero (the
// fixture), as I9 asks of every recursive test. Every case runs on the QEMU
// legs.

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>
#include <limits>
#include <numbers>
#include <span>
#include <type_traits>
#include <vector>

#include <gtest/gtest.h>

#include "support/iir_support.h"
#include "tap/dsp/chain.h"
#include "tap/dsp/decimate.h"
#include "tap/dsp/iir_design.h"
#include "tap/dsp/iir_kernels.h"
#include "tap/dsp/subband.h"

namespace {

    using namespace tap::dsp;
    using tap::dsp::test::flush_to_zero_fixture;

    /// The largest |tau(f) - tau(0)| of a delay function over [0, edge].
    template <typename F>
    double deviation_over(F tau, double edge) {
        const double t0    = tau(0.0);
        double       worst = 0.0;
        for (int g = 0; g <= 1024; ++g) {
            worst = std::max(worst, std::fabs(tau(edge * g / 1024.0) - t0));
        }
        return worst;
    }

    // A nonlinear half-band decimator: output k as input 2k arrives (decimate.h's
    // alignment), y[k] = (H x)[2k]. Its nominal latency is its DC group delay.
    template <typename S>
    class iir_down2 {
      public:
        using sample                        = S;
        static constexpr std::size_t k_up   = 1;
        static constexpr std::size_t k_down = 2;
        static constexpr std::size_t channels() noexcept { return 1; }

        explicit iir_down2(const nonlinear_halfband_design& d, std::size_t window = 2048)
            : m_design(d)
            , m_table(d)
            , m_state(m_table.state_size(), iir_state<S>{})
            , m_window(window)
            , m_deviation(deviation_over([&](double f) { return halfband_group_delay(d, f); }, d.f_pass) / 2.0) {}

        std::size_t process(const S* in, std::size_t n, S* out) noexcept {
            std::size_t made = 0;
            for (std::size_t i = 0; i < n; ++i) {
                if (m_phase == 0) {
                    out[made++] = halfband_decimate(m_table, m_state.data(), in[i], m_odd);
                }
                else {
                    m_odd = in[i];
                }
                m_phase ^= 1u;
            }
            return made;
        }
        std::size_t outputs_for(std::size_t n) const noexcept { return m_phase == 0 ? (n + 1) / 2 : n / 2; }
        void        reset() noexcept {
            std::fill(m_state.begin(), m_state.end(), iir_state<S>{});
            m_phase = 0;
            m_odd   = S{};
        }
        std::size_t window_frames() const noexcept { return m_window; }

        double nominal_latency_output_frames() const noexcept { return halfband_group_delay(m_design, 0.0) / 2.0; }
        double latency_deviation_output_frames() const noexcept { return m_deviation; }
        double group_delay_output_frames(double f) const noexcept {
            return halfband_group_delay(m_design, f / 2.0) / 2.0;
        }

      private:
        nonlinear_halfband_design   m_design;
        nonlinear_halfband_table<S> m_table;
        std::vector<iir_state<S>>   m_state;
        std::size_t                 m_window;
        double                      m_deviation;
        unsigned                    m_phase = 0;
        S                           m_odd{};
    };

    // A near-linear half-band decimator: its delay path, K output frames, is
    // exact (subband.h's delay_path_output_frames()) and is its nominal.
    template <typename S>
    class near_linear_down2 {
      public:
        using sample                        = S;
        static constexpr std::size_t k_up   = 1;
        static constexpr std::size_t k_down = 2;
        static constexpr std::size_t channels() noexcept { return 1; }

        explicit near_linear_down2(const near_linear_halfband_design& d)
            : m_design(&d)
            , m_table(d)
            , m_state(m_table.state_size(), iir_state<S>{})
            , m_evens(m_table.delay + 1, S{}) {}

        std::size_t process(const S* in, std::size_t n, S* out) noexcept {
            std::size_t made = 0;
            for (std::size_t i = 0; i < n; ++i) {
                if (m_phase == 0) {
                    m_evens[m_at] = in[i];
                    m_at          = m_at + 1 == m_evens.size() ? 0 : m_at + 1;
                    // the oldest of the K + 1 held evens is x[2(k - K)]
                    out[made++] = near_linear_decimate(m_table, m_state.data(), m_evens[m_at], m_odd);
                }
                else {
                    m_odd = in[i];
                }
                m_phase ^= 1u;
            }
            return made;
        }
        std::size_t outputs_for(std::size_t n) const noexcept { return m_phase == 0 ? (n + 1) / 2 : n / 2; }
        void        reset() noexcept {
            std::fill(m_state.begin(), m_state.end(), iir_state<S>{});
            std::fill(m_evens.begin(), m_evens.end(), S{});
            m_at    = 0;
            m_phase = 0;
            m_odd   = S{};
        }
        std::size_t window_frames() const noexcept { return 2048; }

        exact_ratio delay_path_output_frames() const noexcept {
            return exact_ratio{static_cast<std::uint64_t>(m_design->delay), 1};
        }
        double nominal_latency_output_frames() const noexcept { return static_cast<double>(m_design->delay); }
        double latency_deviation_output_frames() const noexcept { return m_design->delay_deviation / 2.0; }
        double group_delay_output_frames(double f) const noexcept {
            return halfband_group_delay(*m_design, f / 2.0) / 2.0;
        }

      private:
        const near_linear_halfband_design* m_design;
        near_linear_halfband_table<S>      m_table;
        std::vector<iir_state<S>>          m_state;
        std::vector<S>                     m_evens;
        std::size_t                        m_at    = 0;
        unsigned                           m_phase = 0;
        S                                  m_odd{};
    };

    // A 1:1 stage: the guard band filter's Chebyshev II, one biquad cascade.
    template <typename S>
    class guard_1to1 {
      public:
        using sample                        = S;
        static constexpr std::size_t k_up   = 1;
        static constexpr std::size_t k_down = 1;
        static constexpr std::size_t channels() noexcept { return 1; }

        guard_1to1(const cheby2_design& d, double passband)
            : m_design(d)
            , m_table(d)
            , m_state(m_table.state_size(), iir_state<S>{})
            , m_deviation(deviation_over(
                  [&](double f) {
                      return sections_group_delay(std::span<const iir_section>(d.sections.data(), d.count), f);
                  },
                  passband)) {}

        std::size_t process(const S* in, std::size_t n, S* out) noexcept {
            for (std::size_t i = 0; i < n; ++i) {
                out[i] = biquad_cascade(m_table, m_state.data(), in[i]);
            }
            return n;
        }
        std::size_t outputs_for(std::size_t n) const noexcept { return n; }
        void        reset() noexcept { std::fill(m_state.begin(), m_state.end(), iir_state<S>{}); }
        /// Floating: I6's count under the host's flush-to-zero (the fixture);
        /// fixed point: past the measured decay from full scale (535 frames in
        /// Q31, test_iir_fixed.cpp).
        std::size_t window_frames() const noexcept {
            if constexpr (std::is_floating_point_v<S>) {
                return floating_flush_frames<S>(m_design.pole_radius[m_design.count - 1], m_design.count);
            }
            else {
                return 4096;
            }
        }

        double nominal_latency_output_frames() const noexcept { return group_delay_output_frames(0.0); }
        double latency_deviation_output_frames() const noexcept { return m_deviation; }
        double group_delay_output_frames(double f) const noexcept {
            return sections_group_delay(std::span<const iir_section>(m_design.sections.data(), m_design.count), f);
        }

      private:
        cheby2_design             m_design;
        biquad_cascade_table<S>   m_table;
        std::vector<iir_state<S>> m_state;
        double                    m_deviation;
    };

    static_assert(sync_stage<iir_down2<float>>);
    static_assert(sync_stage<near_linear_down2<std::int16_t>>);
    static_assert(sync_stage<guard_1to1<double>>);
    static_assert(reports_nominal_latency<iir_down2<float>>);
    static_assert(reports_nominal_latency<chain<basic_decimator<float, 2>, iir_down2<float>>>,
                  "a chain reports the nominal triple itself, so chains nest");
    // subband.h's view: the near-linear stage's delay path is exact, the
    // nonlinear stage has none (the scaffold refuses it, X10).
    static_assert(exact_delay_stage<near_linear_down2<float>>);
    static_assert(exact_delay_stage<chain<near_linear_down2<float>, near_linear_down2<float>>>);
    static_assert(!exact_delay_stage<iir_down2<float>>);
    static_assert(!exact_delay_stage<chain<basic_decimator<float, 2>, iir_down2<float>>>);

    template <typename S>
    std::vector<S> signal(std::size_t n) {
        std::vector<S> x(n);
        for (std::size_t i = 0; i < n; ++i) {
            const double t = static_cast<double>(i);
            const double v = 0.4 * std::sin(2.0 * std::numbers::pi * 0.0113 * t)
                             + 0.3 * std::cos(2.0 * std::numbers::pi * 0.21 * t)
                             + 0.1 * std::sin(2.0 * std::numbers::pi * 0.47 * t);
            if constexpr (std::is_floating_point_v<S>) {
                x[i] = static_cast<S>(v);
            }
            else {
                x[i] = static_cast<S>(std::llround(v * (std::is_same_v<S, std::int16_t> ? 32767.0 : 2147483647.0)));
            }
        }
        return x;
    }

    /// The group delay of a sequence at f (cycles per sample): Re(sum n h[n]
    /// e^{-jwn} / sum h[n] e^{-jwn}).
    double measured_group_delay(const std::vector<double>& h, double f) {
        std::complex<double> num = 0.0;
        std::complex<double> den = 0.0;
        for (std::size_t n = 0; n < h.size(); ++n) {
            const auto e = std::polar(1.0, -2.0 * std::numbers::pi * f * static_cast<double>(n));
            den += h[n] * e;
            num += static_cast<double>(n) * h[n] * e;
        }
        return std::real(num / den);
    }

    class chain_iir : public flush_to_zero_fixture {};

    TEST_F(chain_iir, NominalLatencyComposesExactAndRecursiveStages) {
        // decimate.h's by-2 economy (81 taps: 40 input samples, 20 of its
        // outputs, 10 of the chain's) then the nonlinear half-band at 70 dB,
        // 3/8, its DC group delay 2.40 high-rate samples, 1.20 output frames.
        using chain_t = chain<basic_decimator<double, 2>, iir_down2<double>>;
        const auto& d = k_nonlinear_halfband<70, 3, 8, 1>;
        chain_t     c(basic_decimator<double, 2>{}, iir_down2<double>(d));
        EXPECT_DOUBLE_EQ(c.nominal_latency_output_frames(), 10.0 + halfband_group_delay(d, 0.0) / 2.0);
        EXPECT_NEAR(c.nominal_latency_output_frames(), 11.20, 0.005);
        EXPECT_DOUBLE_EQ(c.latency_deviation_output_frames(), c.stage<1>().latency_deviation_output_frames());
        // 2.40 .. 4.96 high-rate samples over the passband: 1.28 output frames.
        EXPECT_NEAR(c.latency_deviation_output_frames(), (4.96 - 2.40) / 2.0, 0.005);
        for (const double f : {0.0, 0.05, 0.1, 0.2, 0.3}) {
            EXPECT_DOUBLE_EQ(c.group_delay_output_frames(f), 10.0 + halfband_group_delay(d, f / 2.0) / 2.0) << f;
        }
        // Measured: the chain's impulse response (its input impulse on both
        // stages' phase 0) delays each frequency by what it reports, within
        // what the decimators' aliasing moves a measurement.
        std::vector<double> x(4096, 0.0);
        x[0] = 1.0;
        std::vector<double> y(c.outputs_for(x.size()));
        y.resize(c.process(x.data(), x.size(), y.data()));
        for (const double f : {0.02, 0.1, 0.2, 0.3}) {
            EXPECT_NEAR(measured_group_delay(y, f), c.group_delay_output_frames(f), 0.01) << f;
        }
    }

    TEST_F(chain_iir, StagesAreConvertedToTheChainsFramesAndFrequencies) {
        // Two recursive stages: the first, at design divisor 2, runs at twice
        // the chain's output rate, so its delay counts half and it sees the
        // chain's frequency f at f / 2 of its own output rate.
        using chain_t  = chain<iir_down2<double>, iir_down2<double>>;
        const auto& d2 = k_nonlinear_halfband<70, 3, 8, 2>;
        const auto& d1 = k_nonlinear_halfband<70, 3, 8, 1>;
        chain_t     c{iir_down2<double>(d2), iir_down2<double>(d1)};
        const auto  stage = [](const nonlinear_halfband_design& d, double f) {
            return halfband_group_delay(d, f / 2.0) / 2.0;
        };
        for (const double f : {0.0, 0.1, 0.25, 0.35}) {
            EXPECT_DOUBLE_EQ(c.group_delay_output_frames(f), stage(d2, f / 2.0) / 2.0 + stage(d1, f)) << f;
        }
        EXPECT_DOUBLE_EQ(c.nominal_latency_output_frames(), stage(d2, 0.0) / 2.0 + stage(d1, 0.0));
        EXPECT_DOUBLE_EQ(c.latency_deviation_output_frames(), c.stage<0>().latency_deviation_output_frames() / 2.0
                                                                  + c.stage<1>().latency_deviation_output_frames());
        std::vector<double> x(4096, 0.0);
        x[0] = 1.0;
        std::vector<double> y(c.outputs_for(x.size()));
        y.resize(c.process(x.data(), x.size(), y.data()));
        for (const double f : {0.02, 0.1, 0.2, 0.3}) {
            EXPECT_NEAR(measured_group_delay(y, f), c.group_delay_output_frames(f), 0.01) << f;
        }
    }

    TEST(ChainIirExact, AnAllExactChainsNominalIsItsExactLatency) {
        chain<basic_decimator<float, 2>, basic_decimator<float, 3>> c(basic_decimator<float, 2>{},
                                                                      basic_decimator<float, 3>{});
        EXPECT_EQ(c.latency_output_frames(), (exact_ratio{80, 3}));
        EXPECT_DOUBLE_EQ(c.nominal_latency_output_frames(), 80.0 / 3.0);
        EXPECT_EQ(c.latency_deviation_output_frames(), 0.0);
        for (const double f : {0.0, 0.1, 0.4}) {
            EXPECT_DOUBLE_EQ(c.group_delay_output_frames(f), 80.0 / 3.0) << "a linear-phase FIR delays every f alike";
        }
    }

    TEST_F(chain_iir, ADelayPathStageIsItsOwnNominal) {
        // The near-linear stages at divisors 2 and 1 (70 dB, 3/8): K = 3 then
        // 7, so the chain's delay path is 3/2 + 7 output frames, exact, and
        // its nominal; subband.h reads the same number.
        using stage_t = near_linear_down2<float>;
        using chain_t = chain<stage_t, stage_t>;
        // Braces: MSVC reads stage_t(k_near_linear_halfband<...>) here as a
        // parameter declaration (the most vexing parse).
        chain_t c(stage_t{k_near_linear_halfband<70, 3, 8, 2>}, stage_t{k_near_linear_halfband<70, 3, 8, 1>});
        EXPECT_DOUBLE_EQ(c.nominal_latency_output_frames(), 1.5 + 7.0);
        EXPECT_EQ(detail::exact_delay_output_frames(c), (exact_ratio{17, 2}));
        // The passband delay strays from it by at most the stages' stated
        // deviations (0.0020 and 0.0175 high-rate samples).
        EXPECT_LE(c.latency_deviation_output_frames(), 0.0175 / 2.0 + 0.0020 / 4.0 + 1e-4);
        for (const double f : {0.0, 0.1, 0.2, 0.3, 0.37}) {
            EXPECT_NEAR(c.group_delay_output_frames(f), 8.5, c.latency_deviation_output_frames() + 1e-12) << f;
        }
    }

    // ------------------------------------------------------------------------
    // A 1:1 stage in a chain (I7's note)
    // ------------------------------------------------------------------------

    template <typename S>
    class chain_one_to_one : public flush_to_zero_fixture {};
    using sample_types = ::testing::Types<float, double, std::int16_t, std::int32_t>;
    TYPED_TEST_SUITE(chain_one_to_one, sample_types);

    TYPED_TEST(chain_one_to_one, AOneToOneStageInAChainIsPinned) {
        // The fixed-point formats run the order-7 Chebyshev II at 3/8; the
        // floating ones the order-5 at 1/3, because a floating cascade of
        // order 7 or more never drains to an exactly zero state (it holds a
        // zero-input limit cycle at a few times the smallest normal under
        // flush-to-zero: iir_kernels.h's honest limit,
        // iir_kernel.TheGuardBandCascadeHoldsALimitCycleAtTheSmallestNormal),
        // so no window would make its flush end in silence.
        using sample         = TypeParam;
        using guard_t        = guard_1to1<sample>;
        using chain_t        = chain<basic_decimator<sample, 2>, guard_t>;
        constexpr bool fixed = !std::is_floating_point_v<sample>;
        const auto&    g     = fixed ? k_guard_band_cheby2<70, 3, 8, 7> : k_guard_band_cheby2<70, 1, 3, 5>;
        const double   p     = fixed ? 3.0 / 8.0 : 1.0 / 3.0;
        const auto     fresh = [&] { return chain_t(basic_decimator<sample, 2>{}, guard_t(g, p)); };
        static_assert(chain_t::k_up == 1 && chain_t::k_down == 2, "a 1:1 stage leaves the ratio alone");
        const auto x = signal<sample>(3001);

        // The two run in sequence, whole: the reference.
        basic_decimator<sample, 2> dec;
        guard_t                    grd(g, p);
        std::vector<sample>        mid(dec.outputs_for(x.size()));
        mid.resize(dec.process(x.data(), x.size(), mid.data()));
        std::vector<sample> ref(grd.outputs_for(mid.size()));
        ref.resize(grd.process(mid.data(), mid.size(), ref.data()));

        // Whole, and in chunks of 1 to 97 frames: bit for bit.
        for (const std::size_t chunk :
             {std::size_t{3001}, std::size_t{1}, std::size_t{2}, std::size_t{7}, std::size_t{64}, std::size_t{97}}) {
            chain_t             c = fresh();
            std::vector<sample> y;
            for (std::size_t at = 0; at < x.size(); at += chunk) {
                const std::size_t   n = std::min(chunk, x.size() - at);
                std::vector<sample> part(c.outputs_for(n));
                ASSERT_EQ(c.process(x.data() + at, n, part.data()), part.size());
                y.insert(y.end(), part.begin(), part.end());
            }
            ASSERT_EQ(y.size(), ref.size());
            EXPECT_TRUE(y == ref) << "chunk " << chunk;
        }

        // outputs_for exact and frames_needed its inverse, from every phase.
        chain_t c = fresh();
        for (std::size_t phase = 0; phase < 4; ++phase) {
            for (std::size_t n = 0; n < 9; ++n) {
                chain_t             probe = c;
                std::vector<sample> out(probe.outputs_for(n) + 1);
                EXPECT_EQ(probe.process(x.data(), n, out.data()), c.outputs_for(n));
            }
            for (std::size_t k = 1; k < 6; ++k) {
                const std::size_t n = c.frames_needed(k);
                EXPECT_GE(c.outputs_for(n), k);
                EXPECT_LT(c.outputs_for(n - 1), k);
            }
            std::vector<sample> one(2);
            c.process(x.data() + phase, 1, one.data());
        }

        // Latency: the FIR's exact 20 of its outputs plus the 1:1 stage's own.
        chain_t lat = fresh();
        EXPECT_DOUBLE_EQ(lat.nominal_latency_output_frames(),
                         20.0 + lat.template stage<1>().nominal_latency_output_frames());

        // flush() equals zero padding, bit for bit, from every phase.
        for (std::size_t fed = 2995; fed <= 3001; ++fed) {
            chain_t             flushed = fresh();
            std::vector<sample> scratch(flushed.outputs_for(fed));
            flushed.process(x.data(), fed, scratch.data());
            chain_t             padded = flushed;
            const std::size_t   expect = flushed.flush_output_frames();
            std::vector<sample> a(expect + 4, sample{1});
            ASSERT_EQ(flushed.flush(a.data()), expect);
            const std::size_t zeros =
                2 * (padded.template stage<0>().window_frames() + 2 * padded.template stage<1>().window_frames());
            std::vector<sample> z(zeros, sample{});
            std::vector<sample> b(padded.outputs_for(zeros));
            const std::size_t   made = padded.process(z.data(), zeros, b.data());
            ASSERT_GT(made, expect);
            for (std::size_t i = 0; i < expect; ++i) {
                ASSERT_EQ(a[i], b[i]) << "fed " << fed << " output " << i;
            }
            for (std::size_t i = expect; i < made; ++i) {
                if constexpr (!std::is_floating_point_v<sample> || k_subnormal_control<sample>) {
                    ASSERT_EQ(b[i], sample{}) << "fed " << fed << " output " << i << ": silence past the flush";
                }
                else {
                    // No bit flushes this type here (double on the M4F and
                    // M33, everything on the soft-float M4): past the flush
                    // is the subnormal residue, the honest limit
                    // (iir_kernels.h; test_iir_subnormals.cpp).
                    ASSERT_LE(std::fabs(static_cast<double>(b[i])),
                              32.0 * static_cast<double>(std::numeric_limits<sample>::denorm_min()))
                        << "fed " << fed << " output " << i << ": a residue past the flush";
                }
            }
        }
    }

} // namespace
