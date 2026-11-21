// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// Locks down the sub-band round trip's contract (tap/dsp/subband.h) on
// test-local stages built on nyquist.h's designs, as chain.h was proven on
// decimate.h (SampleRateTap's docs/X4_X8_PLAN.md, 6.5): an op that changes
// nothing gives the input delayed by T bit for bit; removing the band leaves
// the stages' passband ripple, and only with the exact T (T - 1 and T + 1
// fail by 46 dB and more, which is what pins T against the signal path); T
// is the exact sum of the parts' delays, from whichever form each reports;
// any chunking and in-place running give the same bits; flush is zero
// padding and drains both paths; bands with different T compose by
// addition, each aligned at its low rate; an unmasked change re-injects the
// decimator's transition-band aliasing (X15's honest limit); a stage that
// reports only an exact delay path is accepted, one that reports neither is
// refused at compile time, a non-whole D tau_down + tau_up is refused at
// construction (and, where the types decide it, at compile time:
// tests/compile_fail/subband_non_integer_delay.cpp), and float tracks
// double. Typed over float and double; sized to stream a few thousand frames
// per case, so every case runs on every QEMU leg.
//
// The stages: an M-branch FIR decimator (output k as input kM arrives,
// decimate.h's alignment) and an L-phase FIR interpolator (y = h * the
// zero-stuffed input, L outputs per input), each over a nyquist.h design at
// 70 dB with kaiser_beta(70), accumulating through fir_kernels.h's dot_row.
// The m of each design is the smallest that meets 70 dB with a 1 dB margin
// at its passband edge, found by search_nyquist_m on a 16384-point grid
// (the 1024-point default agrees; measured on the host, 2026-11-21): m = 11
// at a passband of 3/8 of the stage's lower rate, for the 4th-band stage of
// the D = 4 round trip (N = 87) and the half-band stages that end the
// D = 8 chains (N = 43); m = 5 for the 4th-band stage that starts them,
// whose passband is 3/16 of its lower rate (N = 39). So T = 43 + 43 = 86 at
// D = 4 and (19 + 4 x 21) + (4 x 21 + 19) = 206 at D = 8 (D = 8 is X3's
// down by 4 then by 2, and its mirror).

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include "support/db.h"
#include "tap/dsp/chain.h"
#include "tap/dsp/decimate.h"
#include "tap/dsp/detail/expects.h"
#include "tap/dsp/fir_kernels.h"
#include "tap/dsp/kaiser.h"
#include "tap/dsp/nyquist.h"
#include "tap/dsp/subband.h"

namespace {

    using tap::dsp::chain;
    using tap::dsp::exact_delay_stage;
    using tap::dsp::exact_ratio;
    using tap::dsp::subband_op;
    using tap::dsp::subband_output;
    using tap::dsp::sync_stage;
    using tap::dsp::test::amplitude_db;

    // The scaffold's name, spelled twice in this file: here, and in
    // round_trip_type below, which names the class itself rather than this
    // alias (as test_chain.cpp's chainable names chain<T>).
    template <typename Down, typename Up, typename Op>
    using round_trip = tap::dsp::subband_delta<Down, Up, Op>;

    // ------------------------------------------------------------------
    // Test-local stages and op
    // ------------------------------------------------------------------

    /// Which delay member a test stage reports: an exact latency, an exact
    /// delay path (a near-linear recursive stage's form), or neither (a
    /// nonlinear-phase recursive stage's: a group delay that depends on f).
    enum class delay_form : std::uint8_t { latency, delay_path, neither };

    constexpr double      k_atten_db = 70.0;
    constexpr std::size_t k_m_last   = 11; ///< the stage at a passband of 3/8 of its lower rate
    constexpr std::size_t k_m_first  = 5;  ///< the 4th-band stage at 3/16, ahead of a half-band

    /// The L-th-band design at m taps per branch, sum L (nyquist.h).
    std::vector<double> nyquist_design(std::size_t l, std::size_t m) {
        std::vector<double> h(tap::dsp::nyquist_length(l, m));
        tap::dsp::design_nyquist(h, l, tap::dsp::kaiser_beta(k_atten_db));
        return h;
    }

    /// An M-branch FIR decimator over a prototype h at its input rate (DC
    /// gain sum(h)): y[k] = sum_t h[t] x[kM - t], x[< 0] = 0, output k written
    /// as input kM arrives. Group delay (N - 1) / 2 input samples, so
    /// (N - 1) / (2M) output frames.
    template <typename S, std::size_t Ch, std::size_t M, delay_form Form = delay_form::latency>
    class fir_down {
      public:
        using sample                        = S;
        static constexpr std::size_t k_up   = 1;
        static constexpr std::size_t k_down = M;
        static constexpr std::size_t channels() noexcept { return Ch; }

        explicit fir_down(const std::vector<double>& h)
            : m_taps(h.size())
            , m_row(h.size())
            , m_hist(2 * h.size() * Ch, S{0}) {
            for (std::size_t i = 0; i < m_taps; ++i) {
                m_row[i] = static_cast<S>(h[m_taps - 1 - i]); // oldest sample first
            }
        }

        std::size_t process(const S* in, std::size_t n, S* out) noexcept {
            std::size_t made = 0;
            for (std::size_t i = 0; i < n; ++i) {
                for (std::size_t c = 0; c < Ch; ++c) {
                    S* line              = m_hist.data() + c * 2 * m_taps;
                    line[m_pos]          = in[i * Ch + c];
                    line[m_pos + m_taps] = in[i * Ch + c];
                }
                if (m_phase == 0) {
                    for (std::size_t c = 0; c < Ch; ++c) {
                        const S* window    = m_hist.data() + c * 2 * m_taps + m_pos + 1;
                        out[made * Ch + c] = tap::dsp::dot_row<S>(m_row.data(), window, m_taps);
                    }
                    ++made;
                }
                m_pos   = m_pos + 1 == m_taps ? 0 : m_pos + 1;
                m_phase = (m_phase + 1) % M;
            }
            return made;
        }

        std::size_t outputs_for(std::size_t n) const noexcept {
            const std::size_t first = 1 + (M - m_phase) % M;
            return n < first ? 0 : 1 + (n - first) / M;
        }

        void reset() noexcept {
            std::fill(m_hist.begin(), m_hist.end(), S{0});
            m_pos   = 0;
            m_phase = 0;
        }

        std::size_t window_frames() const noexcept { return m_taps - 1; }

        exact_ratio latency_output_frames() const noexcept
            requires(Form == delay_form::latency)
        {
            return group_delay();
        }
        exact_ratio delay_path_output_frames() const noexcept
            requires(Form == delay_form::delay_path)
        {
            return group_delay();
        }
        double group_delay_output_frames(double /*f*/) const noexcept
            requires(Form == delay_form::neither)
        {
            return group_delay().value();
        }

      private:
        exact_ratio group_delay() const noexcept { return exact_ratio{m_taps - 1, 2 * M}; }

        std::size_t    m_taps;
        std::vector<S> m_row;
        std::vector<S> m_hist; ///< per channel: 2N, each sample written twice so a window is contiguous
        std::size_t    m_pos   = 0;
        std::size_t    m_phase = 0;
    };

    /// An L-phase FIR interpolator over a prototype h at its output rate:
    /// y[qL + p] = sum_k x[q - k] h[p + kL], x[< 0] = 0 (h convolved with the
    /// zero-stuffed input), L outputs per input. Group delay (N - 1) / 2
    /// output frames.
    template <typename S, std::size_t Ch, std::size_t L, delay_form Form = delay_form::latency>
    class fir_up {
      public:
        using sample                        = S;
        static constexpr std::size_t k_up   = L;
        static constexpr std::size_t k_down = 1;
        static constexpr std::size_t channels() noexcept { return Ch; }

        explicit fir_up(const std::vector<double>& h)
            : m_taps(h.size())
            , m_per_phase((h.size() + L - 1) / L)
            , m_rows(L * m_per_phase, S{0})
            , m_hist(2 * m_per_phase * Ch, S{0}) {
            for (std::size_t p = 0; p < L; ++p) {
                for (std::size_t i = 0; i < m_per_phase; ++i) {
                    const std::size_t t         = p + (m_per_phase - 1 - i) * L;
                    m_rows[p * m_per_phase + i] = t < m_taps ? static_cast<S>(h[t]) : S{0};
                }
            }
        }

        std::size_t process(const S* in, std::size_t n, S* out) noexcept {
            const std::size_t k = m_per_phase;
            for (std::size_t i = 0; i < n; ++i) {
                for (std::size_t c = 0; c < Ch; ++c) {
                    S* line         = m_hist.data() + c * 2 * k;
                    line[m_pos]     = in[i * Ch + c];
                    line[m_pos + k] = in[i * Ch + c];
                }
                for (std::size_t p = 0; p < L; ++p) {
                    for (std::size_t c = 0; c < Ch; ++c) {
                        const S* window           = m_hist.data() + c * 2 * k + m_pos + 1;
                        out[(i * L + p) * Ch + c] = tap::dsp::dot_row<S>(m_rows.data() + p * k, window, k);
                    }
                }
                m_pos = m_pos + 1 == k ? 0 : m_pos + 1;
            }
            return n * L;
        }

        std::size_t outputs_for(std::size_t n) const noexcept { return n * L; }

        void reset() noexcept {
            std::fill(m_hist.begin(), m_hist.end(), S{0});
            m_pos = 0;
        }

        std::size_t window_frames() const noexcept { return m_per_phase - 1; }

        exact_ratio latency_output_frames() const noexcept
            requires(Form == delay_form::latency)
        {
            return group_delay();
        }
        exact_ratio delay_path_output_frames() const noexcept
            requires(Form == delay_form::delay_path)
        {
            return group_delay();
        }
        double group_delay_output_frames(double /*f*/) const noexcept
            requires(Form == delay_form::neither)
        {
            return group_delay().value();
        }

      private:
        exact_ratio group_delay() const noexcept { return exact_ratio{m_taps - 1, 2}; }

        std::size_t    m_taps;
        std::size_t    m_per_phase;
        std::vector<S> m_rows; ///< phase p's taps, oldest input first
        std::vector<S> m_hist; ///< per channel: 2K, each sample written twice
        std::size_t    m_pos = 0;
    };

    /// The decimator over the L-th-band design, h / M (DC gain 1).
    template <typename S, std::size_t Ch, std::size_t M, delay_form Form = delay_form::latency>
    fir_down<S, Ch, M, Form> nyquist_down(std::size_t m) {
        std::vector<double> h = nyquist_design(M, m);
        for (double& v : h) {
            v /= static_cast<double>(M);
        }
        return fir_down<S, Ch, M, Form>(h);
    }

    /// The interpolator over the L-th-band design (each phase DC gain 1).
    template <typename S, std::size_t Ch, std::size_t L, delay_form Form = delay_form::latency>
    fir_up<S, Ch, L, Form> nyquist_up(std::size_t m) {
        return fir_up<S, Ch, L, Form>(nyquist_design(L, m));
    }

    /// The battery's op. Removing: the change is minus the band, `latency`
    /// frames late (change[k] = -in[k - latency]: the whole low band removed,
    /// with no mask). Otherwise the change is exact zeros (an identity op).
    template <typename S, std::size_t Ch>
    class band_op {
      public:
        using sample = S;
        static constexpr std::size_t channels() noexcept { return Ch; }

        explicit band_op(bool removing, std::size_t latency = 0)
            : m_removing(removing)
            , m_latency(latency)
            , m_line(latency * Ch, S{0}) {}

        std::size_t latency_frames() const noexcept { return m_latency; }
        std::size_t window_frames() const noexcept { return m_latency; }

        void process(const S* in, std::size_t n, S* change) noexcept {
            for (std::size_t i = 0; i < n * Ch; ++i) {
                S held = in[i];
                if (!m_line.empty()) {
                    held          = m_line[m_pos];
                    m_line[m_pos] = in[i];
                    m_pos         = m_pos + 1 == m_line.size() ? 0 : m_pos + 1;
                }
                change[i] = m_removing ? -held : S{0};
            }
        }

        void reset() noexcept {
            std::fill(m_line.begin(), m_line.end(), S{0});
            m_pos = 0;
        }

      private:
        bool           m_removing;
        std::size_t    m_latency;
        std::vector<S> m_line;
        std::size_t    m_pos = 0;
    };

    // The two round trips of the battery and their parts.
    template <typename S, std::size_t Ch = 1>
    using down_4 = fir_down<S, Ch, 4>;
    template <typename S, std::size_t Ch = 1>
    using up_4 = fir_up<S, Ch, 4>;
    template <typename S, std::size_t Ch = 1>
    using down_8 = chain<fir_down<S, Ch, 4>, fir_down<S, Ch, 2>>;
    template <typename S, std::size_t Ch = 1>
    using up_8 = chain<fir_up<S, Ch, 2>, fir_up<S, Ch, 4>>;
    template <typename S, std::size_t Ch = 1>
    using by_4 = round_trip<down_4<S, Ch>, up_4<S, Ch>, band_op<S, Ch>>;
    template <typename S, std::size_t Ch = 1>
    using by_8 = round_trip<down_8<S, Ch>, up_8<S, Ch>, band_op<S, Ch>>;

    template <typename S, std::size_t Ch = 1>
    down_4<S, Ch> make_down_4() {
        return nyquist_down<S, Ch, 4>(k_m_last);
    }
    template <typename S, std::size_t Ch = 1>
    up_4<S, Ch> make_up_4() {
        return nyquist_up<S, Ch, 4>(k_m_last);
    }
    template <typename S, std::size_t Ch = 1>
    down_8<S, Ch> make_down_8() {
        return down_8<S, Ch>(nyquist_down<S, Ch, 4>(k_m_first), nyquist_down<S, Ch, 2>(k_m_last));
    }
    template <typename S, std::size_t Ch = 1>
    up_8<S, Ch> make_up_8() {
        return up_8<S, Ch>(nyquist_up<S, Ch, 2>(k_m_last), nyquist_up<S, Ch, 4>(k_m_first));
    }

    template <typename S, std::size_t Ch = 1>
    by_4<S, Ch> make_by_4(band_op<S, Ch> op, std::size_t aligned = 0,
                          subband_output mode = subband_output::input_plus_change) {
        return by_4<S, Ch>(make_down_4<S, Ch>(), make_up_4<S, Ch>(), std::move(op), aligned, mode);
    }
    template <typename S, std::size_t Ch = 1>
    by_8<S, Ch> make_by_8(band_op<S, Ch> op, std::size_t aligned = 0,
                          subband_output mode = subband_output::input_plus_change) {
        return by_8<S, Ch>(make_down_8<S, Ch>(), make_up_8<S, Ch>(), std::move(op), aligned, mode);
    }

    // ------------------------------------------------------------------
    // Signals and measurement
    // ------------------------------------------------------------------

    constexpr double k_two_pi = 2.0 * std::numbers::pi;

    /// amp sin(2 pi f n), n in [0, frames), by a rotator (a recursion, so
    /// the emulated legs pay no libm per sample).
    std::vector<double> tone(std::size_t frames, double f, double amp) {
        std::vector<double> x(frames);
        const double        c  = std::cos(k_two_pi * f);
        const double        s  = std::sin(k_two_pi * f);
        double              re = amp;
        double              im = 0.0;
        for (double& v : x) {
            v               = im;
            const double nr = re * c - im * s;
            im              = re * s + im * c;
            re              = nr;
        }
        return x;
    }

    void add(std::vector<double>& x, const std::vector<double>& y) {
        for (std::size_t i = 0; i < x.size(); ++i) {
            x[i] += y[i];
        }
    }

    template <typename S>
    std::vector<S> as(const std::vector<double>& x) {
        std::vector<S> y(x.size());
        for (std::size_t i = 0; i < x.size(); ++i) {
            y[i] = static_cast<S>(x[i]);
        }
        return y;
    }

    /// The amplitude of the component at f (cycles per frame) of channel ch
    /// over [at, at + len): one bin of a rectangular DFT, exact for a tone
    /// with a whole number of cycles in len.
    template <typename S>
    double amplitude(const std::vector<S>& y, std::size_t at, std::size_t len, double f, std::size_t channels = 1,
                     std::size_t ch = 0) {
        const double c  = std::cos(k_two_pi * f);
        const double s  = -std::sin(k_two_pi * f);
        double       wr = 1.0;
        double       wi = 0.0;
        double       sr = 0.0;
        double       si = 0.0;
        for (std::size_t n = 0; n < len; ++n) {
            const double v = static_cast<double>(y[(at + n) * channels + ch]);
            sr += v * wr;
            si += v * wi;
            const double nr = wr * c - wi * s;
            wi              = wr * s + wi * c;
            wr              = nr;
        }
        return 2.0 * std::hypot(sr, si) / static_cast<double>(len);
    }

    /// Runs x through a round trip in chunks of `chunk` frames (whole when
    /// 0), in place when asked; returns what it wrote.
    template <typename R>
    std::vector<typename R::sample> run(R& rt, const std::vector<typename R::sample>& x, std::size_t chunk = 0,
                                        bool in_place = false) {
        using sample                 = typename R::sample;
        constexpr std::size_t w      = R::channels();
        const std::size_t     frames = x.size() / w;
        std::vector<sample>   y(x.size(), sample{0});
        if (in_place) {
            y = x;
        }
        const std::size_t step = chunk == 0 ? frames : chunk;
        for (std::size_t at = 0; at < frames; at += step) {
            const std::size_t n   = std::min(step, frames - at);
            const sample*     src = in_place ? y.data() + at * w : x.data() + at * w;
            EXPECT_EQ(rt.process(src, n, y.data() + at * w), n);
        }
        return y;
    }

    template <typename S>
    bool same_bits(const std::vector<S>& a, const std::vector<S>& b) {
        return a.size() == b.size() && (a.empty() || std::memcmp(a.data(), b.data(), a.size() * sizeof(S)) == 0);
    }

    /// x[n - d] + change[n], with the scaffold's select where the change is
    /// exactly zero: the output the contract promises, written by hand.
    template <typename S>
    std::vector<S> delayed_plus(const std::vector<S>& x, const std::vector<S>& change, std::size_t d) {
        std::vector<S> y(x.size());
        for (std::size_t n = 0; n < x.size(); ++n) {
            const S delayed = n >= d ? x[n - d] : S{0};
            y[n]            = change[n] == S{0} ? delayed : delayed + change[n];
        }
        return y;
    }

    template <typename S>
    class subband_test : public ::testing::Test {};
    using sample_types = ::testing::Types<float, double>;
    TYPED_TEST_SUITE(subband_test, sample_types, );

    // ------------------------------------------------------------------
    // The concepts
    // ------------------------------------------------------------------

    template <typename Down, typename Up, typename Op>
    concept round_trip_type = requires { typename tap::dsp::subband_delta<Down, Up, Op>; };

    static_assert(exact_delay_stage<down_4<double>>);
    static_assert(exact_delay_stage<down_8<float>>);
    static_assert(exact_delay_stage<tap::dsp::basic_decimator<double, 2>>); // latency_input_samples()
    static_assert(exact_delay_stage<fir_down<double, 1, 4, delay_form::delay_path>>);
    static_assert(subband_op<band_op<double, 1>>);
    static_assert(round_trip_type<down_4<double>, up_4<double>, band_op<double, 1>>);
    static_assert(round_trip_type<down_8<float, 2>, up_8<float, 2>, band_op<float, 2>>);
    // The round trip is itself a 1/1 sync_stage with an exact delay.
    static_assert(exact_delay_stage<by_4<double>>);
    static_assert(by_4<double>::k_factor == 4 && by_8<float>::k_factor == 8);

    // ------------------------------------------------------------------
    // The battery
    // ------------------------------------------------------------------

    // An op whose change is exactly zero: out[n] == x[n - T], every bit,
    // including negative zero, subnormals and the largest finite values, for
    // one stage and for chains, one channel and two, with an op latency and
    // an aligned latency. Structural: the zero-primed interpolator writes
    // exact zeros and the scaffold writes the delayed input where the change
    // is exactly zero.
    TYPED_TEST(subband_test, IdentityOpIsTheInputDelayedBitForBit) {
        using sample              = TypeParam;
        constexpr std::size_t len = 1500;
        std::vector<double>   xd  = tone(len, 0.0123, 0.4);
        add(xd, tone(len, 0.31, 0.2));
        std::vector<sample> x = as<sample>(xd);
        x[10]                 = -sample{0};
        x[11]                 = std::numeric_limits<sample>::denorm_min();
        x[12]                 = -std::numeric_limits<sample>::denorm_min();
        x[13]                 = std::numeric_limits<sample>::max();
        x[14]                 = std::numeric_limits<sample>::lowest();
        x[15]                 = std::numeric_limits<sample>::min() / sample{4}; // subnormal
        x[len - 1]            = -sample{0};

        const auto expect_delayed = [&x](const std::vector<sample>& y, std::size_t t) {
            std::vector<sample> want(x.size(), sample{0});
            for (std::size_t n = t; n < x.size(); ++n) {
                want[n] = x[n - t];
            }
            EXPECT_TRUE(same_bits(y, want)) << "T " << t;
        };

        auto by4 = make_by_4<sample>(band_op<sample, 1>(false));
        EXPECT_EQ(by4.latency(), 86u);
        expect_delayed(run(by4, x), by4.latency());

        auto by8 = make_by_8<sample>(band_op<sample, 1>(false));
        EXPECT_EQ(by8.latency(), 206u);
        expect_delayed(run(by8, x), by8.latency());

        auto late = make_by_4<sample>(band_op<sample, 1>(false, 7)); // T = 86 + 4 x 7
        EXPECT_EQ(late.latency(), 114u);
        expect_delayed(run(late, x), late.latency());

        auto aligned = make_by_4<sample>(band_op<sample, 1>(false), 301); // 86 + 4 x 53 + 3
        EXPECT_EQ(aligned.latency(), 301u);
        expect_delayed(run(aligned, x), aligned.latency());

        // Two channels: each channel its own stream, delayed alike.
        std::vector<sample> x2(2 * len);
        for (std::size_t n = 0; n < len; ++n) {
            x2[2 * n]     = x[n];
            x2[2 * n + 1] = n % 3 == 0 ? -sample{0} : -x[len - 1 - n];
        }
        auto                stereo = make_by_8<sample, 2>(band_op<sample, 2>(false));
        std::vector<sample> y2     = run(stereo, x2);
        std::vector<sample> want2(2 * len, sample{0});
        for (std::size_t n = stereo.latency(); n < len; ++n) {
            want2[2 * n]     = x2[2 * (n - stereo.latency())];
            want2[2 * n + 1] = x2[2 * (n - stereo.latency()) + 1];
        }
        EXPECT_TRUE(same_bits(y2, want2));
    }

    // The round trip's measured removal floors, in dB re the removed tone
    // (amplitude 0.5, on a bin of the 2000-frame window, after 512 frames
    // of settling): change = -band, so the output at f is 1 - H_down H_up
    // there. Measured on the host (x86-64, GCC 13.3, Release), 2026-11-21,
    // float and double within 0.01 dB of each other but at the near-null
    // (0.02 at D = 4: -116.06 double, -115.87 float):
    //
    //   D = 4 (one stage):  0.02 fs -116.06  0.05 fs  -74.77  0.09 fs  -65.17
    //   D = 4, L_op = 3:                     0.05 fs  -74.77
    //   D = 8 (two stages): 0.01 fs  -70.75  0.025 fs -62.69  0.029 fs -61.25
    //                                        0.045 fs -62.74
    //
    // The same figures follow from the designs' zero-phase responses,
    // |1 - A_down(f) A_up(f)|, to 0.01 dB, and a dense sweep of those over
    // the usable band [0, 3/8 r_low] puts the worst at -65.04 dB (D = 4,
    // 0.0905 fs) and -61.25 dB (D = 8, 0.0289 fs, the 0.029 case). With the
    // delay one sample off either way the residual is 2 sin(pi f): -18.02 /
    // -10.09 / -5.07 dB at D = 4 and -24.04 / -16.08 / -14.80 / -11.00 dB at
    // D = 8, at least 46 dB worse than with T at every frequency. The pins
    // are the measured values rounded up by 0.4 dB or more (6 dB at the
    // near-null), and down by 0.5 dB for the misaligned ones.
    struct removal_case {
        double f;          ///< cycles per full-rate frame
        double floor_db;   ///< the pinned ceiling on the residual with T
        double misaligned; ///< the pinned floor of the residual with T - 1 or T + 1
    };

    template <typename R, typename S>
    void expect_removal(R& with_t, R& change_only, std::size_t settle, std::size_t window, const removal_case& c,
                        double high_f) {
        std::vector<double> xd = tone(settle + window, c.f, 0.5);
        add(xd, tone(settle + window, high_f, 0.25));
        const std::vector<S> x      = as<S>(xd);
        const std::vector<S> y      = run(with_t, x);
        const std::vector<S> change = run(change_only, x);
        const std::size_t    t      = with_t.latency();
        // The scaffold's output is x[n - T] + its change, exactly.
        EXPECT_TRUE(same_bits(y, delayed_plus(x, change, t))) << "f " << c.f;
        for (const std::size_t d : {t - 1, t, t + 1}) {
            const std::vector<S> out      = d == t ? y : delayed_plus(x, change, d);
            const double         residual = amplitude_db(amplitude(out, settle, window, c.f) / 0.5);
            if (d == t) {
                EXPECT_LE(residual, c.floor_db) << "f " << c.f << " with T";
                // A tone above the band reaches the output unchanged: only the
                // round trip's stopband leaks into it (measured within 1e-6 dB).
                const double high = amplitude_db(amplitude(out, settle, window, high_f) / 0.25);
                EXPECT_LT(std::fabs(high), 1e-4) << "f " << high_f;
            }
            else {
                EXPECT_GE(residual, c.misaligned) << "f " << c.f << " with T " << (d < t ? "- 1" : "+ 1");
            }
        }
    }

    TYPED_TEST(subband_test, RemovingTheBandReachesTheStagesRippleFloor) {
        using sample                 = TypeParam;
        constexpr std::size_t settle = 512;
        constexpr std::size_t window = 2000;
        // D = 4 by one stage: the passband is 3/8 of fs / 4, 0.09375 fs.
        for (const removal_case c :
             {removal_case{0.02, -110.0, -18.5}, removal_case{0.05, -74.0, -10.6}, removal_case{0.09, -64.5, -5.6}}) {
            auto with_t      = make_by_4<sample>(band_op<sample, 1>(true));
            auto change_only = make_by_4<sample>(band_op<sample, 1>(true), 0, subband_output::change_only);
            expect_removal<by_4<sample>, sample>(with_t, change_only, settle, window, c, 0.3);
        }
        // An op latency of 3 low-rate frames moves T by 12 and the floor
        // stays: D L_op is in the delay the removal needs.
        {
            auto with_t      = make_by_4<sample>(band_op<sample, 1>(true, 3));
            auto change_only = make_by_4<sample>(band_op<sample, 1>(true, 3), 0, subband_output::change_only);
            EXPECT_EQ(with_t.latency(), 98u);
            expect_removal<by_4<sample>, sample>(with_t, change_only, settle, window, removal_case{0.05, -74.0, -10.6},
                                                 0.3);
        }
        // D = 8 by two stages: the passband is 3/8 of fs / 8, 0.046875 fs.
        for (const removal_case c : {removal_case{0.01, -70.0, -24.5}, removal_case{0.025, -62.0, -16.6},
                                     removal_case{0.029, -60.8, -15.3}, removal_case{0.045, -62.0, -11.5}}) {
            auto with_t      = make_by_8<sample>(band_op<sample, 1>(true));
            auto change_only = make_by_8<sample>(band_op<sample, 1>(true), 0, subband_output::change_only);
            expect_removal<by_8<sample>, sample>(with_t, change_only, settle, window, c, 0.3);
        }
    }

    // T = D tau_down + D L_op + tau_up, from each part's own report: one
    // stage, chains (their delays as chain.h composes them), decimate.h's
    // integer form, a delay path, mixed forms in one chain, an op latency;
    // and latency() = max(T, aligned_latency), the rest held at the low rate
    // and the full rate as the file header says.
    TYPED_TEST(subband_test, LatencyIsTheExactSumOfItsParts) {
        using sample = TypeParam;
        {
            const auto down = make_down_4<sample>();
            const auto up   = make_up_4<sample>();
            EXPECT_EQ(down.latency_output_frames(), (exact_ratio{43, 4})); // (87 - 1) / (2 x 4)
            EXPECT_EQ(up.latency_output_frames(), (exact_ratio{43, 1}));   // (87 - 1) / 2
            EXPECT_EQ(by_4<sample>::stages_delay_frames(down, up), (exact_ratio{86, 1}));
            const band_op<sample, 1> op(true, 5);
            EXPECT_EQ(by_4<sample>::round_trip_latency(down, up, op), 86u + 4u * 5u);
            const auto rt = make_by_4<sample>(band_op<sample, 1>(true, 5));
            EXPECT_EQ(rt.round_trip_latency(), 106u);
            EXPECT_EQ(rt.latency(), 106u);
            EXPECT_EQ(rt.latency_output_frames(), (exact_ratio{106, 1}));
            EXPECT_EQ(rt.extra_low_rate_frames(), 0u);
            EXPECT_EQ(rt.extra_full_rate_frames(), 0u);
            EXPECT_EQ(rt.lead(), 3u); // D - 1: output k by input 4k, four frames per low-rate frame
        }
        {
            const auto down = make_down_8<sample>();
            const auto up   = make_up_8<sample>();
            // chain.h's own composition: 19 at fs over 8, plus 21 at fs / 4 over 2.
            EXPECT_EQ(down.latency_output_frames(), (exact_ratio{103, 8}));
            EXPECT_EQ(up.latency_output_frames(), (exact_ratio{103, 1}));
            EXPECT_EQ(by_8<sample>::stages_delay_frames(down, up),
                      (exact_ratio{8, 1} * down.latency_output_frames() + up.latency_output_frames()));
            const auto rt = make_by_8<sample>(band_op<sample, 1>(false));
            EXPECT_EQ(rt.latency(), 206u);
            EXPECT_EQ(rt.lead(), 7u);
        }
        {
            // decimate.h's by-2 economy (81 taps: 40 input samples, its
            // latency_input_samples() form) and the half-band interpolator.
            using band_2 = round_trip<tap::dsp::basic_decimator<sample, 2>, fir_up<sample, 1, 2>, band_op<sample, 1>>;
            const band_2 rt(tap::dsp::basic_decimator<sample, 2>{}, nyquist_up<sample, 1, 2>(k_m_last),
                            band_op<sample, 1>(false));
            EXPECT_EQ(rt.down_stage().latency_input_samples(), 40u);
            EXPECT_EQ(rt.latency(), 40u + 21u);
        }
        {
            // Mixed forms in one chain each way: decimate.h's integer latency
            // and a delay path down (40 + 2 x 21), a latency and a delay path
            // up (2 x 21 + 13 for the m = 7 half-band, N = 27).
            using down_mixed =
                chain<tap::dsp::basic_decimator<sample, 2>, fir_down<sample, 1, 2, delay_form::delay_path>>;
            using up_mixed = chain<fir_up<sample, 1, 2>, fir_up<sample, 1, 2, delay_form::delay_path>>;
            static_assert(exact_delay_stage<down_mixed> && exact_delay_stage<up_mixed>);
            using mixed = round_trip<down_mixed, up_mixed, band_op<sample, 1>>;
            const mixed rt(
                down_mixed(tap::dsp::basic_decimator<sample, 2>{},
                           nyquist_down<sample, 1, 2, delay_form::delay_path>(k_m_last)),
                up_mixed(nyquist_up<sample, 1, 2>(k_m_last), nyquist_up<sample, 1, 2, delay_form::delay_path>(7)),
                band_op<sample, 1>(false, 2));
            EXPECT_EQ(rt.round_trip_latency(), 82u + 4u * 2u + 55u);
        }
        {
            // Aligned: below T it is T; above, the rest is split at the low
            // rate (53 frames) and the full rate (3).
            const auto below = make_by_4<sample>(band_op<sample, 1>(false), 50);
            EXPECT_EQ(below.latency(), 86u);
            const auto above = make_by_4<sample>(band_op<sample, 1>(false), 301);
            EXPECT_EQ(above.latency(), 301u);
            EXPECT_EQ(above.round_trip_latency(), 86u);
            EXPECT_EQ(above.extra_low_rate_frames(), 53u);
            EXPECT_EQ(above.extra_full_rate_frames(), 3u);
        }
    }

    // Any chunking, and running in place, write the same bits as one call:
    // a chain each way, an op latency, an aligned latency whose rest is held
    // at both rates, both output modes, one channel and two.
    TYPED_TEST(subband_test, ChunkingIsBitIdentical) {
        using sample                 = TypeParam;
        constexpr std::size_t frames = 2000;
        std::vector<double>   xd     = tone(frames, 0.0217, 0.4);
        add(xd, tone(frames, 0.163, 0.3));
        add(xd, tone(frames, 0.41, 0.2));
        const std::vector<sample> x = as<sample>(xd);

        // T = 206 + 8 x 2 = 222; aligned to 251: 3 frames at fs / 8 and 5 at fs.
        const auto fresh = [] { return make_by_8<sample>(band_op<sample, 1>(true, 2), 251); };
        auto       whole = fresh();
        ASSERT_EQ(whole.extra_low_rate_frames(), 3u);
        ASSERT_EQ(whole.extra_full_rate_frames(), 5u);
        const std::vector<sample> ref = run(whole, x);
        for (const std::size_t chunk : {1u, 2u, 3u, 7u, 8u, 63u, 255u, 256u, 257u, 1000u}) {
            auto rt = fresh();
            EXPECT_TRUE(same_bits(run(rt, x, chunk), ref)) << "chunk " << chunk;
        }
        {
            // Irregular chunks.
            auto                  rt = fresh();
            std::vector<sample>   y(frames);
            std::size_t           at        = 0;
            constexpr std::size_t pattern[] = {1, 5, 300, 2, 0, 17, 511};
            for (std::size_t k = 0; at < frames; ++k) {
                const std::size_t n = std::min(pattern[k % 7], frames - at);
                EXPECT_EQ(rt.process(x.data() + at, n, y.data() + at), n);
                at += n;
            }
            EXPECT_TRUE(same_bits(y, ref));
        }
        {
            auto rt = fresh();
            EXPECT_TRUE(same_bits(run(rt, x, 37, true), ref)) << "in place";
        }
        {
            const auto change_fresh = [] {
                return make_by_8<sample>(band_op<sample, 1>(true, 2), 251, subband_output::change_only);
            };
            auto                      a      = change_fresh();
            auto                      b      = change_fresh();
            const std::vector<sample> change = run(a, x);
            EXPECT_TRUE(same_bits(run(b, x, 3), change));
            // and the two modes agree: x[n - A] + change is the scaffold's output
            EXPECT_TRUE(same_bits(delayed_plus(x, change, 251), ref));
        }
        {
            std::vector<sample> x2(2 * frames);
            for (std::size_t n = 0; n < frames; ++n) {
                x2[2 * n]     = x[n];
                x2[2 * n + 1] = x[frames - 1 - n];
            }
            auto                      whole2 = make_by_4<sample, 2>(band_op<sample, 2>(true, 1), 97);
            const std::vector<sample> ref2   = run(whole2, x2);
            for (const std::size_t chunk : {1u, 255u, 256u, 257u}) {
                auto rt = make_by_4<sample, 2>(band_op<sample, 2>(true, 1), 97);
                EXPECT_TRUE(same_bits(run(rt, x2, chunk), ref2)) << "two channels, chunk " << chunk;
            }
            // Each channel is the mono round trip of its own stream.
            auto                mono = make_by_4<sample>(band_op<sample, 1>(true, 1), 97);
            std::vector<sample> left(frames);
            for (std::size_t n = 0; n < frames; ++n) {
                left[n] = ref2[2 * n];
            }
            EXPECT_TRUE(same_bits(left, run(mono, x)));
        }
    }

    // flush() writes window_frames() frames, equal to zero padding bit for
    // bit, and the padded stream is silence after them; the delay path's
    // last input comes out where T puts it and the change path's tail runs
    // past it, so both are drained.
    template <typename R>
    void expect_flush_drains(R rt, const std::vector<typename R::sample>& x) {
        using sample            = typename R::sample;
        constexpr std::size_t w = R::channels();
        std::vector<sample>   scratch(x.size());
        rt.process(x.data(), x.size() / w, scratch.data());
        R padded = rt;

        const std::size_t   flushed = rt.flush_output_frames();
        std::vector<sample> a((flushed + 4) * w, sample{1});
        ASSERT_EQ(rt.flush(a.data()), flushed);
        for (std::size_t i = flushed * w; i < a.size(); ++i) {
            ASSERT_EQ(a[i], sample{1}) << "flush wrote past its count";
        }
        const std::vector<sample> zeros(3 * flushed * w, sample{0});
        std::vector<sample>       b(zeros.size());
        padded.process(zeros.data(), 3 * flushed, b.data());
        std::size_t last_nonzero = 0;
        for (std::size_t i = 0; i < b.size(); ++i) {
            if (i < flushed * w) {
                ASSERT_EQ(std::memcmp(&a[i], &b[i], sizeof(sample)), 0) << "frame " << i / w;
            }
            else {
                ASSERT_EQ(b[i], sample{0}) << "frame " << i / w << " after the window";
            }
            if (b[i] != sample{0}) {
                last_nonzero = i / w;
            }
        }
        if (rt.output_mode() == subband_output::input_plus_change) {
            EXPECT_GE(last_nonzero + 1, rt.latency()) << "the delay path's tail";
        }
        EXPECT_LT(last_nonzero, flushed);
    }

    TYPED_TEST(subband_test, FlushDrainsBothPaths) {
        using sample                 = TypeParam;
        constexpr std::size_t frames = 700;
        std::vector<double>   xd     = tone(frames, 0.03, 0.5);
        add(xd, tone(frames, 0.27, 0.25));
        const std::vector<sample> x = as<sample>(xd);

        // The bound, by hand for D = 4: the decimator's window (86) plus the
        // inputs it needs to write the interpolator's window (21 frames, at
        // most 4 x 21 + 4), plus the lead (3): 177, above T = 86.
        auto removing = make_by_4<sample>(band_op<sample, 1>(true));
        EXPECT_EQ(removing.window_frames(), 86u + (4u * 21u + 4u) + 3u);
        EXPECT_EQ(removing.flush_output_frames(), 177u);
        expect_flush_drains(removing, x);

        // The change path's tail outlasts the delay path's: the flushed
        // output is nonzero after its first T frames.
        {
            auto                rt = make_by_4<sample>(band_op<sample, 1>(true));
            std::vector<sample> scratch(frames);
            rt.process(x.data(), frames, scratch.data());
            std::vector<sample> tail(rt.flush_output_frames());
            rt.flush(tail.data());
            bool after_t = false;
            for (std::size_t i = rt.latency(); i < tail.size(); ++i) {
                after_t = after_t || tail[i] != sample{0};
            }
            EXPECT_TRUE(after_t);
        }

        // An identity op: the flush is the input's last T frames, then silence.
        {
            auto                rt = make_by_4<sample>(band_op<sample, 1>(false));
            std::vector<sample> scratch(frames);
            rt.process(x.data(), frames, scratch.data());
            std::vector<sample> tail(rt.flush_output_frames());
            ASSERT_EQ(rt.flush(tail.data()), 177u);
            for (std::size_t i = 0; i < tail.size(); ++i) {
                const sample want = i < 86 ? x[frames - 86 + i] : sample{0};
                ASSERT_EQ(std::memcmp(&tail[i], &want, sizeof(sample)), 0) << "frame " << i;
            }
            // reset() returns to the fresh state.
            rt.reset();
            auto fresh = make_by_4<sample>(band_op<sample, 1>(false));
            EXPECT_TRUE(same_bits(run(rt, x), run(fresh, x)));
        }

        // Chains, an op latency, both rests of an aligned latency, both modes,
        // and two channels.
        expect_flush_drains(make_by_8<sample>(band_op<sample, 1>(true, 2), 251), x);
        expect_flush_drains(make_by_8<sample>(band_op<sample, 1>(true, 2), 251, subband_output::change_only), x);
        std::vector<sample> x2(2 * frames);
        for (std::size_t n = 0; n < frames; ++n) {
            x2[2 * n]     = x[n];
            x2[2 * n + 1] = -x[n];
        }
        expect_flush_drains(make_by_4<sample, 2>(band_op<sample, 2>(true, 1), 97), x2);
    }

    // Bands compose by addition (6.3): out = x[n - A] + the sum of the
    // bands' changes, every band aligned to A, the largest T. Three bands
    // whose T differ: D = 8 (206), D = 4 (86: 120 more, 30 frames at fs / 4)
    // and D = 2 with decimate.h's decimator (61: 145 more, 72 frames at
    // fs / 2 and 1 at fs). The D = 8 band holds the delay line; the others
    // run change-only.
    template <typename S>
    using band_2 = round_trip<tap::dsp::basic_decimator<S, 2>, fir_up<S, 1, 2>, band_op<S, 1>>;

    template <typename S>
    band_2<S> make_band_2(bool removing, std::size_t aligned, subband_output mode) {
        return band_2<S>(tap::dsp::basic_decimator<S, 2>{}, nyquist_up<S, 1, 2>(k_m_last), band_op<S, 1>(removing),
                         aligned, mode);
    }

    template <typename S>
    std::vector<S> sum_bands(std::vector<S> out, const std::vector<S>& c4, const std::vector<S>& c2) {
        for (std::size_t n = 0; n < out.size(); ++n) {
            out[n] = c4[n] == S{0} ? out[n] : out[n] + c4[n];
            out[n] = c2[n] == S{0} ? out[n] : out[n] + c2[n];
        }
        return out;
    }

    TYPED_TEST(subband_test, BandsComposeByAddition) {
        using sample                 = TypeParam;
        constexpr std::size_t settle = 512;
        constexpr std::size_t window = 2000;
        std::vector<double>   xd     = tone(settle + window, 0.02, 0.5);
        add(xd, tone(settle + window, 0.3, 0.25));
        const std::vector<sample> x = as<sample>(xd);

        // The common latency, from the parts before construction.
        const band_op<sample, 1> idle(false);
        const std::size_t        a =
            std::max({by_8<sample>::round_trip_latency(make_down_8<sample>(), make_up_8<sample>(), idle),
                      by_4<sample>::round_trip_latency(make_down_4<sample>(), make_up_4<sample>(), idle),
                      band_2<sample>::round_trip_latency(tap::dsp::basic_decimator<sample, 2>{},
                                                         nyquist_up<sample, 1, 2>(k_m_last), idle)});
        ASSERT_EQ(a, 206u);

        // Each change-only band at A is its own change delayed by A - T,
        // bit for bit: the alignment, held at the low rate first.
        {
            auto at_a   = make_by_4<sample>(band_op<sample, 1>(true), a, subband_output::change_only);
            auto at_own = make_by_4<sample>(band_op<sample, 1>(true), 0, subband_output::change_only);
            EXPECT_EQ(at_a.latency(), a);
            EXPECT_EQ(at_a.extra_low_rate_frames(), 30u);
            EXPECT_EQ(at_a.extra_full_rate_frames(), 0u);
            const std::vector<sample> shifted =
                delayed_plus(run(at_own, x), std::vector<sample>(x.size(), sample{0}), a - 86);
            EXPECT_TRUE(same_bits(run(at_a, x), shifted));
        }
        {
            auto at_a   = make_band_2<sample>(true, a, subband_output::change_only);
            auto at_own = make_band_2<sample>(true, 0, subband_output::change_only);
            EXPECT_EQ(at_own.latency(), 61u);
            EXPECT_EQ(at_a.extra_low_rate_frames(), 72u);
            EXPECT_EQ(at_a.extra_full_rate_frames(), 1u);
            const std::vector<sample> shifted =
                delayed_plus(run(at_own, x), std::vector<sample>(x.size(), sample{0}), a - 61);
            EXPECT_TRUE(same_bits(run(at_a, x), shifted));
        }

        // The composition, one band removing at a time: its tone goes to that
        // band's floor at A (one sample off, it would stay at -18 dB), and
        // the tone above every band is untouched. Measured (host, float and
        // double within 0.2 dB): D = 8 -61.89, D = 4 -116.06 (-115.87 float),
        // D = 2 -86.99 dB, each that band's own removal floor at 0.02 fs (the
        // D = 2 band's worst over its usable band, by its responses, is
        // -69.91 dB); pinned 0.4 dB up, 6 at the near-null. With no band
        // removing, the composition is the input delayed by A, bit for bit.
        const double floors[] = {-61.4, -110.0, -86.5};
        for (int removing = -1; removing < 3; ++removing) {
            auto band8 = make_by_8<sample>(band_op<sample, 1>(removing == 0), a);
            auto band4 = make_by_4<sample>(band_op<sample, 1>(removing == 1), a, subband_output::change_only);
            auto band2 = make_band_2<sample>(removing == 2, a, subband_output::change_only);
            std::vector<sample> c4  = run(band4, x);
            std::vector<sample> c2  = run(band2, x);
            std::vector<sample> out = sum_bands(run(band8, x), c4, c2);
            if (removing < 0) {
                EXPECT_TRUE(same_bits(out, delayed_plus(x, std::vector<sample>(x.size(), sample{0}), a)));
                continue;
            }
            const double residual = amplitude_db(amplitude(out, settle, window, 0.02) / 0.5);
            EXPECT_LE(residual, floors[removing]) << "band " << removing;
            const double high = amplitude_db(amplitude(out, settle, window, 0.3) / 0.25);
            EXPECT_LT(std::fabs(high), 1e-4) << "band " << removing;
        }
    }

    // X15's honest limit, pinned: an op that changes the bins above p r_low
    // (here, all of the low band: no mask) puts the decimator's
    // transition-band aliasing into its change, and the interpolator writes
    // it back. A tone at 0.55 r_low (0.1375 fs at D = 4) folds to
    // 0.45 r_low (0.1125 fs), where the input has nothing (below -170 dB,
    // the float tone's own rounding); the output there is measured at
    // -21.40 dB re the tone (host, float and double alike), pinned within
    // 1 dB, where a change confined to the bins up to p r_low would leave
    // at most the stages' stopband, 70 dB down. (The tone itself comes out
    // 0.08 dB down: the decimator passes little of it into the change.)
    TYPED_TEST(subband_test, UnmaskedChangeIsReinjected) {
        using sample                     = TypeParam;
        constexpr std::size_t     settle = 512;
        constexpr std::size_t     window = 2000;
        const std::vector<sample> x      = as<sample>(tone(settle + window, 0.1375, 0.5));
        auto                      rt     = make_by_4<sample>(band_op<sample, 1>(true));
        const std::vector<sample> y      = run(rt, x);
        const double              input  = amplitude_db(amplitude(x, settle, window, 0.1125) / 0.5 + 1e-300);
        const double              alias  = amplitude_db(amplitude(y, settle, window, 0.1125) / 0.5);
        EXPECT_LT(input, -150.0);
        EXPECT_GE(alias, -22.4);
        EXPECT_LE(alias, -20.4);
    }

    // A stage that reports only an exact delay path (the near-linear
    // recursive stage's member) enters the scaffold like one that reports a
    // latency: the same T, and the same bits for the same filters.
    TYPED_TEST(subband_test, DelayPathStageIsAccepted) {
        using sample  = TypeParam;
        using path_rt = round_trip<fir_down<sample, 1, 4, delay_form::delay_path>,
                                   fir_up<sample, 1, 4, delay_form::delay_path>, band_op<sample, 1>>;
        static_assert(!tap::dsp::detail::reports_latency_output<fir_down<sample, 1, 4, delay_form::delay_path>>);
        static_assert(tap::dsp::reports_delay_path<fir_down<sample, 1, 4, delay_form::delay_path>>);
        path_rt rt(nyquist_down<sample, 1, 4, delay_form::delay_path>(k_m_last),
                   nyquist_up<sample, 1, 4, delay_form::delay_path>(k_m_last), band_op<sample, 1>(true));
        auto    reference = make_by_4<sample>(band_op<sample, 1>(true));
        EXPECT_EQ(rt.latency(), 86u);
        std::vector<double> xd = tone(1200, 0.04, 0.5);
        add(xd, tone(1200, 0.2, 0.3));
        const std::vector<sample> x = as<sample>(xd);
        EXPECT_TRUE(same_bits(run(rt, x), run(reference, x)));
    }

    // float tracks double: the D = 8 removal of a two-tone input, sample for
    // sample. Measured (host, GCC 13.3): max |float - double| 6.5e-8 over
    // 2512 frames of a 0.75-peak input, the float format's rounding (its
    // ulp at 0.5 is 6e-8); pinned at 1e-7.
    TEST(Subband, FloatTracksDouble) {
        constexpr std::size_t frames = 2512;
        std::vector<double>   xd     = tone(frames, 0.025, 0.5);
        add(xd, tone(frames, 0.3, 0.25));
        auto                      rt_d  = make_by_8<double>(band_op<double, 1>(true));
        auto                      rt_f  = make_by_8<float>(band_op<float, 1>(true));
        const std::vector<double> yd    = run(rt_d, xd);
        const std::vector<float>  yf    = run(rt_f, as<float>(xd));
        double                    worst = 0.0;
        for (std::size_t n = 0; n < frames; ++n) {
            worst = std::max(worst, std::fabs(static_cast<double>(yf[n]) - yd[n]));
        }
        EXPECT_LE(worst, 1e-7);
    }

    // A non-whole D tau_down + tau_up is refused. Here the decimator's tap
    // count is construction geometry: an even-length Kaiser lowpass delays
    // by half a sample at its input rate, (42 - 1) / 2, so with the
    // half-band interpolator's 21 the stages delay 83/2 full-rate samples
    // and the constructor rejects them (raise<std::invalid_argument>, or a
    // terminate without exceptions, where supports() is the check to make
    // first); the same lowpass at 43 taps is accepted. Where the types
    // decide it (static constexpr delays), a non-whole sum is a compile
    // error: tests/compile_fail/subband_non_integer_delay.cpp, run as
    // subband_compile_fail.NonIntegerDelayIsRejected.
    std::vector<double> kaiser_lowpass(std::size_t taps) {
        std::vector<double> h(taps);
        tap::dsp::design_prototype(h, 1, 0.5, tap::dsp::kaiser_beta(k_atten_db));
        return h;
    }

    using fractional_rt = round_trip<fir_down<double, 1, 2>, fir_up<double, 1, 2>, band_op<double, 1>>;

    TEST(Subband, NonIntegerDelayIsRejected) {
        const fir_down<double, 1, 2> even(kaiser_lowpass(42));
        const fir_down<double, 1, 2> odd(kaiser_lowpass(43));
        const fir_up<double, 1, 2>   up = nyquist_up<double, 1, 2>(k_m_last);
        EXPECT_EQ(fractional_rt::stages_delay_frames(even, up), (exact_ratio{83, 2}));
        EXPECT_FALSE(fractional_rt::supports(even, up));
        EXPECT_EQ(fractional_rt::stages_delay_frames(odd, up), (exact_ratio{42, 1}));
        EXPECT_TRUE(fractional_rt::supports(odd, up));
        const fractional_rt accepted(odd, up, band_op<double, 1>(false));
        EXPECT_EQ(accepted.latency(), 42u);
#if TAP_DSP_HAS_EXCEPTIONS
        EXPECT_THROW(fractional_rt(even, up, band_op<double, 1>(false)), std::invalid_argument);
#endif // without exceptions the rejection terminates: SubbandDeathTest below, where gtest has death tests
    }

#if !TAP_DSP_HAS_EXCEPTIONS && defined(GTEST_HAS_DEATH_TEST) && GTEST_HAS_DEATH_TEST
    TEST(SubbandDeathTest, NonIntegerDelayTerminatesWithoutExceptions) {
        const fir_down<double, 1, 2> even(kaiser_lowpass(42));
        const fir_up<double, 1, 2>   up = nyquist_up<double, 1, 2>(k_m_last);
        EXPECT_DEATH(fractional_rt(even, up, band_op<double, 1>(false)), "");
    }
#endif

    // Structure-only stages whose delays are compile-time constants (static
    // constexpr latency_output_frames()), so the types decide whether T is
    // whole: keep_down keeps every M-th frame (output k as input kM arrives,
    // or, Late, as input kM + M - 1 does), hold_up repeats each input L
    // times. They filter nothing; the scaffold's structure does not care.
    template <std::size_t M, std::uint64_t Num, std::uint64_t Den, bool Late = false>
    class keep_down {
      public:
        using sample                        = double;
        static constexpr std::size_t k_up   = 1;
        static constexpr std::size_t k_down = M;
        static constexpr std::size_t channels() noexcept { return 1; }
        static constexpr exact_ratio latency_output_frames() noexcept { return exact_ratio{Num, Den}; }
        std::size_t                  window_frames() const noexcept { return 0; }
        std::size_t                  process(const double* in, std::size_t n, double* out) noexcept {
            std::size_t made = 0;
            for (std::size_t i = 0; i < n; ++i) {
                if (m_phase == (Late ? M - 1 : 0)) {
                    out[made++] = in[i];
                }
                m_phase = (m_phase + 1) % M;
            }
            return made;
        }
        std::size_t outputs_for(std::size_t n) const noexcept {
            std::size_t made  = 0;
            std::size_t phase = m_phase;
            for (std::size_t i = 0; i < n; ++i) {
                made += phase == (Late ? M - 1 : 0) ? 1 : 0;
                phase = (phase + 1) % M;
            }
            return made;
        }
        void reset() noexcept { m_phase = 0; }

      private:
        std::size_t m_phase = 0;
    };

    template <std::size_t L, std::uint64_t Num, std::uint64_t Den>
    class hold_up {
      public:
        using sample                        = double;
        static constexpr std::size_t k_up   = L;
        static constexpr std::size_t k_down = 1;
        static constexpr std::size_t channels() noexcept { return 1; }
        static constexpr exact_ratio latency_output_frames() noexcept { return exact_ratio{Num, Den}; }
        std::size_t                  window_frames() const noexcept { return 0; }
        std::size_t                  process(const double* in, std::size_t n, double* out) noexcept {
            for (std::size_t i = 0; i < n * L; ++i) {
                out[i] = in[i / L];
            }
            return n * L;
        }
        std::size_t outputs_for(std::size_t n) const noexcept { return n * L; }
        void        reset() noexcept {}
    };

    // 2 x 1/4 + 1/2: each part a fraction, the sum whole, decided by the
    // types and accepted. (A non-whole static sum does not compile.)
    TEST(Subband, TheTypesDecideAStaticDelay) {
        using static_rt = round_trip<keep_down<2, 1, 4>, hold_up<2, 1, 2>, band_op<double, 1>>;
        static_assert(tap::dsp::detail::static_delay<keep_down<2, 1, 4>>::k_known);
        static_assert(!tap::dsp::detail::static_delay<down_4<double>>::k_known);
        static_assert((exact_ratio{2, 1} * keep_down<2, 1, 4>::latency_output_frames()
                       + hold_up<2, 1, 2>::latency_output_frames())
                      == exact_ratio{1, 1});
        static_rt                 rt(keep_down<2, 1, 4>{}, hold_up<2, 1, 2>{}, band_op<double, 1>(false));
        const std::vector<double> x = tone(40, 0.1, 1.0);
        EXPECT_EQ(rt.latency(), 1u);
        EXPECT_TRUE(same_bits(run(rt, x), delayed_plus(x, std::vector<double>(x.size(), 0.0), 1)));
    }

    // A change path that falls behind its input cannot be summed frame by
    // frame: a decimator that writes output k at input kM + M - 1 has
    // written nothing after the first input. Refused at construction.
    TEST(Subband, ALaggingChangePathIsRejected) {
        using lagging_rt = round_trip<keep_down<2, 1, 2, true>, hold_up<2, 0, 1>, band_op<double, 1>>;
        EXPECT_FALSE(lagging_rt::supports(keep_down<2, 1, 2, true>{}, hold_up<2, 0, 1>{}));
        EXPECT_TRUE((round_trip<keep_down<2, 1, 2>, hold_up<2, 0, 1>, band_op<double, 1>>::supports(
            keep_down<2, 1, 2>{}, hold_up<2, 0, 1>{})));
#if TAP_DSP_HAS_EXCEPTIONS
        EXPECT_THROW(lagging_rt(keep_down<2, 1, 2, true>{}, hold_up<2, 0, 1>{}, band_op<double, 1>(false)),
                     std::invalid_argument);
#endif // without exceptions: SubbandDeathTest below, where gtest has death tests
    }

#if !TAP_DSP_HAS_EXCEPTIONS && defined(GTEST_HAS_DEATH_TEST) && GTEST_HAS_DEATH_TEST
    TEST(SubbandDeathTest, ALaggingChangePathTerminatesWithoutExceptions) {
        using lagging_rt = round_trip<keep_down<2, 1, 2, true>, hold_up<2, 0, 1>, band_op<double, 1>>;
        EXPECT_DEATH(lagging_rt(keep_down<2, 1, 2, true>{}, hold_up<2, 0, 1>{}, band_op<double, 1>(false)), "");
    }
#endif

    // A stage that reports neither an exact latency nor an exact delay path
    // (the nonlinear-phase recursive stage: a group delay that depends on
    // f) is refused at compile time, alone or inside a chain, either way;
    // concept negatives, so the refusal is pinned without a compile-fail
    // build.
    using neither_down = fir_down<double, 1, 4, delay_form::neither>;
    using neither_up   = fir_up<double, 1, 4, delay_form::neither>;
    static_assert(sync_stage<neither_down> && sync_stage<neither_up>);
    static_assert(!exact_delay_stage<neither_down> && !exact_delay_stage<neither_up>);
    static_assert(!round_trip_type<neither_down, up_4<double>, band_op<double, 1>>);
    static_assert(!round_trip_type<down_4<double>, neither_up, band_op<double, 1>>);
    static_assert(!exact_delay_stage<chain<neither_down, fir_down<double, 1, 2>>>);
    static_assert(!round_trip_type<chain<neither_down, fir_down<double, 1, 2>>, up_8<double>, band_op<double, 1>>);
    static_assert(!round_trip_type<down_8<double>, chain<fir_up<double, 1, 2>, neither_up>, band_op<double, 1>>);
    // An op without its latency is not an op.
    struct no_latency_op {
        using sample = double;
        static constexpr std::size_t channels() noexcept { return 1; }
        void                         process(const double* /*in*/, std::size_t /*n*/, double* /*change*/) noexcept {}
        void                         reset() noexcept {}
    };
    static_assert(!subband_op<no_latency_op>);
    static_assert(!round_trip_type<down_4<double>, up_4<double>, no_latency_op>);

    TEST(Subband, AStageWithNeitherIsRefused) {
        // The static_asserts above are the test; this records it in the run.
        EXPECT_FALSE((exact_delay_stage<neither_down>));
        EXPECT_TRUE((exact_delay_stage<fir_down<double, 1, 4, delay_form::delay_path>>));
    }

} // namespace
