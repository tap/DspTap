/// @file psola.h
/// @brief Pitch-synchronous overlap-add (PSOLA) pitch shifter with a fixed numeric contract.
// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// The real-time TD-PSOLA resynthesis stage from the standard speech-processing
// literature: Hann-windowed grains two source periods long are extracted at
// period-spaced analysis marks and overlap-added at period/ratio-spaced
// synthesis marks, scaled by 1/ratio so the window sum stays unity. The caller
// supplies the period (from tap::dsp::yin or any other tracker) — this class is
// deliberately detection-agnostic so the two concerns stay independently
// testable and reusable (pitch correction, formant work, embedded targets).
//
// Practical notes: analysis marks are a free-running period-synchronous
// scheduler, not glottal-epoch estimates — the standard real-time
// simplification. Synthesis marks are placed with sub-sample precision (the
// grain is resampled through 4-point Hermite interpolation), which keeps the
// period jitter of integer-rounded marks out of the output.
//
// Know what PSOLA is: it resamples the source's SPECTRAL ENVELOPE at the new
// harmonic spacing — which is exactly why it preserves formants on voice, and
// exactly why a PURE TONE far from any new harmonic thins toward silence
// (e.g. a sine shifted up an octave has no harmonic left at the envelope's
// only peak). Feed it harmonic-rich, voice-like material; for pure tones use
// a waveform-preserving shifter (the two-tap engine, or tap::dsp::pvoc). The
// test battery pins both behaviors on purpose.

#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <type_traits>
#include <vector>

namespace tap::dsp {

    /// TD-PSOLA pitch shifter, parameterized over the sample type. Both
    /// instantiations run the identical algorithm in their own precision:
    /// double is the desktop/golden profile, float the embedded profile.
    ///
    /// Geometry is fixed at construction (the deepest period the caller will
    /// ever supply) and every buffer is allocated there; process() is noexcept
    /// and allocation-free, safe on a real-time audio thread.
    ///
    /// Contract points:
    ///   - process(x, period, ratio) consumes one input sample and produces one
    ///     output sample delayed by exactly latency() samples.
    ///   - period is the CURRENT source period in samples (fractional ok),
    ///     clamped to [k_min_period, max_period]. Pass the last known value
    ///     while the source is unpitched — the scheduler keeps running.
    ///   - ratio is the pitch ratio (2 = up an octave), clamped to [1/4, 4].
    ///   - At ratio == 1 the synthesis marks coincide with the analysis
    ///     spacing, the Hann windows sum to exactly one, and the output is the
    ///     input delayed by latency() (plus interpolation error) — pinned by
    ///     the test battery.
    ///   - Run time is unbounded: the sample clock is a fixed-width 32-bit
    ///     count bounded below 2 * clock_wrap() (a multiple of the ring size),
    ///     so no counter overflows on any target. The contract is that bound,
    ///     not an observable. The wrap moves the fractional mark positions in
    ///     magnitude only, so from 2 * clock_wrap() samples on (at most 2^19,
    ///     10.9 s at 48 kHz; 519,552 samples for max_period 900) the output
    ///     diverges from an unbounded-clock reference at the ~1e-8 relative
    ///     level (measured 3.8e-9 absolute in double and one float ulp in
    ///     float, against a 0.35 peak). Pinned by the test battery.
    template <typename Sample>
    class basic_psola {
        static_assert(std::is_same_v<Sample, float> || std::is_same_v<Sample, double>,
                      "basic_psola supports the two Tap numeric profiles: float and double");

      public:
        static constexpr Sample k_min_period = Sample(8);
        static constexpr Sample k_min_ratio  = Sample(0.25);
        static constexpr Sample k_max_ratio  = Sample(4);

        /// @pre 16 <= max_period < 2^26 — the deepest period process() will be
        /// given; the upper bound keeps 2 * clock_wrap() and every grain position
        /// inside the int32 sample clock.
        explicit basic_psola(size_t max_period)
            : m_max_period(static_cast<Sample>(max_period))
            , m_latency(2 * max_period + 2) {
            assert(max_period >= 16);
            assert(max_period < (size_t{1} << 26)); // see @pre
            // One ring size serves both buffers; the clock wrap relies on that.
            // Input history: a grain reaches back to (mark - period) and marks lag the
            // input cursor by up to two periods -> three periods of history plus slack.
            // Output accumulator: emission lags by latency(); grains extend up to one
            // period past their mark -> latency + period ahead of the emit cursor.
            const size_t ring = 4 * max_period + 8;
            m_input.assign(ring, Sample(0));
            m_accum.assign(ring, Sample(0));
            m_wrap = static_cast<std::int32_t>(ring * std::max<size_t>(1, static_cast<size_t>(k_clock_span) / ring));
            clear();
        }

        /// Emission delay of the shifter, in samples (fixed at construction).
        size_t latency() const noexcept { return m_latency; }

        size_t max_period() const noexcept { return static_cast<size_t>(m_max_period); }

        /// Period of the sample clock's wrap, in samples: the largest multiple of
        /// the ring size not above 2^18 (or one ring, if the ring is larger). The
        /// clock first wraps at 2 * clock_wrap() and every clock_wrap() after.
        size_t clock_wrap() const noexcept { return static_cast<size_t>(m_wrap); }

        /// Zero all running state (buffers, marks, counters).
        void clear() noexcept {
            std::fill(m_input.begin(), m_input.end(), Sample(0));
            std::fill(m_accum.begin(), m_accum.end(), Sample(0));
            m_n          = 0;
            m_next_mark  = 0.0;
            m_prev_mark  = 0.0;
            m_have_mark  = false;
            m_next_synth = static_cast<double>(m_latency); // first grain lands at the first emitted sample
        }

        /// Consume one input sample; produce the output sample for time n - latency().
        Sample process(Sample in, Sample period, Sample ratio) noexcept {
            const double t = std::clamp(static_cast<double>(period), static_cast<double>(k_min_period),
                                        static_cast<double>(m_max_period));
            const double r = std::clamp(static_cast<double>(ratio), static_cast<double>(k_min_ratio),
                                        static_cast<double>(k_max_ratio));

            m_input[static_cast<size_t>(m_n % static_cast<std::int32_t>(m_input.size()))] = in;

            // Analysis marks: free-running, one per source period.
            const double now = static_cast<double>(m_n);
            while (m_next_mark <= now) {
                m_prev_mark = m_next_mark;
                m_have_mark = true;
                m_next_mark += t;
            }

            // Synthesis marks: one grain per t/r of output time. A grain centered at
            // synthesis mark s copies the newest fully-received analysis grain — the
            // one centered at m_prev_mark - t (its span [m-t, m+t] ends at m_prev_mark,
            // which the input cursor has already passed).
            while (m_next_synth <= now + t) {
                if (m_have_mark) {
                    const double m = m_prev_mark - t;
                    if (m - t >= now - static_cast<double>(m_input.size()) + 4.0 && m + t <= now) {
                        place_grain(m_next_synth, m, t, static_cast<Sample>(1.0 / r));
                    }
                }
                m_next_synth += t / r;
            }
            // Never let the scheduler fall behind the emit cursor (e.g. after clear()
            // races or extreme ratio jumps).
            const double emit = now - static_cast<double>(m_latency);
            if (m_next_synth < emit) {
                m_next_synth = emit;
            }

            // Emit, then release the slot for reuse.
            Sample y = Sample(0);
            if (m_n >= static_cast<std::int32_t>(m_latency)) {
                const size_t slot = static_cast<size_t>((m_n - static_cast<std::int32_t>(m_latency))
                                                        % static_cast<std::int32_t>(m_accum.size()));
                y                 = m_accum[slot];
                m_accum[slot]     = Sample(0);
            }
            ++m_n;
            if (m_n >= 2 * m_wrap) {
                shift_clock(-m_wrap);
            }
            return y;
        }

        /// Testing seam: advance the sample clock by `samples` (rounded down to a
        /// multiple of the ring size) as if that many samples had elapsed with the
        /// ring contents unchanged, wrapping the clock exactly as process() does.
        /// Lets a test cross the wrap, or an elapsed count past 2^31, in O(1)
        /// instead of processing that many samples. Fractional mark positions that
        /// are not multiples of the new magnitude's ulp may round by that ulp.
        /// Not a contract point; may change or disappear without a version note.
        void advance_clock_for_testing(std::uint64_t samples) noexcept {
            const std::uint64_t ring   = m_input.size();
            const std::uint64_t wrap   = static_cast<std::uint64_t>(m_wrap);
            const std::uint64_t target = static_cast<std::uint64_t>(m_n) + (samples / ring) * ring;
            const std::uint64_t folded = (target < wrap) ? target : wrap + (target - wrap) % wrap;
            shift_clock(static_cast<std::int32_t>(folded) - m_n);
        }

      private:
        /// Move the clock and every absolute position by `by` samples, a multiple of
        /// the ring size, so every ring index is unchanged. process() calls it with
        /// -m_wrap when the clock reaches 2 * m_wrap; m_wrap is at least one ring, so
        /// the clock stays at or above latency() and the warm-up guard stays true.
        /// That subtraction is exact in double: m_wrap is an integer, hence a
        /// multiple of every position's ulp, and the result's magnitude does not
        /// exceed the position's.
        void shift_clock(std::int32_t by) noexcept {
            const double d = static_cast<double>(by);
            m_n += by;
            m_next_mark += d;
            m_prev_mark += d;
            m_next_synth += d;
        }

        /// Overlap-add one Hann grain: output slots o in [s - t, s + t] receive the
        /// source at m + (o - s), read with Hermite interpolation (s is fractional).
        void place_grain(double s, double m, double t, Sample gain) noexcept {
            const long   first = static_cast<long>(std::ceil(s - t));
            const long   last  = static_cast<long>(std::floor(s + t));
            const double inv_t = 1.0 / t;
            const long   an    = static_cast<long>(m_accum.size());
            for (long o = first; o <= last; ++o) {
                const double delta = static_cast<double>(o) - s; // (-t, t)
                const double w     = 0.5 + 0.5 * std::cos(k_pi * delta * inv_t);
                const size_t slot  = static_cast<size_t>(((o % an) + an) % an);
                m_accum[slot] += gain * static_cast<Sample>(w) * read_hermite(m + delta);
            }
        }

        Sample read_hermite(double pos) const noexcept {
            const double fpos = std::floor(pos);
            const double frac = pos - fpos;
            const long   base = static_cast<long>(fpos);
            const long   rn   = static_cast<long>(m_input.size());
            const auto   at   = [&](long i) { return m_input[static_cast<size_t>(((i % rn) + rn) % rn)]; };
            const Sample xm1  = at(base - 1);
            const Sample x0   = at(base);
            const Sample x1   = at(base + 1);
            const Sample x2   = at(base + 2);
            const Sample c    = (x1 - xm1) * Sample(0.5);
            const Sample v    = x0 - x1;
            const Sample w    = c + v;
            const Sample a    = w + v + (x2 - x0) * Sample(0.5);
            const Sample b    = w + a;
            const Sample f    = static_cast<Sample>(frac);
            return (((a * f - b) * f + c) * f + x0);
        }

        static constexpr double k_pi = 3.14159265358979323846;

        /// Span of the sample clock: clock_wrap() is the largest multiple of the
        /// ring size not above this (or one ring, if the ring is larger). The
        /// value is arbitrary within (longest test run, 2^29): 2^18 keeps every
        /// position below 2^19 + 5 * max_period, i.e. with at least 33 fractional
        /// bits, and wraps every 5.4 s at 48 kHz (first at 10.9 s) so the wrap
        /// path is exercised routinely rather than once a shift. Not a contract
        /// point; the tests pin clock_wrap() <= 2^18 as a literal.
        static constexpr std::int32_t k_clock_span = std::int32_t{1} << 18;

        Sample m_max_period;
        size_t m_latency;

        std::vector<Sample> m_input;
        std::vector<Sample> m_accum;
        std::int32_t        m_n{0};    // sample clock, in [0, 2 * m_wrap); see shift_clock()
        std::int32_t        m_wrap{0}; // clock wrap period: a multiple of the ring size
        double              m_next_mark{0.0};
        double              m_prev_mark{0.0};
        bool                m_have_mark{false};
        double              m_next_synth{0.0};
    };

    /// Double-precision shifter — the desktop/golden-model profile.
    using psola = basic_psola<double>;

    /// Single-precision shifter — the embedded real-time profile.
    using psola32 = basic_psola<float>;

} // namespace tap::dsp
