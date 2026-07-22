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
    template <typename Sample>
    class basic_psola {
        static_assert(std::is_same_v<Sample, float> || std::is_same_v<Sample, double>,
                      "basic_psola supports the two Tap numeric profiles: float and double");

      public:
        static constexpr Sample k_min_period = Sample(8);
        static constexpr Sample k_min_ratio  = Sample(0.25);
        static constexpr Sample k_max_ratio  = Sample(4);

        /// @pre max_period >= 16 — the deepest period process() will be given.
        explicit basic_psola(size_t max_period)
            : m_max_period(static_cast<Sample>(max_period))
            , m_latency(2 * max_period + 2) {
            assert(max_period >= 16);
            // Input history: a grain reaches back to (mark - period) and marks lag the
            // input cursor by up to two periods -> three periods of history plus slack.
            m_input.assign(4 * max_period + 8, Sample(0));
            // Output accumulator: emission lags by latency(); grains extend up to one
            // period past their mark -> latency + period ahead of the emit cursor.
            m_accum.assign(4 * max_period + 8, Sample(0));
            clear();
        }

        /// Emission delay of the shifter, in samples (fixed at construction).
        size_t latency() const noexcept { return m_latency; }

        size_t max_period() const noexcept { return static_cast<size_t>(m_max_period); }

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

            m_input[static_cast<size_t>(m_n % static_cast<long>(m_input.size()))] = in;

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
            if (m_n >= static_cast<long>(m_latency)) {
                const size_t slot =
                    static_cast<size_t>((m_n - static_cast<long>(m_latency)) % static_cast<long>(m_accum.size()));
                y             = m_accum[slot];
                m_accum[slot] = Sample(0);
            }
            ++m_n;
            return y;
        }

      private:
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

        Sample m_max_period;
        size_t m_latency;

        std::vector<Sample> m_input;
        std::vector<Sample> m_accum;
        long                m_n{0};
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
