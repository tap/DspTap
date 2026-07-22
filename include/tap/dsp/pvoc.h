/// @file pvoc.h
/// @brief Phase-vocoder pitch shifter with a fixed numeric contract.
// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// An STFT pitch shifter built on the published phase-vocoder literature:
// Hann-windowed frames at 4x overlap, per-bin instantaneous-frequency
// estimation from the frame-to-frame phase increment, and Laroche-Dolson-style
// peak-region shifting — spectral peaks are located, each peak's whole region
// is translated rigidly by an integer bin offset (the Hann mainlobe pattern
// stays intact), and the region is rotated by an accumulated residual phase so
// the overlap-add reconstructs the exact target frequency. Rigid translation
// preserves the phase relationships across each peak, which naive per-bin
// remapping (round(k*ratio) with free-running phase accumulation) destroys —
// that scheme measurably loses half the level on fractional ratios, a failure
// this class's own test battery pins.
//
// Optional formant preservation (set_formant) uses the classic source-filter
// method from the published speech-processing literature: an LPC spectral
// envelope (autocorrelation method + Levinson-Durbin, order 48) is estimated
// from each analysis frame, and every relocated bin is rescaled by
// envelope(target) / envelope(source) — the excitation moves, the envelope
// stays. Implemented from the literature only, per the project's IP policy.
//
// Transient smearing on percussive material remains the known trade of the
// phase-vocoder class. The transform is tap::dsp::basic_real_fft, so the float
// profile rides the vDSP / CMSIS-Helium backends where the build enables them.
// The packed spectrum uses fft.h's conjugated (W = exp(+2*pi*i/N)) convention;
// this class converts to the engineering convention at unpack and back at pack
// so the textbook phase math applies verbatim.

#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <type_traits>
#include <vector>

#include "tap/dsp/fft.h"

namespace tap::dsp {

    /// Phase-vocoder pitch shifter, parameterized over the sample type. Double
    /// is the desktop/golden profile, float the embedded/accelerated profile.
    ///
    /// Geometry (FFT size, 4x overlap) is fixed at construction and every
    /// buffer is allocated there; process() is noexcept and allocation-free,
    /// safe on a real-time audio thread. Latency is exactly latency() samples.
    /// At ratio == 1 every region's bin offset and residual are zero, so the
    /// output reconstructs the input's waveform delayed by exactly one FFT
    /// frame (pinned by the test battery).
    ///
    /// Contract points:
    ///   - process(x, ratio) consumes one input sample and produces one output
    ///     sample; ratio (2 = up an octave) is clamped to [1/4, 4] and sampled
    ///     once per analysis hop.
    ///   - fft_size must be a power of two >= 64. Frequency resolution and
    ///     transient smearing both scale with it; 1024 at 48 kHz is the
    ///     intended desktop operating point.
    ///   - Peak-to-region phase continuity is keyed by the region's target
    ///     bin, so a peak gliding across bins picks up a fresh phase register;
    ///     stationary and slowly-moving material is seamless.
    template <typename Sample>
    class basic_pvoc {
        static_assert(std::is_same_v<Sample, float> || std::is_same_v<Sample, double>,
                      "basic_pvoc supports the two Tap numeric profiles: float and double");

      public:
        static constexpr Sample k_min_ratio = Sample(0.25);
        static constexpr Sample k_max_ratio = Sample(4);
        static constexpr int    k_overlap   = 4;
        static constexpr int    k_lpc_order = 48;   // ~sr/1000 at the intended rates
        static constexpr double k_env_floor = 1e-6; // |A| guard when inverting the envelope
        static constexpr double k_max_boost = 16.0; // per-bin formant-correction gain cap

        /// @pre fft_size is a power of two, >= 64.
        explicit basic_pvoc(size_t fft_size = 1024)
            : m_n_size(static_cast<int>(fft_size))
            , m_hop(static_cast<int>(fft_size) / k_overlap)
            , m_bins(static_cast<int>(fft_size) / 2 + 1)
            , m_fft(fft_size) {
            assert(fft_size >= 64 && (fft_size & (fft_size - 1)) == 0);

            m_window.assign(static_cast<size_t>(m_n_size), Sample(0));
            for (int i = 0; i < m_n_size; ++i) {
                m_window[static_cast<size_t>(i)] = static_cast<Sample>(0.5 - 0.5 * std::cos(2.0 * k_pi * i / m_n_size));
            }

            // COLA normalization: steady-state sum of window^2 at the hop stride
            // (1.5 for Hann at 4x overlap; measured, not assumed). Since the hop
            // divides the window length, the periodic sum equals the steady state.
            double cola = 0.0;
            for (int k = 0; k < k_overlap; ++k) {
                const double w = static_cast<double>(m_window[static_cast<size_t>((m_hop / 2 + k * m_hop) % m_n_size)]);
                cola += w * w;
            }
            m_cola_norm = static_cast<Sample>(1.0 / cola);

            m_input.assign(static_cast<size_t>(m_n_size), Sample(0));
            m_accum.assign(static_cast<size_t>(3 * m_n_size), Sample(0));
            m_frame.assign(static_cast<size_t>(m_n_size), Sample(0));
            m_synth.assign(static_cast<size_t>(m_n_size), Sample(0));
            m_prev_phase.assign(static_cast<size_t>(m_bins), Sample(0));
            m_psi.assign(static_cast<size_t>(m_bins), Sample(0));
            m_mag.assign(static_cast<size_t>(m_bins), 0.0);
            m_true_bin.assign(static_cast<size_t>(m_bins), 0.0);
            m_peaks.reserve(static_cast<size_t>(m_bins));
            m_lpc_time.assign(static_cast<size_t>(m_n_size), Sample(0));
            m_lpc_work.assign(static_cast<size_t>(m_n_size), Sample(0));
            m_lpc_r.assign(static_cast<size_t>(k_lpc_order) + 1, 0.0);
            m_lpc_a.assign(static_cast<size_t>(k_lpc_order) + 1, 0.0);
            m_lpc_tmp.assign(static_cast<size_t>(k_lpc_order) + 1, 0.0);
            m_env.assign(static_cast<size_t>(m_bins), 1.0);
            clear();
        }

        size_t fft_size() const noexcept { return static_cast<size_t>(m_n_size); }
        size_t hop() const noexcept { return static_cast<size_t>(m_hop); }

        /// Emission delay of the shifter, in samples: one FFT frame.
        size_t latency() const noexcept { return static_cast<size_t>(m_n_size); }

        /// Enable LPC formant preservation: relocated bins are rescaled so the
        /// spectral envelope stays where the source put it while the excitation
        /// shifts. Off by default (the shifter then moves envelope and all).
        /// At ratio == 1 the correction is exactly unity, so identity holds
        /// either way. The LPC recursion runs in double in both profiles.
        void set_formant(bool on) noexcept { m_formant = on; }
        bool formant() const noexcept { return m_formant; }

        /// Zero all running state (buffers, phases, counters).
        void clear() noexcept {
            std::fill(m_input.begin(), m_input.end(), Sample(0));
            std::fill(m_accum.begin(), m_accum.end(), Sample(0));
            std::fill(m_prev_phase.begin(), m_prev_phase.end(), Sample(0));
            std::fill(m_psi.begin(), m_psi.end(), Sample(0));
            m_n = 0;
        }

        /// Consume one input sample; produce the output sample for time n - latency().
        Sample process(Sample in, Sample ratio) noexcept {
            const long in_size = static_cast<long>(m_input.size());
            const long an      = static_cast<long>(m_accum.size());

            m_input[static_cast<size_t>(m_n % in_size)] = in;

            if ((m_n + 1) % m_hop == 0 && m_n + 1 >= static_cast<long>(m_n_size)) {
                const double r = std::clamp(static_cast<double>(ratio), static_cast<double>(k_min_ratio),
                                            static_cast<double>(k_max_ratio));
                run_frame(r);
            }

            Sample y = Sample(0);
            if (m_n >= static_cast<long>(m_n_size)) {
                const size_t slot = static_cast<size_t>((m_n - static_cast<long>(m_n_size)) % an);
                y                 = m_accum[slot];
                m_accum[slot]     = Sample(0);
            }
            ++m_n;
            return y;
        }

      private:
        void run_frame(double r) noexcept {
            const long in_size = static_cast<long>(m_input.size());
            const long start   = m_n + 1 - static_cast<long>(m_n_size);

            // analysis: window the newest N samples and transform
            for (int i = 0; i < m_n_size; ++i) {
                m_frame[static_cast<size_t>(i)] =
                    m_input[static_cast<size_t>((start + i) % in_size)] * m_window[static_cast<size_t>(i)];
            }
            if (m_formant) {
                m_lpc_time = m_frame; // keep the windowed time frame; the FFT is in-place
            }
            m_fft.forward_inplace(m_frame.data());
            if (m_formant) {
                compute_envelope();
            }

            // per-bin magnitude and instantaneous frequency (engineering-convention
            // phases: conjugate fft.h's exp(+i) imaginary parts on unpack)
            const double expected = 2.0 * k_pi * static_cast<double>(m_hop) / static_cast<double>(m_n_size);
            for (int k = 1; k < m_bins - 1; ++k) {
                const double re    = static_cast<double>(m_frame[static_cast<size_t>(2 * k)]);
                const double im    = -static_cast<double>(m_frame[static_cast<size_t>(2 * k + 1)]);
                const double phase = std::atan2(im, re);

                double delta = phase - static_cast<double>(m_prev_phase[static_cast<size_t>(k)]) - expected * k;
                m_prev_phase[static_cast<size_t>(k)] = static_cast<Sample>(phase);
                delta -= 2.0 * k_pi * std::round(delta / (2.0 * k_pi));

                m_mag[static_cast<size_t>(k)]      = std::sqrt(re * re + im * im);
                m_true_bin[static_cast<size_t>(k)] = static_cast<double>(k) + delta / expected;
            }

            // peaks: local maxima over +-2 bins (the Laroche-Dolson criterion),
            // gated at -80 dB below the frame's strongest bin so noise-floor
            // maxima neither claim regions nor collide with a real peak's phase
            // register (their bins are simply dropped — they are inaudible)
            double max_mag = 0.0;
            for (int k = 1; k < m_bins - 1; ++k) {
                max_mag = std::max(max_mag, m_mag[static_cast<size_t>(k)]);
            }
            const double floor_mag = max_mag * 1e-4;
            m_peaks.clear();
            for (int k = 3; k < m_bins - 3; ++k) {
                const double m = m_mag[static_cast<size_t>(k)];
                if (m > m_mag[static_cast<size_t>(k - 1)] && m > m_mag[static_cast<size_t>(k - 2)]
                    && m >= m_mag[static_cast<size_t>(k + 1)] && m >= m_mag[static_cast<size_t>(k + 2)]
                    && m > floor_mag) {
                    m_peaks.push_back(k);
                }
            }

            // synthesis: translate each peak's region rigidly by an integer bin
            // offset and rotate it by the accumulated residual phase
            std::fill(m_synth.begin(), m_synth.end(), Sample(0));
            m_synth[0]        = m_frame[0]; // DC and Nyquist pass through untouched: they
            m_synth[1]        = m_frame[1]; // cannot be relocated, and identity stays exact
            const int n_peaks = static_cast<int>(m_peaks.size());
            int       lo      = 1;
            for (int pi = 0; pi < n_peaks; ++pi) {
                const int p         = m_peaks[static_cast<size_t>(pi)];
                const int hi        = (pi == n_peaks - 1) ? m_bins - 2 : (p + m_peaks[static_cast<size_t>(pi + 1)]) / 2;
                const int region_lo = lo;
                lo                  = hi + 1; // next region starts past this one — no bin doubles up

                const double fp    = m_true_bin[static_cast<size_t>(p)];
                const int    shift = static_cast<int>(std::lround(fp * (r - 1.0)));
                const int    c     = p + shift;
                if (c < 1 || c > m_bins - 2) {
                    continue;
                }

                // The integer translation is frame-relative (its implicit modulator
                // e^(2*pi*i*shift*n/N) restarts each frame), so psi must accumulate
                // the FULL per-hop frequency difference fp*(r-1) — subtracting the
                // shift here desynchronizes the overlap-add and guts the level.
                // Keyed by the region's target bin for frame-to-frame continuity.
                const double resid = expected * fp * (r - 1.0);
                double       psi   = static_cast<double>(m_psi[static_cast<size_t>(c)]) + resid;
                psi -= 2.0 * k_pi * std::round(psi / (2.0 * k_pi));
                m_psi[static_cast<size_t>(c)] = static_cast<Sample>(psi);

                const double cs = std::cos(psi);
                const double sn = std::sin(psi);
                for (int k = region_lo; k <= hi; ++k) {
                    const int j = k + shift;
                    if (j < 1 || j > m_bins - 2) {
                        continue;
                    }
                    double re = static_cast<double>(m_frame[static_cast<size_t>(2 * k)]);
                    double im = -static_cast<double>(m_frame[static_cast<size_t>(2 * k + 1)]);
                    if (m_formant) {
                        // keep the envelope in place: excitation from bin k now sits
                        // at bin j, so trade envelope(source) for envelope(target)
                        const double g = std::clamp(m_env[static_cast<size_t>(j)] / m_env[static_cast<size_t>(k)],
                                                    1.0 / k_max_boost, k_max_boost);
                        re *= g;
                        im *= g;
                    }
                    m_synth[static_cast<size_t>(2 * j)] += static_cast<Sample>(re * cs - im * sn);
                    m_synth[static_cast<size_t>(2 * j + 1)] -= static_cast<Sample>(re * sn + im * cs); // conjugate back
                }
            }

            m_fft.inverse_inplace(m_synth.data());

            // synthesis window + COLA-normalized overlap-add (fold in the raw
            // inverse's 2/N normalization here, one multiply per sample)
            const long   an   = static_cast<long>(m_accum.size());
            const Sample norm = m_cola_norm * (Sample(2) / static_cast<Sample>(m_n_size));
            for (int i = 0; i < m_n_size; ++i) {
                const size_t slot = static_cast<size_t>((start + i) % an);
                m_accum[slot] += m_synth[static_cast<size_t>(i)] * m_window[static_cast<size_t>(i)] * norm;
            }
        }

        /// LPC spectral envelope of the current analysis frame (autocorrelation
        /// method + Levinson-Durbin, evaluated over all bins with one extra FFT
        /// of the prediction polynomial). m_env[k] = 1 / max(|A(e^jw_k)|, floor);
        /// only ratios of m_env are used, so the LPC gain term cancels out.
        void compute_envelope() noexcept {
            const int p = k_lpc_order;

            for (int lag = 0; lag <= p; ++lag) {
                double acc = 0.0;
                for (int n = 0; n < m_n_size - lag; ++n) {
                    acc += static_cast<double>(m_lpc_time[static_cast<size_t>(n)])
                           * static_cast<double>(m_lpc_time[static_cast<size_t>(n + lag)]);
                }
                m_lpc_r[static_cast<size_t>(lag)] = acc;
            }
            if (m_lpc_r[0] <= 0.0) {
                std::fill(m_env.begin(), m_env.end(), 1.0); // silent frame: flat envelope
                return;
            }
            m_lpc_r[0] *= 1.0 + 1e-9; // tiny ridge keeps the recursion well-posed

            // Levinson-Durbin, always in double (order-48 float recursion is not safe)
            std::fill(m_lpc_a.begin(), m_lpc_a.end(), 0.0);
            m_lpc_a[0]   = 1.0;
            double error = m_lpc_r[0];
            for (int i = 1; i <= p; ++i) {
                double acc = m_lpc_r[static_cast<size_t>(i)];
                for (int j = 1; j < i; ++j) {
                    acc += m_lpc_a[static_cast<size_t>(j)] * m_lpc_r[static_cast<size_t>(i - j)];
                }
                const double k = -acc / error;
                for (int j = 0; j <= i; ++j) {
                    m_lpc_tmp[static_cast<size_t>(j)] =
                        m_lpc_a[static_cast<size_t>(j)] + k * m_lpc_a[static_cast<size_t>(i - j)];
                }
                std::copy(m_lpc_tmp.begin(), m_lpc_tmp.begin() + i + 1, m_lpc_a.begin());
                error *= (1.0 - k * k);
                if (error <= 0.0) {
                    break;
                }
            }

            // |A| over the bins: FFT of the prediction polynomial (same size, same engine)
            std::fill(m_lpc_work.begin(), m_lpc_work.end(), Sample(0));
            for (int i = 0; i <= p; ++i) {
                m_lpc_work[static_cast<size_t>(i)] = static_cast<Sample>(m_lpc_a[static_cast<size_t>(i)]);
            }
            m_fft.forward_inplace(m_lpc_work.data());
            m_env[0] = 1.0 / std::max(std::abs(static_cast<double>(m_lpc_work[0])), k_env_floor);
            m_env[static_cast<size_t>(m_bins - 1)] =
                1.0 / std::max(std::abs(static_cast<double>(m_lpc_work[1])), k_env_floor);
            for (int k = 1; k < m_bins - 1; ++k) {
                const double re               = static_cast<double>(m_lpc_work[static_cast<size_t>(2 * k)]);
                const double im               = static_cast<double>(m_lpc_work[static_cast<size_t>(2 * k + 1)]);
                m_env[static_cast<size_t>(k)] = 1.0 / std::max(std::sqrt(re * re + im * im), k_env_floor);
            }
        }

        static constexpr double k_pi = 3.14159265358979323846;

        int                    m_n_size;
        int                    m_hop;
        int                    m_bins;
        basic_real_fft<Sample> m_fft;
        Sample                 m_cola_norm{Sample(0)};

        std::vector<Sample> m_window;
        std::vector<Sample> m_input;
        std::vector<Sample> m_accum;
        std::vector<Sample> m_frame;
        std::vector<Sample> m_synth;
        std::vector<Sample> m_prev_phase;
        std::vector<Sample> m_psi;
        std::vector<double> m_mag;
        std::vector<double> m_true_bin;
        std::vector<int>    m_peaks;
        std::vector<Sample> m_lpc_time;
        std::vector<Sample> m_lpc_work;
        std::vector<double> m_lpc_r;
        std::vector<double> m_lpc_a;
        std::vector<double> m_lpc_tmp;
        std::vector<double> m_env;
        bool                m_formant{false};
        long                m_n{0};
    };

    /// Double-precision shifter — the desktop/golden-model profile.
    using pvoc = basic_pvoc<double>;

    /// Single-precision shifter — the embedded/accelerated profile.
    using pvoc32 = basic_pvoc<float>;

} // namespace tap::dsp
