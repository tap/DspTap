/// @file yin.h
/// @brief YIN pitch detector with a fixed numeric contract.
// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// The time-domain YIN estimator of de Cheveigné & Kawahara (2002), steps 1-5:
// squared-difference function, cumulative-mean normalization, absolute
// threshold with local-minimum descent, and parabolic interpolation of the
// selected lag. Step 6 of the paper ("best local estimate") is intentionally
// omitted — it trades latency for a marginal gain this contract does not need.
//
// Grown out of the decimated autocorrelation follower embedded in the TapTools
// pitchaccum kernel; promoted here (full-rate, sub-sample, normalized) so every
// Tap library shares one pitch-detection contract. The difference-function
// inner loop is the hot path and is kept a plain contiguous loop on purpose:
// it is the designated candidate for Helium-MVE / HVX backends behind this
// same interface, with the scalar build remaining the golden model (the
// fft.h backend pattern).

#pragma once

#include <cassert>
#include <cstddef>
#include <type_traits>
#include <vector>

namespace tap::dsp {

    /// One analysis frame's estimate.
    ///
    /// `period` is the detected period in samples (fractional — the lag of the
    /// selected CMND minimum refined by parabolic interpolation), or 0 when the
    /// frame is unvoiced: no candidate lag dipped below the threshold.
    /// `aperiodicity` is the cumulative-mean-normalized difference at the
    /// selected lag — 0 for a perfectly periodic frame, rising toward (and
    /// past) 1 as periodicity disappears. When unvoiced it reports the global
    /// minimum over the lag range, so a caller can still rank near-misses.
    template <typename Sample>
    struct pitch_result {
        Sample period;
        Sample aperiodicity;

        bool voiced() const noexcept { return period > Sample(0); }
    };

    /// YIN pitch detector, parameterized over the sample type. Both
    /// instantiations run the identical algorithm in their own precision:
    /// double is the desktop/golden profile, float the embedded profile, and
    /// the test battery pins their cross-precision agreement.
    ///
    /// Geometry is fixed at construction and every buffer is allocated there;
    /// analyze() is noexcept and allocation-free, safe on a real-time audio
    /// thread. Cost is O(window * tau_max) multiplies per call — callers set
    /// the analysis hop to budget it.
    ///
    /// Contract points:
    ///   - analyze(x) reads frame_size() = window + tau_max samples, oldest
    ///     first: x[j] is compared against x[j + tau] for j in [0, window).
    ///   - Lags searched are [tau_min, tau_max]; bound them from the caller's
    ///     frequency range as tau = sample_rate / frequency.
    ///   - threshold is the paper's absolute threshold on the normalized
    ///     difference (default 0.1); raising it admits noisier frames as
    ///     voiced, lowering it rejects them.
    template <typename Sample>
    class basic_yin {
        static_assert(std::is_same_v<Sample, float> || std::is_same_v<Sample, double>,
                      "basic_yin supports the two Tap numeric profiles: float and double");

      public:
        static constexpr Sample k_default_threshold = Sample(0.1);

        /// @pre 2 <= tau_min < tau_max, window >= tau_max (the paper's W >= tau_max
        /// keeps every searched lag fully supported by the integration window).
        basic_yin(size_t window, size_t tau_min, size_t tau_max)
            : m_window(static_cast<int>(window))
            , m_tau_min(static_cast<int>(tau_min))
            , m_tau_max(static_cast<int>(tau_max)) {
            assert(tau_min >= 2 && tau_min < tau_max && window >= tau_max);
            m_cmnd.assign(static_cast<size_t>(m_tau_max) + 1, Sample(0));
        }

        size_t window() const noexcept { return static_cast<size_t>(m_window); }
        size_t tau_min() const noexcept { return static_cast<size_t>(m_tau_min); }
        size_t tau_max() const noexcept { return static_cast<size_t>(m_tau_max); }

        /// Samples analyze() reads: the integration window plus the deepest lag.
        size_t frame_size() const noexcept { return static_cast<size_t>(m_window + m_tau_max); }

        void   set_threshold(Sample t) noexcept { m_threshold = t; }
        Sample threshold() const noexcept { return m_threshold; }

        /// Analyze one frame of frame_size() samples, oldest first.
        pitch_result<Sample> analyze(const Sample* x) noexcept {
            // Steps 1-3: squared difference per lag, normalized in the same pass
            // by the cumulative mean. m_cmnd[tau] is d'(tau); d'(0) := 1.
            Sample* const cmnd = m_cmnd.data();
            cmnd[0]            = Sample(1);

            Sample cumulative = Sample(0);
            for (int tau = 1; tau <= m_tau_max; ++tau) {
                // Hot loop — keep contiguous and branch-free (MVE/HVX candidate).
                Sample d = Sample(0);
                for (int j = 0; j < m_window; ++j) {
                    const Sample diff = x[j] - x[j + tau];
                    d += diff * diff;
                }
                cumulative += d;
                cmnd[tau] = (cumulative > Sample(0)) ? d * static_cast<Sample>(tau) / cumulative : Sample(1);
            }

            // Step 4: first lag under the threshold, descended to its local
            // minimum. Track the global minimum for the unvoiced report.
            int    best_tau   = 0;
            Sample global_min = cmnd[m_tau_min];
            for (int tau = m_tau_min; tau <= m_tau_max; ++tau) {
                if (cmnd[tau] < global_min) {
                    global_min = cmnd[tau];
                }
                if (best_tau == 0 && cmnd[tau] < m_threshold) {
                    int t = tau;
                    while (t + 1 <= m_tau_max && cmnd[t + 1] < cmnd[t]) {
                        ++t;
                    }
                    best_tau = t;
                }
            }

            if (best_tau == 0) {
                return {Sample(0), global_min};
            }

            // Step 5: parabolic interpolation of the minimum over its immediate
            // neighbors, in the normalized domain, clamped to half a lag.
            Sample period = static_cast<Sample>(best_tau);
            if (best_tau > 1 && best_tau < m_tau_max) {
                const Sample dm    = cmnd[best_tau - 1];
                const Sample d0    = cmnd[best_tau];
                const Sample dp    = cmnd[best_tau + 1];
                const Sample denom = dm - Sample(2) * d0 + dp;
                if (denom > Sample(0)) {
                    Sample offset = (dm - dp) / (Sample(2) * denom);
                    if (offset > Sample(0.5)) {
                        offset = Sample(0.5);
                    }
                    if (offset < Sample(-0.5)) {
                        offset = Sample(-0.5);
                    }
                    period += offset;
                }
            }
            return {period, cmnd[best_tau]};
        }

      private:
        int                 m_window;
        int                 m_tau_min;
        int                 m_tau_max;
        Sample              m_threshold{k_default_threshold};
        std::vector<Sample> m_cmnd;
    };

    /// Double-precision detector — the desktop/golden-model profile.
    using yin = basic_yin<double>;

    /// Single-precision detector — the embedded real-time profile.
    using yin32 = basic_yin<float>;

} // namespace tap::dsp
