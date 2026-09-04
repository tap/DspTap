/// @file log_mel.h
/// @brief Streaming log-mel / PCEN feature front end with a fixed numeric contract.
// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// The analysis front end of a keyword spotter: windowed real FFT, power
// spectrum, triangular mel filterbank, then either a floored, affinely
// normalized log10 (the plain-log path, default) or per-channel energy
// normalization (PCEN, Wang et al. 2017: "Trainable frontend for robust and
// far-field keyword spotting", ICASSP). All formulas are the published
// textbook ones; nothing here is reverse-engineered from a shipping product.
//
// OWNERSHIP RULE (MuTap's wake-word plan, section 5). This header owns the
// FORMULA-level contract: the mel scale (HTK), the filter shape (unit-peak
// triangles), the window (periodic Hann or its square root), the FFT size and
// where the zero padding sits, the frame alignment, the DC policy, the
// pre-emphasis form, the log form and the PCEN recursion. Everything a
// trainer might tune — band count, fmin/fmax, the log floor/shift/scale, every
// PCEN parameter, the pre-emphasis coefficient — is RUNTIME GEOMETRY in
// log_mel_geometry, carried by a trained model's weights and validated at
// load, so retraining never touches this header. k_contract_version changes
// when a formula changes; a model records the version it was trained against.
//
// The formula-level contract, as numbers:
//   - Frame alignment: frame t consists of the `frame` most recent samples
//     when sample (t+1)*hop - 1 has arrived; the history before sample 0 is
//     zero. No centring, no prepended zero frame. Latency = `frame` samples.
//   - Pre-emphasis (when preemphasis != 0): y[n] = x[n] - c * x[n-1] on the
//     continuous stream, x[-1] = 0.
//   - Window: periodic Hann w[n] = 0.5 - 0.5 cos(2 pi n / frame), n in
//     [0, frame), or its square root (mel_window::sqrt_hann). Stated by the
//     geometry, never implied.
//   - FFT: fft_size >= frame, a power of two; the windowed frame occupies
//     [0, frame) and zeros occupy [frame, fft_size). The transform is the
//     unnormalized real DFT of tap::dsp::basic_real_fft (Ooura contract), so
//     the power spectrum is |X_k|^2 with X_k = sum x[n] e^{-i 2 pi k n / N}
//     up to the sign of the imaginary part, which power discards.
//   - Bin frequencies: f_k = k * sample_rate / fft_size, k in [0, fft_size/2].
//   - Mel scale (HTK): mel(f) = 2595 log10(1 + f / 700). Band edges are
//     bands + 2 points equally spaced in mel between fmin_hz and fmax_hz.
//   - Filters: unit-PEAK triangles (weight 1 at the band's centre edge,
//     linear to 0 at its neighbours), evaluated at the bin frequencies. Not
//     Slaney area-normalized. With fmin_hz > 0, bin 0 (DC) has zero weight in
//     every band.
//   - Band energy: E_b = sum_k w_b[k] |X_k|^2 over the band's support.
//   - Plain-log feature: (log10(E_b + log_floor) + log_shift) / log_scale.
//     Defaults 1e-10 / 5 / 5 map E in [1e-10, 1] onto [-1, 1]. log_floor is
//     stated relative to unit full scale: the float profile's power rounding
//     floor (~1e-14 |x|^2) sits well below it for |x| <= 1.
//   - PCEN feature (pcen.enabled): M_b[t] = (1 - s) M_b[t-1] + s E_b[t];
//     out = (E_b / (eps + M_b)^alpha + delta)^r - delta^r. The smoother
//     state after reset() (and at construction) is primed with the FIRST
//     frame's energy, so the first output sees M = E. Defaults are the
//     paper's: s 0.025, alpha 0.98, delta 2, r 0.5, eps 1e-6. The PCEN path
//     replaces the log path; the affine constants do not apply to it.
//
// Profiles: double is the golden model, float the embedded profile; both run
// the identical algorithm in their own precision and the tests pin their
// agreement as a measured number on a signal that excites every band.
// Geometry is fixed at construction and every buffer allocated there;
// process_hop() and process() are noexcept and allocation-free. No double
// arithmetic occurs on the float path: the Cortex-M33 (RP2350) profile has no
// FP64 and this front end must not fall into soft-float there.

#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <numbers>
#include <type_traits>
#include <vector>

#include "tap/dsp/fft.h"

namespace tap::dsp {

    /// Analysis window shape. Both are the periodic form (denominator `frame`).
    enum class mel_window : std::uint8_t {
        hann,     ///< 0.5 - 0.5 cos(2 pi n / frame)
        sqrt_hann ///< sqrt of the above (the suppressor's analysis window)
    };

    /// PCEN parameters (Wang et al. 2017), all runtime geometry.
    struct pcen_params {
        bool   enabled  = false;
        double smoother = 0.025; ///< s: one-pole smoother coefficient per frame
        double alpha    = 0.98;  ///< gain-normalization strength
        double delta    = 2.0;   ///< bias before the root
        double power    = 0.5;   ///< r: compression root
        double epsilon  = 1e-6;  ///< stabilizer added to the smoothed energy

        bool valid() const noexcept {
            return smoother > 0.0 && smoother <= 1.0 && alpha >= 0.0 && alpha <= 1.0 && delta > 0.0 && power > 0.0
                   && power <= 1.0 && epsilon > 0.0;
        }
    };

    /// The runtime geometry of the front end: everything a trained model
    /// carries. Defaults are the reference geometry MuTap's keyword spotter is
    /// developed against (16 kHz, 25 ms frame, 10 ms hop, 512-point FFT, 40
    /// bands 20-7600 Hz).
    struct log_mel_geometry {
        /// Formula-level contract version; bump when any formula in this
        /// header changes. Trained models record the version they used.
        static constexpr std::uint32_t k_contract_version = 1;

        double      sample_rate = 16000.0;
        std::size_t frame       = 400; ///< analysis frame, samples
        std::size_t hop         = 160; ///< frames advance by this many samples
        std::size_t fft_size    = 512; ///< power of two, >= frame
        std::size_t bands       = 40;
        double      fmin_hz     = 20.0;
        double      fmax_hz     = 7600.0;
        mel_window  window      = mel_window::hann;
        double      preemphasis = 0.0; ///< 0 = off
        double      log_floor   = 1e-10;
        double      log_shift   = 5.0;
        double      log_scale   = 5.0;
        pcen_params pcen{};

        std::size_t bins() const noexcept { return fft_size / 2 + 1; }

        bool valid() const noexcept {
            const bool pow2 = fft_size >= 4 && (fft_size & (fft_size - 1)) == 0;
            return sample_rate > 0.0 && hop >= 1 && frame >= hop && pow2 && fft_size >= frame && bands >= 1
                   && fmin_hz >= 0.0 && fmax_hz > fmin_hz && fmax_hz <= 0.5 * sample_rate && log_floor > 0.0
                   && log_scale > 0.0 && pcen.valid();
        }
    };

    /// HTK mel scale.
    inline double hz_to_mel(double hz) noexcept {
        return 2595.0 * std::log10(1.0 + hz / 700.0);
    }
    inline double mel_to_hz(double mel) noexcept {
        return 700.0 * (std::pow(10.0, mel / 2595.0) - 1.0);
    }

    /// Streaming log-mel / PCEN front end, parameterized over the sample type.
    ///
    /// Contract points beyond the file header:
    ///   - process_hop(in, out) consumes exactly hop samples (oldest first) and
    ///     writes bands features for the frame ending at the newest sample.
    ///   - process(in, n, out, max_frames) accepts any n, buffering the
    ///     remainder for the next call; every hop-aligned frame is produced
    ///     exactly as process_hop would, so chunking never changes a feature.
    ///     Frames beyond max_frames are dropped (size out with frames_for(n)).
    ///   - reset() clears the history, the pre-emphasis state, the pending
    ///     input and the PCEN smoother; the next frame is a first frame.
    template <typename Sample>
    class basic_log_mel {
        static_assert(std::is_same_v<Sample, float> || std::is_same_v<Sample, double>,
                      "basic_log_mel supports the two Tap numeric profiles: float and double");

      public:
        /// @pre g.valid()
        explicit basic_log_mel(const log_mel_geometry& g)
            : m_g(g)
            , m_fft(g.fft_size) {
            assert(g.valid());
            m_history.assign(g.frame, Sample(0));
            m_pending.assign(g.hop, Sample(0));
            m_window.assign(g.frame, Sample(0));
            m_spec.assign(g.fft_size, Sample(0));
            m_power.assign(g.bins(), Sample(0));
            m_energy.assign(g.bands, Sample(0));
            m_pcen_state.assign(g.bands, Sample(0));
            build_window();
            build_filterbank();
        }

        const log_mel_geometry& geometry() const noexcept { return m_g; }
        std::size_t             bands() const noexcept { return m_g.bands; }
        std::size_t             hop() const noexcept { return m_g.hop; }
        /// Samples between a sample arriving and the frame containing it
        /// being complete: the frame length.
        std::size_t latency_samples() const noexcept { return m_g.frame; }
        /// Frames process(n) will write from the current pending state.
        std::size_t frames_for(std::size_t n) const noexcept { return (m_pending_count + n) / m_g.hop; }
        std::size_t pending_samples() const noexcept { return m_pending_count; }

        /// Mel-scale band edge i of bands + 2, in Hz (edge b+1 is band b's centre).
        double band_edge_hz(std::size_t i) const noexcept {
            const double lo = hz_to_mel(m_g.fmin_hz);
            const double hi = hz_to_mel(m_g.fmax_hz);
            return mel_to_hz(lo + (hi - lo) * static_cast<double>(i) / static_cast<double>(m_g.bands + 1));
        }
        /// Filter weight of band b at FFT bin k (0 outside the band's support).
        Sample band_weight(std::size_t b, std::size_t k) const noexcept {
            const band& bd = m_bands[b];
            if (k < bd.first_bin || k >= bd.first_bin + bd.weights.size()) {
                return Sample(0);
            }
            return bd.weights[k - bd.first_bin];
        }

        void reset() noexcept {
            std::fill(m_history.begin(), m_history.end(), Sample(0));
            std::fill(m_pcen_state.begin(), m_pcen_state.end(), Sample(0));
            m_pending_count = 0;
            m_preemph_prev  = Sample(0);
            m_pcen_primed   = false;
        }

        /// One hop in, one frame of features out.
        void process_hop(const Sample* in, Sample* features) noexcept {
            const std::size_t frame = m_g.frame;
            const std::size_t hop   = m_g.hop;
            // Shift the history and append the hop, pre-emphasized on the stream.
            std::memmove(m_history.data(), m_history.data() + hop, (frame - hop) * sizeof(Sample));
            Sample* tail = m_history.data() + (frame - hop);
            if (m_g.preemphasis != 0.0) {
                const Sample c = static_cast<Sample>(m_g.preemphasis);
                for (std::size_t i = 0; i < hop; ++i) {
                    tail[i]        = in[i] - c * m_preemph_prev;
                    m_preemph_prev = in[i];
                }
            }
            else {
                std::memcpy(tail, in, hop * sizeof(Sample));
            }
            analyze(features);
        }

        /// Any number of samples in; returns the number of frames written.
        std::size_t process(const Sample* in, std::size_t n, Sample* features, std::size_t max_frames) noexcept {
            std::size_t       written = 0;
            const std::size_t hop     = m_g.hop;
            std::size_t       i       = 0;
            while (i < n) {
                const std::size_t take = std::min(hop - m_pending_count, n - i);
                std::memcpy(m_pending.data() + m_pending_count, in + i, take * sizeof(Sample));
                m_pending_count += take;
                i += take;
                if (m_pending_count == hop) {
                    m_pending_count = 0;
                    if (written < max_frames) {
                        process_hop(m_pending.data(), features + written * m_g.bands);
                        ++written;
                    }
                }
            }
            return written;
        }

      private:
        struct band {
            std::size_t         first_bin = 0;
            std::vector<Sample> weights;
        };

        void build_window() {
            const auto n = static_cast<double>(m_g.frame);
            for (std::size_t i = 0; i < m_g.frame; ++i) {
                const double w = 0.5 - 0.5 * std::cos(2.0 * std::numbers::pi * static_cast<double>(i) / n);
                m_window[i]    = static_cast<Sample>(m_g.window == mel_window::sqrt_hann ? std::sqrt(w) : w);
            }
        }

        void build_filterbank() {
            m_bands.resize(m_g.bands);
            const std::size_t bins = m_g.bins();
            for (std::size_t b = 0; b < m_g.bands; ++b) {
                const double        lo    = band_edge_hz(b);
                const double        mid   = band_edge_hz(b + 1);
                const double        hi    = band_edge_hz(b + 2);
                std::size_t         first = bins;
                std::size_t         last  = 0; // one past
                std::vector<Sample> w;
                for (std::size_t k = 0; k < bins; ++k) {
                    const double f    = static_cast<double>(k) * m_g.sample_rate / static_cast<double>(m_g.fft_size);
                    const double rise = (f - lo) / (mid - lo);
                    const double fall = (hi - f) / (hi - mid);
                    const double v    = std::max(0.0, std::min(rise, fall));
                    if (v > 0.0) {
                        if (first == bins) {
                            first = k;
                        }
                        last = k + 1;
                    }
                }
                if (first == bins) { // no bin falls inside this band (very narrow band / coarse FFT)
                    first = 0;
                    last  = 0;
                }
                w.resize(last - first);
                for (std::size_t k = first; k < last; ++k) {
                    const double f    = static_cast<double>(k) * m_g.sample_rate / static_cast<double>(m_g.fft_size);
                    const double rise = (f - lo) / (mid - lo);
                    const double fall = (hi - f) / (hi - mid);
                    w[k - first]      = static_cast<Sample>(std::max(0.0, std::min(rise, fall)));
                }
                m_bands[b].first_bin = first;
                m_bands[b].weights   = std::move(w);
            }
        }

        void analyze(Sample* features) noexcept {
            const std::size_t frame = m_g.frame;
            const std::size_t n     = m_g.fft_size;
            for (std::size_t i = 0; i < frame; ++i) {
                m_spec[i] = m_history[i] * m_window[i];
            }
            std::fill(m_spec.begin() + static_cast<std::ptrdiff_t>(frame), m_spec.end(), Sample(0));
            m_fft.forward_inplace(m_spec.data());
            // Ooura packing: [0] = DC (real), [1] = Nyquist (real), then (re, im) pairs.
            m_power[0]     = m_spec[0] * m_spec[0];
            m_power[n / 2] = m_spec[1] * m_spec[1];
            for (std::size_t k = 1; k < n / 2; ++k) {
                const Sample re = m_spec[2 * k];
                const Sample im = m_spec[2 * k + 1];
                m_power[k]      = re * re + im * im;
            }
            for (std::size_t b = 0; b < m_g.bands; ++b) {
                const band&   bd  = m_bands[b];
                Sample        acc = Sample(0);
                const Sample* p   = m_power.data() + bd.first_bin;
                for (std::size_t j = 0; j < bd.weights.size(); ++j) {
                    acc += bd.weights[j] * p[j];
                }
                m_energy[b] = acc;
            }
            if (m_g.pcen.enabled) {
                pcen(features);
            }
            else {
                const Sample floor = static_cast<Sample>(m_g.log_floor);
                const Sample shift = static_cast<Sample>(m_g.log_shift);
                const Sample scale = static_cast<Sample>(m_g.log_scale);
                for (std::size_t b = 0; b < m_g.bands; ++b) {
                    features[b] = (std::log10(m_energy[b] + floor) + shift) / scale;
                }
            }
        }

        void pcen(Sample* features) noexcept {
            const Sample s     = static_cast<Sample>(m_g.pcen.smoother);
            const Sample alpha = static_cast<Sample>(m_g.pcen.alpha);
            const Sample delta = static_cast<Sample>(m_g.pcen.delta);
            const Sample r     = static_cast<Sample>(m_g.pcen.power);
            const Sample eps   = static_cast<Sample>(m_g.pcen.epsilon);
            const Sample bias  = std::pow(delta, r);
            if (!m_pcen_primed) {
                std::copy(m_energy.begin(), m_energy.end(), m_pcen_state.begin());
                m_pcen_primed = true;
            }
            for (std::size_t b = 0; b < m_g.bands; ++b) {
                const Sample m  = (Sample(1) - s) * m_pcen_state[b] + s * m_energy[b];
                m_pcen_state[b] = m;
                features[b]     = std::pow(m_energy[b] / std::pow(eps + m, alpha) + delta, r) - bias;
            }
        }

        log_mel_geometry       m_g;
        basic_real_fft<Sample> m_fft;
        std::vector<Sample>    m_history; ///< the last `frame` (pre-emphasized) samples, oldest first
        std::vector<Sample>    m_pending; ///< partial hop buffered by process()
        std::vector<Sample>    m_window;
        std::vector<Sample>    m_spec;
        std::vector<Sample>    m_power;
        std::vector<Sample>    m_energy;
        std::vector<Sample>    m_pcen_state;
        std::vector<band>      m_bands;
        std::size_t            m_pending_count = 0;
        Sample                 m_preemph_prev  = Sample(0);
        bool                   m_pcen_primed   = false;
    };

    using log_mel   = basic_log_mel<double>; ///< golden model
    using log_mel32 = basic_log_mel<float>;  ///< embedded profile

} // namespace tap::dsp
