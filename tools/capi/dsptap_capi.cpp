/// @file dsptap_capi.cpp
/// @brief C ABI over the DspTap primitives — see dsptap_capi.h.
// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.

#include "dsptap_capi.h"

#include <memory>
#include <vector>

#include "tap/dsp/decimate.h"
#include "tap/dsp/log_mel.h"
#include "tap/dsp/psola.h"
#include "tap/dsp/pvoc.h"
#include "tap/dsp/yin.h"

namespace {

    tap::dsp::yin* as_yin(dsptap_yin h) {
        return static_cast<tap::dsp::yin*>(h);
    }
    tap::dsp::psola* as_psola(dsptap_psola h) {
        return static_cast<tap::dsp::psola*>(h);
    }
    tap::dsp::pvoc* as_pvoc(dsptap_pvoc h) {
        return static_cast<tap::dsp::pvoc*>(h);
    }

    // log_mel's geometry is fixed at construction; the C ABI setters rebuild.
    struct log_mel_handle {
        tap::dsp::log_mel_geometry         geometry;
        std::unique_ptr<tap::dsp::log_mel> fe;
    };
    log_mel_handle* as_log_mel(dsptap_log_mel h) {
        return static_cast<log_mel_handle*>(h);
    }

    // One virtual seam over the three ratio types; the boundary converts to the float golden model.
    struct decimator_base {
        virtual ~decimator_base()                                                                      = default;
        virtual int         taps() const                                                               = 0;
        virtual int         latency() const                                                            = 0;
        virtual void        reset()                                                                    = 0;
        virtual std::size_t outputs_for(std::size_t n) const                                           = 0;
        virtual std::size_t process(const double* in, std::size_t n, double* out, std::size_t max_out) = 0;
    };
    template <std::size_t M>
    struct decimator_impl final : decimator_base {
        explicit decimator_impl(const tap::dsp::decimate_profile& p)
            : dec(p) {}
        int         taps() const override { return static_cast<int>(dec.taps()); }
        int         latency() const override { return static_cast<int>(dec.latency_input_samples()); }
        void        reset() override { dec.reset(); }
        std::size_t outputs_for(std::size_t n) const override { return dec.outputs_for(n); }
        std::size_t process(const double* in, std::size_t n, double* out, std::size_t max_out) override {
            fin.resize(n);
            for (std::size_t i = 0; i < n; ++i) {
                fin[i] = static_cast<float>(in[i]);
            }
            fout.resize(dec.outputs_for(n));
            const std::size_t made = dec.process(fin.data(), n, fout.data());
            const std::size_t keep = made < max_out ? made : max_out;
            for (std::size_t i = 0; i < keep; ++i) {
                out[i] = static_cast<double>(fout[i]);
            }
            return keep;
        }
        tap::dsp::basic_decimator<float, M> dec;
        std::vector<float>                  fin, fout;
    };
    decimator_base* as_decimator(dsptap_decimator h) {
        return static_cast<decimator_base*>(h);
    }

} // namespace

extern "C" {

// -- yin --------------------------------------------------------------------------------------

dsptap_yin dsptap_yin_create(int window, int tau_min, int tau_max) {
    if (tau_min < 2 || tau_min >= tau_max || window < tau_max) {
        return nullptr;
    }
    return new tap::dsp::yin(static_cast<size_t>(window), static_cast<size_t>(tau_min), static_cast<size_t>(tau_max));
}

void dsptap_yin_destroy(dsptap_yin h) {
    delete as_yin(h);
}

int dsptap_yin_set_threshold(dsptap_yin h, double threshold) {
    if (h == nullptr) {
        return -1;
    }
    as_yin(h)->set_threshold(threshold);
    return 0;
}

int dsptap_yin_frame_size(dsptap_yin h) {
    if (h == nullptr) {
        return -1;
    }
    return static_cast<int>(as_yin(h)->frame_size());
}

int dsptap_yin_analyze(dsptap_yin h, const double* frame, double* period, double* aperiodicity) {
    if (h == nullptr || frame == nullptr || period == nullptr || aperiodicity == nullptr) {
        return -1;
    }
    const auto r  = as_yin(h)->analyze(frame);
    *period       = r.period;
    *aperiodicity = r.aperiodicity;
    return 0;
}

int dsptap_yin_track(dsptap_yin h, const double* x, int n, int hop, double* periods, int max_out) {
    if (h == nullptr || x == nullptr || periods == nullptr || hop < 1) {
        return -1;
    }
    auto*     det   = as_yin(h);
    const int frame = static_cast<int>(det->frame_size());
    int       count = 0;
    for (int start = 0; start + frame <= n && count < max_out; start += hop) {
        periods[count++] = det->analyze(x + start).period;
    }
    return count;
}

// -- psola ------------------------------------------------------------------------------------

dsptap_psola dsptap_psola_create(int max_period) {
    if (max_period < 16) {
        return nullptr;
    }
    return new tap::dsp::psola(static_cast<size_t>(max_period));
}

void dsptap_psola_destroy(dsptap_psola h) {
    delete as_psola(h);
}

int dsptap_psola_latency(dsptap_psola h) {
    if (h == nullptr) {
        return -1;
    }
    return static_cast<int>(as_psola(h)->latency());
}

int dsptap_psola_clear(dsptap_psola h) {
    if (h == nullptr) {
        return -1;
    }
    as_psola(h)->clear();
    return 0;
}

int dsptap_psola_process(dsptap_psola h, const double* in, double* out, int n, double period, double ratio) {
    if (h == nullptr || in == nullptr || out == nullptr || n < 0) {
        return -1;
    }
    auto* shifter = as_psola(h);
    for (int i = 0; i < n; ++i) {
        out[i] = shifter->process(in[i], period, ratio);
    }
    return 0;
}

// -- pvoc -------------------------------------------------------------------------------------

dsptap_pvoc dsptap_pvoc_create(int fft_size) {
    if (fft_size < 64 || (fft_size & (fft_size - 1)) != 0) {
        return nullptr;
    }
    return new tap::dsp::pvoc(static_cast<size_t>(fft_size));
}

void dsptap_pvoc_destroy(dsptap_pvoc h) {
    delete as_pvoc(h);
}

int dsptap_pvoc_latency(dsptap_pvoc h) {
    if (h == nullptr) {
        return -1;
    }
    return static_cast<int>(as_pvoc(h)->latency());
}

int dsptap_pvoc_set_formant(dsptap_pvoc h, int on) {
    if (h == nullptr) {
        return -1;
    }
    as_pvoc(h)->set_formant(on != 0);
    return 0;
}

int dsptap_pvoc_clear(dsptap_pvoc h) {
    if (h == nullptr) {
        return -1;
    }
    as_pvoc(h)->clear();
    return 0;
}

int dsptap_pvoc_process(dsptap_pvoc h, const double* in, double* out, int n, double ratio) {
    if (h == nullptr || in == nullptr || out == nullptr || n < 0) {
        return -1;
    }
    auto* shifter = as_pvoc(h);
    for (int i = 0; i < n; ++i) {
        out[i] = shifter->process(in[i], ratio);
    }
    return 0;
}

// -- log_mel ----------------------------------------------------------------------------------

dsptap_log_mel dsptap_log_mel_create(double sample_rate, int frame, int hop, int fft_size, int bands, double fmin_hz,
                                     double fmax_hz, int sqrt_window, double preemphasis) {
    if (frame < 1 || hop < 1 || fft_size < 4 || bands < 1) {
        return nullptr;
    }
    auto h                  = std::make_unique<log_mel_handle>();
    h->geometry.sample_rate = sample_rate;
    h->geometry.frame       = static_cast<size_t>(frame);
    h->geometry.hop         = static_cast<size_t>(hop);
    h->geometry.fft_size    = static_cast<size_t>(fft_size);
    h->geometry.bands       = static_cast<size_t>(bands);
    h->geometry.fmin_hz     = fmin_hz;
    h->geometry.fmax_hz     = fmax_hz;
    h->geometry.window      = sqrt_window != 0 ? tap::dsp::mel_window::sqrt_hann : tap::dsp::mel_window::hann;
    h->geometry.preemphasis = preemphasis;
    if (!h->geometry.valid()) {
        return nullptr;
    }
    h->fe = std::make_unique<tap::dsp::log_mel>(h->geometry);
    return h.release();
}

void dsptap_log_mel_destroy(dsptap_log_mel h) {
    delete as_log_mel(h);
}

int dsptap_log_mel_set_log(dsptap_log_mel h, double floor, double shift, double scale) {
    if (h == nullptr) {
        return -1;
    }
    auto* lm    = as_log_mel(h);
    auto  g     = lm->geometry;
    g.log_floor = floor;
    g.log_shift = shift;
    g.log_scale = scale;
    if (!g.valid()) {
        return -1;
    }
    lm->geometry = g;
    lm->fe       = std::make_unique<tap::dsp::log_mel>(g);
    return 0;
}

int dsptap_log_mel_set_pcen(dsptap_log_mel h, int enabled, double smoother, double alpha, double delta, double power,
                            double epsilon) {
    if (h == nullptr) {
        return -1;
    }
    auto* lm        = as_log_mel(h);
    auto  g         = lm->geometry;
    g.pcen.enabled  = enabled != 0;
    g.pcen.smoother = smoother;
    g.pcen.alpha    = alpha;
    g.pcen.delta    = delta;
    g.pcen.power    = power;
    g.pcen.epsilon  = epsilon;
    if (!g.valid()) {
        return -1;
    }
    lm->geometry = g;
    lm->fe       = std::make_unique<tap::dsp::log_mel>(g);
    return 0;
}

int dsptap_log_mel_reset(dsptap_log_mel h) {
    if (h == nullptr) {
        return -1;
    }
    as_log_mel(h)->fe->reset();
    return 0;
}

int dsptap_log_mel_bands(dsptap_log_mel h) {
    if (h == nullptr) {
        return -1;
    }
    return static_cast<int>(as_log_mel(h)->fe->bands());
}

int dsptap_log_mel_latency(dsptap_log_mel h) {
    if (h == nullptr) {
        return -1;
    }
    return static_cast<int>(as_log_mel(h)->fe->latency_samples());
}

int dsptap_log_mel_contract_version(void) {
    return static_cast<int>(tap::dsp::log_mel_geometry::k_contract_version);
}

int dsptap_log_mel_process(dsptap_log_mel h, const double* x, int n, double* features, int max_frames) {
    if (h == nullptr || x == nullptr || features == nullptr || n < 0 || max_frames < 0) {
        return -1;
    }
    return static_cast<int>(
        as_log_mel(h)->fe->process(x, static_cast<size_t>(n), features, static_cast<size_t>(max_frames)));
}

// -- decimate ---------------------------------------------------------------------------------

dsptap_decimator dsptap_decimator_create(int ratio, int transparent) {
    const auto p = transparent != 0 ? tap::dsp::decimate_profile::transparent() : tap::dsp::decimate_profile::economy();
    switch (ratio) {
    case 2:
        return new decimator_impl<2>(p);
    case 3:
        return new decimator_impl<3>(p);
    case 6:
        return new decimator_impl<6>(p);
    default:
        return nullptr;
    }
}

void dsptap_decimator_destroy(dsptap_decimator h) {
    delete as_decimator(h);
}

int dsptap_decimator_taps(dsptap_decimator h) {
    return h == nullptr ? -1 : as_decimator(h)->taps();
}

int dsptap_decimator_latency(dsptap_decimator h) {
    return h == nullptr ? -1 : as_decimator(h)->latency();
}

int dsptap_decimator_reset(dsptap_decimator h) {
    if (h == nullptr) {
        return -1;
    }
    as_decimator(h)->reset();
    return 0;
}

int dsptap_decimator_outputs_for(dsptap_decimator h, int n) {
    if (h == nullptr || n < 0) {
        return -1;
    }
    return static_cast<int>(as_decimator(h)->outputs_for(static_cast<size_t>(n)));
}

int dsptap_decimator_process(dsptap_decimator h, const double* in, int n, double* out, int max_out) {
    if (h == nullptr || in == nullptr || out == nullptr || n < 0 || max_out < 0) {
        return -1;
    }
    return static_cast<int>(as_decimator(h)->process(in, static_cast<size_t>(n), out, static_cast<size_t>(max_out)));
}

} // extern "C"
