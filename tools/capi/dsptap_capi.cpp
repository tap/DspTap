/// @file dsptap_capi.cpp
/// @brief C ABI over the DspTap primitives — see dsptap_capi.h.
// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.

#include "dsptap_capi.h"

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

} // extern "C"
