/// @file dsptap_capi.cpp
/// @brief C ABI over the DspTap primitives — see dsptap_capi.h.
// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.

#include "dsptap_capi.h"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <new>
#include <type_traits>
#include <vector>

#include "tap/dsp/decimate.h"
#include "tap/dsp/fft.h"
#include "tap/dsp/log_mel.h"
#include "tap/dsp/psola.h"
#include "tap/dsp/pvoc.h"
#include "tap/dsp/sample_traits.h"
#include "tap/dsp/yin.h"

// The handle types the header declares incomplete (typedef struct dsptap_*_s* dsptap_*), defined
// here and only here. The two polymorphic seams ARE the handle types, so a create function
// returns its derived object through an implicit upcast and every other entry point uses the
// handle directly: no void*, no cast anywhere at the boundary. The other four wrap their object
// by value.

// One virtual seam over the six FFT profiles (four sample types, two scaling policies for the
// fixed-point ones), implemented by fft_impl below. The double entry points run the double engine
// directly; the float profile stages through a float buffer allocated at construction; the
// fixed-point profiles stage through a native buffer, quantizing at the boundary by the
// substrate's round_sat and reading back as fractions times 2^e, so every transform stays
// allocation-free. The raw entry points hand the native buffer straight to the header's in-place
// transforms, with no conversion, so the notebooks measure the profile's own arithmetic; the seam
// returns the exponent (0 for the floating profiles) and the C entry points decide how to present
// it (the *_exp variants write it out; the legacy raw pair refuses fixed-point handles).
struct dsptap_fft_s {
    dsptap_fft_s()                                               = default;
    dsptap_fft_s(const dsptap_fft_s&)                            = delete;
    dsptap_fft_s& operator=(const dsptap_fft_s&)                 = delete;
    virtual ~dsptap_fft_s()                                      = default;
    virtual int  size() const noexcept                           = 0;
    virtual int  profile() const noexcept                        = 0;
    virtual bool is_fixed_point() const noexcept                 = 0;
    virtual int  sample_bytes() const noexcept                   = 0;
    virtual int  fixed_scaling_exponent() const noexcept         = 0;
    virtual int  forward(const double* in, double* out) noexcept = 0; ///< returns e
    virtual int  inverse(const double* in, double* out) noexcept = 0; ///< returns e
    virtual int  forward_inplace_raw(void* data) noexcept        = 0; ///< returns e
    virtual int  inverse_inplace_raw(void* data) noexcept        = 0; ///< returns e
};

/// One virtual seam over the three decimation ratios (decimator_impl below).
struct dsptap_decimator_s {
    dsptap_decimator_s()                                                                                    = default;
    dsptap_decimator_s(const dsptap_decimator_s&)                                                           = delete;
    dsptap_decimator_s& operator=(const dsptap_decimator_s&)                                                = delete;
    virtual ~dsptap_decimator_s()                                                                           = default;
    virtual int         taps() const noexcept                                                               = 0;
    virtual int         latency() const noexcept                                                            = 0;
    virtual void        reset() noexcept                                                                    = 0;
    virtual std::size_t outputs_for(std::size_t n) const noexcept                                           = 0;
    virtual std::size_t process(const double* in, std::size_t n, double* out, std::size_t max_out) noexcept = 0;
};

struct dsptap_yin_s {
    explicit dsptap_yin_s(std::size_t window, std::size_t tau_min, std::size_t tau_max)
        : impl(window, tau_min, tau_max) {}
    tap::dsp::yin impl;
};

struct dsptap_psola_s {
    explicit dsptap_psola_s(std::size_t max_period)
        : impl(max_period) {}
    tap::dsp::psola impl;
};

struct dsptap_pvoc_s {
    explicit dsptap_pvoc_s(std::size_t fft_size)
        : impl(fft_size) {}
    tap::dsp::pvoc impl;
};

/// log_mel's geometry is fixed at construction; the C ABI setters rebuild.
struct dsptap_log_mel_s {
    tap::dsp::log_mel_geometry         geometry;
    std::unique_ptr<tap::dsp::log_mel> fe;
};

// The point of the typed handles: six distinct pointer types, so passing one primitive's handle to
// another's function is a compile error in C++ and a diagnosed incompatible-pointer conversion in C.
static_assert(!std::is_same_v<dsptap_fft, dsptap_yin> && !std::is_same_v<dsptap_yin, dsptap_psola>
                  && !std::is_same_v<dsptap_psola, dsptap_pvoc> && !std::is_same_v<dsptap_pvoc, dsptap_log_mel>
                  && !std::is_same_v<dsptap_log_mel, dsptap_decimator> && !std::is_same_v<dsptap_decimator, dsptap_fft>
                  && !std::is_convertible_v<dsptap_yin, dsptap_pvoc> && !std::is_convertible_v<void*, dsptap_fft>,
              "each primitive has its own opaque handle type");

namespace {

    template <typename>
    inline constexpr bool k_always_false = false;

    // The DSPTAP_FFT_PROFILE_* a (sample type, scaling policy) pair reports. Exhaustive on
    // purpose: an instantiation that is not mapped here fails to compile instead of reporting a
    // neighbour's profile.
    template <typename Sample, typename Policy>
    constexpr int profile_of() {
        constexpr bool bfp = std::is_same_v<Policy, tap::dsp::scaling::block_floating>;
        if constexpr (std::is_same_v<Sample, double>) {
            static_assert(!bfp, "the floating profiles have no scaling policy");
            return DSPTAP_FFT_PROFILE_DOUBLE;
        }
        else if constexpr (std::is_same_v<Sample, float>) {
            static_assert(!bfp, "the floating profiles have no scaling policy");
            return DSPTAP_FFT_PROFILE_FLOAT;
        }
        else if constexpr (std::is_same_v<Sample, std::int16_t>) {
            return bfp ? DSPTAP_FFT_PROFILE_Q15_BFP : DSPTAP_FFT_PROFILE_Q15;
        }
        else if constexpr (std::is_same_v<Sample, std::int32_t>) {
            return bfp ? DSPTAP_FFT_PROFILE_Q31_BFP : DSPTAP_FFT_PROFILE_Q31;
        }
        else {
            static_assert(k_always_false<Sample>, "no DSPTAP_FFT_PROFILE_* for this sample type");
        }
    }

    // Policy is basic_real_fft's own second argument: the selected engine for the floating
    // profiles (so fft_impl<float> holds exactly tap::dsp::real_fft32, the type pvoc and
    // log_mel in this same library hold — not the pre-Stage-4 <float, scaling::fixed> spelling,
    // a second type with identical code; 35a/F3) and the scaling policy for the fixed-point
    // ones. profile_of<> reads the policy: block_floating names the BFP profiles.
    template <typename Sample, typename Policy = tap::dsp::detail::default_real_fft_policy_t<Sample>>
    struct fft_impl final : dsptap_fft_s {
        using fft_type                      = tap::dsp::basic_real_fft<Sample, Policy>;
        static constexpr bool k_fixed_point = tap::dsp::sample_traits<Sample>::k_is_fixed_point;

        explicit fft_impl(std::size_t n)
            : fft(n)
            , buf(n, Sample(0)) {}
        int  size() const noexcept override { return static_cast<int>(fft.size()); }
        int  profile() const noexcept override { return profile_of<Sample, Policy>(); }
        bool is_fixed_point() const noexcept override { return k_fixed_point; }
        int  sample_bytes() const noexcept override { return static_cast<int>(sizeof(Sample)); }
        int  fixed_scaling_exponent() const noexcept override {
            if constexpr (k_fixed_point) {
                return fft_type::fixed_scaling_exponent(fft.size());
            }
            else {
                return 0;
            }
        }
        int forward(const double* in, double* out) noexcept override {
            if constexpr (std::is_same_v<Sample, double>) {
                fft.forward(in, out);
                return 0;
            }
            else if constexpr (k_fixed_point) {
                quantize(in);
                const int e = fft.forward_inplace(buf.data());
                dequantize(out, e);
                return e;
            }
            else {
                narrow(in);
                fft.forward_inplace(buf.data());
                widen(out, Sample(1));
                return 0;
            }
        }
        int inverse(const double* in, double* out) noexcept override {
            if constexpr (std::is_same_v<Sample, double>) {
                fft.inverse(in, out);
                return 0;
            }
            else if constexpr (k_fixed_point) {
                quantize(in);
                const int e = fft.inverse_inplace(buf.data()); // no 2/N: the exponent is the scale
                dequantize(out, e);
                return e;
            }
            else {
                narrow(in);
                fft.inverse_inplace(buf.data());
                widen(out, Sample(2) / static_cast<Sample>(fft.size())); // inverse()'s scale, in Sample
                return 0;
            }
        }
        int forward_inplace_raw(void* data) noexcept override {
            if constexpr (k_fixed_point) {
                return fft.forward_inplace(static_cast<Sample*>(data));
            }
            else {
                fft.forward_inplace(static_cast<Sample*>(data));
                return 0;
            }
        }
        int inverse_inplace_raw(void* data) noexcept override {
            if constexpr (k_fixed_point) {
                return fft.inverse_inplace(static_cast<Sample*>(data));
            }
            else {
                fft.inverse_inplace(static_cast<Sample*>(data));
                return 0;
            }
        }

        // The floating boundary: plain casts.
        void narrow(const double* in) noexcept {
            for (std::size_t i = 0; i < buf.size(); ++i) {
                buf[i] = static_cast<Sample>(in[i]);
            }
        }
        void widen(double* out, Sample scale) noexcept {
            for (std::size_t i = 0; i < buf.size(); ++i) {
                out[i] = static_cast<double>(buf[i] * scale);
            }
        }
        // The fixed-point boundary: fractions of full scale to the Q format by the substrate's
        // round_sat (half away from zero, saturating), and back as fractions times 2^e (an exact
        // power-of-two scale on an exactly representable integer).
        void quantize(const double* in) noexcept {
            if constexpr (k_fixed_point) {
                constexpr double scale =
                    static_cast<double>(std::int64_t{1} << tap::dsp::sample_traits<Sample>::k_sample_frac_bits);
                for (std::size_t i = 0; i < buf.size(); ++i) {
                    buf[i] = tap::dsp::detail::round_sat<Sample>(in[i] * scale);
                }
            }
        }
        void dequantize(double* out, int exponent) noexcept {
            if constexpr (k_fixed_point) {
                const double scale = std::ldexp(1.0, exponent - tap::dsp::sample_traits<Sample>::k_sample_frac_bits);
                for (std::size_t i = 0; i < buf.size(); ++i) {
                    out[i] = static_cast<double>(buf[i]) * scale;
                }
            }
        }

        fft_type            fft;
        std::vector<Sample> buf;
    };

    // The size gate, per profile, is the header's own: basic_real_fft<...>::supports_size(n),
    // the power-of-two interval [k_min_size, k_max_size] of the engine that profile runs on
    // (fft.h, Stage 4). The header's constructor states the same range as a precondition
    // (TAP_EXPECTS, a debug assertion that a noexcept C entry point could not turn into NULL
    // and that evaluates to nothing in a release build), so this call is the release-mode
    // check the header names as mandatory wherever N comes from outside: the floating
    // profiles' range is the selected engine's (split-radix 4 … 2^30, vDSP 4 … 2^20, CMSIS
    // 32 … 4096), the fixed-point profiles' 4 … 65536.
    template <typename Impl>
    dsptap_fft_s* make_fft_if_supported(int size) {
        if (size <= 0 || !Impl::fft_type::supports_size(static_cast<std::size_t>(size))) {
            return nullptr;
        }
        return new Impl(static_cast<std::size_t>(size));
    }

    // The double boundary stages through float buffers sized here, once: process() walks the
    // input in k_block-sample blocks (decimate.h's process() is chunking-invariant, bit for bit),
    // so no call allocates, whatever n is.
    template <std::size_t M>
    struct decimator_impl final : dsptap_decimator_s {
        static constexpr std::size_t k_block = 4096;

        explicit decimator_impl(const tap::dsp::decimate_profile& p)
            : dec(p)
            , fin(k_block)
            , fout(k_block / M + 1) {}
        int         taps() const noexcept override { return static_cast<int>(dec.taps()); }
        int         latency() const noexcept override { return static_cast<int>(dec.latency_input_samples()); }
        void        reset() noexcept override { dec.reset(); }
        std::size_t outputs_for(std::size_t n) const noexcept override { return dec.outputs_for(n); }
        std::size_t process(const double* in, std::size_t n, double* out, std::size_t max_out) noexcept override {
            std::size_t kept = 0;
            for (std::size_t done = 0; done < n;) {
                const std::size_t len = (n - done) < k_block ? (n - done) : k_block;
                for (std::size_t i = 0; i < len; ++i) {
                    fin[i] = static_cast<float>(in[done + i]);
                }
                // outputs_for(len) <= ceil(len / M) <= k_block / M + 1 == fout.size().
                const std::size_t made = dec.process(fin.data(), len, fout.data());
                for (std::size_t i = 0; i < made && kept < max_out; ++i) {
                    out[kept++] = static_cast<double>(fout[i]);
                }
                done += len;
            }
            return kept;
        }
        tap::dsp::basic_decimator<float, M> dec;
        std::vector<float>                  fin;
        std::vector<float>                  fout;
    };

    /// Rebuild h's front end at geometry g, or report -1 and leave h untouched: validate, build
    /// the new object, and only then commit both fields (neither assignment can throw).
    int rebuild_log_mel(dsptap_log_mel h, const tap::dsp::log_mel_geometry& g) noexcept {
        if (!tap::dsp::log_mel::supports_geometry(g)) { // valid() and the FFT engine's size range
            return -1;
        }
        try {
            auto fe     = std::make_unique<tap::dsp::log_mel>(g);
            h->geometry = g;
            h->fe       = std::move(fe);
            return 0;
        }
        catch (...) {
            return -1;
        }
    }

} // namespace

extern "C" {

// -- fft --------------------------------------------------------------------------------------

dsptap_fft dsptap_fft_create(int size, int profile) DSPTAP_NOEXCEPT {
    try {
        // Each profile is gated by its own class's supports_size (make_fft_if_supported).
        switch (profile) {
        case DSPTAP_FFT_PROFILE_DOUBLE:
            return make_fft_if_supported<fft_impl<double>>(size);
        case DSPTAP_FFT_PROFILE_FLOAT:
            return make_fft_if_supported<fft_impl<float>>(size);
        case DSPTAP_FFT_PROFILE_Q15:
            return make_fft_if_supported<fft_impl<std::int16_t>>(size);
        case DSPTAP_FFT_PROFILE_Q31:
            return make_fft_if_supported<fft_impl<std::int32_t>>(size);
        case DSPTAP_FFT_PROFILE_Q15_BFP:
            return make_fft_if_supported<fft_impl<std::int16_t, tap::dsp::scaling::block_floating>>(size);
        case DSPTAP_FFT_PROFILE_Q31_BFP:
            return make_fft_if_supported<fft_impl<std::int32_t, tap::dsp::scaling::block_floating>>(size);
        default: // anything else is a bad argument
            return nullptr;
        }
    }
    catch (...) {
        return nullptr;
    }
}

void dsptap_fft_destroy(dsptap_fft h) DSPTAP_NOEXCEPT {
    delete h;
}

int dsptap_fft_size(dsptap_fft h) DSPTAP_NOEXCEPT {
    return h == nullptr ? -1 : h->size();
}

int dsptap_fft_num_bins(dsptap_fft h) DSPTAP_NOEXCEPT {
    return h == nullptr ? -1 : h->size() / 2 + 1;
}

int dsptap_fft_profile(dsptap_fft h) DSPTAP_NOEXCEPT {
    return h == nullptr ? -1 : h->profile();
}

int dsptap_fft_sample_bytes(dsptap_fft h) DSPTAP_NOEXCEPT {
    return h == nullptr ? -1 : h->sample_bytes();
}

int dsptap_fft_fixed_scaling_exponent(dsptap_fft h) DSPTAP_NOEXCEPT {
    return h == nullptr ? -1 : h->fixed_scaling_exponent();
}

const char* dsptap_fft_backend(void) DSPTAP_NOEXCEPT {
#if defined(TAP_DSP_FFT_ACCELERATE)
    return "accelerate";
#elif defined(TAP_DSP_FFT_CMSIS)
    return "cmsis";
#else
    // The C++20 split-radix engine (fft/split_radix.h): what basic_real_fft<float>
    // runs where no backend is selected. Named "ooura" until Stage 2c, when the
    // vendored C it is bit-identical to left the shipping tree; the engine
    // is named for what it is (Decision D7), the notices carry the attribution.
    return "split_radix";
#endif
}

int dsptap_fft_forward(dsptap_fft h, const double* in, double* out) DSPTAP_NOEXCEPT {
    if (h == nullptr || in == nullptr || out == nullptr) {
        return -1;
    }
    (void)h->forward(in, out); // out already carries 2^e
    return 0;
}

int dsptap_fft_inverse(dsptap_fft h, const double* in, double* out) DSPTAP_NOEXCEPT {
    if (h == nullptr || in == nullptr || out == nullptr) {
        return -1;
    }
    (void)h->inverse(in, out);
    return 0;
}

int dsptap_fft_forward_exp(dsptap_fft h, const double* in, double* out, int* exponent) DSPTAP_NOEXCEPT {
    if (h == nullptr || in == nullptr || out == nullptr || exponent == nullptr) {
        return -1;
    }
    *exponent = h->forward(in, out);
    return 0;
}

int dsptap_fft_inverse_exp(dsptap_fft h, const double* in, double* out, int* exponent) DSPTAP_NOEXCEPT {
    if (h == nullptr || in == nullptr || out == nullptr || exponent == nullptr) {
        return -1;
    }
    *exponent = h->inverse(in, out);
    return 0;
}

int dsptap_fft_forward_inplace_raw(dsptap_fft h, void* data) DSPTAP_NOEXCEPT {
    if (h == nullptr || data == nullptr || h->is_fixed_point()) {
        return -1; // a fixed-point transform's exponent has nowhere to go here: use the _exp variant
    }
    (void)h->forward_inplace_raw(data); // e == 0 for the floating profiles
    return 0;
}

int dsptap_fft_inverse_inplace_raw(dsptap_fft h, void* data) DSPTAP_NOEXCEPT {
    if (h == nullptr || data == nullptr || h->is_fixed_point()) {
        return -1;
    }
    (void)h->inverse_inplace_raw(data);
    return 0;
}

int dsptap_fft_forward_inplace_raw_exp(dsptap_fft h, void* data, int* exponent) DSPTAP_NOEXCEPT {
    if (h == nullptr || data == nullptr || exponent == nullptr) {
        return -1;
    }
    *exponent = h->forward_inplace_raw(data);
    return 0;
}

int dsptap_fft_inverse_inplace_raw_exp(dsptap_fft h, void* data, int* exponent) DSPTAP_NOEXCEPT {
    if (h == nullptr || data == nullptr || exponent == nullptr) {
        return -1;
    }
    *exponent = h->inverse_inplace_raw(data);
    return 0;
}

// -- yin --------------------------------------------------------------------------------------

dsptap_yin dsptap_yin_create(int window, int tau_min, int tau_max) DSPTAP_NOEXCEPT {
    if (tau_min < 2 || tau_min >= tau_max || window < tau_max) {
        return nullptr; // yin.h's @pre, which the header only asserts
    }
    if (window > std::numeric_limits<int>::max() - tau_max) {
        return nullptr; // yin.h's @pre: window + tau_max (frame_size(), an int sum) must fit an int
    }
    try {
        return new dsptap_yin_s(static_cast<std::size_t>(window), static_cast<std::size_t>(tau_min),
                                static_cast<std::size_t>(tau_max));
    }
    catch (...) {
        return nullptr;
    }
}

void dsptap_yin_destroy(dsptap_yin h) DSPTAP_NOEXCEPT {
    delete h;
}

int dsptap_yin_set_threshold(dsptap_yin h, double threshold) DSPTAP_NOEXCEPT {
    if (h == nullptr) {
        return -1;
    }
    h->impl.set_threshold(threshold);
    return 0;
}

int dsptap_yin_frame_size(dsptap_yin h) DSPTAP_NOEXCEPT {
    if (h == nullptr) {
        return -1;
    }
    return static_cast<int>(h->impl.frame_size());
}

int dsptap_yin_analyze(dsptap_yin h, const double* frame, double* period, double* aperiodicity) DSPTAP_NOEXCEPT {
    if (h == nullptr || frame == nullptr || period == nullptr || aperiodicity == nullptr) {
        return -1;
    }
    const auto r  = h->impl.analyze(frame);
    *period       = r.period;
    *aperiodicity = r.aperiodicity;
    return 0;
}

int dsptap_yin_track(dsptap_yin h, const double* x, int n, int hop, double* periods, int max_out) DSPTAP_NOEXCEPT {
    if (h == nullptr || x == nullptr || periods == nullptr || hop < 1) {
        return -1;
    }
    // Frame starts in 64-bit so start + hop cannot overflow for any int n and hop.
    const auto frame = static_cast<std::int64_t>(h->impl.frame_size());
    int        count = 0;
    for (std::int64_t start = 0; start + frame <= n && count < max_out; start += hop) {
        periods[count++] = h->impl.analyze(x + start).period;
    }
    return count;
}

// -- psola ------------------------------------------------------------------------------------

dsptap_psola dsptap_psola_create(int max_period) DSPTAP_NOEXCEPT {
    if (max_period < 16 || max_period >= (1 << 26)) {
        return nullptr; // psola.h's @pre: 16 <= max_period < 2^26
    }
    try {
        return new dsptap_psola_s(static_cast<std::size_t>(max_period));
    }
    catch (...) {
        return nullptr;
    }
}

void dsptap_psola_destroy(dsptap_psola h) DSPTAP_NOEXCEPT {
    delete h;
}

int dsptap_psola_latency(dsptap_psola h) DSPTAP_NOEXCEPT {
    if (h == nullptr) {
        return -1;
    }
    return static_cast<int>(h->impl.latency());
}

int dsptap_psola_clear(dsptap_psola h) DSPTAP_NOEXCEPT {
    if (h == nullptr) {
        return -1;
    }
    h->impl.clear();
    return 0;
}

int dsptap_psola_process(dsptap_psola h, const double* in, double* out, int n, double period,
                         double ratio) DSPTAP_NOEXCEPT {
    if (h == nullptr || in == nullptr || out == nullptr || n < 0) {
        return -1;
    }
    for (int i = 0; i < n; ++i) {
        out[i] = h->impl.process(in[i], period, ratio);
    }
    return 0;
}

// -- pvoc -------------------------------------------------------------------------------------

dsptap_pvoc dsptap_pvoc_create(int fft_size) DSPTAP_NOEXCEPT {
    // pvoc.h's @pre, as the class states it: a power of two in [k_min_size, k_max_size], the
    // class's own [64, 2^28] intersected with its FFT engine's range (for this double profile,
    // split-radix on every build, the intersection is [64, 2^28]).
    if (fft_size <= 0 || !tap::dsp::pvoc::supports_size(static_cast<std::size_t>(fft_size))) {
        return nullptr;
    }
    try {
        return new dsptap_pvoc_s(static_cast<std::size_t>(fft_size));
    }
    catch (...) {
        return nullptr;
    }
}

void dsptap_pvoc_destroy(dsptap_pvoc h) DSPTAP_NOEXCEPT {
    delete h;
}

int dsptap_pvoc_latency(dsptap_pvoc h) DSPTAP_NOEXCEPT {
    if (h == nullptr) {
        return -1;
    }
    return static_cast<int>(h->impl.latency());
}

int dsptap_pvoc_set_formant(dsptap_pvoc h, int on) DSPTAP_NOEXCEPT {
    if (h == nullptr) {
        return -1;
    }
    h->impl.set_formant(on != 0);
    return 0;
}

int dsptap_pvoc_clear(dsptap_pvoc h) DSPTAP_NOEXCEPT {
    if (h == nullptr) {
        return -1;
    }
    h->impl.clear();
    return 0;
}

int dsptap_pvoc_process(dsptap_pvoc h, const double* in, double* out, int n, double ratio) DSPTAP_NOEXCEPT {
    if (h == nullptr || in == nullptr || out == nullptr || n < 0) {
        return -1;
    }
    for (int i = 0; i < n; ++i) {
        out[i] = h->impl.process(in[i], ratio);
    }
    return 0;
}

// -- log_mel ----------------------------------------------------------------------------------

dsptap_log_mel dsptap_log_mel_create(double sample_rate, int frame, int hop, int fft_size, int bands, double fmin_hz,
                                     double fmax_hz, int sqrt_window, double preemphasis) DSPTAP_NOEXCEPT {
    if (frame < 1 || hop < 1 || fft_size < 4 || bands < 1) {
        return nullptr;
    }
    tap::dsp::log_mel_geometry g;
    g.sample_rate = sample_rate;
    g.frame       = static_cast<std::size_t>(frame);
    g.hop         = static_cast<std::size_t>(hop);
    g.fft_size    = static_cast<std::size_t>(fft_size);
    g.bands       = static_cast<std::size_t>(bands);
    g.fmin_hz     = fmin_hz;
    g.fmax_hz     = fmax_hz;
    g.window      = sqrt_window != 0 ? tap::dsp::mel_window::sqrt_hann : tap::dsp::mel_window::hann;
    g.preemphasis = preemphasis;
    // log_mel.h's @pre: valid() and the FFT engine's size range (for this double profile,
    // split-radix on every build, 4 … 2^30).
    if (!tap::dsp::log_mel::supports_geometry(g)) {
        return nullptr;
    }
    try {
        auto h      = std::make_unique<dsptap_log_mel_s>();
        h->geometry = g;
        h->fe       = std::make_unique<tap::dsp::log_mel>(g);
        return h.release();
    }
    catch (...) {
        return nullptr;
    }
}

void dsptap_log_mel_destroy(dsptap_log_mel h) DSPTAP_NOEXCEPT {
    delete h;
}

int dsptap_log_mel_set_log(dsptap_log_mel h, double floor, double shift, double scale) DSPTAP_NOEXCEPT {
    if (h == nullptr) {
        return -1;
    }
    auto g      = h->geometry;
    g.log_floor = floor;
    g.log_shift = shift;
    g.log_scale = scale;
    return rebuild_log_mel(h, g);
}

int dsptap_log_mel_set_pcen(dsptap_log_mel h, int enabled, double smoother, double alpha, double delta, double power,
                            double epsilon) DSPTAP_NOEXCEPT {
    if (h == nullptr) {
        return -1;
    }
    auto g          = h->geometry;
    g.pcen.enabled  = enabled != 0;
    g.pcen.smoother = smoother;
    g.pcen.alpha    = alpha;
    g.pcen.delta    = delta;
    g.pcen.power    = power;
    g.pcen.epsilon  = epsilon;
    return rebuild_log_mel(h, g);
}

int dsptap_log_mel_reset(dsptap_log_mel h) DSPTAP_NOEXCEPT {
    if (h == nullptr) {
        return -1;
    }
    h->fe->reset();
    return 0;
}

int dsptap_log_mel_bands(dsptap_log_mel h) DSPTAP_NOEXCEPT {
    if (h == nullptr) {
        return -1;
    }
    return static_cast<int>(h->fe->bands());
}

int dsptap_log_mel_latency(dsptap_log_mel h) DSPTAP_NOEXCEPT {
    if (h == nullptr) {
        return -1;
    }
    return static_cast<int>(h->fe->latency_samples());
}

int dsptap_log_mel_contract_version(void) DSPTAP_NOEXCEPT {
    return static_cast<int>(tap::dsp::log_mel_geometry::k_contract_version);
}

int dsptap_log_mel_process(dsptap_log_mel h, const double* x, int n, double* features, int max_frames) DSPTAP_NOEXCEPT {
    if (h == nullptr || x == nullptr || features == nullptr || n < 0 || max_frames < 0) {
        return -1;
    }
    return static_cast<int>(
        h->fe->process(x, static_cast<std::size_t>(n), features, static_cast<std::size_t>(max_frames)));
}

// -- decimate ---------------------------------------------------------------------------------

dsptap_decimator dsptap_decimator_create(int ratio, int transparent) DSPTAP_NOEXCEPT {
    const auto p = transparent != 0 ? tap::dsp::decimate_profile::transparent() : tap::dsp::decimate_profile::economy();
    try {
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
    catch (...) {
        return nullptr;
    }
}

void dsptap_decimator_destroy(dsptap_decimator h) DSPTAP_NOEXCEPT {
    delete h;
}

int dsptap_decimator_taps(dsptap_decimator h) DSPTAP_NOEXCEPT {
    return h == nullptr ? -1 : h->taps();
}

int dsptap_decimator_latency(dsptap_decimator h) DSPTAP_NOEXCEPT {
    return h == nullptr ? -1 : h->latency();
}

int dsptap_decimator_reset(dsptap_decimator h) DSPTAP_NOEXCEPT {
    if (h == nullptr) {
        return -1;
    }
    h->reset();
    return 0;
}

int dsptap_decimator_outputs_for(dsptap_decimator h, int n) DSPTAP_NOEXCEPT {
    if (h == nullptr || n < 0) {
        return -1;
    }
    return static_cast<int>(h->outputs_for(static_cast<std::size_t>(n)));
}

int dsptap_decimator_process(dsptap_decimator h, const double* in, int n, double* out, int max_out) DSPTAP_NOEXCEPT {
    if (h == nullptr || in == nullptr || out == nullptr || n < 0 || max_out < 0) {
        return -1;
    }
    return static_cast<int>(h->process(in, static_cast<std::size_t>(n), out, static_cast<std::size_t>(max_out)));
}

} // extern "C"
