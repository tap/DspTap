// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// ONE translation unit, built twice into two loadable modules by
// tests/CMakeLists.txt (hosted, non-Windows):
//
//   image A  the configured build: tap::dsp's own defines, so the tag is
//            fft_split_radix on linux and fft_vdsp on the macOS leg;
//   image B  -DTAP_DSP_FFT_CMSIS with tests/abi/stub/ first on the include
//            path, so fft.h selects a stub "CMSIS" engine with a different
//            object layout and the tag is fft_cmsis.
//
// tests/test_fft_abi_tag.cpp dlopens A then B into one process with
// RTLD_GLOBAL and asks each image, through the extern "C" probes below,
// what its embedders see. Everything the test needs is exported as plain C
// functions so the host never names a C++ type from either image; the C++
// symbols under test (the embedders' member functions) are the ordinary
// weak, default-visibility template instantiations every consumer image
// exports, which is exactly the shape audit F4 is about.
//
// Two embedders of basic_real_fft<float>, identical except for where they
// are defined: untagged_embedder in plain tap::dsp (what a consumer class
// that has not adopted the tag is — MuTap's today), tagged_embedder inside
// inline namespace TAP_DSP_FFT_ABI (what pvoc.h and log_mel.h do since
// Stage 4). Their member functions are kept out of line and called through
// a pointer read from a volatile, so each call is a real dynamic-symbol
// call the loader resolves, not something the optimizer folded; that is
// what makes the coalescing observable when it happens.

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numbers>
#include <vector>

#include "tap/dsp/fft.h"
#include "tap/dsp/log_mel.h"
#include "tap/dsp/pvoc.h"

#if defined(__clang__)
#define TAP_DSP_ABI_IMAGE_NOINLINE [[gnu::noinline]]
#elif defined(__GNUC__)
#define TAP_DSP_ABI_IMAGE_NOINLINE [[gnu::noinline, gnu::noipa]]
#else
#define TAP_DSP_ABI_IMAGE_NOINLINE
#endif

namespace tap::dsp {

    /// An embedder OUTSIDE the tag.
    template <typename S>
    struct untagged_embedder {
        basic_real_fft<S> fft;
        std::size_t       n;
        explicit untagged_embedder(std::size_t k)
            : fft(k)
            , n(k) {}
        TAP_DSP_ABI_IMAGE_NOINLINE const char* which() const noexcept { return k_real_fft_abi_tag; }
        TAP_DSP_ABI_IMAGE_NOINLINE std::size_t layout() const noexcept { return sizeof(*this); }
    };

    inline namespace TAP_DSP_FFT_ABI {
        /// The same embedder INSIDE the tag.
        template <typename S>
        struct tagged_embedder {
            basic_real_fft<S> fft;
            std::size_t       n;
            explicit tagged_embedder(std::size_t k)
                : fft(k)
                , n(k) {}
            TAP_DSP_ABI_IMAGE_NOINLINE const char* which() const noexcept { return k_real_fft_abi_tag; }
            TAP_DSP_ABI_IMAGE_NOINLINE std::size_t layout() const noexcept { return sizeof(*this); }
        };
    } // namespace TAP_DSP_FFT_ABI

} // namespace tap::dsp

namespace {

    /// Calls a member through a pointer the optimizer cannot see through.
    template <typename Object, typename Result>
    Result call_laundered(const Object& object, Result (Object::*member)() const noexcept) {
        volatile auto laundered = member;
        auto          pointer   = laundered;
        return (object.*pointer)();
    }

    constexpr std::size_t k_n = 512;

} // namespace

extern "C" {

const char* tap_dsp_abi_image_built_as() {
    return TAP_DSP_FFT_ABI_NAME;
}

// The untagged embedder: what() it reports and the layout the bound
// member function sees, against the layout this image really has.
const char* tap_dsp_abi_image_untagged_which() {
    const tap::dsp::untagged_embedder<float> e(k_n);
    return call_laundered(e, &tap::dsp::untagged_embedder<float>::which);
}
std::size_t tap_dsp_abi_image_untagged_layout_seen() {
    const tap::dsp::untagged_embedder<float> e(k_n);
    return call_laundered(e, &tap::dsp::untagged_embedder<float>::layout);
}
std::size_t tap_dsp_abi_image_untagged_layout_real() {
    return sizeof(tap::dsp::untagged_embedder<float>);
}

// The tagged embedder, likewise.
const char* tap_dsp_abi_image_tagged_which() {
    const tap::dsp::tagged_embedder<float> e(k_n);
    return call_laundered(e, &tap::dsp::tagged_embedder<float>::which);
}
std::size_t tap_dsp_abi_image_tagged_layout_seen() {
    const tap::dsp::tagged_embedder<float> e(k_n);
    return call_laundered(e, &tap::dsp::tagged_embedder<float>::layout);
}
std::size_t tap_dsp_abi_image_tagged_layout_real() {
    return sizeof(tap::dsp::tagged_embedder<float>);
}

// The shipping embedders, run for real: pvoc32 shifts a tone and log_mel32
// analyses one, each through its own image's basic_real_fft<float>. With a
// coalesced member (the hazard) the other image's code would run over this
// image's layout; a finite, sane checksum from each image after both are
// loaded is the observation.
double tap_dsp_abi_image_pvoc_checksum() {
    tap::dsp::pvoc32 shifter(k_n);
    double           sum = 0.0;
    for (std::size_t i = 0; i < 8 * k_n; ++i) {
        const float x =
            0.5f * static_cast<float>(std::sin(2.0 * std::numbers::pi * 440.0 * static_cast<double>(i) / 48000.0));
        sum += static_cast<double>(shifter.process(x, 1.5f));
    }
    return sum;
}
double tap_dsp_abi_image_log_mel_checksum() {
    tap::dsp::log_mel_geometry g;
    tap::dsp::log_mel32        front_end(g);
    std::vector<float>         hop(g.hop);
    std::vector<float>         features(g.bands);
    double                     sum = 0.0;
    for (std::size_t frame = 0; frame < 8; ++frame) {
        for (std::size_t i = 0; i < g.hop; ++i) {
            const std::size_t t = frame * g.hop + i;
            hop[i]              = 0.25f
                     * static_cast<float>(
                         std::sin(2.0 * std::numbers::pi * 1000.0 * static_cast<double>(t) / g.sample_rate));
        }
        front_end.process_hop(hop.data(), features.data());
        for (const float f : features) {
            sum += static_cast<double>(f);
        }
    }
    return sum;
}
std::size_t tap_dsp_abi_image_pvoc_layout() {
    return sizeof(tap::dsp::pvoc32);
}
std::size_t tap_dsp_abi_image_log_mel_layout() {
    return sizeof(tap::dsp::log_mel32);
}

} // extern "C"
