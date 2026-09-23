// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// THE TWO-IMAGE TEST OF THE ABI TAG (Stage 4, audit F4; asked for by the
// 35a review, whose experiment this reproduces): two loadable modules built
// from ONE translation unit (tests/abi/abi_tag_image.cpp) with different
// float-default engines and therefore different tags — image A the
// configured build (fft_split_radix on linux, fft_vdsp on macOS), image B a
// stub "CMSIS" engine with a different layout (fft_cmsis) — are dlopen'ed
// into this one process with RTLD_GLOBAL, the loader-visibility mode in
// which an ELF image's weak, default-visibility definitions are candidates
// for every later image's references (on macOS dyld coalesces weak
// definitions across images in either mode). Then each image is asked what
// its embedders see.
//
// What is pinned, and where:
//   - TaggedEmbedderIsImmuneToCrossImageCoalescing: an embedder of
//     basic_real_fft<float> defined inside `inline namespace
//     TAP_DSP_FFT_ABI` (the shape of pvoc.h and log_mel.h) runs its OWN
//     image's code in both images — which() reports the image's own tag and
//     the layout the bound member sees equals the image's real layout — and
//     the shipping embedders pvoc32 and log_mel32 compute finite, sane
//     results in both images with both loaded. This is the F4 property
//     itself, measured in a process rather than read off an nm listing.
//   - UntaggedEmbedderIsWhatTheTagExistsFor: the identical embedder defined
//     in plain tap::dsp (what a consumer class that has not adopted the tag
//     is). Its outcome is platform-defined and is printed as `[ measured ]`;
//     on linux with GCC (the CI host) it is asserted: image B's call binds
//     to image A's definition — which() returns A's tag from inside B and
//     the member sees A's layout against B's real one — so the hazard is
//     real and this test can see it, which is what makes the first test's
//     immunity meaningful. Under RTLD_LOCAL on ELF nothing coalesces (the
//     35a review measured it; not run here, because both modes cannot share
//     one process: the first dlopen decides an image's bindings).
//
// Hosted, non-Windows only (tests/CMakeLists.txt): dlopen and MODULE
// libraries; the QEMU legs have neither, and the design note states the
// Windows situation (an image exports nothing without __declspec, so the
// coalescing has no path there).

#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstring>
#include <string>

#include <dlfcn.h>
#include <gtest/gtest.h>

#include "tap/dsp/fft.h"

#if !defined(TAP_DSP_ABI_IMAGE_A) || !defined(TAP_DSP_ABI_IMAGE_B)
#error "tests/CMakeLists.txt passes the two image paths as TAP_DSP_ABI_IMAGE_A / _B"
#endif

namespace {

    using string_fn = const char* (*)();
    using size_fn   = std::size_t (*)();
    using double_fn = double (*)();

    struct image {
        void*       handle = nullptr;
        const char* name   = "";

        template <typename Fn>
        Fn get(const char* symbol) const {
            // A C-style cast is what dlsym's void* -> function pointer needs; the
            // POSIX-sanctioned form.
            Fn fn = reinterpret_cast<Fn>(dlsym(handle, symbol)); // NOLINT
            EXPECT_NE(fn, nullptr) << name << ": " << symbol << ": " << dlerror();
            return fn;
        }
        const char* str(const char* symbol) const { return get<string_fn>(symbol)(); }
        std::size_t size(const char* symbol) const { return get<size_fn>(symbol)(); }
        double      value(const char* symbol) const { return get<double_fn>(symbol)(); }
    };

    /// Loads A then B with RTLD_NOW | RTLD_GLOBAL, once per process; every
    /// test reads the same handles, so test order cannot change what the
    /// loader bound (the first dlopen of an image decides its bindings).
    struct images {
        image a;
        image b;
    };
    const images& loaded() {
        static const images k_both = [] {
            images r;
            r.a.name   = "image A (" TAP_DSP_ABI_IMAGE_A ")";
            r.b.name   = "image B (" TAP_DSP_ABI_IMAGE_B ")";
            r.a.handle = dlopen(TAP_DSP_ABI_IMAGE_A, RTLD_NOW | RTLD_GLOBAL);
            if (r.a.handle == nullptr) {
                std::fprintf(stderr, "dlopen A: %s\n", dlerror());
            }
            r.b.handle = dlopen(TAP_DSP_ABI_IMAGE_B, RTLD_NOW | RTLD_GLOBAL);
            if (r.b.handle == nullptr) {
                std::fprintf(stderr, "dlopen B: %s\n", dlerror());
            }
            return r;
        }();
        return k_both;
    }

    TEST(fft_abi_tag, TwoImagesWithDifferentDefaultsLoadIntoOneProcess) {
        const images& im = loaded();
        ASSERT_NE(im.a.handle, nullptr);
        ASSERT_NE(im.b.handle, nullptr);
        // A is this build; B is the stub CMSIS default.
        EXPECT_STREQ(im.a.str("tap_dsp_abi_image_built_as"), tap::dsp::k_real_fft_abi_tag);
        EXPECT_STREQ(im.b.str("tap_dsp_abi_image_built_as"), "fft_cmsis");
        // The precondition of the hazard: the two images' layouts differ.
        EXPECT_NE(im.a.size("tap_dsp_abi_image_tagged_layout_real"), im.b.size("tap_dsp_abi_image_tagged_layout_real"));
        EXPECT_NE(im.a.size("tap_dsp_abi_image_pvoc_layout"), im.b.size("tap_dsp_abi_image_pvoc_layout"));
        EXPECT_NE(im.a.size("tap_dsp_abi_image_log_mel_layout"), im.b.size("tap_dsp_abi_image_log_mel_layout"));
    }

    TEST(fft_abi_tag, TaggedEmbedderIsImmuneToCrossImageCoalescing) {
        const images& im = loaded();
        ASSERT_NE(im.a.handle, nullptr);
        ASSERT_NE(im.b.handle, nullptr);
        for (const image* i : {&im.a, &im.b}) {
            const char* built = i->str("tap_dsp_abi_image_built_as");
            EXPECT_STREQ(i->str("tap_dsp_abi_image_tagged_which"), built) << i->name;
            EXPECT_EQ(i->size("tap_dsp_abi_image_tagged_layout_seen"), i->size("tap_dsp_abi_image_tagged_layout_real"))
                << i->name;
        }
        // The shipping embedders, run in both images with both loaded. On
        // linux both images are the split-radix arithmetic underneath (the
        // stub wraps it), so the checksums are bit-identical; on the macOS
        // leg image A is vDSP and the two differ by the engines' float
        // rounding, so the comparison there is an order-of-magnitude
        // crash-and-garbage detector, not a floor (the floors are pinned in
        // test_fft_backend.cpp): 1e-3 relative is three-plus decades above
        // the engines' measured 1e-7-class disagreement and far below what
        // running one image's code over the other's layout produces, if it
        // returns at all.
        const double pv_a = im.a.value("tap_dsp_abi_image_pvoc_checksum");
        const double pv_b = im.b.value("tap_dsp_abi_image_pvoc_checksum");
        const double lm_a = im.a.value("tap_dsp_abi_image_log_mel_checksum");
        const double lm_b = im.b.value("tap_dsp_abi_image_log_mel_checksum");
        std::printf("[ measured ] pvoc checksum A %.9g B %.9g; log_mel checksum A %.9g B %.9g\n", pv_a, pv_b, lm_a,
                    lm_b);
        EXPECT_TRUE(std::isfinite(pv_a) && std::isfinite(pv_b) && std::isfinite(lm_a) && std::isfinite(lm_b));
        EXPECT_NE(pv_a, 0.0);
        EXPECT_NE(lm_a, 0.0);
#if defined(TAP_DSP_FFT_ACCELERATE)
        EXPECT_NEAR(pv_a, pv_b, 1e-3 * std::fabs(pv_a));
        EXPECT_NEAR(lm_a, lm_b, 1e-3 * std::fabs(lm_a));
#else
        EXPECT_EQ(pv_a, pv_b);
        EXPECT_EQ(lm_a, lm_b);
#endif
    }

    TEST(fft_abi_tag, UntaggedEmbedderIsWhatTheTagExistsFor) {
        const images& im = loaded();
        ASSERT_NE(im.a.handle, nullptr);
        ASSERT_NE(im.b.handle, nullptr);
        const char*       a_built = im.a.str("tap_dsp_abi_image_built_as");
        const char*       b_which = im.b.str("tap_dsp_abi_image_untagged_which");
        const std::size_t b_seen  = im.b.size("tap_dsp_abi_image_untagged_layout_seen");
        const std::size_t b_real  = im.b.size("tap_dsp_abi_image_untagged_layout_real");
        std::printf(
            "[ measured ] untagged embedder in image B (fft_cmsis): which() = %s, layout seen/real = %zu/%zu -> %s\n",
            b_which, b_seen, b_real,
            std::strcmp(b_which, a_built) == 0 ? "bound to image A's definition (the F4 hazard)"
                                               : "bound to its own definition");
#if defined(__linux__) && defined(__GNUC__) && !defined(__clang__)
        // The CI host: ELF, GCC, RTLD_GLOBAL — B's references to the untagged
        // member functions bind to A's definitions, loaded first.
        EXPECT_STREQ(b_which, a_built) << "expected image B's untagged embedder to run image A's code";
        EXPECT_NE(b_seen, b_real) << "expected the bound member to see image A's layout";
#endif
        // Image A, loaded first, always sees itself.
        EXPECT_STREQ(im.a.str("tap_dsp_abi_image_untagged_which"), a_built);
    }

} // namespace
