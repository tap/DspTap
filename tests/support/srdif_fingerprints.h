/// @file srdif_fingerprints.h
/// @brief The pinned output fingerprints of the srdif floating FFT engine (tests/test_fft_srdif_fingerprint.cpp).
// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// One row per distinct computation of the engine's output bits. The engine is
// libm-free and the target is built without fp-contraction, so one row is
// expected everywhere (the test's file comment); each row says where it was
// measured.

#pragma once

#include <cstddef>
#include <cstdint>

namespace tap::dsp::test::srdif_fingerprints {

    /// xorshift32 seed of the fixed input (test_fft_srdif_fingerprint.cpp).
    inline constexpr std::uint32_t k_seed = 0x2545F491u;

    /// Sizes 4, 8, … 65536: fifteen fingerprints per profile.
    inline constexpr std::size_t k_sizes = 15;

    struct row {
        const char*   name;
        std::uint64_t float_rows[k_sizes];
        std::uint64_t double_rows[k_sizes];
    };

    // Sizes 4 … 65536 in order. Measured 2026-09-26 on x86-64 Linux (g++
    // 13.3.0 and clang++ 18.1.3, -O3, glibc 2.39 under either libm dispatch)
    // and on the four QEMU legs (arm-none-eabi-gcc 13.2.1, MinSizeRel, newlib:
    // Cortex-M4 soft-float, M4F, M33, M55; they run N <= 4096), every one
    // printing exactly these values; CI's Windows x64 MSVC and macOS arm64
    // AppleClang jobs matched it too (review A of tap/DspTap#42).
    inline constexpr row k_rows[] = {
        {"every compiler and target CI runs (integer tables, -ffp-contract=off)",
         {
             0x5309153b08c0bba3ull,
             0xc527d0fa579ad0ccull,
             0x4a441c08fe488adeull,
             0xd8eca950de16500aull,
             0x9e243c036f0febf4ull,
             0x807bcbd6c7f16e6aull,
             0x0ff0149369e98358ull,
             0x9b34f68a603ea673ull,
             0xee6d855e8bfeccfbull,
             0x7e422b5325560912ull,
             0x25e2afe6da11eb13ull,
             0xf67eb863856172aaull,
             0x8db171597b3f35a3ull,
             0x3f4149c901a6f2feull,
             0x1923fc9f948e764bull,
         },
         {
             0x64c64f6ba81a39c5ull,
             0xea87a5a1806657b9ull,
             0xa815dd94189960b1ull,
             0x9bffba81bb5cb503ull,
             0x6495ef1662b72082ull,
             0xd25bd3850f62c066ull,
             0xa464d0c9d366ba61ull,
             0x79db05f364ec79d6ull,
             0x9b324fb8c8550239ull,
             0x2f1ce7e9491f4841ull,
             0x17077e5dc40796c4ull,
             0xbfef417cfeff5056ull,
             0x82025cf16b0c1d7bull,
             0x165cbd8f43b208fbull,
             0x829f64451960c162ull,
         }},
    };

} // namespace tap::dsp::test::srdif_fingerprints
