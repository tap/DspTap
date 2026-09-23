/// @file split_radix_fingerprints.h
/// @brief The pinned output fingerprints of the split-radix FFT engine (Decision D6/D10).
// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// One procedure and one table of pins, shared by the fingerprint test
// (tests/test_fft_split_radix_fingerprint.cpp) and, while the reference C is
// still in the tree, by the cross-check in tests/test_fft_parity_ooura.cpp
// that computes the same fingerprints from the C.
//
// The procedure. For each N in k_fingerprint_sizes — every power of two from
// 4 to 65536, because Ooura's dispatch branches on N (cftf040 / cftb040 at 8,
// cftf161 and bitrv216 at 32, cftfx41 at 64 and 128, the odd-log4 halves of
// bitrv2 / bitrv2conj, cftleaf's 512-point leaf at odd powers >= 2048) and
// only every size together reaches every statement of the engine — and each
// precision: four materials in a fixed order — xorshift32 broadband noise
// (seed 0x9E3779B9, random_signal, libm-free), a unit impulse, DC (all
// ones), and a full-scale +1/-1 alternation — each taken through
// forward_inplace and then, on that forward output, inverse_inplace
// (unnormalized). Two FNV-1a-64 folds run across all four materials: one
// over every forward output, one over every inverse output, each sample
// folded as its IEEE bit pattern (uint64 for double, uint32 for float),
// least significant byte first. One differing bit anywhere changes the fold,
// and +0 / -0 and NaN payloads count as the bits they are. The on-bin tone
// of the old parity sweep is left out on purpose: tone() synthesizes
// through libm, so its input bits would themselves vary by host.
//
// THE INVARIANT: the float row. split_radix.h has no precision-specific
// branch (no if constexpr, no is_same): every kernel statement is the same
// template code for both precisions. The float instantiation rounds each
// double libm result to float, and on every configuration measured that
// absorbed the libm differences, so the float fingerprints are ONE value
// everywhere (glibc under both CPU dispatches, newlib with and without a
// double-precision FPU, the UCRT, Apple's libm). A float pin that moves is a
// change to the engine's arithmetic, whatever the platform. Double agrees
// everywhere up to N = 64 as well — not because those twiddles are exact
// (only N = 4 reads none; N = 16 already reads cos(pi/4), 0.5 cos(pi/8) and
// 0.5 sin(pi/8)) but because every libm measured rounds them alike.
//
// Why the double rows are per C library build. The engine's tables are
// built from the C library's cos / sin / atan exactly as fftsg.c builds them
// (Decision D10; split_radix.h rule 2), so a last-bit difference between two
// libms moves the double outputs from N = 128 up, and moved the C's outputs
// identically: in one binary the C computed the same double values as the
// port on every row below. A row identifies a libm BUILD, including its
// run-time dispatch: x86-64 glibc selects FMA/AVX2 or SSE2 implementations
// of sin / cos / atan by CPU, and the two differ from N = 8192 up, so glibc
// carries a pair of rows and a run passes when all its cells equal one row.
// A double-only mismatch with the float pins holding is a libm or dispatch
// change, not an engine change: check it against the C (docs/fft-design.md,
// "The bit-identity record after D6") and add a row.

#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <type_traits>
#include <vector>

#include "signals.h"

namespace tap::dsp::test {

    /// The sizes fingerprinted: every power of two from 4 (no real
    /// post-pass) to 65536 (the top of the old parity sweep). Sizes above
    /// TAP_DSP_TEST_MAX_FFT_N (4096 on the emulated legs) are not run.
    inline constexpr std::array<std::size_t, 15> k_fingerprint_sizes{4,    8,    16,   32,   64,    128,   256,  512,
                                                                     1024, 2048, 4096, 8192, 16384, 32768, 65536};

    /// FNV-1a-64 over the IEEE bit patterns of a block, LSB first.
    template <typename Sample>
    void fnv1a64_fold_bits(std::uint64_t& h, const std::vector<Sample>& block) {
        using bits_t = std::conditional_t<sizeof(Sample) == 4, std::uint32_t, std::uint64_t>;
        static_assert(sizeof(bits_t) == sizeof(Sample));
        for (const Sample v : block) {
            bits_t u = 0;
            std::memcpy(&u, &v, sizeof u);
            for (std::size_t b = 0; b < sizeof u; ++b) {
                h ^= static_cast<std::uint64_t>((u >> (8 * b)) & 0xffu);
                h *= 0x100000001b3ull;
            }
        }
    }

    /// The four materials, in fold order.
    template <typename Sample>
    std::vector<std::vector<Sample>> fingerprint_materials(std::size_t n) {
        std::vector<std::vector<Sample>> m;
        m.push_back(random_signal<Sample>(n, 0x9E3779B9u)); // broadband
        std::vector<Sample> impulse(n, Sample(0));
        impulse[0] = Sample(1);
        m.push_back(impulse);
        m.push_back(std::vector<Sample>(n, Sample(1))); // DC
        std::vector<Sample> alternating(n);
        for (std::size_t i = 0; i < n; ++i) {
            alternating[i] = (i % 2 == 0) ? Sample(1) : Sample(-1);
        }
        m.push_back(alternating);
        return m;
    }

    struct fingerprint_pair {
        std::uint64_t forward = 0;
        std::uint64_t inverse = 0;
    };

    /// The procedure above, for any engine with the constructor-from-size /
    /// forward_inplace / inverse_inplace surface (split_radix_rdft; the C's
    /// rdft behind the same surface in the parity TU).
    template <typename Engine, typename Sample>
    fingerprint_pair fingerprint_of(std::size_t n) {
        Engine           engine(n);
        fingerprint_pair f{0xcbf29ce484222325ull, 0xcbf29ce484222325ull};
        for (std::vector<Sample> x : fingerprint_materials<Sample>(n)) {
            engine.forward_inplace(x.data());
            fnv1a64_fold_bits(f.forward, x);
            engine.inverse_inplace(x.data());
            fnv1a64_fold_bits(f.inverse, x);
        }
        return f;
    }

    /// One row of pins, indexed like k_fingerprint_sizes. A zero pair is a
    /// size never measured there (the emulated legs stop at 4096).
    using fingerprint_pins = std::array<fingerprint_pair, k_fingerprint_sizes.size()>;

    /// A double row: the C library build it was measured on, and its pins.
    struct double_row {
        const char*      name;
        fingerprint_pins pins;
    };

    // clang-format off
    /// float: one row, every configuration measured (the invariant above).
    inline constexpr fingerprint_pins k_float_pins{{
        {0xef9b60ebf0e8f587ull, 0x7d9b9e4480f66924ull}, //     4
        {0x96742dab1b577d2cull, 0x02f73648a60330caull}, //     8
        {0xa222e5c5b80bba3bull, 0x3a182210d01399a7ull}, //    16
        {0x06e774c9f71d5af5ull, 0x68592652a0ecee43ull}, //    32
        {0x31eb423adbb65986ull, 0x9d67d026c4f48eccull}, //    64
        {0xc2f53e32cca6945full, 0xd2da49eb3feccae5ull}, //   128
        {0x1dc8ccaa8f8db08dull, 0x50306aa83541329dull}, //   256
        {0x3272fd6659a91c9aull, 0x77384618c779efa6ull}, //   512
        {0xfd0ffea8357b8dddull, 0x526e4272bdaab623ull}, //  1024
        {0x2e8a6f5fb458520dull, 0x23fc4df37b433aecull}, //  2048
        {0x6723317886e7c1f2ull, 0x7a318e19506b1c24ull}, //  4096
        {0xa5edafd7ec823192ull, 0x530ea64869bf67aeull}, //  8192
        {0xe562a21922c4c589ull, 0xed246b24b4c2f3dfull}, // 16384
        {0x5ce3d674258b6035ull, 0x3ad5668f242936a3ull}, // 32768
        {0x141744d674ac1c64ull, 0x8c9d260bedee7129ull}, // 65536
    }};

    /// The rows every platform shares up to N = 64 (see the file comment).
#define TAP_DSP_FINGERPRINT_DOUBLE_TO_64                                                                               \
        {0xf014008109d04382ull, 0xa4e32575aed1c801ull}, /*     4 */                                                    \
        {0x85c71bbb0b4690b6ull, 0x477a3c910df5bca9ull}, /*     8 */                                                    \
        {0xea32e166ad289557ull, 0x4a89b6e54a2b629eull}, /*    16 */                                                    \
        {0xb8aec28ef6d673ebull, 0xe8caa0863bb9b639ull}, /*    32 */                                                    \
        {0x5644ff7c7b0e6419ull, 0x4d1bee09c0ebc556ull}  /*    64 */

    /// double, per C library build. Selected by the macros that identify the
    /// libm the engine's tables are built from; a platform none of them names
    /// has no rows, and the double test fails there after printing its values.
#if defined(__GLIBC__) && defined(__x86_64__)
    inline constexpr const char* k_double_platform = "x86-64 glibc";
    inline constexpr std::array<double_row, 2> k_double_rows{{
        {"x86-64 glibc (2.39), FMA/AVX2 dispatch", {{
            TAP_DSP_FINGERPRINT_DOUBLE_TO_64,
            {0x6a7199cab6054194ull, 0xe27bedfa6fec4919ull}, //   128
            {0x4e70f1158efcfe3dull, 0xb5d51404f905fea6ull}, //   256
            {0x21e35dd0330a809aull, 0xd5c2231d159d4698ull}, //   512
            {0xf24f19a7fe3274a8ull, 0x54dce063fb4161ebull}, //  1024
            {0xf3a7d6278d939926ull, 0x1c95174b8f085c8eull}, //  2048
            {0x753d7cfca461e82bull, 0x8c9e917a44067ff2ull}, //  4096
            {0x6165e9ec3f1269a7ull, 0x32d82388dc2a21fcull}, //  8192
            {0xb54f372489309a93ull, 0x3efc33f7c7e1e0adull}, // 16384
            {0x59c3d563577f7644ull, 0x7d40ace7f287957eull}, // 32768
            {0xb7909fad30a275afull, 0xe72fb50d719ef181ull}, // 65536
        }}},
        {"x86-64 glibc (2.39), SSE2 dispatch", {{
            TAP_DSP_FINGERPRINT_DOUBLE_TO_64,
            {0x6a7199cab6054194ull, 0xe27bedfa6fec4919ull}, //   128
            {0x4e70f1158efcfe3dull, 0xb5d51404f905fea6ull}, //   256
            {0x21e35dd0330a809aull, 0xd5c2231d159d4698ull}, //   512
            {0xf24f19a7fe3274a8ull, 0x54dce063fb4161ebull}, //  1024
            {0xf3a7d6278d939926ull, 0x1c95174b8f085c8eull}, //  2048
            {0x753d7cfca461e82bull, 0x8c9e917a44067ff2ull}, //  4096
            {0xcc1adb16b8385bf9ull, 0x655eced7f8b256f4ull}, //  8192
            {0x22e9355d131725acull, 0xe8b892b022107e49ull}, // 16384
            {0x0b12f1160600a08bull, 0x33cb520f47ba6d5bull}, // 32768
            {0x23f50b8c94aecb5bull, 0xaa919e1c04c85b21ull}, // 65536
        }}},
    }};
#elif defined(__NEWLIB__) && defined(__arm__) && defined(__ARM_FP) && (__ARM_FP & 0x8)
    inline constexpr const char* k_double_platform = "Arm newlib, double-precision FPU";
    inline constexpr std::array<double_row, 1> k_double_rows{{
        {"Arm newlib, double-precision FPU (M55)", {{
            TAP_DSP_FINGERPRINT_DOUBLE_TO_64,
            {0xf94fa452b6b6631bull, 0x73c4adb6801f242aull}, //   128
            {0x5c6c45c5898a558aull, 0x43318e1983699fd3ull}, //   256
            {0x8b67b1235ed466a9ull, 0xb5065c9484181d81ull}, //   512
            {0xf5992e25957ac958ull, 0x6c2e7e635f21adf7ull}, //  1024
            {0xc1c49092ec9a4212ull, 0x6ee7a9907f767089ull}, //  2048
            {0xc42e9fb613c1c2c9ull, 0x81d554f525c82f66ull}, //  4096
            {0, 0}, {0, 0}, {0, 0}, {0, 0},                 // 8192 ... 65536: not run
        }}},
    }};
#elif defined(__NEWLIB__) && defined(__arm__)
    inline constexpr const char* k_double_platform = "Arm newlib, software double";
    inline constexpr std::array<double_row, 1> k_double_rows{{
        {"Arm newlib, software double (M4, M4F, M33)", {{
            TAP_DSP_FINGERPRINT_DOUBLE_TO_64,
            {0xf94fa452b6b6631bull, 0x73c4adb6801f242aull}, //   128
            {0xe59fe4d3e6326eb3ull, 0x43318e1983699fd3ull}, //   256
            {0x5c2479007e97b9c6ull, 0x00306b85c97aee64ull}, //   512
            {0x878a01ec9aa4e1bdull, 0x068c7471dbde8796ull}, //  1024
            {0x4f5075e73328a14aull, 0x731f6158e19aadd7ull}, //  2048
            {0x58c61a3b29d4b55eull, 0x2d0482b12bc4e246ull}, //  4096
            {0, 0}, {0, 0}, {0, 0}, {0, 0},                 // 8192 ... 65536: not run
        }}},
    }};
#elif defined(_MSC_VER) && defined(_M_X64)
    inline constexpr const char* k_double_platform = "x64 MSVC, UCRT";
    inline constexpr std::array<double_row, 1> k_double_rows{{
        {"x64 MSVC, UCRT", {{
            TAP_DSP_FINGERPRINT_DOUBLE_TO_64,
            {0xf94fa452b6b6631bull, 0x73c4adb6801f242aull}, //   128
            {0xcabb356393c2dac7ull, 0x5c0eeb0579868a88ull}, //   256
            {0x5563d637de5e961dull, 0x10e1a0bd60bce082ull}, //   512
            {0xdd44bf873268ba83ull, 0xabc8af10c100bee6ull}, //  1024
            {0xe8a241c435bb2601ull, 0xca0d22bb81356f5aull}, //  2048
            {0x4e0db6e49731d286ull, 0xbcbcc885236c3585ull}, //  4096
            {0x25213b8387ce17dfull, 0x184a6995c63e5507ull}, //  8192
            {0x2e24dc611386f95aull, 0xb6f8aab7dc17962dull}, // 16384
            {0x128e42969938f8d4ull, 0x1fdf8d272bf2b6a1ull}, // 32768
            {0xb02d408a57e2ebdeull, 0x288c58dad9580f4bull}, // 65536
        }}},
    }};
#elif defined(__APPLE__) && defined(__aarch64__)
    inline constexpr const char* k_double_platform = "arm64 macOS";
    inline constexpr std::array<double_row, 1> k_double_rows{{
        {"arm64 macOS", {{
            TAP_DSP_FINGERPRINT_DOUBLE_TO_64,
            {0xf94fa452b6b6631bull, 0x73c4adb6801f242aull}, //   128
            {0x1a682072ef54d745ull, 0xa9e60433894b7013ull}, //   256
            {0xd0a19e500f844c5cull, 0x76d38e82464248afull}, //   512
            {0x5cb137231c0a46e9ull, 0xff6e0741acfebeb3ull}, //  1024
            {0x379d6b9eef0fbe66ull, 0xf80a8d3a5cf07d68ull}, //  2048
            {0x7909ca3d02d78252ull, 0x4a9a5c439158323dull}, //  4096
            {0xdab35c43e9136ed6ull, 0x503f3b059e7ac263ull}, //  8192
            {0xbb8b79e8470c05cfull, 0x0496245fe2a1a0e4ull}, // 16384
            {0x588bd3585488aef2ull, 0x032c156afe537e02ull}, // 32768
            {0xbae4e82f8d79d549ull, 0xbfa3a82a87bbe16dull}, // 65536
        }}},
    }};
#else
    inline constexpr const char*                k_double_platform = "no pinned C library";
    inline constexpr std::array<double_row, 0> k_double_rows{};
#endif
#undef TAP_DSP_FINGERPRINT_DOUBLE_TO_64
    // clang-format on

    /// The number of leading sizes a build runs (those <= max_n).
    constexpr std::size_t fingerprint_size_count(std::size_t max_n) {
        std::size_t count = 0;
        while (count < k_fingerprint_sizes.size() && k_fingerprint_sizes[count] <= max_n) {
            ++count;
        }
        return count;
    }

    /// Whether @p row pins cell @p i and equals @p got there.
    constexpr bool row_matches_cell(const double_row& row, std::size_t i, const fingerprint_pair& got) {
        const fingerprint_pair& pin = row.pins[i];
        return pin.forward != 0 && pin.forward == got.forward && pin.inverse == got.inverse;
    }

    /// The first double row whose pins equal every measured cell of @p got
    /// (got[i] for size index i), or nullptr when no row matches them all.
    inline const double_row* matching_double_row(const std::vector<fingerprint_pair>& got) {
        for (const double_row& row : k_double_rows) {
            bool all = true;
            for (std::size_t i = 0; i < got.size() && all; ++i) {
                all = row_matches_cell(row, i, got[i]);
            }
            if (all) {
                return &row;
            }
        }
        return nullptr;
    }

} // namespace tap::dsp::test
