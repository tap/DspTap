/// @file iir_constexpr_pins.h
/// @brief The pinned bits of the constexpr recursive designs (tests/test_iir_constexpr.cpp).
// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// One fingerprint per design: FNV-1a-64 over its values' bit patterns, in
// the order test_iir_constexpr.cpp hashes them. One row for every compiler
// and target CI runs (that test's file comment says where they were
// measured); a compiler that needs a second row is a finding.

#pragma once

#include <cstdint>

namespace tap::dsp::test::iir_constexpr_pins {

    struct row {
        const char*   name;
        std::uint64_t fingerprint;
    };

    /// The nonlinear half-bands' double coefficients.
    inline constexpr row k_halfband[] = {
        {"70 dB 1/3 d1", 0xdf6e7799f934fb14ull},   {"70 dB 1/3 d2", 0x24fb0306bbe599c8ull},
        {"70 dB 1/3 d4", 0x4dd97b29f692a952ull},   {"70 dB 1/3 d8", 0x196d51279d529ddbull},
        {"70 dB 3/8 d1", 0xa7cefbc947484d09ull},   {"70 dB 3/8 d2", 0x9af13003a0a455b0ull},
        {"70 dB 3/8 d4", 0x2c0485b6e1d3a814ull},   {"70 dB 3/8 d8", 0x5b3bcdbc3e45c53aull},
        {"70 dB 19/48 d1", 0xe6d3c45170533a22ull}, {"70 dB 19/48 d2", 0x5a88c1e7ed669101ull},
        {"70 dB 19/48 d4", 0xcf9f55bbc7b488dbull}, {"70 dB 19/48 d8", 0xeb14d3121acc4c4bull},
        {"120 dB 5/12 d1", 0xa35ef103f21ccc43ull}, {"120 dB 5/12 d2", 0xcfb8af6358f92cfbull},
        {"120 dB 5/12 d4", 0x619118ebcba7d08bull}, {"120 dB 5/12 d8", 0x0b43cfcc5743a176ull},
    };

    /// Their Q0.31 coefficients (Q15 and Q31 share them).
    inline constexpr row k_halfband_fixed[] = {
        {"70 dB 1/3 d1", 0x7f1de3e1ee54273dull},   {"70 dB 1/3 d2", 0xa14ed85c8aa49e47ull},
        {"70 dB 1/3 d4", 0x39e18e14cf591428ull},   {"70 dB 1/3 d8", 0x307977fe25d14241ull},
        {"70 dB 3/8 d1", 0xcce91f88020742a1ull},   {"70 dB 3/8 d2", 0x59d68fe6df6a9355ull},
        {"70 dB 3/8 d4", 0x59345652a5276d0eull},   {"70 dB 3/8 d8", 0xf15295fa74ac8748ull},
        {"70 dB 19/48 d1", 0x5728317a204c0096ull}, {"70 dB 19/48 d2", 0x54c434dcad416480ull},
        {"70 dB 19/48 d4", 0x118351945b8e86d6ull}, {"70 dB 19/48 d8", 0x2eee8f3d7a0a2be8ull},
        {"120 dB 5/12 d1", 0xc72d5fa3d2d3a4c9ull}, {"120 dB 5/12 d2", 0xb4484f64e295fc6dull},
        {"120 dB 5/12 d4", 0x48f6c3e2000da2fdull}, {"120 dB 5/12 d8", 0x5f450ccc2bdf08f5ull},
    };

    /// The Chebyshev II designs' sections and stop edge.
    inline constexpr row k_cheby2[] = {
        {"70 dB 1/3 order 3", 0xd5307bdb64cc35bcull},    {"70 dB 1/3 order 5", 0x4c3ebe62135d4918ull},
        {"70 dB 1/3 order 7", 0x773e7395bae775edull},    {"70 dB 1/3 order 9", 0x790893df06a35405ull},
        {"70 dB 1/3 order 11", 0xa1ef95f5b3578614ull},   {"70 dB 3/8 order 3", 0x3f9f985a01c393efull},
        {"70 dB 3/8 order 5", 0xbe7474df8cf27b30ull},    {"70 dB 3/8 order 7", 0x37476a1ea9497795ull},
        {"70 dB 3/8 order 9", 0x219fb5318fe0e1cfull},    {"70 dB 3/8 order 11", 0x8da362db5be83fd5ull},
        {"70 dB 19/48 order 3", 0x73e9d54439a73400ull},  {"70 dB 19/48 order 5", 0xcce4f2cc489a2e42ull},
        {"70 dB 19/48 order 7", 0x321d3bd5c6f8bf4full},  {"70 dB 19/48 order 9", 0xdcf723dd38995e9bull},
        {"70 dB 19/48 order 11", 0xd7d5827727e1a059ull}, {"120 dB 5/12 order 3", 0x23a866ec091d4962ull},
        {"120 dB 5/12 order 5", 0x421e5212e24f5d6eull},  {"120 dB 5/12 order 7", 0x4d6371895f4699f2ull},
        {"120 dB 5/12 order 9", 0x7005412619e6fb6cull},  {"120 dB 5/12 order 11", 0xee1782e341f2aa98ull},
    };

    /// Their Q2.30 biquad tables ({b0, b1, b2, -a1, -a2} per section).
    inline constexpr row k_cheby2_fixed[] = {
        {"70 dB 1/3 order 3", 0x00798d7eca1e997aull},    {"70 dB 1/3 order 5", 0x4ec732b153eb05f8ull},
        {"70 dB 1/3 order 7", 0xd3ab9868ac9637b3ull},    {"70 dB 1/3 order 9", 0xc64832603f6adac2ull},
        {"70 dB 1/3 order 11", 0x9eb163ae449d8f30ull},   {"70 dB 3/8 order 3", 0x2257bb6229334be7ull},
        {"70 dB 3/8 order 5", 0xd699299b48b58c04ull},    {"70 dB 3/8 order 7", 0x43e11b263ac2eaf1ull},
        {"70 dB 3/8 order 9", 0xf4ba4a9a8af383c3ull},    {"70 dB 3/8 order 11", 0xdde0124922e0fd0full},
        {"70 dB 19/48 order 3", 0x0adad0ab56a0d351ull},  {"70 dB 19/48 order 5", 0x4a5c33c3b3159949ull},
        {"70 dB 19/48 order 7", 0x029d516a2ecba519ull},  {"70 dB 19/48 order 9", 0x439405fa52b584d7ull},
        {"70 dB 19/48 order 11", 0x9b4d8dfe5b38de72ull}, {"120 dB 5/12 order 3", 0x3fef7cf9ad753cffull},
        {"120 dB 5/12 order 5", 0x5cad0b6c79d55ad8ull},  {"120 dB 5/12 order 7", 0xb720eb87d7732f4cull},
        {"120 dB 5/12 order 9", 0x5c9c7e9186f01c84ull},  {"120 dB 5/12 order 11", 0xf33db5f0260fd7dcull},
    };

} // namespace tap::dsp::test::iir_constexpr_pins
