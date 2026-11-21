/// @file iir_fixed_pins.h
/// @brief The pinned output bits of the fixed-point recursive kernels (tests/test_iir_fixed.cpp).
// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// FNV-1a-64 of every tabulated structure's outputs on a fixed input (1024
// xorshift32 samples at 0.9 of full scale, seed 0x5eed), concatenated in the
// order test_iir_fixed.cpp runs them. Integer arithmetic only, so one row
// holds for every compiler and target CI runs; a target that needs a second
// row is a finding.

#pragma once

#include <cstdint>

namespace tap::dsp::test::iir_fixed_pins {

    struct row {
        const char*   name;
        std::uint64_t fingerprint;
    };

    inline constexpr row k_q15[] = {
        {"nonlinear decimators", 0x2e0b8430bcb3d6beull},  {"nonlinear interpolators", 0xb3df8da9d1f93fbaull},
        {"near-linear both ways", 0x34aab8b916cc294eull}, {"four-path both ways", 0x749b0b692717f22cull},
        {"guard band cascades", 0x9e9aca1ce58ca203ull},
    };

    inline constexpr row k_q31[] = {
        {"nonlinear decimators", 0x3b54bc8e678b7093ull},  {"nonlinear interpolators", 0xecc443f6eb3c0433ull},
        {"near-linear both ways", 0x0ad784efa868a662ull}, {"four-path both ways", 0x86e225f30533a6cfull},
        {"guard band cascades", 0x95f675dc93501a9bull},
    };

} // namespace tap::dsp::test::iir_fixed_pins
