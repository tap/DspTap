// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// THE SPLIT-RADIX FINGERPRINTS (Decisions D6 and D10 of
// docs/audit-fft-and-code-smells.md): the bit patterns of
// tap::dsp::detail::split_radix_rdft's double and float outputs, forward and
// inverse, at every power of two from 4 to 65536 (4 to 4096 on the emulated
// legs, TAP_DSP_TEST_MAX_FFT_N), pinned as FNV-1a-64 folds. The procedure,
// the materials and the pins live in tests/support/split_radix_fingerprints.h.
//
// What the pins are. Each value was measured from this engine AND from the
// vendored C it replaced (fftsg.c for double, fftsg_float.c for float, both
// at -ffp-contract=off) in one binary, and the two were equal: the parity
// gate's ReferenceCHasThePinnedFingerprints (tests/test_fft_parity_ooura.cpp)
// computes the same values from the C on every CI leg while the C is in the
// tree. The pins therefore carry the Stage 2a parity gate's promise — the
// engine is the C, bit for bit (D10) — past the gate itself.
//
// How to read a failure. The FLOAT row is the kernel invariant: one value on
// every configuration measured, because the engine has no precision-specific
// branch and float rounding absorbs libm differences; a float pin that moves
// is a change to the engine. The DOUBLE rows are per C library build,
// because the tables come from libm (D10): a double cell that matches no row
// while the float pin at the same N holds points at the platform's libm
// version or CPU dispatch, and the message says which case it is.
//
// Why this target owns -ffp-contract=off. Bit identity with the C held for
// the C and the port compiled to the same sequence of IEEE operations;
// g++ contracts a*b+c into an FMA across statements by default and clang
// within one, so on FMA hardware (Apple arm64, the M55, x86 at -march) a
// default-flags build computes other bits for reasons unrelated to the
// engine's statements. MSVC: nothing is passed; /fp:precise does not
// contract (tests/CMakeLists.txt).
//
// A platform the double rows do not name fails after printing its values:
// its libm is unmeasured, and a new platform is recorded — measured against
// the C — not skipped.

#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <vector>

#include <gtest/gtest.h>

#include "support/split_radix_fingerprints.h"
#include "tap/dsp/fft/split_radix.h"

#ifndef TAP_DSP_TEST_MAX_FFT_N
#define TAP_DSP_TEST_MAX_FFT_N 65536
#endif

namespace {

    namespace t = tap::dsp::test;

    template <typename Sample>
    using engine = tap::dsp::detail::split_radix_rdft<Sample>;

    unsigned long hi32(std::uint64_t v) {
        return static_cast<unsigned long>(v >> 32);
    }
    unsigned long lo32(std::uint64_t v) {
        return static_cast<unsigned long>(v & 0xffffffffu);
    }

    constexpr std::size_t k_count = t::fingerprint_size_count(static_cast<std::size_t>(TAP_DSP_TEST_MAX_FFT_N));

    // Printed before any assertion so a -V log carries every host's values
    // whether or not they match (two 32-bit halves: newlib's printf has no
    // %llx).
    void print_cell(const char* precision, std::size_t n, const t::fingerprint_pair& got, const char* note) {
        std::printf("[ fingerprint ] %s n=%lu forward=%08lx%08lx inverse=%08lx%08lx (%s)\n", precision,
                    static_cast<unsigned long>(n), hi32(got.forward), lo32(got.forward), hi32(got.inverse),
                    lo32(got.inverse), note);
    }

    bool same(const t::fingerprint_pair& a, const t::fingerprint_pair& b) {
        return a.forward == b.forward && a.inverse == b.inverse;
    }

    // The float row is the kernel invariant: one value on every configuration
    // measured, independent of libm. A mismatch here is an engine change.
    TEST(split_radix_fingerprint, FloatOutputIsTheVendoredCBitPattern) {
        for (std::size_t i = 0; i < k_count; ++i) {
            const std::size_t         n   = t::k_fingerprint_sizes[i];
            const t::fingerprint_pair got = t::fingerprint_of<engine<float>, float>(n);
            print_cell("float", n, got, "every platform");
            EXPECT_TRUE(same(got, t::k_float_pins[i]))
                << "float n=" << n << ": the float row is one value on every configuration measured (the engine "
                << "has no precision-specific branch and float rounding absorbs libm differences), so this is a "
                << "change to the kernel's arithmetic, not to the platform: the engine is no longer the vendored C";
        }
    }

    // The double rows are per C library build (the support header). A run
    // passes when every cell equals one row; which one is printed.
    TEST(split_radix_fingerprint, DoubleOutputIsTheVendoredCBitPattern) {
        std::vector<t::fingerprint_pair> got;
        for (std::size_t i = 0; i < k_count; ++i) {
            const std::size_t n = t::k_fingerprint_sizes[i];
            got.push_back(t::fingerprint_of<engine<double>, double>(n));
            print_cell("double", n, got.back(), t::k_double_platform);
        }
        if (const t::double_row* row = t::matching_double_row(got)) {
            std::printf("[ fingerprint ] double matched row: %s\n", row->name);
            return;
        }
        if (t::k_double_rows.empty()) {
            ADD_FAILURE() << "double: no pins for this platform's C library; record the values printed above, "
                          << "measured against the C (docs/fft-design.md, \"The bit-identity record after D6\"), "
                          << "before relying on this build";
            return;
        }
        bool any_cell_mismatch = false;
        for (std::size_t i = 0; i < k_count; ++i) {
            const std::size_t n        = t::k_fingerprint_sizes[i];
            bool              in_a_row = false;
            bool              pinned   = false;
            for (const t::double_row& row : t::k_double_rows) {
                pinned   = pinned || row.pins[i].forward != 0;
                in_a_row = in_a_row || t::row_matches_cell(row, i, got[i]);
            }
            if (in_a_row) {
                continue;
            }
            any_cell_mismatch = true;
            if (!pinned) {
                ADD_FAILURE() << "double n=" << n << ": unmeasured on " << t::k_double_platform
                              << "; record the value printed above after checking it against the C";
                continue;
            }
            const bool float_holds = same(t::fingerprint_of<engine<float>, float>(n), t::k_float_pins[i]);
            if (float_holds) {
                ADD_FAILURE() << "double n=" << n << " matches no " << t::k_double_platform
                              << " row, but the float pin at this N holds, so the kernel is unchanged: suspect "
                              << "this platform's libm version or CPU dispatch, not the engine. Check the value "
                              << "against the C (docs/fft-design.md) and add a row.";
            }
            else {
                ADD_FAILURE() << "double n=" << n << " matches no " << t::k_double_platform
                              << " row and the float pin at this N differs too: the kernel changed, and the "
                              << "engine is no longer the vendored C";
            }
        }
        if (!any_cell_mismatch) {
            ADD_FAILURE() << "double: every cell matches some " << t::k_double_platform
                          << " row, but no single row matches them all";
        }
    }

} // namespace
