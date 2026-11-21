// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// The constexpr designs' bits, pinned: IIR-0's bit test
// (SampleRateTap's rational/tools/spikes/iir/constexpr/, cx.BitsMatchGccAndClang)
// moved into DspTap's battery so every compiler CI runs it — GCC and Clang on
// Linux, MSVC (/W4 /WX) on Windows, AppleClang on macOS, arm-none-eabi-gcc on
// the four QEMU legs — evaluates the same tables at compile time and must
// print the same bits (IIR_PLAN I4).
//
// What is pinned (tests/support/iir_constexpr_pins.h): an FNV-1a-64 over the
// IEEE bit patterns of every tabulated design, in its own row per design so
// a mismatch names it —
//   - the nonlinear half-band (k_nonlinear_halfband) at the four spec points
//     and the divisors 1, 2, 4, 8: its coefficients;
//   - the Chebyshev II (k_guard_band_cheby2, the committed L-inf gains
//     applied) at the four spec points and orders 3 to 11: its sections and
//     stop edge;
//   - the Q15 and Q31 coefficient tables rounded from them at compile time
//     (nonlinear_halfband_table, biquad_cascade_table), whose integers are
//     what the fixed-point outputs are a function of.
// The rule is srdif's (tests/test_fft_srdif_fingerprint.cpp): one row for
// every compiler and target CI runs; a compiler that needs a second row is a
// finding, not a row to add. The arithmetic is IEEE double in the constant
// evaluator (no libm, no contraction, no excess precision), so one row is
// expected. Measured 2026-11-21: g++ 13.3 and clang++ 18.1 on x86-64 Linux,
// arm-none-eabi-gcc 13.2.1 for the M4, M4F, M33 and M55, and the Hexagon
// clang 19.1.5, every one printing these values; MSVC and AppleClang are
// CI's to confirm.
//
// The same functions at run time (a custom design at construction) give the
// same bits here, because this target is built with -ffp-contract=off
// (tests/CMakeLists.txt); with contraction on, a target with an FMA may
// differ in the last bits at run time, never at compile time.
//
// Constant-evaluation cost (I4: keep every evaluation well under MSVC's
// documented default of 100,000 steps), measured by bisecting clang 18's
// -fconstexpr-steps and g++ 13's -fconstexpr-ops-limit over one table each:
// the order-17 half-band (120 dB, 5/12, d = 1, the largest tabulated) 5,835
// clang steps (30,679 g++ operations), the order-11 Chebyshev II (120 dB)
// 2,723 (21,427): IIR-0's prototype took 30,487 and, with the L-inf sweep,
// 414,830.

#include <bit>
#include <cstdint>
#include <cstdio>
#include <span>
#include <string_view>

#include <gtest/gtest.h>

#include "support/iir_constexpr_pins.h"
#include "tap/dsp/iir_design.h"
#include "tap/dsp/iir_kernels.h"

namespace {

    using namespace tap::dsp;
    namespace pins = tap::dsp::test::iir_constexpr_pins;

    struct fnv {
        std::uint64_t h = 0xcbf29ce484222325ull;
        void          add(std::uint64_t v) {
            for (int b = 0; b < 8; ++b) {
                h ^= (v >> (8 * b)) & 0xffu;
                h *= 0x100000001b3ull;
            }
        }
        void add(double v) { add(std::bit_cast<std::uint64_t>(v)); }
        void add(std::int32_t v) { add(static_cast<std::uint64_t>(static_cast<std::uint32_t>(v))); }
    };

    std::uint64_t hash(const nonlinear_halfband_design& d) {
        fnv f;
        for (std::size_t i = 0; i < d.count; ++i) {
            f.add(d.coefficients[i]);
        }
        return f.h;
    }
    std::uint64_t hash(const cheby2_design& d) {
        fnv f;
        for (std::size_t i = 0; i < d.count; ++i) {
            const iir_section& s = d.sections[i];
            for (const double v : {s.b0, s.b1, s.b2, s.a1, s.a2}) {
                f.add(v);
            }
        }
        f.add(d.stop_edge);
        return f.h;
    }
    template <typename S>
    std::uint64_t hash(const nonlinear_halfband_table<S>& t) {
        fnv f;
        for (std::size_t i = 0; i < t.size0; ++i) {
            f.add(t.path0[i]);
        }
        for (std::size_t i = 0; i < t.size1; ++i) {
            f.add(t.path1[i]);
        }
        return f.h;
    }
    template <typename S>
    std::uint64_t hash(const biquad_cascade_table<S>& t) {
        fnv f;
        for (std::size_t i = 0; i < t.size; ++i) {
            for (const auto c : t.sections[i]) {
                f.add(c);
            }
        }
        return f.h;
    }

    std::uint64_t pinned(std::span<const pins::row> rows, std::string_view name) {
        for (const pins::row& r : rows) {
            if (name == r.name) {
                return r.fingerprint;
            }
        }
        return 0;
    }

    void check(std::span<const pins::row> rows, const char* name, std::uint64_t got) {
        // Two 32-bit halves: newlib's printf on the QEMU legs parses neither
        // %zu nor PRIx64 (srdif's fingerprint test does the same).
        std::printf("[ measured ] %-26s 0x%08lx%08lx\n", name, static_cast<unsigned long>(got >> 32),
                    static_cast<unsigned long>(got & 0xffffffffu));
        EXPECT_EQ(got, pinned(rows, name)) << name << ": a constexpr table's bits moved";
    }

    template <int A, int Num, int Den, int D>
    void halfband_row(const char* name) {
        static constexpr nonlinear_halfband_design k_d = k_nonlinear_halfband<A, Num, Den, D>;
        check(pins::k_halfband, name, hash(k_d));
        // The fixed-point tables round from the same compile-time doubles.
        static constexpr nonlinear_halfband_table<std::int16_t> k_q15(k_d);
        static constexpr nonlinear_halfband_table<std::int32_t> k_q31(k_d);
        static_assert(k_q15.multiplies() == k_d.count);
        EXPECT_EQ(hash(k_q15), hash(k_q31)) << "Q15 and Q31 share the Q0.31 one-multiplier coefficients";
        check(pins::k_halfband_fixed, name, hash(k_q15));
        // At run time, the same function gives the same bits (no contraction here).
        const auto runtime = design_nonlinear_halfband(static_cast<double>(A), static_cast<double>(Num) / Den, D);
        EXPECT_EQ(hash(runtime), hash(k_d)) << name;
    }

    template <int A, int Num, int Den, int Order>
    void cheby2_row(const char* name) {
        static constexpr cheby2_design k_d = k_guard_band_cheby2<A, Num, Den, Order>;
        check(pins::k_cheby2, name, hash(k_d));
        static constexpr biquad_cascade_table<std::int16_t> k_q(k_d);
        check(pins::k_cheby2_fixed, name, hash(k_q));
        const auto runtime = design_cheby2_lowpass(Order, static_cast<double>(Num) / Den, static_cast<double>(A));
        EXPECT_EQ(hash(runtime), hash(k_d)) << name << " (its committed gains are all 1.0)";
    }

    TEST(IirConstexpr, HalfbandBitsMatchThePins) {
        halfband_row<70, 1, 3, 1>("70 dB 1/3 d1");
        halfband_row<70, 1, 3, 2>("70 dB 1/3 d2");
        halfband_row<70, 1, 3, 4>("70 dB 1/3 d4");
        halfband_row<70, 1, 3, 8>("70 dB 1/3 d8");
        halfband_row<70, 3, 8, 1>("70 dB 3/8 d1");
        halfband_row<70, 3, 8, 2>("70 dB 3/8 d2");
        halfband_row<70, 3, 8, 4>("70 dB 3/8 d4");
        halfband_row<70, 3, 8, 8>("70 dB 3/8 d8");
        halfband_row<70, 19, 48, 1>("70 dB 19/48 d1");
        halfband_row<70, 19, 48, 2>("70 dB 19/48 d2");
        halfband_row<70, 19, 48, 4>("70 dB 19/48 d4");
        halfband_row<70, 19, 48, 8>("70 dB 19/48 d8");
        halfband_row<120, 5, 12, 1>("120 dB 5/12 d1");
        halfband_row<120, 5, 12, 2>("120 dB 5/12 d2");
        halfband_row<120, 5, 12, 4>("120 dB 5/12 d4");
        halfband_row<120, 5, 12, 8>("120 dB 5/12 d8");
    }

    TEST(IirConstexpr, Cheby2BitsMatchThePins) {
        cheby2_row<70, 1, 3, 3>("70 dB 1/3 order 3");
        cheby2_row<70, 1, 3, 5>("70 dB 1/3 order 5");
        cheby2_row<70, 1, 3, 7>("70 dB 1/3 order 7");
        cheby2_row<70, 1, 3, 9>("70 dB 1/3 order 9");
        cheby2_row<70, 1, 3, 11>("70 dB 1/3 order 11");
        cheby2_row<70, 3, 8, 3>("70 dB 3/8 order 3");
        cheby2_row<70, 3, 8, 5>("70 dB 3/8 order 5");
        cheby2_row<70, 3, 8, 7>("70 dB 3/8 order 7");
        cheby2_row<70, 3, 8, 9>("70 dB 3/8 order 9");
        cheby2_row<70, 3, 8, 11>("70 dB 3/8 order 11");
        cheby2_row<70, 19, 48, 3>("70 dB 19/48 order 3");
        cheby2_row<70, 19, 48, 5>("70 dB 19/48 order 5");
        cheby2_row<70, 19, 48, 7>("70 dB 19/48 order 7");
        cheby2_row<70, 19, 48, 9>("70 dB 19/48 order 9");
        cheby2_row<70, 19, 48, 11>("70 dB 19/48 order 11");
        cheby2_row<120, 5, 12, 3>("120 dB 5/12 order 3");
        cheby2_row<120, 5, 12, 5>("120 dB 5/12 order 5");
        cheby2_row<120, 5, 12, 7>("120 dB 5/12 order 7");
        cheby2_row<120, 5, 12, 9>("120 dB 5/12 order 9");
        cheby2_row<120, 5, 12, 11>("120 dB 5/12 order 11");
    }

    TEST(IirConstexpr, TheTablesAreCompileTimeConstants) {
        // Usable where only a constant expression is: array bounds and
        // static_asserts.
        static_assert(k_nonlinear_halfband<70, 3, 8, 1>.coefficients[0] > 0.0);
        static_assert(k_nonlinear_halfband<70, 3, 8, 1>.coefficients[3] < 1.0);
        static_assert(k_guard_band_cheby2<120, 5, 12, 11>.count == 6);
        static_assert(k_near_linear_halfband<70, 3, 8, 2>.delay == 3);
        constexpr nonlinear_halfband_table<std::int16_t> k_q(k_nonlinear_halfband<70, 3, 8, 1>);
        int                                              sized[k_q.state_size()]{};
        EXPECT_EQ(sizeof(sized) / sizeof(sized[0]), 6u);
    }

} // namespace
