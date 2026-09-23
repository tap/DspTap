// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// THE BIT-IDENTITY GATE FOR THE C++20 PORT OF OOURA'S rdft (Stage 2a of
// docs/audit-fft-and-code-smells.md, Part 3; the policy behind it is Part 4's
// "fp-contraction policy" section and Part 6 items N3/N4).
//
// Two sides, each re-pointable in ONE place below:
//
//   ooura_ref<Sample>      the raw rdft / rdft_f of Takuya Ooura's fftsg.c,
//                          from the REFERENCE copy of the C under
//                          tests/reference/ooura/ (fftsg.c and fftsg_float.c,
//                          declared by tests/reference/ooura_rdft.h) that this
//                          target's CMake block compiles (tests/CMakeLists.txt).
//                          Since Stage 2c that copy is the only C in the repo:
//                          it left the shipping tree (Decision D6) and nothing
//                          but this gate compiles it.
//   engine_under_test      the C++20 port, detail::split_radix_rdft
//                          (include/tap/dsp/fft/split_radix.h), since Stage 2a.
//                          Before the port existed this alias named
//                          basic_real_fft, i.e. the same C, and the suite was
//                          Ooura-vs-Ooura and trivially green; the port was
//                          worked from that red re-point to green statement
//                          by statement. Stage 2b routed basic_real_fft at
//                          the port and the alias STAYS on the engine: on the
//                          M55 and macOS legs basic_real_fft<float> is a
//                          backend (CMSIS, vDSP), not the port, so pointing
//                          the alias at the class would turn the float gate
//                          red there for a reason that has nothing to do
//                          with the port. The class-to-engine identity is
//                          pinned separately by tests/test_fft_routing.cpp
//                          (memcmp, double and, where no backend is
//                          selected, float); the two files together are the
//                          flip's proof.
//
// The comparison is memcmp over the raw output bytes: not EXPECT_EQ (which
// calls +0.0 and -0.0 equal and any NaN unequal to itself), not a tolerance.
// Forward and inverse are both taken in place and unnormalized, exactly as the
// C exposes them, for every power of two from 4 to 65536 plus one run at 2^20
// (the largest size any consumer uses; item N20), on five materials chosen to
// reach different arithmetic: broadband xorshift noise (every butterfly busy),
// an on-bin tone (one live bin, the rest numerically empty), an impulse (flat
// spectrum, no cancellation), DC (maximal cancellation in every non-zero bin)
// and a full-scale alternating +1/-1 (all energy at Nyquist).
//
// Why this target owns its compiler flags. Bit identity between a C function
// and a C++ transliteration of it holds only if both are compiled to the same
// sequence of IEEE operations. Measured for Part 4: gcc -std=c17 does not
// contract a*b+c into an FMA; gcc -std=gnu17 does; g++ contracts in both c++20
// and gnu++20; clang contracts in every mode, statement-scoped. CMake sets none
// of this, so on any ISA with FMA (Apple arm64, Cortex-M55, x86 with -march)
// the two sides fuse differently by default and the identity is false for
// reasons that have nothing to do with the port. The gate therefore compiles
// BOTH the reference C and this TU with -ffp-contract=off (MSVC: /fp:precise
// does not contract, so nothing is passed there). Same binary, same libm:
// cross-platform identity is not claimed and does not hold today either,
// because libm's cos/sin differ in the last bit between glibc, newlib, UCRT
// and Apple (Part 4, "The float twiddles").
//
// To see that the flag is load-bearing rather than take it on faith, count
// the fused instructions in the reference object on an FMA-capable target:
//
//     cc -O3 -march=haswell [-ffp-contract=off] -c tests/reference/ooura/fftsg.c -o probe.o
//     objdump -d probe.o | grep -ciE 'vfmadd|vfmsub|vfnmadd|vfnmsub'
//
// Measured: gcc 13.3 200 / clang 18.1 171 fused instructions by default, 0
// with the flag, for both compilers. On plain x86-64 (no -march) the count is
// 0 either way because the ISA has no FMA, which is why the hosted CI legs
// cannot see a missing flag; the objdump check is the one that can.
//
// The same source builds a SECOND, informational target at default flags
// (TAP_DSP_PARITY_INFORMATIONAL). It measures the max-ulp deviation per N and
// never fails. It is not a pass/fail gate with an assumed bound because a
// different fusion choice per butterfly stage accumulates over log2 N stages
// and "1 ulp" would be a guess, not a measurement (item N3); its job is to
// put the number on the record for each platform so the policy fft.h states
// about exporting -ffp-contract=off can be decided on data. Where the number
// lands: ctest hides the stdout of a passing test under --output-on-failure,
// so the table is visible only when the label is run with -V (CI does that
// as its own step, Part 13) and, as RecordProperty values, in the JUnit XML
// that the CMake block requests with --gtest_output=xml.
//
// TAP_DSP_PARITY_MAX_N caps the sizes run, so the emulated QEMU legs can take
// the suite at N <= 4096 (Part 9, Part 10) where 2^20 does not fit in RAM;
// the CMake block and the toolchain files default it to 4096 when
// cross-compiling, and the 2^20 test is compiled out (not skipped) below it.
//
// Stage 2c moved fftsg.c to tests/reference/ooura/ (Decision D6) and
// fftsg_float.c WITH it, because this gate needs both: the float side
// compares against rdft_f, and without fftsg_float.c it would degenerate to
// port-vs-port. D6 retires the reference copy only after both MuTap and
// MuTap-Max pin a tree containing 2c; when it goes, this file goes with it
// (or its float half first, if only fftsg_float.c does) and the routing
// proof in tests/test_fft_routing.cpp plus the oracle in
// tests/test_fft_oracle.cpp remain.

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <limits>
#include <string>
#include <type_traits>
#include <vector>

#include <gtest/gtest.h>

#include "reference/ooura_rdft.h"
#include "support/signals.h"
#include "tap/dsp/fft/split_radix.h"

#ifndef TAP_DSP_PARITY_MAX_N
#define TAP_DSP_PARITY_MAX_N (1 << 20)
#endif
static_assert(TAP_DSP_PARITY_MAX_N >= 4 && (TAP_DSP_PARITY_MAX_N & (TAP_DSP_PARITY_MAX_N - 1)) == 0,
              "TAP_DSP_PARITY_MAX_N must be a power of two >= 4, or every sweep below is vacuously green");

namespace {

    // ------------------------------------------------------------------------
    // The reference side. Re-point here, nowhere else.
    //
    // Raw Ooura, called exactly as fftsg.c documents: ip[0] = 0 requests table
    // initialization on the first call, and the workspace geometry is the one
    // readme.txt prescribes (ip: 2 + sqrt(n/2), w: n/2). rdft / rdft_f are
    // declared by tests/reference/ooura_rdft.h and come from the reference C
    // library this target links (see the CMake block); the library itself has
    // not linked any C since Stage 2c.
    // ------------------------------------------------------------------------
    inline void ooura_rdft_ref(int n, int isgn, double* a, int* ip, double* w) {
        rdft(n, isgn, a, ip, w);
    }
    inline void ooura_rdft_ref(int n, int isgn, float* a, int* ip, float* w) {
        rdft_f(n, isgn, a, ip, w);
    }

    template <typename Sample>
    class ooura_ref {
      public:
        explicit ooura_ref(std::size_t n)
            : m_n(static_cast<int>(n))
            , m_ip(2 + static_cast<std::size_t>(std::sqrt(static_cast<double>(n) / 2.0)) + 1, 0)
            , m_w(n / 2, Sample(0)) {
            m_ip[0] = 0;
        }
        void forward_inplace(Sample* a) { ooura_rdft_ref(m_n, 1, a, m_ip.data(), m_w.data()); }
        void inverse_inplace(Sample* a) { ooura_rdft_ref(m_n, -1, a, m_ip.data(), m_w.data()); }

      private:
        int                 m_n;
        std::vector<int>    m_ip;
        std::vector<Sample> m_w;
    };

    // ------------------------------------------------------------------------
    // The side under test. Re-point here, nowhere else.
    //
    // Stage 2a (the port lands beside the C, nothing routed): re-pointed at
    //     tap::dsp::detail::split_radix_rdft<Sample>
    // Stage 2b (routing flipped): stays on the engine (see the file comment:
    // basic_real_fft<float> is a backend on two legs), and
    // tests/test_fft_routing.cpp pins basic_real_fft to the engine byte for
    // byte on the profiles that route to it.
    // Both expose the constructor-from-size / forward_inplace / inverse_inplace
    // surface Part 4 fixes, so nothing else in this file changes.
    // ------------------------------------------------------------------------
    template <typename Sample>
    using engine_under_test = tap::dsp::detail::split_radix_rdft<Sample>;

    // ------------------------------------------------------------------------
    // Materials.
    // ------------------------------------------------------------------------
    enum class material { broadband, on_bin_tone, impulse, dc, alternating };

    constexpr material k_materials[] = {material::broadband, material::on_bin_tone, material::impulse, material::dc,
                                        material::alternating};

    const char* name_of(material m) {
        switch (m) {
        case material::broadband:
            return "broadband";
        case material::on_bin_tone:
            return "on-bin tone";
        case material::impulse:
            return "impulse";
        case material::dc:
            return "dc";
        case material::alternating:
            return "alternating +-1";
        }
        return "?";
    }

    template <typename Sample>
    std::vector<Sample> make_material(material m, std::size_t n) {
        switch (m) {
        case material::broadband:
            return tap::dsp::test::random_signal<Sample>(n, 0x9E3779B9u);
        case material::on_bin_tone:
            // Bin n/8 (bin 1 below n = 8): far from DC and Nyquist so the
            // butterflies that combine it see a real twiddle, not +-1 or +-i.
            return tap::dsp::test::tone<Sample>(n, static_cast<double>(std::max<std::size_t>(1, n / 8)), 0.5, 0.3);
        case material::impulse: {
            std::vector<Sample> x(n, Sample(0));
            x[0] = Sample(1);
            return x;
        }
        case material::dc:
            return std::vector<Sample>(n, Sample(1));
        case material::alternating: {
            std::vector<Sample> x(n);
            for (std::size_t i = 0; i < n; ++i) {
                x[i] = (i % 2 == 0) ? Sample(1) : Sample(-1);
            }
            return x;
        }
        }
        return {};
    }

    // ------------------------------------------------------------------------
    // Measurement. ulp distance over the ordered-integer image of the bits
    // (sign-magnitude folded to two's complement), so +0/-0 count as 1 apart,
    // adjacent finite values as 1, and a NaN as "far" rather than as an
    // exception. Used only to make a mismatch report readable and to feed the
    // informational target; the gate itself is memcmp.
    // ------------------------------------------------------------------------
    template <typename Sample>
    std::uint64_t ulp_distance(Sample a, Sample b) {
        using bits_t = std::conditional_t<sizeof(Sample) == 4, std::int32_t, std::int64_t>;
        static_assert(sizeof(bits_t) == sizeof(Sample));
        bits_t ia = 0;
        bits_t ib = 0;
        std::memcpy(&ia, &a, sizeof(Sample));
        std::memcpy(&ib, &b, sizeof(Sample));
        if (ia < 0) {
            ia = std::numeric_limits<bits_t>::min() - ia;
        }
        if (ib < 0) {
            ib = std::numeric_limits<bits_t>::min() - ib;
        }
        // After the fold the images order like the values (negatives below
        // zero), so compare them signed; subtract in uint64, where the true
        // distance always fits (below 2^64 - 1 for the double image, whose
        // signed difference could overflow). Comparing the widened images
        // unsigned put a negative image above every positive one and reported
        // 2^64 - d for an opposite-sign pair (review A of tap/DspTap#32).
        const std::uint64_t ua = static_cast<std::uint64_t>(ia);
        const std::uint64_t ub = static_cast<std::uint64_t>(ib);
        return ia > ib ? ua - ub : ub - ua;
    }

    struct comparison {
        bool          identical  = true;
        std::size_t   first_diff = 0;
        std::uint64_t max_ulp    = 0;
    };

    template <typename Sample>
    comparison compare(const std::vector<Sample>& got, const std::vector<Sample>& ref) {
        comparison c;
        c.identical = std::memcmp(got.data(), ref.data(), got.size() * sizeof(Sample)) == 0;
        if (c.identical) {
            return c;
        }
        bool seen = false;
        for (std::size_t i = 0; i < got.size(); ++i) {
            const std::uint64_t d = ulp_distance(got[i], ref[i]);
            if (d != 0 && !seen) {
                c.first_diff = i;
                seen         = true;
            }
            c.max_ulp = std::max(c.max_ulp, d);
        }
        return c;
    }

    // Runs one (N, material) pair in both directions. The inverse input is the
    // reference forward's output, so both sides receive identical bits and the
    // inverse is judged on its own, not through the forward.
    template <typename Sample>
    struct pair_result {
        comparison forward;
        comparison inverse;
    };

    template <typename Sample>
    pair_result<Sample> run_pair(std::size_t n, material m) {
        const std::vector<Sample> x = make_material<Sample>(m, n);

        ooura_ref<Sample>         ref(n);
        engine_under_test<Sample> dut(n);
        std::vector<Sample>       ref_spec = x;
        std::vector<Sample>       dut_spec = x;
        ref.forward_inplace(ref_spec.data());
        dut.forward_inplace(dut_spec.data());

        std::vector<Sample> ref_time = ref_spec;
        std::vector<Sample> dut_time = ref_spec;
        ref.inverse_inplace(ref_time.data());
        dut.inverse_inplace(dut_time.data());

        return {compare(dut_spec, ref_spec), compare(dut_time, ref_time)};
    }

    template <typename Sample>
    const char* precision_name() {
        return sizeof(Sample) == 4 ? "float" : "double";
    }

    /// The sweep: every power of two from 4 to min(65536, TAP_DSP_PARITY_MAX_N).
    std::vector<std::size_t> sweep_sizes() {
        std::vector<std::size_t> sizes;
        for (std::size_t n = 4; n <= 65536 && n <= static_cast<std::size_t>(TAP_DSP_PARITY_MAX_N); n *= 2) {
            sizes.push_back(n);
        }
        return sizes;
    }

    // Used by the gate only when it is compiled in (see the #if below).
    [[maybe_unused]] constexpr std::size_t k_large_n = std::size_t{1} << 20;

#if !defined(TAP_DSP_PARITY_INFORMATIONAL)

    // ========================================================================
    // THE GATE: -ffp-contract=off on both sides, memcmp identity.
    // ========================================================================

    template <typename Sample>
    void expect_identical(std::size_t n, material m, bool forward) {
        const pair_result<Sample> r = run_pair<Sample>(n, m);
        const comparison&         c = forward ? r.forward : r.inverse;
        ASSERT_TRUE(c.identical) << (forward ? "forward" : "inverse") << " " << precision_name<Sample>() << " N=" << n
                                 << " material=" << name_of(m) << ": first difference at index " << c.first_diff
                                 << ", max " << c.max_ulp << " ulp. The engine under test is not the same sequence "
                                 << "of IEEE operations as Ooura's rdft; see the statement-fidelity rule in Part 4.";
    }

    template <typename Sample>
    void expect_identical_all_sizes(bool forward) {
        for (const std::size_t n : sweep_sizes()) {
            for (const material m : k_materials) {
                expect_identical<Sample>(n, m, forward);
                if (::testing::Test::HasFatalFailure()) {
                    return;
                }
            }
        }
    }

    TEST(fft_parity_ooura, ForwardIsBitIdenticalToOouraDouble) {
        expect_identical_all_sizes<double>(true);
    }

    TEST(fft_parity_ooura, ForwardIsBitIdenticalToOouraFloat) {
        expect_identical_all_sizes<float>(true);
    }

    TEST(fft_parity_ooura, InverseIsBitIdenticalToOouraDouble) {
        expect_identical_all_sizes<double>(false);
    }

    TEST(fft_parity_ooura, InverseIsBitIdenticalToOouraFloat) {
        expect_identical_all_sizes<float>(false);
    }

    // One run at 2^20 — the largest geometry any consumer uses (AmbiTap's
    // long-partition convolution) — kept as its own test so its runtime shows
    // separately. COMPILED OUT, not skipped, when TAP_DSP_PARITY_MAX_N is
    // below 2^20: the bare-metal one-shot main counts a GTEST_SKIP as a failed
    // gate (tests/bare_metal_main.cpp), and the QEMU legs run at 4096.
#if TAP_DSP_PARITY_MAX_N >= (1 << 20)
    template <typename Sample>
    void expect_identical_large() {
        for (const material m : k_materials) {
            expect_identical<Sample>(k_large_n, m, true);
            expect_identical<Sample>(k_large_n, m, false);
            if (::testing::Test::HasFatalFailure()) {
                return;
            }
        }
    }

    TEST(fft_parity_ooura, LargeTransformIsBitIdenticalToOouraDouble) {
        expect_identical_large<double>();
    }

    TEST(fft_parity_ooura, LargeTransformIsBitIdenticalToOouraFloat) {
        expect_identical_large<float>();
    }
#endif // TAP_DSP_PARITY_MAX_N >= (1 << 20)

#else // TAP_DSP_PARITY_INFORMATIONAL

    // ========================================================================
    // INFORMATIONAL: default flags on both sides, measured max ulp per N,
    // never a failure. Printed to stdout (visible under `ctest -V`, NOT under
    // --output-on-failure, which hides a passing test's output) and recorded
    // as gtest properties for the JUnit XML the CMake block requests. One
    // test for both precisions, so one process writes the whole XML.
    // ========================================================================

    template <typename Sample>
    void report_max_ulp() {
        std::vector<std::size_t> sizes = sweep_sizes();
        if (k_large_n <= static_cast<std::size_t>(TAP_DSP_PARITY_MAX_N)) {
            sizes.push_back(k_large_n);
        }
        std::uint64_t overall = 0;
        std::printf("[ ulp ] %s: engine under test vs Ooura, default compiler flags (informational)\n",
                    precision_name<Sample>());
        std::printf("[ ulp ] %10s %14s %14s   %s\n", "N", "forward", "inverse", "worst material (fwd / inv)");
        for (const std::size_t n : sizes) {
            std::uint64_t fwd       = 0;
            std::uint64_t inv       = 0;
            const char*   fwd_worst = name_of(k_materials[0]);
            const char*   inv_worst = name_of(k_materials[0]);
            for (const material m : k_materials) {
                const pair_result<Sample> r = run_pair<Sample>(n, m);
                if (r.forward.max_ulp > fwd) {
                    fwd       = r.forward.max_ulp;
                    fwd_worst = name_of(m);
                }
                if (r.inverse.max_ulp > inv) {
                    inv       = r.inverse.max_ulp;
                    inv_worst = name_of(m);
                }
            }
            overall = std::max({overall, fwd, inv});
            // %lu, not %zu: newlib-nano's printf on the QEMU legs prints "zu".
            std::printf("[ ulp ] %10lu %14llu %14llu   %s / %s\n", static_cast<unsigned long>(n),
                        static_cast<unsigned long long>(fwd), static_cast<unsigned long long>(inv), fwd_worst,
                        inv_worst);
            ::testing::Test::RecordProperty(std::string(precision_name<Sample>()) + "_N" + std::to_string(n) + "_fwd",
                                            std::to_string(fwd));
            ::testing::Test::RecordProperty(std::string(precision_name<Sample>()) + "_N" + std::to_string(n) + "_inv",
                                            std::to_string(inv));
        }
        std::printf("[ ulp ] %s: max over all N and materials = %llu ulp\n", precision_name<Sample>(),
                    static_cast<unsigned long long>(overall));
        std::fflush(stdout);
        SUCCEED() << "informational only: " << overall << " ulp max";
    }

    TEST(fft_parity_ooura_default_flags, ReportsMaxUlpVersusOoura) {
        report_max_ulp<double>();
        report_max_ulp<float>();
    }

#endif // TAP_DSP_PARITY_INFORMATIONAL

} // namespace
