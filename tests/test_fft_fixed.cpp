// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// THE FIXED-POINT CONTRACT BATTERY for tap::dsp::basic_real_fft over the Q15
// (std::int16_t) and Q31 (std::int32_t) profiles under both scaling policies
// (Stage 3b of docs/audit-fft-and-code-smells.md; Part 9 names this file and
// its promises). Written test-first against the API contract; the kernel
// (fft/fixed_point.h) is what has to pass it, and every number below is a
// contract point the kernel's header cites by test name.
//
// The contract, as the tests read it (fft.h and fft/fft_arith.h are the
// specification; this is the summary the assertions are written against):
//
//   - forward_inplace / inverse_inplace / forward / inverse return an
//     exponent e. Read the buffer as fractions of full scale (Q0.15 / Q0.31).
//     With G the double golden model (basic_real_fft<double>) on the SAME
//     quantised input read as fractions:
//         forward:  G.forward_inplace(x)               == out * 2^e
//         inverse:  G.inverse_inplace(a) (UNNORMALISED) == out * 2^e
//     up to the kernel's rounding noise, same packing, same exp(+i) sign.
//     The fixed-point inverse applies NO 2/N: the exponent carries the scale.
//   - scaling::fixed: e is the constant fixed_scaling_exponent(n) =
//     log2(n) + fft_arith<Sample>::k_fixed_scaling_input_pre_shift in both
//     directions: log2 n for Q15 (forward output exactly X/N), log2 n + 1 for
//     Q31 (the one-bit input pre-shift that buys the missing guard bits).
//   - scaling::block_floating: 0 <= e <= that constant, data-dependent; the
//     kernel shifts only when growth requires it.
//   - round trip, both policies: x == out * 2^(e_fwd + e_inv + 1 - log2 n)
//     up to rounding noise (Ooura's unnormalised inverse has gain N/2).
//   - saturation-free for every input under both policies; no alignment
//     requirement; allocation at construction only (test_fft_rt.cpp).
//
// Measured numbers. Every tolerance and ratio here is a number measured on the
// real kernel and pinned at 2x (the log_mel pattern), never a round number;
// the table `pins` carries them all in one place with the host, compiler and
// date they were taken on, and every pinned test prints what it measured so
// a -V run records the current value beside the pin. Fixed seeds, no wall clock, no filesystem, no
// <random>: the battery runs unchanged on the four QEMU legs (Part 10), where
// N <= 2048 fixed point is cheap and the double golden model is the expensive
// part, so sizes are kept modest and TAP_DSP_PARITY_MAX_N caps the sweeps.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <limits>
#include <numbers>
#include <type_traits>
#include <vector>

#include <gtest/gtest.h>

#include "support/signals.h"
#include "tap/dsp/fft.h"
#include "tap/dsp/fft/fft_arith.h"
#include "tap/dsp/fft/tables.h"

#ifndef TAP_DSP_PARITY_MAX_N
#define TAP_DSP_PARITY_MAX_N 65536
#endif

namespace {

    namespace scaling = tap::dsp::scaling;
    using tap::dsp::test::sample_scale;

    // ------------------------------------------------------------------------
    // The four configurations under test.
    // ------------------------------------------------------------------------
    template <typename Sample, typename Scaling>
    struct config {
        using sample = Sample;
        using policy = Scaling;
        using fft    = tap::dsp::basic_real_fft<Sample, Scaling>;
        using arith  = tap::dsp::fft_arith<Sample>;

        static constexpr bool   k_is_q15   = std::is_same_v<Sample, std::int16_t>;
        static constexpr bool   k_is_bfp   = std::is_same_v<Scaling, scaling::block_floating>;
        static constexpr double k_lsb      = sample_scale<Sample>::k_lsb;
        static constexpr double k_full     = sample_scale<Sample>::k_full_scale;
        static constexpr Sample k_max      = std::numeric_limits<Sample>::max();
        static constexpr Sample k_min      = std::numeric_limits<Sample>::min();
        static constexpr int    k_pre      = arith::k_fixed_scaling_input_pre_shift;
        static constexpr int    k_int_frac = k_is_q15 ? 29 : 31; ///< fraction bits of the int32 work format

        static const char* name() {
            if constexpr (k_is_q15) {
                return k_is_bfp ? "Q15/bfp  " : "Q15/fixed";
            }
            else {
                return k_is_bfp ? "Q31/bfp  " : "Q31/fixed";
            }
        }
    };

    using q15_fixed = config<std::int16_t, scaling::fixed>;
    using q31_fixed = config<std::int32_t, scaling::fixed>;
    using q15_bfp   = config<std::int16_t, scaling::block_floating>;
    using q31_bfp   = config<std::int32_t, scaling::block_floating>;

    using all_configs   = ::testing::Types<q15_fixed, q31_fixed, q15_bfp, q31_bfp>;
    using fixed_configs = ::testing::Types<q15_fixed, q31_fixed>;
    using bfp_configs   = ::testing::Types<q15_bfp, q31_bfp>;

    // The aliases are the documented spellings of the four configurations.
    static_assert(std::is_same_v<tap::dsp::real_fft_q15, q15_fixed::fft>);
    static_assert(std::is_same_v<tap::dsp::real_fft_q31, q31_fixed::fft>);
    static_assert(std::is_same_v<tap::dsp::real_fft_q15_bfp, q15_bfp::fft>);
    static_assert(std::is_same_v<tap::dsp::real_fft_q31_bfp, q31_bfp::fft>);
    // The default policy is fixed scaling.
    static_assert(std::is_same_v<tap::dsp::basic_real_fft<std::int16_t>, q15_fixed::fft>);

    // ------------------------------------------------------------------------
    // THE PINS. Measured on the real kernel and pinned at 2x; see each test
    // for what the number bounds. A pin of 0.0 on a field a row reads means
    // "not measured" and fails that test on purpose.
    // ------------------------------------------------------------------------
    struct pin_table {
        double saturation_max_lsb;   ///< SaturationFreeWorstCaseDoesNotWrap: max |out - G/2^e| in output LSB
        double round_trip_k;         ///< RoundTripReconstructsInputPerPolicy: max error in reconstructed LSB
        double bfp_vs_fixed_lsb;     ///< BfpMatchesFixedAfterShift: max |shifted bfp - fixed| in output LSB
        double negation_sum_max_lsb; ///< RoundingBiasOnNegatedInputIsBounded: max |F(x) + F(-x)| in LSB
        double negation_bias_lsb;    ///< ... and the mean of F(x) + F(-x) in LSB
        double noise_ratio_noise;    ///< NoiseFloorTracksWelchModel: max measured/model, white noise
        double noise_ratio_tone;     ///< ... on-bin tone
        double q15_vs_q31_lsb;       ///< Q15AndQ31AgreeToTheQ15Floor: max |v15 - v31| in Q15 output LSB
    };

    // Measured 2026-09-18 on x86-64 Linux (Ubuntu 24.04, glibc 2.39), GCC
    // 13.3.0 and clang 18.1.3 -O3 (identical: the kernel is integer
    // arithmetic and the tables are the same libm), kernel at
    // claude/wave2-stage3b-kernel a55a14f, and pinned at 2x. The measured
    // values, in the pins' order:
    //   Q15/fixed  0.500 / 0.562 /  -  / 1.00 / 0.0005 / 1.055 / 0.022 / 0.500
    //   Q31/fixed  4.250 / 3.382 /  -  / 6.00 / 0.8125 / 1.744 / 0.383 /  -
    //   Q15/bfp    0.750 / 5.000 / 1.0 / 1.00 / 0.0156 / 1.041 / 0.535 / 0.500
    //   Q31/bfp   15.988 / 41.50 / 4.0 / 62.0 / 1.0469 / 2.321 / 0.644 /  -
    // The Q31 block-floating maxima all sit at index 0 or 1: the DC/Nyquist
    // path, where round-half-up biases add coherently (fixed_point.h,
    // "Honest limit"). bfp_vs_fixed_lsb is read by the block-floating rows
    // only and q15_vs_q31_lsb by the Q15 rows only; the others carry 0.0.
    constexpr pin_table k_pins_q15_fixed{1.0, 1.13, 0.0, 2.0, 0.001, 2.11, 0.044, 1.0};
    constexpr pin_table k_pins_q31_fixed{8.5, 6.77, 0.0, 12.0, 1.63, 3.49, 0.77, 0.0};
    constexpr pin_table k_pins_q15_bfp{1.5, 10.0, 2.0, 2.0, 0.032, 2.09, 1.07, 1.0};
    constexpr pin_table k_pins_q31_bfp{32.0, 83.0, 8.0, 124.0, 2.1, 4.65, 1.29, 0.0};

    template <typename Cfg>
    constexpr const pin_table& pins() {
        if constexpr (std::is_same_v<Cfg, q15_fixed>) {
            return k_pins_q15_fixed;
        }
        else if constexpr (std::is_same_v<Cfg, q31_fixed>) {
            return k_pins_q31_fixed;
        }
        else if constexpr (std::is_same_v<Cfg, q15_bfp>) {
            return k_pins_q15_bfp;
        }
        else {
            return k_pins_q31_bfp;
        }
    }

    // ------------------------------------------------------------------------
    // Small helpers.
    // ------------------------------------------------------------------------
    constexpr std::size_t k_max_n = std::min<std::size_t>(65536, TAP_DSP_PARITY_MAX_N);

    constexpr int log2_size(std::size_t n) {
        int l = 0;
        while ((std::size_t{1} << l) < n) {
            ++l;
        }
        return l;
    }
    static_assert(log2_size(4) == 2 && log2_size(65536) == 16);

    double pow2(int e) {
        return std::ldexp(1.0, e);
    }

    double db(double power) {
        return 10.0 * std::log10(std::max(power, 1e-300));
    }

    template <typename Sample>
    std::vector<double> fractions(const std::vector<Sample>& x) {
        std::vector<double> d(x.size());
        for (std::size_t i = 0; i < x.size(); ++i) {
            d[i] = sample_scale<Sample>::to_double(x[i]);
        }
        return d;
    }

    template <typename Sample>
    std::vector<Sample> quantise(const std::vector<double>& x) {
        std::vector<Sample> s(x.size());
        for (std::size_t i = 0; i < x.size(); ++i) {
            s[i] = sample_scale<Sample>::from_double(x[i]);
        }
        return s;
    }

    /// The double golden model's unnormalised forward, on fractions.
    std::vector<double> golden_forward(std::vector<double> x) {
        tap::dsp::basic_real_fft<double> g(x.size());
        g.forward_inplace(x.data());
        return x;
    }

    /// The double golden model's UNNORMALISED inverse, on fractions.
    std::vector<double> golden_inverse(std::vector<double> a) {
        tap::dsp::basic_real_fft<double> g(a.size());
        g.inverse_inplace(a.data());
        return a;
    }

    template <typename Sample>
    struct transform_result {
        std::vector<Sample> out;
        int                 exponent;
    };

    template <typename Cfg>
    transform_result<typename Cfg::sample> run_forward(const std::vector<typename Cfg::sample>& x) {
        typename Cfg::fft                      fft(x.size());
        transform_result<typename Cfg::sample> r{x, 0};
        r.exponent = fft.forward_inplace(r.out.data());
        return r;
    }

    template <typename Cfg>
    transform_result<typename Cfg::sample> run_inverse(const std::vector<typename Cfg::sample>& a) {
        typename Cfg::fft                      fft(a.size());
        transform_result<typename Cfg::sample> r{a, 0};
        r.exponent = fft.inverse_inplace(r.out.data());
        return r;
    }

    /// |out - golden / 2^e| over the whole buffer, in output LSBs.
    struct deviation {
        double      max_lsb;
        double      rms_lsb;
        std::size_t argmax;
    };

    template <typename Cfg>
    deviation deviation_from_golden(const std::vector<typename Cfg::sample>& out, const std::vector<double>& golden,
                                    int e) {
        deviation    d{0.0, 0.0, 0};
        const double scale = pow2(-e) / Cfg::k_lsb; // golden -> output LSB units
        for (std::size_t i = 0; i < out.size(); ++i) {
            const double err = static_cast<double>(out[i]) - golden[i] * scale;
            if (std::fabs(err) > d.max_lsb) {
                d.max_lsb = std::fabs(err);
                d.argmax  = i;
            }
            d.rms_lsb += err * err;
        }
        d.rms_lsb = std::sqrt(d.rms_lsb / static_cast<double>(out.size()));
        return d;
    }

    /// Named test vectors.
    template <typename Sample>
    struct pattern {
        const char*         name;
        std::vector<Sample> x;
    };

    /// Full-scale adversarial time-domain patterns (fft_arith.h, "Scaling"):
    /// every one has |x| at a rail in every sample, so it drives the
    /// magnitude bound; the rotated packed pairs and the square-wave
    /// exponentials are the ones that meet a 45 degree twiddle (Part 6, N6).
    template <typename Sample>
    std::vector<pattern<Sample>> adversarial_time_patterns(std::size_t n) {
        constexpr Sample hi = std::numeric_limits<Sample>::max();
        constexpr Sample lo = std::numeric_limits<Sample>::min();

        std::vector<pattern<Sample>> p;
        p.push_back({"dc +full", std::vector<Sample>(n, hi)});
        p.push_back({"dc INT_MIN", std::vector<Sample>(n, lo)});
        {
            std::vector<Sample> x(n);
            for (std::size_t j = 0; j < n; ++j) {
                x[j] = (j % 2 == 0) ? hi : lo;
            }
            p.push_back({"nyquist +-", x});
            for (std::size_t j = 0; j < n; ++j) {
                x[j] = (j % 2 == 0) ? lo : hi;
            }
            p.push_back({"nyquist -+", x});
        }
        {
            std::vector<Sample> x(n, Sample{0});
            x[0] = hi;
            p.push_back({"impulse +full", x});
            x[0] = lo;
            p.push_back({"impulse INT_MIN", x});
        }
        // Packed pairs z_j = x[2j] + i x[2j+1] at a rail on BOTH components
        // (|z| = sqrt(2) full scale), rotating by +90 / -90 degrees per pair:
        // a complex exponential at bin M/4 of the complex sequence whose whole
        // energy lands in one bin, coherent through every butterfly.
        {
            const std::array<std::array<Sample, 2>, 4> cw{{{hi, hi}, {lo, hi}, {lo, lo}, {hi, lo}}};
            const std::array<std::array<Sample, 2>, 4> ccw{{{hi, hi}, {hi, lo}, {lo, lo}, {lo, hi}}};
            const std::array<std::array<Sample, 2>, 4> cw_neg{{{lo, lo}, {hi, lo}, {hi, hi}, {lo, hi}}};
            std::vector<Sample>                        x(n);
            for (std::size_t j = 0; j < n / 2; ++j) {
                x[2 * j]     = cw[j % 4][0];
                x[2 * j + 1] = cw[j % 4][1];
            }
            p.push_back({"pair rotation +90", x});
            for (std::size_t j = 0; j < n / 2; ++j) {
                x[2 * j]     = ccw[j % 4][0];
                x[2 * j + 1] = ccw[j % 4][1];
            }
            p.push_back({"pair rotation -90", x});
            for (std::size_t j = 0; j < n / 2; ++j) {
                x[2 * j]     = cw_neg[j % 4][0];
                x[2 * j + 1] = cw_neg[j % 4][1];
            }
            p.push_back({"pair rotation +90 from INT_MIN", x});
        }
        // Square-wave complex exponentials: z_j = sgn cos(theta_j + phi) +
        // i sgn sin(theta_j + phi), theta_j = 2 pi j k / M. Every component at
        // a rail; the fundamental sits at complex bin k where the stage
        // twiddles are 45 degrees (k = M/8), 22.5 degrees (M/16) and 135
        // degrees (3M/8). phi keeps the samples off the zero crossings.
        for (const double k_over_m : {1.0 / 8.0, 1.0 / 16.0, 3.0 / 8.0}) {
            std::vector<Sample> x(n);
            for (std::size_t j = 0; j < n / 2; ++j) {
                const double theta =
                    2.0 * std::numbers::pi * k_over_m * static_cast<double>(j) + std::numbers::pi / 8.0;
                x[2 * j]     = std::cos(theta) >= 0.0 ? hi : lo;
                x[2 * j + 1] = std::sin(theta) >= 0.0 ? hi : lo;
            }
            p.push_back({"square exponential", x});
        }
        // Real square waves at bins N/8 and N/4 + 1.
        for (const double k_over_n : {1.0 / 8.0, 1.0 / 4.0 + 1.0 / static_cast<double>(n)}) {
            std::vector<Sample> x(n);
            for (std::size_t j = 0; j < n; ++j) {
                const double theta = 2.0 * std::numbers::pi * k_over_n * static_cast<double>(j) + 0.3;
                x[j]               = std::cos(theta) >= 0.0 ? hi : lo;
            }
            p.push_back({"square wave", x});
        }
        // Full-scale binary noise: the broadband worst case (fixed seeds).
        for (const std::uint32_t seed : {0x2545F491u, 0x9E3779B9u, 0x1D872B41u, 0xC0FFEE01u}) {
            tap::dsp::test::xorshift32 rng(seed);
            std::vector<Sample>        x(n);
            for (std::size_t j = 0; j < n; ++j) {
                x[j] = (rng.next_u32() & 0x10000u) != 0u ? hi : lo;
            }
            p.push_back({"binary noise", x});
        }
        return p;
    }

    /// Full-scale adversarial packed spectra for the inverse direction.
    template <typename Sample>
    std::vector<pattern<Sample>> adversarial_spectrum_patterns(std::size_t n) {
        constexpr Sample hi = std::numeric_limits<Sample>::max();
        constexpr Sample lo = std::numeric_limits<Sample>::min();

        // The time-domain set is a fine set of spectra too (all at a rail),
        // plus the spectral shapes whose inverse concentrates: one full-scale
        // complex bin (a tone at amplitude sqrt(2)), and DC + Nyquist only.
        std::vector<pattern<Sample>> p = adversarial_time_patterns<Sample>(n);
        {
            std::vector<Sample> a(n, Sample{0});
            const std::size_t   k = std::max<std::size_t>(1, n / 8 + 1);
            a[2 * k]              = hi;
            a[2 * k + 1]          = hi;
            p.push_back({"one full complex bin", a});
            a.assign(n, Sample{0});
            a[0] = hi;
            a[1] = hi;
            p.push_back({"dc and nyquist full", a});
            a[0] = lo;
            a[1] = lo;
            p.push_back({"dc and nyquist INT_MIN", a});
        }
        return p;
    }

    /// Sizes for the per-pattern sweeps: the two smallest geometries (M = 2
    /// and M = 4: no radix-4 stage at all, one radix-2), an odd and an even
    /// log2 M, and the two log_mel / pvoc geometries.
    constexpr std::array<std::size_t, 7> k_sweep_sizes{4, 8, 16, 64, 512, 1024, 2048};

    // ========================================================================
    // Exponents.
    // ========================================================================

    template <typename Cfg>
    class fft_fixed_scaling_test : public ::testing::Test {};
    TYPED_TEST_SUITE(fft_fixed_scaling_test, fixed_configs);

    template <typename Cfg>
    class fft_fixed_bfp_test : public ::testing::Test {};
    TYPED_TEST_SUITE(fft_fixed_bfp_test, bfp_configs);

    template <typename Cfg>
    class fft_fixed_point_test : public ::testing::Test {};
    TYPED_TEST_SUITE(fft_fixed_point_test, all_configs);

    // scaling::fixed: e = log2(n) + k_fixed_scaling_input_pre_shift, a
    // compile-time constant, in both directions, from every entry point.
    TYPED_TEST(fft_fixed_scaling_test, FixedExponentIsTheStatedConstant) {
        using cfg = TypeParam;
        using fft = typename cfg::fft;
        static_assert(fft::fixed_scaling_exponent(1024) == 10 + cfg::k_pre);
        static_assert(fft::fixed_scaling_exponent(4) == 2 + cfg::k_pre);
        static_assert(noexcept(fft::fixed_scaling_exponent(4)));
        if constexpr (cfg::k_is_q15) {
            static_assert(fft::fixed_scaling_exponent(512) == 9, "Q15: no pre-shift, output exactly X/N");
        }
        else {
            static_assert(fft::fixed_scaling_exponent(512) == 10, "Q31: the one-bit input pre-shift");
        }
        for (std::size_t n = 4; n <= k_max_n; n *= 2) {
            const int expected = log2_size(n) + cfg::k_pre;
            EXPECT_EQ(fft::fixed_scaling_exponent(n), expected) << "n=" << n;
            fft        f(n);
            const auto x =
                tap::dsp::test::random_signal<typename cfg::sample>(n, 0x2545F491u ^ static_cast<std::uint32_t>(n));
            auto a = x;
            EXPECT_EQ(f.forward_inplace(a.data()), expected) << "forward_inplace n=" << n;
            EXPECT_EQ(f.inverse_inplace(a.data()), expected) << "inverse_inplace n=" << n;
            std::vector<typename cfg::sample> out(n);
            EXPECT_EQ(f.forward(x.data(), out.data()), expected) << "forward n=" << n;
            EXPECT_EQ(f.inverse(out.data(), out.data()), expected) << "inverse (aliased) n=" << n;
            // And a silent block reports the same constant: the exponent is
            // static, not a headroom measurement.
            std::vector<typename cfg::sample> zeros(n, 0);
            EXPECT_EQ(f.forward_inplace(zeros.data()), expected) << "forward of silence n=" << n;
        }
    }

    // scaling::block_floating: 0 <= e <= the fixed constant, in both
    // directions, whatever the input; a louder input never needs a smaller
    // exponent than silence does (0 is the floor).
    TYPED_TEST(fft_fixed_bfp_test, BfpExponentIsWithinRange) {
        using cfg = TypeParam;
        using fft = typename cfg::fft;
        using s   = typename cfg::sample;
        for (const std::size_t n : {std::size_t{4}, std::size_t{16}, std::size_t{512}, std::size_t{2048}}) {
            const int top = fft::fixed_scaling_exponent(n);
            fft       f(n);

            std::vector<std::vector<s>> inputs;
            inputs.push_back(std::vector<s>(n, s{0}));
            inputs.push_back(std::vector<s>(n, s{0}));
            inputs.back()[0] = s{1};
            inputs.push_back(tap::dsp::test::random_signal<s>(n, 0x9E3779B9u, 1e-3));
            inputs.push_back(tap::dsp::test::random_signal<s>(n, 0x9E3779B9u, 0.1));
            inputs.push_back(tap::dsp::test::random_signal<s>(n, 0x9E3779B9u, 1.0));
            inputs.push_back(std::vector<s>(n, cfg::k_max));
            inputs.push_back(std::vector<s>(n, cfg::k_min));
            for (auto& p : adversarial_time_patterns<s>(n)) {
                inputs.push_back(std::move(p.x));
            }

            for (std::size_t i = 0; i < inputs.size(); ++i) {
                auto      a  = inputs[i];
                const int ef = f.forward_inplace(a.data());
                EXPECT_GE(ef, 0) << "forward n=" << n << " input " << i;
                EXPECT_LE(ef, top) << "forward n=" << n << " input " << i;
                const int ei = f.inverse_inplace(a.data());
                EXPECT_GE(ei, 0) << "inverse n=" << n << " input " << i;
                EXPECT_LE(ei, top) << "inverse n=" << n << " input " << i;
            }
        }
    }

    // ========================================================================
    // Exact scale under fixed scaling.
    // ========================================================================

    // An integer-valued impulse of A output LSBs has X_k = A for every k; a
    // constant of c LSBs has X_0 = N c and nothing else. Under fixed scaling
    // the output is X * 2^-e with e = log2 n + pre: when that is an integer
    // number of LSBs no rounding can occur anywhere in the kernel (every
    // stage shift discards zero bits: the impulse's value returns to A / 4^s
    // exactly and the constant's to c at every stage), so the result is
    // bit-exact. A = m * 2^e and c a multiple of 2^(pre + 2) (the Q31
    // pre-shift plus one radix-4 shift of the constant path; Q15's widen
    // supplies fourteen zero bits so any int16 works there).
    TYPED_TEST(fft_fixed_scaling_test, FixedForwardScaleIsExactlyXOverN) {
        using cfg = TypeParam;
        using s   = typename cfg::sample;
        for (const std::size_t n : k_sweep_sizes) {
            const int e = cfg::fft::fixed_scaling_exponent(n);
            // Impulses: A = m * 2^e for several m, the largest that fits.
            const std::int64_t unit  = std::int64_t{1} << e;
            const std::int64_t top_m = static_cast<std::int64_t>(cfg::k_max) / unit;
            for (const std::int64_t m : {std::int64_t{1}, std::int64_t{3}, top_m, -top_m}) {
                if (m == 0) {
                    continue;
                }
                std::vector<s> x(n, s{0});
                x[0]         = static_cast<s>(m * unit);
                const auto r = run_forward<cfg>(x);
                ASSERT_EQ(r.exponent, e);
                // Flat: DC, Nyquist and every real part m; every imaginary part 0.
                for (std::size_t i = 0; i < n; ++i) {
                    const std::int64_t expected = (i < 2 || i % 2 == 0) ? m : 0;
                    ASSERT_EQ(static_cast<std::int64_t>(r.out[i]), expected)
                        << "impulse A=" << m * unit << " n=" << n << " index " << i;
                }
            }
            // Constants: c a multiple of 2^(pre + 2); the extreme values too.
            const std::int64_t        step = std::int64_t{1} << (cfg::k_pre + 2);
            std::vector<std::int64_t> constants{step, -step, 5 * step,
                                                (static_cast<std::int64_t>(cfg::k_max) / step) * step,
                                                static_cast<std::int64_t>(cfg::k_min)};
            if constexpr (cfg::k_is_q15) {
                constants.push_back(cfg::k_max); // 32767: odd, exact through the 14-bit widen
                constants.push_back(1);
                constants.push_back(-1);
            }
            for (const std::int64_t c : constants) {
                const std::vector<s> x(n, static_cast<s>(c));
                const auto           r = run_forward<cfg>(x);
                ASSERT_EQ(r.exponent, e);
                // X_0 = N c; out_0 = N c / 2^e = c / 2^pre.
                const std::int64_t expected_dc = c / (std::int64_t{1} << cfg::k_pre);
                ASSERT_EQ(static_cast<std::int64_t>(r.out[0]), expected_dc) << "dc c=" << c << " n=" << n;
                for (std::size_t i = 1; i < n; ++i) {
                    ASSERT_EQ(r.out[i], s{0}) << "dc c=" << c << " n=" << n << " index " << i;
                }
            }
        }
    }

    // The Q31 profile's fixed forward is X / (2N): the extra bit is the input
    // pre-shift (fft_arith<std::int32_t>::k_fixed_scaling_input_pre_shift = 1),
    // stated here in the profile's own numbers.
    TEST(fft_fixed_q31, FixedForwardScaleIsExactlyXOverTwoN) {
        using fft = tap::dsp::real_fft_q31;
        static_assert(tap::dsp::fft_arith<std::int32_t>::k_fixed_scaling_input_pre_shift == 1);
        static_assert(tap::dsp::fft_arith<std::int16_t>::k_fixed_scaling_input_pre_shift == 0);
        static_assert(fft::fixed_scaling_exponent(512) == 10);
        constexpr std::size_t n = 512;
        fft                   f(n);
        // A constant of 2^30 (0.5): X_0 = 512 * 2^30, out_0 = 2^29 (0.25).
        std::vector<std::int32_t> dc(n, std::int32_t{1} << 30);
        EXPECT_EQ(f.forward_inplace(dc.data()), 10);
        EXPECT_EQ(dc[0], std::int32_t{1} << 29);
        // An impulse of 2^20: out = 2^20 / 1024 = 2^10 in every slot.
        std::vector<std::int32_t> imp(n, 0);
        imp[0] = std::int32_t{1} << 20;
        EXPECT_EQ(f.forward_inplace(imp.data()), 10);
        for (std::size_t i = 0; i < n; ++i) {
            ASSERT_EQ(imp[i], (i < 2 || i % 2 == 0) ? std::int32_t{1} << 10 : 0) << i;
        }
        // The Q15 profile, for contrast, is X / N: the same 0.5 constant
        // comes back as 0.5.
        tap::dsp::real_fft_q15    f15(n);
        std::vector<std::int16_t> dc15(n, std::int16_t{1} << 14);
        EXPECT_EQ(f15.forward_inplace(dc15.data()), 9);
        EXPECT_EQ(dc15[0], std::int16_t{1} << 14);
    }

    // The inverse under fixed scaling carries the same constant and is
    // exact on the same kind of input: the unnormalised inverse of a
    // DC-only spectrum a[0] = c is c/2 everywhere (fft.h's packing), so the
    // output is c / 2^(e+1); of a flat spectrum (the impulse response) it is
    // N c / 2 at k = 0 and 0 elsewhere, so out_0 = c / 2^(pre + 1).
    TYPED_TEST(fft_fixed_scaling_test, FixedInverseCarriesTheSameExponent) {
        using cfg = TypeParam;
        using s   = typename cfg::sample;
        for (const std::size_t n : k_sweep_sizes) {
            const int          e    = cfg::fft::fixed_scaling_exponent(n);
            const std::int64_t unit = std::int64_t{1} << (e + 1);
            for (const std::int64_t m :
                 {std::int64_t{1}, std::int64_t{-3}, static_cast<std::int64_t>(cfg::k_max) / unit}) {
                if (m == 0) {
                    continue;
                }
                std::vector<s> a(n, s{0});
                a[0]         = static_cast<s>(m * unit);
                const auto r = run_inverse<cfg>(a);
                ASSERT_EQ(r.exponent, e);
                for (std::size_t i = 0; i < n; ++i) {
                    ASSERT_EQ(static_cast<std::int64_t>(r.out[i]), m) << "dc spectrum n=" << n << " index " << i;
                }
            }
            // The flat spectrum: a[0] = a[1] = a[2k] = c, a[2k+1] = 0. The
            // pre-pass takes pre + 1 bits, then the constant c / 2^(pre+1)
            // rides the unrotated path through every stage's 2-bit shift, so
            // c must carry pre + 3 zero bits for the trip to be exact (with
            // fewer, round-half-up doubles a 2-LSB constant at every stage:
            // that is the kernel's honest bias, not a scale error).
            const std::int64_t step = std::int64_t{1} << (cfg::k_pre + 3);
            for (const std::int64_t c : {step, -step, (static_cast<std::int64_t>(cfg::k_max) / step) * step}) {
                std::vector<s> a(n, s{0});
                a[0] = static_cast<s>(c);
                a[1] = static_cast<s>(c);
                for (std::size_t k = 1; k < n / 2; ++k) {
                    a[2 * k] = static_cast<s>(c);
                }
                const auto r = run_inverse<cfg>(a);
                ASSERT_EQ(r.exponent, e);
                // x_0 = N c / 2, out_0 = N c / 2^(e + 1) = c / 2^(pre + 1).
                ASSERT_EQ(static_cast<std::int64_t>(r.out[0]), c / (std::int64_t{1} << (cfg::k_pre + 1)))
                    << "flat spectrum c=" << c << " n=" << n;
                for (std::size_t i = 1; i < n; ++i) {
                    ASSERT_EQ(r.out[i], s{0}) << "flat spectrum c=" << c << " n=" << n << " index " << i;
                }
            }
        }
    }

    // ========================================================================
    // Saturation.
    // ========================================================================

    // Every adversarial full-scale pattern, both directions, both policies:
    // the output tracks the golden model at the returned exponent to within
    // the pinned number of LSBs (a wrap is 2^15 / 2^31 LSBs off, a clamp at
    // least the amount clamped), and no output sample sits on a rail unless
    // the golden model puts it within the pin of that rail.
    struct worst_case {
        double      value = 0.0;
        const char* what  = "";
        const char* name  = "";
        std::size_t n     = 0;
        int         e     = 0;
        std::size_t index = 0;
        void        note(double v, const char* w, const char* nm, std::size_t size, int exponent, std::size_t i) {
            if (v > value) {
                value = v;
                what  = w;
                name  = nm;
                n     = size;
                e     = exponent;
                index = i;
            }
        }
    };

    TYPED_TEST(fft_fixed_point_test, SaturationFreeWorstCaseDoesNotWrap) {
        using cfg      = TypeParam;
        using s        = typename cfg::sample;
        const auto pin = pins<cfg>().saturation_max_lsb;
        worst_case worst;
        worst_case rail; // largest shortfall of the golden model below a rail the output sits on
        for (const std::size_t n : k_sweep_sizes) {
            for (const bool inverse : {false, true}) {
                const auto patterns = inverse ? adversarial_spectrum_patterns<s>(n) : adversarial_time_patterns<s>(n);
                for (const auto& p : patterns) {
                    const auto r      = inverse ? run_inverse<cfg>(p.x) : run_forward<cfg>(p.x);
                    const auto golden = inverse ? golden_inverse(fractions(p.x)) : golden_forward(fractions(p.x));
                    const auto d      = deviation_from_golden<cfg>(r.out, golden, r.exponent);
                    worst.note(d.max_lsb, inverse ? "inverse" : "forward", p.name, n, r.exponent, d.argmax);
                    for (std::size_t i = 0; i < n; ++i) {
                        if (r.out[i] == cfg::k_max || r.out[i] == cfg::k_min) {
                            const double g = std::fabs(golden[i]) * pow2(-r.exponent) / cfg::k_lsb;
                            rail.note(static_cast<double>(cfg::k_max) - g, inverse ? "inverse" : "forward", p.name, n,
                                      r.exponent, i);
                        }
                    }
                }
            }
        }
        std::printf("[ measured ] %s saturation sweep: max |out - G/2^e| = %.3f LSB (pin %.3f) at %s %s n=%zu e=%d "
                    "index %zu; largest rail shortfall %.3f LSB\n",
                    cfg::name(), worst.value, pin, worst.what, worst.name, worst.n, worst.e, worst.index, rail.value);
        EXPECT_GT(pin, 0.0) << "unmeasured pin";
        EXPECT_LE(worst.value, pin) << cfg::name() << " " << worst.what << " " << worst.name << " n=" << worst.n
                                    << " e=" << worst.e << " index " << worst.index;
        // No output sits on a rail unless the golden model is within the pin
        // of that rail: a clamp would show as a shortfall of at least what was
        // clamped, a wrap as 2^15 / 2^31 LSB of deviation above.
        EXPECT_LE(rail.value, pin) << cfg::name() << " " << rail.what << " " << rail.name << " n=" << rail.n
                                   << ": sample " << rail.index << " sits on a rail the golden model does not reach";
    }

    // A silent block is silent out, in both directions.
    TYPED_TEST(fft_fixed_point_test, SilenceIsSilence) {
        using cfg = TypeParam;
        using s   = typename cfg::sample;
        for (const std::size_t n : k_sweep_sizes) {
            const std::vector<s> zeros(n, s{0});
            const auto           f = run_forward<cfg>(zeros);
            const auto           i = run_inverse<cfg>(zeros);
            for (std::size_t j = 0; j < n; ++j) {
                ASSERT_EQ(f.out[j], s{0}) << "forward n=" << n << " index " << j;
                ASSERT_EQ(i.out[j], s{0}) << "inverse n=" << n << " index " << j;
            }
            EXPECT_GE(f.exponent, 0);
            EXPECT_LE(f.exponent, cfg::fft::fixed_scaling_exponent(n));
        }
    }

    // ========================================================================
    // Round trip.
    // ========================================================================

    // x == out * 2^(e_fwd + e_inv + 1 - log2 n) up to rounding noise. The
    // error is measured in units of the reconstructed LSB, i.e. the output
    // LSB times that same power of two: under fixed scaling the round trip
    // discards log2(N) bits twice and the unit is coarse (Q15, N = 1024: 2^-4
    // of full scale, an honest number the header states), under block
    // floating point the exponents are smaller and the unit finer. The pin
    // is the largest error in that unit over the sizes and levels below.
    TYPED_TEST(fft_fixed_point_test, RoundTripReconstructsInputPerPolicy) {
        using cfg      = TypeParam;
        using s        = typename cfg::sample;
        const auto pin = pins<cfg>().round_trip_k;
        worst_case worst;
        for (const std::size_t n : {std::size_t{16}, std::size_t{64}, std::size_t{256}, std::size_t{1024}}) {
            for (const double amplitude : {1.0, 0.25, 0.01}) {
                const auto x =
                    tap::dsp::test::random_signal<s>(n, 0xC0FFEE01u ^ static_cast<std::uint32_t>(n), amplitude);
                const auto   xf    = fractions(x);
                const auto   fwd   = run_forward<cfg>(x);
                const auto   inv   = run_inverse<cfg>(fwd.out);
                const int    shift = fwd.exponent + inv.exponent + 1 - log2_size(n);
                const double unit  = cfg::k_lsb * pow2(shift);
                for (std::size_t i = 0; i < n; ++i) {
                    const double back = sample_scale<s>::to_double(inv.out[i]) * pow2(shift);
                    worst.note(std::fabs(back - xf[i]) / unit,
                               amplitude == 1.0 ? "0 dBFS" : (amplitude == 0.25 ? "-12 dBFS" : "-40 dBFS"),
                               "round trip", n, fwd.exponent * 100 + inv.exponent, i);
                }
            }
        }
        std::printf("[ measured ] %s round trip: max error %.3f reconstructed LSB (pin %.3f) at %s n=%zu "
                    "e_fwd=%d e_inv=%d index %zu\n",
                    cfg::name(), worst.value, pin, worst.what, worst.n, worst.e / 100, worst.e % 100, worst.index);
        EXPECT_GT(pin, 0.0) << "unmeasured pin";
        EXPECT_LE(worst.value, pin) << cfg::name() << " " << worst.what << " n=" << worst.n
                                    << " e_fwd=" << worst.e / 100 << " e_inv=" << worst.e % 100 << " index "
                                    << worst.index;
    }

    // forward()/inverse() are copy-then-in-place: same bits, same exponent,
    // with and without aliasing; and the fixed-point inverse() applies no
    // 2/N (the exponent carries the scale).
    TYPED_TEST(fft_fixed_point_test, OutOfPlaceIsCopyThenInPlace) {
        using cfg = TypeParam;
        using s   = typename cfg::sample;
        for (const std::size_t n : {std::size_t{8}, std::size_t{512}}) {
            typename cfg::fft fft(n);
            const auto        x = tap::dsp::test::random_signal<s>(n, 0x1D872B41u, 0.9);

            auto           inplace   = x;
            const int      e_inplace = fft.forward_inplace(inplace.data());
            std::vector<s> out(n, s{0});
            const int      e_out   = fft.forward(x.data(), out.data());
            auto           alias   = x;
            const int      e_alias = fft.forward(alias.data(), alias.data());
            EXPECT_EQ(e_out, e_inplace);
            EXPECT_EQ(e_alias, e_inplace);
            EXPECT_EQ(std::memcmp(out.data(), inplace.data(), n * sizeof(s)), 0) << "forward out-of-place n=" << n;
            EXPECT_EQ(std::memcmp(alias.data(), inplace.data(), n * sizeof(s)), 0) << "forward aliased n=" << n;

            auto           inv_inplace = inplace;
            const int      ei_inplace  = fft.inverse_inplace(inv_inplace.data());
            std::vector<s> inv_out(n, s{0});
            const int      ei_out    = fft.inverse(inplace.data(), inv_out.data());
            auto           inv_alias = inplace;
            const int      ei_alias  = fft.inverse(inv_alias.data(), inv_alias.data());
            EXPECT_EQ(ei_out, ei_inplace);
            EXPECT_EQ(ei_alias, ei_inplace);
            EXPECT_EQ(std::memcmp(inv_out.data(), inv_inplace.data(), n * sizeof(s)), 0)
                << "inverse out-of-place n=" << n;
            EXPECT_EQ(std::memcmp(inv_alias.data(), inv_inplace.data(), n * sizeof(s)), 0) << "inverse aliased n=" << n;
        }
    }

    // ========================================================================
    // Block floating point against fixed scaling.
    // ========================================================================

    // The BFP output, shifted right (round-half-up) by e_fixed - e_bfp,
    // agrees with the fixed output to within the pinned number of LSBs (the
    // difference is the fixed path's extra rounding noise plus one rounding
    // of the shift); and on an input where BFP reports the full constant it
    // shifted like fixed at every stage, so the two are bit-identical.
    TYPED_TEST(fft_fixed_bfp_test, BfpMatchesFixedAfterShift) {
        using cfg      = TypeParam;
        using s        = typename cfg::sample;
        using fixed    = config<s, scaling::fixed>;
        const auto pin = pins<cfg>().bfp_vs_fixed_lsb;
        worst_case worst;
        for (const std::size_t n : k_sweep_sizes) {
            std::vector<pattern<s>> inputs;
            inputs.push_back({"noise 0 dBFS", tap::dsp::test::random_signal<s>(n, 0x2545F491u, 1.0)});
            inputs.push_back({"noise -20 dBFS", tap::dsp::test::random_signal<s>(n, 0x2545F491u, 0.1)});
            inputs.push_back({"noise -40 dBFS", tap::dsp::test::random_signal<s>(n, 0x2545F491u, 0.01)});
            inputs.push_back({"tone -6 dBFS", tap::dsp::test::tone<s>(n, static_cast<double>(n / 8 + 1), 0.5, 0.3)});
            inputs.push_back({"dc +full", std::vector<s>(n, cfg::k_max)});
            for (const auto& p : inputs) {
                for (const bool inverse : {false, true}) {
                    const auto f = inverse ? run_inverse<fixed>(p.x) : run_forward<fixed>(p.x);
                    const auto b = inverse ? run_inverse<cfg>(p.x) : run_forward<cfg>(p.x);
                    ASSERT_LE(b.exponent, f.exponent) << p.name << " n=" << n;
                    const int shift = f.exponent - b.exponent;
                    for (std::size_t i = 0; i < n; ++i) {
                        const std::int64_t bv = static_cast<std::int64_t>(b.out[i]);
                        const std::int64_t shifted =
                            shift == 0 ? bv : ((bv + (std::int64_t{1} << (shift - 1))) >> shift);
                        const double diff =
                            std::fabs(static_cast<double>(shifted - static_cast<std::int64_t>(f.out[i])));
                        worst.note(diff, inverse ? "inverse" : "forward", p.name, n, f.exponent * 100 + b.exponent, i);
                    }
                }
            }
        }
        std::printf("[ measured ] %s vs fixed after shift: max %.1f LSB (pin %.1f) at %s %s n=%zu e_fixed=%d e_bfp=%d "
                    "index %zu\n",
                    cfg::name(), worst.value, pin, worst.what, worst.name, worst.n, worst.e / 100, worst.e % 100,
                    worst.index);
        EXPECT_GT(pin, 0.0) << "unmeasured pin";
        EXPECT_LE(worst.value, pin) << cfg::name() << " " << worst.what << " " << worst.name << " n=" << worst.n
                                    << " e_fixed=" << worst.e / 100 << " e_bfp=" << worst.e % 100 << " index "
                                    << worst.index;
    }

    TYPED_TEST(fft_fixed_bfp_test, BfpAtTheFullExponentIsBitIdenticalToFixed) {
        using cfg            = TypeParam;
        using s              = typename cfg::sample;
        using fixed          = config<s, scaling::fixed>;
        std::size_t attained = 0;
        std::size_t tried    = 0;
        for (const std::size_t n : k_sweep_sizes) {
            for (const auto& p : adversarial_time_patterns<s>(n)) {
                ++tried;
                const auto f = run_forward<fixed>(p.x);
                const auto b = run_forward<cfg>(p.x);
                if (b.exponent != f.exponent) {
                    continue;
                }
                ++attained;
                EXPECT_EQ(std::memcmp(b.out.data(), f.out.data(), n * sizeof(s)), 0)
                    << cfg::name() << " " << p.name << " n=" << n << ": BFP reported the full exponent " << f.exponent
                    << " but its output is not the fixed output";
            }
        }
        // The contract's upper bound is attainable: some full-scale input
        // makes BFP shift like fixed at every stage.
        EXPECT_GT(attained, 0u) << cfg::name() << ": no full-scale pattern out of " << tried
                                << " reached the fixed exponent";
        std::printf("[ measured ] %s: %zu of %zu full-scale patterns reach the fixed exponent\n", cfg::name(), attained,
                    tried);
    }

    // ========================================================================
    // Rounding symmetry.
    // ========================================================================

    // Every rounding in the kernel is round-half-up (fft_arith.h), which is
    // not odd-symmetric: F(-x) is not exactly -F(x). The sum F(x) + F(-x) is
    // the asymmetry, at most one LSB per rounding that lands on a tie and
    // reaches the output; its maximum and its mean (the bias) are pinned.
    TYPED_TEST(fft_fixed_point_test, RoundingBiasOnNegatedInputIsBounded) {
        using cfg           = TypeParam;
        using s             = typename cfg::sample;
        const auto max_pin  = pins<cfg>().negation_sum_max_lsb;
        const auto bias_pin = pins<cfg>().negation_bias_lsb;
        worst_case worst_max;
        worst_case worst_bias;
        for (const std::size_t n : k_sweep_sizes) {
            for (const double amplitude : {0.999, 0.1}) {
                // 0.999 keeps -x representable (-INT_MIN would saturate).
                const auto x =
                    tap::dsp::test::random_signal<s>(n, 0x9E3779B9u ^ static_cast<std::uint32_t>(n), amplitude);
                std::vector<s> neg(n);
                for (std::size_t i = 0; i < n; ++i) {
                    neg[i] = static_cast<s>(-x[i]);
                }
                for (const bool inverse : {false, true}) {
                    const auto a = inverse ? run_inverse<cfg>(x) : run_forward<cfg>(x);
                    const auto b = inverse ? run_inverse<cfg>(neg) : run_forward<cfg>(neg);
                    ASSERT_EQ(a.exponent, b.exponent) << "negation changed the exponent, n=" << n;
                    double mean = 0.0;
                    for (std::size_t i = 0; i < n; ++i) {
                        const double sum = static_cast<double>(a.out[i]) + static_cast<double>(b.out[i]);
                        worst_max.note(std::fabs(sum), inverse ? "inverse" : "forward",
                                       amplitude > 0.5 ? "-0 dBFS" : "-20 dBFS", n, a.exponent, i);
                        mean += sum;
                    }
                    // The bias is a mean over the block: read it where the
                    // block is large enough for a mean to say something.
                    if (n >= 64) {
                        worst_bias.note(std::fabs(mean / static_cast<double>(n)), inverse ? "inverse" : "forward",
                                        amplitude > 0.5 ? "-0 dBFS" : "-20 dBFS", n, a.exponent, 0);
                    }
                }
            }
        }
        std::printf("[ measured ] %s F(x)+F(-x): max %.2f LSB (pin %.2f) at %s %s n=%zu index %zu; bias %.4f LSB "
                    "(pin %.4f) at %s %s n=%zu\n",
                    cfg::name(), worst_max.value, max_pin, worst_max.what, worst_max.name, worst_max.n, worst_max.index,
                    worst_bias.value, bias_pin, worst_bias.what, worst_bias.name, worst_bias.n);
        EXPECT_GT(max_pin, 0.0) << "unmeasured pin";
        EXPECT_GT(bias_pin, 0.0) << "unmeasured pin";
        EXPECT_LE(worst_max.value, max_pin) << cfg::name() << " " << worst_max.what << " " << worst_max.name
                                            << " n=" << worst_max.n << " index " << worst_max.index;
        EXPECT_LE(worst_bias.value, bias_pin)
            << cfg::name() << " " << worst_bias.what << " " << worst_bias.name << " n=" << worst_bias.n;
    }

    // ========================================================================
    // Q15 against Q31: the one sibling-profile comparison (Part 9, Rules).
    // ========================================================================

    // The same signal (a Q15 vector, widened to Q31 exactly) through both
    // profiles under the same policy: brought to a common exponent, the two
    // spectra differ by the Q15 output quantisation plus the Q15 profile's
    // (coarser) internal noise, pinned in Q15 output LSBs.
    template <typename Scaling>
    void q15_vs_q31(const char* what, double pin) {
        using c15 = config<std::int16_t, Scaling>;
        using c31 = config<std::int32_t, Scaling>;
        worst_case worst;
        for (const std::size_t n : k_sweep_sizes) {
            for (const double amplitude : {1.0, 0.05}) {
                const auto x15 = tap::dsp::test::random_signal<std::int16_t>(
                    n, 0xC0FFEE01u ^ static_cast<std::uint32_t>(n), amplitude);
                std::vector<std::int32_t> x31(n);
                for (std::size_t i = 0; i < n; ++i) {
                    x31[i] = static_cast<std::int32_t>(x15[i]) << 16; // exact
                }
                for (const bool inverse : {false, true}) {
                    const auto r15 = inverse ? run_inverse<c15>(x15) : run_forward<c15>(x15);
                    const auto r31 = inverse ? run_inverse<c31>(x31) : run_forward<c31>(x31);
                    // Compare in the unnormalised frame, in Q15 output LSBs at e15.
                    const double unit = c15::k_lsb * pow2(r15.exponent);
                    for (std::size_t i = 0; i < n; ++i) {
                        const double v15 = sample_scale<std::int16_t>::to_double(r15.out[i]) * pow2(r15.exponent);
                        const double v31 = sample_scale<std::int32_t>::to_double(r31.out[i]) * pow2(r31.exponent);
                        worst.note(std::fabs(v15 - v31) / unit, inverse ? "inverse" : "forward",
                                   amplitude == 1.0 ? "0 dBFS" : "-26 dBFS", n, r15.exponent * 100 + r31.exponent, i);
                    }
                }
            }
        }
        std::printf(
            "[ measured ] Q15 vs Q31 (%s): max %.3f Q15 LSB (pin %.3f) at %s %s n=%zu e15=%d e31=%d index %zu\n", what,
            worst.value, pin, worst.what, worst.name, worst.n, worst.e / 100, worst.e % 100, worst.index);
        EXPECT_GT(pin, 0.0) << "unmeasured pin";
        EXPECT_LE(worst.value, pin) << what << " " << worst.what << " " << worst.name << " n=" << worst.n
                                    << " e15=" << worst.e / 100 << " e31=" << worst.e % 100 << " index " << worst.index;
    }

    TEST(fft_fixed_profiles, Q15AndQ31AgreeToTheQ15Floor) {
        q15_vs_q31<scaling::fixed>("fixed", k_pins_q15_fixed.q15_vs_q31_lsb);
        q15_vs_q31<scaling::block_floating>("bfp", k_pins_q15_bfp.q15_vs_q31_lsb);
    }

    // ========================================================================
    // The noise floor against Welch's model.
    // ========================================================================
    //
    // THE MODEL (Welch 1969, "A fixed-point fast Fourier transform error
    // analysis"; Oppenheim & Weinstein 1972), specialised to this kernel's
    // arithmetic as fft_arith.h states it and to the structure fixed_point.h
    // documents (the model is written here, independently, from those two
    // headers; the kernel's design note carries its own derivation):
    //
    //   - Internal LSB q: 2^-31 for Q31, 2^-29 for the widened Q15 (Q2.29),
    //     both as fractions of full scale. Every rounding is additive noise,
    //     independent of every other: a mul_coeff rounding is uniform on a
    //     continuous interval, variance q^2/12; an s-bit shr_round of an
    //     integer takes one of 2^s equally likely residues, variance
    //     (1 - 4^-s) q^2/12 (3/4 of q^2/12 for one bit, 15/16 for two).
    //   - Structure: a radix-4 DIF complex FFT of length M = N/2 (one stage
    //     per factor of 4, spans M, M/4, ..., down to 4; a radix-2 stage when
    //     log2 M is odd), then Ooura's real post-pass pairing bins k and
    //     M - k with one complex product per pair (bins 1 .. N/4 - 1 and
    //     their mirrors: every interior bin but N/4). Under fixed scaling the
    //     Q31 profile folds its one-bit input pre-shift into the first
    //     stage's shift (a single 3-bit rounding).
    //   - Per stage, per output component. Shift-before-butterfly: each of
    //     the stage's sum_gain inputs (4 for radix-4, 2 for radix-2) is
    //     rounded once per component and reaches every output component with
    //     unit weight, so the shift injects sum_gain * v(s) -- when the stage
    //     shifts at all (shr_round(x, 0) is exact, so a BFP stage shifting by
    //     0 injects nothing). The twiddle rotation is the two-rounding
    //     complex multiply, 2 q^2/12 per component, on 3 of the 4 outputs of
    //     every butterfly whose index q is not 0: averaged over a stage of
    //     span L that is 2 * (3/4) * (1 - 4/L) q^2/12 per component, zero for
    //     the last radix-4 stage (L = 4) and for the radix-2 stage. The
    //     twiddle's own quantisation, |dw| <= 2^-31 per component (half an
    //     LSB of Q1.30), uniform, adds 2 * P * 2^-62 / 3 on the rotated
    //     outputs, P being the signal variance per component entering the
    //     rotation.
    //   - Post-pass: its one-bit shift injects v(1) (the two paired values'
    //     roundings reach the output with weights |1 - w|^2 + |w|^2 = 1),
    //     its rotation 2 q^2/12; it carries earlier noise and the signal with
    //     gain 1 * 4^-shift.
    //   - Propagation: noise present after a stage reaches the output through
    //     each later stage u with power gain sum_gain_u * 4^-shift_u (a sum of
    //     sum_gain_u uncorrelated terms through unit-magnitude twiddles, then
    //     the shift). Signal: a white input of variance sigma^2 per sample
    //     enters as M complex values of variance sigma^2 per component and
    //     follows the same gains without the injections; through the
    //     post-pass it lands at sigma^2 / (2N) * 4^(log2 N - e), the DFT's
    //     N sigma^2 / 2 per component at exponent e.
    //   - The Q15 narrowing adds (2^-15)^2 / 12 per output value; under block
    //     floating point one more shift (growth 1) precedes it.
    //   - Shift schedule. scaling::fixed: every stage shifts its full amount
    //     and the sum is fixed_scaling_exponent(N). scaling::block_floating:
    //     the kernel returns the total e and the per-stage schedule is not
    //     observable, so the model brackets it: the EARLY schedule (e
    //     assigned to the first stages, each up to its fixed amount) is the
    //     worst case for noise, since a bit shifted early attenuates nothing
    //     injected after it, and is what the pin is measured against; the
    //     LATE schedule is the best case and the table prints it beside.
    //
    // The model's number is the noise variance per output component in the
    // output's own units (fractions of full scale at exponent e). The
    // measurement is the mean over the interior bins (indices 2 .. N-1; DC
    // and Nyquist take the glue path) of |out - G(x)/2^e|^2 per component,
    // with G the double golden model on the SAME quantised input, so only
    // the kernel's roundings and its twiddle quantisation are in it. What
    // the variance model leaves out is the round-half-up bias: an s-bit
    // shr_round has mean +2^-(s+1) LSB, which the final shifts put on every
    // bin and the earlier ones spread with rotated phases; it shows as a
    // measured/model ratio above one and is pinned as such. The pin is the
    // largest measured/model ratio over the sizes and levels, per
    // configuration and material, at 2x; the table is printed for a -V run
    // to record (Part 13).
    struct stage_spec {
        int    sum_gain;          ///< inputs summed into each output component
        int    max_shift;         ///< the fixed-scaling shift, and the BFP maximum
        double rotation_fraction; ///< fraction of output components that see a twiddle rotation
    };

    std::vector<stage_spec> kernel_stages(std::size_t n, bool fold_pre_shift, bool q15_final_shift) {
        std::vector<stage_spec> stages;
        const std::size_t       m = n / 2;
        for (std::size_t span = m; span >= 4; span /= 4) {
            const double rotated = 0.75 * (1.0 - 4.0 / static_cast<double>(span)); // 3 of 4 outputs, q != 0
            stages.push_back({4, 2, rotated});
        }
        if (stages.empty() || (m >> (2 * stages.size())) == 2) {
            stages.push_back({2, 1, 0.0}); // the radix-2 final stage (log2 M odd)
        }
        if (fold_pre_shift) {
            stages.front().max_shift += 1; // Q31 fixed: pre-shift and first stage, one rounding
        }
        stages.push_back({1, 1, 1.0}); // the real post-pass
        if (q15_final_shift) {
            stages.push_back({1, 31, 0.0}); // Q15 BFP: the shift before the narrow, unbounded
        }
        return stages;
    }

    enum class schedule { fixed, early, late };

    std::vector<int> shifts_for(const std::vector<stage_spec>& stages, int e, schedule which) {
        std::vector<int> shifts(stages.size(), 0);
        if (which == schedule::fixed) {
            for (std::size_t t = 0; t < stages.size(); ++t) {
                shifts[t] = stages[t].max_shift;
            }
            return shifts;
        }
        int remaining = e;
        if (which == schedule::early) {
            for (std::size_t t = 0; t < stages.size() && remaining > 0; ++t) {
                shifts[t] = std::min(stages[t].max_shift, remaining);
                remaining -= shifts[t];
            }
        }
        else {
            for (std::size_t t = stages.size(); t-- > 0 && remaining > 0;) {
                shifts[t] = std::min(stages[t].max_shift, remaining);
                remaining -= shifts[t];
            }
        }
        return shifts;
    }

    struct welch_prediction {
        double noise;  ///< variance per output component, output units
        double signal; ///< the model's signal variance per output component
    };

    welch_prediction welch_model(const std::vector<stage_spec>& stages, const std::vector<int>& shifts, double q_int,
                                 double q_out_extra, double input_variance) {
        constexpr double twiddle_lsb = 0x1p-31; // half an LSB of Q1.30
        const double     rounding    = q_int * q_int / 12.0;
        double           noise       = 0.0;
        double           signal      = input_variance;
        for (std::size_t t = 0; t < stages.size(); ++t) {
            const stage_spec& st         = stages[t];
            const double      shift_gain = pow2(-2 * shifts[t]);
            const double      stage_gain = static_cast<double>(st.sum_gain) * shift_gain;
            noise                        = noise * stage_gain;
            signal                       = signal * stage_gain;
            if (shifts[t] > 0) {
                noise += static_cast<double>(st.sum_gain) * (1.0 - pow2(-2 * shifts[t])) * rounding;
            }
            if (st.rotation_fraction > 0.0) {
                noise += st.rotation_fraction * (2.0 * rounding + 2.0 * signal * twiddle_lsb * twiddle_lsb / 3.0);
            }
        }
        noise += q_out_extra * q_out_extra / 12.0;
        return {noise, signal};
    }

    struct noise_row {
        const char* material;
        std::size_t n;
        double      level_db;
        int         exponent;
        double      signal_db;
        double      floor_db;
        double      model_db;
        double      model_late_db;
        double      ratio;
    };

    template <typename Cfg>
    noise_row measure_noise_floor(const char* material, std::size_t n, double level_db,
                                  const std::vector<typename Cfg::sample>& x, double input_variance) {
        const auto   r      = run_forward<Cfg>(x);
        const auto   golden = golden_forward(fractions(x));
        const double scale  = pow2(-r.exponent);

        double noise  = 0.0;
        double signal = 0.0;
        for (std::size_t i = 2; i < n; ++i) {
            const double g   = golden[i] * scale;
            const double err = sample_scale<typename Cfg::sample>::to_double(r.out[i]) - g;
            noise += err * err;
            signal += g * g;
        }
        const auto count = static_cast<double>(n - 2);
        noise /= count;
        signal /= count;

        const double q_int        = pow2(-Cfg::k_int_frac);
        const double q_out_extra  = Cfg::k_is_q15 ? Cfg::k_lsb : 0.0;
        const auto   stages       = kernel_stages(n, !Cfg::k_is_bfp && Cfg::k_pre > 0, Cfg::k_is_bfp && Cfg::k_is_q15);
        const auto   model_shifts = shifts_for(stages, r.exponent, Cfg::k_is_bfp ? schedule::early : schedule::fixed);
        const auto   late_shifts  = shifts_for(stages, r.exponent, Cfg::k_is_bfp ? schedule::late : schedule::fixed);
        const auto   model        = welch_model(stages, model_shifts, q_int, q_out_extra, input_variance);
        const auto   late         = welch_model(stages, late_shifts, q_int, q_out_extra, input_variance);

        noise_row row{material,           n,         level_db,        r.exponent,
                      db(signal),         db(noise), db(model.noise), db(late.noise),
                      noise / model.noise};
        std::printf("[ floor ] %s %-5s N=%5zu %4.0f dBFS e=%2d  signal %7.2f  floor %7.2f  model %7.2f (late %7.2f) "
                    "dBFS/component  ratio %.3f  snr %6.2f dB\n",
                    Cfg::name(), material, n, level_db, r.exponent, row.signal_db, row.floor_db, row.model_db,
                    row.model_late_db, row.ratio, row.signal_db - row.floor_db);
        return row;
    }

    TYPED_TEST(fft_fixed_point_test, NoiseFloorTracksWelchModel) {
        using cfg              = TypeParam;
        using s                = typename cfg::sample;
        const auto noise_pin   = pins<cfg>().noise_ratio_noise;
        const auto tone_pin    = pins<cfg>().noise_ratio_tone;
        double     worst_noise = 0.0;
        double     worst_tone  = 0.0;
        for (const std::size_t n : {std::size_t{256}, std::size_t{512}, std::size_t{2048}}) {
            for (const double level_db : {0.0, -20.0, -40.0, -60.0}) {
                const double amplitude = std::pow(10.0, level_db / 20.0);
                // White noise, uniform in [-A, A): variance A^2 / 3 per sample.
                const auto noise =
                    tap::dsp::test::random_signal<s>(n, 0x2545F491u ^ static_cast<std::uint32_t>(n), amplitude);
                const auto nrow = measure_noise_floor<cfg>("noise", n, level_db, noise, amplitude * amplitude / 3.0);
                worst_noise     = std::max(worst_noise, nrow.ratio);
                // On-bin tone at bin N/8 + 3, peak A * full scale: variance A^2 / 2.
                const double a    = amplitude * cfg::k_full;
                const auto   tone = tap::dsp::test::tone<s>(n, static_cast<double>(n / 8 + 3), a, 0.3);
                const auto   trow = measure_noise_floor<cfg>("tone", n, level_db, tone, a * a / 2.0);
                worst_tone        = std::max(worst_tone, trow.ratio);
            }
        }
        std::printf("[ measured ] %s noise/model ratio: noise %.3f (pin %.3f), tone %.3f (pin %.3f)\n", cfg::name(),
                    worst_noise, noise_pin, worst_tone, tone_pin);
        EXPECT_GT(noise_pin, 0.0) << "unmeasured pin";
        EXPECT_GT(tone_pin, 0.0) << "unmeasured pin";
        EXPECT_LE(worst_noise, noise_pin) << cfg::name() << " white-noise floor above the pinned ratio to the model";
        EXPECT_LE(worst_tone, tone_pin) << cfg::name() << " on-bin tone floor above the pinned ratio to the model";
        // A floor far BELOW the model would mean a different arithmetic (a
        // fused complex multiply, an exact shift) or a broken measurement,
        // not a better kernel: the two-rounding form is the contract.
        EXPECT_GE(worst_noise, noise_pin / 8.0) << cfg::name() << " white-noise floor implausibly far below the model";
    }

    // ========================================================================
    // The tables (fft/tables.h).
    // ========================================================================

    /// FNV-1a 64 over a table as little-endian int32 bytes, written here
    /// independently of detail::table_checksum so that helper is itself
    /// pinned: the fold is over bit patterns, host-endianness independent,
    /// and one differing last bit anywhere changes it (Part 13: an integer
    /// fold, never a float).
    std::uint64_t fnv1a64(const std::vector<std::int32_t>& table) {
        std::uint64_t h = 0xcbf29ce484222325ull;
        for (const std::int32_t v : table) {
            const auto u = static_cast<std::uint32_t>(v);
            for (int b = 0; b < 4; ++b) {
                h ^= static_cast<std::uint64_t>((u >> (8 * b)) & 0xffu);
                h *= 0x100000001b3ull;
            }
        }
        return h;
    }

    // |w_q - w| <= 0.5 LSB of Q1.30 on every host, for the kernel twiddles
    // (W_M^k, M = N/2, interleaved cos/sin) and the real post-pass pairs
    // (0.5 - 0.5 sin, 0.5 cos over N): each is cos/sin from libm rounded
    // once by make_coeff (half away from zero), never a recurrence. The
    // reference is libm too, on the exact turn fraction; two double
    // evaluations of the same angle differ by ~1e-16, far inside the 2^-20
    // LSB slack, so a last-bit libm difference passes here (the checksum
    // below is where it is detected) while a recurrence-drifted or
    // Q1.14-derived table fails. The exact entries (1, 0, i) and the
    // k <-> M-k symmetry are checked bit-exactly.
    TEST(fft_fixed_tables, TwiddleTableIsWithinHalfLsb) {
        constexpr double one   = 0x1p30;
        constexpr double slack = 0.5 + 0x1p-20;
        for (std::size_t n = 4; n <= k_max_n; n *= 2) {
            const std::size_t m     = n / 2;
            const auto        table = tap::dsp::detail::make_twiddle_table(m);
            ASSERT_EQ(table.size(), 2 * m) << "n=" << n;
            double worst = 0.0;
            for (std::size_t k = 0; k < m; ++k) {
                const double turns = static_cast<double>(k) / static_cast<double>(m); // exact
                const double c     = std::cos(2.0 * std::numbers::pi * turns) * one;
                const double s     = std::sin(2.0 * std::numbers::pi * turns) * one;
                const double ec    = std::fabs(static_cast<double>(table[2 * k]) - c);
                const double es    = std::fabs(static_cast<double>(table[2 * k + 1]) - s);
                worst              = std::max({worst, ec, es});
                ASSERT_LE(ec, slack) << "cos m=" << m << " k=" << k;
                ASSERT_LE(es, slack) << "sin m=" << m << " k=" << k;
            }
            EXPECT_EQ(table[0], std::int32_t{1} << 30) << "m=" << m;
            EXPECT_EQ(table[1], 0) << "m=" << m;
            if (m >= 4) {
                EXPECT_EQ(table[2 * (m / 4)], 0) << "m=" << m;
                EXPECT_EQ(table[2 * (m / 4) + 1], std::int32_t{1} << 30) << "m=" << m;
            }
            for (std::size_t k = 1; k < m / 2; ++k) {
                ASSERT_EQ(table[2 * k], table[2 * (m - k)]) << "cos symmetry m=" << m << " k=" << k;
                ASSERT_EQ(table[2 * k + 1], -table[2 * (m - k) + 1]) << "sin symmetry m=" << m << " k=" << k;
            }

            const auto post = tap::dsp::detail::make_real_post_pass_table(n);
            ASSERT_EQ(post.size(), 2 * (n / 4)) << "n=" << n;
            EXPECT_EQ(post[0], std::int32_t{1} << 29) << "n=" << n; // 0.5 - 0.5 sin 0
            EXPECT_EQ(post[1], std::int32_t{1} << 29) << "n=" << n; // 0.5 cos 0
            for (std::size_t k = 0; k < n / 4; ++k) {
                const double turns = static_cast<double>(k) / static_cast<double>(n);
                const double wkr   = (0.5 - 0.5 * std::sin(2.0 * std::numbers::pi * turns)) * one;
                const double wki   = 0.5 * std::cos(2.0 * std::numbers::pi * turns) * one;
                const double er    = std::fabs(static_cast<double>(post[2 * k]) - wkr);
                const double ei    = std::fabs(static_cast<double>(post[2 * k + 1]) - wki);
                worst              = std::max({worst, er, ei});
                ASSERT_LE(er, slack) << "wkr n=" << n << " k=" << k;
                ASSERT_LE(ei, slack) << "wki n=" << n << " k=" << k;
            }
            if (n <= 2048) {
                std::printf("[ measured ] tables n=%zu: max |w_q - w| = %.6f LSB\n", n, worst);
            }
        }
    }

    // Host libm last-bit differences can move a double lying within 2^-31 of
    // a Q1.30 rounding boundary to the other side (fft_arith.h, "Twiddles"),
    // and then fixed-point outputs differ between hosts by design. The
    // checksum makes that visible instead of absorbed: a disagreement between
    // the CI hosts is a finding to record (which host, which N), not a skip.
    // Both tables a transform of N builds are pinned, for the three
    // certified geometries.
    TEST(fft_fixed_tables, TwiddleTableChecksumIsPinned) {
        struct pinned {
            std::size_t   n;
            std::uint64_t twiddles;  ///< make_twiddle_table(n / 2)
            std::uint64_t post_pass; ///< make_real_post_pass_table(n)
        };
        // Taken 2026-09-18 on x86-64 Linux, glibc 2.39 (GCC 13.3.0 and clang
        // 18.1.3 agree), kernel a55a14f. The other CI hosts (macOS arm64,
        // Windows UCRT, the newlib QEMU legs) either reproduce these or the
        // difference is a recorded finding per host and N.
        constexpr std::array<pinned, 3> expected{{{256, 0x95f5c68afe494835ull, 0x66c84a75861eafb6ull},
                                                  {512, 0x6df6ff99a3ed3c85ull, 0x4b1374200abed27cull},
                                                  {2048, 0xe42c528f3ae88b45ull, 0xa0f40e80bbf4efb9ull}}};
        for (const auto& p : expected) {
            const auto twiddles = tap::dsp::detail::make_twiddle_table(p.n / 2);
            const auto post     = tap::dsp::detail::make_real_post_pass_table(p.n);
            const auto tw_sum   = fnv1a64(twiddles);
            const auto post_sum = fnv1a64(post);
            // Printed before any assertion so a -V log carries every host's
            // values whether or not they match.
            std::printf("[ checksum ] n=%zu twiddles fnv1a64=%016llx post-pass fnv1a64=%016llx\n", p.n,
                        static_cast<unsigned long long>(tw_sum), static_cast<unsigned long long>(post_sum));
            // The kernel's own checksum helper is the same fold.
            EXPECT_EQ(tap::dsp::detail::table_checksum(twiddles.data(), twiddles.size()), tw_sum);
            EXPECT_EQ(tap::dsp::detail::table_checksum(post.data(), post.size()), post_sum);
            EXPECT_NE(p.twiddles, 0u) << "unmeasured checksum n=" << p.n;
            EXPECT_EQ(tw_sum, p.twiddles) << "n=" << p.n << ": this host's libm produced a different twiddle table";
            EXPECT_EQ(post_sum, p.post_pass)
                << "n=" << p.n << ": this host's libm produced a different post-pass table";
        }
    }

} // namespace
