// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// THE INDEPENDENT ORACLE for tap::dsp::basic_real_fft (Stage 2a of
// docs/audit-fft-and-code-smells.md; Part 9 names this file).
//
// test_fft_parity_ooura.cpp proves the port IS Ooura, bit for bit. It cannot
// prove Ooura is a Fourier transform: if the reference C and the port shared a
// defect, parity would be green. This file answers that with two references
// that share nothing with fftsg.c:
//
//   1. Closed-form vectors whose transforms are known exactly — an impulse
//      (flat spectrum), DC (N in slot 0), the Nyquist alternation (N in slot
//      1), an on-bin cosine and sine at bin k (+N/2 in the real, resp. IMAG,
//      part of bin k — the plus sign in the imaginary part is the documented
//      W = exp(+2*pi*i/N) convention, the conjugate of the engineering DFT),
//      and a two-tone superposition — and their unnormalized inverses, whose
//      gain is N/2 per the contract in fft.h. Checked to a tolerance DERIVED
//      from the sample type's epsilon and log2 N (see "Tolerance" below).
//
//   2. A compensated-summation DFT for N <= 256: the O(N^2) definition with
//      double-double accumulation (TwoSum / TwoProd, Dekker 1971 and Knuth
//      TAOCP vol. 2 4.2.2; Shewchuk 1997 for the error-free transformations)
//      and double-double twiddles from a Taylor series after exact octant
//      reduction on the integer bin*sample index. Not `long double`: that is
//      `double` on MSVC and on Apple arm64 (Part 6, item N19), so it would be
//      no oracle at all on two of the three hosted CI legs.
//
// Relation to test_fft.cpp. That file's contract battery already pins the
// impulse, DC/Nyquist packing and the +i sign convention at one size each
// (N = 64 and 128) with a round tolerance; Part 9 lists the closed forms in
// both files on purpose. What this file adds is the sweep to N = 65536, the
// inverse closed forms, the derived (not round) tolerance, and the DFT
// reference; the test names are distinct so the ctest listing carries each
// promise once.
//
// The suite is typed over float and double now. The fixed-point stage
// (Stage 3, Part 7) extends `profile<Sample>` below with int16_t and int32_t;
// the tests are written against that trait (scale as a function of N, the
// profile's own tolerance, a full-scale amplitude below saturation), so the
// change is confined to the trait and the type list.

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numbers>
#include <vector>

#include <gtest/gtest.h>

#include "support/signals.h"
#include "tap/dsp/fft.h"

#ifndef TAP_DSP_PARITY_MAX_N
#define TAP_DSP_PARITY_MAX_N (1 << 20)
#endif

namespace {

    // ------------------------------------------------------------------------
    // Tolerance.
    //
    // Higham, Accuracy and Stability of Numerical Algorithms (2nd ed., 2002),
    // Theorem 24.2: a radix-2 FFT of length N = 2^t whose twiddles carry
    // relative error at most mu computes y_hat with
    //
    //     || y_hat - y ||_2  <=  t * eta / (1 - t * eta) * || y ||_2,
    //     eta = mu + gamma_4 (1 + mu),  gamma_4 ~ 4u,
    //
    // where u is the unit roundoff (epsilon / 2). mu is NOT one rounding here:
    // Ooura's makewt takes cos/sin from libm for a quarter of the table and
    // derives the rest arithmetically (w[2] = 0.5 / cos(2 delta), the halving
    // recurrences 0.5 / wk1r, ...), and the float instantiation forms
    // delta * j in float, so mu is a small multiple of u. With mu ~ 2-3u,
    // eta ~ 6-7u and, for t * eta << 1, the 2-norm error is below about
    // 7 * u * log2(N) * ||y||_2. The largest single-element error cannot exceed
    // the 2-norm, and ||y||_2 comes free with the exact answer. Casting the
    // closed-form input to Sample adds at most u * ||y||_2 (Parseval), and
    // Ooura is split-radix rather than the theorem's radix-2, which changes the
    // constant by a factor of order one. The bound is taken as
    //
    //     k_higham_constant * epsilon * log2(N) * ||y||_2,  k_higham_constant = 4
    //
    // (epsilon = 2u, so 8u: the ~7u above plus the input cast), i.e. derived
    // with slack, not fitted. Measured against it on x86-64 Linux (GCC 13, -O3,
    // glibc, Ooura as the engine): the largest |error| / (epsilon * log2(N) *
    // ||y||_2) over every test in this file is 0.25 for double and 0.17 for
    // float, 16x and 23x inside the constant.
    //
    // What the bound cannot see. A 2-norm used per element is loose by up to
    // sqrt(N): at float, N = 65536, DC input, the tolerance is ~0.5 on a
    // spectrum whose one live bin is 65536, so a single empty bin that is off
    // by 0.4 would pass. That is acceptable for the closed forms (their job is
    // the sweep over sizes and the sign convention) and is exactly why the DFT
    // comparison is confined to N <= 256, where sqrt(N) <= 16 and the reference
    // is exact to ~1e-30: there the check is tight enough to catch a single
    // wrong bin. Precision pins as measured numbers live in test_fft.cpp and
    // test_fft_backend.cpp; this file is an oracle for correctness.
    // ------------------------------------------------------------------------
    constexpr double k_higham_constant = 4.0;

    double higham_tolerance(double epsilon, std::size_t n, double spectrum_norm2) {
        const double log2n = std::log2(static_cast<double>(n));
        return k_higham_constant * epsilon * log2n * spectrum_norm2;
    }

    // ------------------------------------------------------------------------
    // Per-profile traits — THE EXTENSION POINT FOR THE FIXED-POINT STAGE.
    //
    // A profile says how a sample crosses into the double domain, the amplitude
    // the closed forms are driven at, what scale the engine's forward and
    // unnormalized inverse carry relative to the mathematical DFT as a function
    // of N, and its own tolerance. For the float profiles: full scale 1.0, both
    // scales 1, and the Higham bound above. Stage 3 (Part 7) adds
    //     template <> struct profile<std::int16_t> { ... };
    //     template <> struct profile<std::int32_t> { ... };
    // with k_full_scale below 1 - 2^-15 (a full-scale 1.0 input whose X/N
    // expectation is exactly 1.0 saturates in Q15), forward_scale(n) = 1/n
    // under `fixed` scaling (BFP reports its exponent alongside), the profile's
    // inverse scale, and a tolerance built from Part 7's quantization-noise
    // numbers rather than from epsilon. Then append the types to oracle_types.
    // ------------------------------------------------------------------------
    template <typename Sample>
    struct profile;

    // Size range a profile's engine accepts. fft.h promises any power of two
    // >= 4 for Ooura, and the double profile is always Ooura. The float
    // profile is whatever backend the build selected: under TAP_DSP_FFT_CMSIS
    // (the M55 leg) CMSIS-DSP's arm_rfft_fast_init_f32 accepts 32..4096 only
    // (arm_rfft_fast_init_f32.c, the switch at the end) and fft.h's wrapper
    // does not check its return status, so a size outside that range is
    // undefined behaviour — a HardFault at N = 4 on the QEMU M55 leg is how
    // this was found. Until fft.h rejects or falls back on those sizes (a
    // finding for the fft.h owner, not this file), the float sweeps on that
    // backend run over the range the backend supports. vDSP (macOS) takes the
    // full range.
    constexpr std::size_t k_ooura_min_n = 4;
    constexpr std::size_t k_ooura_max_n = std::size_t{1} << 20;
#if defined(TAP_DSP_FFT_CMSIS)
    constexpr std::size_t k_float_backend_min_n = 32;
    constexpr std::size_t k_float_backend_max_n = 4096;
#else
    constexpr std::size_t k_float_backend_min_n = k_ooura_min_n;
    constexpr std::size_t k_float_backend_max_n = k_ooura_max_n;
#endif

    template <>
    struct profile<double> {
        static constexpr double      k_epsilon    = std::numeric_limits<double>::epsilon();
        static constexpr double      k_full_scale = 1.0; ///< closed-form drive amplitude
        static constexpr std::size_t k_min_n      = k_ooura_min_n;
        static constexpr std::size_t k_max_n      = k_ooura_max_n;
        static double                to_double(double v) { return v; }
        static double                from_double(double v) { return v; }
        static double                forward_scale(std::size_t) { return 1.0; } ///< engine forward = scale * DFT
        static double inverse_scale(std::size_t) { return 1.0; } ///< engine inverse = scale * unnormalized
        static double tolerance(std::size_t n, double norm2) { return higham_tolerance(k_epsilon, n, norm2); }
    };

    template <>
    struct profile<float> {
        static constexpr double      k_epsilon    = std::numeric_limits<float>::epsilon();
        static constexpr double      k_full_scale = 1.0;
        static constexpr std::size_t k_min_n      = k_float_backend_min_n;
        static constexpr std::size_t k_max_n      = k_float_backend_max_n;
        static double                to_double(float v) { return static_cast<double>(v); }
        static float                 from_double(double v) { return static_cast<float>(v); }
        static double                forward_scale(std::size_t) { return 1.0; }
        static double                inverse_scale(std::size_t) { return 1.0; }
        static double tolerance(std::size_t n, double norm2) { return higham_tolerance(k_epsilon, n, norm2); }
    };

    using oracle_types = ::testing::Types<float, double>;

    /// ||y||_2 of the full complex spectrum of x, from Parseval: sqrt(N) * ||x||_2.
    double spectrum_norm2(const std::vector<double>& x) {
        double e = 0.0;
        for (const double v : x) {
            e += v * v;
        }
        return std::sqrt(static_cast<double>(x.size()) * e);
    }

    std::vector<double> scale_by(std::vector<double> x, double factor) {
        for (double& v : x) {
            v *= factor;
        }
        return x;
    }

    template <typename Sample>
    std::vector<double> to_doubles(const std::vector<Sample>& x) {
        std::vector<double> d(x.size());
        for (std::size_t i = 0; i < x.size(); ++i) {
            d[i] = profile<Sample>::to_double(x[i]);
        }
        return d;
    }

    template <typename Sample>
    std::vector<Sample> from_doubles(const std::vector<double>& x) {
        std::vector<Sample> s(x.size());
        for (std::size_t i = 0; i < x.size(); ++i) {
            s[i] = profile<Sample>::from_double(x[i]);
        }
        return s;
    }

    // ------------------------------------------------------------------------
    // Double-double arithmetic: an unevaluated sum hi + lo with |lo| <= ulp(hi)/2,
    // ~106 bits of significand. Only what the DFT below needs.
    // ------------------------------------------------------------------------
    struct dd {
        double hi;
        double lo;
    };

    /// Knuth TwoSum: s + e == a + b exactly, no ordering assumption.
    inline dd two_sum(double a, double b) {
        const double s  = a + b;
        const double bb = s - a;
        const double e  = (a - (s - bb)) + (b - bb);
        return {s, e};
    }

    /// Dekker FastTwoSum, valid when |a| >= |b|.
    inline dd quick_two_sum(double a, double b) {
        const double s = a + b;
        const double e = b - (s - a);
        return {s, e};
    }

    /// TwoProd via a correctly rounded fma: p + e == a * b exactly. std::fma
    /// is required to round once, whatever the hardware.
    inline dd two_prod(double a, double b) {
        const double p = a * b;
        const double e = std::fma(a, b, -p);
        return {p, e};
    }

    inline dd dd_add(dd a, dd b) {
        dd s = two_sum(a.hi, b.hi);
        dd t = two_sum(a.lo, b.lo);
        s.lo += t.hi;
        s = quick_two_sum(s.hi, s.lo);
        s.lo += t.lo;
        return quick_two_sum(s.hi, s.lo);
    }

    inline dd dd_neg(dd a) {
        return {-a.hi, -a.lo};
    }

    inline dd dd_mul(dd a, dd b) {
        dd p = two_prod(a.hi, b.hi);
        p.lo += a.hi * b.lo + a.lo * b.hi;
        return quick_two_sum(p.hi, p.lo);
    }

    inline dd dd_mul_d(dd a, double b) {
        dd p = two_prod(a.hi, b);
        p.lo += a.lo * b;
        return quick_two_sum(p.hi, p.lo);
    }

    inline dd dd_div_d(dd a, double b) {
        const double q1 = a.hi / b;
        const dd     r  = dd_add(a, dd_neg(dd_mul_d({b, 0.0}, q1)));
        const double q2 = r.hi / b;
        return quick_two_sum(q1, q2);
    }

    inline double dd_round(dd a) {
        return a.hi + a.lo;
    }

    // pi as a double-double: the standard 3.141592653589793 + 1.2246467991473532e-16.
    constexpr dd k_pi_dd{3.141592653589793, 1.2246467991473532e-16};

    /// sin and cos of |theta| <= pi/4 by Taylor series in double-double. Terms
    /// fall below 1e-34 relative by k = 15 for |theta| <= pi/4; 20 terms are
    /// summed for margin, largest first, which double-double handles exactly
    /// enough for a 1e-30 result. The oracle needs ~1e-20, so this is ample.
    struct dd_sincos {
        dd sin;
        dd cos;
    };

    dd_sincos sincos_small(dd theta) {
        const dd theta2 = dd_mul(theta, theta);
        dd       s_term = theta;
        dd       c_term = {1.0, 0.0};
        dd       s      = s_term;
        dd       c      = c_term;
        for (int k = 1; k <= 20; ++k) {
            // sin: term_k = -term_{k-1} * theta^2 / ((2k)(2k+1))
            // cos: term_k = -term_{k-1} * theta^2 / ((2k-1)(2k))
            s_term = dd_neg(dd_div_d(dd_mul(s_term, theta2), static_cast<double>((2 * k) * (2 * k + 1))));
            c_term = dd_neg(dd_div_d(dd_mul(c_term, theta2), static_cast<double>((2 * k - 1) * (2 * k))));
            s      = dd_add(s, s_term);
            c      = dd_add(c, c_term);
        }
        return {s, c};
    }

    /// (cos, sin)(2*pi*m/n) for integer m in [0, n) and n a power of two >= 4,
    /// by exact octant reduction on the integers (m and n are exact, m/n is
    /// exact because n is a power of two) followed by the small-angle series.
    dd_sincos twiddle(std::size_t m, std::size_t n) {
        const std::size_t half    = n / 2;
        const std::size_t quarter = n / 4;
        bool              negate  = false;
        bool              rotate  = false; // (cos, sin) -> (-sin, cos)
        m %= n;
        if (m >= half) {
            m -= half;
            negate = true;
        }
        if (m >= quarter) {
            m -= quarter;
            rotate = true;
        }
        // m in [0, quarter): theta in [0, pi/2). Fold to [0, pi/4] via
        // cos(theta) = sin(pi/2 - theta).
        bool swap = false;
        if (2 * m > quarter) {
            m    = quarter - m;
            swap = true;
        }
        const double    fraction = static_cast<double>(m) / static_cast<double>(n); // exact
        const dd        theta    = dd_mul_d({2.0 * k_pi_dd.hi, 2.0 * k_pi_dd.lo}, fraction);
        const dd_sincos sc       = sincos_small(theta);
        dd              c        = swap ? sc.sin : sc.cos;
        dd              s        = swap ? sc.cos : sc.sin;
        if (rotate) {
            const dd t = c;
            c          = dd_neg(s);
            s          = t;
        }
        if (negate) {
            c = dd_neg(c);
            s = dd_neg(s);
        }
        return {s, c};
    }

    // ------------------------------------------------------------------------
    // The compensated DFT, in Ooura's packing and sign convention (fft.h):
    //   forward: a[0] = Re X[0], a[1] = Re X[N/2], a[2k] = Re X[k], a[2k+1] = Im X[k],
    //            X[k] = sum_j x[j] * exp(+2*pi*i*j*k/N)
    //   inverse (unnormalized, gain N/2 on a round trip):
    //            x[k] = (R[0] + R[N/2] * (-1)^k) / 2
    //                 + sum_{j=1}^{N/2-1} (R[j] cos(2*pi*j*k/N) + I[j] sin(2*pi*j*k/N))
    // The inverse formula is fftsg.c's own statement of what rdft(n, -1, ...)
    // computes; both are restated from the definition, not from the code path.
    // ------------------------------------------------------------------------
    class compensated_dft {
      public:
        explicit compensated_dft(std::size_t n)
            : m_n(n)
            , m_twiddles(n) {
            for (std::size_t m = 0; m < n; ++m) {
                m_twiddles[m] = twiddle(m, n);
            }
        }

        std::vector<double> forward(const std::vector<double>& x) const {
            std::vector<double> a(m_n, 0.0);
            for (std::size_t k = 0; k <= m_n / 2; ++k) {
                dd re{0.0, 0.0};
                dd im{0.0, 0.0};
                for (std::size_t j = 0; j < m_n; ++j) {
                    const dd_sincos& w = m_twiddles[(j * k) % m_n];
                    re                 = dd_add(re, dd_mul_d(w.cos, x[j]));
                    im                 = dd_add(im, dd_mul_d(w.sin, x[j]));
                }
                if (k == 0) {
                    a[0] = dd_round(re);
                }
                else if (k == m_n / 2) {
                    a[1] = dd_round(re);
                }
                else {
                    a[2 * k]     = dd_round(re);
                    a[2 * k + 1] = dd_round(im);
                }
            }
            return a;
        }

        std::vector<double> inverse_unnormalized(const std::vector<double>& a) const {
            std::vector<double> x(m_n, 0.0);
            for (std::size_t k = 0; k < m_n; ++k) {
                // (R[0] +- R[N/2]) / 2: the sum error-free via TwoSum, the
                // halving exact.
                dd acc = dd_mul_d(two_sum(a[0], (k % 2 == 0) ? a[1] : -a[1]), 0.5);
                for (std::size_t j = 1; j < m_n / 2; ++j) {
                    const dd_sincos& w = m_twiddles[(j * k) % m_n];
                    acc                = dd_add(acc, dd_mul_d(w.cos, a[2 * j]));
                    acc                = dd_add(acc, dd_mul_d(w.sin, a[2 * j + 1]));
                }
                x[k] = dd_round(acc);
            }
            return x;
        }

        const dd_sincos& twiddle_at(std::size_t m) const { return m_twiddles[m % m_n]; }

      private:
        std::size_t            m_n;
        std::vector<dd_sincos> m_twiddles;
    };

    // ------------------------------------------------------------------------
    // Shared checkers.
    // ------------------------------------------------------------------------
    template <typename Sample>
    void expect_close(const std::vector<Sample>& got, const std::vector<double>& expected, double tol, const char* what,
                      std::size_t n) {
        ASSERT_EQ(got.size(), expected.size());
        for (std::size_t i = 0; i < got.size(); ++i) {
            const double g = profile<Sample>::to_double(got[i]);
            ASSERT_NEAR(g, expected[i], tol) << what << " N=" << n << " index " << i;
        }
    }

    /// Powers of two from the profile's minimum up to min(limit,
    /// TAP_DSP_PARITY_MAX_N, the profile's maximum). The cap is the same knob
    /// the parity gate uses (tests/CMakeLists.txt); the QEMU legs set 4096
    /// because a 65536-point double sweep needs several 512 KB buffers that
    /// the MPS2 data region does not have.
    template <typename Sample>
    std::vector<std::size_t> sizes_up_to(std::size_t limit) {
        const std::size_t top =
            std::min({limit, static_cast<std::size_t>(TAP_DSP_PARITY_MAX_N), profile<Sample>::k_max_n});
        std::vector<std::size_t> s;
        for (std::size_t n = profile<Sample>::k_min_n; n <= top; n *= 2) {
            s.push_back(n);
        }
        return s;
    }

    /// The closed-form sweep: up to 65536.
    template <typename Sample>
    std::vector<std::size_t> closed_form_sizes() {
        return sizes_up_to<Sample>(65536);
    }

    /// The compensated-DFT sizes: up to 256 (see "What the bound cannot see").
    template <typename Sample>
    std::vector<std::size_t> dft_sizes() {
        return sizes_up_to<Sample>(256);
    }

    /// Bin used for the on-bin materials: n/8, or 1 below n = 8. Always in
    /// [1, n/2), i.e. a genuinely complex bin, never DC or Nyquist.
    std::size_t tone_bin(std::size_t n) {
        return std::max<std::size_t>(1, n / 8);
    }

    // Runs the engine forward on x (given in the double domain, cast to the
    // profile) and checks it against the exact spectrum `expected`, with the
    // profile's forward scale applied.
    template <typename Sample>
    void check_forward(const std::vector<double>& x, const std::vector<double>& expected, const char* what) {
        const std::size_t                n = x.size();
        tap::dsp::basic_real_fft<Sample> fft(n);
        std::vector<Sample>              buf = from_doubles<Sample>(x);
        fft.forward_inplace(buf.data());
        const std::vector<double> scaled = scale_by(expected, profile<Sample>::forward_scale(n));
        expect_close(buf, scaled, profile<Sample>::tolerance(n, spectrum_norm2(x)), what, n);
    }

    // Runs the engine's UNNORMALIZED inverse on the packed spectrum `a` and
    // checks it against `expected` (which already carries the N/2 gain), with
    // the profile's inverse scale applied.
    template <typename Sample>
    void check_inverse(const std::vector<double>& a, const std::vector<double>& expected, const char* what) {
        const std::size_t                n = a.size();
        tap::dsp::basic_real_fft<Sample> fft(n);
        std::vector<Sample>              buf = from_doubles<Sample>(a);
        fft.inverse_inplace(buf.data());
        const std::vector<double> scaled = scale_by(expected, profile<Sample>::inverse_scale(n));
        // ||expected||_2 plays the role of ||y||_2 for the inverse direction;
        // the packed spectrum's 2-norm is within sqrt(2) of the true one, and
        // the N/2 gain means ||expected||_2 >= ||a||_2 * sqrt(N)/2 for these
        // vectors, so the bound stays of the same form.
        double e = 0.0;
        for (const double v : expected) {
            e += v * v;
        }
        expect_close(buf, scaled, profile<Sample>::tolerance(n, std::sqrt(e)), what, n);
    }

    template <typename Sample>
    class fft_oracle_test : public ::testing::Test {};
    TYPED_TEST_SUITE(fft_oracle_test, oracle_types);

    // ========================================================================
    // Closed forms, forward. Every one is driven at the profile's full-scale
    // amplitude (1.0 for float and double); the exact answers scale with it.
    // ========================================================================

    TYPED_TEST(fft_oracle_test, ImpulseIsFlatAtEverySize) {
        const double a = profile<TypeParam>::k_full_scale;
        for (const std::size_t n : closed_form_sizes<TypeParam>()) {
            std::vector<double> x(n, 0.0);
            x[0] = a;
            std::vector<double> expected(n, 0.0);
            expected[0] = a; // DC
            expected[1] = a; // Nyquist
            for (std::size_t k = 1; k < n / 2; ++k) {
                expected[2 * k] = a;
            }
            check_forward<TypeParam>(x, expected, "impulse");
        }
    }

    TYPED_TEST(fft_oracle_test, DcLandsAsNInSlotZero) {
        const double a = profile<TypeParam>::k_full_scale;
        for (const std::size_t n : closed_form_sizes<TypeParam>()) {
            std::vector<double> x(n, a);
            std::vector<double> expected(n, 0.0);
            expected[0] = a * static_cast<double>(n);
            check_forward<TypeParam>(x, expected, "dc");
        }
    }

    TYPED_TEST(fft_oracle_test, NyquistAlternationLandsAsNInSlotOne) {
        const double a = profile<TypeParam>::k_full_scale;
        for (const std::size_t n : closed_form_sizes<TypeParam>()) {
            std::vector<double> x(n);
            for (std::size_t j = 0; j < n; ++j) {
                x[j] = (j % 2 == 0) ? a : -a;
            }
            std::vector<double> expected(n, 0.0);
            expected[1] = a * static_cast<double>(n);
            check_forward<TypeParam>(x, expected, "nyquist");
        }
    }

    TYPED_TEST(fft_oracle_test, OnBinCosineIsPlusHalfNInTheRealSlot) {
        const double a = profile<TypeParam>::k_full_scale;
        for (const std::size_t n : closed_form_sizes<TypeParam>()) {
            const std::size_t   k = tone_bin(n);
            const auto          x = tap::dsp::test::tone<double>(n, static_cast<double>(k), a, 0.0);
            std::vector<double> expected(n, 0.0);
            expected[2 * k] = a * static_cast<double>(n) / 2.0;
            check_forward<TypeParam>(x, expected, "cosine");
        }
    }

    // The sign convention: a sine at bin k lands at +N/2 in the imaginary slot
    // (it would be -N/2 in the engineering convention exp(-2*pi*i/N)).
    TYPED_TEST(fft_oracle_test, OnBinSineIsPlusHalfNInTheImaginarySlot) {
        const double a = profile<TypeParam>::k_full_scale;
        for (const std::size_t n : closed_form_sizes<TypeParam>()) {
            const std::size_t k = tone_bin(n);
            // sin(w j) = cos(w j - pi/2)
            const auto          x = tap::dsp::test::tone<double>(n, static_cast<double>(k), a, -std::numbers::pi / 2.0);
            std::vector<double> expected(n, 0.0);
            expected[2 * k + 1] = a * static_cast<double>(n) / 2.0;
            check_forward<TypeParam>(x, expected, "sine");
        }
    }

    TYPED_TEST(fft_oracle_test, TwoToneSuperposesLinearly) {
        const double a = profile<TypeParam>::k_full_scale;
        for (const std::size_t n : closed_form_sizes<TypeParam>()) {
            if (n < 8) {
                continue; // needs two distinct complex bins
            }
            const std::size_t k1 = tone_bin(n);
            const std::size_t k2 = k1 + 1;
            const double      a1 = 0.75 * a;
            const double      a2 = -0.25 * a;
            const auto        c  = tap::dsp::test::tone<double>(n, static_cast<double>(k1), a1, 0.0);
            const auto        s = tap::dsp::test::tone<double>(n, static_cast<double>(k2), a2, -std::numbers::pi / 2.0);
            std::vector<double> x(n);
            for (std::size_t j = 0; j < n; ++j) {
                x[j] = c[j] + s[j];
            }
            std::vector<double> expected(n, 0.0);
            expected[2 * k1]     = a1 * static_cast<double>(n) / 2.0;
            expected[2 * k2 + 1] = a2 * static_cast<double>(n) / 2.0;
            check_forward<TypeParam>(x, expected, "two-tone");
        }
    }

    // ========================================================================
    // Closed forms, inverse (unnormalized: gain N/2).
    // ========================================================================

    TYPED_TEST(fft_oracle_test, InverseOfFlatSpectrumIsHalfNImpulse) {
        const double amp = profile<TypeParam>::k_full_scale;
        for (const std::size_t n : closed_form_sizes<TypeParam>()) {
            std::vector<double> a(n, 0.0);
            a[0] = amp;
            a[1] = amp;
            for (std::size_t k = 1; k < n / 2; ++k) {
                a[2 * k] = amp;
            }
            std::vector<double> expected(n, 0.0);
            expected[0] = amp * static_cast<double>(n) / 2.0;
            check_inverse<TypeParam>(a, expected, "inverse impulse");
        }
    }

    TYPED_TEST(fft_oracle_test, InverseOfDcSpectrumIsHalfNConstant) {
        const double amp = profile<TypeParam>::k_full_scale;
        for (const std::size_t n : closed_form_sizes<TypeParam>()) {
            std::vector<double> a(n, 0.0);
            a[0] = amp;
            std::vector<double> expected(n, 0.5 * amp);
            check_inverse<TypeParam>(a, expected, "inverse dc");
        }
    }

    TYPED_TEST(fft_oracle_test, InverseOfNyquistSpectrumIsHalfNAlternation) {
        const double amp = profile<TypeParam>::k_full_scale;
        for (const std::size_t n : closed_form_sizes<TypeParam>()) {
            std::vector<double> a(n, 0.0);
            a[1] = amp;
            std::vector<double> expected(n);
            for (std::size_t k = 0; k < n; ++k) {
                expected[k] = (k % 2 == 0) ? 0.5 * amp : -0.5 * amp;
            }
            check_inverse<TypeParam>(a, expected, "inverse nyquist");
        }
    }

    TYPED_TEST(fft_oracle_test, InverseOfOnBinSpectrumIsHalfNTone) {
        const double amp = profile<TypeParam>::k_full_scale;
        for (const std::size_t n : closed_form_sizes<TypeParam>()) {
            const std::size_t   k = tone_bin(n);
            std::vector<double> a(n, 0.0);
            a[2 * k]     = amp; // cosine part
            a[2 * k + 1] = amp; // sine part, +i convention
            // x[j] = amp * (cos(w j) + sin(w j)): the unnormalized inverse has
            // gain N/2 on the round trip, so a bin of amp comes back with
            // amplitude amp.
            const auto c = tap::dsp::test::tone<double>(n, static_cast<double>(k), amp, 0.0);
            const auto s = tap::dsp::test::tone<double>(n, static_cast<double>(k), amp, -std::numbers::pi / 2.0);
            std::vector<double> expected(n);
            for (std::size_t j = 0; j < n; ++j) {
                expected[j] = c[j] + s[j];
            }
            check_inverse<TypeParam>(a, expected, "inverse tone");
        }
    }

    // ========================================================================
    // The compensated DFT, N <= 256.
    // ========================================================================

    // The oracle checks itself first: its twiddles agree with libm to double
    // rounding, and its own forward followed by its own inverse is the
    // identity to double rounding, so a disagreement below is the engine's.
    // The two constants are derived bounds, not fitted numbers. Twiddles: the
    // comparison forms theta = 2*pi*m/n in double, and rounding theta (at most
    // half an ulp of 2*pi, 4.4e-16 = 2 eps) moves cos/sin by up to that much;
    // libm adds at most 1 ulp (1 eps for values in [0.5, 1]); dd_round adds
    // half an ulp. Bound 3.5 eps, constant 4 eps; measured maximum 6.4e-16 =
    // 2.9 eps (x86-64 Linux, glibc 2.39), i.e. the dd twiddles are as good as
    // libm and the argument rounding dominates. Round trip: forward then
    // inverse re-rounds through two dd_round steps plus the final 2/N
    // multiply (exact, power of two), so ~1.5 eps relative on unit-scale data;
    // constant 8 eps; measured maximum 1.1e-16 = 0.5 eps.
    TEST(fft_oracle_self_check, TwiddlesMatchLibmToDoubleRounding) {
        for (const std::size_t n : dft_sizes<double>()) {
            const compensated_dft oracle(n);
            for (std::size_t m = 0; m < n; ++m) {
                const double theta = 2.0 * std::numbers::pi * static_cast<double>(m) / static_cast<double>(n);
                const auto&  w     = oracle.twiddle_at(m);
                EXPECT_NEAR(dd_round(w.cos), std::cos(theta), 4.0 * std::numeric_limits<double>::epsilon())
                    << "N=" << n << " m=" << m;
                EXPECT_NEAR(dd_round(w.sin), std::sin(theta), 4.0 * std::numeric_limits<double>::epsilon())
                    << "N=" << n << " m=" << m;
                // The low word is a genuine correction, not noise: it is at
                // most half an ulp of the high word.
                EXPECT_LE(std::fabs(w.cos.lo), std::fabs(w.cos.hi) * std::numeric_limits<double>::epsilon() + 1e-300);
            }
        }
    }

    TEST(fft_oracle_self_check, RoundTripIsIdentityToDoubleRounding) {
        for (const std::size_t n : dft_sizes<double>()) {
            const compensated_dft     oracle(n);
            const std::vector<double> x    = tap::dsp::test::random_signal<double>(n, 0xC0FFEEu);
            const std::vector<double> a    = oracle.forward(x);
            const std::vector<double> back = oracle.inverse_unnormalized(a);
            const double              gain = 2.0 / static_cast<double>(n);
            for (std::size_t i = 0; i < n; ++i) {
                EXPECT_NEAR(back[i] * gain, x[i], 8.0 * std::numeric_limits<double>::epsilon())
                    << "N=" << n << " " << i;
            }
        }
    }

    template <typename Sample>
    void expect_forward_matches_dft(std::size_t n, const std::vector<double>& x, const char* what) {
        const compensated_dft     oracle(n);
        const std::vector<double> expected = oracle.forward(x);
        check_forward<Sample>(x, expected, what);
    }

    template <typename Sample>
    void expect_inverse_matches_dft(std::size_t n, const std::vector<double>& a, const char* what) {
        const compensated_dft     oracle(n);
        const std::vector<double> expected = oracle.inverse_unnormalized(a);
        check_inverse<Sample>(a, expected, what);
    }

    TYPED_TEST(fft_oracle_test, ForwardMatchesCompensatedDftOnBroadband) {
        for (const std::size_t n : dft_sizes<TypeParam>()) {
            // Drawn in the profile first so the oracle sees exactly the bits
            // the engine sees; from_doubles is then the identity.
            const auto x = to_doubles(tap::dsp::test::random_signal<TypeParam>(n, 0x2545F491u));
            expect_forward_matches_dft<TypeParam>(n, x, "dft forward broadband");
        }
    }

    TYPED_TEST(fft_oracle_test, ForwardMatchesCompensatedDftOnOffBinTone) {
        for (const std::size_t n : dft_sizes<TypeParam>()) {
            // Off-bin (k + 0.37) so every bin carries leakage and none is
            // numerically empty: the opposite regime from the closed forms.
            const auto x =
                to_doubles(tap::dsp::test::tone<TypeParam>(n, static_cast<double>(tone_bin(n)) + 0.37, 0.8, 1.1));
            expect_forward_matches_dft<TypeParam>(n, x, "dft forward off-bin tone");
        }
    }

    TYPED_TEST(fft_oracle_test, InverseMatchesCompensatedDftOnBroadbandSpectrum) {
        for (const std::size_t n : dft_sizes<TypeParam>()) {
            // An arbitrary packed spectrum, not one produced by a forward, so
            // the inverse is judged on its own definition.
            const auto a = to_doubles(tap::dsp::test::random_signal<TypeParam>(n, 0x1D872B41u));
            expect_inverse_matches_dft<TypeParam>(n, a, "dft inverse broadband");
        }
    }

    TYPED_TEST(fft_oracle_test, InverseMatchesCompensatedDftOnToneSpectrum) {
        for (const std::size_t n : dft_sizes<TypeParam>()) {
            const compensated_dft oracle(n);
            const auto            x =
                to_doubles(tap::dsp::test::tone<TypeParam>(n, static_cast<double>(tone_bin(n)) + 0.37, 0.8, 1.1));
            // The spectrum is rounded to the profile before both sides see it.
            const auto a = to_doubles(from_doubles<TypeParam>(oracle.forward(x)));
            expect_inverse_matches_dft<TypeParam>(n, a, "dft inverse tone spectrum");
        }
    }

} // namespace
