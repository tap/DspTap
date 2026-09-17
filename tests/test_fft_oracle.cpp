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
//      from the sample type's epsilon and log2 N (see tolerance() below).
//
//   2. A compensated-summation DFT for N <= 256: the O(N^2) definition with
//      double-double accumulation (TwoSum / TwoProd, Dekker 1971 and Knuth
//      TAOCP vol. 2 4.2.2; Shewchuk 1997 for the error-free transformations)
//      and double-double twiddles from a Taylor series after exact octant
//      reduction on the integer bin*sample index. Not `long double`: that is
//      `double` on MSVC and on Apple arm64 (Part 6, item N19), so it would be
//      no oracle at all on two of the three hosted CI legs.
//
// The suite is typed over float and double now. The fixed-point stage
// (Stage 3, Part 7) extends `profile<Sample>` below with int16_t and int32_t,
// supplying the profile's documented forward scale and the conversion into
// the double domain; the tests themselves are written against that trait so
// nothing else changes.

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

namespace {

    // ------------------------------------------------------------------------
    // Per-profile traits — THE EXTENSION POINT FOR THE FIXED-POINT STAGE.
    //
    // A profile says how to get a sample into the double domain, what its
    // unit roundoff is, and what scale the engine's forward / unnormalized
    // inverse carry relative to the mathematical DFT (1 and 1 for the float
    // profiles; the Q15/Q31 profiles document X/N under `fixed` scaling and
    // will state their own inverse scale). Stage 3 adds:
    //     template <> struct profile<std::int16_t> { ... };
    //     template <> struct profile<std::int32_t> { ... };
    // and appends the two types to oracle_types. The tolerance for those
    // profiles is a quantization-noise number from Part 7, not epsilon-based;
    // tolerance() dispatches on the trait so that is a local change too.
    // ------------------------------------------------------------------------
    template <typename Sample>
    struct profile;

    template <>
    struct profile<double> {
        static constexpr double k_epsilon       = std::numeric_limits<double>::epsilon();
        static constexpr double k_forward_scale = 1.0; ///< engine forward = k * DFT
        static constexpr double k_inverse_scale = 1.0; ///< engine inverse = k * Ooura's unnormalized inverse
        static double           to_double(double v) { return v; }
        static double           from_double(double v) { return v; }
    };

    template <>
    struct profile<float> {
        static constexpr double k_epsilon       = std::numeric_limits<float>::epsilon();
        static constexpr double k_forward_scale = 1.0;
        static constexpr double k_inverse_scale = 1.0;
        static double           to_double(float v) { return static_cast<double>(v); }
        static float            from_double(double v) { return static_cast<float>(v); }
    };

    using oracle_types = ::testing::Types<float, double>;

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
    // where u is the unit roundoff (epsilon / 2). Ooura's tables are computed
    // from libm and stored in Sample, so mu <= u, hence eta <= 5u + O(u^2)
    // and, for t*eta << 1, the 2-norm error is below 5 * u * log2(N) * ||y||_2.
    // The largest single-element error cannot exceed the 2-norm, and a
    // measured (not remembered) value for || y ||_2 comes free with the exact
    // answer. Two further terms ride on top: casting the closed-form input to
    // Sample perturbs each x_j by at most u|x_j|, which by Parseval moves y by
    // at most u * ||y||_2 in 2-norm; and Ooura is split-radix, not the radix-2
    // of the theorem, so its constant differs in a factor of order one. The
    // bound is therefore taken as k_higham_constant * epsilon * log2(N) *
    // ||y||_2 with k_higham_constant = 4 (epsilon = 2u, so 4 * epsilon = 8u:
    // the 5u of the theorem, the u of the input cast, and 2u of margin for
    // the split-radix constant).
    //
    // Measured against that bound on x86-64 Linux (GCC 13, -O3, glibc libm,
    // Ooura as the engine): the largest ratio |error| / (epsilon * log2(N) *
    // ||y||_2) over every test in this file was 0.25 for double and 0.17 for
    // float, i.e. the engine sits 16x (double) and 23x (float) inside the
    // derived bound of 4. That margin is deliberate: this is an oracle for
    // correctness, not a precision ratchet; the measured-number pins live in
    // test_fft.cpp and test_fft_backend.cpp.
    // ------------------------------------------------------------------------
    constexpr double k_higham_constant = 4.0;

    template <typename Sample>
    double tolerance(std::size_t n, double spectrum_norm2) {
        const double log2n = std::log2(static_cast<double>(n));
        return k_higham_constant * profile<Sample>::k_epsilon * log2n * spectrum_norm2;
    }

    /// ||y||_2 of the full complex spectrum of x, from Parseval: sqrt(N) * ||x||_2.
    double spectrum_norm2(const std::vector<double>& x) {
        double e = 0.0;
        for (const double v : x) {
            e += v * v;
        }
        return std::sqrt(static_cast<double>(x.size()) * e);
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
                const double edge = (k % 2 == 0) ? a[0] + a[1] : a[0] - a[1];
                dd           acc  = dd_mul_d({edge, 0.0}, 0.5);
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

    std::vector<std::size_t> closed_form_sizes() {
        std::vector<std::size_t> s;
        for (std::size_t n = 4; n <= 65536; n *= 2) {
            s.push_back(n);
        }
        return s;
    }

    std::vector<std::size_t> dft_sizes() {
        std::vector<std::size_t> s;
        for (std::size_t n = 4; n <= 256; n *= 2) {
            s.push_back(n);
        }
        return s;
    }

    /// Bin used for the on-bin materials: n/8, or 1 below n = 8. Always in
    /// [1, n/2), i.e. a genuinely complex bin, never DC or Nyquist.
    std::size_t tone_bin(std::size_t n) {
        return std::max<std::size_t>(1, n / 8);
    }

    // Runs the engine forward on x (given in the double domain, cast to the
    // profile) and checks it against the exact spectrum `expected`.
    template <typename Sample>
    void check_forward(const std::vector<double>& x, const std::vector<double>& expected, const char* what) {
        const std::size_t                n = x.size();
        tap::dsp::basic_real_fft<Sample> fft(n);
        std::vector<Sample>              buf = from_doubles<Sample>(x);
        fft.forward_inplace(buf.data());
        std::vector<double> scaled = expected;
        for (double& v : scaled) {
            v *= profile<Sample>::k_forward_scale;
        }
        expect_close(buf, scaled, tolerance<Sample>(n, spectrum_norm2(x)), what, n);
    }

    // Runs the engine's UNNORMALIZED inverse on the packed spectrum `a` and
    // checks it against `expected` (which already carries the N/2 gain).
    template <typename Sample>
    void check_inverse(const std::vector<double>& a, const std::vector<double>& expected, const char* what) {
        const std::size_t                n = a.size();
        tap::dsp::basic_real_fft<Sample> fft(n);
        std::vector<Sample>              buf = from_doubles<Sample>(a);
        fft.inverse_inplace(buf.data());
        std::vector<double> scaled = expected;
        for (double& v : scaled) {
            v *= profile<Sample>::k_inverse_scale;
        }
        // ||expected||_2 plays the role of ||y||_2 for the inverse direction;
        // the packed spectrum's 2-norm is within sqrt(2) of the true one, and
        // the N/2 gain means ||expected||_2 >= ||a||_2 * sqrt(N)/2 for these
        // vectors, so the bound stays of the same form.
        double e = 0.0;
        for (const double v : expected) {
            e += v * v;
        }
        expect_close(buf, scaled, tolerance<Sample>(n, std::sqrt(e)), what, n);
    }

    template <typename Sample>
    class fft_oracle_test : public ::testing::Test {};
    TYPED_TEST_SUITE(fft_oracle_test, oracle_types);

    // ========================================================================
    // Closed forms, forward.
    // ========================================================================

    TYPED_TEST(fft_oracle_test, ImpulseHasFlatSpectrum) {
        for (const std::size_t n : closed_form_sizes()) {
            std::vector<double> x(n, 0.0);
            x[0] = 1.0;
            std::vector<double> expected(n, 0.0);
            expected[0] = 1.0; // DC
            expected[1] = 1.0; // Nyquist
            for (std::size_t k = 1; k < n / 2; ++k) {
                expected[2 * k] = 1.0;
            }
            check_forward<TypeParam>(x, expected, "impulse");
        }
    }

    TYPED_TEST(fft_oracle_test, DcLandsAsNInSlotZero) {
        for (const std::size_t n : closed_form_sizes()) {
            std::vector<double> x(n, 1.0);
            std::vector<double> expected(n, 0.0);
            expected[0] = static_cast<double>(n);
            check_forward<TypeParam>(x, expected, "dc");
        }
    }

    TYPED_TEST(fft_oracle_test, NyquistAlternationLandsAsNInSlotOne) {
        for (const std::size_t n : closed_form_sizes()) {
            std::vector<double> x(n);
            for (std::size_t j = 0; j < n; ++j) {
                x[j] = (j % 2 == 0) ? 1.0 : -1.0;
            }
            std::vector<double> expected(n, 0.0);
            expected[1] = static_cast<double>(n);
            check_forward<TypeParam>(x, expected, "nyquist");
        }
    }

    TYPED_TEST(fft_oracle_test, OnBinCosineIsPlusHalfNInTheRealSlot) {
        for (const std::size_t n : closed_form_sizes()) {
            const std::size_t   k = tone_bin(n);
            const auto          x = tap::dsp::test::tone<double>(n, static_cast<double>(k), 1.0, 0.0);
            std::vector<double> expected(n, 0.0);
            expected[2 * k] = static_cast<double>(n) / 2.0;
            check_forward<TypeParam>(x, expected, "cosine");
        }
    }

    // The sign convention: a sine at bin k lands at +N/2 in the imaginary slot
    // (it would be -N/2 in the engineering convention exp(-2*pi*i/N)).
    TYPED_TEST(fft_oracle_test, OnBinSineIsPlusHalfNInTheImaginarySlot) {
        for (const std::size_t n : closed_form_sizes()) {
            const std::size_t k = tone_bin(n);
            // sin(w j) = cos(w j - pi/2)
            const auto x = tap::dsp::test::tone<double>(n, static_cast<double>(k), 1.0, -std::numbers::pi / 2.0);
            std::vector<double> expected(n, 0.0);
            expected[2 * k + 1] = static_cast<double>(n) / 2.0;
            check_forward<TypeParam>(x, expected, "sine");
        }
    }

    TYPED_TEST(fft_oracle_test, TwoToneSuperposesLinearly) {
        for (const std::size_t n : closed_form_sizes()) {
            if (n < 8) {
                continue; // needs two distinct complex bins
            }
            const std::size_t k1 = tone_bin(n);
            const std::size_t k2 = k1 + 1;
            const double      a1 = 0.75;
            const double      a2 = -0.25;
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
        for (const std::size_t n : closed_form_sizes()) {
            std::vector<double> a(n, 0.0);
            a[0] = 1.0;
            a[1] = 1.0;
            for (std::size_t k = 1; k < n / 2; ++k) {
                a[2 * k] = 1.0;
            }
            std::vector<double> expected(n, 0.0);
            expected[0] = static_cast<double>(n) / 2.0;
            check_inverse<TypeParam>(a, expected, "inverse impulse");
        }
    }

    TYPED_TEST(fft_oracle_test, InverseOfDcSpectrumIsHalfNConstant) {
        for (const std::size_t n : closed_form_sizes()) {
            std::vector<double> a(n, 0.0);
            a[0] = 1.0;
            std::vector<double> expected(n, 0.5);
            check_inverse<TypeParam>(a, expected, "inverse dc");
        }
    }

    TYPED_TEST(fft_oracle_test, InverseOfNyquistSpectrumIsHalfNAlternation) {
        for (const std::size_t n : closed_form_sizes()) {
            std::vector<double> a(n, 0.0);
            a[1] = 1.0;
            std::vector<double> expected(n);
            for (std::size_t k = 0; k < n; ++k) {
                expected[k] = (k % 2 == 0) ? 0.5 : -0.5;
            }
            check_inverse<TypeParam>(a, expected, "inverse nyquist");
        }
    }

    TYPED_TEST(fft_oracle_test, InverseOfOnBinSpectrumIsHalfNTone) {
        for (const std::size_t n : closed_form_sizes()) {
            const std::size_t   k = tone_bin(n);
            std::vector<double> a(n, 0.0);
            a[2 * k]     = 1.0; // cosine part
            a[2 * k + 1] = 1.0; // sine part, +i convention
            // x[j] = cos(w j) + sin(w j), unnormalized inverse gain is N/2 on
            // the round trip, i.e. a unit bin comes back with amplitude 1.
            const auto c = tap::dsp::test::tone<double>(n, static_cast<double>(k), 1.0, 0.0);
            const auto s = tap::dsp::test::tone<double>(n, static_cast<double>(k), 1.0, -std::numbers::pi / 2.0);
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
    // identity to ~1e-15 relative, so a disagreement below is the engine's.
    TEST(fft_oracle_self_check, TwiddlesMatchLibmToDoubleRounding) {
        for (const std::size_t n : dft_sizes()) {
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
        for (const std::size_t n : dft_sizes()) {
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
        for (const std::size_t n : dft_sizes()) {
            // Drawn in the profile first so the oracle sees exactly the bits
            // the engine sees; from_doubles is then the identity.
            const auto x = to_doubles(tap::dsp::test::random_signal<TypeParam>(n, 0x2545F491u));
            expect_forward_matches_dft<TypeParam>(n, x, "dft forward broadband");
        }
    }

    TYPED_TEST(fft_oracle_test, ForwardMatchesCompensatedDftOnOffBinTone) {
        for (const std::size_t n : dft_sizes()) {
            // Off-bin (k + 0.37) so every bin carries leakage and none is
            // numerically empty: the opposite regime from the closed forms.
            const auto x =
                to_doubles(tap::dsp::test::tone<TypeParam>(n, static_cast<double>(tone_bin(n)) + 0.37, 0.8, 1.1));
            expect_forward_matches_dft<TypeParam>(n, x, "dft forward off-bin tone");
        }
    }

    TYPED_TEST(fft_oracle_test, InverseMatchesCompensatedDftOnBroadbandSpectrum) {
        for (const std::size_t n : dft_sizes()) {
            // An arbitrary packed spectrum, not one produced by a forward, so
            // the inverse is judged on its own definition.
            const auto a = to_doubles(tap::dsp::test::random_signal<TypeParam>(n, 0x1D872B41u));
            expect_inverse_matches_dft<TypeParam>(n, a, "dft inverse broadband");
        }
    }

    TYPED_TEST(fft_oracle_test, InverseMatchesCompensatedDftOnToneSpectrum) {
        for (const std::size_t n : dft_sizes()) {
            const compensated_dft oracle(n);
            const auto            x =
                to_doubles(tap::dsp::test::tone<TypeParam>(n, static_cast<double>(tone_bin(n)) + 0.37, 0.8, 1.1));
            // The spectrum is rounded to the profile before both sides see it.
            const auto a = to_doubles(from_doubles<TypeParam>(oracle.forward(x)));
            expect_inverse_matches_dft<TypeParam>(n, a, "dft inverse tone spectrum");
        }
    }

} // namespace
