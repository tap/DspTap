/// @file iir_design.h
/// @brief Recursive design math: the two-path allpass half-band (constexpr), the odd-order Chebyshev II, the committed
/// tables.
// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// The design half of SampleRateTap's recursive polyphase stages
// (docs/IIR_PLAN.md, I1 and I4), landed here first per the substrate rule:
// rational holds only the stage types. Published literature only; NOTICE.md
// keeps the clean-room record and traces every formula below to its source.
//
// The structures (Valenzuela and Constantinides 1983; Vaidyanathan 1993, the
// power-symmetric IIR filter and its allpass decomposition):
//
//   H(z) = 1/2 [A0(z^2) + z^-1 A1(z^2)]
//
// a half-band low-pass of the higher rate F_hi of a by-2 stage. The nonlinear
// variant's paths are cascades of first-order allpass sections in w = z^2,
// (a + w^-1) / (1 + a w^-1), 0 < a < 1; the near-linear variant's A0 is a
// pure delay w^-K and A1 a real allpass of order N in w that approximates
// the delay w^-(K - 1/2) (Renfors and Saramaki 1987's structure). The
// four-path stage is H(z) = 1/4 sum_{k=0}^{3} z^-k A_k(z^4) (harris 2004's
// M-path recursive polyphase filter; Renfors and Saramaki 1987).
//
// A spec point is (A, p, d): the stopband attenuation A in dB, the passband
// edge p as a fraction of the chain's lowest rate r_min, and the design
// divisor d, the stage's lower rate over r_min (rational's PLAN.md 3.1;
// IIR_PLAN I13). A stage at divisor d has its passband edge at p / d of its
// lower rate, so f_pass = p / (2 d) of its higher rate, and its stopband
// mirrored about its own low Nyquist: from 1/2 - f_pass to 1/2 of F_hi.
// Every design here meets -(A + margin) dB there, margin 1 dB
// (k_iir_design_margin_db: rational's design margin and nyquist.h's search
// criterion), and the four spec points of the family's profiles, at the
// divisors 1, 2, 4 and 8 a power-of-two cascade uses, are tabulated:
//
//   | A      | p     | (rational's profile name, for orientation only) |
//   |--------|-------|-------------------------------------------------|
//   | 70 dB  | 1/3   | super_economy                                   |
//   | 70 dB  | 3/8   | economy                                         |
//   | 70 dB  | 19/48 | balanced                                        |
//   | 120 dB | 5/12  | transparent                                     |
//
// keyed by the numbers, never by those names (I1: the names are rational's).
//
// Contract, as numbers:
//   - The nonlinear half-band (design_nonlinear_halfband): the elliptic
//     (bireciprocal) half-band of odd order N, (N - 1) / 2 allpass
//     coefficients in (0, 1), ascending, alternating between the paths
//     (coefficients 0, 2, 4, ... in A0; 1, 3, 5, ... in A1), from Gazsi 1985's
//     explicit formulas in the bireciprocal case. N is the smallest odd order
//     whose closed-form attenuation, 10 log10(1 + q^(-N/2) / 4) with q the
//     nome, reaches A + margin (halfband_order). Power complementary:
//     |H(f)|^2 + |H(1/2 - f)|^2 = 1, so the passband deviates from 0 dB by
//     about 10^(-A/10) / 2 in power (4e-7 dB at 70 dB). constexpr, no libm
//     (detail/cx_math.h); the named spec points are compile-time constants
//     (k_nonlinear_halfband<A, PassNum, PassDen, D>), and a custom spec calls
//     the same function at run time. Accuracy: every tabulated coefficient
//     within 13.2 ulps of an mpmath leg at 60 digits (tools/iir/
//     gen_cx_reference.py; tests/test_iir_design.cpp pins 13.5), under the
//     15 ulps IIR-0's amendment promised; bits pinned per compiler in
//     tests/test_iir_constexpr.cpp (one row for every compiler CI runs).
//     Three steps are rewritten by algebra on the published formulas, each
//     removing a cancellation: the nome step's 1 - k'^(1/2) = k^2 / ((1 +
//     k'^(1/2)) (1 + k')) (IIR-0, SampleRateTap#104; without it the relaxed
//     designs lost up to 94,049 ulps); the last step a = (1 - x) / (1 + x)
//     as a = (1 - x^2) / (1 + x)^2 with 1 - x^2 = w^2 (1 + k)^2 / (k (1 +
//     w^2)^2); and 1 - w^2 / k as (t - w)(t + w) / k, t = tan(pi f_pass),
//     whose difference is exact where w nears t at the band edge (15.4 ulps
//     at worst without the last, 13.2 with it).
//   - Orders at the tabulated spec points (the closed form, reproducing
//     IIR-0's grid search): N9 / N5 / N5 / N3 at the three 70 dB points (d =
//     1 / 2 / 4 / 8), N17 / N9 / N7 / N5 at 120 dB, 5/12.
//   - The near-linear half-band: an iterative design (minimax over the
//     passband phase error, which equals the mirrored stopband), so its
//     tables are committed (iir_tables.h, tools/iir/gen_iir_tables.py) and
//     pinned by tests: order N, delay K (K = N at every tabulated point), the
//     sections of A1 (second-order (c1, c2) pairs ordered by increasing pole
//     radius, then one first-order c1 when N is odd), the attained stopband,
//     and the passband delay's deviation from 2K high-rate samples.
//   - The four-path stage: committed tables (the same generator) at the
//     placements SampleRateTap's IIR-0 found to pass I12's bar, design
//     divisor 2 at the three 70 dB spec points: six sections, laid out
//     [2, 2, 1, 1] over the paths.
//   - The guard band filter's Chebyshev II (design_cheby2_lowpass): odd
//     order n, flat within ripple_db (0.01 dB) to p of its rate, A dB down
//     from the stopband edge that follows, Os = Op cosh(acosh(D) / n): the
//     analogue prototype's poles the scaled reciprocals of Chebyshev I poles,
//     its zeros at j Os / cos(theta_k), the bilinear transform with
//     prewarping (Antoniou 2006),
//     (n - 1) / 2 biquads each pairing a pole pair with the zero pair of the
//     same theta, and one first-order section with its zero at z = -1;
//     sections ordered by increasing pole radius, each at unity DC gain.
//     constexpr: poles, zeros, pairing, order and the stop edge. The L-inf
//     gain distribution (Jackson 1970: each partial cascade peaks at 0 dB,
//     the last section restoring unity DC gain) is NOT constexpr: its sweep
//     takes about 415k constant-evaluation steps at order 11, four times
//     MSVC's default limit (IIR-0), so its factors are committed per spec
//     point and order (orders 3 to 11; tools/iir/gen_iir_tables.py) and
//     applied in constexpr (k_guard_band_cheby2<A, PassNum, PassDen, Order>);
//     cheby2_linf_gains() computes them at run time for a custom design.
//     Every committed factor is exactly 1.0: with the sections in this order,
//     each at unity DC gain, every partial cascade of every committed design
//     peaks at DC, so there the L-inf distribution is unity DC gain per
//     section (tests/test_iir_design.cpp pins it; a custom design need not
//     share it). Accuracy against the mpmath leg: 3.4 ulps at worst over
//     every committed section, the stop edges within 0.6.
//   - The linear guard band form: even-length symmetric (type II) FIR,
//     linear phase, a zero at Nyquist by its symmetry. Committed tables by
//     scipy.signal.remez at every (spec point, Chebyshev order) where it
//     converges, flat within 0.01 dB peak to peak to p and A dB down from the
//     Chebyshev's stop edge at that order (the same generator); a custom
//     setting designs at run time through kaiser.h
//     (design_guard_band_linear), pinned by tolerance. Q15 has rows of its
//     own (k_guard_band_fir_q15<A, PassNum, PassDen, Order>, since 0.4.0) at
//     the three 70 dB spec points and orders 3 to 11: designs that keep the
//     same spec from the same stop edge once quantized to Q1.14 by
//     quantize_symmetric_row_preserving_sum<std::int16_t> (quantize.h), the
//     same length as the double row or up to 12 taps longer; none at 120 dB,
//     which Q1.14 cannot carry (iir_tables.h states both).
//   - Responses and group delays (the *_response and *_group_delay
//     functions): runtime double over libm, for the stages' I5 latency
//     members and the tests.
//
// The constexpr designers are a handful of coefficients per design, not
// long FIR tables (kaiser.h's design note): each evaluation takes at most a
// few thousand of Clang's constant-evaluation steps (the largest tabulated,
// tests/test_iir_constexpr.cpp states them), far below MSVC's documented
// default of 100,000.

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numbers>
#include <span>
#include <vector>

#include "tap/dsp/detail/cx_math.h"
#include "tap/dsp/iir_tables.h"
#include "tap/dsp/kaiser.h"

namespace tap::dsp {

    /// Every design here meets -(A + this) dB on its stopband.
    inline constexpr double k_iir_design_margin_db = 1.0;

    /// The most allpass coefficients a nonlinear half-band design holds:
    /// orders up to 33.
    inline constexpr std::size_t k_max_halfband_coefficients = 16;

    // ------------------------------------------------------------------------
    // The nonlinear (elliptic) two-path half-band
    // ------------------------------------------------------------------------

    /// The selectivity k = tan^2(pi f_pass) and the nome q of the elliptic
    /// half-band whose passband ends at f_pass of its higher rate (Gazsi
    /// 1985: k' = sqrt(1 - k^2), q0 = (1 - k'^(1/2)) / (2 (1 + k'^(1/2))), q =
    /// q0 + 2 q0^5 + 15 q0^9 + 150 q0^13), with 1 - k'^(1/2) computed as
    /// k^2 / ((1 + k'^(1/2)) (1 + k')), which has no cancellation.
    struct halfband_modulus {
        double selectivity = 0.0; ///< k
        double nome        = 0.0; ///< q
        double root        = 0.0; ///< k^(1/2) = tan(pi f_pass), unrounded by the square
    };

    /// @pre 0 < f_pass < 1/4.
    constexpr halfband_modulus halfband_nome(double f_pass) noexcept {
        const double t  = detail::cx::tan_pi(f_pass);
        const double k  = t * t;
        const double kp = detail::cx::sqrt((1.0 - k) * (1.0 + k)); // k'
        const double kk = detail::cx::sqrt(kp);                    // k'^(1/2)
        const double e  = 0.5 * (k * k / ((1.0 + kk) * (1.0 + kp))) / (1.0 + kk);
        const double e4 = e * e * e * e;
        return {k, e * (1.0 + e4 * (2.0 + e4 * (15.0 + 150.0 * e4))), t};
    }

    /// The closed-form stopband attenuation, in dB, of the order-N elliptic
    /// half-band with passband edge f_pass: 10 log10(1 + q^(-N/2) / 4)
    /// (the degree equation q1 = q^N with the bireciprocal ripple tie).
    constexpr double halfband_attenuation_db(int order, double f_pass) noexcept {
        const double lq = -detail::cx::log(halfband_nome(f_pass).nome);
        const double t  = detail::cx::exp(0.5 * static_cast<double>(order) * lq - 2.0 * detail::cx::k_ln2);
        return 10.0 * detail::cx::log1p(t) / detail::cx::k_ln10;
    }

    /// The smallest odd order N >= 3 whose attenuation reaches stop_db at
    /// passband edge f_pass of the higher rate; 0 when it exceeds 2
    /// k_max_halfband_coefficients + 1.
    constexpr int halfband_order(double stop_db, double f_pass) noexcept {
        const double lq   = -detail::cx::log(halfband_nome(f_pass).nome);
        const double need = 2.0 * (2.0 * detail::cx::k_ln2 + detail::cx::log(detail::cx::exp10m1(stop_db / 10.0))) / lq;
        int          n    = 3;
        while (static_cast<double>(n) < need) {
            n += 2;
        }
        return n <= 2 * static_cast<int>(k_max_halfband_coefficients) + 1 ? n : 0;
    }

    /// A nonlinear two-path half-band: (order - 1) / 2 allpass coefficients,
    /// ascending; even indices are path 0's sections, odd indices path 1's.
    struct nonlinear_halfband_design {
        std::array<double, k_max_halfband_coefficients> coefficients{};
        std::size_t                                     count       = 0; ///< (order - 1) / 2
        int                                             order       = 0; ///< N, odd; 0: no design
        double                                          f_pass      = 0; ///< passband edge over the higher rate
        double                                          stopband_db = 0; ///< closed-form attenuation, dB

        /// Sections in path p (0 or 1).
        constexpr std::size_t path_size(std::size_t p) const noexcept { return (count + 1 - p) / 2; }
        /// Path p's i-th section coefficient.
        constexpr double path_coefficient(std::size_t p, std::size_t i) const noexcept {
            return coefficients[2 * i + p];
        }
    };

    /// The order-N elliptic half-band with passband edge f_pass of its higher
    /// rate. @pre order odd, 3 <= order <= 2 k_max_halfband_coefficients + 1,
    /// 0 < f_pass < 1/4.
    constexpr nonlinear_halfband_design design_nonlinear_halfband_order(int order, double f_pass) noexcept {
        nonlinear_halfband_design out{};
        if (order < 3 || order % 2 == 0 || order > 2 * static_cast<int>(k_max_halfband_coefficients) + 1) {
            return out;
        }
        const halfband_modulus m  = halfband_nome(f_pass);
        const double           k  = m.selectivity;
        const double           t  = m.root;
        const double           q  = m.nome;
        const double           q2 = q * q;
        const double           q4 = detail::cx::sqrt(detail::cx::sqrt(q)); // q^(1/4)
        const double           n  = static_cast<double>(order);
        out.order                 = order;
        out.count                 = static_cast<std::size_t>((order - 1) / 2);
        out.f_pass                = f_pass;
        out.stopband_db           = halfband_attenuation_db(order, f_pass);
        for (std::size_t c = 1; c <= out.count; ++c) {
            const double cd = static_cast<double>(c);
            // The critical frequency Omega_c by its theta-function series
            // (Gazsi 1985): 2 q^(1/4) sum_m (-1)^m q^(m(m+1)) sin((2m+1) c pi / n)
            // over 1 + 2 sum_{m>=1} (-1)^m q^(m^2) cos(2 m c pi / n). The
            // powers step by q^(2m) and q^(2m-1); a term below 2^-60 of the
            // leading one ends the sum.
            double num     = detail::cx::sin_pi(cd / n);
            double den     = 1.0;
            double p_num   = 1.0; // q^(m(m+1))
            double p_den   = 1.0; // q^(m^2)
            double inc_num = 1.0; // q^(2m)
            double inc_den = q;   // q^(2m - 1)
            double sign    = 1.0;
            for (int i = 1; i < 12; ++i) {
                const double id = static_cast<double>(i);
                inc_num *= q2;
                p_num *= inc_num;
                p_den *= inc_den;
                inc_den *= q2;
                sign = -sign;
                if (p_den < 8.673617379884035e-19) { // 2^-60, and p_num < p_den
                    break;
                }
                num += sign * p_num * detail::cx::sin_pi((2.0 * id + 1.0) * cd / n);
                den += 2.0 * sign * p_den * detail::cx::cos_pi(2.0 * id * cd / n);
            }
            const double w  = 2.0 * q4 * num / den;
            const double w2 = w * w;
            // a = (1 - x) / (1 + x), x = sqrt((1 - w^2 k)(1 - w^2 / k)) / (1 + w^2)
            // (Gazsi 1985, the bireciprocal case), as (1 - x^2) / (1 + x)^2,
            // with 1 - w^2 / k = (t - w)(t + w) / k (t = k^(1/2) = tan(pi
            // f_pass)): near the band edge w approaches t, and the difference
            // of two doubles is exact there, where 1 - w^2 / k would carry the
            // roundings of both squares, amplified.
            const double x            = detail::cx::sqrt((1.0 - w2 * k) * ((t - w) * (t + w) / k)) / (1.0 + w2);
            const double one_minus_x2 = w2 * ((1.0 + k) * (1.0 + k)) / (k * ((1.0 + w2) * (1.0 + w2)));
            out.coefficients[c - 1]   = one_minus_x2 / ((1.0 + x) * (1.0 + x));
        }
        // Ascending (they come out so; the insertion sort makes it the contract).
        for (std::size_t i = 1; i < out.count; ++i) {
            for (std::size_t j = i; j > 0 && out.coefficients[j] < out.coefficients[j - 1]; --j) {
                const double swap       = out.coefficients[j];
                out.coefficients[j]     = out.coefficients[j - 1];
                out.coefficients[j - 1] = swap;
            }
        }
        return out;
    }

    /// The design at a spec point: stopband A (stop_db), passband edge p of
    /// the chain's lowest rate, design divisor d; the smallest order meeting
    /// A + margin_db. order == 0 when no order up to 33 does or the spec is
    /// out of range (0 < p / d < 1/2).
    constexpr nonlinear_halfband_design design_nonlinear_halfband(double stop_db, double passband, int divisor = 1,
                                                                  double margin_db = k_iir_design_margin_db) noexcept {
        if (!(divisor >= 1) || !(passband > 0.0) || !(passband / static_cast<double>(divisor) < 0.5)) {
            return {};
        }
        const double f_pass = passband / (2.0 * static_cast<double>(divisor));
        const int    order  = halfband_order(stop_db + margin_db, f_pass);
        if (order == 0) {
            return {};
        }
        return design_nonlinear_halfband_order(order, f_pass);
    }

    /// The nonlinear half-band at a tabulated or custom integer spec point,
    /// a compile-time constant: A dB, p = PassNum / PassDen, divisor D.
    template <int StopDb, int PassNum, int PassDen, int Divisor>
    inline constexpr nonlinear_halfband_design k_nonlinear_halfband = design_nonlinear_halfband(
        static_cast<double>(StopDb), static_cast<double>(PassNum) / static_cast<double>(PassDen), Divisor);

    // ------------------------------------------------------------------------
    // The near-linear two-path half-band and the four-path stage: committed
    // ------------------------------------------------------------------------

    namespace detail {

        template <typename Table>
        constexpr const Table* find_spec(std::span<const Table> tables, int stop_db, int pass_num, int pass_den,
                                         int key) noexcept {
            for (const Table& t : tables) {
                if (t.stop_db == stop_db && t.pass_num == pass_num && t.pass_den == pass_den && t.key() == key) {
                    return &t;
                }
            }
            return nullptr;
        }

        template <int StopDb, int PassNum, int PassDen, int Divisor>
        struct committed_near_linear {
            static constexpr const near_linear_halfband_design* k_found =
                find_spec(std::span<const near_linear_halfband_design>(k_near_linear_halfband_tables), StopDb, PassNum,
                          PassDen, Divisor);
            static_assert(k_found != nullptr,
                          "tap::dsp: no committed near-linear half-band at this spec point (iir_tables.h lists them; "
                          "tools/iir/gen_iir_tables.py generates more)");
        };

        template <int StopDb, int PassNum, int PassDen, int Divisor>
        struct committed_fourpath {
            static constexpr const fourpath_design* k_found =
                find_spec(std::span<const fourpath_design>(k_fourpath_tables), StopDb, PassNum, PassDen, Divisor);
            static_assert(k_found != nullptr,
                          "tap::dsp: no committed four-path design at this spec point (iir_tables.h lists them)");
        };

    } // namespace detail

    /// The committed near-linear half-band at a spec point (a static_assert
    /// names a missing one).
    template <int StopDb, int PassNum, int PassDen, int Divisor>
    inline constexpr const near_linear_halfband_design& k_near_linear_halfband =
        *detail::committed_near_linear<StopDb, PassNum, PassDen, Divisor>::k_found;

    /// The committed four-path design at a spec point.
    template <int StopDb, int PassNum, int PassDen, int Divisor>
    inline constexpr const fourpath_design& k_fourpath =
        *detail::committed_fourpath<StopDb, PassNum, PassDen, Divisor>::k_found;

    // ------------------------------------------------------------------------
    // The guard band filter: the odd-order Chebyshev II
    // ------------------------------------------------------------------------

    /// One second-order section (b0 + b1 z^-1 + b2 z^-2) / (1 + a1 z^-1 + a2
    /// z^-2); a first-order section has b2 == a2 == 0.
    struct iir_section {
        double b0 = 0.0;
        double b1 = 0.0;
        double b2 = 0.0;
        double a1 = 0.0;
        double a2 = 0.0;
    };

    /// An odd-order Chebyshev II low-pass as sections ordered by increasing
    /// pole radius (the first-order section among them, where its radius
    /// falls).
    struct cheby2_design {
        std::array<iir_section, k_max_cheby2_sections> sections{};
        std::array<double, k_max_cheby2_sections>      pole_radius{};
        std::size_t                                    count     = 0; ///< (order + 1) / 2
        int                                            order     = 0; ///< odd; 0: no design
        double                                         stop_edge = 0; ///< where it reaches A, a fraction of its rate
    };

    namespace detail {

        struct cplx {
            double re = 0.0;
            double im = 0.0;
        };
        constexpr cplx operator+(cplx a, cplx b) noexcept {
            return {a.re + b.re, a.im + b.im};
        }
        constexpr cplx operator-(cplx a, cplx b) noexcept {
            return {a.re - b.re, a.im - b.im};
        }
        constexpr cplx operator/(cplx a, cplx b) noexcept {
            const double d = b.re * b.re + b.im * b.im;
            return {(a.re * b.re + a.im * b.im) / d, (a.im * b.re - a.re * b.im) / d};
        }
        constexpr double norm(cplx a) noexcept {
            return a.re * a.re + a.im * a.im;
        }

    } // namespace detail

    /// The order-n Chebyshev II low-pass flat within ripple_db to passband p
    /// (a fraction of its rate) and stop_db down from its stopband edge, every
    /// section at unity DC gain (the L-inf gains apply on top:
    /// apply_section_gains). @pre order odd, 3 <= order <= 2
    /// k_max_cheby2_sections - 1, 0 < p < 1/2, stop_db > ripple_db > 0.
    constexpr cheby2_design design_cheby2_lowpass(int order, double passband, double stop_db,
                                                  double ripple_db = 0.01) noexcept {
        cheby2_design out{};
        if (order < 3 || order % 2 == 0 || order > 2 * static_cast<int>(k_max_cheby2_sections) - 1
            || !(passband > 0.0 && passband < 0.5) || !(stop_db > ripple_db && ripple_db > 0.0)) {
            return out;
        }
        namespace cx           = detail::cx;
        const double      n    = static_cast<double>(order);
        const std::size_t m    = static_cast<std::size_t>((order - 1) / 2);
        const double      a10  = cx::exp10m1(stop_db / 10.0);   // 10^(A/10) - 1 = 1 / eps^2
        const double      r10  = cx::exp10m1(ripple_db / 10.0); // 10^(rp/10) - 1
        const double      op   = cx::tan_pi(passband);          // prewarped passband edge
        const double      dd   = cx::sqrt(a10 / r10);
        const double      os   = op * cx::cosh(cx::acosh(dd) / n); // stopband edge
        const double      mu   = cx::asinh(cx::sqrt(a10)) / n;     // asinh(1 / eps) / n
        const double      sh   = cx::sinh(mu);
        const double      ch   = cx::cosh(mu);
        out.order              = order;
        out.count              = m + 1;
        out.stop_edge          = cx::atan_pi(os);
        const detail::cplx one = {1.0, 0.0};
        for (std::size_t k = 1; k <= m; ++k) {
            const double       th = static_cast<double>(2 * k - 1) / (2.0 * n); // theta_k / pi
            const double       st = cx::sin_pi(th);
            const double       ct = cx::cos_pi(th);
            const detail::cplx cheb{-sh * st, ch * ct};                 // Chebyshev I pole
            const detail::cplx pz  = detail::cplx{os, 0.0} / cheb;      // the analogue pole
            const detail::cplx zp  = (one + pz) / (one - pz);           // bilinear
            const double       oz  = os / ct;                           // the analogue zero j oz
            const double       cz  = (1.0 - oz * oz) / (1.0 + oz * oz); // Re of its image on |z| = 1
            iir_section&       s   = out.sections[k - 1];
            const double       g   = detail::norm(one - zp) / (2.0 * (1.0 - cz)); // unity DC
            s.b0                   = g;
            s.b1                   = -2.0 * cz * g;
            s.b2                   = g;
            s.a1                   = -2.0 * zp.re;
            s.a2                   = detail::norm(zp);
            out.pole_radius[k - 1] = cx::sqrt(s.a2);
        }
        {
            const double pr    = -os / sh; // the real pole, theta = pi / 2
            const double zr    = (1.0 + pr) / (1.0 - pr);
            iir_section& s     = out.sections[m];
            const double g     = 0.5 * (1.0 - zr);
            s.b0               = g;
            s.b1               = g;
            s.b2               = 0.0;
            s.a1               = -zr;
            s.a2               = 0.0;
            out.pole_radius[m] = cx::abs(zr);
        }
        for (std::size_t i = 1; i < out.count; ++i) {
            for (std::size_t j = i; j > 0 && out.pole_radius[j] < out.pole_radius[j - 1]; --j) {
                const iir_section ts   = out.sections[j];
                out.sections[j]        = out.sections[j - 1];
                out.sections[j - 1]    = ts;
                const double tr        = out.pole_radius[j];
                out.pole_radius[j]     = out.pole_radius[j - 1];
                out.pole_radius[j - 1] = tr;
            }
        }
        return out;
    }

    /// Scales section i's numerator by gains[i] for i < count - 1 and the
    /// last section's by 1 / prod(gains): the L-inf distribution (Jackson
    /// 1970) when gains are each partial cascade's 1 / peak; the DC gain of
    /// the whole cascade is unchanged.
    constexpr cheby2_design apply_section_gains(cheby2_design d, std::span<const double> gains) noexcept {
        if (d.count == 0) {
            return d;
        }
        double carried = 1.0;
        for (std::size_t i = 0; i + 1 < d.count && i < gains.size(); ++i) {
            d.sections[i].b0 *= gains[i];
            d.sections[i].b1 *= gains[i];
            d.sections[i].b2 *= gains[i];
            carried *= gains[i];
        }
        iir_section& last = d.sections[d.count - 1];
        last.b0 /= carried;
        last.b1 /= carried;
        last.b2 /= carried;
        return d;
    }

    namespace detail {

        template <int StopDb, int PassNum, int PassDen, int Order>
        struct committed_cheby2 {
            static constexpr const cheby2_gains* k_found =
                find_spec(std::span<const cheby2_gains>(k_cheby2_gain_tables), StopDb, PassNum, PassDen, Order);
            static_assert(k_found != nullptr, "tap::dsp: no committed L-inf gains for this Chebyshev II (iir_tables.h "
                                              "lists them; cheby2_linf_gains() computes a custom design's)");
            static constexpr cheby2_design k_design = apply_section_gains(
                design_cheby2_lowpass(Order, static_cast<double>(PassNum) / static_cast<double>(PassDen),
                                      static_cast<double>(StopDb)),
                std::span<const double>(k_found->gains.data(), static_cast<std::size_t>((Order - 1) / 2)));
        };

        template <int StopDb, int PassNum, int PassDen, int Order>
        struct committed_guard_fir {
            static constexpr const guard_band_fir_design* k_found = find_spec(
                std::span<const guard_band_fir_design>(k_guard_band_fir_tables), StopDb, PassNum, PassDen, Order);
            static_assert(k_found != nullptr, "tap::dsp: no committed linear guard band filter at this spec point "
                                              "(iir_tables.h lists them; design_guard_band_linear() designs a custom "
                                              "one)");
        };

        template <int StopDb, int PassNum, int PassDen, int Order>
        struct committed_guard_fir_q15 {
            static constexpr const guard_band_fir_q15_design* k_found =
                find_spec(std::span<const guard_band_fir_q15_design>(k_guard_band_fir_q15_tables), StopDb, PassNum,
                          PassDen, Order);
            static_assert(k_found != nullptr, "tap::dsp: no committed Q15 linear guard band filter at this spec point "
                                              "(iir_tables.h lists them: the three 70 dB points, orders 3 to 11; none "
                                              "at 120 dB)");
        };

    } // namespace detail

    /// The guard band filter's Chebyshev II at a spec point (A dB, passband
    /// PassNum / PassDen of its rate, the odd Order), its L-inf gains the
    /// committed ones: a compile-time constant.
    template <int StopDb, int PassNum, int PassDen, int Order>
    inline constexpr cheby2_design k_guard_band_cheby2 =
        detail::committed_cheby2<StopDb, PassNum, PassDen, Order>::k_design;

    /// The committed linear (type II) guard band filter whose stop edge is
    /// the order-Order Chebyshev II's at the same spec point.
    template <int StopDb, int PassNum, int PassDen, int Order>
    inline constexpr const guard_band_fir_design& k_guard_band_fir =
        *detail::committed_guard_fir<StopDb, PassNum, PassDen, Order>::k_found;

    /// The committed linear (type II) guard band filter for Q15 at the same
    /// spec point and stop edge: its taps, quantized to Q1.14 by
    /// quantize_symmetric_row_preserving_sum<std::int16_t>, keep the spec
    /// k_guard_band_fir's double taps lose there (iir_tables.h).
    template <int StopDb, int PassNum, int PassDen, int Order>
    inline constexpr const guard_band_fir_q15_design& k_guard_band_fir_q15 =
        *detail::committed_guard_fir_q15<StopDb, PassNum, PassDen, Order>::k_found;

    // ------------------------------------------------------------------------
    // Responses and group delays (runtime, libm)
    // ------------------------------------------------------------------------

    /// prod (a + w^-1) / (1 + a w^-1) at w = e^{j theta}.
    inline std::complex<double> allpass1_cascade_response(std::span<const double> a, double theta) noexcept {
        const std::complex<double> wi = std::polar(1.0, -theta);
        std::complex<double>       r  = 1.0;
        for (const double c : a) {
            r *= (c + wi) / (1.0 + c * wi);
        }
        return r;
    }

    /// The group delay of that cascade, in path samples (theta's unit).
    inline double allpass1_cascade_group_delay(std::span<const double> a, double theta) noexcept {
        double tau = 0.0;
        for (const double c : a) {
            tau += (1.0 - c * c) / (1.0 + 2.0 * c * std::cos(theta) + c * c);
        }
        return tau;
    }

    /// The near-linear path A1 (second-order sections (c2 + c1 w^-1 + w^-2) /
    /// (1 + c1 w^-1 + c2 w^-2), then one first-order c1 when N is odd).
    inline std::complex<double> allpass2_cascade_response(std::span<const double> sections, std::size_t order,
                                                          double theta) noexcept {
        const std::complex<double> wi  = std::polar(1.0, -theta);
        const std::complex<double> wi2 = wi * wi;
        std::complex<double>       r   = 1.0;
        std::size_t                i   = 0;
        for (std::size_t left = order; left > 0;) {
            if (left >= 2) {
                const double c1 = sections[i];
                const double c2 = sections[i + 1];
                r *= (c2 + c1 * wi + wi2) / (1.0 + c1 * wi + c2 * wi2);
                i += 2;
                left -= 2;
            }
            else {
                const double c1 = sections[i];
                r *= (c1 + wi) / (1.0 + c1 * wi);
                i += 1;
                left -= 1;
            }
        }
        return r;
    }

    /// Its group delay in path samples, summed over the sections: an allpass
    /// w^-n D(1/w) / D(w), D = sum_k d_k w^-k, delays by n - 2 Re(S / D) with
    /// S = sum_k k d_k w^-k (the derivative of its phase, -n theta - 2 arg D).
    inline double allpass2_cascade_group_delay(std::span<const double> sections, std::size_t order,
                                               double theta) noexcept {
        const std::complex<double> wi  = std::polar(1.0, -theta);
        const std::complex<double> wi2 = wi * wi;
        double                     tau = 0.0;
        std::size_t                i   = 0;
        for (std::size_t left = order; left > 0;) {
            if (left >= 2) {
                const double               c1 = sections[i];
                const double               c2 = sections[i + 1];
                const std::complex<double> d  = 1.0 + c1 * wi + c2 * wi2;
                const std::complex<double> dd = c1 * wi + 2.0 * c2 * wi2;
                tau += 2.0 - 2.0 * std::real(dd / d);
                i += 2;
                left -= 2;
            }
            else {
                const double               c1 = sections[i];
                const std::complex<double> d  = 1.0 + c1 * wi;
                const std::complex<double> dd = c1 * wi;
                tau += 1.0 - 2.0 * std::real(dd / d);
                i += 1;
                left -= 1;
            }
        }
        return tau;
    }

    /// The nonlinear half-band's response at f, cycles per high-rate sample.
    inline std::complex<double> halfband_response(const nonlinear_halfband_design& d, double f) noexcept {
        std::array<double, k_max_halfband_coefficients> p0{};
        std::array<double, k_max_halfband_coefficients> p1{};
        const std::size_t                               n0 = d.path_size(0);
        const std::size_t                               n1 = d.path_size(1);
        for (std::size_t i = 0; i < n0; ++i) {
            p0[i] = d.path_coefficient(0, i);
        }
        for (std::size_t i = 0; i < n1; ++i) {
            p1[i] = d.path_coefficient(1, i);
        }
        const double theta = 4.0 * std::numbers::pi * f;
        const auto   a0    = allpass1_cascade_response(std::span<const double>(p0.data(), n0), theta);
        const auto   a1    = allpass1_cascade_response(std::span<const double>(p1.data(), n1), theta);
        return 0.5 * (a0 + a1 * std::polar(1.0, -2.0 * std::numbers::pi * f));
    }

    /// Its group delay in high-rate samples: -d arg H / d omega, exact from
    /// the paths' delays (H' / H).
    inline double halfband_group_delay(const nonlinear_halfband_design& d, double f) noexcept {
        std::array<double, k_max_halfband_coefficients> p0{};
        std::array<double, k_max_halfband_coefficients> p1{};
        const std::size_t                               n0 = d.path_size(0);
        const std::size_t                               n1 = d.path_size(1);
        for (std::size_t i = 0; i < n0; ++i) {
            p0[i] = d.path_coefficient(0, i);
        }
        for (std::size_t i = 0; i < n1; ++i) {
            p1[i] = d.path_coefficient(1, i);
        }
        const double theta = 4.0 * std::numbers::pi * f;
        const auto   a0    = allpass1_cascade_response(std::span<const double>(p0.data(), n0), theta);
        const auto   a1    = allpass1_cascade_response(std::span<const double>(p1.data(), n1), theta);
        const double t0    = 2.0 * allpass1_cascade_group_delay(std::span<const double>(p0.data(), n0), theta);
        const double t1    = 2.0 * allpass1_cascade_group_delay(std::span<const double>(p1.data(), n1), theta);
        const std::complex<double> zi = std::polar(1.0, -2.0 * std::numbers::pi * f);
        const std::complex<double> j{0.0, 1.0};
        const std::complex<double> h  = 0.5 * (a0 + a1 * zi);
        const std::complex<double> dh = 0.5 * (-j * t0 * a0 - j * (1.0 + t1) * a1 * zi);
        return -std::imag(dh / h);
    }

    /// The near-linear half-band's response at f (cycles per high-rate sample).
    inline std::complex<double> halfband_response(const near_linear_halfband_design& d, double f) noexcept {
        const double theta = 4.0 * std::numbers::pi * f;
        const auto   a0    = std::polar(1.0, -static_cast<double>(d.delay) * theta);
        const auto   a1    = allpass2_cascade_response(d.sections, static_cast<std::size_t>(d.order), theta);
        return 0.5 * (a0 + a1 * std::polar(1.0, -2.0 * std::numbers::pi * f));
    }

    /// Its group delay in high-rate samples.
    inline double halfband_group_delay(const near_linear_halfband_design& d, double f) noexcept {
        const double theta = 4.0 * std::numbers::pi * f;
        const auto   a0    = std::polar(1.0, -static_cast<double>(d.delay) * theta);
        const auto   a1    = allpass2_cascade_response(d.sections, static_cast<std::size_t>(d.order), theta);
        const double t0    = 2.0 * static_cast<double>(d.delay);
        const double t1    = 2.0 * allpass2_cascade_group_delay(d.sections, static_cast<std::size_t>(d.order), theta);
        const std::complex<double> zi = std::polar(1.0, -2.0 * std::numbers::pi * f);
        const std::complex<double> j{0.0, 1.0};
        const std::complex<double> h  = 0.5 * (a0 + a1 * zi);
        const std::complex<double> dh = 0.5 * (-j * t0 * a0 - j * (1.0 + t1) * a1 * zi);
        return -std::imag(dh / h);
    }

    /// The four-path stage's response at f, cycles per high-rate (input) sample.
    inline std::complex<double> fourpath_response(const fourpath_design& d, double f) noexcept {
        const double         theta = 8.0 * std::numbers::pi * f;
        std::complex<double> h     = 0.0;
        for (std::size_t k = 0; k < 4; ++k) {
            h += allpass1_cascade_response(d.paths[k], theta)
                 * std::polar(1.0, -2.0 * std::numbers::pi * f * static_cast<double>(k));
        }
        return 0.25 * h;
    }

    /// Its group delay in high-rate samples.
    inline double fourpath_group_delay(const fourpath_design& d, double f) noexcept {
        const double               theta = 8.0 * std::numbers::pi * f;
        const std::complex<double> j{0.0, 1.0};
        std::complex<double>       h  = 0.0;
        std::complex<double>       dh = 0.0;
        for (std::size_t k = 0; k < 4; ++k) {
            const auto a = allpass1_cascade_response(d.paths[k], theta)
                           * std::polar(1.0, -2.0 * std::numbers::pi * f * static_cast<double>(k));
            const double tau = static_cast<double>(k) + 4.0 * allpass1_cascade_group_delay(d.paths[k], theta);
            h += a;
            dh -= j * tau * a;
        }
        return -std::imag(dh / h);
    }

    /// A cascade of sections' response at f, cycles per sample.
    inline std::complex<double> sections_response(std::span<const iir_section> s, double f) noexcept {
        const std::complex<double> zi  = std::polar(1.0, -2.0 * std::numbers::pi * f);
        const std::complex<double> zi2 = zi * zi;
        std::complex<double>       r   = 1.0;
        for (const iir_section& c : s) {
            r *= (c.b0 + c.b1 * zi + c.b2 * zi2) / (1.0 + c.a1 * zi + c.a2 * zi2);
        }
        return r;
    }

    /// A cascade of sections' group delay at f, in samples: per section the
    /// numerator's delay less the denominator's, a polynomial P = sum_k p_k
    /// z^-k delaying by Re(sum_k k p_k z^-k / P). Away from the numerators'
    /// zeros (in the passband).
    inline double sections_group_delay(std::span<const iir_section> s, double f) noexcept {
        const std::complex<double> zi  = std::polar(1.0, -2.0 * std::numbers::pi * f);
        const std::complex<double> zi2 = zi * zi;
        double                     tau = 0.0;
        for (const iir_section& c : s) {
            const std::complex<double> b = c.b0 + c.b1 * zi + c.b2 * zi2;
            const std::complex<double> a = 1.0 + c.a1 * zi + c.a2 * zi2;
            tau += std::real((c.b1 * zi + 2.0 * c.b2 * zi2) / b) - std::real((c.a1 * zi + 2.0 * c.a2 * zi2) / a);
        }
        return tau;
    }

    /// The Chebyshev II design's response at f (cycles per sample).
    inline std::complex<double> cheby2_response(const cheby2_design& d, double f) noexcept {
        return sections_response(std::span<const iir_section>(d.sections.data(), d.count), f);
    }

    /// The L-inf gains of a Chebyshev II design at unity DC gain per section
    /// (design_cheby2_lowpass), at run time: gains[i] = 1 / (the peak of the
    /// partial cascade 0 .. i, its earlier sections already scaled), for i <
    /// count - 1, the peak found on a grid of grid_points + 1 frequencies over
    /// [0, 1/2] and refined by golden-section search around the largest. The
    /// committed tables were computed the same way at 60 digits
    /// (tools/iir/gen_iir_tables.py). @pre gains.size() >= count - 1.
    inline void cheby2_linf_gains(const cheby2_design& d, std::span<double> gains,
                                  std::size_t grid_points = 4096) noexcept {
        cheby2_design scaled = d;
        for (std::size_t i = 0; i + 1 < d.count; ++i) {
            const auto partial = [&](double f) noexcept {
                return std::abs(sections_response(std::span<const iir_section>(scaled.sections.data(), i + 1), f));
            };
            std::size_t best = 0;
            double      peak = 0.0;
            for (std::size_t g = 0; g <= grid_points; ++g) {
                const double v = partial(0.5 * static_cast<double>(g) / static_cast<double>(grid_points));
                if (v > peak) {
                    peak = v;
                    best = g;
                }
            }
            // Golden-section refinement over the neighbouring grid cells.
            double lo = 0.5 * static_cast<double>(best > 0 ? best - 1 : 0) / static_cast<double>(grid_points);
            double hi = 0.5 * static_cast<double>(best < grid_points ? best + 1 : grid_points)
                        / static_cast<double>(grid_points);
            constexpr double inv_phi = 0.6180339887498949;
            for (int it = 0; it < 80; ++it) {
                const double m1 = hi - inv_phi * (hi - lo);
                const double m2 = lo + inv_phi * (hi - lo);
                if (partial(m1) < partial(m2)) {
                    lo = m1;
                }
                else {
                    hi = m2;
                }
            }
            peak            = std::max(peak, partial(0.5 * (lo + hi)));
            const double sc = 1.0 / peak;
            gains[i]        = sc;
            scaled.sections[i].b0 *= sc;
            scaled.sections[i].b1 *= sc;
            scaled.sections[i].b2 *= sc;
        }
    }

    // ------------------------------------------------------------------------
    // The linear guard band form, custom: a Kaiser-windowed type II FIR
    // ------------------------------------------------------------------------

    /// Writes the even-length (type II) Kaiser-windowed sinc low-pass with
    /// its cutoff midway between passband p and stop edge f_stop (fractions of
    /// its rate), its beta for stop_db + margin (kaiser.h), DC gain 1. Linear
    /// phase, symmetric (h[i] == h[N - 1 - i] bit for bit: the second half is
    /// copied), a zero at Nyquist by its symmetry; latency (N - 1) / 2
    /// samples, a half-integer. @pre h.size() even, >= 2.
    inline void design_guard_band_linear(std::span<double> h, double passband, double f_stop, double stop_db,
                                         double margin_db = k_iir_design_margin_db) noexcept {
        const std::size_t n = h.size();
        // kaiser.h's prototype with one phase: cutoff_norm is 2 fc in units of
        // the rate's Nyquist.
        design_prototype(h, 1, passband + f_stop, kaiser_beta(stop_db + margin_db));
        for (std::size_t i = 0; i < n / 2; ++i) {
            h[n - 1 - i] = h[i];
        }
        double sum = 0.0;
        for (const double v : h) {
            sum += v;
        }
        for (double& v : h) {
            v /= sum;
        }
    }

    /// The worst stopband of a FIR over [f_stop, 1/2] on grid_points + 1
    /// frequencies, in dB, and its peak-to-peak passband ripple over [0, p].
    struct fir_band_measure {
        double stopband_db = 0.0;
        double ripple_db   = 0.0;
    };
    inline fir_band_measure measure_fir_bands(std::span<const double> h, double passband, double f_stop,
                                              std::size_t grid_points = 4096) noexcept {
        const auto mag = [&](double f) noexcept {
            std::complex<double> acc = 0.0;
            for (std::size_t i = 0; i < h.size(); ++i) {
                acc += h[i] * std::polar(1.0, -2.0 * std::numbers::pi * f * static_cast<double>(i));
            }
            return std::abs(acc);
        };
        double worst = 0.0;
        double lo    = 1e300;
        double hi    = 0.0;
        for (std::size_t g = 0; g <= grid_points; ++g) {
            const double t = static_cast<double>(g) / static_cast<double>(grid_points);
            worst          = std::max(worst, mag(f_stop + (0.5 - f_stop) * t));
            const double p = mag(passband * t);
            lo             = std::min(lo, p);
            hi             = std::max(hi, p);
        }
        return {20.0 * std::log10(std::max(worst, 1e-300)), 20.0 * std::log10(hi / lo)};
    }

    /// The shortest even length (from 2) whose design_guard_band_linear meets
    /// stop_db + margin_db on grid_points; allocates (design time only).
    inline std::size_t search_guard_band_linear_length(double passband, double f_stop, double stop_db,
                                                       double      margin_db   = k_iir_design_margin_db,
                                                       std::size_t grid_points = 4096, std::size_t max_taps = 512) {
        std::vector<double> h;
        for (std::size_t n = 2; n <= max_taps; n += 2) {
            h.assign(n, 0.0);
            design_guard_band_linear(h, passband, f_stop, stop_db, margin_db);
            if (measure_fir_bands(h, passband, f_stop, grid_points).stopband_db <= -(stop_db + margin_db)) {
                return n;
            }
        }
        return 0;
    }

    // ------------------------------------------------------------------------
    // The flush count of a floating-point state (IIR_PLAN I6, recommended)
    // ------------------------------------------------------------------------

    /// Path frames from the format's largest finite value to its smallest
    /// normal value through the largest pole radius r (in path samples), plus
    /// one frame per section: ceil(ln(max / min_normal) / ln(1 / r)) +
    /// sections. Below the smallest normal the host's flush-to-zero
    /// (subnormals.h) zeroes the state; without it the state sticks (the
    /// header's honest limit). @pre 0 < r < 1.
    template <typename Float>
    inline std::size_t floating_flush_frames(double r, std::size_t sections) noexcept {
        const double span = std::log(static_cast<double>(std::numeric_limits<Float>::max()))
                            - std::log(static_cast<double>(std::numeric_limits<Float>::min()));
        return static_cast<std::size_t>(std::ceil(span / -std::log(r))) + sections;
    }

} // namespace tap::dsp
