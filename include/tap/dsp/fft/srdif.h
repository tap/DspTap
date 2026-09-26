/// @file srdif.h
/// @brief The floating real FFT engine: a split-radix DIF kernel of length N/2 and the real post-pass.
// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// detail::srdif_rdft<Sample> is the engine basic_real_fft<double> runs on
// every build and basic_real_fft<float> runs wherever the build selects no
// accelerated backend (fft.h). It replaced, at the same contract, a port of a
// third-party split-radix engine, and was written without reference to that
// engine or to any other FFT library: from the literature below and from
// DspTap's own fixed-point engine (fft/fixed_point.h), whose structure it
// generalizes to floating point (docs/fft-design.md, "The floating engine
// (srdif)").
//
// The contract it presents (fft.h is the authority): N a power of two in
// [4, 2^30]; packing a[0] = Re X[0], a[1] = Re X[N/2], a[2k] + i a[2k+1] =
// X[k] for 1 <= k < N/2; forward X[k] = sum_j x[j] W^jk with W = exp(+2 pi i/N);
// an unnormalized inverse (forward then inverse is N/2 times the input).
//
// Implemented from the published literature only (house IP policy), by part:
//
//   The real transform through a half-length complex one, and its inverse:
//     J. W. Cooley, P. A. W. Lewis and P. D. Welch, "The fast Fourier
//       transform algorithm: Programming considerations in the calculation
//       of sine, cosine and Laplace transforms," J. Sound Vib. 12(3),
//       315-337, 1970.
//     H. V. Sorensen, D. L. Jones, M. T. Heideman and C. S. Burrus,
//       "Real-valued fast Fourier transform algorithms," IEEE Trans. ASSP
//       35(6), 849-863, 1987.
//     The floating statements are fixed_point.h's real_post_pass (which
//     derives them from those equations, tap/DspTap#39) with the scaling
//     removed: the floating contract has no exponent, so no one-bit shift
//     precedes the pass and nothing is halved except the inverse's
//     DC / Nyquist pair.
//
//   The complex kernel, split-radix decimation in frequency:
//     P. Duhamel and H. Hollmann, "Split radix FFT algorithm," Electron.
//       Lett. 20(1), 14-16, 1984 (the decomposition: X[2k] from a half-length
//       transform, X[4k+1] and X[4k+3] from two quarter-length ones).
//     H. V. Sorensen, M. T. Heideman and C. S. Burrus, "On computing the
//       split-radix FFT," IEEE Trans. ASSP 34(1), 152-156, 1986 (the
//       decimation-in-frequency "L-shaped" butterfly, the special butterflies
//       at the trivial and the eighth-turn twiddles, and the operation count
//       4M log2 M - 6M + 8 real additions and multiplications that this
//       kernel attains exactly).
//     P. Duhamel and M. Vetterli, "Fast Fourier transforms: a tutorial
//       review and a state of the art," Signal Processing 19, 259-299, 1990
//       (the index maps behind the fused pass and the bit-reversed output
//       order of in-place decimation in frequency).
//
//   The permutation: A. H. Karp, "Bit reversal on uniprocessors," SIAM
//     Review 38(1), 1-26, 1996 (the table-driven swap method).
//
//   The tables: the Maclaurin series of cos and sin evaluated by Horner's
//     rule in 64-bit fixed point (D. E. Knuth, The Art of Computer
//     Programming vol. 2, 3rd ed., Sec. 4.3.1 for the multiple-precision
//     products and Sec. 4.6.4 for Horner's rule), and the exact symmetries
//     of the unit circle for everything outside the first octant, as
//     tables.h does for the fixed-point tables.
//
// The arrangement of these pieces (the fused two-level pass, the twiddle
// pairing, the compile-time blocks and register leaves, the table layouts)
// is derived in the docstrings below.

#pragma once

#include <array>
#include <bit>
#include <cassert>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <type_traits>
#include <utility>
#include <vector>

// The transforms are kept out of line (a performance attribute only: each is
// one call per transform, and inlined into a caller's loop the permutation
// and the post-pass lose their registers to the caller's live values).
#if defined(__GNUC__) || defined(__clang__)
#define TAP_DSP_SRDIF_NOINLINE [[gnu::noinline]]
#elif defined(_MSC_VER)
#define TAP_DSP_SRDIF_NOINLINE __declspec(noinline)
#else
#define TAP_DSP_SRDIF_NOINLINE
#endif

namespace tap::dsp::detail {

    /// The engine's trigonometry: cos and sin of 2 pi k / 2^L in integer
    /// arithmetic, so the tables are the same bits on every host.
    ///
    /// Why not libm. A table built from std::cos / std::sin is only as
    /// portable as the last bit of each libm (glibc, newlib, UCRT and Apple
    /// differ there, and x86-64 glibc differs from itself between its CPU
    /// dispatches), so every output bit of a transform would be pinned per C
    /// library. Evaluated here on unsigned 64-bit integers, the table is a
    /// function of (k, L) alone; the one floating step is to_sample's exact
    /// conversion of an at most 53-bit integer and its exact power-of-two
    /// scaling. Together with arithmetic in a fixed order this makes the
    /// engine's output bits a function of the source and the compiler's
    /// fp-contraction choice only (tests/test_fft_srdif_fingerprint.cpp).
    /// It also keeps libm out of the float profile's image: a Cortex-M
    /// build of basic_real_fft<float> links no sin / cos at all.
    ///
    /// Method (first_octant). theta = 2 pi k / 2^L, 0 < k <= 2^L / 8, is
    /// formed exactly as the 128-bit product P = k * floor(2 pi 2^61) (a
    /// relative error of 2^-63 from the truncated constant), read two ways:
    /// as Q0.64 (absolute; theta <= pi/4 < 1) for y = theta^2, and
    /// normalized (relative) for the products that must keep precision at
    /// small angles. Horner's rule in Q0.64 gives tc = (1 - cos) / y and
    /// ts = (1 - sin/theta) / y from the Maclaurin coefficients 1/(2j)!
    /// (j <= 10) and 1/(2j+1)! (j <= 9), each floor(2^64 / n!) (the
    /// truncated terms are below 2^-66 at y <= (pi/4)^2); every
    /// intermediate is positive and below 1. Then cos = 1 - y tc and
    /// 1 - sin = 1 - theta (1 - y ts) in Q0.64 (both at least 0.29 in the
    /// first octant, so absolute precision is relative precision there),
    /// sin = theta (1 - y ts) and 1 - cos = theta^2 tc with the normalized
    /// theta, so a small angle keeps full relative precision in both.
    ///
    /// Accuracy, measured against __float128 (libquadmath) over every k of
    /// every 2^L <= 2^18 and 200,000 random k per L = 19 … 30 (2026-09-26):
    /// the 64-bit mantissa is within 12.7 units of its last place of the
    /// exact value (1 - cos is the worst of the four; cos 2.7, sin 4.6,
    /// 1 - sin 4.8). Rounded by to_sample: every float value is the correctly
    /// rounded one (no exact value lies within 2.2e4 units of a float
    /// rounding tie); a double value is within 0.5 + 2^-7.3 ulp of the exact
    /// one, and 0.08 % of them are the neighbour of the correctly rounded
    /// double. Deterministic either way. `GeneratorAgreesWithLibmToTheLastBits`,
    /// `OctantEndpointsAreExact`, `KernelLiteralsAreTheGeneratorsValues`,
    /// `ToSampleRoundsHalfToEven` (tests/test_fft_srdif.cpp).
    ///
    /// Cost: about 25 64 x 64 -> 128 high products per angle (four 32-bit
    /// multiplies each); the engine evaluates about 3N/16 angles at
    /// construction. Measured on the Cortex-M4F (QEMU, -O3), construction
    /// whole, allocation included: 103 k instructions at N = 512, 396 k at
    /// N = 2048 (0.11 % / 0.37 % of the ratchet scenarios).
    namespace srdif_trig {

        /// A non-negative value m * 2^e, m normalized (bit 63 set) or zero.
        struct wide_real {
            std::uint64_t m;
            int           e;
        };

        /// floor(a * b / 2^64), exact, from four 32 x 32 -> 64 products.
        constexpr std::uint64_t mul_hi(std::uint64_t a, std::uint64_t b) noexcept {
            const std::uint64_t a_lo = a & 0xffffffffu;
            const std::uint64_t a_hi = a >> 32;
            const std::uint64_t b_lo = b & 0xffffffffu;
            const std::uint64_t b_hi = b >> 32;
            const std::uint64_t ll   = a_lo * b_lo;
            const std::uint64_t lh   = a_lo * b_hi;
            const std::uint64_t hl   = a_hi * b_lo;
            const std::uint64_t hh   = a_hi * b_hi;
            const std::uint64_t mid  = (ll >> 32) + (lh & 0xffffffffu) + (hl & 0xffffffffu);
            return hh + (lh >> 32) + (hl >> 32) + (mid >> 32);
        }

        /// m * 2^e normalized so that bit 63 of m is set (m != 0).
        constexpr wide_real normalized(std::uint64_t m, int e) noexcept {
            if (m == 0) {
                return {0, 0};
            }
            const int z = std::countl_zero(m);
            return {m << z, e - z};
        }

        /// floor(2^64 / d) for d >= 2.
        constexpr std::uint64_t reciprocal_q64(std::uint64_t d) noexcept {
            const std::uint64_t all = ~std::uint64_t{0};
            std::uint64_t       q   = all / d;
            if (all % d == d - 1) {
                ++q; // d divides 2^64
            }
            return q;
        }

        constexpr std::uint64_t factorial(int n) noexcept {
            std::uint64_t f = 1;
            for (int i = 2; i <= n; ++i) {
                f *= static_cast<std::uint64_t>(i);
            }
            return f;
        }

        /// floor(2^64 / (2j)!) for j = 0 .. 10 (entry 0 unused) and floor(2^64 / (2j+1)!) for
        /// j = 0 .. 9 (entry 0 unused): the Taylor coefficients of cos and sin in Q0.64.
        inline constexpr std::array<std::uint64_t, 11> k_cos_coeff = [] {
            std::array<std::uint64_t, 11> c{};
            for (int j = 1; j <= 10; ++j) {
                c[static_cast<std::size_t>(j)] = reciprocal_q64(factorial(2 * j));
            }
            return c;
        }();
        inline constexpr std::array<std::uint64_t, 10> k_sin_coeff = [] {
            std::array<std::uint64_t, 10> c{};
            for (int j = 1; j <= 9; ++j) {
                c[static_cast<std::size_t>(j)] = reciprocal_q64(factorial(2 * j + 1));
            }
            return c;
        }();

        /// 2*pi in Q3.61, truncated: floor(2 pi 2^61).
        inline constexpr std::uint64_t k_two_pi_q61 = 0xc90fdaa22168c234ull;

        /// The four values the tables need at one angle theta = 2 pi k / 2^log2n,
        /// 0 <= k <= 2^log2n / 8 (the first octant), each with relative precision.
        struct octant_values {
            wide_real cos;
            wide_real sin;
            wide_real one_minus_cos;
            wide_real one_minus_sin;
        };

        /// cos, sin, 1 - cos and 1 - sin of theta = 2 pi k / 2^log2n.
        /// @pre 0 <= k <= 2^log2n / 8, 3 <= log2n <= 33.
        constexpr octant_values first_octant(std::uint64_t k, int log2n) noexcept {
            if (k == 0) {
                return {{std::uint64_t{1} << 63, -63}, {0, 0}, {0, 0}, {std::uint64_t{1} << 63, -63}};
            }
            // P = k * (2 pi 2^61), exact in 128 bits; theta = P 2^-(61 + log2n).
            const std::uint64_t p_hi = mul_hi(k, k_two_pi_q61);
            const std::uint64_t p_lo = k * k_two_pi_q61;
            // theta in Q0.64 (absolute): P >> (log2n - 3); k <= 2^(log2n - 3) keeps it below 2^64.
            const int           s       = log2n - 3;
            const std::uint64_t theta_q = s == 0 ? p_lo : (p_lo >> s) | (p_hi << (64 - s));
            // theta with relative precision: P's top 64 bits.
            wide_real theta{};
            if (p_hi != 0) {
                const int z = std::countl_zero(p_hi);
                theta       = {z == 0 ? p_hi : (p_hi << z) | (p_lo >> (64 - z)), 64 - z - 61 - log2n};
            }
            else {
                theta = normalized(p_lo, -61 - log2n);
            }
            // y = theta^2 in Q0.64 (theta <= pi/4, y < 0.62).
            const std::uint64_t y = mul_hi(theta_q, theta_q);
            // (1 - cos) / y = 1/2! - y/4! + y^2/6! - ... ; (1 - sin/theta) / y = 1/3! - y/5! + ...
            std::uint64_t tc = k_cos_coeff[10];
            for (std::size_t j = 9; j >= 1; --j) {
                tc = k_cos_coeff[j] - mul_hi(y, tc);
            }
            std::uint64_t ts = k_sin_coeff[9];
            for (std::size_t j = 8; j >= 1; --j) {
                ts = k_sin_coeff[j] - mul_hi(y, ts);
            }
            const std::uint64_t cos_q   = std::uint64_t{0} - mul_hi(y, tc); // 1 - y tc, in (0.707, 1)
            const std::uint64_t sinc_q  = std::uint64_t{0} - mul_hi(y, ts); // sin(theta)/theta, in (0.9, 1)
            const std::uint64_t sin_abs = mul_hi(theta_q, sinc_q);          // sin theta, Q0.64 absolute
            const wide_real     y_rel   = normalized(mul_hi(theta.m, theta.m), 2 * theta.e + 64);
            octant_values       r{};
            r.cos           = {cos_q, -64};
            r.sin           = normalized(mul_hi(theta.m, sinc_q), theta.e);
            r.one_minus_cos = normalized(mul_hi(y_rel.m, tc), y_rel.e);
            r.one_minus_sin = normalized(std::uint64_t{0} - sin_abs, -64);
            return r;
        }

        /// The Sample nearest m * 2^e (ties to even on the 64-bit mantissa).
        template <std::floating_point Sample>
        Sample to_sample(wide_real v) noexcept {
            static_assert(std::numeric_limits<Sample>::is_iec559, "IEEE 754 binary32 / binary64");
            if (v.m == 0) {
                return Sample(0);
            }
            constexpr int       p    = std::numeric_limits<Sample>::digits;
            std::uint64_t       keep = v.m >> (64 - p);
            const std::uint64_t rest = v.m << p;
            constexpr auto      half = std::uint64_t{1} << 63;
            if (rest > half || (rest == half && (keep & 1u) != 0)) {
                ++keep;
            }
            const int e = v.e + 64 - p;
            // 2^e as a Sample, built from its bit pattern (normal range here).
            const Sample scale = [e] {
                if constexpr (std::is_same_v<Sample, float>) {
                    return std::bit_cast<float>(static_cast<std::uint32_t>(127 + e) << 23);
                }
                else {
                    return std::bit_cast<double>(static_cast<std::uint64_t>(1023 + e) << 52);
                }
            }();
            return static_cast<Sample>(keep) * scale;
        }

        /// v / 2 (exact).
        constexpr wide_real halved(wide_real v) noexcept {
            return v.m == 0 ? v : wide_real{v.m, v.e - 1};
        }

    } // namespace srdif_trig

    /// The floating real FFT engine: float or double, N a power of two in
    /// [4, 2^30].
    ///
    /// Structure (fixed_point_rdft's, generalized). The N real samples are
    /// read as M = N/2 complex values z[n] = x[2n] + i x[2n+1]; a split-radix
    /// decimation-in-frequency kernel of length M leaves Z in bit-reversed
    /// order; the permutation puts it in natural order; the real post-pass
    /// turns Z into the packed spectrum. The inverse runs the pre-pass, the
    /// conjugate kernel and the permutation. Operation count, exactly
    /// (counted through an instrumented Sample type): the kernel's
    /// 4M log2 M - 6M + 8 real additions and multiplications (the
    /// split-radix count; M >= 2), 12 per bin pair of the post-pass
    /// (N/4 - 1 pairs) and 2 for DC / Nyquist, i.e. 2N log2 N - 2N - 2 per
    /// forward transform and two more (the halving of the DC / Nyquist pair)
    /// per inverse: 8,190 / 8,192 at N = 512, 40,958 / 40,960 at 2048.
    ///
    /// The kernel. A split-radix step on a block of l values computes, for
    /// each j < l/4, the butterfly on x[j], x[j + l/4], x[j + l/2],
    /// x[j + 3l/4]: the sums x[j] + x[j + l/2] and x[j + l/4] + x[j + 3l/4]
    /// feed the half-length block (the even outputs), and with t1 = x[j] -
    /// x[j + l/2], t2 = x[j + l/4] - x[j + 3l/4] the values (t1 + i t2) W^j
    /// and (t1 - i t2) W^3j feed the two quarter-length blocks (outputs
    /// 4k + 1 and 4k + 3), W = exp(+2 pi i/l) (the conjugates, and -i for +i,
    /// in the inverse). Recursing on [0, l/2), [l/2, 3l/4) and [3l/4, l)
    /// leaves the outputs where a radix-2 decimation in frequency leaves
    /// them: X[bitrev(p)] at position p. Arranged here as:
    ///   - The fused pass (fused_pass, block<32 | 64>): the half-length
    ///     block's own butterfly at j reads exactly the half-length outputs
    ///     of the level-l butterflies at j and j + l/8 (its inputs sit at j,
    ///     j + l/8, j + l/4, j + 3l/8), so one group of eight values
    ///     x[j + r l/8], r = 0 … 7, runs both levels between one load and one
    ///     store: three butterflies per 16 loads and 16 stores instead of
    ///     48 of each. Every block of l >= 32 is processed this way and
    ///     recurses on the five blocks the two levels leave: l/4 (the first
    ///     half's half), l/8 twice (its quarters) and l/4 twice (the
    ///     quarters). Groups run two per loop step in the run-time pass
    ///     (l >= 128) so its eight row pointers serve both.
    ///   - Special butterflies, as Sorensen, Heideman and Burrus count them:
    ///     j = 0 (W = 1, additions only) and j = l/8 (W^j = (1 + i)/sqrt 2,
    ///     W^3j = (-1 + i)/sqrt 2: two additions and two multiplications per
    ///     product instead of two and four); in a group, the level-l/2
    ///     butterfly at j = l/16 is its eighth-turn one, and its level-l
    ///     twiddles are W_16^1, W_16^3, W_16^3 and W_16^9, spelled as
    ///     literals (KernelLiteralsAreTheGeneratorsValues).
    ///   - Twiddle pairing: the groups j and l/8 - j use the same six
    ///     twiddles, by W^(l/4 - j) = i conj W^j and W^3(l/4 - j) =
    ///     -i conj W^3j (an exchange of cos and sin with sign flips, exact),
    ///     so the table holds j < l/16 only and the second half of the pass
    ///     walks it backwards (group_high).
    ///   - Compile-time blocks: blocks of 32 and 64 (block<32>, block<64>)
    ///     run with every offset a constant and their own small tables, and
    ///     blocks of 16 and fewer are register leaves (leaf<L>: loaded,
    ///     transformed through every remaining level and stored once). The
    ///     run-time recursion (kernel) handles l >= 128 only.
    ///
    /// The tables (built in the constructor, never touched by a transform
    /// except to read):
    ///   - m_twiddles: for the fused pass at run time, twelve Samples per
    ///     jj in [1, M/16): W_M^jj, W_M^3jj, W_M^(jj + M/8),
    ///     W_M^(3jj + 3M/8), W_M^2jj, W_M^6jj (cos, sin each); a block of
    ///     length l reads the entry j (M/l) for its group j. 12 (M/16 - 1)
    ///     Samples (none below M = 32).
    ///   - m_small: the same layout at M = 32 (one entry) and 64 (three), for
    ///     block<32> and block<64>: 12 + 36 Samples where used.
    ///   - m_post: C_k = (1 + i W_N^k) / 2 for k < N/4 (fixed_point.h's
    ///     real_post_pass derives it), stored as ((1 - sin theta_k) / 2,
    ///     cos theta_k / 2) with 1 - sin and, past the octant, 1 - cos
    ///     evaluated without cancellation by srdif_trig, each rounded once:
    ///     N/2 Samples.
    ///   - m_swaps: the bit-reversal permutation as the list of pairs
    ///     (i, bitrev i), i < bitrev i, as Sample offsets: M - 2^ceil(log2 M
    ///     / 2) uint32 words.
    ///   Every table value is srdif_trig's, taken to the other octants by
    ///   exact symmetry (negation and exchange). Heap per object, measured
    ///   (x86-64, libstdc++): 2,896 B float / 4,832 B double at N = 512,
    ///   11,280 B / 18,592 B at 2048, 359,568 B / 589,088 B at 65536, 2.6 to
    ///   2.8 times what the engine it replaced held (docs/fft-design.md has
    ///   the comparison); sizeof 104 B on LP64.
    ///
    /// Contract numbers (fft.h re-exports the first three):
    ///   - k_min_size = 4, k_max_size = 2^30 (the uint32 swap offsets reach
    ///     2^30 - 2 there; at 2^30 the float tables alone are 5.5 GiB beside
    ///     the caller's 4 GiB, and no test constructs it).
    ///   - k_is_shareable = true: the transforms are const, read the tables
    ///     and write the caller's buffer only, so two threads may transform
    ///     through one object at once.
    ///   - Transforms noexcept and allocation-free (tests/test_fft_rt.cpp),
    ///     the first as cheap as every later one; no alignment requirement;
    ///     no data-dependent branch, so NaN propagates like any value;
    ///     copyable, and a copy is bit-identical to its source.
    ///   - Output bits: a function of the input, of this source and of
    ///     fp-contraction alone (srdif_trig). Built without contraction they
    ///     are the same on every host and QEMU leg, pinned at every power of
    ///     two 4 … 65536 (tests/test_fft_srdif_fingerprint.cpp). A build that
    ///     fuses a*b + c (GCC's default on FMA targets such as the M4F, M33
    ///     and M55 legs; clang's, within a statement, on Apple arm64) moves
    ///     last bits and not the error statistics (fft.h, "fp-contraction
    ///     policy").
    ///
    /// Accuracy (the method of the maintainer's targets sheet: against a
    /// __float128 reference, white noise at full scale, eight trials pooled;
    /// x86-64, GCC 13.3 -O3, no FMA), relative rms error of the forward /
    /// unnormalized inverse / round trip:
    ///   float  N = 512   1.06e-7 / 1.02e-7 / 1.44e-7   (1.8 u / 1.7 u / 2.4 u)
    ///          N = 2048  1.17e-7 / 1.16e-7 / 1.62e-7
    ///          N = 65536 1.43e-7 / 1.42e-7 / 2.00e-7
    ///   double N = 512   1.89e-16 / 1.89e-16 / 2.61e-16
    ///          N = 2048  2.18e-16 / 2.17e-16 / 3.00e-16
    ///          N = 65536 2.66e-16 / 2.65e-16 / 3.69e-16
    /// (u the unit roundoff): 4 to 14 % below the engine it replaced at every
    /// N >= 128 in both profiles and every direction, and at most 4 % above
    /// it below that (float round trip at N = 64; docs/fft-design.md has the
    /// full tables, the per-signal cells and the maxima). Against the
    /// battery's own pins: double forward vs the compensated DFT at N = 256
    /// 1.60e-16 (pin 7.4e-16), float vs double at N = 512 9.73e-8 (pin
    /// 2.25e-7), the oracle's worst ratio to 4 eps log2 N ||y||_2 0.155
    /// (float) and 0.25 (double).
    ///
    /// Speed: the instruction-count ratchet's float scenarios (forward plus
    /// scaled inverse over 2^20 samples, harness included; arm-none-eabi-gcc
    /// 13.2.1 -O3, QEMU 8.2.2), against the engine it replaced, N = 512 /
    /// 2048: Cortex-M4 soft-float -2.67 % / -2.18 %, M4F -3.23 % / -2.92 %,
    /// M33 -6.50 % / -6.06 %, M55 without CMSIS -3.46 % / -2.62 %
    /// (bench/README.md, docs/fft-design.md, "The floating engine (srdif)").
    /// On x86-64 without -march it is 3-10 % slower in float and 6-20 % in
    /// double than the engine it replaced (informational; the design note
    /// has the table).
    template <std::floating_point Sample>
    class srdif_rdft {
        static_assert(std::is_same_v<Sample, float> || std::is_same_v<Sample, double>, "srdif_rdft: float or double");

      public:
        static constexpr std::size_t k_min_size     = 4;
        static constexpr std::size_t k_max_size     = std::size_t{1} << 30;
        static constexpr bool        k_is_shareable = true;

        explicit srdif_rdft(std::size_t n)
            : m_n(n)
            , m_twiddles(make_kernel_table(n / 2))
            , m_small(make_small_tables(n / 2))
            , m_post(make_post_table(n))
            , m_swaps(make_swap_table(n / 2)) {
            assert(n >= k_min_size && n <= k_max_size && (n & (n - 1)) == 0);
        }

        [[nodiscard]] std::size_t size() const noexcept { return m_n; }

        TAP_DSP_SRDIF_NOINLINE void forward_inplace(Sample* a) const noexcept {
            kernel<false>(a, m_n / 2);
            permute(a);
            post_pass<false>(a);
        }

        TAP_DSP_SRDIF_NOINLINE void inverse_inplace(Sample* a) const noexcept {
            post_pass<true>(a);
            kernel<true>(a, m_n / 2);
            permute(a);
        }

      private:
        static constexpr Sample k_sqrt_half =
            std::is_same_v<Sample, float> ? Sample(0x1.6a09e6p-1) : Sample(0x1.6a09e667f3bcdp-1);

        static int log2_of(std::size_t n) noexcept { return std::bit_width(n) - 1; }

        /// The fused pass's twiddles at resolution m: for jj in [1, m/16), twelve
        /// Samples W^jj, W^3jj, W^(jj + m/8), W^(3jj + 3m/8), W^2jj, W^6jj
        /// (W = exp(+2 pi i/m), each as cos, sin). Entry 0 is not stored.
        static std::vector<Sample> make_kernel_table(std::size_t m) {
            const std::size_t   count = m >= 32 ? m / 16 - 1 : 0; // m == 0: empty
            std::vector<Sample> t(12 * count);
            if (count == 0) {
                return t;
            }
            // The first octant, b in [0, m/8], once; every other angle by exact symmetry.
            const std::size_t   octant = m / 8;
            std::vector<Sample> fo(2 * (octant + 1));
            const int           lg = log2_of(m);
            for (std::size_t b = 0; b <= octant; ++b) {
                const auto v  = srdif_trig::first_octant(b, lg);
                fo[2 * b]     = srdif_trig::to_sample<Sample>(v.cos);
                fo[2 * b + 1] = srdif_trig::to_sample<Sample>(v.sin);
            }
            const std::size_t quarter = m / 4;
            const auto        unit    = [&](std::size_t angle, Sample* out) {
                const std::size_t a  = angle % m;
                const std::size_t qd = a / quarter;
                const std::size_t r  = a % quarter;
                Sample            c;
                Sample            sn;
                if (r <= octant) {
                    c  = fo[2 * r];
                    sn = fo[2 * r + 1];
                }
                else { // W^r = i conj W^(m/4 - r)
                    c  = fo[2 * (quarter - r) + 1];
                    sn = fo[2 * (quarter - r)];
                }
                for (std::size_t k = 0; k < qd; ++k) { // times i per quarter turn
                    const Sample t0 = c;
                    c               = -sn;
                    sn              = t0;
                }
                out[0] = c;
                out[1] = sn;
            };
            for (std::size_t jj = 1; jj <= count; ++jj) {
                Sample* const e = t.data() + 12 * (jj - 1);
                unit(jj, e);
                unit(3 * jj, e + 2);
                unit(jj + m / 8, e + 4);
                unit(3 * jj + 3 * (m / 8), e + 6);
                unit(2 * jj, e + 8);
                unit(6 * jj, e + 10);
            }
            return t;
        }

        /// The fused pass's twiddles for the blocks of 32 and 64 (the kernel
        /// tables of those lengths, concatenated: 12 + 36 Samples), read by the
        /// compile-time blocks with constant offsets.
        static std::vector<Sample> make_small_tables(std::size_t m) {
            std::vector<Sample> t = make_kernel_table(m >= 32 ? 32 : 0);
            if (m >= 64) {
                const std::vector<Sample> t64 = make_kernel_table(64);
                t.insert(t.end(), t64.begin(), t64.end());
            }
            return t;
        }

        /// C_k = (1 + i W_n^k) / 2 for k in [0, n/4), interleaved (re, im).
        static std::vector<Sample> make_post_table(std::size_t n) {
            const std::size_t   quarter = n / 4;
            const std::size_t   octant  = n / 8;
            std::vector<Sample> t(2 * quarter);
            if (quarter == 0) {
                return t;
            }
            const int lg = log2_of(n);
            t[0]         = Sample(0.5);
            t[1]         = Sample(0.5);
            for (std::size_t b = 1; b <= octant; ++b) {
                const auto v = srdif_trig::first_octant(b, lg);
                // k = b: re = (1 - sin) / 2, im = cos / 2.
                t[2 * b]     = srdif_trig::to_sample<Sample>(srdif_trig::halved(v.one_minus_sin));
                t[2 * b + 1] = srdif_trig::to_sample<Sample>(srdif_trig::halved(v.cos));
                // k = n/4 - b: sin theta_k = cos theta_b, cos theta_k = sin theta_b.
                const std::size_t k = quarter - b;
                if (k != b && k > 0) {
                    t[2 * k]     = srdif_trig::to_sample<Sample>(srdif_trig::halved(v.one_minus_cos));
                    t[2 * k + 1] = srdif_trig::to_sample<Sample>(srdif_trig::halved(v.sin));
                }
            }
            return t;
        }

        /// The bit-reversal permutation of m complex values as swap pairs (i, j), i < j.
        static std::vector<std::uint32_t> make_swap_table(std::size_t m) {
            // m - 2^ceil(bits/2) indices are not their own reversal: half as many pairs.
            const int                  bits = log2_of(m);
            std::vector<std::uint32_t> t;
            t.reserve(m - (std::size_t{1} << ((bits + 1) / 2)));
            for (std::size_t i = 0; i < m; ++i) {
                std::size_t r = 0;
                for (int b = 0; b < bits; ++b) {
                    r |= ((i >> b) & 1u) << (bits - 1 - b);
                }
                if (i < r) { // Sample offsets of the two complex values
                    t.push_back(static_cast<std::uint32_t>(2 * i));
                    t.push_back(static_cast<std::uint32_t>(2 * r));
                }
            }
            return t;
        }

        void permute(Sample* a) const noexcept {
            const std::uint32_t* s   = m_swaps.data();
            const std::uint32_t* end = s + m_swaps.size();
            for (; s != end; s += 2) {
                Sample* const p  = a + s[0];
                Sample* const q  = a + s[1];
                const Sample  p0 = p[0];
                const Sample  p1 = p[1];
                const Sample  q0 = q[0];
                const Sample  q1 = q[1];
                p[0]             = q0;
                p[1]             = q1;
                q[0]             = p0;
                q[1]             = p1;
            }
        }

        static constexpr Sample k_cos_pi8 =
            std::is_same_v<Sample, float> ? Sample(0x1.d906bcp-1) : Sample(0x1.d906bcf328d46p-1);
        static constexpr Sample k_sin_pi8 =
            std::is_same_v<Sample, float> ? Sample(0x1.87de2ap-2) : Sample(0x1.87de2a6aea963p-2);

        /// A complex value in registers.
        struct cv {
            Sample r;
            Sample i;
        };

        static cv   ld(const Sample* p) noexcept { return {p[0], p[1]}; }
        static void st(Sample* p, cv v) noexcept {
            p[0] = v.r;
            p[1] = v.i;
        }

        /// The three kinds of split-radix butterfly, by twiddle.
        enum class tw { one, eighth, general };

        /// One split-radix DIF butterfly on the four values of a length-l block
        /// at j, j + l/4, j + l/2, j + 3l/4, in place: x0 <- x0 + x2,
        /// x1 <- x1 + x3 (the half-length block's inputs), x2 <- (t1 + i t2) W^j,
        /// x3 <- (t1 - i t2) W^3j with t1 = x0 - x2, t2 = x1 - x3 (forward; the
        /// inverse conjugates i and both twiddles). Kind one: W^j = W^3j = 1;
        /// eighth: j = l/8, W^j = (1 + i)/sqrt 2, W^3j = (-1 + i)/sqrt 2;
        /// general: w1 = W^j, w3 = W^3j as given.
        template <bool Inverse, tw Kind>
        static void bf(cv& x0, cv& x1, cv& x2, cv& x3, cv w1 = {}, cv w3 = {}) noexcept {
            const Sample t1r = x0.r - x2.r;
            const Sample t1i = x0.i - x2.i;
            const Sample t2r = x1.r - x3.r;
            const Sample t2i = x1.i - x3.i;
            x0               = {x0.r + x2.r, x0.i + x2.i};
            x1               = {x1.r + x3.r, x1.i + x3.i};
            // u = t1 + i t2, v = t1 - i t2 (forward); swapped for the inverse.
            Sample ur;
            Sample ui;
            Sample vr;
            Sample vi;
            if constexpr (Inverse) {
                ur = t1r + t2i;
                ui = t1i - t2r;
                vr = t1r - t2i;
                vi = t1i + t2r;
            }
            else {
                ur = t1r - t2i;
                ui = t1i + t2r;
                vr = t1r + t2i;
                vi = t1i - t2r;
            }
            if constexpr (Kind == tw::one) {
                (void)w1; // read by the general kind only
                (void)w3;
                x2 = {ur, ui};
                x3 = {vr, vi};
            }
            else if constexpr (Kind == tw::eighth) {
                (void)w1;
                (void)w3;
                const Sample h = k_sqrt_half;
                if constexpr (Inverse) { // u (1 - i) h, v (-1 - i) h
                    x2 = {h * (ur + ui), h * (ui - ur)};
                    x3 = {h * (vi - vr), -h * (vr + vi)};
                }
                else { // u (1 + i) h, v (-1 + i) h
                    x2 = {h * (ur - ui), h * (ur + ui)};
                    x3 = {-h * (vr + vi), h * (vr - vi)};
                }
            }
            else {
                if constexpr (Inverse) {
                    x2 = {ur * w1.r + ui * w1.i, ui * w1.r - ur * w1.i};
                    x3 = {vr * w3.r + vi * w3.i, vi * w3.r - vr * w3.i};
                }
                else {
                    x2 = {ur * w1.r - ui * w1.i, ur * w1.i + ui * w1.r};
                    x3 = {vr * w3.r - vi * w3.i, vr * w3.i + vi * w3.r};
                }
            }
        }

        /// One group of the fused pass: the eight values at p + r e (r = 0..7,
        /// complex units) of a length-l block with e = l/8. The level-l
        /// butterflies at j and j + e (twiddles a1, a3 and b1, b3), then the
        /// level-l/2 butterfly at j of the first half (c1, c3), whose four
        /// inputs are exactly the first-half outputs of the other two.
        template <bool Inverse, tw KA, tw KB, tw KC>
        static void group(Sample* p, std::size_t e, cv a1, cv a3, cv b1, cv b3, cv c1, cv c3) noexcept {
            const std::size_t d  = 2 * e;
            cv                x0 = ld(p);
            cv                x2 = ld(p + 2 * d);
            cv                x4 = ld(p + 4 * d);
            cv                x6 = ld(p + 6 * d);
            bf<Inverse, KA>(x0, x2, x4, x6, a1, a3);
            st(p + 4 * d, x4);
            st(p + 6 * d, x6);
            cv x1 = ld(p + d);
            cv x3 = ld(p + 3 * d);
            cv x5 = ld(p + 5 * d);
            cv x7 = ld(p + 7 * d);
            bf<Inverse, KB>(x1, x3, x5, x7, b1, b3);
            st(p + 5 * d, x5);
            st(p + 7 * d, x7);
            bf<Inverse, KC>(x0, x1, x2, x3, c1, c3);
            st(p, x0);
            st(p + d, x1);
            st(p + 2 * d, x2);
            st(p + 3 * d, x3);
        }

        /// The four kinds of group in a fused pass over a block of length l = 8e.
        template <bool Inverse>
        static void group_first(Sample* a, std::size_t e) noexcept {
            // j = 0: W = 1, then the eighth-turn at j + e, then W = 1 at level l/2.
            group<Inverse, tw::one, tw::eighth, tw::one>(a, e, {}, {}, {}, {}, {}, {});
        }

        /// j in [1, e/2), twelve twiddles at w.
        template <bool Inverse>
        static void group_low(Sample* a, std::size_t e, const Sample* w) noexcept {
            group<Inverse, tw::general, tw::general, tw::general>(a, e, {w[0], w[1]}, {w[2], w[3]}, {w[4], w[5]},
                                                                  {w[6], w[7]}, {w[8], w[9]}, {w[10], w[11]});
        }

        /// j = e/2 = l/16: W_16^1, W_16^3; W_16^3, W_16^9 = -W_16^1; the eighth-turn at level l/2.
        template <bool Inverse>
        static void group_mid(Sample* a, std::size_t e) noexcept {
            const cv c16{k_cos_pi8, k_sin_pi8};
            const cv s16{k_sin_pi8, k_cos_pi8};
            group<Inverse, tw::general, tw::general, tw::eighth>(a, e, c16, s16, s16, {-k_cos_pi8, -k_sin_pi8}, {}, {});
        }

        /// j in (e/2, e), reading the entry of j' = e - j: W^j = i conj W^(j'+e),
        /// W^3j = -i conj W^3(j'+e), W^(j+e) = i conj W^j', W^3(j+e) = -i conj W^3j',
        /// and at level l/2 W^j = i conj W^j', W^3j = -i conj W^3j'.
        template <bool Inverse>
        static void group_high(Sample* a, std::size_t e, const Sample* w) noexcept {
            group<Inverse, tw::general, tw::general, tw::general>(a, e, {w[5], w[4]}, {-w[7], -w[6]}, {w[1], w[0]},
                                                                  {-w[3], -w[2]}, {w[9], w[8]}, {-w[11], -w[10]});
        }

        /// The fused pass over a block of l >= 128 complex values: level l and
        /// level l/2 of its first half; the table entry of j is j (m/l) steps.
        template <bool Inverse>
        void fused_pass(Sample* a, std::size_t l) const noexcept {
            const std::size_t e      = l / 8;
            const std::size_t h      = e / 2;
            const std::size_t stride = 12 * ((m_n / 2) / l); // Samples per table step at this length
            group_first<Inverse>(a, e);
            // j = 1 alone, then two groups per step (e/2 - 1 is odd): the
            // eight row pointers serve both groups at constant offsets.
            const Sample* w = m_twiddles.data() + stride - 12;
            group_low<Inverse>(a + 2, e, w);
            w += stride;
            for (std::size_t j = 2; j < h; j += 2, w += 2 * stride) {
                group_low<Inverse>(a + 2 * j, e, w);
                group_low<Inverse>(a + 2 * j + 2, e, w + stride);
            }
            group_mid<Inverse>(a + 2 * h, e);
            // j in (e/2, e) descending through the table: j' = e - j.
            w -= stride;
            group_high<Inverse>(a + 2 * (h + 1), e, w);
            for (std::size_t j = h + 2; j < e; j += 2) {
                w -= 2 * stride;
                group_high<Inverse>(a + 2 * j, e, w + stride);
                group_high<Inverse>(a + 2 * j + 2, e, w);
            }
        }

        /// A block of L <= 64 complex values with everything known at compile
        /// time: the fused pass (L = 32, 64) over the small tables, or a leaf.
        template <bool Inverse, std::size_t L>
        static void block(Sample* a, const Sample* small) noexcept {
            if constexpr (L >= 32) {
                constexpr std::size_t e = L / 8;
                constexpr std::size_t h = e / 2;
                const Sample* const   w = small + (L == 32 ? 0 : 12); // this length's table
                group_first<Inverse>(a, e);
                [&]<std::size_t... J>(std::index_sequence<J...>) {
                    (group_low<Inverse>(a + 2 * (J + 1), e, w + 12 * J), ...);
                    (group_high<Inverse>(a + 2 * (e - 1 - J), e, w + 12 * J), ...);
                }(std::make_index_sequence<h - 1>{});
                group_mid<Inverse>(a + 2 * h, e);
                block<Inverse, L / 4>(a, small);
                block<Inverse, L / 8>(a + L / 2, small);
                block<Inverse, L / 8>(a + L / 2 + L / 4, small);
                block<Inverse, L / 4>(a + L, small);
                block<Inverse, L / 4>(a + L + L / 2, small);
            }
            else {
                (void)small; // the leaves need no table
                leaf<Inverse, L>(a);
            }
        }

        /// The split-radix recursion over a block held in registers (a leaf of
        /// 2 ... 16 complex values), block by block.
        template <bool Inverse, std::size_t L>
        static void leaf_blocks(cv* v) noexcept {
            if constexpr (L == 1) {
                (void)v; // a block of one value is its own transform
            }
            else if constexpr (L == 2) {
                const cv x0 = v[0];
                const cv x1 = v[1];
                v[0]        = {x0.r + x1.r, x0.i + x1.i};
                v[1]        = {x0.r - x1.r, x0.i - x1.i};
            }
            else if constexpr (L >= 4) {
                static_assert(L <= 16, "leaves are 2 ... 16 complex values");
                constexpr std::size_t q = L / 4;
                bf<Inverse, tw::one>(v[0], v[q], v[2 * q], v[3 * q]);
                if constexpr (L >= 8) {
                    constexpr std::size_t e = L / 8;
                    bf<Inverse, tw::eighth>(v[e], v[q + e], v[2 * q + e], v[3 * q + e]);
                }
                if constexpr (L == 16) { // j = 1: W16^1, W16^3; j = 3: W16^3, W16^9 = -W16^1
                    const cv c16{k_cos_pi8, k_sin_pi8};
                    const cv s16{k_sin_pi8, k_cos_pi8};
                    bf<Inverse, tw::general>(v[1], v[5], v[9], v[13], c16, s16);
                    bf<Inverse, tw::general>(v[3], v[7], v[11], v[15], s16, {-k_cos_pi8, -k_sin_pi8});
                }
                leaf_blocks<Inverse, L / 4>(v + L / 2);
                leaf_blocks<Inverse, L / 4>(v + 3 * q);
                leaf_blocks<Inverse, L / 2>(v);
            }
        }

        template <bool Inverse, std::size_t L>
        static void leaf(Sample* a) noexcept {
            [&]<std::size_t... I>(std::index_sequence<I...>) {
                cv v[L] = {ld(a + 2 * I)...};
                leaf_blocks<Inverse, L>(v);
                (st(a + 2 * I, v[I]), ...);
            }(std::make_index_sequence<L>{});
        }

        /// Split-radix DIF over l complex values at a (bit-reversed output).
        template <bool Inverse>
        void kernel(Sample* a, std::size_t l) const noexcept {
            const Sample* const small = m_small.data();
            switch (l) {
            case 64:
                block<Inverse, 64>(a, small);
                return;
            case 32:
                block<Inverse, 32>(a, small);
                return;
            case 16:
                block<Inverse, 16>(a, small);
                return;
            case 8:
                block<Inverse, 8>(a, small);
                return;
            case 4:
                block<Inverse, 4>(a, small);
                return;
            case 2:
                block<Inverse, 2>(a, small);
                return;
            default:
                break;
            }
            if (l >= 128) {
                fused_pass<Inverse>(a, l);
                kernel<Inverse>(a, l / 4);         // the first half's half
                kernel<Inverse>(a + l / 2, l / 8); // the first half's quarters
                kernel<Inverse>(a + l / 2 + l / 4, l / 8);
                kernel<Inverse>(a + l, l / 4); // the quarters
                kernel<Inverse>(a + l + l / 2, l / 4);
            }
        }

        /// The real post-pass (forward) / pre-pass (Inverse); fixed_point.h derives it.
        template <bool Inverse>
        void post_pass(Sample* a) const noexcept {
            const std::size_t m = m_n / 2;
            const Sample      r = a[0];
            const Sample      i = a[1];
            if constexpr (Inverse) {
                a[0] = Sample(0.5) * (r + i);
                a[1] = Sample(0.5) * (r - i);
            }
            else {
                a[0] = r + i;
                a[1] = r - i;
            }
            if (m < 4) {
                return;
            }
            // Pairs (k, m - k) for k in [1, m/2): k = 1 alone, then two per step
            // (m/2 - 1 is odd for every m >= 4).
            Sample*       pk  = a + 2;
            Sample*       pj  = a + 2 * (m - 1);
            const Sample* c   = m_post.data() + 2;
            Sample* const end = a + m;
            post_pair<Inverse>(pk, pj, c);
            for (pk += 2, pj -= 2, c += 2; pk < end; pk += 4, pj -= 4, c += 4) {
                post_pair<Inverse>(pk, pj, c);
                post_pair<Inverse>(pk + 2, pj - 2, c + 2);
            }
        }

        /// One bin pair of the post-pass: u = a[k], v = conj a[m - k],
        /// G = C_k (u - v) (conj C_k for the inverse), a[k] <- u - G,
        /// a[m - k] <- conj(v + G).
        template <bool Inverse>
        static void post_pair(Sample* pk, Sample* pj, const Sample* c) noexcept {
            const Sample ur = pk[0];
            const Sample ui = pk[1];
            const Sample vr = pj[0];
            const Sample vi = pj[1];
            const Sample dr = ur - vr;
            const Sample di = ui + vi;
            const Sample cr = c[0];
            const Sample ci = c[1];
            Sample       gr;
            Sample       gi;
            if constexpr (Inverse) {
                gr = dr * cr + di * ci;
                gi = di * cr - dr * ci;
            }
            else {
                gr = dr * cr - di * ci;
                gi = dr * ci + di * cr;
            }
            pk[0] = ur - gr;
            pk[1] = ui - gi;
            pj[0] = vr + gr;
            pj[1] = vi - gi;
        }

        std::size_t                m_n;
        std::vector<Sample>        m_twiddles;
        std::vector<Sample>        m_small;
        std::vector<Sample>        m_post;
        std::vector<std::uint32_t> m_swaps;
    };

} // namespace tap::dsp::detail
