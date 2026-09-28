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
// DspTap's own fixed-point engine (fft/fixed_point.h), whose real post-pass
// tap/DspTap#39 records as the same arrangement of the published method as
// the third-party package's, and whose structure it generalizes to floating
// point (docs/fft-design.md, "The floating engine (srdif)"; that this engine
// is not derived from the package is the maintainer's judgement, NOTICE.md).
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
//   The permutation: B. Gold and C. M. Rader, Digital Processing of
//     Signals, McGraw-Hill, 1969 (the reverse-carry counter: a bit-reversed
//     index advanced beside the natural one, swapping when i < rev(i));
//     A. H. Karp, "Bit reversal on uniprocessors," SIAM Review 38(1), 1-26,
//     1996 (the survey of such methods, among them splitting the index
//     into its end bits and its middle so that one counter serves several
//     swaps).
//
//   The tables: the Maclaurin series of cos and sin evaluated by Horner's
//     rule in 64-bit fixed point (D. E. Knuth, The Art of Computer
//     Programming vol. 2, 3rd ed., Sec. 4.3.1 for the multiple-precision
//     products and Sec. 4.6.4 for Horner's rule), and the exact symmetries
//     of the unit circle for everything outside the first octant, as
//     tables.h does for the fixed-point tables.
//
// The arrangement of these pieces (the fused two-level pass, the twiddle
// pairing, the compile-time blocks and their in-memory levels, the table layouts)
// is derived in the docstrings below.

#pragma once

#include <array>
#include <bit>
#include <cassert>
#include <concepts>
#include <cstddef>
#include <cstdint>
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

// Under clang, the kernel's building blocks (butterflies, groups, the
// compile-time blocks and their levels, the post-pass pairs, the permutation's
// swaps) are forced inline into the pass that calls them (a performance
// attribute only; the arithmetic and its order are the same either way).
// clang's inliner, left to itself, kept them out of line on Hexagon (clang
// 19, v68): the register leaves this engine had then went through a stack
// array and every group through a call. Forcing them was the largest single
// step of the Hexagon tuning (rfft_f32_512 57.66 -> 49.67 M instructions
// with the leaves of the time; the tree as it stands reads 54.13 M without
// the attribute and 45.94 M with it: docs/fft-design.md, "Hexagon (clang)
// tuning"). GCC is left to its own choices: forced there, the Cortex-M
// keys with an FPU (M4F, M33, M55 without CMSIS) cost 4.5 - 9.1 % more
// instructions (spills in the larger bodies; the soft-float M4 0.1 - 0.2 %
// less). Not forced when optimizing for size (-Os / -Oz define
// __OPTIMIZE_SIZE__): forced, the Hexagon float size probe's .text grows
// from 238,628 to 263,204 bytes.
#if defined(__clang__) && !defined(__OPTIMIZE_SIZE__)
#define TAP_DSP_SRDIF_INLINE __attribute__((always_inline))
#else
#define TAP_DSP_SRDIF_INLINE
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
    /// Accuracy, against __float128 (libquadmath), exhaustive over every k
    /// of every 2^L, L = 3 … 30 (every angle any N <= 2^30 builds; review A
    /// of tap/DspTap#42, 2026-09-26): the 64-bit mantissa is within 13.44
    /// units of its last place of the exact value (1 - cos is the worst of
    /// the four, at L = 30; 12.34 at L = 19, 13.07 at L = 24). Rounded by
    /// to_sample: every float value is the correctly rounded one, with no
    /// misrounding at any L <= 30 (the nearest exact value to a float
    /// rounding tie is 547 units away, first at L = 25; 21,970 units up to
    /// L = 24); a double value is within 0.5 + 13.44/2048 = 0.5 + 2^-7.25 ulp
    /// of the exact one (0.5052 ulp measured), and 0.04 - 0.17 % of them,
    /// by function, are the neighbour of the correctly rounded double
    /// (0.08 % of the 2^20 post table). Deterministic either way.
    /// `GeneratorAgreesWithLibmToTheLastBits`, `OctantEndpointsAreExact`,
    /// `KernelLiteralsAreTheGeneratorsValues`, `ToSampleRoundsHalfToEven`
    /// (tests/test_fft_srdif.cpp).
    ///
    /// Cost: about 25 64 x 64 -> 128 high products per angle (four 32-bit
    /// multiplies each); the engine evaluates N/8 angles at construction
    /// (the post table's; the kernel tables are read out of it). Measured on
    /// the Cortex-M4F (QEMU, -O3), construction whole, allocation included:
    /// 52 k instructions at N = 512, 204 k at N = 2048 (0.06 % / 0.19 % of
    /// the ratchet scenarios).
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
    ///     quarters). In float, groups run two per loop step in the
    ///     run-time pass (l >= 128) so its eight row pointers serve both; in
    ///     double one (k_one_group_per_step).
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
    ///     blocks of 16 and fewer run one split-radix level at a time in
    ///     memory (level<L>, one butterfly at a time) down to blocks of two.
    ///     The run-time recursion (kernel) handles l >= 128 only.
    ///
    /// The tables: one allocation, sized once and built in place by the
    /// constructor, never touched by a transform except to read. In Samples,
    /// in this order:
    ///   - post: C_k = (1 + i W_N^k) / 2 for k < N/4 (fixed_point.h's
    ///     real_post_pass derives it), stored as ((1 - sin theta_k) / 2,
    ///     cos theta_k / 2) with 1 - sin and, past the octant, 1 - cos
    ///     evaluated without cancellation by srdif_trig, each rounded once:
    ///     N/2 Samples.
    ///   - small: the fused pass's twiddles at M = 32 (one entry) and 64
    ///     (three), for block<32> and block<64>: 12 Samples at M = 32, 12 + 36
    ///     from M = 64.
    ///   - kernel: for the run-time fused pass (M >= 128), twelve Samples per
    ///     jj in [1, M/16): W_M^jj, W_M^3jj, W_M^(jj + M/8),
    ///     W_M^(3jj + 3M/8), W_M^2jj, W_M^6jj (cos, sin each); a block of
    ///     length l reads the entry j (M/l) for its group j. 12 (M/16 - 1)
    ///     Samples.
    ///   The post table is srdif_trig's; small and kernel are read out of
    ///   its imaginary parts (fill_kernel_table), which gives the same bits
    ///   as evaluating srdif_trig at their own resolution. Every value
    ///   outside the first octant is taken there by exact symmetry (negation
    ///   and exchange). The permutation needs no table (permute).
    ///   Heap per object, exactly (heap_bytes; one allocation, no temporary
    ///   at construction): sizeof(Sample) x (N/2 for N <= 32, 44 at N = 64,
    ///   112 at N = 128, 7N/8 + 36 from N = 256). Float 1,936 B at N = 512,
    ///   7,312 B at 2048, 229,520 B at 65536, double twice that: 1.73 - 1.82
    ///   times what the engine it replaced held (docs/fft-design.md has the
    ///   comparison). sizeof 32 B on LP64.
    ///
    /// Contract numbers (fft.h re-exports the first three):
    ///   - k_min_size = 4, k_max_size = 2^30 (the contract's range, kept from
    ///     the engine this one replaced; every index is a size_t, and at 2^30
    ///     the float tables are 3.5 GiB beside the caller's 4 GiB). The
    ///     constructor checks the precondition (assert) before any table is
    ///     sized from n.
    ///   - k_is_shareable = true: the transforms are const, read the tables
    ///     and write the caller's buffer only, so two threads may transform
    ///     through one object at once.
    ///   - Transforms noexcept and allocation-free (tests/test_fft_rt.cpp),
    ///     the first as cheap as every later one; no alignment requirement;
    ///     no data-dependent branch, so a NaN or Inf in the input reaches
    ///     every output bin (per bin, not per word: one component of bin N/4,
    ///     whose weight is structurally zero, can stay finite, as it did in
    ///     the engine this one replaced; review A of tap/DspTap#42); copyable,
    ///     and a copy is bit-identical to its source.
    ///   - Output bits: a function of the input, of this source and of
    ///     fp-contraction alone (srdif_trig), given arithmetic that rounds
    ///     every operation to its format (FLT_EVAL_METHOD == 0: not an x87
    ///     build) and no flush-to-zero or -ffast-math. Built without
    ///     contraction they are one row for every compiler and target CI
    ///     runs (x86-64 g++ and clang++, Windows x64 MSVC, macOS arm64
    ///     AppleClang, the four QEMU legs; MSVC on ARM64 has no leg), pinned
    ///     at every power of two 4 … 65536
    ///     (tests/test_fft_srdif_fingerprint.cpp). A build that fuses a*b + c
    ///     (GCC's default on FMA targets such as the M4F, M33 and M55 legs;
    ///     clang's, within a statement, on Apple arm64) moves last bits and
    ///     not the error statistics (fft.h, "fp-contraction policy").
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
    /// 2048: Cortex-M4 soft-float -2.78 % / -2.28 %, M4F -5.76 % / -5.51 %,
    /// M33 -8.24 % / -7.98 %, M55 without CMSIS -6.33 % / -5.52 % (at #42,
    /// before the Hexagon tuning: -2.71 / -2.23, -4.70 / -4.49, -7.89 /
    /// -7.53, -5.98 / -5.04 %); Hexagon (clang 19.1.5 -O3, -mv68 -mhvx,
    /// qemu-hexagon 8.2.2, which counts packets) float -7.64 % / -5.88 %,
    /// double at N = 512 -3.20 %, where #42 read +15.9 % / +19.8 % and
    /// +7.2 % (bench/README.md; docs/fft-design.md, "The floating engine
    /// (srdif)" and its "Hexagon (clang) tuning"). On x86-64 (a same-machine
    /// A/B of the targets sheet's hostbench, g++ 13.3, pinned core, 11 runs,
    /// at #42, before the Hexagon tuning; the predecessor's medians are review A's
    /// of tap/DspTap#42, this engine's measured the same way on the same
    /// machine, within 2 % of review A's own re-run of the previous tree),
    /// N = 256 … 4096, forward / inverse against the engine it replaced:
    /// -O3 float +0.5 … +5.7 % / +1.9 … +5.9 %, double +0.1 … +3.6 % /
    /// -3.4 … +1.0 %; -O3 -march=native float -9.0 … -13.4 % /
    /// -14.5 … -18.1 %, double -0.2 … +3.5 % / -2.2 … -7.5 %
    /// (informational; the design note has the table).
    ///
    /// GCC's basic-block (SLP) vectorizer is off for this class (the pragma
    /// below; a performance setting only). It packs the (re, im) halves of
    /// the complex values into vector registers and then pays a lane
    /// shuffle for every multiplication by i and every complex product,
    /// which costs more than it saves: with it on, the double forward at
    /// N = 512 took 1,403 ns at -O3 and 1,441 ns at -march=native (review A
    /// measured +14 % and +28 % against the predecessor), and the M55 float
    /// scenario (MVE) 1.9 - 2.3 % more instructions. The M4, M4F and M33 have no
    /// vector unit and compile the same code either way; clang's SLP
    /// vectorizer helps here and is left alone. The output bits do not
    /// depend on it (it reassociates nothing; the fingerprints are the same
    /// with and without).
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC push_options
#pragma GCC optimize("no-tree-slp-vectorize")
#endif
    template <std::floating_point Sample>
    class srdif_rdft {
        static_assert(std::is_same_v<Sample, float> || std::is_same_v<Sample, double>, "srdif_rdft: float or double");

      public:
        static constexpr std::size_t k_min_size     = 4;
        static constexpr std::size_t k_max_size     = std::size_t{1} << 30;
        static constexpr bool        k_is_shareable = true;

        explicit srdif_rdft(std::size_t n)
            : m_n(checked_size(n))
            , m_tables(make_tables(n)) {}

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

        /// Heap bytes one engine of size n holds: its single table allocation.
        [[nodiscard]] static constexpr std::size_t heap_bytes(std::size_t n) noexcept {
            return sizeof(Sample) * table_samples(n);
        }

      private:
        static constexpr Sample k_sqrt_half =
            std::is_same_v<Sample, float> ? Sample(0x1.6a09e6p-1) : Sample(0x1.6a09e667f3bcdp-1);

        static int log2_of(std::size_t n) noexcept { return std::bit_width(n) - 1; }

        /// The precondition, checked before any table is sized from n.
        static std::size_t checked_size(std::size_t n) noexcept {
            assert(n >= k_min_size && n <= k_max_size && (n & (n - 1)) == 0);
            return n;
        }

        // The one table allocation, in Samples: [post | small | kernel].
        //   post:   C_k for k < N/4, N/2 Samples (make_post_table).
        //   small:  the fused-pass twiddles at M = 32 (12 Samples) and M = 64
        //           (36), for block<32> and block<64>, where the kernel uses them.
        //   kernel: the fused-pass twiddles at M, twelve Samples per jj in
        //           [1, M/16), where the run-time pass runs (M >= 128).
        static constexpr std::size_t small_samples(std::size_t m) noexcept { return m >= 64 ? 48 : (m == 32 ? 12 : 0); }
        static constexpr std::size_t kernel_samples(std::size_t m) noexcept { return m >= 128 ? 12 * (m / 16 - 1) : 0; }
        static constexpr std::size_t table_samples(std::size_t n) noexcept {
            return n / 2 + small_samples(n / 2) + kernel_samples(n / 2);
        }

        static std::vector<Sample> make_tables(std::size_t n) {
            const std::size_t   m = n / 2;
            std::vector<Sample> t(table_samples(n));
            Sample* const       post = t.data();
            fill_post_table(post, n);
            Sample* const small = post + n / 2;
            if (m >= 32) {
                fill_kernel_table(small, 32, post, n);
            }
            if (m >= 64) {
                fill_kernel_table(small + 12, 64, post, n);
            }
            if (m >= 128) {
                fill_kernel_table(small + 48, m, post, n);
            }
            return t;
        }

        /// C_k = (1 + i W_n^k) / 2 for k in [0, n/4), interleaved (re, im):
        /// n/2 Samples at t.
        static void fill_post_table(Sample* t, std::size_t n) {
            const std::size_t quarter = n / 4;
            const std::size_t octant  = n / 8;
            if (quarter == 0) { // n < 4: outside the contract; write nothing
                return;
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
        }

        /// The fused pass's twiddles at resolution res (a power of two with
        /// 32 <= res <= n/2), twelve Samples per jj in [1, res/16): W^jj,
        /// W^3jj, W^(jj + res/8), W^(3jj + 3res/8), W^2jj, W^6jj (W =
        /// exp(+2 pi i/res), each as cos, sin), at t. Every value is read
        /// from the post table's imaginary parts, which hold cos theta_k / 2
        /// for every k < n/4 (theta_k = 2 pi k/n; srdif_trig's cos or sin of
        /// the first-octant angle, halved exactly): cos theta_r = 2 im_r and,
        /// for r > 0, sin theta_r = cos theta_(n/4 - r) = 2 im_(n/4 - r).
        /// Doubling undoes the halving exactly, so these are srdif_trig's
        /// values, the same bits a table built from srdif_trig at resolution
        /// res directly would hold (srdif_trig's value at (k, L) and at
        /// (2k, L + 1) is the same computation).
        static void fill_kernel_table(Sample* t, std::size_t res, const Sample* post, std::size_t n) {
            const std::size_t quarter = n / 4;
            const std::size_t step    = n / res; // index k at resolution n of angle 1 at res
            const auto        unit    = [&](std::size_t angle, Sample* out) {
                const std::size_t k  = (angle % res) * step;
                const std::size_t qd = k / quarter;
                const std::size_t r  = k % quarter;
                Sample            c  = Sample(2) * post[2 * r + 1];
                Sample            sn = r == 0 ? Sample(0) : Sample(2) * post[2 * (quarter - r) + 1];
                for (std::size_t q = 0; q < qd; ++q) { // times i per quarter turn
                    const Sample t0 = c;
                    c               = -sn;
                    sn              = t0;
                }
                out[0] = c;
                out[1] = sn;
            };
            for (std::size_t jj = 1; jj < res / 16; ++jj) {
                Sample* const e = t + 12 * (jj - 1);
                unit(jj, e);
                unit(3 * jj, e + 2);
                unit(jj + res / 8, e + 4);
                unit(3 * jj + 3 * (res / 8), e + 6);
                unit(2 * jj, e + 8);
                unit(6 * jj, e + 10);
            }
        }

        /// The bit-reversal permutation of the m = N/2 complex values, with no
        /// index table: a reversed counter walked beside the forward index
        /// (Gold and Rader's reverse-carry increment). Adding one to an index
        /// flips its trailing ones and the zero above them; the reversed
        /// counter flips the mirror image of those bits. With m = 2^b >= 16,
        /// write i = A m/4 + 4y + B (A, B < 4, y < m/16); then rev(i) =
        /// rev2(B) m/4 + 4 rev'(y) + rev2(A), rev2 reversing two bits and rev'
        /// the middle b - 4. One counter over y therefore serves sixteen
        /// indices: the six with A < rev2(B) are below their partners
        /// whatever y is and swap unconditionally, the four with
        /// A = rev2(B) swap when y < rev'(y), and the other six are the
        /// partners. Below m = 16 the counter walks every index.
        void permute(Sample* a) const noexcept {
            const std::size_t m = m_n / 2;
            if (m < 16) {
                std::size_t j = 0; // rev(i)
                for (std::size_t i = 0; i < m; ++i) {
                    if (i < j) {
                        swap2(a + 2 * i, a + 2 * j);
                    }
                    j ^= m - (m >> (std::countr_one(i) + 1));
                }
                return;
            }
            // Sample offsets: a quarter (A) is m/2 Samples, and i, j are the
            // offsets of complex indices 4y and 4 rev'(y) within a quarter.
            const std::size_t q = m / 2;
            Sample*           p = a; // row 0 at i; rows 1 ... 3 at p + q, p + 2q, p + 3q
            std::size_t       j = 0;
            for (std::size_t i = 0; i < q; i += 8, p += 8) {
                Sample* const s0 = a + j;
                Sample* const s1 = s0 + q;
                Sample* const s2 = s1 + q;
                Sample* const s3 = s2 + q;
                Sample* const p1 = p + q;
                Sample* const p2 = p1 + q;
                // (A, B) <-> (rev2 B, rev2 A), A < rev2 B.
                swap2(p + 2, s2);      // (0, 1) <-> (2, 0)
                swap2(p1 + 2, s2 + 4); // (1, 1) <-> (2, 2)
                swap2(p + 4, s1);      // (0, 2) <-> (1, 0)
                swap2(p + 6, s3);      // (0, 3) <-> (3, 0)
                swap2(p1 + 6, s3 + 4); // (1, 3) <-> (3, 2)
                swap2(p2 + 6, s3 + 2); // (2, 3) <-> (3, 1)
                if (i < j) {           // A = rev2 B: (0, 0), (2, 1), (1, 2), (3, 3)
                    swap2(p, s0);
                    swap2(p2 + 2, s2 + 2);
                    swap2(p1 + 4, s1 + 4);
                    swap2(p2 + q + 6, s3 + 6);
                }
                j ^= q - (q >> (std::countr_one(i >> 3) + 1));
            }
        }

        TAP_DSP_SRDIF_INLINE static void swap2(Sample* p, Sample* q) noexcept {
            const Sample p0 = p[0];
            const Sample p1 = p[1];
            const Sample q0 = q[0];
            const Sample q1 = q[1];
            p[0]            = q0;
            p[1]            = q1;
            q[0]            = p0;
            q[1]            = p1;
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

        TAP_DSP_SRDIF_INLINE static cv   ld(const Sample* p) noexcept { return {p[0], p[1]}; }
        TAP_DSP_SRDIF_INLINE static void st(Sample* p, cv v) noexcept {
            p[0] = v.r;
            p[1] = v.i;
        }

        /// a b - c d, the first product rounded in a statement of its own.
        ///
        /// Without contraction this is the same two products and one
        /// subtraction, rounded the same way, as the plain expression (the
        /// output bits and the fingerprints do not move). The float spelling
        /// is for the compilers that contract within a statement (clang,
        /// AppleClang): on `a * b - c * d` clang fuses the left product and
        /// negates the right one, fma(a, b, -(c d)), which costs Hexagon a
        /// separate negation in the two slots the float arithmetic itself
        /// needs; here the statement is `ab - c * d`, which it fuses as
        /// fma(-c, d, ab), one Hexagon multiply-subtract (Rx -= sfmpy(Rs,
        /// Rt)). Hexagon float scenarios -0.65 % / -0.43 % (N = 512 / 2048).
        /// GCC contracts across statements after SSA and compiles both
        /// spellings the same (the Cortex-M counts are identical). Double
        /// keeps the plain expression: Hexagon has no double fused
        /// multiply-add, and the split spelling there only moved the schedule
        /// (rfft_f64_512 +0.14 %).
        TAP_DSP_SRDIF_INLINE static Sample mul_sub(Sample a, Sample b, Sample c, Sample d) noexcept {
            if constexpr (std::is_same_v<Sample, float>) {
                const Sample ab = a * b;
                return ab - c * d;
            }
            else {
                return a * b - c * d;
            }
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
        TAP_DSP_SRDIF_INLINE static void bf(cv& x0, cv& x1, cv& x2, cv& x3, cv w1 = {}, cv w3 = {}) noexcept {
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
                    x2 = {ur * w1.r + ui * w1.i, mul_sub(ui, w1.r, ur, w1.i)};
                    x3 = {vr * w3.r + vi * w3.i, mul_sub(vi, w3.r, vr, w3.i)};
                }
                else {
                    x2 = {mul_sub(ur, w1.r, ui, w1.i), ur * w1.i + ui * w1.r};
                    x3 = {mul_sub(vr, w3.r, vi, w3.i), vr * w3.i + vi * w3.r};
                }
            }
        }

        /// One group of the fused pass: the eight values at p + r e (r = 0..7,
        /// complex units) of a length-l block with e = l/8. The level-l
        /// butterflies at j and j + e (twiddles a1, a3 and b1, b3), then the
        /// level-l/2 butterfly at j of the first half (c1, c3), whose four
        /// inputs are exactly the first-half outputs of the other two.
        template <bool Inverse, tw KA, tw KB, tw KC>
        TAP_DSP_SRDIF_INLINE static void group(Sample* p, std::size_t e, cv a1, cv a3, cv b1, cv b3, cv c1,
                                               cv c3) noexcept {
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
        TAP_DSP_SRDIF_INLINE static void group_first(Sample* a, std::size_t e) noexcept {
            // j = 0: W = 1, then the eighth-turn at j + e, then W = 1 at level l/2.
            group<Inverse, tw::one, tw::eighth, tw::one>(a, e, {}, {}, {}, {}, {}, {});
        }

        /// j in [1, e/2), twelve twiddles at w.
        template <bool Inverse>
        TAP_DSP_SRDIF_INLINE static void group_low(Sample* a, std::size_t e, const Sample* w) noexcept {
            group<Inverse, tw::general, tw::general, tw::general>(a, e, {w[0], w[1]}, {w[2], w[3]}, {w[4], w[5]},
                                                                  {w[6], w[7]}, {w[8], w[9]}, {w[10], w[11]});
        }

        /// j = e/2 = l/16: W_16^1, W_16^3; W_16^3, W_16^9 = -W_16^1; the eighth-turn at level l/2.
        template <bool Inverse>
        TAP_DSP_SRDIF_INLINE static void group_mid(Sample* a, std::size_t e) noexcept {
            const cv c16{k_cos_pi8, k_sin_pi8};
            const cv s16{k_sin_pi8, k_cos_pi8};
            group<Inverse, tw::general, tw::general, tw::eighth>(a, e, c16, s16, s16, {-k_cos_pi8, -k_sin_pi8}, {}, {});
        }

        /// j in (e/2, e), reading the entry of j' = e - j: W^j = i conj W^(j'+e),
        /// W^3j = -i conj W^3(j'+e), W^(j+e) = i conj W^j', W^3(j+e) = -i conj W^3j',
        /// and at level l/2 W^j = i conj W^j', W^3j = -i conj W^3j'.
        template <bool Inverse>
        TAP_DSP_SRDIF_INLINE static void group_high(Sample* a, std::size_t e, const Sample* w) noexcept {
            group<Inverse, tw::general, tw::general, tw::general>(a, e, {w[5], w[4]}, {-w[7], -w[6]}, {w[1], w[0]},
                                                                  {-w[3], -w[2]}, {w[9], w[8]}, {-w[11], -w[10]});
        }

        /// Whether the run-time fused pass runs one group per loop step
        /// (double) or two (float, whose eight row pointers then serve both
        /// groups at constant offsets). Two groups of doubles hold 32 values
        /// and 24 twiddles, twice Hexagon's register file, and spill;
        /// measured on Hexagon, one per step: rfft_f64_512 -0.38 %, the float
        /// scenarios +0.48 % / +0.99 % (N = 512 / 2048), so each keeps its
        /// best. The Cortex-M keys measure float only.
        static constexpr bool k_one_group_per_step = sizeof(Sample) == 8;

        /// The fused pass over a block of l >= 128 complex values: level l and
        /// level l/2 of its first half; the table entry of j is j (m/l) steps.
        template <bool Inverse>
        void fused_pass(Sample* a, std::size_t l) const noexcept {
            const std::size_t e      = l / 8;
            const std::size_t h      = e / 2;
            const std::size_t stride = 12 * ((m_n / 2) / l); // Samples per table step at this length
            group_first<Inverse>(a, e);
            const Sample* w = m_tables.data() + m_n / 2 + 48 + stride - 12;
            if constexpr (k_one_group_per_step) {
                for (std::size_t j = 1; j < h; ++j, w += stride) {
                    group_low<Inverse>(a + 2 * j, e, w);
                }
                group_mid<Inverse>(a + 2 * h, e);
                // j in (e/2, e) descending through the table: j' = e - j.
                for (std::size_t j = h + 1; j < e; ++j) {
                    w -= stride;
                    group_high<Inverse>(a + 2 * j, e, w);
                }
            }
            else {
                // j = 1 alone, then two groups per step (e/2 - 1 is odd): the
                // eight row pointers serve both groups at constant offsets.
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
        }

        /// A block of L <= 64 complex values with everything known at compile
        /// time: the fused pass over the small tables (L = 32, 64), one
        /// split-radix level in memory (L = 4, 8, 16) or the final pair.
        ///
        /// Below 32 every level runs in memory, one butterfly at a time,
        /// down to blocks of two; there are no register leaves. Until the
        /// Hexagon tuning blocks of 16 and fewer were register leaves (loaded
        /// once, transformed through every remaining level in registers,
        /// stored once), which on the ratchet's targets costs more than it
        /// saves: a leaf of 16 holds 32 values, the M4F's whole float register
        /// file and all of Hexagon's general registers (a double takes two),
        /// and spills. The arithmetic of every output is the same either way
        /// (the output bits do not move). Measured in this tree with register
        /// leaves of at most 16 / 8 / 4 values restored, against the pairs
        /// (docs/fft-design.md, "Hexagon (clang) tuning"), millions of
        /// instructions: Hexagon rfft_f32_512 46.10 / 45.88 / 45.89 / 45.94,
        /// rfft_f32_2048 51.79 / 51.58 / 51.55 / 51.61, rfft_f64_512 101.30 /
        /// 100.63 / 100.12 / 99.94; M4F rfft_f32_512 92.54 / 92.00 / 91.85 /
        /// 91.55, rfft_f32_2048 106.30 / 105.83 / 105.58 / 105.15; the M33 and
        /// the M55 without CMSIS 0.3 - 0.5 % below leaves of 8 as well, the
        /// soft-float M4 flat. Pairs cost Hexagon float 0.1 % against its best
        /// (leaves of 4 or 8) and are the best everywhere else, in one code
        /// path for both profiles.
        template <bool Inverse, std::size_t L>
        TAP_DSP_SRDIF_INLINE static void block(Sample* a, const Sample* small) noexcept {
            if constexpr (L >= 32) {
                constexpr std::size_t e = L / 8;
                constexpr std::size_t h = e / 2;
                const Sample* const   w = small + (L == 32 ? 0 : 12); // this length's table
                group_first<Inverse>(a, e);
                [&]<std::size_t... J>(std::index_sequence<J...>) TAP_DSP_SRDIF_INLINE {
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
            else if constexpr (L >= 4) {
                level<Inverse, L>(a);
                block<Inverse, L / 2>(a, small);
                if constexpr (L >= 8) { // blocks of one value are their own transforms
                    block<Inverse, L / 4>(a + L, small);
                    block<Inverse, L / 4>(a + L + L / 2, small);
                }
            }
            else {
                static_assert(L == 2, "blocks of 2 ... 64 values");
                (void)small; // the pairs need no table
                const cv x0 = ld(a);
                const cv x1 = ld(a + 2);
                st(a, {x0.r + x1.r, x0.i + x1.i});
                st(a + 2, {x0.r - x1.r, x0.i - x1.i});
            }
        }

        /// The first split-radix level of a block of L = 4, 8 or 16 values, in
        /// memory, one butterfly at a time: j = 0 (W = 1), j = L/8 (the
        /// eighth turn) and, at L = 16, j = 1 and 3 (W_16^1, W_16^3; W_16^3,
        /// W_16^9 = -W_16^1).
        template <bool Inverse, std::size_t L>
        TAP_DSP_SRDIF_INLINE static void level(Sample* a) noexcept {
            [&]<std::size_t... J>(std::index_sequence<J...>)
                TAP_DSP_SRDIF_INLINE { (level_bf<Inverse, L, J>(a), ...); }(std::make_index_sequence<L / 4>{});
        }

        template <bool Inverse, std::size_t L, std::size_t J>
        TAP_DSP_SRDIF_INLINE static void level_bf(Sample* a) noexcept {
            constexpr std::size_t q  = L / 4;
            cv                    x0 = ld(a + 2 * J);
            cv                    x1 = ld(a + 2 * (J + q));
            cv                    x2 = ld(a + 2 * (J + 2 * q));
            cv                    x3 = ld(a + 2 * (J + 3 * q));
            if constexpr (J == 0) {
                bf<Inverse, tw::one>(x0, x1, x2, x3);
            }
            else if constexpr (J == L / 8) {
                bf<Inverse, tw::eighth>(x0, x1, x2, x3);
            }
            else if constexpr (J == 1) { // L = 16: W_16^1, W_16^3
                bf<Inverse, tw::general>(x0, x1, x2, x3, {k_cos_pi8, k_sin_pi8}, {k_sin_pi8, k_cos_pi8});
            }
            else { // L = 16, J = 3: W_16^3, W_16^9 = -W_16^1
                bf<Inverse, tw::general>(x0, x1, x2, x3, {k_sin_pi8, k_cos_pi8}, {-k_cos_pi8, -k_sin_pi8});
            }
            st(a + 2 * J, x0);
            st(a + 2 * (J + q), x1);
            st(a + 2 * (J + 2 * q), x2);
            st(a + 2 * (J + 3 * q), x3);
        }

        /// Split-radix DIF over l complex values at a (bit-reversed output).
        template <bool Inverse>
        void kernel(Sample* a, std::size_t l) const noexcept {
            const Sample* const small = m_tables.data() + m_n / 2;
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
            const Sample* c   = m_tables.data() + 2;
            Sample* const end = a + m;
            post_pair<Inverse>(pk, pj, c);
            for (pk += 2, pj -= 2, c += 2; pk < end; pk += 4, pj -= 4, c += 4) {
                post_two<Inverse>(pk, pj, c);
            }
        }

        /// The four outputs of one bin pair of the post-pass.
        struct post_out {
            Sample kr;
            Sample ki;
            Sample jr;
            Sample ji;
        };

        /// One bin pair of the post-pass: u = a[k], v = conj a[m - k],
        /// G = C_k (u - v) (conj C_k for the inverse), a[k] <- u - G,
        /// a[m - k] <- conj(v + G); u, v and C_k given, the new a[k] and
        /// conj a[m - k] returned.
        template <bool Inverse>
        TAP_DSP_SRDIF_INLINE static post_out post_math(Sample ur, Sample ui, Sample vr, Sample vi, Sample cr,
                                                       Sample ci) noexcept {
            const Sample dr = ur - vr;
            const Sample di = ui + vi;
            Sample       gr;
            Sample       gi;
            if constexpr (Inverse) {
                gr = dr * cr + di * ci;
                gi = mul_sub(di, cr, dr, ci);
            }
            else {
                gr = mul_sub(dr, cr, di, ci);
                gi = dr * ci + di * cr;
            }
            return {ur - gr, ui - gi, vr + gr, vi - gi};
        }

        template <bool Inverse>
        TAP_DSP_SRDIF_INLINE static void post_pair(Sample* pk, Sample* pj, const Sample* c) noexcept {
            const post_out o = post_math<Inverse>(pk[0], pk[1], pj[0], pj[1], c[0], c[1]);
            pk[0]            = o.kr;
            pk[1]            = o.ki;
            pj[0]            = o.jr;
            pj[1]            = o.ji;
        }

        /// Two bin pairs, (k, m - k) and (k + 1, m - k - 1), every load before
        /// any store. The compiler cannot tell a store through pk from a later
        /// load through pj, so pair by pair the second pair's loads wait for
        /// the first pair's stores, and on Hexagon, which issues up to four
        /// instructions per packet and counts packets, the two pairs' arithmetic
        /// then cannot share them: without this, the Hexagon scenarios read
        /// +3.3 % / +3.0 % (float, N = 512 / 2048) and +1.0 % (double). The
        /// Cortex-M keys pay at most 0.02 % for it (the soft-float M4).
        template <bool Inverse>
        TAP_DSP_SRDIF_INLINE static void post_two(Sample* pk, Sample* pj, const Sample* c) noexcept {
            const Sample   u0r = pk[0];
            const Sample   u0i = pk[1];
            const Sample   u1r = pk[2];
            const Sample   u1i = pk[3];
            const Sample   v1r = pj[-2];
            const Sample   v1i = pj[-1];
            const Sample   v0r = pj[0];
            const Sample   v0i = pj[1];
            const post_out o0  = post_math<Inverse>(u0r, u0i, v0r, v0i, c[0], c[1]);
            const post_out o1  = post_math<Inverse>(u1r, u1i, v1r, v1i, c[2], c[3]);
            pk[0]              = o0.kr;
            pk[1]              = o0.ki;
            pk[2]              = o1.kr;
            pk[3]              = o1.ki;
            pj[-2]             = o1.jr;
            pj[-1]             = o1.ji;
            pj[0]              = o0.jr;
            pj[1]              = o0.ji;
        }

        std::size_t         m_n;
        std::vector<Sample> m_tables;
    };
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC pop_options
#endif

} // namespace tap::dsp::detail
