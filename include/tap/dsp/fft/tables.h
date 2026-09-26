/// @file tables.h
/// @brief Table generators for the fixed-point real FFT: bit reversal and Q1.30 twiddles.
// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// The three tables detail::fixed_point_rdft (fft/fixed_point.h) builds at
// construction, as free functions returning std::vector so the test battery
// can pin them independently of any transform: the bit-reversal permutation
// of the complex kernel, the kernel's Q1.30 twiddles, and the Q1.30
// coefficients of the real post-pass. Shared with nothing else by design
// (audit Part 7): the floating engine keeps its own tables, and these tables
// are what the fixed-point contract is pinned against instead.
//
// Every coefficient is generated in double through std::cos / std::sin and
// rounded ONCE by fft_arith<std::int32_t>::make_coeff (round half away from
// zero, saturating), so |w_q - w| <= 0.5 LSB of Q1.30 on every host
// (`TwiddleTableIsWithinHalfLsb` in the battery); the kernel twiddles go
// through libm for the first octant only and reach the rest of the circle by
// exact symmetry. Two things can still move a coefficient by one Q1.30 LSB
// between hosts: a libm last-bit difference (glibc, newlib, UCRT, Apple),
// and fp-contraction of a generator expression (Decision D9: a compiler that
// fuses an a - b*c into one rounding, the default on Apple arm64 and on the
// M55 leg), either of which can carry a double that lies within 2^-31 of a
// rounding boundary onto the other side. The generators therefore contain no
// contractible expression (<<CLEAN-ROOM: the post-pass generator's example>>),
// and fixed-point transform outputs
// are host-identical only if the table is: the battery pins each certified
// N's table checksum (FNV-1a-64 over the int32 bit patterns, in index order)
// so a libm difference is detected rather than silently absorbed
// (`TwiddleTableChecksumIsPinned`), and pins the transforms' own output
// fingerprints on top (`OutputFingerprintIsPinned`).

#pragma once

#include <bit>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numbers>
#include <vector>

#include "tap/dsp/fft/fft_arith.h"

namespace tap::dsp::detail {

    /// Bit-reversal permutation for a complex transform of length m.
    ///
    /// Layout: m entries; table[i] is i with its log2(m) low bits reversed,
    /// so a decimation-in-frequency kernel that leaves element X[table[i]] at
    /// position i is put in natural order by swapping every pair (i, table[i])
    /// with i < table[i]. table[0] == 0 and table[m - 1] == m - 1 always.
    ///
    /// @pre m is a power of two (m == 1 yields the one-entry identity).
    inline std::vector<std::uint32_t> make_bit_reversal_table(std::size_t m) {
        assert(m >= 1 && (m & (m - 1)) == 0);
        const int                  bits = std::bit_width(m) - 1;
        std::vector<std::uint32_t> table(m);
        for (std::size_t i = 0; i < m; ++i) {
            std::uint32_t r = 0;
            for (int b = 0; b < bits; ++b) {
                r |= static_cast<std::uint32_t>((i >> b) & 1u) << (bits - 1 - b);
            }
            table[i] = r;
        }
        return table;
    }

    /// Q1.30 twiddles of the complex kernel: the m roots of unity in the
    /// library's exp(+i) convention.
    ///
    /// Layout: 2m entries, interleaved; for k in [0, m):
    ///   table[2k]     = make_coeff(cos(2*pi*k/m))
    ///   table[2k + 1] = make_coeff(sin(2*pi*k/m))
    /// i.e. W_m^k = table[2k] + i*table[2k+1] with W_m = exp(+2*pi*i/m). The
    /// forward kernel multiplies by W, the inverse by its conjugate (same
    /// table, the sine term subtracted instead of added). k = 0 is exactly
    /// 1.0 (k_coeff_one = 2^30) and k = m/4 is exactly i, since make_coeff
    /// rounds the double 1.0 / 0.0 exactly; the kernel relies on neither.
    ///
    /// Only the first octant, k in [0, m/8], goes through libm: cos and sin
    /// of the angle 2*pi*k/m, formed in double from the integer k (one
    /// rounding per angle; for the kernel's m = N/2 <= 2^15 the angle error,
    /// < 2^-49 rad, is five orders of magnitude below the Q1.30 quantum, so
    /// the 0.5 LSB bound holds with margin), each rounded once by make_coeff.
    /// At k = m/8 the sine is set equal to the cosine (the exact value is
    /// sqrt(1/2) for both). The other seven octants are the first one's
    /// quantized values under the exact symmetries of the circle,
    ///   W^(m/4 - k) = i conj(W^k),  W^(m/2 - k) = -conj(W^k),  W^(m - k) = conj(W^k),
    /// applied to the Q1.30 integers (a negation is exact in Q1.30, and
    /// make_coeff rounds half away from zero, so a negated coefficient is the
    /// coefficient of the negated double). Consequences: the k <-> m - k
    /// conjugate symmetry and the k <-> m/4 - k swap hold bit-exactly on every
    /// host by construction, and a libm difference can enter only through the
    /// m/8 + 1 first-octant evaluations rather than through 2m of them
    /// (`TwiddleTableIsWithinHalfLsb` checks the symmetries and the bound).
    ///
    /// @pre m is a power of two >= 2 (the kernel's length, N/2 for a real
    ///      transform of N).
    inline std::vector<fft_arith<std::int32_t>::coeff> make_twiddle_table(std::size_t m) {
        assert(m >= 2 && (m & (m - 1)) == 0);
        using arith = fft_arith<std::int32_t>;
        std::vector<arith::coeff> table(2 * m);
        const std::size_t         octant  = m / 8; // 0 for m < 8: only k = 0 is evaluated
        const std::size_t         quarter = m / 4;
        const std::size_t         half    = m / 2;
        const double              step    = 2.0 * std::numbers::pi / static_cast<double>(m);
        for (std::size_t k = 0; k <= octant; ++k) {
            const double angle = step * static_cast<double>(k);
            table[2 * k]       = arith::make_coeff(std::cos(angle));
            table[2 * k + 1]   = (k == octant && octant > 0) ? table[2 * k] : arith::make_coeff(std::sin(angle));
        }
        for (std::size_t k = octant + 1; k <= quarter; ++k) { // W^k = i conj(W^(m/4 - k)): (c, s) -> (s, c)
            const std::size_t j = quarter - k;
            table[2 * k]        = table[2 * j + 1];
            table[2 * k + 1]    = table[2 * j];
        }
        for (std::size_t k = quarter + 1; k <= half; ++k) { // W^k = -conj(W^(m/2 - k)): (c, s) -> (-c, s)
            const std::size_t j = half - k;
            table[2 * k]        = static_cast<arith::coeff>(-table[2 * j]);
            table[2 * k + 1]    = table[2 * j + 1];
        }
        for (std::size_t k = half + 1; k < m; ++k) { // W^k = conj(W^(m - k)): (c, s) -> (c, -s)
            const std::size_t j = m - k;
            table[2 * k]        = table[2 * j];
            table[2 * k + 1]    = static_cast<arith::coeff>(-table[2 * j + 1]);
        }
        return table;
    }

    /// <<CLEAN-ROOM: make_real_post_pass_table (or whatever the derived
    /// post-pass needs), its layout, its accuracy statement, and its
    /// no-contractible-expression note (file header, D9).>>

    /// FNV-1a-64 over a table's bit patterns in index order: the checksum the
    /// battery pins per certified N so a host libm difference in the table is
    /// detected rather than absorbed. Each int32 is folded as its unsigned
    /// 32-bit pattern, least significant byte first.
    inline std::uint64_t table_checksum(const std::int32_t* table, std::size_t count) noexcept {
        std::uint64_t h = 0xcbf29ce484222325ull;
        for (std::size_t i = 0; i < count; ++i) {
            const auto bits = static_cast<std::uint32_t>(table[i]);
            for (int b = 0; b < 4; ++b) {
                h ^= (bits >> (8 * b)) & 0xffu;
                h *= 0x100000001b3ull;
            }
        }
        return h;
    }

} // namespace tap::dsp::detail
