/// @file tables.h
/// @brief Table generators for the fixed-point real FFT: bit reversal and Q1.30 twiddles.
// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// The three tables detail::fixed_point_rdft (fft/fixed_point.h) builds at
// construction, as free functions returning std::vector so the test battery
// can pin them independently of any transform: the bit-reversal permutation
// of the complex kernel, the kernel's Q1.30 twiddles, and the Q1.30
// coefficients of Ooura's real post-pass. Shared with nothing else by design
// (audit Part 7): the floating profiles keep Ooura's table semantics for bit
// identity with the vendored C, and these tables are what the fixed-point
// contract is pinned against instead.
//
// Every coefficient is generated in double through std::cos / std::sin and
// rounded ONCE by fft_arith<std::int32_t>::make_coeff (round half away from
// zero, saturating), so |w_q - w| <= 0.5 LSB of Q1.30 on every host
// (`TwiddleTableIsWithinHalfLsb` in the battery). Host libm last-bit
// differences can move a double that lies within 2^-31 of a rounding
// boundary onto the other side, so fixed-point transform outputs are
// host-identical only if the table is; the battery pins each certified N's
// table checksum (FNV-1a-64 over the int32 bit patterns, in index order) so
// a libm difference is detected rather than silently absorbed
// (`TwiddleTableChecksumIsPinned`).

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
    /// The angle is formed as 2*pi*k/m in double from the integer k, one
    /// rounding per angle; for m <= 2^16 the angle error (< 2^-48 rad) is
    /// three orders of magnitude below the Q1.30 quantum, so the 0.5 LSB
    /// bound holds with margin.
    ///
    /// @pre m is a power of two >= 2 (the kernel's length, N/2 for a real
    ///      transform of N).
    inline std::vector<fft_arith<std::int32_t>::coeff> make_twiddle_table(std::size_t m) {
        assert(m >= 2 && (m & (m - 1)) == 0);
        using arith = fft_arith<std::int32_t>;
        std::vector<arith::coeff> table(2 * m);
        const double              step = 2.0 * std::numbers::pi / static_cast<double>(m);
        for (std::size_t k = 0; k < m; ++k) {
            const double angle = step * static_cast<double>(k);
            table[2 * k]       = arith::make_coeff(std::cos(angle));
            table[2 * k + 1]   = arith::make_coeff(std::sin(angle));
        }
        return table;
    }

    /// Q1.30 coefficients of the real post-pass (Ooura's rftfsub / rftbsub
    /// formulas, third_party/ooura/fftsg.c), for a real transform of length n.
    ///
    /// Ooura's makect stores 0.5*cos(2*pi*j/n) and reads the pair
    ///   wkr = 0.5 - 0.5*sin(2*pi*k/n),   wki = 0.5*cos(2*pi*k/n)
    /// for bin k; this table stores that pair directly so each value is one
    /// rounding from its double. Layout: 2*(n/4) entries, interleaved; for
    /// k in [0, n/4):
    ///   table[2k]     = make_coeff(0.5 - 0.5*sin(2*pi*k/n))
    ///   table[2k + 1] = make_coeff(0.5*cos(2*pi*k/n))
    /// The post-pass runs k over [1, n/4) (bins 1 .. n/4 - 1 paired with their
    /// mirrors n/2 - k); entry 0 is (0.5, 0.5) and unused, kept so that the
    /// index is the bin number. |wkr + i*wki| = sqrt(0.5*(1 - sin)) <= 1/sqrt(2),
    /// so the products never approach the coefficient's own range.
    ///
    /// @pre n is a power of two >= 4 (n == 4 yields the one unused entry).
    inline std::vector<fft_arith<std::int32_t>::coeff> make_real_post_pass_table(std::size_t n) {
        assert(n >= 4 && (n & (n - 1)) == 0);
        using arith                       = fft_arith<std::int32_t>;
        const auto                quarter = n / 4;
        std::vector<arith::coeff> table(2 * quarter);
        const double              step = 2.0 * std::numbers::pi / static_cast<double>(n);
        for (std::size_t k = 0; k < quarter; ++k) {
            const double angle = step * static_cast<double>(k);
            table[2 * k]       = arith::make_coeff(0.5 - 0.5 * std::sin(angle));
            table[2 * k + 1]   = arith::make_coeff(0.5 * std::cos(angle));
        }
        return table;
    }

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
