/// @file split_radix.h
/// @brief Split-radix real FFT engine: a C++20 transliteration of Ooura's rdft.
// SPDX-License-Identifier: LicenseRef-Ooura AND MIT
// Copyright(C) 1996-2001 Takuya OOURA
// Copyright 2026 Timothy Place and the DspTap contributors (modifications).
// Derived from fftsg.c in Takuya Ooura's General Purpose FFT Package; this is
// a derivative work, not the ORIGINAL package. Modifications (the C++20
// template form, the constructor-time tables, the explicit conversions that
// spell out the C's implicit ones) September 2026.
//
// The original package's notice, which governs the derived portion (the
// transform itself), verbatim from its readme.txt:
//
//     Copyright(C) 1996-2001 Takuya OOURA
//         email: ooura@mmm.t.u-tokyo.ac.jp
//         download: http://momonga.t.u-tokyo.ac.jp/~ooura/fft.html
//         You may use, copy, modify this code for any purpose and
//         without fee. You may distribute this ORIGINAL package.
//
// The wrapper and DspTap's additions are MIT (see LICENSE). The full
// statement is NOTICE.md; the original package's readme.txt is kept in-tree
// at third_party/ooura/readme.txt, and LICENSES/LicenseRef-Ooura.txt carries
// the notice text the SPDX reference denotes.
//
// ---------------------------------------------------------------------------
// TRANSLITERATION RULES — read before editing anything below the class
// docstring. The engine is BIT-IDENTICAL to the vendored C (fftsg.c for
// double, fftsg_float.c for float) and tests/test_fft_parity_ooura.cpp is
// the gate that holds it there: memcmp identity, forward and inverse, at
// every power of two from 4 to 65536 plus 2^20, both precisions, with both
// sides compiled at -ffp-contract=off. These rules are what the gate depends
// on (docs/audit-fft-and-code-smells.md, Part 4; docs/fft-design.md).
//
// 1. STATEMENT FIDELITY. Every Ooura arithmetic statement stays textually
//    intact: same operands, same order, same grouping, one statement per
//    Ooura statement. No cmul helpers, no lambdas, no std::complex, no
//    reassociation, no hoisting of a repeated product into a local. Clang
//    contracts a*b+c into an FMA within a statement only; GCC contracts
//    across statements after inlining. Any restructuring changes which
//    products fuse at default flags and, worse, changes the operation
//    sequence itself, which breaks the identity on one compiler or the
//    other without a location to point at.
//
// 2. TABLE SEMANTICS PER PRECISION (Decision D10). The C's float build is
//    fftsg.c compiled under `#define double float`, which retargets every
//    local in makewt/makect to float — delta, wn4r, wk1r, ... — and the
//    products fed to the trig calls (`delta * j`, `3 * delta * j`) to float
//    arithmetic, while cos/sin/atan remain the DOUBLE libm functions, called
//    with the float argument promoted and their double result converted to
//    float on assignment; `0.5 / cos(...)`, `0.5 / wk1r` and `0.5 * cos(...)`
//    are double arithmetic converted on assignment. In C++ `std::cos(float)`
//    would call cosf and the tables would differ from the C. So: the locals
//    are typed Sample, every libm call takes an explicit
//    `static_cast<double>(...)` argument, and the result is converted to
//    Sample with an explicit cast exactly where the C's assignment converted
//    it. Each explicit cast below reproduces an implicit conversion of the C
//    (int operand to Sample, double result to Sample); none introduces one.
//
// 3. TABLES ARE BUILT ONCE, IN THE CONSTRUCTOR, with the first-call protocol
//    of the C's rdft (ip[0] = 0 requests the build) executed there and
//    nowhere else. The lazy `nw`/`nc` re-initialization is gone, which is
//    what makes the transforms const and the object shareable across threads
//    once constructed.
//
// 4. INDEX TYPES STAY THE C's `int`; the tables are addressed through raw
//    pointers loaded into locals once per transform, so the helpers receive
//    `a`, `ip` and `w` by parameter exactly as the C functions do and the
//    aliasing analysis is the same. The one std::size_t -> int narrowing is
//    the constructor's.
//
// Scope: rdft and what it reaches — makewt, makeipt, makect, bitrv2,
// bitrv2conj, bitrv216, bitrv216neg, bitrv208, bitrv208neg, cftfsub, cftbsub,
// cftf1st, cftb1st, cftrec4, cfttree, cftleaf, cftmdl1, cftmdl2, cftfx41,
// cftf161, cftf162, cftf081, cftf082, cftf040, cftb040, cftx020, rftfsub,
// rftbsub. Not ported: cdft, the DCT/DST family (ddct, ddst, dfct, dfst,
// dctsub, dstsub) and the USE_CDFT_PTHREADS / USE_CDFT_WINTHREADS
// scaffolding. Nothing in the reachable set self-recurses (cftrec4 is a
// while plus a for over cfttree/cftleaf).
// ---------------------------------------------------------------------------

#pragma once

#include <cassert>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <vector>

namespace tap::dsp::detail {

    /// Ooura's split-radix real discrete Fourier transform (rdft), as a
    /// C++20 class template over the sample type: `double` is the golden
    /// model and `float` the embedded profile, and each instantiation is
    /// bit-identical to the corresponding build of the vendored C (see the
    /// transliteration rules in the file banner).
    ///
    /// Contract (the numbers `basic_real_fft` re-presents; this class is the
    /// engine, not the consumer-facing surface):
    ///   - Size N is a power of two, N >= 4, fixed at construction. The
    ///     bit-reversal table (2 + sqrt(N/2) + 1 ints) and the trig table
    ///     (N/2 samples) are allocated and built in the constructor; the
    ///     transforms allocate nothing and are `noexcept` and `const`.
    ///   - Packing after forward_inplace of N real samples: a[0] = DC,
    ///     a[1] = Nyquist (both real), a[2k] / a[2k+1] = re / im of bin k
    ///     for 1 <= k < N/2.
    ///   - Sign: A[k] = sum_j a[j] W^(jk) with W = exp(+2*pi*i/N), i.e. the
    ///     imaginary parts are conjugated relative to the engineering DFT.
    ///   - Forward scale 1; inverse_inplace is UNNORMALIZED: the caller
    ///     applies 2/N for a round trip.
    ///   - No alignment requirement on the data pointer. NaN propagates to
    ///     every bin (no data-dependent branches). Latency 0.
    ///   - Copyable; a copy is bit-identical to its source. Two threads may
    ///     run transforms on one object concurrently (nothing is mutated
    ///     after construction).
    ///
    /// The transforms' arithmetic is Ooura's, statement for statement; the
    /// class adds the geometry check, the table build and the const surface.
    template <std::floating_point Sample>
    class split_radix_rdft {
      public:
        /// Builds the tables for a transform of @p n real samples, exactly as
        /// the C's rdft builds them on its first call (makewt for the
        /// complex stages, makect for the real post-pass; makeipt inside
        /// makewt for the bit-reversal permutation).
        /// @param n transform size: a power of two, n >= 4 (asserted).
        explicit split_radix_rdft(std::size_t n)
            : m_size(static_cast<int>(n))
            , m_ip(2 + static_cast<std::size_t>(std::sqrt(static_cast<double>(n) / 2.0)) + 1, 0)
            , m_w(n / 2, Sample(0)) {
            assert(n >= 4 && (n & (n - 1)) == 0);
            // rdft's first-call protocol (ip[0] = 0, ip[1] = 0 on entry), run
            // once here. makewt sets ip[1] = 1, so for n = 4 makect is never
            // called and nc stays 1: the real post-pass is not run at n = 4.
            int*    ip = m_ip.data();
            Sample* w  = m_w.data();
            int     nw = ip[0];
            if (m_size > (nw << 2)) {
                nw = m_size >> 2;
                makewt(nw, ip, w);
            }
            int nc = ip[1];
            if (m_size > (nc << 2)) {
                nc = m_size >> 2;
                makect(nc, ip, w + nw);
            }
            m_nw = nw;
            m_nc = nc;
        }

        /// @return the transform size N.
        [[nodiscard]] std::size_t size() const noexcept { return static_cast<std::size_t>(m_size); }

      private:
        // -------- initializing routines (rule 2 applies to every line) --------

        static void makewt(int nw, int* ip, Sample* w) noexcept {
            int    j, nwh, nw0, nw1;
            Sample delta, wn4r, wk1r, wk1i, wk3r, wk3i;

            ip[0] = nw;
            ip[1] = 1;
            if (nw > 2) {
                nwh   = nw >> 1;
                delta = static_cast<Sample>(std::atan(1.0) / nwh);
                wn4r  = static_cast<Sample>(std::cos(static_cast<double>(delta * static_cast<Sample>(nwh))));
                w[0]  = 1;
                w[1]  = wn4r;
                if (nwh == 4) {
                    w[2] = static_cast<Sample>(std::cos(static_cast<double>(delta * 2)));
                    w[3] = static_cast<Sample>(std::sin(static_cast<double>(delta * 2)));
                }
                else if (nwh > 4) {
                    makeipt(nw, ip);
                    w[2] = static_cast<Sample>(0.5 / std::cos(static_cast<double>(delta * 2)));
                    w[3] = static_cast<Sample>(0.5 / std::cos(static_cast<double>(delta * 6)));
                    for (j = 4; j < nwh; j += 4) {
                        w[j]     = static_cast<Sample>(std::cos(static_cast<double>(delta * static_cast<Sample>(j))));
                        w[j + 1] = static_cast<Sample>(std::sin(static_cast<double>(delta * static_cast<Sample>(j))));
                        w[j + 2] =
                            static_cast<Sample>(std::cos(static_cast<double>(3 * delta * static_cast<Sample>(j))));
                        w[j + 3] =
                            static_cast<Sample>(-std::sin(static_cast<double>(3 * delta * static_cast<Sample>(j))));
                    }
                }
                nw0 = 0;
                while (nwh > 2) {
                    nw1 = nw0 + nwh;
                    nwh >>= 1;
                    w[nw1]     = 1;
                    w[nw1 + 1] = wn4r;
                    if (nwh == 4) {
                        wk1r       = w[nw0 + 4];
                        wk1i       = w[nw0 + 5];
                        w[nw1 + 2] = wk1r;
                        w[nw1 + 3] = wk1i;
                    }
                    else if (nwh > 4) {
                        wk1r       = w[nw0 + 4];
                        wk3r       = w[nw0 + 6];
                        w[nw1 + 2] = static_cast<Sample>(0.5 / wk1r);
                        w[nw1 + 3] = static_cast<Sample>(0.5 / wk3r);
                        for (j = 4; j < nwh; j += 4) {
                            wk1r           = w[nw0 + 2 * j];
                            wk1i           = w[nw0 + 2 * j + 1];
                            wk3r           = w[nw0 + 2 * j + 2];
                            wk3i           = w[nw0 + 2 * j + 3];
                            w[nw1 + j]     = wk1r;
                            w[nw1 + j + 1] = wk1i;
                            w[nw1 + j + 2] = wk3r;
                            w[nw1 + j + 3] = wk3i;
                        }
                    }
                    nw0 = nw1;
                }
            }
        }

        static void makeipt(int nw, int* ip) noexcept {
            int j, l, m, m2, p, q;

            ip[2] = 0;
            ip[3] = 16;
            m     = 2;
            for (l = nw; l > 32; l >>= 2) {
                m2 = m << 1;
                q  = m2 << 3;
                for (j = m; j < m2; j++) {
                    p          = ip[j] << 2;
                    ip[m + j]  = p;
                    ip[m2 + j] = p + q;
                }
                m = m2;
            }
        }

        static void makect(int nc, int* ip, Sample* c) noexcept {
            int    j, nch;
            Sample delta;

            ip[1] = nc;
            if (nc > 1) {
                nch    = nc >> 1;
                delta  = static_cast<Sample>(std::atan(1.0) / nch);
                c[0]   = static_cast<Sample>(std::cos(static_cast<double>(delta * static_cast<Sample>(nch))));
                c[nch] = static_cast<Sample>(0.5 * c[0]);
                for (j = 1; j < nch; j++) {
                    c[j] = static_cast<Sample>(0.5 * std::cos(static_cast<double>(delta * static_cast<Sample>(j))));
                    c[nc - j] =
                        static_cast<Sample>(0.5 * std::sin(static_cast<double>(delta * static_cast<Sample>(j))));
                }
            }
        }

        int              m_size; ///< N, the transform size
        int              m_nw;   ///< length of the complex-stage twiddle table (N/4), as rdft computes it
        int              m_nc;   ///< length of the real post-pass table (N/4; 1 at N = 4, where it is unused)
        std::vector<int> m_ip;   ///< bit-reversal work table: 2 + sqrt(N/2) + 1 entries
        std::vector<Sample>
            m_w; ///< trig tables: m_w[0, m_nw) for the complex stages, m_w[m_nw, N/2) for rftfsub/rftbsub
    };

} // namespace tap::dsp::detail
