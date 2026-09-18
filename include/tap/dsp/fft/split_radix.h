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

        // -------- child routines: Ooura's statements, verbatim (rule 1) --------

        static void bitrv2(int n, const int* ip, Sample* a) noexcept {
            int    j, j1, k, k1, l, m, nh, nm;
            Sample xr, xi, yr, yi;

            m = 1;
            for (l = n >> 2; l > 8; l >>= 2) {
                m <<= 1;
            }
            nh = n >> 1;
            nm = 4 * m;
            if (l == 8) {
                for (k = 0; k < m; k++) {
                    for (j = 0; j < k; j++) {
                        j1        = 4 * j + 2 * ip[m + k];
                        k1        = 4 * k + 2 * ip[m + j];
                        xr        = a[j1];
                        xi        = a[j1 + 1];
                        yr        = a[k1];
                        yi        = a[k1 + 1];
                        a[j1]     = yr;
                        a[j1 + 1] = yi;
                        a[k1]     = xr;
                        a[k1 + 1] = xi;
                        j1 += nm;
                        k1 += 2 * nm;
                        xr        = a[j1];
                        xi        = a[j1 + 1];
                        yr        = a[k1];
                        yi        = a[k1 + 1];
                        a[j1]     = yr;
                        a[j1 + 1] = yi;
                        a[k1]     = xr;
                        a[k1 + 1] = xi;
                        j1 += nm;
                        k1 -= nm;
                        xr        = a[j1];
                        xi        = a[j1 + 1];
                        yr        = a[k1];
                        yi        = a[k1 + 1];
                        a[j1]     = yr;
                        a[j1 + 1] = yi;
                        a[k1]     = xr;
                        a[k1 + 1] = xi;
                        j1 += nm;
                        k1 += 2 * nm;
                        xr        = a[j1];
                        xi        = a[j1 + 1];
                        yr        = a[k1];
                        yi        = a[k1 + 1];
                        a[j1]     = yr;
                        a[j1 + 1] = yi;
                        a[k1]     = xr;
                        a[k1 + 1] = xi;
                        j1 += nh;
                        k1 += 2;
                        xr        = a[j1];
                        xi        = a[j1 + 1];
                        yr        = a[k1];
                        yi        = a[k1 + 1];
                        a[j1]     = yr;
                        a[j1 + 1] = yi;
                        a[k1]     = xr;
                        a[k1 + 1] = xi;
                        j1 -= nm;
                        k1 -= 2 * nm;
                        xr        = a[j1];
                        xi        = a[j1 + 1];
                        yr        = a[k1];
                        yi        = a[k1 + 1];
                        a[j1]     = yr;
                        a[j1 + 1] = yi;
                        a[k1]     = xr;
                        a[k1 + 1] = xi;
                        j1 -= nm;
                        k1 += nm;
                        xr        = a[j1];
                        xi        = a[j1 + 1];
                        yr        = a[k1];
                        yi        = a[k1 + 1];
                        a[j1]     = yr;
                        a[j1 + 1] = yi;
                        a[k1]     = xr;
                        a[k1 + 1] = xi;
                        j1 -= nm;
                        k1 -= 2 * nm;
                        xr        = a[j1];
                        xi        = a[j1 + 1];
                        yr        = a[k1];
                        yi        = a[k1 + 1];
                        a[j1]     = yr;
                        a[j1 + 1] = yi;
                        a[k1]     = xr;
                        a[k1 + 1] = xi;
                        j1 += 2;
                        k1 += nh;
                        xr        = a[j1];
                        xi        = a[j1 + 1];
                        yr        = a[k1];
                        yi        = a[k1 + 1];
                        a[j1]     = yr;
                        a[j1 + 1] = yi;
                        a[k1]     = xr;
                        a[k1 + 1] = xi;
                        j1 += nm;
                        k1 += 2 * nm;
                        xr        = a[j1];
                        xi        = a[j1 + 1];
                        yr        = a[k1];
                        yi        = a[k1 + 1];
                        a[j1]     = yr;
                        a[j1 + 1] = yi;
                        a[k1]     = xr;
                        a[k1 + 1] = xi;
                        j1 += nm;
                        k1 -= nm;
                        xr        = a[j1];
                        xi        = a[j1 + 1];
                        yr        = a[k1];
                        yi        = a[k1 + 1];
                        a[j1]     = yr;
                        a[j1 + 1] = yi;
                        a[k1]     = xr;
                        a[k1 + 1] = xi;
                        j1 += nm;
                        k1 += 2 * nm;
                        xr        = a[j1];
                        xi        = a[j1 + 1];
                        yr        = a[k1];
                        yi        = a[k1 + 1];
                        a[j1]     = yr;
                        a[j1 + 1] = yi;
                        a[k1]     = xr;
                        a[k1 + 1] = xi;
                        j1 -= nh;
                        k1 -= 2;
                        xr        = a[j1];
                        xi        = a[j1 + 1];
                        yr        = a[k1];
                        yi        = a[k1 + 1];
                        a[j1]     = yr;
                        a[j1 + 1] = yi;
                        a[k1]     = xr;
                        a[k1 + 1] = xi;
                        j1 -= nm;
                        k1 -= 2 * nm;
                        xr        = a[j1];
                        xi        = a[j1 + 1];
                        yr        = a[k1];
                        yi        = a[k1 + 1];
                        a[j1]     = yr;
                        a[j1 + 1] = yi;
                        a[k1]     = xr;
                        a[k1 + 1] = xi;
                        j1 -= nm;
                        k1 += nm;
                        xr        = a[j1];
                        xi        = a[j1 + 1];
                        yr        = a[k1];
                        yi        = a[k1 + 1];
                        a[j1]     = yr;
                        a[j1 + 1] = yi;
                        a[k1]     = xr;
                        a[k1 + 1] = xi;
                        j1 -= nm;
                        k1 -= 2 * nm;
                        xr        = a[j1];
                        xi        = a[j1 + 1];
                        yr        = a[k1];
                        yi        = a[k1 + 1];
                        a[j1]     = yr;
                        a[j1 + 1] = yi;
                        a[k1]     = xr;
                        a[k1 + 1] = xi;
                    }
                    k1 = 4 * k + 2 * ip[m + k];
                    j1 = k1 + 2;
                    k1 += nh;
                    xr        = a[j1];
                    xi        = a[j1 + 1];
                    yr        = a[k1];
                    yi        = a[k1 + 1];
                    a[j1]     = yr;
                    a[j1 + 1] = yi;
                    a[k1]     = xr;
                    a[k1 + 1] = xi;
                    j1 += nm;
                    k1 += 2 * nm;
                    xr        = a[j1];
                    xi        = a[j1 + 1];
                    yr        = a[k1];
                    yi        = a[k1 + 1];
                    a[j1]     = yr;
                    a[j1 + 1] = yi;
                    a[k1]     = xr;
                    a[k1 + 1] = xi;
                    j1 += nm;
                    k1 -= nm;
                    xr        = a[j1];
                    xi        = a[j1 + 1];
                    yr        = a[k1];
                    yi        = a[k1 + 1];
                    a[j1]     = yr;
                    a[j1 + 1] = yi;
                    a[k1]     = xr;
                    a[k1 + 1] = xi;
                    j1 -= 2;
                    k1 -= nh;
                    xr        = a[j1];
                    xi        = a[j1 + 1];
                    yr        = a[k1];
                    yi        = a[k1 + 1];
                    a[j1]     = yr;
                    a[j1 + 1] = yi;
                    a[k1]     = xr;
                    a[k1 + 1] = xi;
                    j1 += nh + 2;
                    k1 += nh + 2;
                    xr        = a[j1];
                    xi        = a[j1 + 1];
                    yr        = a[k1];
                    yi        = a[k1 + 1];
                    a[j1]     = yr;
                    a[j1 + 1] = yi;
                    a[k1]     = xr;
                    a[k1 + 1] = xi;
                    j1 -= nh - nm;
                    k1 += 2 * nm - 2;
                    xr        = a[j1];
                    xi        = a[j1 + 1];
                    yr        = a[k1];
                    yi        = a[k1 + 1];
                    a[j1]     = yr;
                    a[j1 + 1] = yi;
                    a[k1]     = xr;
                    a[k1 + 1] = xi;
                }
            }
            else {
                for (k = 0; k < m; k++) {
                    for (j = 0; j < k; j++) {
                        j1        = 4 * j + ip[m + k];
                        k1        = 4 * k + ip[m + j];
                        xr        = a[j1];
                        xi        = a[j1 + 1];
                        yr        = a[k1];
                        yi        = a[k1 + 1];
                        a[j1]     = yr;
                        a[j1 + 1] = yi;
                        a[k1]     = xr;
                        a[k1 + 1] = xi;
                        j1 += nm;
                        k1 += nm;
                        xr        = a[j1];
                        xi        = a[j1 + 1];
                        yr        = a[k1];
                        yi        = a[k1 + 1];
                        a[j1]     = yr;
                        a[j1 + 1] = yi;
                        a[k1]     = xr;
                        a[k1 + 1] = xi;
                        j1 += nh;
                        k1 += 2;
                        xr        = a[j1];
                        xi        = a[j1 + 1];
                        yr        = a[k1];
                        yi        = a[k1 + 1];
                        a[j1]     = yr;
                        a[j1 + 1] = yi;
                        a[k1]     = xr;
                        a[k1 + 1] = xi;
                        j1 -= nm;
                        k1 -= nm;
                        xr        = a[j1];
                        xi        = a[j1 + 1];
                        yr        = a[k1];
                        yi        = a[k1 + 1];
                        a[j1]     = yr;
                        a[j1 + 1] = yi;
                        a[k1]     = xr;
                        a[k1 + 1] = xi;
                        j1 += 2;
                        k1 += nh;
                        xr        = a[j1];
                        xi        = a[j1 + 1];
                        yr        = a[k1];
                        yi        = a[k1 + 1];
                        a[j1]     = yr;
                        a[j1 + 1] = yi;
                        a[k1]     = xr;
                        a[k1 + 1] = xi;
                        j1 += nm;
                        k1 += nm;
                        xr        = a[j1];
                        xi        = a[j1 + 1];
                        yr        = a[k1];
                        yi        = a[k1 + 1];
                        a[j1]     = yr;
                        a[j1 + 1] = yi;
                        a[k1]     = xr;
                        a[k1 + 1] = xi;
                        j1 -= nh;
                        k1 -= 2;
                        xr        = a[j1];
                        xi        = a[j1 + 1];
                        yr        = a[k1];
                        yi        = a[k1 + 1];
                        a[j1]     = yr;
                        a[j1 + 1] = yi;
                        a[k1]     = xr;
                        a[k1 + 1] = xi;
                        j1 -= nm;
                        k1 -= nm;
                        xr        = a[j1];
                        xi        = a[j1 + 1];
                        yr        = a[k1];
                        yi        = a[k1 + 1];
                        a[j1]     = yr;
                        a[j1 + 1] = yi;
                        a[k1]     = xr;
                        a[k1 + 1] = xi;
                    }
                    k1 = 4 * k + ip[m + k];
                    j1 = k1 + 2;
                    k1 += nh;
                    xr        = a[j1];
                    xi        = a[j1 + 1];
                    yr        = a[k1];
                    yi        = a[k1 + 1];
                    a[j1]     = yr;
                    a[j1 + 1] = yi;
                    a[k1]     = xr;
                    a[k1 + 1] = xi;
                    j1 += nm;
                    k1 += nm;
                    xr        = a[j1];
                    xi        = a[j1 + 1];
                    yr        = a[k1];
                    yi        = a[k1 + 1];
                    a[j1]     = yr;
                    a[j1 + 1] = yi;
                    a[k1]     = xr;
                    a[k1 + 1] = xi;
                }
            }
        }

        static void bitrv2conj(int n, const int* ip, Sample* a) noexcept {
            int    j, j1, k, k1, l, m, nh, nm;
            Sample xr, xi, yr, yi;

            m = 1;
            for (l = n >> 2; l > 8; l >>= 2) {
                m <<= 1;
            }
            nh = n >> 1;
            nm = 4 * m;
            if (l == 8) {
                for (k = 0; k < m; k++) {
                    for (j = 0; j < k; j++) {
                        j1        = 4 * j + 2 * ip[m + k];
                        k1        = 4 * k + 2 * ip[m + j];
                        xr        = a[j1];
                        xi        = -a[j1 + 1];
                        yr        = a[k1];
                        yi        = -a[k1 + 1];
                        a[j1]     = yr;
                        a[j1 + 1] = yi;
                        a[k1]     = xr;
                        a[k1 + 1] = xi;
                        j1 += nm;
                        k1 += 2 * nm;
                        xr        = a[j1];
                        xi        = -a[j1 + 1];
                        yr        = a[k1];
                        yi        = -a[k1 + 1];
                        a[j1]     = yr;
                        a[j1 + 1] = yi;
                        a[k1]     = xr;
                        a[k1 + 1] = xi;
                        j1 += nm;
                        k1 -= nm;
                        xr        = a[j1];
                        xi        = -a[j1 + 1];
                        yr        = a[k1];
                        yi        = -a[k1 + 1];
                        a[j1]     = yr;
                        a[j1 + 1] = yi;
                        a[k1]     = xr;
                        a[k1 + 1] = xi;
                        j1 += nm;
                        k1 += 2 * nm;
                        xr        = a[j1];
                        xi        = -a[j1 + 1];
                        yr        = a[k1];
                        yi        = -a[k1 + 1];
                        a[j1]     = yr;
                        a[j1 + 1] = yi;
                        a[k1]     = xr;
                        a[k1 + 1] = xi;
                        j1 += nh;
                        k1 += 2;
                        xr        = a[j1];
                        xi        = -a[j1 + 1];
                        yr        = a[k1];
                        yi        = -a[k1 + 1];
                        a[j1]     = yr;
                        a[j1 + 1] = yi;
                        a[k1]     = xr;
                        a[k1 + 1] = xi;
                        j1 -= nm;
                        k1 -= 2 * nm;
                        xr        = a[j1];
                        xi        = -a[j1 + 1];
                        yr        = a[k1];
                        yi        = -a[k1 + 1];
                        a[j1]     = yr;
                        a[j1 + 1] = yi;
                        a[k1]     = xr;
                        a[k1 + 1] = xi;
                        j1 -= nm;
                        k1 += nm;
                        xr        = a[j1];
                        xi        = -a[j1 + 1];
                        yr        = a[k1];
                        yi        = -a[k1 + 1];
                        a[j1]     = yr;
                        a[j1 + 1] = yi;
                        a[k1]     = xr;
                        a[k1 + 1] = xi;
                        j1 -= nm;
                        k1 -= 2 * nm;
                        xr        = a[j1];
                        xi        = -a[j1 + 1];
                        yr        = a[k1];
                        yi        = -a[k1 + 1];
                        a[j1]     = yr;
                        a[j1 + 1] = yi;
                        a[k1]     = xr;
                        a[k1 + 1] = xi;
                        j1 += 2;
                        k1 += nh;
                        xr        = a[j1];
                        xi        = -a[j1 + 1];
                        yr        = a[k1];
                        yi        = -a[k1 + 1];
                        a[j1]     = yr;
                        a[j1 + 1] = yi;
                        a[k1]     = xr;
                        a[k1 + 1] = xi;
                        j1 += nm;
                        k1 += 2 * nm;
                        xr        = a[j1];
                        xi        = -a[j1 + 1];
                        yr        = a[k1];
                        yi        = -a[k1 + 1];
                        a[j1]     = yr;
                        a[j1 + 1] = yi;
                        a[k1]     = xr;
                        a[k1 + 1] = xi;
                        j1 += nm;
                        k1 -= nm;
                        xr        = a[j1];
                        xi        = -a[j1 + 1];
                        yr        = a[k1];
                        yi        = -a[k1 + 1];
                        a[j1]     = yr;
                        a[j1 + 1] = yi;
                        a[k1]     = xr;
                        a[k1 + 1] = xi;
                        j1 += nm;
                        k1 += 2 * nm;
                        xr        = a[j1];
                        xi        = -a[j1 + 1];
                        yr        = a[k1];
                        yi        = -a[k1 + 1];
                        a[j1]     = yr;
                        a[j1 + 1] = yi;
                        a[k1]     = xr;
                        a[k1 + 1] = xi;
                        j1 -= nh;
                        k1 -= 2;
                        xr        = a[j1];
                        xi        = -a[j1 + 1];
                        yr        = a[k1];
                        yi        = -a[k1 + 1];
                        a[j1]     = yr;
                        a[j1 + 1] = yi;
                        a[k1]     = xr;
                        a[k1 + 1] = xi;
                        j1 -= nm;
                        k1 -= 2 * nm;
                        xr        = a[j1];
                        xi        = -a[j1 + 1];
                        yr        = a[k1];
                        yi        = -a[k1 + 1];
                        a[j1]     = yr;
                        a[j1 + 1] = yi;
                        a[k1]     = xr;
                        a[k1 + 1] = xi;
                        j1 -= nm;
                        k1 += nm;
                        xr        = a[j1];
                        xi        = -a[j1 + 1];
                        yr        = a[k1];
                        yi        = -a[k1 + 1];
                        a[j1]     = yr;
                        a[j1 + 1] = yi;
                        a[k1]     = xr;
                        a[k1 + 1] = xi;
                        j1 -= nm;
                        k1 -= 2 * nm;
                        xr        = a[j1];
                        xi        = -a[j1 + 1];
                        yr        = a[k1];
                        yi        = -a[k1 + 1];
                        a[j1]     = yr;
                        a[j1 + 1] = yi;
                        a[k1]     = xr;
                        a[k1 + 1] = xi;
                    }
                    k1 = 4 * k + 2 * ip[m + k];
                    j1 = k1 + 2;
                    k1 += nh;
                    a[j1 - 1] = -a[j1 - 1];
                    xr        = a[j1];
                    xi        = -a[j1 + 1];
                    yr        = a[k1];
                    yi        = -a[k1 + 1];
                    a[j1]     = yr;
                    a[j1 + 1] = yi;
                    a[k1]     = xr;
                    a[k1 + 1] = xi;
                    a[k1 + 3] = -a[k1 + 3];
                    j1 += nm;
                    k1 += 2 * nm;
                    xr        = a[j1];
                    xi        = -a[j1 + 1];
                    yr        = a[k1];
                    yi        = -a[k1 + 1];
                    a[j1]     = yr;
                    a[j1 + 1] = yi;
                    a[k1]     = xr;
                    a[k1 + 1] = xi;
                    j1 += nm;
                    k1 -= nm;
                    xr        = a[j1];
                    xi        = -a[j1 + 1];
                    yr        = a[k1];
                    yi        = -a[k1 + 1];
                    a[j1]     = yr;
                    a[j1 + 1] = yi;
                    a[k1]     = xr;
                    a[k1 + 1] = xi;
                    j1 -= 2;
                    k1 -= nh;
                    xr        = a[j1];
                    xi        = -a[j1 + 1];
                    yr        = a[k1];
                    yi        = -a[k1 + 1];
                    a[j1]     = yr;
                    a[j1 + 1] = yi;
                    a[k1]     = xr;
                    a[k1 + 1] = xi;
                    j1 += nh + 2;
                    k1 += nh + 2;
                    xr        = a[j1];
                    xi        = -a[j1 + 1];
                    yr        = a[k1];
                    yi        = -a[k1 + 1];
                    a[j1]     = yr;
                    a[j1 + 1] = yi;
                    a[k1]     = xr;
                    a[k1 + 1] = xi;
                    j1 -= nh - nm;
                    k1 += 2 * nm - 2;
                    a[j1 - 1] = -a[j1 - 1];
                    xr        = a[j1];
                    xi        = -a[j1 + 1];
                    yr        = a[k1];
                    yi        = -a[k1 + 1];
                    a[j1]     = yr;
                    a[j1 + 1] = yi;
                    a[k1]     = xr;
                    a[k1 + 1] = xi;
                    a[k1 + 3] = -a[k1 + 3];
                }
            }
            else {
                for (k = 0; k < m; k++) {
                    for (j = 0; j < k; j++) {
                        j1        = 4 * j + ip[m + k];
                        k1        = 4 * k + ip[m + j];
                        xr        = a[j1];
                        xi        = -a[j1 + 1];
                        yr        = a[k1];
                        yi        = -a[k1 + 1];
                        a[j1]     = yr;
                        a[j1 + 1] = yi;
                        a[k1]     = xr;
                        a[k1 + 1] = xi;
                        j1 += nm;
                        k1 += nm;
                        xr        = a[j1];
                        xi        = -a[j1 + 1];
                        yr        = a[k1];
                        yi        = -a[k1 + 1];
                        a[j1]     = yr;
                        a[j1 + 1] = yi;
                        a[k1]     = xr;
                        a[k1 + 1] = xi;
                        j1 += nh;
                        k1 += 2;
                        xr        = a[j1];
                        xi        = -a[j1 + 1];
                        yr        = a[k1];
                        yi        = -a[k1 + 1];
                        a[j1]     = yr;
                        a[j1 + 1] = yi;
                        a[k1]     = xr;
                        a[k1 + 1] = xi;
                        j1 -= nm;
                        k1 -= nm;
                        xr        = a[j1];
                        xi        = -a[j1 + 1];
                        yr        = a[k1];
                        yi        = -a[k1 + 1];
                        a[j1]     = yr;
                        a[j1 + 1] = yi;
                        a[k1]     = xr;
                        a[k1 + 1] = xi;
                        j1 += 2;
                        k1 += nh;
                        xr        = a[j1];
                        xi        = -a[j1 + 1];
                        yr        = a[k1];
                        yi        = -a[k1 + 1];
                        a[j1]     = yr;
                        a[j1 + 1] = yi;
                        a[k1]     = xr;
                        a[k1 + 1] = xi;
                        j1 += nm;
                        k1 += nm;
                        xr        = a[j1];
                        xi        = -a[j1 + 1];
                        yr        = a[k1];
                        yi        = -a[k1 + 1];
                        a[j1]     = yr;
                        a[j1 + 1] = yi;
                        a[k1]     = xr;
                        a[k1 + 1] = xi;
                        j1 -= nh;
                        k1 -= 2;
                        xr        = a[j1];
                        xi        = -a[j1 + 1];
                        yr        = a[k1];
                        yi        = -a[k1 + 1];
                        a[j1]     = yr;
                        a[j1 + 1] = yi;
                        a[k1]     = xr;
                        a[k1 + 1] = xi;
                        j1 -= nm;
                        k1 -= nm;
                        xr        = a[j1];
                        xi        = -a[j1 + 1];
                        yr        = a[k1];
                        yi        = -a[k1 + 1];
                        a[j1]     = yr;
                        a[j1 + 1] = yi;
                        a[k1]     = xr;
                        a[k1 + 1] = xi;
                    }
                    k1 = 4 * k + ip[m + k];
                    j1 = k1 + 2;
                    k1 += nh;
                    a[j1 - 1] = -a[j1 - 1];
                    xr        = a[j1];
                    xi        = -a[j1 + 1];
                    yr        = a[k1];
                    yi        = -a[k1 + 1];
                    a[j1]     = yr;
                    a[j1 + 1] = yi;
                    a[k1]     = xr;
                    a[k1 + 1] = xi;
                    a[k1 + 3] = -a[k1 + 3];
                    j1 += nm;
                    k1 += nm;
                    a[j1 - 1] = -a[j1 - 1];
                    xr        = a[j1];
                    xi        = -a[j1 + 1];
                    yr        = a[k1];
                    yi        = -a[k1 + 1];
                    a[j1]     = yr;
                    a[j1 + 1] = yi;
                    a[k1]     = xr;
                    a[k1 + 1] = xi;
                    a[k1 + 3] = -a[k1 + 3];
                }
            }
        }

        static void bitrv216(Sample* a) noexcept {
            Sample x1r, x1i, x2r, x2i, x3r, x3i, x4r, x4i, x5r, x5i, x7r, x7i, x8r, x8i, x10r, x10i, x11r, x11i, x12r,
                x12i, x13r, x13i, x14r, x14i;

            x1r   = a[2];
            x1i   = a[3];
            x2r   = a[4];
            x2i   = a[5];
            x3r   = a[6];
            x3i   = a[7];
            x4r   = a[8];
            x4i   = a[9];
            x5r   = a[10];
            x5i   = a[11];
            x7r   = a[14];
            x7i   = a[15];
            x8r   = a[16];
            x8i   = a[17];
            x10r  = a[20];
            x10i  = a[21];
            x11r  = a[22];
            x11i  = a[23];
            x12r  = a[24];
            x12i  = a[25];
            x13r  = a[26];
            x13i  = a[27];
            x14r  = a[28];
            x14i  = a[29];
            a[2]  = x8r;
            a[3]  = x8i;
            a[4]  = x4r;
            a[5]  = x4i;
            a[6]  = x12r;
            a[7]  = x12i;
            a[8]  = x2r;
            a[9]  = x2i;
            a[10] = x10r;
            a[11] = x10i;
            a[14] = x14r;
            a[15] = x14i;
            a[16] = x1r;
            a[17] = x1i;
            a[20] = x5r;
            a[21] = x5i;
            a[22] = x13r;
            a[23] = x13i;
            a[24] = x3r;
            a[25] = x3i;
            a[26] = x11r;
            a[27] = x11i;
            a[28] = x7r;
            a[29] = x7i;
        }

        static void bitrv216neg(Sample* a) noexcept {
            Sample x1r, x1i, x2r, x2i, x3r, x3i, x4r, x4i, x5r, x5i, x6r, x6i, x7r, x7i, x8r, x8i, x9r, x9i, x10r, x10i,
                x11r, x11i, x12r, x12i, x13r, x13i, x14r, x14i, x15r, x15i;

            x1r   = a[2];
            x1i   = a[3];
            x2r   = a[4];
            x2i   = a[5];
            x3r   = a[6];
            x3i   = a[7];
            x4r   = a[8];
            x4i   = a[9];
            x5r   = a[10];
            x5i   = a[11];
            x6r   = a[12];
            x6i   = a[13];
            x7r   = a[14];
            x7i   = a[15];
            x8r   = a[16];
            x8i   = a[17];
            x9r   = a[18];
            x9i   = a[19];
            x10r  = a[20];
            x10i  = a[21];
            x11r  = a[22];
            x11i  = a[23];
            x12r  = a[24];
            x12i  = a[25];
            x13r  = a[26];
            x13i  = a[27];
            x14r  = a[28];
            x14i  = a[29];
            x15r  = a[30];
            x15i  = a[31];
            a[2]  = x15r;
            a[3]  = x15i;
            a[4]  = x7r;
            a[5]  = x7i;
            a[6]  = x11r;
            a[7]  = x11i;
            a[8]  = x3r;
            a[9]  = x3i;
            a[10] = x13r;
            a[11] = x13i;
            a[12] = x5r;
            a[13] = x5i;
            a[14] = x9r;
            a[15] = x9i;
            a[16] = x1r;
            a[17] = x1i;
            a[18] = x14r;
            a[19] = x14i;
            a[20] = x6r;
            a[21] = x6i;
            a[22] = x10r;
            a[23] = x10i;
            a[24] = x2r;
            a[25] = x2i;
            a[26] = x12r;
            a[27] = x12i;
            a[28] = x4r;
            a[29] = x4i;
            a[30] = x8r;
            a[31] = x8i;
        }

        static void bitrv208(Sample* a) noexcept {
            Sample x1r, x1i, x3r, x3i, x4r, x4i, x6r, x6i;

            x1r   = a[2];
            x1i   = a[3];
            x3r   = a[6];
            x3i   = a[7];
            x4r   = a[8];
            x4i   = a[9];
            x6r   = a[12];
            x6i   = a[13];
            a[2]  = x4r;
            a[3]  = x4i;
            a[6]  = x6r;
            a[7]  = x6i;
            a[8]  = x1r;
            a[9]  = x1i;
            a[12] = x3r;
            a[13] = x3i;
        }

        static void bitrv208neg(Sample* a) noexcept {
            Sample x1r, x1i, x2r, x2i, x3r, x3i, x4r, x4i, x5r, x5i, x6r, x6i, x7r, x7i;

            x1r   = a[2];
            x1i   = a[3];
            x2r   = a[4];
            x2i   = a[5];
            x3r   = a[6];
            x3i   = a[7];
            x4r   = a[8];
            x4i   = a[9];
            x5r   = a[10];
            x5i   = a[11];
            x6r   = a[12];
            x6i   = a[13];
            x7r   = a[14];
            x7i   = a[15];
            a[2]  = x7r;
            a[3]  = x7i;
            a[4]  = x3r;
            a[5]  = x3i;
            a[6]  = x5r;
            a[7]  = x5i;
            a[8]  = x1r;
            a[9]  = x1i;
            a[10] = x6r;
            a[11] = x6i;
            a[12] = x2r;
            a[13] = x2i;
            a[14] = x4r;
            a[15] = x4i;
        }

        static void cftfsub(int n, Sample* a, const int* ip, int nw, const Sample* w) noexcept {
            if (n > 8) {
                if (n > 32) {
                    cftf1st(n, a, &w[nw - (n >> 2)]);
                    if (n > 512) {
                        cftrec4(n, a, nw, w);
                    }
                    else if (n > 128) {
                        cftleaf(n, 1, a, nw, w);
                    }
                    else {
                        cftfx41(n, a, nw, w);
                    }
                    bitrv2(n, ip, a);
                }
                else if (n == 32) {
                    cftf161(a, &w[nw - 8]);
                    bitrv216(a);
                }
                else {
                    cftf081(a, w);
                    bitrv208(a);
                }
            }
            else if (n == 8) {
                cftf040(a);
            }
            else if (n == 4) {
                cftx020(a);
            }
        }

        static void cftbsub(int n, Sample* a, const int* ip, int nw, const Sample* w) noexcept {
            if (n > 8) {
                if (n > 32) {
                    cftb1st(n, a, &w[nw - (n >> 2)]);
                    if (n > 512) {
                        cftrec4(n, a, nw, w);
                    }
                    else if (n > 128) {
                        cftleaf(n, 1, a, nw, w);
                    }
                    else {
                        cftfx41(n, a, nw, w);
                    }
                    bitrv2conj(n, ip, a);
                }
                else if (n == 32) {
                    cftf161(a, &w[nw - 8]);
                    bitrv216neg(a);
                }
                else {
                    cftf081(a, w);
                    bitrv208neg(a);
                }
            }
            else if (n == 8) {
                cftb040(a);
            }
            else if (n == 4) {
                cftx020(a);
            }
        }

        static void cftf1st(int n, Sample* a, const Sample* w) noexcept {
            int    j, j0, j1, j2, j3, k, m, mh;
            Sample wn4r, csc1, csc3, wk1r, wk1i, wk3r, wk3i, wd1r, wd1i, wd3r, wd3i;
            Sample x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i, y0r, y0i, y1r, y1i, y2r, y2i, y3r, y3i;

            mh        = n >> 3;
            m         = 2 * mh;
            j1        = m;
            j2        = j1 + m;
            j3        = j2 + m;
            x0r       = a[0] + a[j2];
            x0i       = a[1] + a[j2 + 1];
            x1r       = a[0] - a[j2];
            x1i       = a[1] - a[j2 + 1];
            x2r       = a[j1] + a[j3];
            x2i       = a[j1 + 1] + a[j3 + 1];
            x3r       = a[j1] - a[j3];
            x3i       = a[j1 + 1] - a[j3 + 1];
            a[0]      = x0r + x2r;
            a[1]      = x0i + x2i;
            a[j1]     = x0r - x2r;
            a[j1 + 1] = x0i - x2i;
            a[j2]     = x1r - x3i;
            a[j2 + 1] = x1i + x3r;
            a[j3]     = x1r + x3i;
            a[j3 + 1] = x1i - x3r;
            wn4r      = w[1];
            csc1      = w[2];
            csc3      = w[3];
            wd1r      = 1;
            wd1i      = 0;
            wd3r      = 1;
            wd3i      = 0;
            k         = 0;
            for (j = 2; j < mh - 2; j += 4) {
                k += 4;
                wk1r      = csc1 * (wd1r + w[k]);
                wk1i      = csc1 * (wd1i + w[k + 1]);
                wk3r      = csc3 * (wd3r + w[k + 2]);
                wk3i      = csc3 * (wd3i + w[k + 3]);
                wd1r      = w[k];
                wd1i      = w[k + 1];
                wd3r      = w[k + 2];
                wd3i      = w[k + 3];
                j1        = j + m;
                j2        = j1 + m;
                j3        = j2 + m;
                x0r       = a[j] + a[j2];
                x0i       = a[j + 1] + a[j2 + 1];
                x1r       = a[j] - a[j2];
                x1i       = a[j + 1] - a[j2 + 1];
                y0r       = a[j + 2] + a[j2 + 2];
                y0i       = a[j + 3] + a[j2 + 3];
                y1r       = a[j + 2] - a[j2 + 2];
                y1i       = a[j + 3] - a[j2 + 3];
                x2r       = a[j1] + a[j3];
                x2i       = a[j1 + 1] + a[j3 + 1];
                x3r       = a[j1] - a[j3];
                x3i       = a[j1 + 1] - a[j3 + 1];
                y2r       = a[j1 + 2] + a[j3 + 2];
                y2i       = a[j1 + 3] + a[j3 + 3];
                y3r       = a[j1 + 2] - a[j3 + 2];
                y3i       = a[j1 + 3] - a[j3 + 3];
                a[j]      = x0r + x2r;
                a[j + 1]  = x0i + x2i;
                a[j + 2]  = y0r + y2r;
                a[j + 3]  = y0i + y2i;
                a[j1]     = x0r - x2r;
                a[j1 + 1] = x0i - x2i;
                a[j1 + 2] = y0r - y2r;
                a[j1 + 3] = y0i - y2i;
                x0r       = x1r - x3i;
                x0i       = x1i + x3r;
                a[j2]     = wk1r * x0r - wk1i * x0i;
                a[j2 + 1] = wk1r * x0i + wk1i * x0r;
                x0r       = y1r - y3i;
                x0i       = y1i + y3r;
                a[j2 + 2] = wd1r * x0r - wd1i * x0i;
                a[j2 + 3] = wd1r * x0i + wd1i * x0r;
                x0r       = x1r + x3i;
                x0i       = x1i - x3r;
                a[j3]     = wk3r * x0r + wk3i * x0i;
                a[j3 + 1] = wk3r * x0i - wk3i * x0r;
                x0r       = y1r + y3i;
                x0i       = y1i - y3r;
                a[j3 + 2] = wd3r * x0r + wd3i * x0i;
                a[j3 + 3] = wd3r * x0i - wd3i * x0r;
                j0        = m - j;
                j1        = j0 + m;
                j2        = j1 + m;
                j3        = j2 + m;
                x0r       = a[j0] + a[j2];
                x0i       = a[j0 + 1] + a[j2 + 1];
                x1r       = a[j0] - a[j2];
                x1i       = a[j0 + 1] - a[j2 + 1];
                y0r       = a[j0 - 2] + a[j2 - 2];
                y0i       = a[j0 - 1] + a[j2 - 1];
                y1r       = a[j0 - 2] - a[j2 - 2];
                y1i       = a[j0 - 1] - a[j2 - 1];
                x2r       = a[j1] + a[j3];
                x2i       = a[j1 + 1] + a[j3 + 1];
                x3r       = a[j1] - a[j3];
                x3i       = a[j1 + 1] - a[j3 + 1];
                y2r       = a[j1 - 2] + a[j3 - 2];
                y2i       = a[j1 - 1] + a[j3 - 1];
                y3r       = a[j1 - 2] - a[j3 - 2];
                y3i       = a[j1 - 1] - a[j3 - 1];
                a[j0]     = x0r + x2r;
                a[j0 + 1] = x0i + x2i;
                a[j0 - 2] = y0r + y2r;
                a[j0 - 1] = y0i + y2i;
                a[j1]     = x0r - x2r;
                a[j1 + 1] = x0i - x2i;
                a[j1 - 2] = y0r - y2r;
                a[j1 - 1] = y0i - y2i;
                x0r       = x1r - x3i;
                x0i       = x1i + x3r;
                a[j2]     = wk1i * x0r - wk1r * x0i;
                a[j2 + 1] = wk1i * x0i + wk1r * x0r;
                x0r       = y1r - y3i;
                x0i       = y1i + y3r;
                a[j2 - 2] = wd1i * x0r - wd1r * x0i;
                a[j2 - 1] = wd1i * x0i + wd1r * x0r;
                x0r       = x1r + x3i;
                x0i       = x1i - x3r;
                a[j3]     = wk3i * x0r + wk3r * x0i;
                a[j3 + 1] = wk3i * x0i - wk3r * x0r;
                x0r       = y1r + y3i;
                x0i       = y1i - y3r;
                a[j3 - 2] = wd3i * x0r + wd3r * x0i;
                a[j3 - 1] = wd3i * x0i - wd3r * x0r;
            }
            wk1r      = csc1 * (wd1r + wn4r);
            wk1i      = csc1 * (wd1i + wn4r);
            wk3r      = csc3 * (wd3r - wn4r);
            wk3i      = csc3 * (wd3i - wn4r);
            j0        = mh;
            j1        = j0 + m;
            j2        = j1 + m;
            j3        = j2 + m;
            x0r       = a[j0 - 2] + a[j2 - 2];
            x0i       = a[j0 - 1] + a[j2 - 1];
            x1r       = a[j0 - 2] - a[j2 - 2];
            x1i       = a[j0 - 1] - a[j2 - 1];
            x2r       = a[j1 - 2] + a[j3 - 2];
            x2i       = a[j1 - 1] + a[j3 - 1];
            x3r       = a[j1 - 2] - a[j3 - 2];
            x3i       = a[j1 - 1] - a[j3 - 1];
            a[j0 - 2] = x0r + x2r;
            a[j0 - 1] = x0i + x2i;
            a[j1 - 2] = x0r - x2r;
            a[j1 - 1] = x0i - x2i;
            x0r       = x1r - x3i;
            x0i       = x1i + x3r;
            a[j2 - 2] = wk1r * x0r - wk1i * x0i;
            a[j2 - 1] = wk1r * x0i + wk1i * x0r;
            x0r       = x1r + x3i;
            x0i       = x1i - x3r;
            a[j3 - 2] = wk3r * x0r + wk3i * x0i;
            a[j3 - 1] = wk3r * x0i - wk3i * x0r;
            x0r       = a[j0] + a[j2];
            x0i       = a[j0 + 1] + a[j2 + 1];
            x1r       = a[j0] - a[j2];
            x1i       = a[j0 + 1] - a[j2 + 1];
            x2r       = a[j1] + a[j3];
            x2i       = a[j1 + 1] + a[j3 + 1];
            x3r       = a[j1] - a[j3];
            x3i       = a[j1 + 1] - a[j3 + 1];
            a[j0]     = x0r + x2r;
            a[j0 + 1] = x0i + x2i;
            a[j1]     = x0r - x2r;
            a[j1 + 1] = x0i - x2i;
            x0r       = x1r - x3i;
            x0i       = x1i + x3r;
            a[j2]     = wn4r * (x0r - x0i);
            a[j2 + 1] = wn4r * (x0i + x0r);
            x0r       = x1r + x3i;
            x0i       = x1i - x3r;
            a[j3]     = -wn4r * (x0r + x0i);
            a[j3 + 1] = -wn4r * (x0i - x0r);
            x0r       = a[j0 + 2] + a[j2 + 2];
            x0i       = a[j0 + 3] + a[j2 + 3];
            x1r       = a[j0 + 2] - a[j2 + 2];
            x1i       = a[j0 + 3] - a[j2 + 3];
            x2r       = a[j1 + 2] + a[j3 + 2];
            x2i       = a[j1 + 3] + a[j3 + 3];
            x3r       = a[j1 + 2] - a[j3 + 2];
            x3i       = a[j1 + 3] - a[j3 + 3];
            a[j0 + 2] = x0r + x2r;
            a[j0 + 3] = x0i + x2i;
            a[j1 + 2] = x0r - x2r;
            a[j1 + 3] = x0i - x2i;
            x0r       = x1r - x3i;
            x0i       = x1i + x3r;
            a[j2 + 2] = wk1i * x0r - wk1r * x0i;
            a[j2 + 3] = wk1i * x0i + wk1r * x0r;
            x0r       = x1r + x3i;
            x0i       = x1i - x3r;
            a[j3 + 2] = wk3i * x0r + wk3r * x0i;
            a[j3 + 3] = wk3i * x0i - wk3r * x0r;
        }

        static void cftb1st(int n, Sample* a, const Sample* w) noexcept {
            int    j, j0, j1, j2, j3, k, m, mh;
            Sample wn4r, csc1, csc3, wk1r, wk1i, wk3r, wk3i, wd1r, wd1i, wd3r, wd3i;
            Sample x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i, y0r, y0i, y1r, y1i, y2r, y2i, y3r, y3i;

            mh        = n >> 3;
            m         = 2 * mh;
            j1        = m;
            j2        = j1 + m;
            j3        = j2 + m;
            x0r       = a[0] + a[j2];
            x0i       = -a[1] - a[j2 + 1];
            x1r       = a[0] - a[j2];
            x1i       = -a[1] + a[j2 + 1];
            x2r       = a[j1] + a[j3];
            x2i       = a[j1 + 1] + a[j3 + 1];
            x3r       = a[j1] - a[j3];
            x3i       = a[j1 + 1] - a[j3 + 1];
            a[0]      = x0r + x2r;
            a[1]      = x0i - x2i;
            a[j1]     = x0r - x2r;
            a[j1 + 1] = x0i + x2i;
            a[j2]     = x1r + x3i;
            a[j2 + 1] = x1i + x3r;
            a[j3]     = x1r - x3i;
            a[j3 + 1] = x1i - x3r;
            wn4r      = w[1];
            csc1      = w[2];
            csc3      = w[3];
            wd1r      = 1;
            wd1i      = 0;
            wd3r      = 1;
            wd3i      = 0;
            k         = 0;
            for (j = 2; j < mh - 2; j += 4) {
                k += 4;
                wk1r      = csc1 * (wd1r + w[k]);
                wk1i      = csc1 * (wd1i + w[k + 1]);
                wk3r      = csc3 * (wd3r + w[k + 2]);
                wk3i      = csc3 * (wd3i + w[k + 3]);
                wd1r      = w[k];
                wd1i      = w[k + 1];
                wd3r      = w[k + 2];
                wd3i      = w[k + 3];
                j1        = j + m;
                j2        = j1 + m;
                j3        = j2 + m;
                x0r       = a[j] + a[j2];
                x0i       = -a[j + 1] - a[j2 + 1];
                x1r       = a[j] - a[j2];
                x1i       = -a[j + 1] + a[j2 + 1];
                y0r       = a[j + 2] + a[j2 + 2];
                y0i       = -a[j + 3] - a[j2 + 3];
                y1r       = a[j + 2] - a[j2 + 2];
                y1i       = -a[j + 3] + a[j2 + 3];
                x2r       = a[j1] + a[j3];
                x2i       = a[j1 + 1] + a[j3 + 1];
                x3r       = a[j1] - a[j3];
                x3i       = a[j1 + 1] - a[j3 + 1];
                y2r       = a[j1 + 2] + a[j3 + 2];
                y2i       = a[j1 + 3] + a[j3 + 3];
                y3r       = a[j1 + 2] - a[j3 + 2];
                y3i       = a[j1 + 3] - a[j3 + 3];
                a[j]      = x0r + x2r;
                a[j + 1]  = x0i - x2i;
                a[j + 2]  = y0r + y2r;
                a[j + 3]  = y0i - y2i;
                a[j1]     = x0r - x2r;
                a[j1 + 1] = x0i + x2i;
                a[j1 + 2] = y0r - y2r;
                a[j1 + 3] = y0i + y2i;
                x0r       = x1r + x3i;
                x0i       = x1i + x3r;
                a[j2]     = wk1r * x0r - wk1i * x0i;
                a[j2 + 1] = wk1r * x0i + wk1i * x0r;
                x0r       = y1r + y3i;
                x0i       = y1i + y3r;
                a[j2 + 2] = wd1r * x0r - wd1i * x0i;
                a[j2 + 3] = wd1r * x0i + wd1i * x0r;
                x0r       = x1r - x3i;
                x0i       = x1i - x3r;
                a[j3]     = wk3r * x0r + wk3i * x0i;
                a[j3 + 1] = wk3r * x0i - wk3i * x0r;
                x0r       = y1r - y3i;
                x0i       = y1i - y3r;
                a[j3 + 2] = wd3r * x0r + wd3i * x0i;
                a[j3 + 3] = wd3r * x0i - wd3i * x0r;
                j0        = m - j;
                j1        = j0 + m;
                j2        = j1 + m;
                j3        = j2 + m;
                x0r       = a[j0] + a[j2];
                x0i       = -a[j0 + 1] - a[j2 + 1];
                x1r       = a[j0] - a[j2];
                x1i       = -a[j0 + 1] + a[j2 + 1];
                y0r       = a[j0 - 2] + a[j2 - 2];
                y0i       = -a[j0 - 1] - a[j2 - 1];
                y1r       = a[j0 - 2] - a[j2 - 2];
                y1i       = -a[j0 - 1] + a[j2 - 1];
                x2r       = a[j1] + a[j3];
                x2i       = a[j1 + 1] + a[j3 + 1];
                x3r       = a[j1] - a[j3];
                x3i       = a[j1 + 1] - a[j3 + 1];
                y2r       = a[j1 - 2] + a[j3 - 2];
                y2i       = a[j1 - 1] + a[j3 - 1];
                y3r       = a[j1 - 2] - a[j3 - 2];
                y3i       = a[j1 - 1] - a[j3 - 1];
                a[j0]     = x0r + x2r;
                a[j0 + 1] = x0i - x2i;
                a[j0 - 2] = y0r + y2r;
                a[j0 - 1] = y0i - y2i;
                a[j1]     = x0r - x2r;
                a[j1 + 1] = x0i + x2i;
                a[j1 - 2] = y0r - y2r;
                a[j1 - 1] = y0i + y2i;
                x0r       = x1r + x3i;
                x0i       = x1i + x3r;
                a[j2]     = wk1i * x0r - wk1r * x0i;
                a[j2 + 1] = wk1i * x0i + wk1r * x0r;
                x0r       = y1r + y3i;
                x0i       = y1i + y3r;
                a[j2 - 2] = wd1i * x0r - wd1r * x0i;
                a[j2 - 1] = wd1i * x0i + wd1r * x0r;
                x0r       = x1r - x3i;
                x0i       = x1i - x3r;
                a[j3]     = wk3i * x0r + wk3r * x0i;
                a[j3 + 1] = wk3i * x0i - wk3r * x0r;
                x0r       = y1r - y3i;
                x0i       = y1i - y3r;
                a[j3 - 2] = wd3i * x0r + wd3r * x0i;
                a[j3 - 1] = wd3i * x0i - wd3r * x0r;
            }
            wk1r      = csc1 * (wd1r + wn4r);
            wk1i      = csc1 * (wd1i + wn4r);
            wk3r      = csc3 * (wd3r - wn4r);
            wk3i      = csc3 * (wd3i - wn4r);
            j0        = mh;
            j1        = j0 + m;
            j2        = j1 + m;
            j3        = j2 + m;
            x0r       = a[j0 - 2] + a[j2 - 2];
            x0i       = -a[j0 - 1] - a[j2 - 1];
            x1r       = a[j0 - 2] - a[j2 - 2];
            x1i       = -a[j0 - 1] + a[j2 - 1];
            x2r       = a[j1 - 2] + a[j3 - 2];
            x2i       = a[j1 - 1] + a[j3 - 1];
            x3r       = a[j1 - 2] - a[j3 - 2];
            x3i       = a[j1 - 1] - a[j3 - 1];
            a[j0 - 2] = x0r + x2r;
            a[j0 - 1] = x0i - x2i;
            a[j1 - 2] = x0r - x2r;
            a[j1 - 1] = x0i + x2i;
            x0r       = x1r + x3i;
            x0i       = x1i + x3r;
            a[j2 - 2] = wk1r * x0r - wk1i * x0i;
            a[j2 - 1] = wk1r * x0i + wk1i * x0r;
            x0r       = x1r - x3i;
            x0i       = x1i - x3r;
            a[j3 - 2] = wk3r * x0r + wk3i * x0i;
            a[j3 - 1] = wk3r * x0i - wk3i * x0r;
            x0r       = a[j0] + a[j2];
            x0i       = -a[j0 + 1] - a[j2 + 1];
            x1r       = a[j0] - a[j2];
            x1i       = -a[j0 + 1] + a[j2 + 1];
            x2r       = a[j1] + a[j3];
            x2i       = a[j1 + 1] + a[j3 + 1];
            x3r       = a[j1] - a[j3];
            x3i       = a[j1 + 1] - a[j3 + 1];
            a[j0]     = x0r + x2r;
            a[j0 + 1] = x0i - x2i;
            a[j1]     = x0r - x2r;
            a[j1 + 1] = x0i + x2i;
            x0r       = x1r + x3i;
            x0i       = x1i + x3r;
            a[j2]     = wn4r * (x0r - x0i);
            a[j2 + 1] = wn4r * (x0i + x0r);
            x0r       = x1r - x3i;
            x0i       = x1i - x3r;
            a[j3]     = -wn4r * (x0r + x0i);
            a[j3 + 1] = -wn4r * (x0i - x0r);
            x0r       = a[j0 + 2] + a[j2 + 2];
            x0i       = -a[j0 + 3] - a[j2 + 3];
            x1r       = a[j0 + 2] - a[j2 + 2];
            x1i       = -a[j0 + 3] + a[j2 + 3];
            x2r       = a[j1 + 2] + a[j3 + 2];
            x2i       = a[j1 + 3] + a[j3 + 3];
            x3r       = a[j1 + 2] - a[j3 + 2];
            x3i       = a[j1 + 3] - a[j3 + 3];
            a[j0 + 2] = x0r + x2r;
            a[j0 + 3] = x0i - x2i;
            a[j1 + 2] = x0r - x2r;
            a[j1 + 3] = x0i + x2i;
            x0r       = x1r + x3i;
            x0i       = x1i + x3r;
            a[j2 + 2] = wk1i * x0r - wk1r * x0i;
            a[j2 + 3] = wk1i * x0i + wk1r * x0r;
            x0r       = x1r - x3i;
            x0i       = x1i - x3r;
            a[j3 + 2] = wk3i * x0r + wk3r * x0i;
            a[j3 + 3] = wk3i * x0i - wk3r * x0r;
        }

        static void cftrec4(int n, Sample* a, int nw, const Sample* w) noexcept {
            int isplt, j, k, m;

            m = n;
            while (m > 512) {
                m >>= 2;
                cftmdl1(m, &a[n - m], &w[nw - (m >> 1)]);
            }
            cftleaf(m, 1, &a[n - m], nw, w);
            k = 0;
            for (j = n - m; j > 0; j -= m) {
                k++;
                isplt = cfttree(m, j, k, a, nw, w);
                cftleaf(m, isplt, &a[j - m], nw, w);
            }
        }

        static int cfttree(int n, int j, int k, Sample* a, int nw, const Sample* w) noexcept {
            int i, isplt, m;

            if ((k & 3) != 0) {
                isplt = k & 1;
                if (isplt != 0) {
                    cftmdl1(n, &a[j - n], &w[nw - (n >> 1)]);
                }
                else {
                    cftmdl2(n, &a[j - n], &w[nw - n]);
                }
            }
            else {
                m = n;
                for (i = k; (i & 3) == 0; i >>= 2) {
                    m <<= 2;
                }
                isplt = i & 1;
                if (isplt != 0) {
                    while (m > 128) {
                        cftmdl1(m, &a[j - m], &w[nw - (m >> 1)]);
                        m >>= 2;
                    }
                }
                else {
                    while (m > 128) {
                        cftmdl2(m, &a[j - m], &w[nw - m]);
                        m >>= 2;
                    }
                }
            }
            return isplt;
        }

        static void cftleaf(int n, int isplt, Sample* a, int nw, const Sample* w) noexcept {
            if (n == 512) {
                cftmdl1(128, a, &w[nw - 64]);
                cftf161(a, &w[nw - 8]);
                cftf162(&a[32], &w[nw - 32]);
                cftf161(&a[64], &w[nw - 8]);
                cftf161(&a[96], &w[nw - 8]);
                cftmdl2(128, &a[128], &w[nw - 128]);
                cftf161(&a[128], &w[nw - 8]);
                cftf162(&a[160], &w[nw - 32]);
                cftf161(&a[192], &w[nw - 8]);
                cftf162(&a[224], &w[nw - 32]);
                cftmdl1(128, &a[256], &w[nw - 64]);
                cftf161(&a[256], &w[nw - 8]);
                cftf162(&a[288], &w[nw - 32]);
                cftf161(&a[320], &w[nw - 8]);
                cftf161(&a[352], &w[nw - 8]);
                if (isplt != 0) {
                    cftmdl1(128, &a[384], &w[nw - 64]);
                    cftf161(&a[480], &w[nw - 8]);
                }
                else {
                    cftmdl2(128, &a[384], &w[nw - 128]);
                    cftf162(&a[480], &w[nw - 32]);
                }
                cftf161(&a[384], &w[nw - 8]);
                cftf162(&a[416], &w[nw - 32]);
                cftf161(&a[448], &w[nw - 8]);
            }
            else {
                cftmdl1(64, a, &w[nw - 32]);
                cftf081(a, &w[nw - 8]);
                cftf082(&a[16], &w[nw - 8]);
                cftf081(&a[32], &w[nw - 8]);
                cftf081(&a[48], &w[nw - 8]);
                cftmdl2(64, &a[64], &w[nw - 64]);
                cftf081(&a[64], &w[nw - 8]);
                cftf082(&a[80], &w[nw - 8]);
                cftf081(&a[96], &w[nw - 8]);
                cftf082(&a[112], &w[nw - 8]);
                cftmdl1(64, &a[128], &w[nw - 32]);
                cftf081(&a[128], &w[nw - 8]);
                cftf082(&a[144], &w[nw - 8]);
                cftf081(&a[160], &w[nw - 8]);
                cftf081(&a[176], &w[nw - 8]);
                if (isplt != 0) {
                    cftmdl1(64, &a[192], &w[nw - 32]);
                    cftf081(&a[240], &w[nw - 8]);
                }
                else {
                    cftmdl2(64, &a[192], &w[nw - 64]);
                    cftf082(&a[240], &w[nw - 8]);
                }
                cftf081(&a[192], &w[nw - 8]);
                cftf082(&a[208], &w[nw - 8]);
                cftf081(&a[224], &w[nw - 8]);
            }
        }

        static void cftmdl1(int n, Sample* a, const Sample* w) noexcept {
            int    j, j0, j1, j2, j3, k, m, mh;
            Sample wn4r, wk1r, wk1i, wk3r, wk3i;
            Sample x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i;

            mh        = n >> 3;
            m         = 2 * mh;
            j1        = m;
            j2        = j1 + m;
            j3        = j2 + m;
            x0r       = a[0] + a[j2];
            x0i       = a[1] + a[j2 + 1];
            x1r       = a[0] - a[j2];
            x1i       = a[1] - a[j2 + 1];
            x2r       = a[j1] + a[j3];
            x2i       = a[j1 + 1] + a[j3 + 1];
            x3r       = a[j1] - a[j3];
            x3i       = a[j1 + 1] - a[j3 + 1];
            a[0]      = x0r + x2r;
            a[1]      = x0i + x2i;
            a[j1]     = x0r - x2r;
            a[j1 + 1] = x0i - x2i;
            a[j2]     = x1r - x3i;
            a[j2 + 1] = x1i + x3r;
            a[j3]     = x1r + x3i;
            a[j3 + 1] = x1i - x3r;
            wn4r      = w[1];
            k         = 0;
            for (j = 2; j < mh; j += 2) {
                k += 4;
                wk1r      = w[k];
                wk1i      = w[k + 1];
                wk3r      = w[k + 2];
                wk3i      = w[k + 3];
                j1        = j + m;
                j2        = j1 + m;
                j3        = j2 + m;
                x0r       = a[j] + a[j2];
                x0i       = a[j + 1] + a[j2 + 1];
                x1r       = a[j] - a[j2];
                x1i       = a[j + 1] - a[j2 + 1];
                x2r       = a[j1] + a[j3];
                x2i       = a[j1 + 1] + a[j3 + 1];
                x3r       = a[j1] - a[j3];
                x3i       = a[j1 + 1] - a[j3 + 1];
                a[j]      = x0r + x2r;
                a[j + 1]  = x0i + x2i;
                a[j1]     = x0r - x2r;
                a[j1 + 1] = x0i - x2i;
                x0r       = x1r - x3i;
                x0i       = x1i + x3r;
                a[j2]     = wk1r * x0r - wk1i * x0i;
                a[j2 + 1] = wk1r * x0i + wk1i * x0r;
                x0r       = x1r + x3i;
                x0i       = x1i - x3r;
                a[j3]     = wk3r * x0r + wk3i * x0i;
                a[j3 + 1] = wk3r * x0i - wk3i * x0r;
                j0        = m - j;
                j1        = j0 + m;
                j2        = j1 + m;
                j3        = j2 + m;
                x0r       = a[j0] + a[j2];
                x0i       = a[j0 + 1] + a[j2 + 1];
                x1r       = a[j0] - a[j2];
                x1i       = a[j0 + 1] - a[j2 + 1];
                x2r       = a[j1] + a[j3];
                x2i       = a[j1 + 1] + a[j3 + 1];
                x3r       = a[j1] - a[j3];
                x3i       = a[j1 + 1] - a[j3 + 1];
                a[j0]     = x0r + x2r;
                a[j0 + 1] = x0i + x2i;
                a[j1]     = x0r - x2r;
                a[j1 + 1] = x0i - x2i;
                x0r       = x1r - x3i;
                x0i       = x1i + x3r;
                a[j2]     = wk1i * x0r - wk1r * x0i;
                a[j2 + 1] = wk1i * x0i + wk1r * x0r;
                x0r       = x1r + x3i;
                x0i       = x1i - x3r;
                a[j3]     = wk3i * x0r + wk3r * x0i;
                a[j3 + 1] = wk3i * x0i - wk3r * x0r;
            }
            j0        = mh;
            j1        = j0 + m;
            j2        = j1 + m;
            j3        = j2 + m;
            x0r       = a[j0] + a[j2];
            x0i       = a[j0 + 1] + a[j2 + 1];
            x1r       = a[j0] - a[j2];
            x1i       = a[j0 + 1] - a[j2 + 1];
            x2r       = a[j1] + a[j3];
            x2i       = a[j1 + 1] + a[j3 + 1];
            x3r       = a[j1] - a[j3];
            x3i       = a[j1 + 1] - a[j3 + 1];
            a[j0]     = x0r + x2r;
            a[j0 + 1] = x0i + x2i;
            a[j1]     = x0r - x2r;
            a[j1 + 1] = x0i - x2i;
            x0r       = x1r - x3i;
            x0i       = x1i + x3r;
            a[j2]     = wn4r * (x0r - x0i);
            a[j2 + 1] = wn4r * (x0i + x0r);
            x0r       = x1r + x3i;
            x0i       = x1i - x3r;
            a[j3]     = -wn4r * (x0r + x0i);
            a[j3 + 1] = -wn4r * (x0i - x0r);
        }

        static void cftmdl2(int n, Sample* a, const Sample* w) noexcept {
            int    j, j0, j1, j2, j3, k, kr, m, mh;
            Sample wn4r, wk1r, wk1i, wk3r, wk3i, wd1r, wd1i, wd3r, wd3i;
            Sample x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i, y0r, y0i, y2r, y2i;

            mh        = n >> 3;
            m         = 2 * mh;
            wn4r      = w[1];
            j1        = m;
            j2        = j1 + m;
            j3        = j2 + m;
            x0r       = a[0] - a[j2 + 1];
            x0i       = a[1] + a[j2];
            x1r       = a[0] + a[j2 + 1];
            x1i       = a[1] - a[j2];
            x2r       = a[j1] - a[j3 + 1];
            x2i       = a[j1 + 1] + a[j3];
            x3r       = a[j1] + a[j3 + 1];
            x3i       = a[j1 + 1] - a[j3];
            y0r       = wn4r * (x2r - x2i);
            y0i       = wn4r * (x2i + x2r);
            a[0]      = x0r + y0r;
            a[1]      = x0i + y0i;
            a[j1]     = x0r - y0r;
            a[j1 + 1] = x0i - y0i;
            y0r       = wn4r * (x3r - x3i);
            y0i       = wn4r * (x3i + x3r);
            a[j2]     = x1r - y0i;
            a[j2 + 1] = x1i + y0r;
            a[j3]     = x1r + y0i;
            a[j3 + 1] = x1i - y0r;
            k         = 0;
            kr        = 2 * m;
            for (j = 2; j < mh; j += 2) {
                k += 4;
                wk1r = w[k];
                wk1i = w[k + 1];
                wk3r = w[k + 2];
                wk3i = w[k + 3];
                kr -= 4;
                wd1i      = w[kr];
                wd1r      = w[kr + 1];
                wd3i      = w[kr + 2];
                wd3r      = w[kr + 3];
                j1        = j + m;
                j2        = j1 + m;
                j3        = j2 + m;
                x0r       = a[j] - a[j2 + 1];
                x0i       = a[j + 1] + a[j2];
                x1r       = a[j] + a[j2 + 1];
                x1i       = a[j + 1] - a[j2];
                x2r       = a[j1] - a[j3 + 1];
                x2i       = a[j1 + 1] + a[j3];
                x3r       = a[j1] + a[j3 + 1];
                x3i       = a[j1 + 1] - a[j3];
                y0r       = wk1r * x0r - wk1i * x0i;
                y0i       = wk1r * x0i + wk1i * x0r;
                y2r       = wd1r * x2r - wd1i * x2i;
                y2i       = wd1r * x2i + wd1i * x2r;
                a[j]      = y0r + y2r;
                a[j + 1]  = y0i + y2i;
                a[j1]     = y0r - y2r;
                a[j1 + 1] = y0i - y2i;
                y0r       = wk3r * x1r + wk3i * x1i;
                y0i       = wk3r * x1i - wk3i * x1r;
                y2r       = wd3r * x3r + wd3i * x3i;
                y2i       = wd3r * x3i - wd3i * x3r;
                a[j2]     = y0r + y2r;
                a[j2 + 1] = y0i + y2i;
                a[j3]     = y0r - y2r;
                a[j3 + 1] = y0i - y2i;
                j0        = m - j;
                j1        = j0 + m;
                j2        = j1 + m;
                j3        = j2 + m;
                x0r       = a[j0] - a[j2 + 1];
                x0i       = a[j0 + 1] + a[j2];
                x1r       = a[j0] + a[j2 + 1];
                x1i       = a[j0 + 1] - a[j2];
                x2r       = a[j1] - a[j3 + 1];
                x2i       = a[j1 + 1] + a[j3];
                x3r       = a[j1] + a[j3 + 1];
                x3i       = a[j1 + 1] - a[j3];
                y0r       = wd1i * x0r - wd1r * x0i;
                y0i       = wd1i * x0i + wd1r * x0r;
                y2r       = wk1i * x2r - wk1r * x2i;
                y2i       = wk1i * x2i + wk1r * x2r;
                a[j0]     = y0r + y2r;
                a[j0 + 1] = y0i + y2i;
                a[j1]     = y0r - y2r;
                a[j1 + 1] = y0i - y2i;
                y0r       = wd3i * x1r + wd3r * x1i;
                y0i       = wd3i * x1i - wd3r * x1r;
                y2r       = wk3i * x3r + wk3r * x3i;
                y2i       = wk3i * x3i - wk3r * x3r;
                a[j2]     = y0r + y2r;
                a[j2 + 1] = y0i + y2i;
                a[j3]     = y0r - y2r;
                a[j3 + 1] = y0i - y2i;
            }
            wk1r      = w[m];
            wk1i      = w[m + 1];
            j0        = mh;
            j1        = j0 + m;
            j2        = j1 + m;
            j3        = j2 + m;
            x0r       = a[j0] - a[j2 + 1];
            x0i       = a[j0 + 1] + a[j2];
            x1r       = a[j0] + a[j2 + 1];
            x1i       = a[j0 + 1] - a[j2];
            x2r       = a[j1] - a[j3 + 1];
            x2i       = a[j1 + 1] + a[j3];
            x3r       = a[j1] + a[j3 + 1];
            x3i       = a[j1 + 1] - a[j3];
            y0r       = wk1r * x0r - wk1i * x0i;
            y0i       = wk1r * x0i + wk1i * x0r;
            y2r       = wk1i * x2r - wk1r * x2i;
            y2i       = wk1i * x2i + wk1r * x2r;
            a[j0]     = y0r + y2r;
            a[j0 + 1] = y0i + y2i;
            a[j1]     = y0r - y2r;
            a[j1 + 1] = y0i - y2i;
            y0r       = wk1i * x1r - wk1r * x1i;
            y0i       = wk1i * x1i + wk1r * x1r;
            y2r       = wk1r * x3r - wk1i * x3i;
            y2i       = wk1r * x3i + wk1i * x3r;
            a[j2]     = y0r - y2r;
            a[j2 + 1] = y0i - y2i;
            a[j3]     = y0r + y2r;
            a[j3 + 1] = y0i + y2i;
        }

        static void cftfx41(int n, Sample* a, int nw, const Sample* w) noexcept {
            if (n == 128) {
                cftf161(a, &w[nw - 8]);
                cftf162(&a[32], &w[nw - 32]);
                cftf161(&a[64], &w[nw - 8]);
                cftf161(&a[96], &w[nw - 8]);
            }
            else {
                cftf081(a, &w[nw - 8]);
                cftf082(&a[16], &w[nw - 8]);
                cftf081(&a[32], &w[nw - 8]);
                cftf081(&a[48], &w[nw - 8]);
            }
        }

        static void cftf161(Sample* a, const Sample* w) noexcept {
            Sample wn4r, wk1r, wk1i, x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i, y0r, y0i, y1r, y1i, y2r, y2i, y3r, y3i,
                y4r, y4i, y5r, y5i, y6r, y6i, y7r, y7i, y8r, y8i, y9r, y9i, y10r, y10i, y11r, y11i, y12r, y12i, y13r,
                y13i, y14r, y14i, y15r, y15i;

            wn4r  = w[1];
            wk1r  = w[2];
            wk1i  = w[3];
            x0r   = a[0] + a[16];
            x0i   = a[1] + a[17];
            x1r   = a[0] - a[16];
            x1i   = a[1] - a[17];
            x2r   = a[8] + a[24];
            x2i   = a[9] + a[25];
            x3r   = a[8] - a[24];
            x3i   = a[9] - a[25];
            y0r   = x0r + x2r;
            y0i   = x0i + x2i;
            y4r   = x0r - x2r;
            y4i   = x0i - x2i;
            y8r   = x1r - x3i;
            y8i   = x1i + x3r;
            y12r  = x1r + x3i;
            y12i  = x1i - x3r;
            x0r   = a[2] + a[18];
            x0i   = a[3] + a[19];
            x1r   = a[2] - a[18];
            x1i   = a[3] - a[19];
            x2r   = a[10] + a[26];
            x2i   = a[11] + a[27];
            x3r   = a[10] - a[26];
            x3i   = a[11] - a[27];
            y1r   = x0r + x2r;
            y1i   = x0i + x2i;
            y5r   = x0r - x2r;
            y5i   = x0i - x2i;
            x0r   = x1r - x3i;
            x0i   = x1i + x3r;
            y9r   = wk1r * x0r - wk1i * x0i;
            y9i   = wk1r * x0i + wk1i * x0r;
            x0r   = x1r + x3i;
            x0i   = x1i - x3r;
            y13r  = wk1i * x0r - wk1r * x0i;
            y13i  = wk1i * x0i + wk1r * x0r;
            x0r   = a[4] + a[20];
            x0i   = a[5] + a[21];
            x1r   = a[4] - a[20];
            x1i   = a[5] - a[21];
            x2r   = a[12] + a[28];
            x2i   = a[13] + a[29];
            x3r   = a[12] - a[28];
            x3i   = a[13] - a[29];
            y2r   = x0r + x2r;
            y2i   = x0i + x2i;
            y6r   = x0r - x2r;
            y6i   = x0i - x2i;
            x0r   = x1r - x3i;
            x0i   = x1i + x3r;
            y10r  = wn4r * (x0r - x0i);
            y10i  = wn4r * (x0i + x0r);
            x0r   = x1r + x3i;
            x0i   = x1i - x3r;
            y14r  = wn4r * (x0r + x0i);
            y14i  = wn4r * (x0i - x0r);
            x0r   = a[6] + a[22];
            x0i   = a[7] + a[23];
            x1r   = a[6] - a[22];
            x1i   = a[7] - a[23];
            x2r   = a[14] + a[30];
            x2i   = a[15] + a[31];
            x3r   = a[14] - a[30];
            x3i   = a[15] - a[31];
            y3r   = x0r + x2r;
            y3i   = x0i + x2i;
            y7r   = x0r - x2r;
            y7i   = x0i - x2i;
            x0r   = x1r - x3i;
            x0i   = x1i + x3r;
            y11r  = wk1i * x0r - wk1r * x0i;
            y11i  = wk1i * x0i + wk1r * x0r;
            x0r   = x1r + x3i;
            x0i   = x1i - x3r;
            y15r  = wk1r * x0r - wk1i * x0i;
            y15i  = wk1r * x0i + wk1i * x0r;
            x0r   = y12r - y14r;
            x0i   = y12i - y14i;
            x1r   = y12r + y14r;
            x1i   = y12i + y14i;
            x2r   = y13r - y15r;
            x2i   = y13i - y15i;
            x3r   = y13r + y15r;
            x3i   = y13i + y15i;
            a[24] = x0r + x2r;
            a[25] = x0i + x2i;
            a[26] = x0r - x2r;
            a[27] = x0i - x2i;
            a[28] = x1r - x3i;
            a[29] = x1i + x3r;
            a[30] = x1r + x3i;
            a[31] = x1i - x3r;
            x0r   = y8r + y10r;
            x0i   = y8i + y10i;
            x1r   = y8r - y10r;
            x1i   = y8i - y10i;
            x2r   = y9r + y11r;
            x2i   = y9i + y11i;
            x3r   = y9r - y11r;
            x3i   = y9i - y11i;
            a[16] = x0r + x2r;
            a[17] = x0i + x2i;
            a[18] = x0r - x2r;
            a[19] = x0i - x2i;
            a[20] = x1r - x3i;
            a[21] = x1i + x3r;
            a[22] = x1r + x3i;
            a[23] = x1i - x3r;
            x0r   = y5r - y7i;
            x0i   = y5i + y7r;
            x2r   = wn4r * (x0r - x0i);
            x2i   = wn4r * (x0i + x0r);
            x0r   = y5r + y7i;
            x0i   = y5i - y7r;
            x3r   = wn4r * (x0r - x0i);
            x3i   = wn4r * (x0i + x0r);
            x0r   = y4r - y6i;
            x0i   = y4i + y6r;
            x1r   = y4r + y6i;
            x1i   = y4i - y6r;
            a[8]  = x0r + x2r;
            a[9]  = x0i + x2i;
            a[10] = x0r - x2r;
            a[11] = x0i - x2i;
            a[12] = x1r - x3i;
            a[13] = x1i + x3r;
            a[14] = x1r + x3i;
            a[15] = x1i - x3r;
            x0r   = y0r + y2r;
            x0i   = y0i + y2i;
            x1r   = y0r - y2r;
            x1i   = y0i - y2i;
            x2r   = y1r + y3r;
            x2i   = y1i + y3i;
            x3r   = y1r - y3r;
            x3i   = y1i - y3i;
            a[0]  = x0r + x2r;
            a[1]  = x0i + x2i;
            a[2]  = x0r - x2r;
            a[3]  = x0i - x2i;
            a[4]  = x1r - x3i;
            a[5]  = x1i + x3r;
            a[6]  = x1r + x3i;
            a[7]  = x1i - x3r;
        }

        static void cftf162(Sample* a, const Sample* w) noexcept {
            Sample wn4r, wk1r, wk1i, wk2r, wk2i, wk3r, wk3i, x0r, x0i, x1r, x1i, x2r, x2i, y0r, y0i, y1r, y1i, y2r, y2i,
                y3r, y3i, y4r, y4i, y5r, y5i, y6r, y6i, y7r, y7i, y8r, y8i, y9r, y9i, y10r, y10i, y11r, y11i, y12r,
                y12i, y13r, y13i, y14r, y14i, y15r, y15i;

            wn4r  = w[1];
            wk1r  = w[4];
            wk1i  = w[5];
            wk3r  = w[6];
            wk3i  = -w[7];
            wk2r  = w[8];
            wk2i  = w[9];
            x1r   = a[0] - a[17];
            x1i   = a[1] + a[16];
            x0r   = a[8] - a[25];
            x0i   = a[9] + a[24];
            x2r   = wn4r * (x0r - x0i);
            x2i   = wn4r * (x0i + x0r);
            y0r   = x1r + x2r;
            y0i   = x1i + x2i;
            y4r   = x1r - x2r;
            y4i   = x1i - x2i;
            x1r   = a[0] + a[17];
            x1i   = a[1] - a[16];
            x0r   = a[8] + a[25];
            x0i   = a[9] - a[24];
            x2r   = wn4r * (x0r - x0i);
            x2i   = wn4r * (x0i + x0r);
            y8r   = x1r - x2i;
            y8i   = x1i + x2r;
            y12r  = x1r + x2i;
            y12i  = x1i - x2r;
            x0r   = a[2] - a[19];
            x0i   = a[3] + a[18];
            x1r   = wk1r * x0r - wk1i * x0i;
            x1i   = wk1r * x0i + wk1i * x0r;
            x0r   = a[10] - a[27];
            x0i   = a[11] + a[26];
            x2r   = wk3i * x0r - wk3r * x0i;
            x2i   = wk3i * x0i + wk3r * x0r;
            y1r   = x1r + x2r;
            y1i   = x1i + x2i;
            y5r   = x1r - x2r;
            y5i   = x1i - x2i;
            x0r   = a[2] + a[19];
            x0i   = a[3] - a[18];
            x1r   = wk3r * x0r - wk3i * x0i;
            x1i   = wk3r * x0i + wk3i * x0r;
            x0r   = a[10] + a[27];
            x0i   = a[11] - a[26];
            x2r   = wk1r * x0r + wk1i * x0i;
            x2i   = wk1r * x0i - wk1i * x0r;
            y9r   = x1r - x2r;
            y9i   = x1i - x2i;
            y13r  = x1r + x2r;
            y13i  = x1i + x2i;
            x0r   = a[4] - a[21];
            x0i   = a[5] + a[20];
            x1r   = wk2r * x0r - wk2i * x0i;
            x1i   = wk2r * x0i + wk2i * x0r;
            x0r   = a[12] - a[29];
            x0i   = a[13] + a[28];
            x2r   = wk2i * x0r - wk2r * x0i;
            x2i   = wk2i * x0i + wk2r * x0r;
            y2r   = x1r + x2r;
            y2i   = x1i + x2i;
            y6r   = x1r - x2r;
            y6i   = x1i - x2i;
            x0r   = a[4] + a[21];
            x0i   = a[5] - a[20];
            x1r   = wk2i * x0r - wk2r * x0i;
            x1i   = wk2i * x0i + wk2r * x0r;
            x0r   = a[12] + a[29];
            x0i   = a[13] - a[28];
            x2r   = wk2r * x0r - wk2i * x0i;
            x2i   = wk2r * x0i + wk2i * x0r;
            y10r  = x1r - x2r;
            y10i  = x1i - x2i;
            y14r  = x1r + x2r;
            y14i  = x1i + x2i;
            x0r   = a[6] - a[23];
            x0i   = a[7] + a[22];
            x1r   = wk3r * x0r - wk3i * x0i;
            x1i   = wk3r * x0i + wk3i * x0r;
            x0r   = a[14] - a[31];
            x0i   = a[15] + a[30];
            x2r   = wk1i * x0r - wk1r * x0i;
            x2i   = wk1i * x0i + wk1r * x0r;
            y3r   = x1r + x2r;
            y3i   = x1i + x2i;
            y7r   = x1r - x2r;
            y7i   = x1i - x2i;
            x0r   = a[6] + a[23];
            x0i   = a[7] - a[22];
            x1r   = wk1i * x0r + wk1r * x0i;
            x1i   = wk1i * x0i - wk1r * x0r;
            x0r   = a[14] + a[31];
            x0i   = a[15] - a[30];
            x2r   = wk3i * x0r - wk3r * x0i;
            x2i   = wk3i * x0i + wk3r * x0r;
            y11r  = x1r + x2r;
            y11i  = x1i + x2i;
            y15r  = x1r - x2r;
            y15i  = x1i - x2i;
            x1r   = y0r + y2r;
            x1i   = y0i + y2i;
            x2r   = y1r + y3r;
            x2i   = y1i + y3i;
            a[0]  = x1r + x2r;
            a[1]  = x1i + x2i;
            a[2]  = x1r - x2r;
            a[3]  = x1i - x2i;
            x1r   = y0r - y2r;
            x1i   = y0i - y2i;
            x2r   = y1r - y3r;
            x2i   = y1i - y3i;
            a[4]  = x1r - x2i;
            a[5]  = x1i + x2r;
            a[6]  = x1r + x2i;
            a[7]  = x1i - x2r;
            x1r   = y4r - y6i;
            x1i   = y4i + y6r;
            x0r   = y5r - y7i;
            x0i   = y5i + y7r;
            x2r   = wn4r * (x0r - x0i);
            x2i   = wn4r * (x0i + x0r);
            a[8]  = x1r + x2r;
            a[9]  = x1i + x2i;
            a[10] = x1r - x2r;
            a[11] = x1i - x2i;
            x1r   = y4r + y6i;
            x1i   = y4i - y6r;
            x0r   = y5r + y7i;
            x0i   = y5i - y7r;
            x2r   = wn4r * (x0r - x0i);
            x2i   = wn4r * (x0i + x0r);
            a[12] = x1r - x2i;
            a[13] = x1i + x2r;
            a[14] = x1r + x2i;
            a[15] = x1i - x2r;
            x1r   = y8r + y10r;
            x1i   = y8i + y10i;
            x2r   = y9r - y11r;
            x2i   = y9i - y11i;
            a[16] = x1r + x2r;
            a[17] = x1i + x2i;
            a[18] = x1r - x2r;
            a[19] = x1i - x2i;
            x1r   = y8r - y10r;
            x1i   = y8i - y10i;
            x2r   = y9r + y11r;
            x2i   = y9i + y11i;
            a[20] = x1r - x2i;
            a[21] = x1i + x2r;
            a[22] = x1r + x2i;
            a[23] = x1i - x2r;
            x1r   = y12r - y14i;
            x1i   = y12i + y14r;
            x0r   = y13r + y15i;
            x0i   = y13i - y15r;
            x2r   = wn4r * (x0r - x0i);
            x2i   = wn4r * (x0i + x0r);
            a[24] = x1r + x2r;
            a[25] = x1i + x2i;
            a[26] = x1r - x2r;
            a[27] = x1i - x2i;
            x1r   = y12r + y14i;
            x1i   = y12i - y14r;
            x0r   = y13r - y15i;
            x0i   = y13i + y15r;
            x2r   = wn4r * (x0r - x0i);
            x2i   = wn4r * (x0i + x0r);
            a[28] = x1r - x2i;
            a[29] = x1i + x2r;
            a[30] = x1r + x2i;
            a[31] = x1i - x2r;
        }

        static void cftf081(Sample* a, const Sample* w) noexcept {
            Sample wn4r, x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i, y0r, y0i, y1r, y1i, y2r, y2i, y3r, y3i, y4r, y4i, y5r,
                y5i, y6r, y6i, y7r, y7i;

            wn4r  = w[1];
            x0r   = a[0] + a[8];
            x0i   = a[1] + a[9];
            x1r   = a[0] - a[8];
            x1i   = a[1] - a[9];
            x2r   = a[4] + a[12];
            x2i   = a[5] + a[13];
            x3r   = a[4] - a[12];
            x3i   = a[5] - a[13];
            y0r   = x0r + x2r;
            y0i   = x0i + x2i;
            y2r   = x0r - x2r;
            y2i   = x0i - x2i;
            y1r   = x1r - x3i;
            y1i   = x1i + x3r;
            y3r   = x1r + x3i;
            y3i   = x1i - x3r;
            x0r   = a[2] + a[10];
            x0i   = a[3] + a[11];
            x1r   = a[2] - a[10];
            x1i   = a[3] - a[11];
            x2r   = a[6] + a[14];
            x2i   = a[7] + a[15];
            x3r   = a[6] - a[14];
            x3i   = a[7] - a[15];
            y4r   = x0r + x2r;
            y4i   = x0i + x2i;
            y6r   = x0r - x2r;
            y6i   = x0i - x2i;
            x0r   = x1r - x3i;
            x0i   = x1i + x3r;
            x2r   = x1r + x3i;
            x2i   = x1i - x3r;
            y5r   = wn4r * (x0r - x0i);
            y5i   = wn4r * (x0r + x0i);
            y7r   = wn4r * (x2r - x2i);
            y7i   = wn4r * (x2r + x2i);
            a[8]  = y1r + y5r;
            a[9]  = y1i + y5i;
            a[10] = y1r - y5r;
            a[11] = y1i - y5i;
            a[12] = y3r - y7i;
            a[13] = y3i + y7r;
            a[14] = y3r + y7i;
            a[15] = y3i - y7r;
            a[0]  = y0r + y4r;
            a[1]  = y0i + y4i;
            a[2]  = y0r - y4r;
            a[3]  = y0i - y4i;
            a[4]  = y2r - y6i;
            a[5]  = y2i + y6r;
            a[6]  = y2r + y6i;
            a[7]  = y2i - y6r;
        }

        static void cftf082(Sample* a, const Sample* w) noexcept {
            Sample wn4r, wk1r, wk1i, x0r, x0i, x1r, x1i, y0r, y0i, y1r, y1i, y2r, y2i, y3r, y3i, y4r, y4i, y5r, y5i,
                y6r, y6i, y7r, y7i;

            wn4r  = w[1];
            wk1r  = w[2];
            wk1i  = w[3];
            y0r   = a[0] - a[9];
            y0i   = a[1] + a[8];
            y1r   = a[0] + a[9];
            y1i   = a[1] - a[8];
            x0r   = a[4] - a[13];
            x0i   = a[5] + a[12];
            y2r   = wn4r * (x0r - x0i);
            y2i   = wn4r * (x0i + x0r);
            x0r   = a[4] + a[13];
            x0i   = a[5] - a[12];
            y3r   = wn4r * (x0r - x0i);
            y3i   = wn4r * (x0i + x0r);
            x0r   = a[2] - a[11];
            x0i   = a[3] + a[10];
            y4r   = wk1r * x0r - wk1i * x0i;
            y4i   = wk1r * x0i + wk1i * x0r;
            x0r   = a[2] + a[11];
            x0i   = a[3] - a[10];
            y5r   = wk1i * x0r - wk1r * x0i;
            y5i   = wk1i * x0i + wk1r * x0r;
            x0r   = a[6] - a[15];
            x0i   = a[7] + a[14];
            y6r   = wk1i * x0r - wk1r * x0i;
            y6i   = wk1i * x0i + wk1r * x0r;
            x0r   = a[6] + a[15];
            x0i   = a[7] - a[14];
            y7r   = wk1r * x0r - wk1i * x0i;
            y7i   = wk1r * x0i + wk1i * x0r;
            x0r   = y0r + y2r;
            x0i   = y0i + y2i;
            x1r   = y4r + y6r;
            x1i   = y4i + y6i;
            a[0]  = x0r + x1r;
            a[1]  = x0i + x1i;
            a[2]  = x0r - x1r;
            a[3]  = x0i - x1i;
            x0r   = y0r - y2r;
            x0i   = y0i - y2i;
            x1r   = y4r - y6r;
            x1i   = y4i - y6i;
            a[4]  = x0r - x1i;
            a[5]  = x0i + x1r;
            a[6]  = x0r + x1i;
            a[7]  = x0i - x1r;
            x0r   = y1r - y3i;
            x0i   = y1i + y3r;
            x1r   = y5r - y7r;
            x1i   = y5i - y7i;
            a[8]  = x0r + x1r;
            a[9]  = x0i + x1i;
            a[10] = x0r - x1r;
            a[11] = x0i - x1i;
            x0r   = y1r + y3i;
            x0i   = y1i - y3r;
            x1r   = y5r + y7r;
            x1i   = y5i + y7i;
            a[12] = x0r - x1i;
            a[13] = x0i + x1r;
            a[14] = x0r + x1i;
            a[15] = x0i - x1r;
        }

        static void cftf040(Sample* a) noexcept {
            Sample x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i;

            x0r  = a[0] + a[4];
            x0i  = a[1] + a[5];
            x1r  = a[0] - a[4];
            x1i  = a[1] - a[5];
            x2r  = a[2] + a[6];
            x2i  = a[3] + a[7];
            x3r  = a[2] - a[6];
            x3i  = a[3] - a[7];
            a[0] = x0r + x2r;
            a[1] = x0i + x2i;
            a[2] = x1r - x3i;
            a[3] = x1i + x3r;
            a[4] = x0r - x2r;
            a[5] = x0i - x2i;
            a[6] = x1r + x3i;
            a[7] = x1i - x3r;
        }

        static void cftb040(Sample* a) noexcept {
            Sample x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i;

            x0r  = a[0] + a[4];
            x0i  = a[1] + a[5];
            x1r  = a[0] - a[4];
            x1i  = a[1] - a[5];
            x2r  = a[2] + a[6];
            x2i  = a[3] + a[7];
            x3r  = a[2] - a[6];
            x3i  = a[3] - a[7];
            a[0] = x0r + x2r;
            a[1] = x0i + x2i;
            a[2] = x1r + x3i;
            a[3] = x1i - x3r;
            a[4] = x0r - x2r;
            a[5] = x0i - x2i;
            a[6] = x1r - x3i;
            a[7] = x1i + x3r;
        }

        static void cftx020(Sample* a) noexcept {
            Sample x0r, x0i;

            x0r = a[0] - a[2];
            x0i = a[1] - a[3];
            a[0] += a[2];
            a[1] += a[3];
            a[2] = x0r;
            a[3] = x0i;
        }

        int              m_size; ///< N, the transform size
        int              m_nw;   ///< length of the complex-stage twiddle table (N/4), as rdft computes it
        int              m_nc;   ///< length of the real post-pass table (N/4; 1 at N = 4, where it is unused)
        std::vector<int> m_ip;   ///< bit-reversal work table: 2 + sqrt(N/2) + 1 entries
        std::vector<Sample>
            m_w; ///< trig tables: m_w[0, m_nw) for the complex stages, m_w[m_nw, N/2) for rftfsub/rftbsub
    };

} // namespace tap::dsp::detail
