/// @file butterworth.h
/// @brief Butterworth low-pass: the even-order design (bilinear, pre-warped) and a 1:1 chain stage over it.
// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// A recursive 1:1 chain stage over the IIR substrate (iir_kernels.h's biquad
// cascade, iir_design.h's sections and group delay, chain.h's nominal
// latency), under SampleRateTap's IIR plan's contracts for a recursive stage
// in a chain: latency I5 and flush I6. Subnormals are the host's (I9,
// subnormals.h): the host sets flush-to-zero; the kernel only reads it, as
// a debug-build precondition.
//
// Design (Oppenheim & Schafer, "Discrete-Time Signal Processing", 3rd ed.,
// 2010, sec. 7.2: the bilinear transform with frequency pre-warping, and
// the Butterworth prototype's pole placement). For order N (even) and cutoff
// fc in cycles per sample (the -3 dB point), K = tan(pi fc); the analogue
// prototype's poles pair into N/2 sections s^2 + alpha_k s + 1 with
// alpha_k = 2 sin((2k + 1) pi / (2N)), and the bilinear map gives each
//
//     b = K^2 [1, 2, 1] / a0,   a1 = 2 (K^2 - 1) / a0,   a2 = (1 - alpha_k K + K^2) / a0,
//     a0 = 1 + alpha_k K + K^2.
//
// Each section has unity gain at DC and a double zero at Nyquist. In exact
// arithmetic each is the cookbook low-pass biquad (Bristow-Johnson, "Audio
// EQ Cookbook") at Q_k = 1 / alpha_k, the standard Butterworth Q set (0.5412
// and 1.3066 at order 4): the same filter, so a cookbook cascade at those Qs
// is this design within its rounding
// (`butterworth.SectionsAreTheCookbookBiquadsAtTheButterworthQs`). The
// sections run in increasing pole radius (lowest Q first), the order the
// plan's I8 requires for fixed point. The design is runtime double, like
// every DspTap designer (kaiser.h, nyquist.h): a few coefficients through
// libm's tan and sin. The plan's constexpr designers (I4) are for its named
// tables; a Butterworth with a run-time cutoff has none.
//
// Contract, as numbers:
//   - Order N even, 2 ... 16. Cutoff k_butterworth_min_cutoff<S> <= fc <=
//     1/2 - k_butterworth_min_cutoff<S> cycles per sample: 1e-5 for double,
//     1e-4 for float (butterworth_supports_cutoff<S>; the constructor
//     rejects anything else through raise<std::invalid_argument>,
//     set_cutoff() states it as its precondition). Float's bound is its
//     poles': each section's stability margin, 1 + a2 - |a1|, is 4 K^2 / a0
//     at the low end (4 / a0 at the high end, symmetrically), and rounding a1
//     and a2 to float can take up to 9e-8 of it, so the float design loses a
//     pole to the unit circle below fc = 4.7e-5 and its slowest pole moves
//     by more than 0.1 % below 6.0e-5 (measured, every even order, either
//     end); at 1e-4 the designed margin, 3.9e-7, is 4.4 times what rounding
//     can take. Double's keeps the flush window below in a 32-bit size_t:
//     the slowest pole decays by about alpha_min pi d per frame, d the
//     cutoff's distance from 0 or 1/2, so the window grows as 1 / d (the
//     largest, order 16 at 1/2 - 1e-5: 232,710,082 frames).
//     Magnitude |H(fc)| = 1/sqrt(2) (-3.0103 dB), |H(0)| = 1, H(1/2) = 0;
//     monotone in between.
//   - A 1:1 sync_stage (chain.h): k_up = k_down = 1, channels() = Ch, any
//     Ch >= 1, frames interleaved; process() writes exactly n frames
//     (outputs_for(n) = n), noexcept, allocation-free; in and out the same
//     buffer or disjoint ones. Bit-identical for any chunking (the state is
//     per channel and per sample), across channel counts and across chain.h's
//     flush and process paths, at fft.h's fp-contraction policy (contraction
//     off: tests/CMakeLists.txt's comment on tap_dsp_chain_stages says why).
//     With contraction on, a compiler may fuse a call's expressions by the
//     path it takes (a single frame or a block, a pair of channels or one:
//     iir_kernels.h's biquad_cascade_block), so there the bits may depend on
//     the chunking and the channel count on a target with a fused
//     multiply-add, within the float-against-double bound below. Until 0.6.0
//     the stage ran the per-sample kernel on every path and its chunking
//     could not move a bit under any flags. The arithmetic is iir_kernels.h's
//     floating biquad cascade, transposed direct form II, one state pair per
//     section and channel, run a block at a time (biquad_cascade_block: the
//     per-sample kernel's expressions with the coefficients and state held
//     in locals for the call; 0.26 to 0.52 of the per-sample stage's
//     instructions at orders 2 to 8 on the family's embedded cores with an
//     FPU, bench/README.md).
//   - set_cutoff(fc): the constructor's design, in place. The sections, the
//     kernel's table, the flush window and the nominal latency are bit for
//     bit what a stage constructed at fc holds, and the deviation it reports
//     is the same number (the constructor runs the same function; bit for
//     bit at fft.h's fp-contraction policy, as above). noexcept and
//     allocation-free; its precondition is butterworth_supports_cutoff<S>(fc),
//     a TAP_EXPECTS (a debug-build assertion), never a throw. It leaves the
//     state alone: the recursion continues from where it was under the new
//     coefficients, a step the caller chooses (reset() first for a clean
//     start). It costs one design, nearly all of it the flush window:
//     measured on x86-64 (a shared 2.1 GHz Xeon VM, GCC 13.3 -O3, the
//     median of 401 calls, three runs), 1.1 to 2.8 us at order 4 and 5.5 to
//     14 us at order 8, float or double, at the cutoffs 0.45 / 2 ... 0.45 /
//     16 an oversampling kernel sets, and 77 us at order 16 in double at its
//     lowest cutoff. The design is double, so a core without a double FPU
//     pays libgcc for it: 2.3 M instructions a call at order 8 on a
//     Cortex-M33 or M4F over the cutoffs 0.45 / 2, 4, 8, 16, 2.8 M at most
//     (bench/README.md; 0.1 M on a Cortex-M55), more than a block of audio
//     can spare there. For occasional
//     changes (an oversampling factor from a setter), not for per-sample or
//     per-block modulation. Through 0.6.0 set_cutoff() and the constructor
//     also ran the deviation's group-delay scan, 80 to 93 % of their cost;
//     the deviation is now computed when it is read (below). Its stack is
//     the window's matrices (below): -fstack-usage gives the window's frame
//     as 1,432, 5,272 and 20,640 bytes at orders 4, 8 and 16 on the
//     Cortex-M33 (arm-none-eabi-gcc 13.2 -Os), the design around it 0.3 to
//     0.4 KiB more.
//   - Latency (the plan's I5, chain.h's reports_nominal_latency): a
//     recursive stage has no single delay, so it never provides
//     latency_output_frames(). It reports three numbers, in frames (x 1 /
//     rate for seconds):
//       group_delay_output_frames(f): -d(phase)/d(omega) at f cycles per
//         sample, 0 <= f < 1/2 (iir_design.h's sections_group_delay, double);
//       nominal_latency_output_frames(): the group delay at DC, the delay a
//         caller gives a parallel path (a dry signal, a sibling channel) to
//         line it up with this one;
//       latency_deviation_output_frames(): the passband group-delay
//         deviation, the largest |group_delay(f) - nominal| over 0 <= f <=
//         fc, on a 1,024-point grid. It is the alignment error a caller
//         accepts by compensating the nominal alone: a passband component
//         still arrives up to this many frames off the delayed parallel
//         path, so a sum of the two cancels where that offset is half the
//         component's period. A Butterworth's group delay rises toward fc,
//         so the deviation is its excess at the band edge (order 4: 2.51
//         frames on a nominal of 4.02 at fc = 0.1, 2.39 on 1.31 at 0.25).
//         Computed on every read, from the sections and the cutoff
//         (butterworth_latency_deviation, the same scan as a free function
//         for a caller that keeps the number): 1,025 group-delay evaluations
//         in double, 40 to 72 us at orders 4 and 8 on the x86-64 above and
//         about 28 M instructions at order 8 on a Cortex-M33. A design-time
//         query, not for the audio thread; const and without a cache, so two
//         threads may read it at once. A chain's
//         latency_deviation_output_frames() reads each stage's.
//     All three are the double design's.
//   - Flush (the plan's I6): window_frames() is a design-time constant, the
//     zero input frames after which the exact recursion of the coefficients
//     the kernel runs (the design rounded to S) leaves a state below the
//     smallest normal value of S from any state S can hold. With zero input
//     the cascade is x' = A x, A its 2(N/2) x 2(N/2) state-transition matrix;
//     the window is the smallest n with
//
//         ln(max_finite) + ln ||A^n||_inf  <  ln(min_normal),
//
//     plus a margin of one frame per section, found by binary lifting over
//     A^(2^j) kept normalized with a separate log scale. Each power is
//     recomputed from A by squaring rather than kept, so the stack holds a
//     few N x N matrices whatever the window (an array of the 63 powers a
//     64-bit size_t allows would take 32 KiB at order 8 and 126 KiB at
//     order 16, past the 64 KiB stack of the M33 legs): O(log^2 n) products
//     instead of O(log n), the same bits.
//     What the arithmetic leaves after it is the host's FTZ's
//     (iir_kernels.h's honest limit): usually exact silence; without FTZ a
//     subnormal residue; under FTZ, for some designs, a limit cycle just
//     above the smallest normal value (FTZ flushes a2 y while s1 survives,
//     and |a1| > 1 lifts it again), growing about as 1 / fc as the cutoff
//     falls (measured on x86-64, from full scale, in multiples of
//     std::numeric_limits<S>::min(): 124 at order 8, fc = 0.002, float, 243
//     at 0.001 and 2,351 at the float bound 1e-4; 99 at order 4, fc =
//     0.002, double, and 12,474 at double's bound 1e-5; none for most
//     designs between 0.01 and 0.3). Pinned bound at the battery's cutoffs,
//     0.005 to 0.45: 256 times the smallest normal value (3e-36 in float,
//     about -700 dB re full scale), orders 4 and 8, either FTZ setting
//     (`butterworth_test.TheTailAfterTheWindowIsBounded` under FTZ,
//     `butterworth_subnormals.TheTailAfterTheWindowIsBounded` without, in
//     the tap_dsp_iir_subnormals executable). In a chain, flush() therefore
//     equals zero padding up to that residue rather than bit for bit (chain.h
//     promises bit identity for stages that drain to an exactly zero
//     state), and what follows it is the residue, not exact silence.
//   - Non-finite input poisons the state until reset() (I9: no per-sample
//     check).
//   - Formats: double (golden) and float (float states, coefficients and
//     products, I14). Float against double, the largest difference over
//     48,000 frames of white noise at 0.9 of full scale (x86-64): order 4,
//     2.0e-7 at fc = 0.25, 4.5e-7 at 0.07, 9.5e-6 at 0.01, 1.1e-4 at 0.001
//     and 2.4e-3 at the float bound 1e-4; order 8, 2.6e-7, 8.9e-7, 1.2e-5,
//     8.3e-5 and 3.4e-3. Transposed DF II's float error grows as the poles
//     approach z = 1; below about 0.001 the float stage is stable but far
//     from the double one, and double is the profile for low cutoffs. The
//     battery bounds fc = 0.07 at 2e-6 (`TracksTheDoubleRecursion`).
//     Q15 / Q31: the kernel has them (direct form I, iir_kernels.h); this
//     stage does not offer them yet (a Butterworth at a low cutoff needs the
//     L-infinity section scaling iir_design.h does for the Chebyshev II).

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <span>
#include <stdexcept>
#include <type_traits>

#include "tap/dsp/detail/expects.h"
#include "tap/dsp/iir_design.h"
#include "tap/dsp/iir_kernels.h"

namespace tap::dsp {

    /// The lowest cutoff of a stage in format S, and the highest's distance
    /// from 1/2, in cycles per sample (the file header says why each).
    template <typename S>
    inline constexpr double k_butterworth_min_cutoff = std::is_same_v<S, float> ? 1e-4 : 1e-5;

    /// Whether a Butterworth stage in format S is designable at fc cycles per
    /// sample: k_butterworth_min_cutoff<S> <= fc <= 1/2 -
    /// k_butterworth_min_cutoff<S> (false for NaN). The design itself
    /// (design_butterworth_lowpass, double) asks for S = double's range.
    template <typename S>
    constexpr bool butterworth_supports_cutoff(double fc) noexcept {
        return fc >= k_butterworth_min_cutoff<S> && fc <= 0.5 - k_butterworth_min_cutoff<S>;
    }

    namespace detail {

        /// A cascade's state-transition matrix or a power of it, kept as
        /// M * e^log_scale with M normalized to ||M||_inf = 1. Block lower
        /// triangular in 2 x 2 blocks, one per section (a section's next state
        /// depends on its own state and the sections' before it, never a
        /// later one's), and so is every product of two such: row r is zero
        /// past column row_end(r). The loops stop there. The terms they leave
        /// out are exact zeros added to a sum that is never a negative zero,
        /// so every element has the bits the full loops give.
        template <std::size_t D>
        struct scaled_matrix {
            std::array<double, D * D> m{};
            double                    log_scale = 0.0;

            /// One past the last column of row r that can be nonzero.
            static constexpr std::size_t row_end(std::size_t r) noexcept { return 2 * (r / 2) + 2; }

            void normalize() noexcept {
                double n = 0.0;
                for (std::size_t r = 0; r < D; ++r) {
                    double row = 0.0;
                    for (std::size_t c = 0; c < row_end(r); ++c) {
                        row += std::fabs(m[r * D + c]);
                    }
                    n = std::max(n, row);
                }
                if (n > 0.0) {
                    for (std::size_t r = 0; r < D; ++r) {
                        for (std::size_t c = 0; c < row_end(r); ++c) {
                            m[r * D + c] /= n;
                        }
                    }
                    log_scale += std::log(n);
                }
                else {
                    log_scale = -std::numeric_limits<double>::infinity();
                }
            }

            friend scaled_matrix operator*(const scaled_matrix& a, const scaled_matrix& b) noexcept {
                scaled_matrix r;
                for (std::size_t i = 0; i < D; ++i) {
                    for (std::size_t k = 0; k < row_end(i); ++k) {
                        const double aik = a.m[i * D + k];
                        for (std::size_t j = 0; j < row_end(k); ++j) {
                            r.m[i * D + j] += aik * b.m[k * D + j];
                        }
                    }
                }
                r.log_scale = a.log_scale + b.log_scale;
                r.normalize();
                return r;
            }
        };

        /// The zero-input state-transition matrix of a transposed direct form
        /// II cascade of K sections, state (s1, s2) per section in processing
        /// order.
        template <std::size_t K>
        scaled_matrix<2 * K> cascade_transition(const std::array<iir_section, K>& sec) noexcept {
            constexpr std::size_t D = 2 * K;
            scaled_matrix<D>      a;
            std::array<double, D> x{}; // section 0's input: the zero input, as a row over the old state
            for (std::size_t s = 0; s < K; ++s) {
                const iir_section&    c = sec[s];
                std::array<double, D> y = x; // y = b0 x + s1 (old)
                for (double& v : y) {
                    v *= c.b0;
                }
                y[2 * s] += 1.0;
                for (std::size_t j = 0; j < D; ++j) { // s1' = b1 x - a1 y + s2; s2' = b2 x - a2 y
                    a.m[(2 * s) * D + j]     = c.b1 * x[j] - c.a1 * y[j];
                    a.m[(2 * s + 1) * D + j] = c.b2 * x[j] - c.a2 * y[j];
                }
                a.m[(2 * s) * D + 2 * s + 1] += 1.0;
                x = y; // the next section's input is this one's output
            }
            a.normalize();
            return a;
        }

        /// The flush window (the file header's formula) of the cascade the
        /// sections describe: zero input frames after which its exact
        /// recursion's state is below S's smallest normal value from any
        /// state S can hold, plus one frame per section. noexcept and
        /// allocation-free: binary lifting over A^(2^j). The descent needs
        /// every power below the top one, highest first; k_kept of them, evenly
        /// spaced, are kept from the ascent and the others are squared up from
        /// the nearest kept one below (the same products in the same order as
        /// squaring up from A, so the same bits, as keeping all of them would
        /// give), so the stack holds k_kept matrices and a few more whatever
        /// the window. A cascade whose window does not fit half of
        /// size_t's range, or that does not decay (none that
        /// butterworth_supports_cutoff admits), returns size_t's largest value.
        template <typename S, std::size_t K>
        std::size_t cascade_flush_window(const std::array<iir_section, K>& sections) noexcept {
            using mat                         = scaled_matrix<2 * K>;
            constexpr std::size_t k_levels    = std::numeric_limits<std::size_t>::digits - 1;
            constexpr std::size_t k_saturated = std::numeric_limits<std::size_t>::max();
            const mat             a           = cascade_transition(sections);
            const double          start       = std::log(static_cast<double>(std::numeric_limits<S>::max()));
            const double          goal        = std::log(static_cast<double>(std::numeric_limits<S>::min()));
            const auto            done        = [&](const mat& p) noexcept { return start + p.log_scale < goal; };

            // The powers A^(2^j), j < levels, up to the first that alone is done.
            std::size_t levels = 1;
            for (mat top = a; !done(top); ++levels) {
                if (levels == k_levels) {
                    return k_saturated;
                }
                top = top * top;
            }
            // A^(2^(i * stride)), i < k_kept: the checkpoints the descent squares up from.
            constexpr std::size_t   k_kept = 4;
            const std::size_t       stride = (levels + k_kept - 1) / k_kept;
            std::array<mat, k_kept> kept;
            kept[0] = a;
            for (std::size_t i = 1; i < k_kept && i * stride < levels; ++i) {
                kept[i] = kept[i - 1];
                for (std::size_t q = 0; q < stride; ++q) {
                    kept[i] = kept[i] * kept[i];
                }
            }
            const auto power = [&](std::size_t j) noexcept { // A^(2^j), j < levels
                mat p = kept[j / stride];
                for (std::size_t i = j % stride; i > 0; --i) {
                    p = p * p;
                }
                return p;
            };
            // Binary lifting: the largest n from those powers not yet done.
            // The top power is done alone, so the descent starts below it.
            mat         acc;
            bool        have = false;
            std::size_t n    = 0;
            for (std::size_t j = levels - 1; j-- > 0;) {
                const mat pj   = power(j);
                const mat cand = have ? acc * pj : pj;
                if (!done(cand)) {
                    acc  = cand;
                    have = true;
                    n += std::size_t{1} << j;
                }
            }
            mat p = have ? acc * a : a;
            ++n;
            while (!done(p)) {
                p = p * a;
                ++n;
            }
            return n + K;
        }

    } // namespace detail

    /// The N/2 sections of an order-N Butterworth low-pass at fc cycles per
    /// sample, in increasing pole radius. Precondition:
    /// butterworth_supports_cutoff<double>(fc).
    template <std::size_t N>
    std::array<iir_section, N / 2> design_butterworth_lowpass(double fc) noexcept {
        static_assert(N >= 2 && N <= 16 && N % 2 == 0, "design_butterworth_lowpass: even orders 2 ... 16");
        TAP_EXPECTS(butterworth_supports_cutoff<double>(fc));
        constexpr double               pi = std::numbers::pi;
        const double                   k  = std::tan(pi * fc);
        const double                   k2 = k * k;
        std::array<iir_section, N / 2> sec{};
        // k = N/2 - 1 has the largest alpha (lowest Q, smallest radius): run it first.
        for (std::size_t i = 0; i < N / 2; ++i) {
            const std::size_t idx = N / 2 - 1 - i;
            const double alpha = 2.0 * std::sin(static_cast<double>(2 * idx + 1) * pi / (2.0 * static_cast<double>(N)));
            const double a0    = 1.0 + alpha * k + k2;
            sec[i].b0          = k2 / a0;
            sec[i].b1          = 2.0 * k2 / a0;
            sec[i].b2          = k2 / a0;
            sec[i].a1          = 2.0 * (k2 - 1.0) / a0;
            sec[i].a2          = (1.0 - alpha * k + k2) / a0;
        }
        return sec;
    }

    /// The passband group-delay deviation of a low-pass cascade, in frames:
    /// max over 0 <= f <= cutoff of |group delay(f) - group delay(0)|, on a
    /// 1,024-point grid (1,025 evaluations of iir_design.h's
    /// sections_group_delay, double). What
    /// basic_butterworth_lowpass::latency_deviation_output_frames() returns,
    /// as a function of the design alone, for a caller that wants the number
    /// once and keeps it. A design-time computation, not for the audio
    /// thread: the scan is double, 28 M instructions at order 8 on a
    /// Cortex-M33 (soft double), 0.7 M on a Cortex-M55 (bench/README.md).
    template <std::size_t K>
    double butterworth_latency_deviation(const std::array<iir_section, K>& sections, double cutoff) noexcept {
        const double  nominal = sections_group_delay(sections, 0.0);
        double        dev     = 0.0;
        constexpr int k_grid  = 1024;
        for (int i = 0; i <= k_grid; ++i) {
            const double f = cutoff * static_cast<double>(i) / k_grid;
            dev            = std::max(dev, std::fabs(sections_group_delay(sections, f) - nominal));
        }
        return dev;
    }

    /// An order-N Butterworth low-pass as a 1:1 chain stage over Ch
    /// interleaved channels, float or double (the file header's contract).
    template <typename S, std::size_t N, std::size_t Ch = 1>
    class basic_butterworth_lowpass {
        static_assert(N >= 2 && N <= 16 && N % 2 == 0, "basic_butterworth_lowpass: even orders 2 ... 16");
        static_assert(N / 2 <= k_max_cheby2_sections, "basic_butterworth_lowpass: the kernel's section capacity");
        static_assert(Ch >= 1, "basic_butterworth_lowpass: at least one channel");
        static_assert(std::is_same_v<S, float> || std::is_same_v<S, double>,
                      "basic_butterworth_lowpass: float and double (the file header: Q15 / Q31 not offered yet)");

      public:
        using sample                            = S;
        static constexpr std::size_t k_up       = 1;
        static constexpr std::size_t k_down     = 1;
        static constexpr std::size_t k_order    = N;
        static constexpr std::size_t k_sections = N / 2;
        static constexpr std::size_t channels() noexcept { return Ch; }

        /// The construction-time check, as a value: empty when cutoff is
        /// designable (butterworth_supports_cutoff<S>).
        static constexpr const char* validate(double cutoff) noexcept {
            if (butterworth_supports_cutoff<S>(cutoff)) {
                return nullptr;
            }
            return std::is_same_v<S, float>
                       ? "basic_butterworth_lowpass: a float stage's cutoff must be in [1e-4, 1/2 - 1e-4] cycles per "
                         "sample"
                       : "basic_butterworth_lowpass: a double stage's cutoff must be in [1e-5, 1/2 - 1e-5] cycles per "
                         "sample";
        }

        /// cutoff: the -3 dB frequency in cycles per sample
        /// (butterworth_supports_cutoff<S>); the state starts at zero.
        explicit basic_butterworth_lowpass(double cutoff) {
            if (const char* why = validate(cutoff)) {
                raise<std::invalid_argument>(why);
            }
            design(cutoff);
        }

        /// Redesigns at cutoff, in place: bit for bit what the constructor
        /// computes there (the file header's set_cutoff()). Leaves the state
        /// alone; reset() for a clean start. noexcept, allocation-free, one
        /// design's cost: for occasional changes, not modulation.
        /// Precondition: butterworth_supports_cutoff<S>(cutoff).
        void set_cutoff(double cutoff) noexcept {
            TAP_EXPECTS(butterworth_supports_cutoff<S>(cutoff));
            design(cutoff);
        }

        /// n interleaved frames in, n out; in and out may be the same buffer.
        /// Precondition (debug builds, where a bit governs S): the calling
        /// thread flushes subnormals (subnormals.h).
        std::size_t process(const S* in, std::size_t n, S* out) noexcept {
            biquad_cascade_block<S, N / 2, Ch>(m_table, m_state.data(), in, n, out);
            return n;
        }

        std::size_t outputs_for(std::size_t n) const noexcept { return n; }
        void        reset() noexcept { m_state.fill(S(0)); }

        /// Zero input frames after which the exact recursion's state is below
        /// the smallest normal value of S (the residue after it: the file
        /// header's Flush).
        std::size_t window_frames() const noexcept { return m_window; }

        /// The group delay at f cycles per sample, 0 <= f < 1/2, in frames:
        /// -d(phase)/d(omega) of the double design.
        double group_delay_output_frames(double f) const noexcept { return sections_group_delay(m_sections, f); }

        /// The nominal latency: the group delay at DC, in frames. Delay a
        /// parallel path by this to align it with the stage's output.
        double nominal_latency_output_frames() const noexcept { return m_nominal; }

        /// The passband group-delay deviation, in frames: max over 0 <= f <=
        /// cutoff of |group_delay_output_frames(f) - nominal|, on a
        /// 1,024-point grid. The alignment error left after compensating the
        /// nominal alone (chain.h's reports_nominal_latency); 0 would mean the
        /// passband is delayed alike. A Butterworth's peaks at the band edge.
        /// Computed on each call (butterworth_latency_deviation: a design-time
        /// cost, the file header's numbers), not stored.
        double latency_deviation_output_frames() const noexcept {
            return butterworth_latency_deviation(m_sections, m_cutoff);
        }

        double                                cutoff() const noexcept { return m_cutoff; }
        const std::array<iir_section, N / 2>& sections() const noexcept { return m_sections; }

        /// True when every channel's state is exactly zero (silence from here on).
        bool is_silent() const noexcept {
            for (const S v : m_state) {
                if (v != S(0)) {
                    return false;
                }
            }
            return true;
        }

      private:
        /// Everything that depends on the cutoff: the constructor's and
        /// set_cutoff()'s one path. Leaves m_state alone.
        void design(double cutoff) noexcept {
            m_cutoff   = cutoff;
            m_sections = design_butterworth_lowpass<N>(cutoff);
            m_table    = biquad_cascade_table<S>(std::span<const iir_section>(m_sections));
            // The window of the recursion the kernel runs: the design rounded to S.
            const auto in_s = [](double v) noexcept { return static_cast<double>(static_cast<S>(v)); };
            std::array<iir_section, N / 2> run{};
            for (std::size_t i = 0; i < N / 2; ++i) {
                const iir_section& c = m_sections[i];
                run[i]               = {in_s(c.b0), in_s(c.b1), in_s(c.b2), in_s(c.a1), in_s(c.a2)};
            }
            m_window  = detail::cascade_flush_window<S>(run);
            m_nominal = sections_group_delay(m_sections, 0.0);
        }

        double                           m_cutoff = 0.0;
        std::array<iir_section, N / 2>   m_sections{};
        biquad_cascade_table<S>          m_table{};
        std::array<iir_state<S>, N * Ch> m_state{}; ///< per channel, z1 then z2 per section (iir_kernels.h)
        std::size_t                      m_window  = 0;
        double                           m_nominal = 0.0;
    };

    template <std::size_t Ch = 1>
    using butterworth4 = basic_butterworth_lowpass<double, 4, Ch>; ///< order 4, the golden model
    template <std::size_t Ch = 1>
    using butterworth4_32 = basic_butterworth_lowpass<float, 4, Ch>; ///< order 4, the embedded profile
    template <std::size_t Ch = 1>
    using butterworth8 = basic_butterworth_lowpass<double, 8, Ch>; ///< order 8, the golden model
    template <std::size_t Ch = 1>
    using butterworth8_32 = basic_butterworth_lowpass<float, 8, Ch>; ///< order 8, the embedded profile

} // namespace tap::dsp
