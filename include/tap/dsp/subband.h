/// @file subband.h
/// @brief The sub-band round trip: a low-rate change, interpolated back onto the input delayed by its exact latency.
// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// A process that needs fine frequency resolution in a low band (a 4096-point
// FFT whose bins are D times narrower after decimating by D) runs at the low
// rate, and only what it CHANGES comes back: SampleRateTap's
// docs/X4_X8_PLAN.md, sections 6.3 and 6.4 (X9, X10, X12, X15), landed here
// as the substrate rule asks, generic over chain.h's sync_stage and proven on
// test-local stages built on nyquist.h, as chain.h was proven on decimate.h.
//
//   out[n] = x[n - T] + up(change(down(x)))[n],   T = D tau_down + D L_op + tau_up
//
//   D         the factor: Down converts by 1/D and Up by D/1, inverse integer
//             factors (a static_assert);
//   tau_down  Down's exact delay in its own output frames (the low rate), so
//             D tau_down is in full-rate samples;
//   L_op      the op's latency in whole low-rate frames (subband_op);
//   tau_up    Up's exact delay in its own output frames (the full rate).
//
// Only the change travels down and back. Everything the op leaves alone,
// inside the low band and above it, reaches the output through the delay
// line untouched, so the decimator's aliasing in the bins the op does not
// touch is never re-injected, no guard band filter is needed, and bands
// compose by addition with no crossover filter (Crochiere and Rabiner 1983,
// sub-band structures; Vaidyanathan 1993, the delay-complementary pair
// z^-T - H(z) of a linear-phase H; Allen 1977 for the overlap-add
// modification an STFT op makes at the low rate). Published literature only.
//
// Contract, as numbers:
//   - Factor and formats: Down's ratio is 1/D and Up's D/1 for one integer
//     D >= 1, both as reduced compile-time ratios; the sample type and the
//     channel count Ch come from the stages, and Down, Up and Op must agree
//     on both (static_asserts). float and double only: fixed point is out of
//     scope (a static_assert), because the sum and the delay line need a
//     Q-format design of their own (headroom for x + change, where the
//     single rounding sits), and per-primitive fixed-point adoption is a
//     documented design each time.
//   - Where T comes from (exact_delay_stage): each stage reports an exact
//     latency, latency_output_frames() (an exact_ratio in its output frames)
//     or latency_input_samples() (an integer at its input rate), as chain.h
//     reads them, or an exact delay path, delay_path_output_frames()
//     (reports_delay_path: the member a recursive stage with a near-linear
//     phase provides in place of a latency; SampleRateTap's docs/IIR_PLAN.md,
//     I5). A tap::dsp::chain is walked stage by stage, so a chain qualifies
//     when every stage does, whichever form each reports. A stage with none,
//     such as a recursive stage with a nonlinear phase, is not an
//     exact_delay_stage and the scaffold does not compile with it.
//   - T is an integer number of full-rate frames, exactly, or the scaffold is
//     refused: D tau_down + tau_up must be whole (D L_op always is). When both
//     stages' delays are compile-time constants (static constexpr members,
//     or a chain of such stages), the types decide it and a non-whole sum is
//     a static_assert; when either is decided at run time (a stage's tap
//     count is construction geometry) the constructor rejects it through
//     tap::dsp::raise<std::invalid_argument>, and supports(down, up) is the
//     predicate to check first where the build has no exceptions. A chain of
//     integer-factor L-th-band FIR stages is always whole: a stage of N =
//     2mL - 1 taps delays (N - 1) / 2 = mL - 1 samples at its higher rate.
//     A chain with a mixed-ratio stage may not be (a 3/8 stage puts a 3 in
//     the denominator).
//   - The change path leads: from the reset state, after j input frames the
//     change path has written at least j full-rate frames (Down writes its
//     output k by input k D, decimate.h's and nyquist.h's alignment, and Up
//     writes D frames per input), checked for j = 1 .. 2D at construction
//     (stages are synchronous, so the count repeats with period D) and
//     rejected through raise<std::invalid_argument> otherwise. The surplus,
//     at most lead() = max_j (written - j) frames (D - 1 for those stages),
//     waits in a small queue for the input it belongs to.
//   - Identity, by construction: an op whose change is exactly zero gives
//     out[n] == x[n - T] bit for bit, for every bit pattern of x (negative
//     zero and subnormals included), whatever T and whatever the stages. A
//     zero-primed FIR turns zeros into exact zeros, and where the change is
//     exactly zero the delayed input is written as it is: a select, not an
//     addition, because IEEE 754 makes -0 + +0 = +0. It is a structural
//     guarantee and does not test T; only a cancellation does.
//   - Latency: latency() = max(T, aligned_latency) full-rate frames, an
//     integer; round_trip_latency() = T. latency_output_frames() reports it as
//     an exact_ratio, so the scaffold is itself an exact_delay_stage and a
//     sync_stage of ratio 1/1 (k_up = k_down = 1), and composes in a chain.
//   - Composition (bands compose by addition): out = x[n - A] + sum over bands
//     of each band's change, every band aligned to the common latency A, the
//     largest T. A band constructed with aligned_latency = A holds its extra
//     delay A - T as floor((A - T) / D) frames at its low rate, before Up, and
//     the (A - T) mod D left over at the full rate, after it. One band runs
//     subband_output::input_plus_change and holds the delay line; every other
//     band runs subband_output::change_only and holds none, so the input
//     path is counted once; the caller adds the change-only bands' outputs.
//     Disjoint ownership of the bins, and complementary tapers where the
//     bands meet, are the ops' business (the X15 precondition below).
//   - Memory (X12), samples of the stage type, fixed at construction:
//     the delay line A Ch (input_plus_change only); the extra low-rate delay
//     floor((A - T) / D) Ch; the change queue (lead() + (A - T) mod D +
//     B_up(B_down(k_block))) Ch, with B the frames a stage writes for n
//     inputs, chain.h's scratch bound n L / M + 2 per stage; two low-rate
//     blocks of B_down(k_block) Ch each, and one k_block Ch block of zeros
//     for flush; plus whatever the stages and the op hold. X12's figures,
//     for an op of one 4096-point frame's latency on SampleRateTap's economy
//     chains: T = 110 + 4 x 4096 = 16494 at D = 4 (129 KiB of double delay
//     line per channel), 206 + 8 x 4096 = 32974 at D = 8 (258 KiB) and
//     378 + 16 x 4096 = 65914 at D = 16 (515 KiB); a D = 4 band beside a
//     D = 8 band holds its 16480 extra frames as 4120 at 12 kHz (32 KiB in
//     double).
//   - process(in, n, out) reads n interleaved frames and writes exactly n,
//     for any n, through internal blocks of k_block = 256 frames: noexcept,
//     allocation-free, and bit-identical for any chunking of the same stream
//     (every part is). in == out is allowed (in place); a partial overlap is
//     not.
//   - flush(out): end of stream. Feeds window_frames() zero frames through
//     process() and writes as many, so it equals zero padding by
//     construction. window_frames() is a bound, valid from any state, on the
//     zero input frames after which the output is silence:
//       max(A, W_down + F_down(W_op + E + W_up) + lead() + (A - T) mod D)
//     (the first term input_plus_change only), with W a stage's
//     window_frames() (a chain's composed from its stages, each stage's
//     window plus the inputs the stages before it need to write it), E the
//     extra low-rate frames, and F_down(k) = ceil(k M / L) + M per stage,
//     composed backward, the inputs a stage on its ratio's schedule needs
//     for k outputs from any phase. Down, Up and Op must each report
//     window_frames() for flush to compile (a static_assert names the
//     rule). The scaffold is left in the zero-fed state; reset() before
//     reuse.
//   - reset(): every stage, the op, both delay lines and the change queue
//     back to the zero-primed state (the constructor resets them too:
//     stages arrive in whatever state their previous owner left them in).
//   - Construction allocates and may reject (raise); process(), flush() and
//     reset() are noexcept and allocation-free. The stages and the op are
//     driven only through the scaffold (read-only accessors), as chain.h's
//     scratch bound requires.
//
// Honest limits:
//   - Cancellation needs linear phase and an exact delay. A round trip whose
//     phase is off by phi at a frequency leaves 2 sin(phi / 2) of a removed
//     component there: -40 dB needs phi <= 0.57 degrees and -60 dB phi <=
//     0.057 degrees, summed over every stage down and up. The scaffold
//     cannot check phase; a stage's exact delay is a promise that every
//     frequency of the usable band is delayed by it (linear-phase FIR stages
//     keep it by construction).
//   - Removing a low-band component (change = -down(x)) leaves the round
//     trip's passband ripple, |1 - H_down H_up| at the component, not
//     silence. With the battery's test stages (nyquist.h at 70 dB, passband
//     3/8 of the low rate) it is -61.3 to -116 dB at the tones the battery
//     removes, and at worst -65.0 dB (D = 4, one stage each way) and
//     -61.3 dB (D = 8, two each way) over the whole usable band, by a dense
//     sweep of the designs' responses, which the streamed removals match to
//     0.01 dB (test_subband.cpp, RemovingTheBandReachesTheStagesRippleFloor;
//     float and double alike). With the delay one sample off either way the
//     same removal leaves 2 sin(pi f): -5.1 to -24.0 dB at those tones.
//   - The change must be confined to the bins up to p r_low (X15, a
//     precondition on Op, not enforced): p is the stages' passband edge as a
//     fraction of the low rate r_low (3/8 at the family's economy profile:
//     75 % of the low band). Above it, from p r_low to (1 - p) r_low, the
//     low band carries the decimator's transition-band aliasing; an op that
//     changes those bins (an unmasked time-domain op, say) puts the aliasing
//     into the change and Up writes it back into the output: a tone at
//     0.55 r_low comes back with an alias at 0.45 r_low, where the input has
//     nothing, at -21.4 dB re the tone with the test stages at D = 4
//     (UnmaskedChangeIsReinjected), where a confined change would leave at
//     most the stopband, 70 dB down. Masking and tapering belong to the op
//     (an STFT's change mode).
//   - A non-whole D tau_down + tau_up is refused, so a chain with a factor of
//     3 in its latency's denominator cannot enter the scaffold.
//   - Resolution is bought with latency and memory: T grows by D L_op, and
//     the delay line by as much (the memory figures above).

#pragma once

#include <algorithm>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include "tap/dsp/chain.h"
#include "tap/dsp/detail/expects.h"
#include "tap/dsp/sample_traits.h"

namespace tap::dsp {

    // ANCHOR: subband_delay_path
    /// A stage's exact delay path: the member a recursive stage with a
    /// near-linear phase reports in place of an exact latency (it delays every
    /// passband frequency by about this many of its output frames, a whole
    /// number of samples at its higher rate, with a stated deviation). The
    /// sub-band round trip takes T from it; a chain.h latency it is not, so a
    /// chain's latency_output_frames() never reads it.
    template <typename T>
    concept reports_delay_path = requires(const T& s) {
        { s.delay_path_output_frames() } -> std::same_as<exact_ratio>;
    };
    // ANCHOR_END: subband_delay_path

    namespace detail {

        /// Whether T is a tap::dsp::chain, and its stage types when it is.
        template <typename T>
        struct chain_parts {
            static constexpr bool k_is_chain = false;
        };
        template <sync_stage First, sync_stage... Rest>
        struct chain_parts<chain<First, Rest...>> {
            static constexpr bool k_is_chain = true;
            using stages                     = std::tuple<First, Rest...>;
        };
        template <typename T>
        inline constexpr bool k_is_chain = chain_parts<T>::k_is_chain;

        /// A stage that reports an exact delay in one of the three forms; a
        /// chain when every stage does.
        template <typename T>
        struct has_exact_delay {
            static constexpr bool k_value =
                reports_latency_output<T> || reports_latency_input<T> || reports_delay_path<T>;
        };
        template <sync_stage First, sync_stage... Rest>
        struct has_exact_delay<chain<First, Rest...>> {
            static constexpr bool k_value = has_exact_delay<First>::k_value && (has_exact_delay<Rest>::k_value && ...);
        };

    } // namespace detail

    // ANCHOR: subband_exact_delay_stage
    /// A sync_stage whose delay is exact: latency_output_frames(),
    /// latency_input_samples() or delay_path_output_frames(), or a chain whose
    /// every stage reports one of them.
    template <typename T>
    concept exact_delay_stage = sync_stage<T> && detail::has_exact_delay<T>::k_value;
    // ANCHOR_END: subband_exact_delay_stage

    // ANCHOR: subband_op
    /// The low-rate process of the sub-band round trip. Contract, as numbers:
    ///   - sample and channels() as a sync_stage declares them (a floating
    ///     sample type; a static constexpr count of at least 1), equal to the
    ///     stages';
    ///   - latency_frames(): L_op, the op's latency in whole low-rate frames,
    ///     constant for the op's life;
    ///   - process(in, n, change): reads n interleaved low-rate frames and
    ///     writes n frames of change, change[k] = result[k] - in[k - L_op] (in
    ///     before the stream is zero): what the op does to the band, not the
    ///     band. n may be 0. Bit-identical for any chunking of the stream;
    ///     in and change never alias; noexcept and allocation-free. An op that
    ///     changes nothing writes exact zeros;
    ///   - reset(): back to the zero-primed state;
    ///   - optional, for flush: window_frames(), the zero input frames after
    ///     which the change is silence (L_op for an op that holds L_op frames).
    /// Precondition, not enforced (X15): the change is confined to the bins up
    /// to p r_low, the stages' passband (the file header's honest limits).
    template <typename T>
    concept subband_op =
        requires(T op, const T& cop, const typename T::sample* in, typename T::sample* change, std::size_t n) {
            typename T::sample;
            requires sample_type<typename T::sample>;
            typename std::integral_constant<std::size_t, T::channels()>;
            requires T::channels() >= 1;
            { cop.latency_frames() } noexcept -> std::same_as<std::size_t>;
            { op.process(in, n, change) } noexcept;
            { op.reset() } noexcept;
        };
    // ANCHOR_END: subband_op

    /// What a sub-band round trip writes (the file header's composition).
    enum class subband_output : std::uint8_t {
        input_plus_change, ///< x[n - A] + change: holds the delay line, A frames per channel
        change_only,       ///< the change alone, aligned to A: one band of a composition, no delay line
    };

    namespace detail {

        /// The product of the ratios of the stages after I in Stages (a
        /// std::tuple of stage types): converts stage I's frames to the last's.
        template <typename Stages, std::size_t I>
        constexpr exact_ratio ratio_after_stage() noexcept {
            return []<std::size_t... J>(std::index_sequence<J...>) {
                exact_ratio r{1, 1};
                ((r = J > I ? r
                                  * exact_ratio{std::tuple_element_t<J, Stages>::k_up,
                                                std::tuple_element_t<J, Stages>::k_down}
                            : r),
                 ...);
                return r;
            }(std::make_index_sequence<std::tuple_size_v<Stages>>{});
        }

        /// A stage's exact delay in its own output frames, from whichever
        /// form it reports; a chain's composed over its stages.
        template <typename T>
        exact_ratio exact_delay_output_frames(const T& s) noexcept {
            if constexpr (k_is_chain<T>) {
                using stages      = typename chain_parts<T>::stages;
                exact_ratio total = {0, 1};
                [&]<std::size_t... I>(std::index_sequence<I...>) {
                    ((total =
                          total + exact_delay_output_frames(s.template stage<I>()) * ratio_after_stage<stages, I>()),
                     ...);
                }(std::make_index_sequence<std::tuple_size_v<stages>>{});
                return total;
            }
            else if constexpr (reports_latency_output<T> || reports_latency_input<T>) {
                return stage_latency_output_frames(s);
            }
            else {
                static_assert(reports_delay_path<T>,
                              "the sub-band round trip needs every stage to report an exact latency or an exact "
                              "delay path");
                return s.delay_path_output_frames();
            }
        }

        // The same delay when the type decides it: a static constexpr member
        // in the form the stage reports (the form exact_delay_output_frames
        // reads), or a chain of such stages. Each concept first asks that the
        // member be static (its address a plain function pointer), and only
        // then calls it without an object in a constant expression: the
        // conjunction never forms that call for a non-static member.
        template <typename F>
        concept static_member_function = !std::is_member_function_pointer_v<F>;
        template <typename T>
        concept static_latency_output = static_member_function<decltype(&T::latency_output_frames)> && requires {
            typename std::integral_constant<std::uint64_t, T::latency_output_frames().den>;
        };
        template <typename T>
        concept static_latency_input = static_member_function<decltype(&T::latency_input_samples)> && requires {
            typename std::integral_constant<std::uint64_t, static_cast<std::uint64_t>(T::latency_input_samples())>;
        };
        template <typename T>
        concept static_delay_path = static_member_function<decltype(&T::delay_path_output_frames)> && requires {
            typename std::integral_constant<std::uint64_t, T::delay_path_output_frames().den>;
        };

        template <typename T>
        struct static_delay {
            static constexpr bool k_known = reports_latency_output<T>  ? static_latency_output<T>
                                            : reports_latency_input<T> ? static_latency_input<T>
                                                                       : static_delay_path<T>;
        };
        template <sync_stage First, sync_stage... Rest>
        struct static_delay<chain<First, Rest...>> {
            static constexpr bool k_known = static_delay<First>::k_known && (static_delay<Rest>::k_known && ...);
        };

        /// exact_delay_output_frames for a type whose delay is a constant.
        template <typename T>
        constexpr exact_ratio static_delay_output_frames() noexcept {
            if constexpr (k_is_chain<T>) {
                using stages = typename chain_parts<T>::stages;
                return []<std::size_t... I>(std::index_sequence<I...>) {
                    exact_ratio total = {0, 1};
                    ((total = total
                              + static_delay_output_frames<std::tuple_element_t<I, stages>>()
                                    * ratio_after_stage<stages, I>()),
                     ...);
                    return total;
                }(std::make_index_sequence<std::tuple_size_v<stages>>{});
            }
            else if constexpr (reports_latency_output<T>) {
                return T::latency_output_frames();
            }
            else if constexpr (reports_latency_input<T>) {
                return exact_ratio{static_cast<std::uint64_t>(T::latency_input_samples()), 1}
                       * exact_ratio{static_cast<std::uint64_t>(T::k_up), static_cast<std::uint64_t>(T::k_down)};
            }
            else {
                return T::delay_path_output_frames();
            }
        }

        /// D tau_down + tau_up is whole, where the types decide it; true
        /// (decided at construction) where they do not.
        template <typename Down, typename Up, std::uint64_t D>
        constexpr bool stages_delay_whole_if_static() noexcept {
            if constexpr (static_delay<Down>::k_known && static_delay<Up>::k_known) {
                return (exact_ratio{D, 1} * static_delay_output_frames<Down>() + static_delay_output_frames<Up>()).den
                       == 1;
            }
            else {
                return true;
            }
        }

        /// Frames a stage writes for n input frames, at most, from any state
        /// its process() leaves it in: chain.h's scratch bound, composed
        /// through a chain's stages.
        template <typename T>
        constexpr std::size_t max_frames_out(std::size_t n) noexcept {
            if constexpr (k_is_chain<T>) {
                using stages = typename chain_parts<T>::stages;
                return [n]<std::size_t... I>(std::index_sequence<I...>) {
                    std::size_t c = n;
                    ((c = max_outputs_for<std::tuple_element_t<I, stages>>(c)), ...);
                    return c;
                }(std::make_index_sequence<std::tuple_size_v<stages>>{});
            }
            else {
                return max_outputs_for<T>(n);
            }
        }

        /// Input frames a stage needs, at most, to write k frames from any
        /// phase: ceil(k M / L) + M for a stage on its ratio's schedule,
        /// composed backward through a chain's stages.
        template <typename T>
        constexpr std::size_t max_frames_in(std::size_t k) noexcept {
            if constexpr (k_is_chain<T>) {
                using stages = typename chain_parts<T>::stages;
                return [k]<std::size_t... I>(std::index_sequence<I...>) {
                    std::size_t c = k;
                    ((c = max_frames_in<std::tuple_element_t<std::tuple_size_v<stages> - 1 - I, stages>>(c)), ...);
                    return c;
                }(std::make_index_sequence<std::tuple_size_v<stages>>{});
            }
            else {
                if (k == 0) {
                    return 0;
                }
                const std::size_t l = T::k_up;
                const std::size_t m = T::k_down;
                return (k * m + l - 1) / l + m;
            }
        }

        /// Whether a stage reports the zero input frames that drain it: its
        /// window_frames(), or every stage's for a chain.
        template <typename T>
        struct drains {
            static constexpr bool k_value = reports_window<T>;
        };
        template <sync_stage First, sync_stage... Rest>
        struct drains<chain<First, Rest...>> {
            static constexpr bool k_value = drains<First>::k_value && (drains<Rest>::k_value && ...);
        };

        /// Zero input frames after which a stage's output is silence, from any
        /// state: its window_frames(); for a chain, the first stage's window
        /// plus the inputs it needs to write the rest of the chain's window,
        /// composed backward.
        template <typename T>
        std::size_t drain_frames(const T& s) noexcept {
            if constexpr (k_is_chain<T>) {
                using stages  = typename chain_parts<T>::stages;
                std::size_t w = 0;
                [&]<std::size_t... I>(std::index_sequence<I...>) {
                    ((w = drain_frames(s.template stage<std::tuple_size_v<stages> - 1 - I>())
                          + max_frames_in<std::tuple_element_t<std::tuple_size_v<stages> - 1 - I, stages>>(w)),
                     ...);
                }(std::make_index_sequence<std::tuple_size_v<stages>>{});
                return w;
            }
            else {
                return static_cast<std::size_t>(s.window_frames());
            }
        }

    } // namespace detail

    /// The sub-band round trip over Down (1/D), Up (D/1) and the low-rate Op;
    /// see the file header.
    template <exact_delay_stage Down, exact_delay_stage Up, subband_op Op>
    class subband_delta {
      public:
        using sample = typename Down::sample;
        static_assert(std::is_same_v<sample, typename Up::sample> && std::is_same_v<sample, typename Op::sample>,
                      "the sub-band round trip's stages and op must use the same sample type");
        static_assert(std::is_floating_point_v<sample>,
                      "the sub-band round trip is float and double only: fixed point is out of scope");

        /// The stages' ratios, reduced.
        static constexpr exact_ratio k_down_ratio = {Down::k_up, Down::k_down};
        static constexpr exact_ratio k_up_ratio   = {Up::k_up, Up::k_down};
        static_assert(k_down_ratio.num == 1 && k_up_ratio.den == 1 && k_up_ratio.num == k_down_ratio.den,
                      "the sub-band round trip needs inverse integer factors: Down converts by 1/D, Up by D/1");

        /// D, the factor.
        static constexpr std::size_t k_factor = static_cast<std::size_t>(k_down_ratio.den);

        /// The frame width every part runs, interleaved.
        static constexpr std::size_t k_channels = Down::channels();
        static_assert(Up::channels() == k_channels && Op::channels() == k_channels,
                      "the sub-band round trip's stages and op must run the same channel count");

        static_assert(detail::stages_delay_whole_if_static<Down, Up, k_factor>(),
                      "the sub-band round trip needs D * tau_down + tau_up to be a whole number of full-rate "
                      "samples: these stages' delays are compile-time constants and their sum is not");

        /// A sync_stage of ratio 1/1.
        static constexpr std::size_t k_up    = 1;
        static constexpr std::size_t k_down  = 1;
        static constexpr std::size_t k_block = 256; ///< full-rate frames per internal block

        /// Takes the stages and the op by value (move them in) and resets them.
        /// aligned_latency: the latency to align to (a composition's largest
        /// T); the instance's latency is the larger of it and its own T, so 0
        /// means its own. Rejects stages whose D tau_down + tau_up is not whole
        /// or whose change path falls behind (raise<std::invalid_argument>;
        /// supports() is the predicate).
        explicit subband_delta(Down down, Up up, Op op, std::size_t aligned_latency = 0,
                               subband_output output = subband_output::input_plus_change)
            : m_down(std::move(down))
            , m_up(std::move(up))
            , m_op(std::move(op))
            , m_output(output) {
            m_down.reset();
            m_up.reset();
            m_op.reset();
            const exact_ratio stages = stages_delay_frames(m_down, m_up);
            if (stages.den != 1) {
                raise<std::invalid_argument>("the sub-band round trip needs D * tau_down + tau_up to be a whole "
                                             "number of full-rate samples");
            }
            if (!lead_of(m_down, m_up, m_lead)) {
                raise<std::invalid_argument>("the sub-band round trip's change path falls behind its input");
            }
            m_round_trip            = static_cast<std::size_t>(stages.num) + k_factor * m_op.latency_frames();
            m_latency               = std::max(m_round_trip, aligned_latency);
            const std::size_t extra = m_latency - m_round_trip;
            m_extra_frames          = extra / k_factor;
            m_remainder             = extra % k_factor;

            const std::size_t low_cap = detail::max_frames_out<Down>(k_block);
            m_low.assign(low_cap * k_channels, sample{0});
            m_change.assign(low_cap * k_channels, sample{0});
            m_queue.assign((m_lead + m_remainder + detail::max_frames_out<Up>(low_cap)) * k_channels, sample{0});
            m_extra.assign(m_extra_frames * k_channels, sample{0});
            if (m_output == subband_output::input_plus_change) {
                m_delay.assign(m_latency * k_channels, sample{0});
            }
            m_zeros.assign(k_block * k_channels, sample{0});
            reset_lines();
        }

        /// D tau_down + tau_up in full-rate frames, exact: the stages' share
        /// of T.
        [[nodiscard]] static exact_ratio stages_delay_frames(const Down& down, const Up& up) noexcept {
            return exact_ratio{k_factor, 1} * detail::exact_delay_output_frames(down)
                   + detail::exact_delay_output_frames(up);
        }

        /// Whether the constructor accepts these stages, read from the state
        /// they are in (the reset state, as constructed): D tau_down + tau_up
        /// whole, and the change path leading.
        [[nodiscard]] static bool supports(const Down& down, const Up& up) noexcept {
            std::size_t surplus = 0;
            return stages_delay_frames(down, up).den == 1 && lead_of(down, up, surplus);
        }

        /// T = D tau_down + D L_op + tau_up from the parts, before construction
        /// (to align a composition's bands). @pre supports(down, up).
        [[nodiscard]] static std::size_t round_trip_latency(const Down& down, const Up& up, const Op& op) noexcept {
            const exact_ratio stages = stages_delay_frames(down, up);
            TAP_EXPECTS(stages.den == 1);
            return static_cast<std::size_t>(stages.num / stages.den) + k_factor * op.latency_frames();
        }

        /// Writes n frames: x[n - A] + change, or the change alone (the
        /// output mode). noexcept and allocation-free; in == out allowed.
        std::size_t process(const sample* in, std::size_t n, sample* out) noexcept {
            std::size_t done = 0;
            while (done < n) {
                const std::size_t take = std::min(n - done, k_block);
                run_block(in + done * k_channels, take, out + done * k_channels);
                done += take;
            }
            return n;
        }

        /// Frames process(n) writes: n.
        std::size_t outputs_for(std::size_t n) const noexcept { return n; }

        /// The instance's latency in full-rate frames: max(T, aligned_latency).
        std::size_t latency() const noexcept { return m_latency; }
        /// T = D tau_down + D L_op + tau_up, this band's own round trip.
        std::size_t round_trip_latency() const noexcept { return m_round_trip; }
        /// latency() as an exact_ratio of output frames (chain.h's form).
        exact_ratio latency_output_frames() const noexcept { return exact_ratio{m_latency, 1}; }
        /// The change queue's surplus bound: max over j of the frames the
        /// change path has written after j inputs, minus j (D - 1 for
        /// decimators on decimate.h's alignment).
        std::size_t lead() const noexcept { return m_lead; }
        /// The extra delay held at the low rate, floor((latency() - T) / D) frames.
        std::size_t extra_low_rate_frames() const noexcept { return m_extra_frames; }
        /// The rest of the extra delay, (latency() - T) mod D full-rate frames.
        std::size_t extra_full_rate_frames() const noexcept { return m_remainder; }
        /// What process() writes.
        subband_output output_mode() const noexcept { return m_output; }

        /// Zero input frames after which the output is silence, from any state
        /// (the file header's bound).
        std::size_t window_frames() const noexcept {
            static_assert(k_flushable, "the sub-band round trip's flush needs Down, Up and Op to report "
                                       "window_frames(), the zero input frames that drain each");
            const std::size_t low =
                static_cast<std::size_t>(m_op.window_frames()) + m_extra_frames + detail::drain_frames(m_up);
            const std::size_t change =
                detail::drain_frames(m_down) + detail::max_frames_in<Down>(low) + m_lead + m_remainder;
            return m_output == subband_output::input_plus_change ? std::max(m_latency, change) : change;
        }

        /// Frames flush() writes: window_frames().
        std::size_t flush_output_frames() const noexcept { return window_frames(); }

        /// End of stream: feeds window_frames() zero frames through process()
        /// and writes as many; zero padding by construction. reset() before
        /// reuse.
        std::size_t flush(sample* out) noexcept {
            const std::size_t w    = window_frames();
            std::size_t       done = 0;
            while (done < w) {
                const std::size_t take = std::min(w - done, k_block);
                run_block(m_zeros.data(), take, out + done * k_channels);
                done += take;
            }
            return w;
        }

        /// Back to the zero-primed state.
        void reset() noexcept {
            m_down.reset();
            m_up.reset();
            m_op.reset();
            reset_lines();
        }

        /// The channel count every part runs (k_channels).
        static constexpr std::size_t channels() noexcept { return k_channels; }

        /// The parts, read-only: they are driven only through the scaffold.
        const Down& down_stage() const noexcept { return m_down; }
        const Up&   up_stage() const noexcept { return m_up; }
        const Op&   low_rate_op() const noexcept { return m_op; }

      private:
        static constexpr bool k_flushable =
            detail::drains<Down>::k_value && detail::drains<Up>::k_value && detail::reports_window<Op>;

        /// The change path's lead over its input from the stages' state: false
        /// if it falls behind for some j in 1 .. 2D, else the largest surplus.
        static bool lead_of(const Down& down, const Up& up, std::size_t& surplus) noexcept {
            surplus = 0;
            for (std::size_t j = 1; j <= 2 * k_factor; ++j) {
                const std::size_t written = up.outputs_for(down.outputs_for(j));
                if (written < j) {
                    return false;
                }
                surplus = std::max(surplus, written - j);
            }
            return true;
        }

        /// Where the change is exactly zero the delayed input is written as it
        /// is (the identity's select; -0 + +0 would write +0).
        static sample add_change(sample delayed, sample change) noexcept {
            return change == sample{0} ? delayed : delayed + change;
        }

        void reset_lines() noexcept {
            std::fill(m_delay.begin(), m_delay.end(), sample{0});
            std::fill(m_extra.begin(), m_extra.end(), sample{0});
            m_delay_pos = 0;
            m_extra_pos = 0;
            // The full-rate rest of the extra delay: the change stream starts
            // with that many zero frames.
            std::fill(m_queue.begin(), m_queue.begin() + static_cast<std::ptrdiff_t>(m_remainder * k_channels),
                      sample{0});
            m_queued = m_remainder;
        }

        /// One block of at most k_block frames through both paths.
        void run_block(const sample* in, std::size_t frames, sample* out) noexcept {
            const std::size_t low = m_down.process(in, frames, m_low.data());
            m_op.process(m_low.data(), low, m_change.data());
            delay_low_rate(m_change.data(), low);
            TAP_EXPECTS((m_queued + detail::max_frames_out<Up>(low)) * k_channels <= m_queue.size());
            m_queued += m_up.process(m_change.data(), low, m_queue.data() + m_queued * k_channels);
            TAP_EXPECTS(m_queued >= frames); // the change path leads (checked at construction)
            const sample* change = m_queue.data();
            if (m_output == subband_output::input_plus_change) {
                add_to_delayed_input(in, frames, change, out);
            }
            else {
                std::copy(change, change + frames * k_channels, out);
            }
            // Keep the surplus (at most lead() + the remainder) for the inputs
            // it belongs to.
            m_queued -= frames;
            std::copy(m_queue.begin() + static_cast<std::ptrdiff_t>(frames * k_channels),
                      m_queue.begin() + static_cast<std::ptrdiff_t>((frames + m_queued) * k_channels), m_queue.begin());
        }

        /// The extra low-rate delay, in place on the change block.
        void delay_low_rate(sample* change, std::size_t frames) noexcept {
            std::size_t i = 0;
            while (i < frames && m_extra_frames > 0) {
                const std::size_t seg  = std::min(frames - i, m_extra_frames - m_extra_pos);
                sample*           line = m_extra.data() + m_extra_pos * k_channels;
                sample*           at   = change + i * k_channels;
                for (std::size_t j = 0; j < seg * k_channels; ++j) {
                    const sample held = line[j];
                    line[j]           = at[j];
                    at[j]             = held;
                }
                i += seg;
                m_extra_pos += seg;
                if (m_extra_pos == m_extra_frames) {
                    m_extra_pos = 0;
                }
            }
        }

        /// out = x[n - A] + change through the delay line; reads in[j] before
        /// writing out[j], so in == out is safe.
        void add_to_delayed_input(const sample* in, std::size_t frames, const sample* change, sample* out) noexcept {
            if (m_latency == 0) {
                for (std::size_t j = 0; j < frames * k_channels; ++j) {
                    out[j] = add_change(in[j], change[j]);
                }
                return;
            }
            std::size_t i = 0;
            while (i < frames) {
                const std::size_t seg  = std::min(frames - i, m_latency - m_delay_pos);
                sample*           line = m_delay.data() + m_delay_pos * k_channels;
                const std::size_t base = i * k_channels;
                for (std::size_t j = 0; j < seg * k_channels; ++j) {
                    const sample delayed = line[j];
                    line[j]              = in[base + j];
                    out[base + j]        = add_change(delayed, change[base + j]);
                }
                i += seg;
                m_delay_pos += seg;
                if (m_delay_pos == m_latency) {
                    m_delay_pos = 0;
                }
            }
        }

        Down                m_down;
        Up                  m_up;
        Op                  m_op;
        subband_output      m_output;
        std::size_t         m_round_trip   = 0; ///< T
        std::size_t         m_latency      = 0; ///< A = max(T, aligned_latency)
        std::size_t         m_extra_frames = 0; ///< floor((A - T) / D), at the low rate
        std::size_t         m_remainder    = 0; ///< (A - T) mod D, at the full rate
        std::size_t         m_lead         = 0;
        std::vector<sample> m_low;    ///< one block's low-rate frames from Down
        std::vector<sample> m_change; ///< the op's change for them
        std::vector<sample> m_queue;  ///< the change at the full rate, ahead of the input
        std::vector<sample> m_extra;  ///< the extra low-rate delay line
        std::vector<sample> m_delay;  ///< the input's delay line (input_plus_change)
        std::vector<sample> m_zeros;  ///< one block of silence, what flush feeds
        std::size_t         m_queued    = 0;
        std::size_t         m_delay_pos = 0;
        std::size_t         m_extra_pos = 0;
    };

} // namespace tap::dsp
