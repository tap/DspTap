/// @file chain.h
/// @brief The synchronous-stage concept and chain<>, the composition point below the engines.
// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// The SampleRateTap family's sync engines each convert by one ratio; a rate
// pair that needs several (96 -> 44.1 is a by-2 stage then the 147/160
// bridge) is a chain of stages, and the family rule is that a chain is
// written by the caller, as a type, never looked up from a rate pair. This
// header holds that composition so it can sit below every engine: a chain of
// a bridge converter and rational stages needs neither engine to name the
// other. Landed here first (rational's PLAN.md, R6) and proven on
// decimate.h, an existing primitive, before any new engine depends on it.
//
// The sync_stage concept is structural and minimal: a sample type, the
// rate ratio k_up / k_down as compile-time numbers, process(), outputs_for()
// and reset(), which is what decimate.h and the bridge converter already
// have. Latency composes when every stage reports its own.
//
// Contract, as numbers:
//   - process(in, n, out) feeds n input frames through every stage in order
//     and returns the frames the last stage wrote; it is noexcept and
//     allocation-free, working in chunks of k_block input frames through
//     scratch buffers sized at construction (every intermediate count is
//     bounded by n * k_up / k_down + 2 per stage, the stages' contract).
//     Bit-identical for any chunking of the same stream, because every
//     stage is.
//   - outputs_for(n) composes forward from the current position of every
//     stage: exactly what the next process(n) call returns.
//   - frames_needed(k): the smallest n with outputs_for(n) >= k, composed
//     backward stage by stage; a pull that delivers exactly frames_needed(k)
//     frames yields exactly k outputs (outputs_for(frames_needed(k)) == k,
//     and one frame fewer yields fewer).
//   - latency_output_frames(): the sum over stages of each stage's group
//     delay converted to the chain's output rate, an exact rational
//     (exact_ratio, reduced); latency_seconds(out_rate_hz) is its double.
//     Each stage reports latency_output_frames() (an exact_ratio in its own
//     output frames) or latency_input_samples() (an integer at its input
//     rate, decimate.h's form); the chain converts either.
//   - reset() resets every stage.
//
// Channels: every stage of a chain must agree on the frame width; the chain
// is told the count only to size its scratch (decimate.h is mono).

#pragma once

#include <concepts>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include "tap/dsp/sample_traits.h"

namespace tap::dsp {

    /// A reduced non-negative rational: numerator / denominator, denominator >= 1.
    struct exact_ratio {
        std::uint64_t num = 0;
        std::uint64_t den = 1;

        constexpr exact_ratio() = default;
        constexpr exact_ratio(std::uint64_t n, std::uint64_t d) noexcept
            : num(n)
            , den(d) {
            reduce();
        }

        constexpr exact_ratio operator+(exact_ratio o) const noexcept {
            return exact_ratio{num * o.den + o.num * den, den * o.den};
        }
        constexpr exact_ratio operator*(exact_ratio o) const noexcept { return exact_ratio{num * o.num, den * o.den}; }
        constexpr bool        operator==(const exact_ratio&) const noexcept = default;
        constexpr double      value() const noexcept { return static_cast<double>(num) / static_cast<double>(den); }

      private:
        constexpr void reduce() noexcept {
            const std::uint64_t g = std::gcd(num, den);
            if (g > 1) {
                num /= g;
                den /= g;
            }
            if (den == 0) {
                den = 1;
            }
        }
    };

    /// A synchronous rate-changing stage: the structural minimum a chain
    /// needs. k_up / k_down is the output / input frame ratio (L / M).
    template <typename T>
    concept sync_stage =
        requires(T s, const T& cs, const typename T::sample* in, typename T::sample* out, std::size_t n) {
            typename T::sample;
            requires sample_type<typename T::sample>;
            { T::k_up } -> std::convertible_to<std::size_t>;
            { T::k_down } -> std::convertible_to<std::size_t>;
            { s.process(in, n, out) } noexcept -> std::same_as<std::size_t>;
            { cs.outputs_for(n) } noexcept -> std::same_as<std::size_t>;
            { s.reset() } noexcept;
        };

    namespace detail {

        template <typename T>
        concept reports_latency_output = requires(const T& s) {
            { s.latency_output_frames() } -> std::same_as<exact_ratio>;
        };
        template <typename T>
        concept reports_latency_input = requires(const T& s) {
            { s.latency_input_samples() } -> std::convertible_to<std::uint64_t>;
        };

        /// A stage's group delay in its own output frames, from whichever
        /// form it reports.
        template <typename T>
        constexpr exact_ratio stage_latency_output_frames(const T& s) noexcept {
            if constexpr (reports_latency_output<T>) {
                return s.latency_output_frames();
            }
            else {
                static_assert(
                    reports_latency_input<T>,
                    "a chain's latency needs every stage to report latency_output_frames() or latency_input_samples()");
                return exact_ratio{static_cast<std::uint64_t>(s.latency_input_samples()), 1}
                       * exact_ratio{static_cast<std::uint64_t>(T::k_up), static_cast<std::uint64_t>(T::k_down)};
            }
        }

        /// Output frames a stage can produce for n input frames: floor(nL/M)
        /// plus one for the ceiling and one for the phase it starts in (the
        /// bound every stage's contract keeps: a decimator emits at most
        /// ceil(n/M), an interpolator exactly nL, a mixed L/M stage at most
        /// ceil(nL/M) + 1 from mid-superblock).
        template <typename T>
        constexpr std::size_t max_outputs_for(std::size_t n) noexcept {
            return (n * static_cast<std::size_t>(T::k_up)) / static_cast<std::size_t>(T::k_down) + 2;
        }

    } // namespace detail

    /// A chain of sync stages run in order; see the file header.
    template <sync_stage First, sync_stage... Rest>
    class chain {
      public:
        using sample = typename First::sample;
        static_assert((std::is_same_v<sample, typename Rest::sample> && ...),
                      "every stage of a chain must use the same sample type");

        static constexpr std::size_t k_stages = 1 + sizeof...(Rest);
        static constexpr std::size_t k_block  = 64; ///< input frames per internal chunk

        /// The chain's own ratio: the product of its stages', reduced.
        static constexpr exact_ratio k_ratio =
            (exact_ratio{First::k_up, First::k_down} * ... * exact_ratio{Rest::k_up, Rest::k_down});
        static constexpr std::size_t k_up   = static_cast<std::size_t>(k_ratio.num);
        static constexpr std::size_t k_down = static_cast<std::size_t>(k_ratio.den);

        /// Takes the stages by value (move them in); mono.
        explicit chain(First first, Rest... rest)
            : chain(std::size_t{1}, std::move(first), std::move(rest)...) {}

        /// Takes the stages by value; channels sizes the scratch only.
        chain(std::size_t channels, First first, Rest... rest)
            : m_stages(std::move(first), std::move(rest)...)
            , m_channels(channels) {
            // Scratch i holds stage i's output for a k_block-frame input chunk
            // (the last stage writes to the caller's buffer).
            std::size_t cap = k_block;
            [&]<std::size_t... I>(std::index_sequence<I...>) {
                ((cap = detail::max_outputs_for<std::tuple_element_t<I, stages_tuple>>(cap),
                  I + 1 < k_stages ? m_scratch[I].assign(cap * m_channels, sample_traits<sample>::silence()) : void()),
                 ...);
            }(std::make_index_sequence<k_stages>{});
        }

        /// Feeds n input frames through every stage; returns the frames the
        /// last stage wrote. noexcept and allocation-free.
        std::size_t process(const sample* in, std::size_t n, sample* out) noexcept {
            std::size_t produced = 0;
            while (n > 0) {
                const std::size_t take = n < k_block ? n : k_block;
                produced += run_chunk<0>(in, take, out + produced * m_channels);
                in += take * m_channels;
                n -= take;
            }
            return produced;
        }

        /// Frames process(n) will write from the current position.
        std::size_t outputs_for(std::size_t n) const noexcept {
            return [&]<std::size_t... I>(std::index_sequence<I...>) {
                ((n = std::get<I>(m_stages).outputs_for(n)), ...);
                return n;
            }(std::make_index_sequence<k_stages>{});
        }

        /// The smallest n with outputs_for(n) >= k, composed backward.
        std::size_t frames_needed(std::size_t k) const noexcept {
            return [&]<std::size_t... I>(std::index_sequence<I...>) {
                ((k = stage_frames_needed(std::get<k_stages - 1 - I>(m_stages), k)), ...);
                return k;
            }(std::make_index_sequence<k_stages>{});
        }

        /// The chain's group delay at its output rate, exact.
        exact_ratio latency_output_frames() const noexcept {
            exact_ratio total{0, 1};
            [&]<std::size_t... I>(std::index_sequence<I...>) {
                ((total = total + detail::stage_latency_output_frames(std::get<I>(m_stages)) * ratio_after<I>()), ...);
            }(std::make_index_sequence<k_stages>{});
            return total;
        }

        /// latency_output_frames() in seconds at the chain's output rate.
        double latency_seconds(double out_rate_hz) const noexcept {
            return latency_output_frames().value() / out_rate_hz;
        }

        void reset() noexcept {
            std::apply([](auto&... s) { (s.reset(), ...); }, m_stages);
        }

        std::size_t channels() const noexcept { return m_channels; }

        template <std::size_t I>
        auto& stage() noexcept {
            return std::get<I>(m_stages);
        }
        template <std::size_t I>
        const auto& stage() const noexcept {
            return std::get<I>(m_stages);
        }

      private:
        using stages_tuple = std::tuple<First, Rest...>;

        /// The product of the ratios of the stages after I: converts stage
        /// I's output frames to the chain's.
        template <std::size_t I>
        static constexpr exact_ratio ratio_after() noexcept {
            return [&]<std::size_t... J>(std::index_sequence<J...>) {
                exact_ratio r{1, 1};
                ((r = J > I ? r
                                  * exact_ratio{std::tuple_element_t<J, stages_tuple>::k_up,
                                                std::tuple_element_t<J, stages_tuple>::k_down}
                            : r),
                 ...);
                return r;
            }(std::make_index_sequence<k_stages>{});
        }

        template <std::size_t I>
        std::size_t run_chunk(const sample* in, std::size_t n, sample* out) noexcept {
            if constexpr (I + 1 == k_stages) {
                return std::get<I>(m_stages).process(in, n, out);
            }
            else {
                const std::size_t made = std::get<I>(m_stages).process(in, n, m_scratch[I].data());
                return run_chunk<I + 1>(m_scratch[I].data(), made, out);
            }
        }

        /// The smallest n with s.outputs_for(n) >= k for one stage: outputs_for
        /// is exact and monotone, so a bisection over [0, bound] finds it;
        /// the bound comes from the stage's ratio and is widened if the phase
        /// the stage is in needs more.
        template <typename T>
        static std::size_t stage_frames_needed(const T& s, std::size_t k) noexcept {
            if (k == 0) {
                return 0;
            }
            std::size_t hi = (k * static_cast<std::size_t>(T::k_down)) / static_cast<std::size_t>(T::k_up)
                             + static_cast<std::size_t>(T::k_down) + 2;
            while (s.outputs_for(hi) < k) {
                hi *= 2;
            }
            std::size_t lo = 0; // outputs_for(lo) < k
            while (hi - lo > 1) {
                const std::size_t mid = lo + (hi - lo) / 2;
                if (s.outputs_for(mid) >= k) {
                    hi = mid;
                }
                else {
                    lo = mid;
                }
            }
            return hi;
        }

        stages_tuple        m_stages;
        std::size_t         m_channels;
        std::vector<sample> m_scratch[k_stages > 1 ? k_stages - 1 : 1];
    };

} // namespace tap::dsp
