/// @file signals.h
/// @brief Shared test-signal generators for the DspTap test battery.
// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// One xorshift32, one random_signal and one tone synthesizer, so a new test
// file does not grow its own copy (three already exist in test_fft.cpp,
// test_fft_backend.cpp and test_nn.cpp; migrating them, and adding the dB
// helper Part 9 lists once a caller exists, is the Stage 6 hygiene item in
// docs/audit-fft-and-code-smells.md, not this file's job). Everything here is deterministic — fixed seeds, no wall
// clock, no filesystem — so it can run unchanged on the bare-metal QEMU legs.
//
// sample_scale<Sample> is how a generator lands in a profile's sample type:
// the identity cast for float and double, and round-to-nearest with saturation
// into Q0.15 / Q0.31 for the fixed-point profiles (Stage 3b), so the same
// seed produces the same signal in every profile up to the profile's own
// quantisation. Four files share it (test_fft.cpp, test_fft_fixed.cpp,
// test_fft_oracle.cpp, test_fft_rt.cpp), which is the bar for a helper here.

#pragma once

#include <cmath>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numbers>
#include <type_traits>
#include <vector>

#include "tap/dsp/sample_traits.h"

namespace tap::dsp::test {

    /// How a value in the double domain (fractions of full scale, 1.0 = full
    /// scale) crosses into and out of a profile's sample type.
    ///   - float, double : the plain cast; k_lsb is 0 and k_full_scale is 1.
    ///   - int16_t (Q0.15), int32_t (Q0.31): from_double rounds half away
    ///     from zero and saturates (the substrate's round_sat, the same
    ///     rounding the coefficient generators use); k_lsb is 2^-15 / 2^-31
    ///     and k_full_scale the largest representable positive value,
    ///     1 - k_lsb. Full-scale negative (-1.0) is INT_MIN exactly.
    template <typename Sample>
    struct sample_scale;

    template <std::floating_point F>
    struct sample_scale<F> {
        static constexpr int    k_frac_bits  = 0;
        static constexpr double k_lsb        = 0.0;
        static constexpr double k_full_scale = 1.0;
        static constexpr double to_double(F v) noexcept { return static_cast<double>(v); }
        static constexpr F      from_double(double v) noexcept { return static_cast<F>(v); }
    };

    template <std::integral I>
    struct sample_scale<I> {
        static_assert(std::is_same_v<I, std::int16_t> || std::is_same_v<I, std::int32_t>,
                      "the fixed-point profiles are Q0.15 (int16_t) and Q0.31 (int32_t)");
        static constexpr int    k_frac_bits  = std::numeric_limits<I>::digits; // 15 or 31
        static constexpr double k_scale      = static_cast<double>(std::int64_t{1} << k_frac_bits);
        static constexpr double k_lsb        = 1.0 / k_scale;
        static constexpr double k_full_scale = 1.0 - k_lsb;
        static constexpr double to_double(I v) noexcept { return static_cast<double>(v) / k_scale; }
        static constexpr I      from_double(double v) noexcept { return tap::dsp::detail::round_sat<I>(v * k_scale); }
    };

    static_assert(sample_scale<std::int16_t>::from_double(1.0) == std::numeric_limits<std::int16_t>::max());
    static_assert(sample_scale<std::int16_t>::from_double(-1.0) == std::numeric_limits<std::int16_t>::min());
    static_assert(sample_scale<std::int32_t>::from_double(0.5) == std::int32_t{1} << 30);
    static_assert(sample_scale<std::int32_t>::to_double(std::int32_t{1} << 30) == 0.5);
    static_assert(sample_scale<float>::from_double(0.25) == 0.25f);

    /// Marsaglia's xorshift32 (Journal of Statistical Software 8(14), 2003),
    /// the generator every existing copy in this repo already uses. The state
    /// must be non-zero; a zero seed is replaced by a fixed non-zero constant.
    class xorshift32 {
      public:
        explicit xorshift32(std::uint32_t seed) noexcept
            : m_s(seed != 0u ? seed : 0x9E3779B9u) {}

        /// Next raw 32-bit state.
        std::uint32_t next_u32() noexcept {
            m_s ^= m_s << 13;
            m_s ^= m_s >> 17;
            m_s ^= m_s << 5;
            return m_s;
        }

        /// Uniform in [-1, 1), computed in double so float and double signals
        /// drawn from the same seed are the same values up to the final cast.
        double next_unit() noexcept { return static_cast<double>(next_u32()) / 2147483648.0 - 1.0; }

      private:
        std::uint32_t m_s;
    };

    /// n samples uniform in [-amplitude, amplitude), from a fixed seed, landed
    /// in the profile through sample_scale (a plain cast for float/double).
    template <typename Sample>
    std::vector<Sample> random_signal(std::size_t n, std::uint32_t seed, double amplitude = 1.0) {
        xorshift32          rng(seed);
        std::vector<Sample> x(n);
        for (auto& v : x) {
            v = sample_scale<Sample>::from_double(amplitude * rng.next_unit());
        }
        return x;
    }

    /// amplitude * cos(2*pi*bin*j/n + phase) for j in [0, n).
    ///
    /// The angle is formed as 2*pi * fmod(bin*j, n) / n, not as (2*pi*bin/n) * j:
    /// with a rounded per-sample increment the phase error grows like j*ulp,
    /// and at n = 65536 a supposedly on-bin cosine is ~1e-11 off, a thousand
    /// times the engine's own error — which is what a closed-form oracle would
    /// then wrongly report. fmod is exact, bin*j is exact for integer bins below
    /// 2^53, and the one division by a power-of-two n is exact, so an integer-bin
    /// tone here is the sampled cosine to one rounding of cos() per sample.
    template <typename Sample>
    std::vector<Sample> tone(std::size_t n, double bin, double amplitude, double phase = 0.0) {
        std::vector<Sample> x(n);
        const double        period = static_cast<double>(n);
        for (std::size_t j = 0; j < n; ++j) {
            const double turns = std::fmod(bin * static_cast<double>(j), period) / period;
            x[j] = sample_scale<Sample>::from_double(amplitude * std::cos(2.0 * std::numbers::pi * turns + phase));
        }
        return x;
    }

} // namespace tap::dsp::test
