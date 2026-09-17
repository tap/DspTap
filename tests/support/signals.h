/// @file signals.h
/// @brief Shared test-signal generators for the DspTap test battery.
// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// One xorshift32, one random_signal, one tone synthesizer and one dB helper,
// so a new test file does not grow its own copy (three already exist in
// test_fft.cpp, test_fft_backend.cpp and test_nn.cpp; migrating them is the
// Stage 6 hygiene item in docs/audit-fft-and-code-smells.md, not this file's
// job). Everything here is deterministic — fixed seeds, no wall clock, no
// filesystem — so it can run unchanged on the bare-metal QEMU legs.

#pragma once

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numbers>
#include <vector>

namespace tap::dsp::test {

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

    /// n samples uniform in [-amplitude, amplitude), from a fixed seed.
    template <typename Sample>
    std::vector<Sample> random_signal(std::size_t n, std::uint32_t seed, double amplitude = 1.0) {
        xorshift32          rng(seed);
        std::vector<Sample> x(n);
        for (auto& v : x) {
            v = static_cast<Sample>(amplitude * rng.next_unit());
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
            x[j]               = static_cast<Sample>(amplitude * std::cos(2.0 * std::numbers::pi * turns + phase));
        }
        return x;
    }

    /// 20*log10(ratio) of an amplitude ratio; -inf for zero, as std::log10 gives.
    inline double db_from_ratio(double ratio) noexcept {
        return 20.0 * std::log10(ratio);
    }

} // namespace tap::dsp::test
