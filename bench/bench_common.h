/// @file bench_common.h
/// @brief Deterministic input and the bit-exact fold shared by the FFT benches.
// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// Shared by bench/icount/icount_main.cpp (the instruction-count ratchet
// workloads), bench/size_probe.cpp (the .text probe) and bench/bench_fft.cpp
// (the host wall-clock microbenchmark) so all measure the same transform
// over the same input (bench/README.md).

#pragma once

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <type_traits>

#include "tap/dsp/fft.h"

namespace tap::dsp::bench {

    /// What every bench binary measures: tap::dsp::basic_real_fft<Sample> as
    /// built, i.e. over the engine the build selected as that profile's
    /// default (Stage 4: default_real_fft_engine_t<Sample> for the floating
    /// profiles) — the split-radix engine (include/tap/dsp/fft/split_radix.h)
    /// for double and, where no backend define is active, for float; the
    /// accelerated engine behind the same class where the build routes
    /// float32 through one (TAP_DSP_FFT_CMSIS on the `m55` key, vDSP on
    /// Apple); the int32 fixed-point kernel (fft/fixed_point.h) for the Q15
    /// and Q31 profiles. The printed `backend=` field (backend_name below)
    /// says which.
    ///
    /// From Stage 2a until Stage 2c this header also carried an engine
    /// selector (TAP_DSP_BENCH_ENGINE) and an adapter that presented the
    /// vendored Ooura C behind the class's surface, so every float scenario
    /// could be counted twice and the job could print the C/port ratio and
    /// whether the two output checksums agreed. Stage 2c retired the C from
    /// the shipping tree (docs/audit-fft-and-code-smells.md, Part 3; the
    /// reference copy lives under tests/reference/ooura/ for the parity gate
    /// alone) and the pair with it: the bench measures only what ships. The
    /// engine name below is kept as the constant it always printed for the
    /// shipping class, so the DONE line — and with it the recorded counts,
    /// which include the print — did not move at 2c, nor at Stage 4: the
    /// engine became an explicit class parameter there, but the printed
    /// field stays the class's name (the `backend=` field already names the
    /// engine the default resolved to: split_radix / cmsis / accelerate /
    /// fixed_point) because the recorded baselines include the print and a
    /// longer string is a counted change (the 2b re-record measured +75
    /// instructions for three characters, bench/README.md).
    constexpr const char* k_engine_name = "basic_real_fft";

    /// The transform under test.
    template <typename Sample>
    using fft_under_test = tap::dsp::basic_real_fft<Sample>;

    /// Which backend the transform under test runs on: what the class was
    /// built over, from the macros fft.h switches on (the accelerated
    /// backends apply to float only; double is always the split-radix
    /// engine, Q15 and Q31 always the fixed-point kernel).
    template <typename Sample>
    constexpr const char* backend_name() noexcept {
        if constexpr (std::is_same_v<Sample, float>) {
#if defined(TAP_DSP_FFT_CMSIS)
            return "cmsis";
#elif defined(TAP_DSP_FFT_ACCELERATE)
            return "accelerate";
#else
            return "split_radix";
#endif
        }
        else if constexpr (std::is_integral_v<Sample>) {
            return "fixed_point";
        }
        else {
            return "split_radix";
        }
    }

    /// xorshift32 (Marsaglia 2003), the family's deterministic source: no
    /// <random>, so counts cannot drift with a toolchain's libstdc++, and no
    /// double anywhere, so the float and fixed-point scenarios never touch
    /// libgcc's soft double on the M4 soft-float leg. For the floating
    /// profiles fill() lands in the CLOSED interval [-0.5, +0.5]:
    /// static_cast<float>(0xFFFFFFFF) rounds up to 2^32, so the top value is
    /// reached exactly. For the fixed-point profiles it is the raw state
    /// narrowed to the sample width, i.e. uniform over the whole Q0.15 /
    /// Q0.31 range including the rails — the saturation-free contract holds
    /// for any input under scaling::fixed (fft.h), so full scale is the
    /// honest workload — with no floating-point conversion anywhere.
    class xorshift32 {
      public:
        constexpr explicit xorshift32(std::uint32_t seed) noexcept
            : m_state(seed) {}

        constexpr std::uint32_t next() noexcept {
            m_state ^= m_state << 13;
            m_state ^= m_state >> 17;
            m_state ^= m_state << 5;
            return m_state;
        }

        template <typename Sample>
        void fill(Sample* out, std::size_t n) noexcept {
            for (std::size_t i = 0; i < n; ++i) {
                if constexpr (std::is_integral_v<Sample>) {
                    // The top sizeof(Sample) bytes of the state, reinterpreted
                    // as two's complement: -2^(b-1) .. 2^(b-1)-1, uniform.
                    using bits_type = std::make_unsigned_t<Sample>;
                    const auto bits = static_cast<bits_type>(next() >> (32 - 8 * sizeof(Sample)));
                    out[i]          = static_cast<Sample>(bits);
                }
                else {
                    out[i] = (static_cast<Sample>(next()) / Sample(4294967296.0)) - Sample(0.5);
                }
            }
        }

      private:
        std::uint32_t m_state;
    };

    /// FNV-1a-64 (Fowler–Noll–Vo; offset basis and prime as published) over
    /// the bit pattern of every sample, one 16-, 32- or 64-bit word per step
    /// rather than one octet: h ^= bits(v[i]); h *= prime. Integer only, so
    /// it is exact and order-sensitive — a 1-ulp change in any single output
    /// changes the printed value, which a floating running sum cannot promise
    /// (it absorbs differences below the accumulator's ulp). This is what
    /// makes the printed checksum a fingerprint of the engine's output (the
    /// C-vs-port comparison Stage 2b was judged on, on the QEMU legs), and
    /// not only a dead-code-elimination guard. The fixed-point scenarios
    /// fold the returned exponent through the same function (as one int
    /// word) after each block, so a transform that returned the right bits
    /// under the wrong scale would change the checksum too.
    inline constexpr std::uint64_t k_fnv1a64_offset = 0xcbf29ce484222325ull;
    inline constexpr std::uint64_t k_fnv1a64_prime  = 0x100000001b3ull;

    template <typename Sample>
    std::uint64_t fold(std::uint64_t h, const Sample* v, std::size_t n) noexcept {
        using bits_type = std::conditional_t<
            sizeof(Sample) == sizeof(std::uint16_t), std::uint16_t,
            std::conditional_t<sizeof(Sample) == sizeof(std::uint32_t), std::uint32_t, std::uint64_t>>;
        static_assert(sizeof(bits_type) == sizeof(Sample), "fold expects a 16-, 32- or 64-bit sample");
        for (std::size_t i = 0; i < n; ++i) {
            bits_type bits = 0;
            std::memcpy(&bits, &v[i], sizeof bits);
            h ^= bits;
            h *= k_fnv1a64_prime;
        }
        return h;
    }

} // namespace tap::dsp::bench
