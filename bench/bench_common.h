/// @file bench_common.h
/// @brief Engine selector, deterministic input and the bit-exact fold shared by the FFT benches.
// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// Shared by bench/icount/icount_main.cpp (the instruction-count ratchet
// workloads), bench/size_probe.cpp (the .text probe) and bench/bench_fft.cpp
// (the host wall-clock microbenchmark) so all measure the same engine over
// the same input (bench/README.md).

#pragma once

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <string_view>
#include <type_traits>

#include "tap/dsp/fft.h"

// TAP_DSP_BENCH_ENGINE — the engine a scenario measures. A bare identifier,
// set by bench/CMakeLists.txt (cache variable of the same name; default
// reference_c) and stringized here so it is printed with every result:
//
//   reference_c  tap::dsp::basic_real_fft<Sample> as built today: the vendored
//                Ooura C (third_party/ooura/fftsg.c, fftsg_float.c). Where the
//                build routes float32 through an accelerated backend behind the
//                same class (TAP_DSP_FFT_CMSIS on the `m55` key, vDSP on Apple),
//                that backend is what this value measures; the printed
//                `backend=` field (backend_name below) says which.
//   split_radix  TODO(Stage 2a): the C++20 port, detail::split_radix_rdft<Sample>,
//                built beside the C. The port agent replaces the static_assert
//                below with the second alias and makes bench/icount/CMakeLists.txt
//                build every float scenario twice (docs/audit-fft-and-code-smells.md,
//                Part 11: the C-vs-port ratio is the Stage 2b gate).
//
// A documented macro, on purpose — not a registry. Once Stage 4 makes the
// engine an explicit class parameter this selector becomes that parameter.
#ifndef TAP_DSP_BENCH_ENGINE
#define TAP_DSP_BENCH_ENGINE reference_c
#endif
#define TAP_DSP_BENCH_STRINGIZE_(x) #x
#define TAP_DSP_BENCH_STRINGIZE(x) TAP_DSP_BENCH_STRINGIZE_(x)

namespace tap::dsp::bench {

    constexpr const char* k_engine_name = TAP_DSP_BENCH_STRINGIZE(TAP_DSP_BENCH_ENGINE);
    static_assert(std::string_view{k_engine_name} == "reference_c",
                  "TAP_DSP_BENCH_ENGINE: only reference_c exists until the Stage 2a port lands split_radix");

    /// The transform under test for the selected engine.
    template <typename Sample>
    using fft_under_test = tap::dsp::basic_real_fft<Sample>;

    /// Which backend basic_real_fft<Sample> was built over, from the macros
    /// fft.h switches on. The accelerated backends apply to float only;
    /// double is always the Ooura C.
    template <typename Sample>
    constexpr const char* backend_name() noexcept {
        if constexpr (std::is_same_v<Sample, float>) {
#if defined(TAP_DSP_FFT_CMSIS)
            return "cmsis";
#elif defined(TAP_DSP_FFT_ACCELERATE)
            return "accelerate";
#else
            return "ooura";
#endif
        }
        else {
            return "ooura";
        }
    }

    /// xorshift32 (Marsaglia 2003), the family's deterministic source: no
    /// <random>, so counts cannot drift with a toolchain's libstdc++, and no
    /// double anywhere, so the float scenarios never touch libgcc's soft
    /// double on the M4 soft-float leg. fill() lands in the CLOSED interval
    /// [-0.5, +0.5]: static_cast<float>(0xFFFFFFFF) rounds up to 2^32, so
    /// the top value is reached exactly.
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
                out[i] = (static_cast<Sample>(next()) / Sample(4294967296.0)) - Sample(0.5);
            }
        }

      private:
        std::uint32_t m_state;
    };

    /// FNV-1a-64 (Fowler–Noll–Vo; offset basis and prime as published) over
    /// the bit pattern of every sample, one 32- or 64-bit word per step
    /// rather than one octet: h ^= bits(v[i]); h *= prime. Integer only, so
    /// it is exact and order-sensitive — a 1-ulp change in any single output
    /// changes the printed value, which a floating running sum cannot promise
    /// (it absorbs differences below the accumulator's ulp). This is what
    /// makes the printed checksum a fingerprint of the engine's output, the
    /// C-vs-port comparison Stage 2b needs on the QEMU legs, and not only a
    /// dead-code-elimination guard.
    inline constexpr std::uint64_t k_fnv1a64_offset = 0xcbf29ce484222325ull;
    inline constexpr std::uint64_t k_fnv1a64_prime  = 0x100000001b3ull;

    template <typename Sample>
    std::uint64_t fold(std::uint64_t h, const Sample* v, std::size_t n) noexcept {
        using bits_type = std::conditional_t<sizeof(Sample) == sizeof(std::uint32_t), std::uint32_t, std::uint64_t>;
        static_assert(sizeof(bits_type) == sizeof(Sample), "fold expects a 32- or 64-bit sample");
        for (std::size_t i = 0; i < n; ++i) {
            bits_type bits = 0;
            std::memcpy(&bits, &v[i], sizeof bits);
            h ^= bits;
            h *= k_fnv1a64_prime;
        }
        return h;
    }

} // namespace tap::dsp::bench
