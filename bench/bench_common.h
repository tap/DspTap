/// @file bench_common.h
/// @brief Engine selector and deterministic input shared by the two FFT benches.
// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// Shared by bench/icount/icount_main.cpp (the instruction-count ratchet
// workloads) and bench/bench_fft.cpp (the host wall-clock microbenchmark) so
// both measure the same engine over the same input (bench/README.md).

#pragma once

#include <cstddef>
#include <cstdint>
#include <string_view>

#include "tap/dsp/fft.h"

// TAP_DSP_BENCH_ENGINE — the engine a scenario measures. A bare identifier,
// set by bench/CMakeLists.txt (cache variable of the same name; default
// reference_c) and stringized here so it is printed with every result:
//
//   reference_c  tap::dsp::basic_real_fft<Sample> as built today: the vendored
//                Ooura C (third_party/ooura/fftsg.c, fftsg_float.c). Where the
//                build routes float32 through an accelerated backend behind the
//                same class (TAP_DSP_FFT_CMSIS on the `m55` key, vDSP on Apple),
//                that backend is what this value measures; the ratchet key says
//                which (bench/README.md).
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

    constexpr std::string_view k_engine_name = TAP_DSP_BENCH_STRINGIZE(TAP_DSP_BENCH_ENGINE);
    static_assert(k_engine_name == "reference_c",
                  "TAP_DSP_BENCH_ENGINE: only reference_c exists until the Stage 2a port lands split_radix");

    /// The transform under test for the selected engine.
    template <typename Sample>
    using fft_under_test = tap::dsp::basic_real_fft<Sample>;

    /// xorshift32 (Marsaglia 2003), the family's deterministic source: no
    /// <random>, so counts cannot drift with a toolchain's libstdc++, and no
    /// double anywhere, so the float scenarios never touch libgcc's soft
    /// double on the M4 soft-float leg. Uniform in [-0.5, 0.5).
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

    /// Fold every element into a running sum of the same type: defeats
    /// dead-code elimination and pins cross-run determinism (the printed
    /// checksum must be bit-identical between two runs of the same binary).
    template <typename Sample>
    Sample fold(Sample acc, const Sample* v, std::size_t n) noexcept {
        for (std::size_t i = 0; i < n; ++i) {
            acc += v[i];
        }
        return acc;
    }

} // namespace tap::dsp::bench
