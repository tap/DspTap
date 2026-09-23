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

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <string_view>
#include <type_traits>
#include <vector>

#include "tap/dsp/fft.h"

// TAP_DSP_BENCH_ENGINE — the engine a scenario measures. A bare identifier,
// set by bench/CMakeLists.txt (cache variable of the same name; default
// basic_real_fft) and stringized here so it is printed with every result:
//
//   basic_real_fft  tap::dsp::basic_real_fft<Sample> as built: what ships.
//                   Since Stage 2b that is the C++20 split-radix engine
//                   (include/tap/dsp/fft/split_radix.h) for double and, where
//                   no backend define is active, for float; where the build
//                   routes float32 through an accelerated backend behind the
//                   same class (TAP_DSP_FFT_CMSIS on the `m55` key, vDSP on
//                   Apple), that backend is what this value measures. The
//                   printed `backend=` field (backend_name below) says which.
//                   The bare-keyed scenarios, i.e. the gated ones, measure this.
//   reference_c     the vendored Ooura C called directly (rdft / rdft_f from
//                   third_party/ooura/fftsg.c and fftsg_float.c, through the
//                   reference_c_bench_adapter below, which presents
//                   basic_real_fft's out-of-place forward()/inverse() surface
//                   with the same copy loop and the same 2/N arithmetic, so
//                   the two engines' counts differ only by the transform).
//                   INFORMATIONAL from Stage 2b until Stage 2c retires the C:
//                   bench/icount/CMakeLists.txt builds every float scenario a
//                   second time against it under a `_c` key suffix, which
//                   scripts/icount.py reports beside the gated sibling with
//                   the ratio and whether the two output checksums agree, and
//                   never gates (bench/README.md). Before 2b this value named
//                   basic_real_fft itself, which then WAS the C; after the
//                   flip that spelling would have been a lie, so the C is
//                   measured by name and the shipping class by its own.
//
// A documented macro, on purpose — not a registry. Once Stage 4 makes the
// engine an explicit class parameter this selector becomes that parameter.
#ifndef TAP_DSP_BENCH_ENGINE
#define TAP_DSP_BENCH_ENGINE basic_real_fft
#endif
#define TAP_DSP_BENCH_STRINGIZE_(x) #x
#define TAP_DSP_BENCH_STRINGIZE(x) TAP_DSP_BENCH_STRINGIZE_(x)

namespace tap::dsp::bench {

    constexpr const char* k_engine_name           = TAP_DSP_BENCH_STRINGIZE(TAP_DSP_BENCH_ENGINE);
    constexpr bool        k_engine_is_reference_c = std::string_view{k_engine_name} == "reference_c";
    static_assert(k_engine_is_reference_c || std::string_view{k_engine_name} == "basic_real_fft",
                  "TAP_DSP_BENCH_ENGINE must be basic_real_fft (what ships) or reference_c (the vendored C, "
                  "informational until Stage 2c)");

    namespace detail {
        inline void reference_rdft(int n, int isgn, double* a, int* ip, double* w) {
            rdft(n, isgn, a, ip, w);
        }
        inline void reference_rdft(int n, int isgn, float* a, int* ip, float* w) {
            rdft_f(n, isgn, a, ip, w);
        }
    } // namespace detail

    /// The vendored C behind basic_real_fft's out-of-place surface, for
    /// like-for-like counts: the workspace geometry readme.txt prescribes
    /// (ip: 2 + sqrt(n/2), w: n/2), the first-call protocol (ip[0] = 0, so
    /// the C builds its tables inside the first transform, as it always did),
    /// and copy() and the 2/N scaling exactly as basic_real_fft::forward /
    /// ::inverse write them, so a C-vs-shipping delta is the transform's
    /// alone. A bench adapter, not a consumer surface; it goes with the C at
    /// Stage 2c. rdft / rdft_f are the extern "C" declarations fft.h still
    /// carries, resolved from the tap_dsp_fft library tap::dsp links.
    template <typename Sample>
    class reference_c_bench_adapter {
      public:
        explicit reference_c_bench_adapter(std::size_t n)
            : m_size(static_cast<int>(n))
            , m_ip(2 + static_cast<std::size_t>(std::sqrt(static_cast<double>(n) / 2.0)) + 1, 0)
            , m_w(n / 2, Sample(0)) {
            m_ip[0] = 0;
        }

        std::size_t size() const noexcept { return static_cast<std::size_t>(m_size); }

        void forward_inplace(Sample* data) noexcept {
            detail::reference_rdft(m_size, 1, data, m_ip.data(), m_w.data());
        }
        void inverse_inplace(Sample* data) noexcept {
            detail::reference_rdft(m_size, -1, data, m_ip.data(), m_w.data());
        }

        void forward(const Sample* input, Sample* output) noexcept {
            copy(input, output);
            forward_inplace(output);
        }

        void inverse(const Sample* input, Sample* output) noexcept {
            copy(input, output);
            inverse_inplace(output);
            const Sample scale = Sample(2) / static_cast<Sample>(m_size);
            for (int i = 0; i < m_size; ++i) {
                output[i] *= scale;
            }
        }

      private:
        void copy(const Sample* input, Sample* output) noexcept {
            if (input != output) {
                for (int i = 0; i < m_size; ++i) {
                    output[i] = input[i];
                }
            }
        }

        int                 m_size;
        std::vector<int>    m_ip;
        std::vector<Sample> m_w;
    };

    /// The transform under test for the selected engine.
    template <typename Sample>
    using fft_under_test = std::conditional_t<k_engine_is_reference_c, reference_c_bench_adapter<Sample>,
                                              tap::dsp::basic_real_fft<Sample>>;

    /// Which backend the transform under test runs on. For basic_real_fft:
    /// what the class was built over, from the macros fft.h switches on (the
    /// accelerated backends apply to float only; double is always the
    /// split-radix engine). For reference_c: the vendored C, which has no
    /// backend.
    template <typename Sample>
    constexpr const char* backend_name() noexcept {
        if constexpr (k_engine_is_reference_c) {
            return "ooura_c";
        }
        else if constexpr (std::is_same_v<Sample, float>) {
#if defined(TAP_DSP_FFT_CMSIS)
            return "cmsis";
#elif defined(TAP_DSP_FFT_ACCELERATE)
            return "accelerate";
#else
            return "split_radix";
#endif
        }
        else {
            return "split_radix";
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
