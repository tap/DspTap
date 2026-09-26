/// @file cmsis.h
/// @brief CMSIS-DSP Helium float32 real FFT engine, re-presenting the split-radix contract.
// SPDX-License-Identifier: MIT
// Copyright 2025-2026 Timothy Place and the DspTap contributors.
//
// One of the two accelerated float32 engines behind tap::dsp::basic_real_fft
// (Stage 4 of docs/audit-fft-and-code-smells.md: the engine is a template
// parameter of the class; this header is included by whoever selects this
// engine — fft.h under TAP_DSP_FFT_CMSIS, where it is the float default, or a
// translation unit that names detail::cmsis_real_fft_f32 explicitly). It
// needs "arm_math.h" on the include path and the CMSIS-DSP objects at link
// time (the tap_dsp_fft library the root CMakeLists.txt builds under
// TAP_DSP_FFT_CMSIS), i.e. an Arm target with Helium (the root CMakeLists.txt
// defaults that option ON only when the compiler defines __ARM_FEATURE_MVE
// with the floating-point bit); it is not a host engine and nothing hosted
// compiles it.

#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include "arm_math.h"
#include "tap/dsp/detail/expects.h"

namespace tap::dsp::detail {

    /// Wraps CMSIS-DSP's radix-4/8 Helium real FFT to reproduce the
    /// split-radix engine's exact float32 contract (Ooura's). CMSIS uses the
    /// engineering convention exp(-i2*pi/N) and a 1/N-normalized inverse;
    /// the contract is exp(+i2*pi/N) and an unnormalized inverse (caller
    /// applies 2/N). Reconciled by conjugating the imaginary bins on every
    /// transform and scaling the inverse by N/2. The pre-Stage-4 de-risk
    /// measurement of the reconciliation was <2e-7 relative error at N = 512
    /// and 2048; what the battery PINS on the Cortex-M55 QEMU leg is
    /// 5e-6 x peak per bin against the split-radix engine and a 2e-5
    /// absolute round trip, at the certified geometries
    /// (`fft_backend_parity/cmsis.ForwardMatchesOoura`,
    /// `RoundTripReproducesInput`, tests/test_fft_backend.cpp), plus
    /// bit-stability across buffer addresses and tonal accuracy at 512 /
    /// 2048 / 4096.
    ///
    /// Engine contract numbers (what basic_real_fft reads from every engine):
    ///   - Size range k_min_size = 32 … k_max_size = 4096, a power of two.
    ///     This is CMSIS-DSP's own table: arm_rfft_fast_init_f32 dispatches
    ///     on N through a switch over exactly {32, 64, …, 4096} and returns
    ///     ARM_MATH_ARGUMENT_ERROR for anything else, leaving the instance
    ///     uninitialized. Until Stage 4 the wrapper ignored that status and
    ///     fft.h documented "a power of two >= 4" for every profile, so N < 32
    ///     or N > 4096 under TAP_DSP_FFT_CMSIS was undefined behaviour — a
    ///     HardFault at N = 4 in the first QEMU run of the Stage 2a oracle on
    ///     the M55 leg is how it was found (tap/DspTap#24; audit Part 13). The constructor states the range as a
    ///     precondition (TAP_EXPECTS: a debug assertion, STYLE.md §4), never truncates a size above k_max_size into the
    ///     library's uint16_t argument (65536 would arrive as 0), and checks the init status the same way. RELEASE-MODE
    ///     BEHAVIOUR, DECIDED: the check evaluates to nothing, so a size outside the range remains a precondition
    ///     violation — the instance stays zero-initialized and every transform is undefined behaviour. No fault is
    ///     promised, and mostly none happens: per the final audit's reading of the library (2026-09-26), it computes
    ///     fftLen - 1 = -1 from the zeroed instance, skips its loops and reads
    ///     its coefficient table through a null pointer, and address 0 is
    ///     readable memory on the M55 (ITCM on mps3-an547; the vector table on
    ///     Cortex-M parts generally). Measured under QEMU in Release and
    ///     MinSizeRel: N = 16 and 8192 return a wrong spectrum and nothing
    ///     faults; N = 4 returns a wrong spectrum and corrupts the heap, so a
    ///     later, unrelated call fails — a HardFault in one probe, newlib's
    ///     "Balloc succeeded" assertion in the next printf in another. Real
    ///     silicon is not measured. Deliberately no defined fallback: a branch
    ///     in the transforms would cost the hot path on every call for a case
    ///     the precondition excludes, and would hand a consumer an object
    ///     that silently transforms nothing.
    ///     basic_real_fft<float, cmsis_real_fft_f32>::supports_size(n) is the
    ///     mandatory gate wherever N comes from configuration (the capi's
    ///     dsptap_fft_create applies it per profile; MuTap's config path must,
    ///     docs/fft-design.md "MuTap bump checklist"); its values for this
    ///     engine are pinned on the M55 leg (tests/test_fft_engine.cpp).
    ///   - k_is_shareable = false: the transform runs through the per-object
    ///     scratch buffer (arm_rfft_fast_f32 is out of place), so two threads
    ///     may not transform through one object concurrently; the transforms
    ///     are non-const accordingly.
    ///   - Everything is allocated in the constructor; the transforms are
    ///     noexcept and allocate nothing (the vendor layer is not audited by
    ///     tests/test_fft_rt.cpp's operator-new guard, only the wrapper).
    class cmsis_real_fft_f32 {
      public:
        static constexpr std::size_t k_min_size     = 32;
        static constexpr std::size_t k_max_size     = 4096;
        static constexpr bool        k_is_shareable = false;

        /// @pre n is a power of two in [k_min_size, k_max_size] (TAP_EXPECTS).
        explicit cmsis_real_fft_f32(std::size_t n)
            : m_scratch(n, 0.0f)
            , m_n(static_cast<int>(n)) {
            TAP_EXPECTS(n >= k_min_size && n <= k_max_size && (n & (n - 1)) == 0);
            // The range check guards the narrowing: the library takes a uint16_t and
            // 65536 would arrive as 0, so an out-of-range n never reaches it.
            const arm_status status =
                n <= k_max_size ? arm_rfft_fast_init_f32(&m_inst, static_cast<uint16_t>(n)) : ARM_MATH_ARGUMENT_ERROR;
            TAP_EXPECTS(status == ARM_MATH_SUCCESS);
            static_cast<void>(status); // the release build evaluates TAP_EXPECTS to nothing
        }

        [[nodiscard]] std::size_t size() const noexcept { return static_cast<std::size_t>(m_n); }

        void forward_inplace(float* a) noexcept {
            arm_rfft_fast_f32(&m_inst, a, m_scratch.data(), 0);
            // Copy back, conjugating imaginary bins into the exp(+i) convention.
            a[0] = m_scratch[0]; // DC (real)
            a[1] = m_scratch[1]; // Nyquist (real)
            for (int k = 2; k < m_n; ++k) {
                a[k] = (k & 1) ? -m_scratch[k] : m_scratch[k];
            }
        }

        void inverse_inplace(float* a) noexcept {
            // a is a packed exp(+i) spectrum; conjugate back to the CMSIS
            // convention, invert, and rescale to the unnormalized inverse
            // (the caller's 2/N then normalizes the round trip).
            for (int k = 3; k < m_n; k += 2) {
                a[k] = -a[k];
            }
            arm_rfft_fast_f32(&m_inst, a, m_scratch.data(), 1);
            const float s = 0.5f * static_cast<float>(m_n);
            for (int i = 0; i < m_n; ++i) {
                a[i] = m_scratch[i] * s;
            }
        }

      private:
        arm_rfft_fast_instance_f32 m_inst{};
        std::vector<float>         m_scratch;
        int                        m_n = 0;
    };

} // namespace tap::dsp::detail
