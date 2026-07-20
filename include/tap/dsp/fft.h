/// @file fft.h
/// @brief Real FFT with a fixed numeric contract and pluggable float32 backends.
// SPDX-License-Identifier: MIT
// Copyright 2025-2026 Timothy Place and the DspTap contributors.
//
// Extracted from the Tap family DSP libraries (MuTap's adaptive-filtering FFT
// and AmbiTap's binaural convolution FFT), which each carried a byte-identical
// copy of the vendored Ooura transform under a diverging wrapper. This is the
// consolidated wrapper: one Ooura numeric contract, one place to add a faster
// backend. See README.md for the provenance and migration notes.

#pragma once

#include <cassert>
#include <cmath>
#include <cstddef>
#include <memory>
#include <type_traits>
#include <vector>

// Ooura C functions. rdft comes from third_party/ooura/fftsg.c (double);
// rdft_f is the same source instantiated for float (fftsg_float.c) so both
// precisions can link into one binary. They are built into the DspTap::fft
// static library; the tap::dsp INTERFACE target links it automatically.
// cdft/cdft_f (the complex transform) are declared for consumers that need raw
// Ooura access; the wrappers below use only rdft/rdft_f.
extern "C" {
void rdft(int n, int isgn, double* a, int* ip, double* w);
void cdft(int n, int isgn, double* a, int* ip, double* w);
void rdft_f(int n, int isgn, float* a, int* ip, float* w);
void cdft_f(int n, int isgn, float* a, int* ip, float* w);
}

// Optional per-platform float32 FFT backends, chosen by the build. AT MOST ONE
// may be defined (they are mutually exclusive), and each applies ONLY to float
// — double always stays Ooura, the golden model, so the double-precision
// reference battery is unaffected:
//   TAP_DSP_FFT_CMSIS       CMSIS-DSP Helium on the bare-metal Cortex-M55
//   TAP_DSP_FFT_ACCELERATE  Apple's vDSP (Accelerate) on macOS
// Each wrapper below re-presents its backend in Ooura's EXACT numeric contract
// (same packed layout, exp(+i) sign convention, unnormalized inverse), so every
// intermediate spectrum matches the Ooura build to float epsilon and the whole
// float32 test battery stays a valid oracle. Measured vs autovectorized Ooura:
// ~3x fewer instructions on the M55, ~3x faster on Apple Silicon per transform.
#if defined(TAP_DSP_FFT_CMSIS) && defined(TAP_DSP_FFT_ACCELERATE)
#error "TAP_DSP_FFT_CMSIS and TAP_DSP_FFT_ACCELERATE are mutually exclusive"
#endif
#if defined(TAP_DSP_FFT_CMSIS)
#include "arm_math.h"
#endif
#if defined(TAP_DSP_FFT_ACCELERATE)
#include <Accelerate/Accelerate.h>
#endif
#if defined(TAP_DSP_FFT_CMSIS) || defined(TAP_DSP_FFT_ACCELERATE)
#define TAP_DSP_FFT_FLOAT_BACKEND 1
#endif

namespace tap::dsp {

    namespace detail {
        inline void ooura_rdft(int n, int isgn, double* a, int* ip, double* w) {
            rdft(n, isgn, a, ip, w);
        }
        inline void ooura_rdft(int n, int isgn, float* a, int* ip, float* w) {
            rdft_f(n, isgn, a, ip, w);
        }

#if defined(TAP_DSP_FFT_CMSIS)
        // Wraps CMSIS-DSP's radix-4/8 Helium real FFT to reproduce Ooura's
        // exact float32 contract. CMSIS uses the engineering convention
        // exp(-i2*pi/N) and a 1/N-normalized inverse; Ooura uses exp(+i2*pi/N)
        // and an unnormalized inverse (caller applies 2/N). We reconcile by
        // conjugating the imaginary bins on every transform and scaling the
        // inverse by N/2 — both verified against Ooura to <2e-7 relative error
        // at N=512 and N=2048 (the certified geometries).
        class cmsis_real_fft_f32 {
          public:
            void init(int n) {
                m_n = n;
                m_scratch.assign(static_cast<size_t>(n), 0.0f);
                arm_rfft_fast_init_f32(&m_inst, static_cast<uint16_t>(n));
            }

            void forward_inplace(float* a) noexcept {
                arm_rfft_fast_f32(&m_inst, a, m_scratch.data(), 0);
                // Copy back, conjugating imaginary bins into Ooura convention.
                a[0] = m_scratch[0]; // DC (real)
                a[1] = m_scratch[1]; // Nyquist (real)
                for (int k = 2; k < m_n; ++k) {
                    a[k] = (k & 1) ? -m_scratch[k] : m_scratch[k];
                }
            }

            void inverse_inplace(float* a) noexcept {
                // a is an Ooura-convention packed spectrum; conjugate back to
                // CMSIS convention, invert, and rescale to Ooura's unnormalized
                // inverse (caller's 2/N then normalizes the round trip).
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

#endif // TAP_DSP_FFT_CMSIS

#if defined(TAP_DSP_FFT_ACCELERATE)
        // Wraps Apple's vDSP real FFT (Accelerate) to reproduce Ooura's exact
        // float32 contract. vDSP works in split-complex form, in the
        // engineering convention exp(-i2*pi/N), with a 2x-scaled forward; we
        // deinterleave/reinterleave (ctoz/ztoc), conjugate the imaginary bins,
        // and apply the measured scales — x0.5 on the forward, x0.25 on the
        // inverse, which lands on Ooura's UNNORMALIZED inverse (the caller's
        // 2/N then normalizes the round trip). Constants verified on Apple
        // Silicon to <4e-7 relative error at N=512 and N=2048 (bench/vdsp).
        class accelerate_real_fft_f32 {
          public:
            void init(int n) {
                m_n     = n;
                m_log2n = static_cast<int>(std::lround(std::log2(static_cast<double>(n))));
                // Shared, read-only twiddle tables: copyable value semantics
                // (basic_real_fft is held by value in the chain) with a single
                // owner-managed lifetime, and safe to share across transforms.
                m_setup = std::shared_ptr<std::remove_pointer_t<FFTSetup>>(
                    vDSP_create_fftsetup(static_cast<vDSP_Length>(m_log2n), kFFTRadix2), vDSP_destroy_fftsetup);
                m_rp.assign(static_cast<size_t>(n) / 2, 0.0f);
                m_ip.assign(static_cast<size_t>(n) / 2, 0.0f);
            }

            void forward_inplace(float* a) noexcept {
                // Index the split halves through raw pointers (pointer + int is
                // warning-free; a std::vector subscript would be int->size_t).
                float* const    rp = m_rp.data();
                float* const    ip = m_ip.data();
                DSPSplitComplex sp{rp, ip};
                vDSP_ctoz(reinterpret_cast<const DSPComplex*>(a), 2, &sp, 1, static_cast<vDSP_Length>(m_n / 2));
                vDSP_fft_zrip(m_setup.get(), &sp, 1, static_cast<vDSP_Length>(m_log2n), kFFTDirection_Forward);
                a[0] = rp[0] * 0.5f; // DC (real)
                a[1] = ip[0] * 0.5f; // Nyquist (real)
                for (int k = 1; k < m_n / 2; ++k) {
                    a[2 * k]     = rp[k] * 0.5f;
                    a[2 * k + 1] = -ip[k] * 0.5f; // conjugate into Ooura's exp(+i)
                }
            }

            void inverse_inplace(float* a) noexcept {
                // a is an Ooura-packed spectrum: rebuild vDSP's split form (undo
                // the 0.5, conjugate back to exp(-i)), invert, interleave, and
                // rescale to Ooura's unnormalized inverse.
                float* const    rp = m_rp.data();
                float* const    ip = m_ip.data();
                DSPSplitComplex sp{rp, ip};
                rp[0] = 2.0f * a[0];
                ip[0] = 2.0f * a[1];
                for (int k = 1; k < m_n / 2; ++k) {
                    rp[k] = 2.0f * a[2 * k];
                    ip[k] = -2.0f * a[2 * k + 1];
                }
                vDSP_fft_zrip(m_setup.get(), &sp, 1, static_cast<vDSP_Length>(m_log2n), kFFTDirection_Inverse);
                vDSP_ztoc(&sp, 1, reinterpret_cast<DSPComplex*>(a), 2, static_cast<vDSP_Length>(m_n / 2));
                for (int i = 0; i < m_n; ++i) {
                    a[i] *= 0.25f;
                }
            }

          private:
            std::shared_ptr<std::remove_pointer_t<FFTSetup>> m_setup;
            std::vector<float>                               m_rp, m_ip;
            int                                              m_n     = 0;
            int                                              m_log2n = 0;
        };
#endif // TAP_DSP_FFT_ACCELERATE

#if defined(TAP_DSP_FFT_FLOAT_BACKEND)
#if defined(TAP_DSP_FFT_CMSIS)
        using float_fft_engine = cmsis_real_fft_f32;
#else
        using float_fft_engine = accelerate_real_fft_f32;
#endif
        // Empty stand-in so basic_real_fft<double> carries no backend state.
        struct fft_engine_noop {
            void init(int) noexcept {}
        };
        template <typename Sample>
        using float_engine_t = std::conditional_t<std::is_same_v<Sample, float>, float_fft_engine, fft_engine_noop>;
#endif
    } // namespace detail

    /// Real FFT using the Ooura split-radix algorithm, parameterized over the
    /// sample type (float or double — the two Ooura instantiations).
    ///
    /// FFT size must be a power of 2 (>= 4), fixed at construction. Workspace
    /// (bit-reversal and trig tables) is allocated once in the constructor;
    /// the transforms themselves are noexcept and allocation-free, so they
    /// are safe on a real-time audio thread.
    ///
    /// Packing after a forward transform of N real samples (N/2 + 1 bins):
    ///   - bin[0].real    = data[0]   (DC;      imag is zero, not stored)
    ///   - bin[N/2].real  = data[1]   (Nyquist; imag is zero, not stored)
    ///   - bin[k].real    = data[2k], bin[k].imag = data[2k+1]  for 1 <= k < N/2
    ///
    /// Sign convention: the forward transform is A[k] = sum_j a[j]*W^(jk) with
    /// W = exp(+2*pi*i/N) — the imaginary parts are CONJUGATED relative to the
    /// engineering-convention DFT (exp(-2*pi*i/N)). Spectral products (e.g.
    /// fast convolution, FDAF regressor accumulation) are unaffected as long
    /// as every operand uses this class; conjugate if importing spectra
    /// computed elsewhere.
    ///
    /// The raw inverse is unnormalized: inverse_inplace() must be followed by
    /// a 2/N scaling for a round trip, which inverse() applies for you.
    template <typename Sample>
    class basic_real_fft {
        static_assert(std::is_same_v<Sample, float> || std::is_same_v<Sample, double>,
                      "basic_real_fft supports the two Ooura instantiations: float and double");

      public:
        explicit basic_real_fft(size_t size)
            : m_size(static_cast<int>(size)) {
            assert(size >= 4 && (size & (size - 1)) == 0);
#if defined(TAP_DSP_FFT_FLOAT_BACKEND)
            if constexpr (std::is_same_v<Sample, float>) {
                m_engine.init(m_size); // Ooura workspace stays unallocated
                return;
            }
#endif
            m_ip.assign(2 + static_cast<size_t>(std::sqrt(static_cast<double>(size) / 2.0)) + 1, 0);
            m_w.assign(size / 2, Sample(0));
            m_ip[0] = 0; // triggers Ooura table init on first call
        }

        size_t size() const noexcept { return static_cast<size_t>(m_size); }
        size_t num_bins() const noexcept { return static_cast<size_t>(m_size / 2 + 1); }

        /// In-place forward FFT: time-domain Sample[size] -> packed spectrum Sample[size].
        void forward_inplace(Sample* data) noexcept {
#if defined(TAP_DSP_FFT_FLOAT_BACKEND)
            if constexpr (std::is_same_v<Sample, float>) {
                m_engine.forward_inplace(data);
                return;
            }
#endif
            detail::ooura_rdft(m_size, 1, data, m_ip.data(), m_w.data());
        }

        /// In-place inverse FFT: packed spectrum Sample[size] -> time-domain Sample[size],
        /// UNSCALED — multiply by 2/size for a normalized round trip.
        void inverse_inplace(Sample* data) noexcept {
#if defined(TAP_DSP_FFT_FLOAT_BACKEND)
            if constexpr (std::is_same_v<Sample, float>) {
                m_engine.inverse_inplace(data);
                return;
            }
#endif
            detail::ooura_rdft(m_size, -1, data, m_ip.data(), m_w.data());
        }

        /// Out-of-place forward FFT. Output may alias input.
        void forward(const Sample* input, Sample* output) noexcept {
            copy(input, output);
            forward_inplace(output);
        }

        /// Out-of-place inverse FFT, scaled by 2/size so forward() -> inverse()
        /// reproduces the input. Output may alias input.
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
#if defined(TAP_DSP_FFT_FLOAT_BACKEND)
        detail::float_engine_t<Sample> m_engine;
#endif
    };

    /// Double-precision real FFT — the desktop/golden-model profile.
    using real_fft = basic_real_fft<double>;

    /// Single-precision real FFT — the embedded real-time profile (Cortex-M55,
    /// Hexagon HVX), where hardware floating point is single-precision only.
    using real_fft32 = basic_real_fft<float>;

} // namespace tap::dsp
