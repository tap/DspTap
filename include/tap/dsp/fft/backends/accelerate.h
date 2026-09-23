/// @file accelerate.h
/// @brief Apple vDSP (Accelerate) float32 real FFT engine, re-presenting the split-radix contract.
// SPDX-License-Identifier: MIT
// Copyright 2025-2026 Timothy Place and the DspTap contributors.
//
// One of the two accelerated float32 engines behind tap::dsp::basic_real_fft
// (Stage 4 of docs/audit-fft-and-code-smells.md: the engine is a template
// parameter of the class; this header is included by whoever selects this
// engine — fft.h under TAP_DSP_FFT_ACCELERATE, where it is the float
// default, or a translation unit that names detail::accelerate_real_fft_f32
// explicitly, e.g. the same-binary parity test on the macOS leg). It needs
// the Accelerate framework at link time, so it is Apple-only.

#pragma once

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <new>
#include <type_traits>
#include <vector>

#include <Accelerate/Accelerate.h>

#include "tap/dsp/detail/expects.h"

namespace tap::dsp::detail {

    /// Wraps Apple's vDSP real FFT (Accelerate) to reproduce the split-radix
    /// engine's exact float32 contract (Ooura's). vDSP works in split-complex
    /// form, in the engineering convention exp(-i2*pi/N), with a 2x-scaled
    /// forward; the wrapper deinterleaves/reinterleaves (ctoz/ztoc),
    /// conjugates the imaginary bins, and applies the measured scales — x0.5
    /// on the forward, x0.25 on the inverse, which lands on the UNNORMALIZED
    /// inverse (the caller's 2/N then normalizes the round trip). Constants
    /// verified on Apple Silicon to <4e-7 relative error at N = 512 and
    /// N = 2048 (bench/vdsp; tests/test_fft_backend.cpp on the macOS leg).
    ///
    /// Engine contract numbers (what basic_real_fft reads from every engine):
    ///   - Size range k_min_size = 4 … k_max_size = 2^20, a power of two.
    ///     vDSP documents no maximum; 2^20 is the bound the float oracle
    ///     (tests/test_fft_oracle.cpp) has always stated for this engine and
    ///     is the size the split-radix gate itself runs to. The wrapper is
    ///     swept through the class by the oracle's closed forms to 65536 on
    ///     the macOS leg and pinned bin-for-bin at 512 / 2048; sizes between
    ///     65536 and 2^20 are inside the stated range on vDSP's word, not on
    ///     a measurement in this repo.
    ///   - k_is_shareable = false: the split-complex halves are per-object
    ///     scratch, so two threads may not transform through one object
    ///     concurrently; the transforms are non-const accordingly.
    ///   - Everything is allocated in the constructor; the transforms are
    ///     noexcept and allocate nothing in the wrapper (the framework is not
    ///     audited by tests/test_fft_rt.cpp's operator-new guard).
    ///
    /// REPRODUCIBILITY — vDSP dispatches on buffer alignment, and the paths
    /// it dispatches to do not agree bit-for-bit. Diagnosed in tap/MuTap#31;
    /// the buffers below are aligned to compensate, and the detail is
    /// documented at k_align_bytes.
    ///
    /// What this engine does and does not promise, stated precisely, because
    /// a compliance battery was built on the wrong reading of it:
    ///
    ///   IT DOES  return a fixed function of its input, now that the split
    ///            buffers are aligned. Verified by fft_alignment_stability,
    ///            which is the gate; fft_backend_parity cannot see this
    ///            class of bug because it runs a single process against a
    ///            single allocation.
    ///   IT DOES  agree with the split-radix engine to <4e-7 measured as
    ///            peak-normalized absolute error, which is what "relative"
    ///            meant when the bound was recorded.
    ///   IT DOES NOT agree with the split-radix engine per bin. On material
    ///            where most bins are numerically empty — an on-bin tone,
    ///            say — per-bin relative error against a double-precision
    ///            reference runs to 1e6 and beyond, for BOTH engines. The
    ///            split-radix engine is not the accurate one there; measured
    ///            on Intel it is ~4x worse than vDSP against double. Any
    ///            consumer whose behaviour depends on the contents of empty
    ///            bins is depending on rounding noise, and no engine choice
    ///            fixes that.
    class accelerate_real_fft_f32 {
      public:
        static constexpr std::size_t k_min_size     = 4;
        static constexpr std::size_t k_max_size     = std::size_t{1} << 20;
        static constexpr bool        k_is_shareable = false;

        /// @pre n is a power of two in [k_min_size, k_max_size] (TAP_EXPECTS).
        /// @throws std::bad_alloc if vDSP cannot allocate its twiddle tables.
        explicit accelerate_real_fft_f32(std::size_t n)
            : m_n(static_cast<int>(n))
            , m_log2n(static_cast<int>(std::lround(std::log2(static_cast<double>(n))))) {
            TAP_EXPECTS(n >= k_min_size && n <= k_max_size && (n & (n - 1)) == 0);
            // Shared, read-only twiddle tables: copyable value semantics
            // (basic_real_fft is held by value in the chain) with a single
            // owner-managed lifetime, and safe to share across transforms.
            //
            // vDSP_create_fftsetup returns NULL if it cannot allocate the
            // tables. Unchecked, that NULL flows straight into vDSP_fft_zrip
            // below, which is undefined behaviour — and an allocation failure
            // on a loaded machine is exactly the shape of fault that presents
            // as an intermittent one. Fail here instead, where construction
            // is already allowed to throw and the caller can fall back.
            FFTSetup setup = vDSP_create_fftsetup(static_cast<vDSP_Length>(m_log2n), kFFTRadix2);
            if (setup == nullptr) {
                throw std::bad_alloc();
            }
            m_setup = std::shared_ptr<std::remove_pointer_t<FFTSetup>>(setup, vDSP_destroy_fftsetup);
            // Over-allocated by k_align_pad so the halves can be handed to
            // vDSP 64-byte aligned regardless of where the allocator put
            // them — see rp()/ip() and the k_align comment.
            m_rp.assign(n / 2 + k_align_pad, 0.0f);
            m_ip.assign(n / 2 + k_align_pad, 0.0f);
        }

        [[nodiscard]] std::size_t size() const noexcept { return static_cast<std::size_t>(m_n); }

        void forward_inplace(float* a) noexcept {
            // Index the split halves through raw pointers (pointer + int is
            // warning-free; a std::vector subscript would be int->size_t).
            float* const    rp = this->rp();
            float* const    ip = this->ip();
            DSPSplitComplex sp{rp, ip};
            vDSP_ctoz(reinterpret_cast<const DSPComplex*>(a), 2, &sp, 1, static_cast<vDSP_Length>(m_n / 2));
            vDSP_fft_zrip(m_setup.get(), &sp, 1, static_cast<vDSP_Length>(m_log2n), kFFTDirection_Forward);
            a[0] = rp[0] * 0.5f; // DC (real)
            a[1] = ip[0] * 0.5f; // Nyquist (real)
            for (int k = 1; k < m_n / 2; ++k) {
                a[2 * k]     = rp[k] * 0.5f;
                a[2 * k + 1] = -ip[k] * 0.5f; // conjugate into exp(+i)
            }
        }

        void inverse_inplace(float* a) noexcept {
            // a is a packed exp(+i) spectrum: rebuild vDSP's split form (undo
            // the 0.5, conjugate back to exp(-i)), invert, interleave, and
            // rescale to the unnormalized inverse.
            float* const    rp = this->rp();
            float* const    ip = this->ip();
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
        // vDSP DISPATCHES ON BUFFER ALIGNMENT, and the two paths do not
        // agree bit-for-bit. Measured on Apple M1 / macOS 26.5.2 /
        // Xcode 26.6 at N=2048 and N=4096: the same input through the same
        // FFTSetup produces one output when the split-complex halves are
        // 64-byte aligned and a different one when they are not, with the
        // two differing by ~2e-7 peak-normalized. std::vector's allocator
        // hands out 16-byte-aligned storage, so which path a process took
        // was decided by wherever the heap happened to land — the
        // per-process nondeterminism in tap/MuTap#31 (140/60 over 200
        // processes at N=2048, 110/90 at N=4096).
        //
        // Placing the halves at a fixed offset removes that variable. WHICH
        // offset matters as much as fixing it — the two kernels differ in
        // accuracy, not merely in rounding:
        //
        // Two kernels exist and they are not equally good. Measured on
        // Apple M1 / macOS 26.5.2 / Xcode 26.6, median per-bin relative
        // error against a double-precision reference at N=2048:
        //
        //   material         64-byte aligned   NOT 64-byte aligned   split-radix
        //   broadband        1.6e-07           1.2e-07               1.2e-07
        //   tone, off-bin    1.3e-06           1.3e-06               6.4e-07
        //   tone, ON-bin     0.65              1.2e-07               1.1e-07
        //
        // On material whose spectrum has exactly-empty bins — an on-bin
        // tone is the clean case — the 64-byte-aligned kernel puts the
        // MEDIAN bin 65% away from truth, while the other kernel and the
        // split-radix engine both track it to float epsilon. Any leakage
        // that lifts those bins above the noise floor hides the difference,
        // which is why a peak-normalized parity check cannot see it.
        //
        // So we pin the alignment to a fixed value that is 32-byte aligned
        // (Apple's vDSP.h asks for "preferably 16-byte aligned or better",
        // and unaligned vector loads would cost performance) but NOT
        // 64-byte aligned, which selects the accurate kernel. Fixed, so the
        // transform is a function of its input; skewed, so it is the better
        // of the two functions available.
        //
        // This steers around an UNDOCUMENTED dispatch rule. Apple's vDSP.h
        // says the routines are "free to rearrange calculations for better
        // performance" and are "not expected to conform to IEEE 754", so
        // nothing stops a future SDK from dispatching differently. If that
        // happens this stops helping — silently, which is the real risk —
        // so fft_tonal_accuracy asserts the property directly rather than
        // trusting the skew. A consumer who cannot accept that exposure at
        // all should build with TAP_DSP_FFT_ACCELERATE=OFF; MuTap's
        // compliance battery does exactly that.
        //
        // The placement is computed at use (rp()/ip()) rather than cached,
        // because basic_real_fft is copyable and held by value: a copy's
        // storage lands wherever the allocator puts it, so a cached pointer
        // would silently lose the placement the copy still needs.
        static constexpr std::size_t k_align_bytes = 64; ///< boundary vDSP dispatches on
        static constexpr std::size_t k_skew_bytes  = 32; ///< offset past it: 32-byte aligned, not 64
        static constexpr std::size_t k_align_pad   = (k_align_bytes + k_skew_bytes) / sizeof(float) + 1;

        static float* place(float* p) noexcept {
            const auto addr = reinterpret_cast<std::uintptr_t>(p);
            const auto up   = (addr + (k_align_bytes - 1)) & ~static_cast<std::uintptr_t>(k_align_bytes - 1);
            return reinterpret_cast<float*>(up + k_skew_bytes);
        }

        float* rp() noexcept { return place(m_rp.data()); }
        float* ip() noexcept { return place(m_ip.data()); }

        std::shared_ptr<std::remove_pointer_t<FFTSetup>> m_setup;
        std::vector<float>                               m_rp, m_ip;
        int                                              m_n     = 0;
        int                                              m_log2n = 0;
    };

} // namespace tap::dsp::detail
