/// @file fft.h
/// @brief Real FFT with a fixed numeric contract: double, float, Q15 and Q31 profiles.
// SPDX-License-Identifier: MIT
// Copyright 2025-2026 Timothy Place and the DspTap contributors.
//
// Extracted from the Tap family DSP libraries (MuTap's adaptive-filtering FFT
// and AmbiTap's binaural convolution FFT), which each carried a byte-identical
// copy of the vendored Ooura transform under a diverging wrapper. This is the
// consolidated wrapper: one numeric contract (Ooura's, carried since Stage 2b
// by the bit-identical C++20 port), one place to add a faster backend. See
// README.md for the provenance and migration notes.
//
// Four profiles share the contract (packing, exp(+i), unnormalized inverse):
// double (the golden model) and float (the embedded floating profile) run the
// split-radix engine in fft/split_radix.h, the C++20 transliteration of
// Ooura's rdft that Stage 2a landed bit-identical to the vendored C and
// Stage 2b routed here (docs/audit-fft-and-code-smells.md, Part 3); the
// float profile may instead be routed to an accelerated backend by the build
// (below). std::int16_t (Q15) and std::int32_t (Q31) run the int32
// fixed-point kernel in fft/fixed_point.h, whose transforms return an
// exponent in place of a floating scale (Stage 3b).

#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <new>
#include <type_traits>
#include <vector>

#include "tap/dsp/fft/fixed_point.h"
#include "tap/dsp/fft/split_radix.h"

// Header-only. No vendored C is compiled into what ships: the split-radix
// engine replaced Ooura's rdft at Stage 2b, and Stage 2c moved the C out of
// the shipping tree (docs/audit-fft-and-code-smells.md, Part 3; Decision
// D6) together with the extern "C" rdft / rdft_f declarations this header
// carried for it. The reference copy lives under tests/reference/ooura/ and
// is compiled only by tests/test_fft_parity_ooura.cpp, the bit-identity gate
// for the engine (which declares it through tests/reference/ooura_rdft.h).
// tap::dsp is a pure INTERFACE target unless TAP_DSP_FFT_CMSIS is on, in
// which case it links the CMSIS-DSP objects (root CMakeLists.txt).

// Optional per-platform float32 FFT backends, chosen by the build. AT MOST ONE
// may be defined (they are mutually exclusive), and each applies ONLY to float
// — double always runs the split-radix engine, the golden model, so the
// double-precision reference battery is unaffected:
//   TAP_DSP_FFT_CMSIS       CMSIS-DSP Helium on the bare-metal Cortex-M55
//   TAP_DSP_FFT_ACCELERATE  Apple's vDSP (Accelerate) on macOS
// Each wrapper below re-presents its backend in the split-radix engine's
// EXACT numeric contract (Ooura's: same packed layout, exp(+i) sign
// convention, unnormalized inverse), so every intermediate spectrum matches
// the default build to float epsilon and the whole float32 test battery stays
// a valid oracle. Against the vendored C, to which the engine is
// bit-identical: on the M55 icount key the C executes 1.81x (N = 512) /
// 1.97x (N = 2048) the instructions of the CMSIS build over the whole
// ratchet scenario (bench/README.md, run 35844483811); the "~3x faster on
// Apple Silicon" figure is MuTap's transform-only measurement on the C
// (tap/MuTap#31), not re-measured here (Stage 4, docs/fft-design.md).
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
        //
        // REPRODUCIBILITY — vDSP dispatches on buffer alignment, and the paths
        // it dispatches to do not agree bit-for-bit. Diagnosed in
        // tap/MuTap#31; the buffers below are aligned to compensate, and the
        // detail is documented at k_align_bytes.
        //
        // What this backend does and does not promise, stated precisely,
        // because a compliance battery was built on the wrong reading of it:
        //
        //   IT DOES  return a fixed function of its input, now that the split
        //            buffers are aligned. Verified by fft_alignment_stability,
        //            which is the gate; fft_backend_parity cannot see this
        //            class of bug because it runs a single process against a
        //            single allocation.
        //   IT DOES  agree with Ooura to <4e-7 measured as peak-normalized
        //            absolute error, which is what "relative" meant when the
        //            bound was recorded.
        //   IT DOES NOT agree with Ooura per bin. On material where most bins
        //            are numerically empty — an on-bin tone, say — per-bin
        //            relative error against a double-precision reference runs
        //            to 1e6 and beyond, for BOTH backends. Ooura is not the
        //            accurate one there; measured on Intel it is ~4x worse
        //            than vDSP against double. Any consumer whose behaviour
        //            depends on the contents of empty bins is depending on
        //            rounding noise, and no backend choice fixes that.
        class accelerate_real_fft_f32 {
          public:
            void init(int n) {
                m_n     = n;
                m_log2n = static_cast<int>(std::lround(std::log2(static_cast<double>(n))));
                // Shared, read-only twiddle tables: copyable value semantics
                // (basic_real_fft is held by value in the chain) with a single
                // owner-managed lifetime, and safe to share across transforms.
                //
                // vDSP_create_fftsetup returns NULL if it cannot allocate the
                // tables. Unchecked, that NULL flows straight into
                // vDSP_fft_zrip below, which is undefined behaviour — and an
                // allocation failure on a loaded machine is exactly the shape
                // of fault that presents as an intermittent one. Fail here
                // instead, where construction is already allowed to throw and
                // the caller can fall back.
                FFTSetup setup = vDSP_create_fftsetup(static_cast<vDSP_Length>(m_log2n), kFFTRadix2);
                if (setup == nullptr) {
                    throw std::bad_alloc();
                }
                m_setup = std::shared_ptr<std::remove_pointer_t<FFTSetup>>(setup, vDSP_destroy_fftsetup);
                // Over-allocated by k_align_pad so the halves can be handed to
                // vDSP 64-byte aligned regardless of where the allocator put
                // them — see rp()/ip() and the k_align comment.
                m_rp.assign(static_cast<size_t>(n) / 2 + k_align_pad, 0.0f);
                m_ip.assign(static_cast<size_t>(n) / 2 + k_align_pad, 0.0f);
            }

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
                    a[2 * k + 1] = -ip[k] * 0.5f; // conjugate into Ooura's exp(+i)
                }
            }

            void inverse_inplace(float* a) noexcept {
                // a is an Ooura-packed spectrum: rebuild vDSP's split form (undo
                // the 0.5, conjugate back to exp(-i)), invert, interleave, and
                // rescale to Ooura's unnormalized inverse.
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
            // accuracy, not merely in rounding. See k_skew_bytes.
            //
            // Computed at use rather than cached, because basic_real_fft is
            // copyable and held by value: a copy's storage lands wherever the
            // allocator puts it, so a cached pointer would silently lose the
            // alignment the copy still needs.
            // WHICH alignment, and why it is deliberately not the "best" one.
            //
            // Two kernels exist and they are not equally good. Measured on
            // Apple M1 / macOS 26.5.2 / Xcode 26.6, median per-bin relative
            // error against a double-precision reference at N=2048:
            //
            //   material         64-byte aligned   NOT 64-byte aligned   Ooura
            //   broadband        1.6e-07           1.2e-07               1.2e-07
            //   tone, off-bin    1.3e-06           1.3e-06               6.4e-07
            //   tone, ON-bin     0.65              1.2e-07               1.1e-07
            //
            // On material whose spectrum has exactly-empty bins — an on-bin
            // tone is the clean case — the 64-byte-aligned kernel puts the
            // MEDIAN bin 65% away from truth, while the other kernel and Ooura
            // both track it to float epsilon. Any leakage that lifts those bins
            // above the noise floor hides the difference, which is why a
            // peak-normalized parity check cannot see it.
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
            // Computed at use rather than cached, because basic_real_fft is
            // copyable and held by value: a copy's storage lands wherever the
            // allocator puts it, so a cached pointer would silently lose the
            // placement the copy still needs.
            static constexpr size_t k_align_bytes = 64; ///< boundary vDSP dispatches on
            static constexpr size_t k_skew_bytes  = 32; ///< offset past it: 32-byte aligned, not 64
            static constexpr size_t k_align_pad   = (k_align_bytes + k_skew_bytes) / sizeof(float) + 1;

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
#endif // TAP_DSP_FFT_ACCELERATE

        // The engine behind basic_real_fft's floating profiles: the split-radix
        // engine for double always, and for float unless the build selected an
        // accelerated backend above. The two shapes differ only in how they are
        // constructed (the engine takes its size in the constructor; the
        // backend wrappers are default-constructed and init()ed, unchanged
        // from before the Stage 2b flip), which make_floating_engine hides.
        template <typename Sample>
        struct floating_engine {
            using type = split_radix_rdft<Sample>;
        };
#if defined(TAP_DSP_FFT_FLOAT_BACKEND)
        template <>
        struct floating_engine<float> {
#if defined(TAP_DSP_FFT_CMSIS)
            using type = cmsis_real_fft_f32;
#else
            using type = accelerate_real_fft_f32;
#endif
        };
#endif
        template <typename Sample>
        using floating_engine_t = typename floating_engine<Sample>::type;

        template <typename Sample>
        floating_engine_t<Sample> make_floating_engine(std::size_t n) {
            if constexpr (std::is_same_v<floating_engine_t<Sample>, split_radix_rdft<Sample>>) {
                return split_radix_rdft<Sample>(n);
            }
            else {
                floating_engine_t<Sample> engine;
                engine.init(static_cast<int>(n));
                return engine;
            }
        }
    } // namespace detail

    /// Real FFT with a fixed numeric contract, parameterized over the sample
    /// type: float or double (the two floating profiles, this primary
    /// template) and std::int16_t or std::int32_t (the Q15 and Q31
    /// fixed-point profiles, the specialization below, which adds the Scaling
    /// policy parameter and returns an exponent from every transform).
    /// Scaling is meaningful for the fixed-point profiles only; the floating
    /// profiles accept scaling::fixed (the default) and nothing else.
    ///
    /// Routing (Stage 2b of docs/audit-fft-and-code-smells.md, Part 3):
    ///   - double  -> detail::split_radix_rdft<double> (fft/split_radix.h),
    ///                always. The golden model.
    ///   - float   -> detail::split_radix_rdft<float>, unless the build
    ///                defines TAP_DSP_FFT_CMSIS (CMSIS-DSP Helium, bare-metal
    ///                Cortex-M55) or TAP_DSP_FFT_ACCELERATE (Apple vDSP), in
    ///                which case the backend wrapper above re-presents the
    ///                same contract to float epsilon (tests/test_fft_backend.cpp).
    ///   - Q15/Q31 -> detail::fixed_point_rdft (the specialization below).
    /// The split-radix engine is the C++20 transliteration of Ooura's rdft
    /// and is BIT-IDENTICAL to the C it replaced for both precisions (the
    /// reference copy under tests/reference/ooura/, compiled by the gate
    /// tests/test_fft_parity_ooura.cpp alone since Stage 2c; Decision D10),
    /// so the flip changed no
    /// output bit of any consumer built at default fp-contraction (measured
    /// on every CI platform) or with clang and FMA; the one measured
    /// exception is a g++ x86-64 build with -march (FMA), which moves the
    /// float profile's last bits at N >= 1024 (the fp-contraction policy
    /// below, D9). tests/test_fft_routing.cpp pins that this class's output
    /// is byte-identical to the engine's.
    ///
    /// FFT size must be a power of 2 (>= 4), fixed at construction. Workspace
    /// (bit-reversal and trig tables) is allocated AND BUILT in the constructor
    /// (the vendored C built its tables lazily on the first transform; the
    /// engine does not, audit item F6), so the first transform costs what
    /// every later one costs; the transforms themselves are noexcept and
    /// allocation-free, so they are safe on a real-time audio thread
    /// (tests/test_fft_rt.cpp). No alignment requirement on the data pointer.
    /// NaN propagates to every bin (no data-dependent branches). Latency 0.
    /// Copyable, and a copy is bit-identical to its source; the object is held
    /// by value in every consumer.
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
    ///
    /// Noise floor (docs/fft-design.md, contract table; measured 2026-09-23
    /// on x86-64 Linux, GCC 13.3.0 and Clang 18.1.3 -O3, glibc 2.39, the
    /// engine as routed here): double against the compensated-summation DFT
    /// oracle, relative 2-norm error of the forward on full-scale uniform
    /// noise at N = 256, 1.85e-16, pinned at 4x (7.4e-16) for the libm and
    /// fp-contraction spread across hosts
    /// (`fft_oracle_floor.DoubleForwardTracksCompensatedDft`); float against
    /// double, same metric at N = 512 on the split-radix engine itself,
    /// 1.12e-7 (the audit's Part 6 N2 probe value of 1.105e-7, re-measured on
    /// what ships; the two runs differ in material, not in engine), pinned at
    /// 2x (2.25e-7) (`RealFftCrossPrecision.FloatEngineTracksDoubleAtN512`),
    /// and < 1e-6 at N = 1024 through basic_real_fft (`FloatTracksDouble`).
    ///
    /// fp-contraction policy (Decision D9; measurements in docs/fft-design.md,
    /// "The fp-contraction policy"). tap::dsp does NOT export an
    /// -ffp-contract setting to its consumers and this header sets none: the
    /// engine's arithmetic is compiled under whatever contraction the
    /// consumer's compiler applies, alike for every instantiation in that
    /// build. The bit-identity gate against the reference C runs at
    /// -ffp-contract=off on both sides, and at default flags the identity was
    /// MEASURED to hold as well on every CI platform (0 ulp on linux, windows
    /// and macOS arm64 and on the Cortex-M4, M4F, M33 and M55 legs, three of
    /// them VFMA; identical bench checksums between the C and the engine on
    /// every Ooura bench key; MuTap's fingerprint rows byte-identical through
    /// the flip, float rows included, with no flag set anywhere), because the
    /// engine's statements are textually the C's and each compiler fuses both
    /// sides alike. That is an observation, not a guarantee: g++ on x86-64
    /// built with -march (FMA) is measured to fuse the two sides differently
    /// (float moves by a few ulp at N >= 1024, double does not; clang with
    /// the same flags is identical) and is not claimed. Why no export: an
    /// INTERFACE -ffp-contract=off would reach
    /// every consumer translation unit that includes this header and would
    /// pessimize the VFMA / FMA targets the float profile exists for (the
    /// M55, Apple arm64) for the whole of that code, in exchange for a
    /// cross-compiler bit reproducibility that libm's last-bit cos/sin
    /// differences between glibc, newlib, UCRT and Apple already deny across
    /// hosts. A consumer that needs bit reproducibility across its own
    /// compilers sets -ffp-contract=off on its own targets; the contract this
    /// header makes is the one above.
    ///
    /// Thread rule: one transform at a time per object, as before the flip.
    /// The split-radix engine's transforms are const and it has no mutable
    /// state after construction, so a double object (or a float one with no
    /// backend) could be shared; this class's transforms stay non-const
    /// because the CMSIS and vDSP wrappers carry scratch, and Stage 4 states
    /// shareability as an engine trait rather than per build define.
    template <typename Sample, typename Scaling = scaling::fixed>
    class basic_real_fft {
        static_assert(std::is_same_v<Sample, float> || std::is_same_v<Sample, double>,
                      "basic_real_fft supports float and double (split radix) and std::int16_t / std::int32_t "
                      "(Q15 / Q31)");
        static_assert(std::is_same_v<Scaling, scaling::fixed>, "scaling policies apply to the fixed-point profiles");

      public:
        /// The engine this instantiation routes to (see the class docstring).
        using engine = detail::floating_engine_t<Sample>;

        /// @pre size is a power of two, size >= 4 (asserted).
        explicit basic_real_fft(size_t size)
            : m_size(static_cast<int>(size))
            , m_engine(detail::make_floating_engine<Sample>(size)) {
            assert(size >= 4 && (size & (size - 1)) == 0);
        }

        size_t size() const noexcept { return static_cast<size_t>(m_size); }
        size_t num_bins() const noexcept { return static_cast<size_t>(m_size / 2 + 1); }

        /// In-place forward FFT: time-domain Sample[size] -> packed spectrum Sample[size].
        void forward_inplace(Sample* data) noexcept { m_engine.forward_inplace(data); }

        /// In-place inverse FFT: packed spectrum Sample[size] -> time-domain Sample[size],
        /// UNSCALED — multiply by 2/size for a normalized round trip.
        void inverse_inplace(Sample* data) noexcept { m_engine.inverse_inplace(data); }

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

        /// Float-I/O convenience on the DOUBLE engine: run the double-precision
        /// transform over float buffers (copy in, transform, copy out), so a
        /// caller holding float data can use the double FFT without maintaining
        /// its own double staging buffer — e.g. AmbiTap's binaural HRTF analysis,
        /// which wants double-precision spectra from float impulse responses.
        /// Only the double instantiation offers these; basic_real_fft<float>
        /// already takes float in the same-type forward()/inverse() above. These
        /// allocate a staging buffer (a setup-time path), unlike the noexcept,
        /// allocation-free in-place transforms.
        ///
        /// DEPRECATED (Decision D5, audit item F7): removed after one consumer
        /// cycle, i.e. once MuTap and MuTap-Max have pinned a tree containing
        /// this deprecation. No consumer on disk calls them (grepped at Stage
        /// 2b: DspTap's tools/capi, MuTap, MuTap-Max); the replacement is the
        /// caller's own double staging buffer and the same-type transforms.
        [[deprecated("basic_real_fft<double>::forward(const float*, float*) allocates per call and is removed "
                     "after one consumer cycle (Decision D5); stage through a double buffer instead")]]
        void forward(const float* input, float* output)
            requires std::is_same_v<Sample, double>
        {
            std::vector<double> buf(static_cast<size_t>(m_size));
            for (int i = 0; i < m_size; ++i) {
                buf[static_cast<size_t>(i)] = static_cast<double>(input[i]);
            }
            forward_inplace(buf.data());
            for (int i = 0; i < m_size; ++i) {
                output[i] = static_cast<float>(buf[static_cast<size_t>(i)]);
            }
        }

        /// Float-I/O inverse on the double engine, scaled by 2/size like inverse().
        /// DEPRECATED with forward(const float*, float*) above (Decision D5).
        [[deprecated("basic_real_fft<double>::inverse(const float*, float*) allocates per call and is removed "
                     "after one consumer cycle (Decision D5); stage through a double buffer instead")]]
        void inverse(const float* input, float* output)
            requires std::is_same_v<Sample, double>
        {
            std::vector<double> buf(static_cast<size_t>(m_size));
            for (int i = 0; i < m_size; ++i) {
                buf[static_cast<size_t>(i)] = static_cast<double>(input[i]);
            }
            inverse_inplace(buf.data());
            const double scale = 2.0 / static_cast<double>(m_size);
            for (int i = 0; i < m_size; ++i) {
                output[i] = static_cast<float>(buf[static_cast<size_t>(i)] * scale);
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

        int    m_size;
        engine m_engine;
    };

    /// The fixed-point profiles: Q15 (std::int16_t) and Q31 (std::int32_t)
    /// I/O over one int32 radix-4 kernel (fft/fixed_point.h), same packing and
    /// exp(+i) convention as the floating profiles, with the scale carried by
    /// an exponent instead of a floating mantissa. Contract, as numbers (each
    /// pinned by the named test of the Stage 3b battery, tests/test_fft_fixed.cpp
    /// and the widened tests/test_fft.cpp):
    ///
    ///  - Exponent. Every transform returns e. Read the buffer as fractions of
    ///    full scale (Q0.15 / Q0.31) and let G be basic_real_fft<double> on the
    ///    same input read the same way: after forward_inplace, G's forward
    ///    result == data * 2^e (same packing); after inverse_inplace, G's
    ///    UNNORMALIZED inverse result == data * 2^e. The fixed-point inverse()
    ///    applies NO 2/N (unlike the floating profiles); a round trip
    ///    reconstructs x == out * 2^(e_fwd + e_inv + 1 - log2 N) under both
    ///    policies. `RoundTripReproducesInput`, `RoundTripReconstructsInputPerPolicy`,
    ///    `FixedForwardScaleIsExactlyXOverN`, `FixedForwardScaleIsExactlyXOverTwoN`.
    ///  - scaling::fixed (default): e == fixed_scaling_exponent(N)
    ///    == log2 N + fft_arith<Sample>::k_fixed_scaling_input_pre_shift in
    ///    both directions: log2 N for Q15 (output exactly X / N), log2 N + 1
    ///    for Q31 (X / 2N; the one-bit input pre-shift is the price of no
    ///    guard bits, fft_arith.h). `FixedExponentIsTheStatedConstant`,
    ///    `FixedInverseCarriesTheSameExponent`.
    ///  - scaling::block_floating: 0 <= e <= fixed_scaling_exponent(N),
    ///    data-dependent; the kernel never shifts more than the stage's
    ///    growth requires beyond the headroom present, and never consumes
    ///    the last bit. Brought to a common exponent, the block-floating
    ///    output agrees with the fixed one within a pinned bound (Q15 1.0
    ///    LSB, Q31 4.0 LSB measured; pinned 2 / 8). At the full exponent it
    ///    is bit-identical to the fixed output only when every stage shifted
    ///    the fixed amount (the full-scale patterns, pinned); in general the
    ///    two schedules reach the same exponent through different rounding
    ///    histories (an input with headroom at entry shifts less at the
    ///    first stage and catches up later) and differ by up to 4 LSB (Q31,
    ///    at the DC index; 2 LSB elsewhere) or 1 LSB (Q15, where the int32
    ///    history's few LSB32 survive the narrow only on a rounding tie),
    ///    measured and pinned at 2x. `BfpExponentIsWithinRange`,
    ///    `BfpMatchesFixedAfterShift`, `BfpAtTheFullExponentIsBitIdenticalToFixed`,
    ///    `BfpAtTheFullExponentAgreesWithFixedWithinPin`, `SilenceIsSilence`.
    ///  - Saturation. The int32 kernel performs no saturating operation for
    ///    any input under either policy: Q15 through the two guard bits of
    ///    the widened Q2.29 data, Q31 through the pre-shift (fixed) or the
    ///    headroom rule (block floating); the worst case is the packed pair
    ///    at full scale under a 45-degree twiddle. The Q15 output narrowing
    ///    (round-half-up, saturating) clamps at the rail exactly when the
    ///    true value is the rail: a full-scale Nyquist alternation lands on
    ///    32767.5 LSB and comes back as 32767, the pinned 0.5 LSB rail
    ///    shortfall. Largest deviation from the golden model over the
    ///    adversarial full-scale sweep, both directions: Q15 0.50 / 0.75 LSB
    ///    (fixed / block floating), Q31 4.25 / 15.99 LSB, at index 0 for
    ///    Q31; pinned at 1.0 / 1.5 / 8.5 / 32. `SaturationFreeWorstCaseDoesNotWrap`.
    ///  - Host identity. The fixed-point output is integer arithmetic over a
    ///    checksum-pinned table: for a fixed input it is one bit pattern on
    ///    every host, pinned per profile, policy, direction and N = 512 /
    ///    2048. `OutputFingerprintIsPinned`, `TwiddleTableChecksumIsPinned`.
    ///  - Noise floor (output-referred, against the double golden model on
    ///    the same quantized input; N = 256 / 512 / 2048, 0 to -60 dBFS; the
    ///    numbers are the `[ floor ]` rows `NoiseFloorTracksWelchModel`
    ///    prints, forward, white noise): Q15 fixed 0.26 - 0.30 LSB rms at
    ///    every N and level (the narrow's own rounding; per-bin SNR 69.1 /
    ///    66.5 / 60.2 dB at 0 dBFS); Q31 fixed 0.67 - 0.75 LSB rms (152.0 /
    ///    149.4 / 142.8 dB), level-independent, i.e. SNR falls 20 dB per
    ///    20 dB of level; block floating point keeps 84 - 91 dB (Q15) and
    ///    157 - 161 dB (Q31) at 0 dBFS and does not lose the low-level
    ///    signal (78 - 86 dB Q15, 154 - 156 dB Q31 at -40 dBFS). Welch's
    ///    variance model predicts 0.57 - 0.59 LSB32 for the kernel; the
    ///    measured/model ratio of 1.31 - 1.74 (Q31 fixed) is the
    ///    round-half-up bias, largest at index 0 under block floating point
    ///    (fft/fixed_point.h, "Honest limit"). `NoiseFloorTracksWelchModel`,
    ///    `RoundingBiasOnNegatedInputIsBounded`, `Q15TracksDouble`, `Q31TracksDouble`,
    ///    `Q15AndQ31AgreeToTheQ15Floor`.
    ///  - CMSIS-DSP compatibility (Decision D3): CMSIS documents its q15 /
    ///    q31 real FFTs as "downscaled by 2 for every stage", i.e. a forward
    ///    output of X / N with log2 N bits to upscale. The Q15 fixed forward
    ///    here is that same X / N; the Q31 fixed forward is X / 2N, one bit
    ///    below, because of the pre-shift; the inverse here is Ooura's
    ///    unnormalized inverse over 2^e, i.e. (1/N) sum X W^-jk divided by 2
    ///    (Q15) or 4 (Q31), which is not CMSIS's inverse table; block
    ///    floating point has no CMSIS analogue. The exponent, not a fixed
    ///    Q format per N, is the contract here; nothing of CMSIS's
    ///    behaviour was measured or reproduced, only its documentation read.
    ///  - Size 4 <= N <= 65536, a power of two, fixed at construction. Q15
    ///    allocates an int32 work buffer of N at construction (in-place API
    ///    preserved at the caller's int16 buffer); Q31 transforms in place.
    ///    Transforms are noexcept and allocation-free, the object is copyable,
    ///    there is no alignment requirement. `TransformsAreNoexcept`,
    ///    `ForwardInplaceAllocatesNothing` (and the inverse and out-of-place
    ///    forms), `CopyProducesBitIdenticalOutput`, `OutOfPlaceIsCopyThenInPlace`.
    ///  - Thread rule: one transform at a time per object. This deviates
    ///    from the design record (audit Part 7 lists `is_shareable` true for
    ///    both fixed profiles): Q15 is NOT shareable, since the in-place
    ///    int16 API needs the per-object int32 work buffer; Q31 has no
    ///    mutable state (its transforms touch the caller's buffer and the
    ///    tables only) but its transforms are non-const in this stage, as
    ///    the floating profiles' are. Stage 4 states shareability as an
    ///    engine trait (Q31 true, Q15 false) and decides transform constness
    ///    across engines (the split-radix port's transforms are const).
    ///  - Every transform returns the exponent and is [[nodiscard]]: under
    ///    block floating point a discarded e is a silent scale error.
    ///  - Latency 0; no NaN or denormal behaviour to state (integer data).
    ///
    /// Per-profile noise floors and the Welch-model derivation are in
    /// docs/fft-design.md ("The fixed-point profiles") and the README
    /// profiles table.
    template <typename Sample, typename Scaling>
        requires fft_arith<Sample>::k_is_fixed_point
    class basic_real_fft<Sample, Scaling> {
      public:
        using engine = detail::fixed_point_rdft<Sample, Scaling>;

        /// The constant exponent of scaling::fixed (and the upper bound of
        /// scaling::block_floating) for size n; the tests use it in place of
        /// magic numbers.
        [[nodiscard]] static constexpr int fixed_scaling_exponent(std::size_t n) noexcept {
            return engine::fixed_scaling_exponent(n);
        }

        /// @pre size is a power of two in [4, 65536] (asserted).
        explicit basic_real_fft(std::size_t size)
            : m_engine(size) {}

        [[nodiscard]] std::size_t size() const noexcept { return m_engine.size(); }
        [[nodiscard]] std::size_t num_bins() const noexcept { return m_engine.size() / 2 + 1; }

        /// In-place forward FFT: Sample[size] -> packed spectrum Sample[size],
        /// scaled by 2^-e. @return e, the output's scale: under block floating
        /// point a discarded exponent is a silent scale error, hence nodiscard.
        [[nodiscard]] int forward_inplace(Sample* data) noexcept { return m_engine.forward_inplace(data); }

        /// In-place inverse FFT: packed spectrum -> Sample[size], the
        /// UNNORMALIZED inverse scaled by 2^-e. @return e.
        [[nodiscard]] int inverse_inplace(Sample* data) noexcept { return m_engine.inverse_inplace(data); }

        /// Out-of-place forward FFT: copy, then forward_inplace. Output may
        /// alias input. @return e.
        [[nodiscard]] int forward(const Sample* input, Sample* output) noexcept {
            copy(input, output);
            return forward_inplace(output);
        }

        /// Out-of-place inverse FFT: copy, then inverse_inplace. NO 2/N is
        /// applied (the exponent carries the scale, unlike the floating
        /// profiles' inverse()). Output may alias input. @return e.
        [[nodiscard]] int inverse(const Sample* input, Sample* output) noexcept {
            copy(input, output);
            return inverse_inplace(output);
        }

      private:
        void copy(const Sample* input, Sample* output) noexcept {
            if (input != output) {
                std::copy_n(input, m_engine.size(), output);
            }
        }

        engine m_engine;
    };

    /// Double-precision real FFT — the desktop/golden-model profile.
    using real_fft = basic_real_fft<double>;

    /// Single-precision real FFT — the embedded real-time profile (Cortex-M55,
    /// Hexagon HVX), where hardware floating point is single-precision only.
    using real_fft32 = basic_real_fft<float>;

    /// Q15 real FFT (Q0.15 I/O), fixed scaling: e == log2 N, output exactly X / N.
    using real_fft_q15 = basic_real_fft<std::int16_t>;
    /// Q31 real FFT (Q0.31 I/O), fixed scaling: e == log2 N + 1, output exactly X / 2N.
    using real_fft_q31 = basic_real_fft<std::int32_t>;
    /// Q15 real FFT, block floating point: 0 <= e <= log2 N, data-dependent.
    using real_fft_q15_bfp = basic_real_fft<std::int16_t, scaling::block_floating>;
    /// Q31 real FFT, block floating point: 0 <= e <= log2 N + 1, data-dependent.
    using real_fft_q31_bfp = basic_real_fft<std::int32_t, scaling::block_floating>;

} // namespace tap::dsp
