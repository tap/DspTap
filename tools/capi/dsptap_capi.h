/// @file dsptap_capi.h
/// @brief Minimal C ABI over the DspTap primitives (fft, yin, psola, pvoc, log_mel, decimate), for language
///        bindings and the verification notebooks (notebooks/ drive it via ctypes).
///
///        Conventions: plain C types only; the caller owns all arrays and sizes them. Handle-based
///        functions return 0 on success and -1 on any error (bad argument, bad handle). No global
///        state. Everything runs the double-precision golden profile — the notebooks verify the
///        same code the consuming libraries compile — except where a function documents a profile
///        selector (the FFT), so the notebooks can measure the embedded profiles as well.
// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.

#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#if defined(_WIN32)
#define DSPTAP_API __declspec(dllexport)
#else
#define DSPTAP_API __attribute__((visibility("default")))
#endif

/// The FFT entry points are declared noexcept: dsptap_fft_create catches and returns NULL; the
/// others call the header's noexcept transforms and plain copy loops and cannot throw (if that
/// ever changed they would std::terminate rather than return -1). A C translation unit sees
/// nothing. The older entry points predate this and are not annotated; making the whole ABI
/// uniform (typed opaque handles, exception safety, no allocation in process) is Stage 6 hygiene
/// in docs/audit-fft-and-code-smells.md.
#ifdef __cplusplus
#define DSPTAP_NOEXCEPT noexcept
#else
#define DSPTAP_NOEXCEPT
#endif

typedef void* dsptap_fft;
typedef void* dsptap_yin;
typedef void* dsptap_psola;
typedef void* dsptap_pvoc;
typedef void* dsptap_log_mel;
typedef void* dsptap_decimator;

/// -- fft ------------------------------------------------------------------------------------

/// Numeric profile of a real FFT handle: which tap::dsp::basic_real_fft<Sample> instantiation
/// it runs. DOUBLE is the golden model; FLOAT is the embedded profile (Ooura, or the platform
/// backend the build selected — see dsptap_fft_backend). Q15 and Q31 are reserved for the
/// fixed-point profiles of Stage 3 (docs/audit-fft-and-code-smells.md, Stage 3c wires them
/// here); until then dsptap_fft_create returns NULL for them.
#define DSPTAP_FFT_PROFILE_DOUBLE 0
#define DSPTAP_FFT_PROFILE_FLOAT 1
#define DSPTAP_FFT_PROFILE_Q15 2
#define DSPTAP_FFT_PROFILE_Q31 3

/// Create a real FFT of `size` points (a power of two >= 4) in the given profile, or NULL on a
/// bad size, an unavailable profile, or allocation failure. All workspace is allocated here; the
/// transforms below are allocation-free.
DSPTAP_API dsptap_fft dsptap_fft_create(int size, int profile) DSPTAP_NOEXCEPT;
DSPTAP_API void       dsptap_fft_destroy(dsptap_fft h) DSPTAP_NOEXCEPT;
DSPTAP_API int        dsptap_fft_size(dsptap_fft h) DSPTAP_NOEXCEPT;
/// size/2 + 1 (DC .. Nyquist).
DSPTAP_API int dsptap_fft_num_bins(dsptap_fft h) DSPTAP_NOEXCEPT;
/// The DSPTAP_FFT_PROFILE_* the handle was created with.
DSPTAP_API int dsptap_fft_profile(dsptap_fft h) DSPTAP_NOEXCEPT;
/// sizeof the profile's native sample (8 for DOUBLE, 4 for FLOAT) — the element width of the
/// buffers the *_inplace_raw transforms take.
DSPTAP_API int dsptap_fft_sample_bytes(dsptap_fft h) DSPTAP_NOEXCEPT;
/// The float32 engine this build compiled: "ooura", "accelerate" (Apple vDSP) or "cmsis"
/// (CMSIS-DSP Helium). The double profile is always Ooura.
DSPTAP_API const char* dsptap_fft_backend(void) DSPTAP_NOEXCEPT;

/// Forward transform of size() real samples into the packed spectrum, in double regardless of
/// profile (the FLOAT profile converts at the boundary, so `in` is rounded to float32 first and
/// the result is widened back). `out` may alias `in`. Packing and sign convention are the
/// header's contract (fft.h): out[0] = DC, out[1] = Nyquist (both real), out[2k] + i*out[2k+1]
/// = bin k for 1 <= k < size/2, with W = exp(+2*pi*i/size) — CONJUGATE to the engineering
/// convention numpy.fft.rfft uses.
///
/// Fixed point (Stage 3c, forward-looking so it extends rather than reinterprets): the Q15 and
/// Q31 profiles convert at this boundary as x / 2^15 and x / 2^31 (full scale = 1.0) on the way
/// in and the inverse on the way out, with the profile's fixed scaling policy; the block-
/// floating-point policy and its per-transform exponent arrive through a separate create
/// variant and a separate raw-path entry point, not through a change to these signatures.
DSPTAP_API int dsptap_fft_forward(dsptap_fft h, const double* in, double* out) DSPTAP_NOEXCEPT;
/// Inverse of dsptap_fft_forward, scaled by 2/size like basic_real_fft::inverse() so
/// forward -> inverse reproduces the input. `out` may alias `in`.
DSPTAP_API int dsptap_fft_inverse(dsptap_fft h, const double* in, double* out) DSPTAP_NOEXCEPT;

/// In-place transforms on the profile's NATIVE sample type: `data` points at size() samples of
/// dsptap_fft_sample_bytes() each (double for DOUBLE, float for FLOAT), in the same packing,
/// and must be aligned for that type (a double* or float* the caller obtained as such; not an
/// offset into a byte buffer). These are the header's forward_inplace()/inverse_inplace() with
/// no conversion at all, so the notebooks measure the embedded profile's own arithmetic rather
/// than a double round trip. The inverse is UNSCALED (multiply by 2/size for a round trip),
/// exactly as in fft.h.
DSPTAP_API int dsptap_fft_forward_inplace_raw(dsptap_fft h, void* data) DSPTAP_NOEXCEPT;
DSPTAP_API int dsptap_fft_inverse_inplace_raw(dsptap_fft h, void* data) DSPTAP_NOEXCEPT;

/// -- yin ------------------------------------------------------------------------------------

/// Create a detector (window, tau_min, tau_max as in tap::dsp::basic_yin), or NULL on bad geometry.
DSPTAP_API dsptap_yin dsptap_yin_create(int window, int tau_min, int tau_max);
DSPTAP_API void       dsptap_yin_destroy(dsptap_yin h);
DSPTAP_API int        dsptap_yin_set_threshold(dsptap_yin h, double threshold);
DSPTAP_API int        dsptap_yin_frame_size(dsptap_yin h);

/// Analyze one frame of frame_size() samples (oldest first). Writes the fractional period in
/// samples (0 = unvoiced) and the normalized aperiodicity. Returns 0, or -1 on a bad handle.
DSPTAP_API int dsptap_yin_analyze(dsptap_yin h, const double* frame, double* period, double* aperiodicity);

/// Track a whole signal: analyze every `hop` samples and write up to max_out periods (0 where
/// unvoiced). Returns the number of analyses written, or -1 on error.
DSPTAP_API int dsptap_yin_track(dsptap_yin h, const double* x, int n, int hop, double* periods, int max_out);

/// -- psola ----------------------------------------------------------------------------------

DSPTAP_API dsptap_psola dsptap_psola_create(int max_period);
DSPTAP_API void         dsptap_psola_destroy(dsptap_psola h);
DSPTAP_API int          dsptap_psola_latency(dsptap_psola h);
DSPTAP_API int          dsptap_psola_clear(dsptap_psola h);

/// Shift n samples at a fixed source period and ratio (state persists across calls).
DSPTAP_API int dsptap_psola_process(dsptap_psola h, const double* in, double* out, int n, double period, double ratio);

/// -- pvoc -----------------------------------------------------------------------------------

DSPTAP_API dsptap_pvoc dsptap_pvoc_create(int fft_size);
DSPTAP_API void        dsptap_pvoc_destroy(dsptap_pvoc h);
DSPTAP_API int         dsptap_pvoc_latency(dsptap_pvoc h);
DSPTAP_API int         dsptap_pvoc_set_formant(dsptap_pvoc h, int on);
DSPTAP_API int         dsptap_pvoc_clear(dsptap_pvoc h);

/// Shift n samples at a fixed ratio (state persists across calls).
DSPTAP_API int dsptap_pvoc_process(dsptap_pvoc h, const double* in, double* out, int n, double ratio);

/// -- log_mel --------------------------------------------------------------------------------

/// Create a log-mel front end at the given geometry (tap::dsp::log_mel_geometry; sqrt_window
/// selects the sqrt-Hann window, preemphasis 0 = off), or NULL on invalid geometry. The log
/// constants and PCEN default per the header; the setters below rebuild the object (geometry is
/// fixed at construction), which also resets it.
DSPTAP_API dsptap_log_mel dsptap_log_mel_create(double sample_rate, int frame, int hop, int fft_size, int bands,
                                                double fmin_hz, double fmax_hz, int sqrt_window, double preemphasis);
DSPTAP_API void           dsptap_log_mel_destroy(dsptap_log_mel h);
DSPTAP_API int            dsptap_log_mel_set_log(dsptap_log_mel h, double floor, double shift, double scale);
DSPTAP_API int dsptap_log_mel_set_pcen(dsptap_log_mel h, int enabled, double smoother, double alpha, double delta,
                                       double power, double epsilon);
DSPTAP_API int dsptap_log_mel_reset(dsptap_log_mel h);
DSPTAP_API int dsptap_log_mel_bands(dsptap_log_mel h);
DSPTAP_API int dsptap_log_mel_latency(dsptap_log_mel h);
/// The formula-level contract version of log_mel.h (trained models record it).
DSPTAP_API int dsptap_log_mel_contract_version(void);

/// Stream n samples; writes up to max_frames frames of bands() features (row-major) and returns the
/// number written, or -1 on error. Partial hops carry over to the next call.
DSPTAP_API int dsptap_log_mel_process(dsptap_log_mel h, const double* x, int n, double* features, int max_frames);

/// -- decimate -------------------------------------------------------------------------------

/// Create a decimator by ratio 2, 3 or 6 (tap::dsp::basic_decimator<float, M>; transparent selects
/// that profile over economy), or NULL on a bad ratio. The float profile IS the golden model on
/// the FIR substrate; the double arrays here are converted at the boundary.
DSPTAP_API dsptap_decimator dsptap_decimator_create(int ratio, int transparent);
DSPTAP_API void             dsptap_decimator_destroy(dsptap_decimator h);
DSPTAP_API int              dsptap_decimator_taps(dsptap_decimator h);
DSPTAP_API int              dsptap_decimator_latency(dsptap_decimator h);
DSPTAP_API int              dsptap_decimator_reset(dsptap_decimator h);
/// Outputs the next call with n inputs will produce.
DSPTAP_API int dsptap_decimator_outputs_for(dsptap_decimator h, int n);

/// Decimate n samples; writes up to max_out outputs and returns the number written (all of
/// outputs_for(n) when max_out allows; excess outputs are dropped), or -1 on error.
DSPTAP_API int dsptap_decimator_process(dsptap_decimator h, const double* in, int n, double* out, int max_out);

#ifdef __cplusplus
}
#endif
