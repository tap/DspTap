/// @file dsptap_capi.h
/// @brief Minimal C ABI over the DspTap primitives (yin, psola, pvoc, log_mel, decimate), for language bindings and
///        the verification notebooks (notebooks/ drive it via ctypes).
///
///        Conventions: plain C types only; the caller owns all arrays and sizes them. Handle-based
///        functions return 0 on success and -1 on any error (bad argument, bad handle). No global
///        state. Everything runs the double-precision golden profile — the notebooks verify the
///        same code the consuming libraries compile.
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

typedef void* dsptap_yin;
typedef void* dsptap_psola;
typedef void* dsptap_pvoc;
typedef void* dsptap_log_mel;
typedef void* dsptap_decimator;

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
