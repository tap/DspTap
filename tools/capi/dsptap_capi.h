/// @file dsptap_capi.h
/// @brief Minimal C ABI over the DspTap primitives (yin, psola, pvoc), for language bindings and
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

#ifdef __cplusplus
}
#endif
