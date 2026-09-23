/// @file test_capi.c
/// @brief Smoke test of the DspTap C ABI, compiled as C: every FFT profile constructs, transforms and round-trips
///        through the exponent contract dsptap_capi.h states, and the header is C-clean.
// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// Built by the root project when its tests are on (tools/capi/CMakeLists.txt) and registered with
// CTest as dsptap_capi_smoke; not a gtest because it exists to compile the header with a C front
// end and to exercise the ABI the way a binding does (the notebooks' ctypes bridge is such a
// binding). Returns the number of failed checks.

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "dsptap_capi.h"

static int failures = 0;

#define CHECK(cond)                                                                                                    \
    do {                                                                                                               \
        if (!(cond)) {                                                                                                 \
            ++failures;                                                                                                \
            printf("FAILED %s:%d: %s\n", __FILE__, __LINE__, #cond);                                                   \
        }                                                                                                              \
    } while (0)

static const int         k_profiles[] = {DSPTAP_FFT_PROFILE_DOUBLE, DSPTAP_FFT_PROFILE_FLOAT,   DSPTAP_FFT_PROFILE_Q15,
                                         DSPTAP_FFT_PROFILE_Q31,    DSPTAP_FFT_PROFILE_Q15_BFP, DSPTAP_FFT_PROFILE_Q31_BFP};
static const char* const k_names[]    = {"double", "float", "q15", "q31", "q15_bfp", "q31_bfp"};

static int is_fixed_point(int profile) {
    return profile >= DSPTAP_FFT_PROFILE_Q15;
}

static int log2_int(int n) {
    int l = 0;
    while ((1 << l) < n) {
        ++l;
    }
    return l;
}

/// The test signal in fractions of full scale: a cosine at bin 5 plus one at bin 37, peak 0.75.
static void test_signal(double* x, int n) {
    for (int j = 0; j < n; ++j) {
        const double t = 2.0 * 3.14159265358979323846 * (double)j / (double)n;
        x[j]           = 0.5 * cos(5.0 * t + 0.3) + 0.25 * cos(37.0 * t - 1.1);
    }
}

/// Quantize a fraction to a Q format the way the capi documents (round half away from zero,
/// saturating), for the raw-path buffers.
static int64_t quantize(double v, int frac_bits) {
    const double scale = (double)((int64_t)1 << frac_bits);
    const double r     = v < 0.0 ? v * scale - 0.5 : v * scale + 0.5;
    const double lo    = -scale;
    const double hi    = scale - 1.0;
    if (r <= lo) {
        return (int64_t)lo;
    }
    if (r >= hi) {
        return (int64_t)hi;
    }
    return (int64_t)r;
}

/// The native-buffer form of the test signal: the profile's own dtype, quantized as the boundary
/// quantizes.
static void fill_raw(void* raw, int profile, int bytes, const double* x, const double* xq, int n) {
    for (int j = 0; j < n; ++j) {
        if (profile == DSPTAP_FFT_PROFILE_DOUBLE) {
            ((double*)raw)[j] = xq[j];
        }
        else if (profile == DSPTAP_FFT_PROFILE_FLOAT) {
            ((float*)raw)[j] = (float)xq[j];
        }
        else if (bytes == 2) {
            ((int16_t*)raw)[j] = (int16_t)quantize(x[j], 15);
        }
        else {
            ((int32_t*)raw)[j] = (int32_t)quantize(x[j], 31);
        }
    }
}

static void check_profile(int profile, const char* name, int n) {
    dsptap_fft h = dsptap_fft_create(n, profile);
    CHECK(h != NULL);
    if (h == NULL) {
        return;
    }
    CHECK(dsptap_fft_size(h) == n);
    CHECK(dsptap_fft_num_bins(h) == n / 2 + 1);
    CHECK(dsptap_fft_profile(h) == profile);
    const int bytes = dsptap_fft_sample_bytes(h);
    const int bound = dsptap_fft_fixed_scaling_exponent(h);
    switch (profile) {
    case DSPTAP_FFT_PROFILE_DOUBLE:
        CHECK(bytes == 8 && bound == 0);
        break;
    case DSPTAP_FFT_PROFILE_FLOAT:
        CHECK(bytes == 4 && bound == 0);
        break;
    case DSPTAP_FFT_PROFILE_Q15:
    case DSPTAP_FFT_PROFILE_Q15_BFP:
        CHECK(bytes == 2 && bound == log2_int(n));
        break;
    default:
        CHECK(bytes == 4 && bound == log2_int(n) + 1);
        break;
    }

    // The double boundary: forward against the DOUBLE profile on the same quantized input, and
    // the exponent the _exp variant reports.
    double* x      = (double*)malloc((size_t)n * sizeof(double));
    double* xq     = (double*)malloc((size_t)n * sizeof(double));
    double* out    = (double*)malloc((size_t)n * sizeof(double));
    double* golden = (double*)malloc((size_t)n * sizeof(double));
    test_signal(x, n);
    const int frac_bits = bytes == 2 ? 15 : 31;
    for (int j = 0; j < n; ++j) {
        xq[j] = is_fixed_point(profile) ? (double)quantize(x[j], frac_bits) / (double)((int64_t)1 << frac_bits)
                                        : (profile == DSPTAP_FFT_PROFILE_FLOAT ? (double)(float)x[j] : x[j]);
    }
    dsptap_fft g = dsptap_fft_create(n, DSPTAP_FFT_PROFILE_DOUBLE);
    CHECK(g != NULL && dsptap_fft_forward(g, xq, golden) == 0);
    int e = -99;
    CHECK(dsptap_fft_forward_exp(h, x, out, &e) == 0);
    if (is_fixed_point(profile)) {
        CHECK(e >= 0 && e <= bound);
        if (profile == DSPTAP_FFT_PROFILE_Q15 || profile == DSPTAP_FFT_PROFILE_Q31) {
            CHECK(e == bound);
        }
    }
    else {
        CHECK(e == 0);
    }
    // Tolerance: one output LSB at the returned exponent, times a small factor for the kernel's
    // rounding (the battery pins the exact floors; this is a smoke test of the plumbing).
    double max_err = 0.0;
    double peak    = 0.0;
    for (int j = 0; j < n; ++j) {
        const double err = fabs(out[j] - golden[j]);
        if (err > max_err) {
            max_err = err;
        }
        if (fabs(golden[j]) > peak) {
            peak = fabs(golden[j]);
        }
    }
    double tol = 1e-12 * peak;
    if (profile == DSPTAP_FFT_PROFILE_FLOAT) {
        tol = 4e-7 * peak; // float epsilon, peak-normalized (fft.h's parity metric)
    }
    else if (is_fixed_point(profile)) {
        tol = 4.0 * ldexp(1.0, e - frac_bits);
    }
    CHECK(max_err <= tol);
    // dsptap_fft_forward (no exponent out) is the same result.
    double* out2 = (double*)malloc((size_t)n * sizeof(double));
    CHECK(dsptap_fft_forward(h, x, out2) == 0);
    CHECK(memcmp(out, out2, (size_t)n * sizeof(double)) == 0);
    printf("  %-8s n=%d  forward e=%d (bound %d)  max |out - golden| %.3e (tol %.3e)\n", name, n, e, bound, max_err,
           tol);

    // The raw path on the native buffer: forward then inverse, and the round-trip identity
    // x == out * 2^(e_fwd + e_inv + 1 - log2 n) (fixed point) or out * 2/n (floating).
    void* raw = malloc((size_t)n * (size_t)bytes);
    fill_raw(raw, profile, bytes, x, xq, n);
    int e_fwd = -99;
    int e_inv = -99;
    CHECK(dsptap_fft_forward_inplace_raw_exp(h, raw, &e_fwd) == 0);
    CHECK(dsptap_fft_inverse_inplace_raw_exp(h, raw, &e_inv) == 0);
    if (is_fixed_point(profile)) {
        CHECK(e_fwd == e);
        CHECK(e_inv >= 0 && e_inv <= bound);
    }
    else {
        CHECK(e_fwd == 0 && e_inv == 0);
    }
    const double scale =
        is_fixed_point(profile) ? ldexp(1.0, e_fwd + e_inv + 1 - log2_int(n) - frac_bits) : 2.0 / (double)n;
    double rt_err = 0.0;
    for (int j = 0; j < n; ++j) {
        double v;
        if (profile == DSPTAP_FFT_PROFILE_DOUBLE) {
            v = ((double*)raw)[j];
        }
        else if (profile == DSPTAP_FFT_PROFILE_FLOAT) {
            v = (double)((float*)raw)[j];
        }
        else if (bytes == 2) {
            v = (double)((int16_t*)raw)[j];
        }
        else {
            v = (double)((int32_t*)raw)[j];
        }
        const double err = fabs(v * scale - xq[j]);
        if (err > rt_err) {
            rt_err = err;
        }
    }
    // The reconstruction's resolution is one native LSB times the round-trip scale; the battery
    // (tests/test_fft_fixed.cpp, RoundTripReconstructsInputPerPolicy) pins the exact numbers
    // (up to 41.5 reconstructed LSB for Q31 block floating, at the DC index); this is plumbing.
    const double rt_tol =
        is_fixed_point(profile) ? 128.0 * scale : (profile == DSPTAP_FFT_PROFILE_FLOAT ? 1e-6 : 1e-14);
    CHECK(rt_err <= rt_tol);
    printf("  %-8s n=%d  raw round trip e_fwd=%d e_inv=%d  max |x - out * 2^k| %.3e (tol %.3e)\n", name, n, e_fwd,
           e_inv, rt_err, rt_tol);

    // The status-returning raw functions return the same exponent (0 for the floating profiles,
    // which is what they returned before the fixed-point profiles existed).
    fill_raw(raw, profile, bytes, x, xq, n);
    const int status = dsptap_fft_forward_inplace_raw(h, raw);
    CHECK(status == e_fwd);
    const int status_inv = dsptap_fft_inverse_inplace_raw(h, raw);
    CHECK(status_inv == e_inv);

    free(raw);
    free(out2);
    free(golden);
    free(out);
    free(xq);
    free(x);
    dsptap_fft_destroy(g);
    dsptap_fft_destroy(h);
}

int main(void) {
    printf("dsptap_capi smoke: float32 backend %s\n", dsptap_fft_backend());

    // Bad arguments: NULL on failure, -1 from every handle function, never a crash.
    CHECK(dsptap_fft_create(3, DSPTAP_FFT_PROFILE_DOUBLE) == NULL);
    CHECK(dsptap_fft_create(0, DSPTAP_FFT_PROFILE_Q15) == NULL);
    CHECK(dsptap_fft_create(512, 6) == NULL);
    CHECK(dsptap_fft_create(512, -1) == NULL);
    CHECK(dsptap_fft_size(NULL) == -1);
    CHECK(dsptap_fft_num_bins(NULL) == -1);
    CHECK(dsptap_fft_profile(NULL) == -1);
    CHECK(dsptap_fft_sample_bytes(NULL) == -1);
    CHECK(dsptap_fft_fixed_scaling_exponent(NULL) == -1);
    {
        double     buf[16] = {0};
        int        e       = 0;
        dsptap_fft h       = dsptap_fft_create(16, DSPTAP_FFT_PROFILE_Q31_BFP);
        CHECK(h != NULL);
        CHECK(dsptap_fft_forward(NULL, buf, buf) == -1);
        CHECK(dsptap_fft_forward(h, NULL, buf) == -1);
        CHECK(dsptap_fft_inverse(h, buf, NULL) == -1);
        CHECK(dsptap_fft_forward_exp(h, buf, buf, NULL) == -1);
        CHECK(dsptap_fft_inverse_exp(h, buf, buf, NULL) == -1);
        CHECK(dsptap_fft_forward_inplace_raw(h, NULL) == -1);
        CHECK(dsptap_fft_inverse_inplace_raw(NULL, buf) == -1);
        CHECK(dsptap_fft_forward_inplace_raw_exp(h, buf, NULL) == -1);
        CHECK(dsptap_fft_inverse_inplace_raw_exp(NULL, buf, &e) == -1);
        // Silence under block floating point: e == 0, output all zero.
        CHECK(dsptap_fft_forward_exp(h, buf, buf, &e) == 0 && e == 0);
        for (int j = 0; j < 16; ++j) {
            CHECK(buf[j] == 0.0);
        }
        dsptap_fft_destroy(h);
        dsptap_fft_destroy(NULL);
    }

    for (size_t p = 0; p < sizeof(k_profiles) / sizeof(k_profiles[0]); ++p) {
        check_profile(k_profiles[p], k_names[p], 512);
        check_profile(k_profiles[p], k_names[p], 64);
    }

    if (failures == 0) {
        printf("dsptap_capi smoke: all checks passed\n");
    }
    else {
        printf("dsptap_capi smoke: %d check(s) FAILED\n", failures);
    }
    return failures == 0 ? 0 : 1;
}
