// Reproducer: vDSP_fft_zrip selects a substantially less accurate kernel when
// the DSPSplitComplex halves are 64-byte aligned.
//
// Build and run:
//   cc -O2 vdsp_fft_alignment_repro.c -framework Accelerate -o repro && ./repro
//
// The program performs the SAME forward real FFT of the SAME input through the
// SAME FFTSetup twice, changing only where the split-complex buffers sit:
// once with both halves 64-byte aligned, once 32 bytes past that boundary.
// Each result is compared against vDSP's OWN double-precision real FFT
// (vDSP_fft_zripD) of the same input, so the reference is Accelerate itself.
//
// A naive O(N^2) DFT is NOT usable as the reference here: its summation error
// grows as N*eps*peak (~2e-10 at this scale), which is a thousand times larger
// than the ~1e-13 bins being measured. That mistake makes both placements look
// equally wrong. vDSP_fft_zripD is O(N log N), so its error stays around
// 1e-12, and it is the natural double-precision counterpart of the routine
// under test.
//
// Input is a pure sinusoid at exactly bin N/16, so all but one bin of the true
// spectrum is at the level of the input's own float32 quantisation — the
// regime where the two kernels differ. Broadband input shows no difference,
// which is why this is easy to miss.

#include <Accelerate/Accelerate.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

enum { kN = 2048, kBins = kN / 2 + 1 };

static int cmp_double(const void* a, const void* b) {
    const double x = *(const double*)a, y = *(const double*)b;
    return (x > y) - (x < y);
}

// Reference: vDSP's own double-precision real FFT of the same input.
static void reference_magnitudes(const float* x, double* mag) {
    static double rpd[kN / 2], ipd[kN / 2], in[kN];
    for (int i = 0; i < kN; ++i) {
        in[i] = (double)x[i];
    }
    FFTSetupD setup = vDSP_create_fftsetupD(11, kFFTRadix2);
    DSPDoubleSplitComplex sc = {rpd, ipd};
    vDSP_ctozD((const DSPDoubleComplex*)in, 2, &sc, 1, kN / 2);
    vDSP_fft_zripD(setup, &sc, 1, 11, kFFTDirection_Forward);
    mag[0]      = 0.5 * fabs(rpd[0]);
    mag[kN / 2] = 0.5 * fabs(ipd[0]);
    for (int k = 1; k < kN / 2; ++k) {
        mag[k] = 0.5 * hypot(rpd[k], ipd[k]);
    }
    vDSP_destroy_fftsetupD(setup);
}

// One forward transform with the split halves placed at `skew` bytes past a
// 64-byte boundary, writing per-bin magnitudes.
static void transform_at(FFTSetup setup, unsigned char* slab, size_t skew, const float* input, double* out_mag) {
    float* in = (float*)(slab + skew);
    float* rp = (float*)(slab + skew + kN * sizeof(float) + 4096);
    float* ip = rp + kN / 2;
    memcpy(in, input, kN * sizeof(float));

    DSPSplitComplex sc = {rp, ip};
    vDSP_ctoz((const DSPComplex*)in, 2, &sc, 1, kN / 2);
    vDSP_fft_zrip(setup, &sc, 1, 11 /* log2(2048) */, kFFTDirection_Forward);

    // vDSP packs DC in rp[0], Nyquist in ip[0], and carries a factor of 2.
    out_mag[0]        = 0.5 * fabs((double)rp[0]);
    out_mag[kN / 2]   = 0.5 * fabs((double)ip[0]);
    for (int k = 1; k < kN / 2; ++k) {
        out_mag[k] = 0.5 * hypot((double)rp[k], (double)ip[k]);
    }
    printf("  rp = %p  (mod 64 = %2lu)\n", (void*)rp, (unsigned long)((uintptr_t)rp % 64));
}

static double median_relative_error(const double* got, const double* ref) {
    double* rel = malloc(kBins * sizeof(double));
    int     n   = 0;
    for (int k = 0; k < kBins; ++k) {
        if (ref[k] > 0.0) {
            rel[n++] = fabs(got[k] - ref[k]) / ref[k];
        }
    }
    qsort(rel, n, sizeof(double), cmp_double);
    const double m = rel[n / 2];
    free(rel);
    return m;
}

int main(void) {
    static float  input[kN];
    static double ref[kBins], got[kBins];

    // Exactly on bin N/16 — no leakage, so the other bins stay at the
    // quantisation floor.
    for (int i = 0; i < kN; ++i) {
        input[i] = (float)(0.5 * sin(2.0 * M_PI * (kN / 16.0) * i / (double)kN));
    }

    reference_magnitudes(input, ref);

    FFTSetup setup = vDSP_create_fftsetup(11, kFFTRadix2);
    if (!setup) {
        fprintf(stderr, "vDSP_create_fftsetup failed\n");
        return 1;
    }
    void* slab = NULL;
    if (posix_memalign(&slab, 4096, kN * sizeof(float) * 8 + 3 * 4096) != 0) {
        fprintf(stderr, "posix_memalign failed\n");
        return 1;
    }

    printf("vDSP_fft_zrip, N=%d, identical input and FFTSetup, only buffer placement differs\n\n", kN);

    printf("split halves 64-byte aligned:\n");
    transform_at(setup, slab, 0, input, got);
    const double aligned = median_relative_error(got, ref);
    printf("  median per-bin relative error vs double DFT: %.6g\n\n", aligned);

    printf("split halves 32 bytes past the boundary:\n");
    transform_at(setup, slab, 32, input, got);
    const double skewed = median_relative_error(got, ref);
    printf("  median per-bin relative error vs double DFT: %.6g\n\n", skewed);

    printf("ratio (aligned / skewed): %.6g\n", aligned / skewed);

    free(slab);
    vDSP_destroy_fftsetup(setup);
    return 0;
}
