// Test runner main for the bare-metal emulated targets (Cortex-M4 soft-float
// and M4F on qemu mps2-an386, Cortex-M33 on mps2-an505, Cortex-M55 on
// mps3-an547): there is no argv on the target, so the emulation-appropriate
// selection is baked in. This is a POSITIVE filter — a suite is on the target
// because it is named here, never by default.
//
// Selected:
//   - every float32 typed suite (/0 in sample_types order): real_fft, yin,
//     psola, pvoc, log_mel, nn — the embedded profile these legs exist for;
//   - the double real_fft contract suite (/1) and RealFftCrossPrecision: small
//     N, and on the M4/M33 legs a soft-float double correctness check;
//   - the FFT backend suites (parity, alignment stability, tonal accuracy):
//     Ooura-vs-Ooura identity on the M4/M33 legs, and on the M55 leg the
//     first place the CMSIS-DSP Helium backend actually runs against Ooura;
//   - the fir_kernels typed suite over float/Q15/Q31 — the SMLALD dual-MAC
//     path is compiled and run for the first time (both M4 and M33 define
//     __ARM_FEATURE_DSP) — plus the sample_traits and quantize batteries;
//   - the decimate battery (float and Q15 profiles; the numpy pins are
//     compiled-in arrays, no filesystem needed);
//   - every float-tracks-double oracle check (each costs one double instance
//     at a small geometry);
//   - the cheap Kaiser design tests and the sine-fit instrument floors.
//
// Excluded: the double-typed psola/pvoc/yin/log_mel/nn suites (target-
// independent math already covered on every host, and soft-float double on
// three of the four legs), the double Kaiser prototype searches (~0.7 s of
// double on a desktop host, i.e. minutes of soft-float emulation) and the
// double multitone instrument fits.
//
// Budget: under 5 minutes of emulation per leg. The whole hosted battery is
// ~1.4 s on a desktop core; the selection here drops most of the double
// work, and the M4 soft-float leg (every float op a library call) is the one
// that sets the ceiling. Measured times per leg are in the CI job logs
// (ctest prints the wall time of tap_dsp_tests_emulated).
// SPDX-License-Identifier: MIT
// Copyright 2025-2026 Timothy Place and the DspTap contributors.
#include <cstdio>

#include <gtest/gtest.h>

int main() {
    // Typed-suite naming: /0 = float, /1 = double (sample_types order);
    // fir_kernels_test/0,1,2 = float, int16 (Q15), int32 (Q31).
    ::testing::GTEST_FLAG(filter) = "real_fft_test/0.*:real_fft_test/1.*:RealFftCrossPrecision.*:RealFftFloatIO.*:"
                                    "CertifiedGeometries/fft_backend_parity.*:"
                                    "CertifiedGeometries/fft_alignment_stability.*:"
                                    "CertifiedGeometries/fft_tonal_accuracy.*:"
                                    "fir_kernels_test/*.*:SampleTraits.*:Quantize.*:"
                                    "Kaiser.BesselI0ReferenceValues:Kaiser.BetaReferenceValues:"
                                    "Kaiser.TapEstimateMatchesHarrisFormula:Kaiser.CompensatedBranchSumsAreUniform:"
                                    "Kaiser.SolveDenseSolvesKnownSystem:"
                                    "SineAnalysis.*:"
                                    "Decimate.*:"
                                    "log_mel_test/0.*:LogMelCrossPrecision.*:"
                                    "nn_test/0.*:NnCrossPrecision.*:"
                                    "yin_test/0.*:yin_cross_precision.*:"
                                    "psola_test/0.*:psola_cross_precision.*:"
                                    "pvoc_test/0.*:pvoc_cross_precision.*";
    ::testing::InitGoogleTest();
    const int rc = RUN_ALL_TESTS();
    // A filter typo selects zero tests and RUN_ALL_TESTS() returns 0 — an
    // empty run must not pass green. Checked after the run because gtest
    // only applies the filter inside RUN_ALL_TESTS. The on-target selection
    // is ~100 tests; 60 leaves headroom for legitimate removals without
    // masking a typo.
    const int selected = ::testing::UnitTest::GetInstance()->test_to_run_count();
    if (selected < 60) {
        std::printf("only %d tests selected (expected >= 60): filter is broken\n", selected);
        std::printf("TAP_DSP_TESTS_COMPLETE rc=1\n");
        return 1;
    }
    // CTest's pass criterion: printed only if we get all the way here, so a
    // crash after gtest's summary cannot register as a pass.
    std::printf("TAP_DSP_TESTS_COMPLETE rc=%d\n", rc);
    return rc;
}
