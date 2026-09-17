// SPDX-License-Identifier: MIT
// Copyright 2025-2026 Timothy Place and the DspTap contributors.
//
// One-shot gtest main for the bare-metal emulated targets (Cortex-M4 soft-float
// and M4F on qemu mps2-an386, Cortex-M33 on mps2-an505, Cortex-M55 on
// mps3-an547). tap_dsp_add_gtest_executable() in tests/CMakeLists.txt links
// this into every test executable when TAP_DSP_BARE_METAL is set. There is no
// argv on the target, so the selection is baked in at build time:
//
//   TAP_DSP_BARE_METAL_FILTER  the gtest filter (the helper's MAIN_FILTER);
//                              unset means "*": everything in the binary runs.
//
// The filter is NEGATIVE by policy (docs/audit-fft-and-code-smells.md, Part
// 13): a suite compiled into a test executable runs on the target unless it is
// excluded by name, and every exclusion is a written budget decision next to
// the source list. The positive filter this file had first would have silently
// dropped every suite added afterwards, and no count check can tell an
// omission from a small suite.
//
// Completion contract, read by CTest (see the PASS/FAIL regular expressions in
// tests/CMakeLists.txt):
//   TAP_DSP_TESTS_COMPLETE rc=<n> selected=<n> skipped=<n>
// printed only if the run reached the end, so a crash after gtest's own summary
// cannot register as a pass. rc is 1 when any test failed, when nothing was
// selected (a filter typo must not pass green), and when any test was skipped:
// gtest 1.14 cannot fail on GTEST_SKIP by itself, and a skipped gate is not a
// passed gate — exclude the test by name in the filter, or compile it out,
// instead of skipping at run time. A fault before this line prints
// TAP_DSP_TESTS_FAULT from platform/cortexm_startup.c and exits, which CTest
// also treats as failure.
//
// Budget: under 5 minutes of emulation per leg. The measured wall time of each
// leg is in its CI job log (ctest -V, test tap_dsp_tests_emulated); it moves by
// +-50 % between identical runs on shared runners (TCG wall is noise), so it
// is a budget signal, never a number to pin. Instruction counts (Stage 1b) are
// the ratchet.
#include <cstdio>

#include <gtest/gtest.h>

#ifndef TAP_DSP_BARE_METAL_FILTER
#define TAP_DSP_BARE_METAL_FILTER "*"
#endif

int main() {
    ::testing::GTEST_FLAG(filter) = TAP_DSP_BARE_METAL_FILTER;
    ::testing::InitGoogleTest();
    int rc = RUN_ALL_TESTS();
    // gtest applies the filter inside RUN_ALL_TESTS, so the counts are read
    // afterwards.
    const ::testing::UnitTest* const unit     = ::testing::UnitTest::GetInstance();
    const int                        selected = unit->test_to_run_count();
    const int                        skipped  = unit->skipped_test_count();
    if (selected == 0) {
        std::printf("no tests selected: filter \"%s\" is broken\n", TAP_DSP_BARE_METAL_FILTER);
        rc = 1;
    }
    if (skipped > 0) {
        std::printf("%d test(s) skipped: a skipped gate is a failed gate on the target; "
                    "exclude it in MAIN_FILTER or compile it out\n",
                    skipped);
        rc = 1;
    }
    std::printf("TAP_DSP_TESTS_COMPLETE rc=%d selected=%d skipped=%d\n", rc, selected, skipped);
    return rc;
}
