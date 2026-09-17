include(FetchContent)
include(GoogleTest)

# ------------------------------------------------------------------------------
# Two test build modes, one registration helper.
#
# Hosted (default): each test executable links GTest::gtest_main and its tests
# are registered one by one through gtest_discover_tests.
#
# TAP_DSP_BARE_METAL (set by the four cmake/arm-cortex-*.cmake toolchains):
# one-shot mode for Cortex-M under qemu-system-arm. There is no argv on the
# target and semihosting does not reliably propagate exit codes, so each test
# executable gets tests/bare_metal_main.cpp as its main, runs its whole
# selection in one process, and is registered as ONE CTest test judged on the
# text it prints: PASS on "TAP_DSP_TESTS_COMPLETE rc=0 ...", FAIL on any gtest
# "[  FAILED  ]" line or on the "TAP_DSP_TESTS_FAULT" marker the startup file
# (platform/cortexm_startup.c) prints from a fault handler before exiting.
#
#   tap_dsp_add_gtest_executable(<target>
#       SOURCES <src>...
#       [MAIN_FILTER <gtest filter>]
#       [LABELS <label>...])
#
# creates <target> from SOURCES and does the mode-specific wiring. Add include
# directories, libraries, options and definitions afterwards with the usual
# target_* commands; do NOT link gtest_main and do NOT call
# gtest_discover_tests yourself (the helper does, in hosted mode only —
# discovery would run the cross binary under the emulator at build time with
# argv it cannot see). MAIN_FILTER is a gtest filter applied only in
# bare-metal mode, default "*": everything compiled into the binary runs on
# the target. Use a NEGATIVE filter ("*-suite_a.*:suite_b.Test") to exclude
# what the emulation budget cannot afford, never a positive one — a positive
# filter silently drops every suite added later, and no count check can tell
# an omission from a small suite. LABELS go on the discovered tests (hosted)
# and on the one-shot test (bare metal). The GTEST_HAS_* feature macros bare
# metal needs are directory-wide, so a caller needs nothing else.
#
# Skips are failures on the target: bare_metal_main.cpp reports rc=1 when any
# test was GTEST_SKIP'ped, so a suite that would skip at a size the target
# cannot afford must be excluded by name in MAIN_FILTER or compiled out
# (e.g. under TAP_DSP_PARITY_MAX_N, which the toolchains default to 4096).
# ------------------------------------------------------------------------------
if(NOT TAP_DSP_BARE_METAL)
    find_package(Threads QUIET)
endif()
if(NOT Threads_FOUND)
    set(gtest_disable_pthreads ON CACHE BOOL "" FORCE)
endif()
if(TAP_DSP_BARE_METAL)
    # Value-checked feature macros only: GTEST_HAS_DEATH_TEST is tested with
    # #ifdef in gtest 1.14, so it must stay undefined (no GTEST_OS_* is
    # detected on bare metal, which disables it).
    add_compile_definitions(GTEST_HAS_POSIX_RE=0 GTEST_HAS_STREAM_REDIRECTION=0
                            GTEST_HAS_FILE_SYSTEM=0)
endif()

FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    # Commit pin, not the movable tag: tags can be re-pointed upstream.
    GIT_TAG f8d7d77c06936315286eb55f8de22cd23c188571 # v1.14.0
    FIND_PACKAGE_ARGS NAMES GTest)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

function(tap_dsp_add_gtest_executable target)
    cmake_parse_arguments(PARSE_ARGV 1 arg "" "MAIN_FILTER" "SOURCES;LABELS")
    if(arg_UNPARSED_ARGUMENTS)
        message(FATAL_ERROR
            "tap_dsp_add_gtest_executable(${target}): unknown arguments '${arg_UNPARSED_ARGUMENTS}'")
    endif()
    if(NOT arg_SOURCES)
        message(FATAL_ERROR "tap_dsp_add_gtest_executable(${target}): SOURCES is required")
    endif()
    add_executable(${target} ${arg_SOURCES})
    if(TAP_DSP_BARE_METAL)
        target_sources(${target} PRIVATE ${CMAKE_CURRENT_FUNCTION_LIST_DIR}/bare_metal_main.cpp)
        target_link_libraries(${target} PRIVATE GTest::gtest)
        if(arg_MAIN_FILTER)
            target_compile_definitions(${target} PRIVATE
                "TAP_DSP_BARE_METAL_FILTER=\"${arg_MAIN_FILTER}\"")
        endif()
        # CTest prepends CMAKE_CROSSCOMPILING_EMULATOR (the toolchain's qemu
        # line) to a target-name COMMAND automatically. TIMEOUT is ~6x the
        # slowest measured leg; a fault exits in seconds via the marker.
        add_test(NAME ${target}_emulated COMMAND ${target})
        set_tests_properties(${target}_emulated PROPERTIES
            PASS_REGULAR_EXPRESSION "TAP_DSP_TESTS_COMPLETE rc=0 "
            FAIL_REGULAR_EXPRESSION "\\[  FAILED  \\];TAP_DSP_TESTS_FAULT"
            TIMEOUT 600)
        if(arg_LABELS)
            set_tests_properties(${target}_emulated PROPERTIES LABELS "${arg_LABELS}")
        endif()
    else()
        target_link_libraries(${target} PRIVATE GTest::gtest_main)
        set(_tap_dsp_props TIMEOUT 900)
        if(arg_LABELS)
            list(APPEND _tap_dsp_props LABELS "${arg_LABELS}")
        endif()
        gtest_discover_tests(${target} DISCOVERY_TIMEOUT 120 PROPERTIES ${_tap_dsp_props})
    endif()
endfunction()

# ------------------------------------------------------------------------------
# The main battery. Adding a source here puts its suites on the four QEMU legs
# as well; the MAIN_FILTER below is the complete list of what does not run
# there, each entry a budget decision:
#   - the double-typed psola/pvoc/yin/log_mel/nn suites (/1): target-
#     independent math already covered on every host, and soft-float double on
#     three of the four legs;
#   - the six Kaiser design searches (the four prototype specs plus the
#     rational-phase-count and compensated-at-16k specs): 30-340 ms of double
#     each on a desktop core, i.e. tens of seconds of soft-float emulation
#     apiece, pinning target-independent design numbers;
#   - MultitoneAnalysis: the double joint-fit instrument floor, same reason.
# Everything else — every float32 suite, the double real_fft contract suite,
# every float-tracks-double oracle check, the FFT backend suites (CMSIS-DSP
# Helium against Ooura on the M55), fir_kernels over float/Q15/Q31, decimate,
# sample_traits, quantize, the cheap Kaiser tests, SineAnalysis, and whatever
# a later PR appends — runs on target.
# ------------------------------------------------------------------------------
tap_dsp_add_gtest_executable(tap_dsp_tests
    SOURCES
        test_analysis.cpp
        test_decimate.cpp
        test_fft.cpp
        test_fft_arith.cpp
        test_fft_backend.cpp
        test_fir_kernels.cpp
        test_kaiser.cpp
        test_log_mel.cpp
        test_nn.cpp
        test_psola.cpp
        test_pvoc.cpp
        test_quantize.cpp
        test_sample_traits.cpp
        test_yin.cpp
    MAIN_FILTER
        "*-psola_test/1.*:pvoc_test/1.*:yin_test/1.*:log_mel_test/1.*:nn_test/1.*:Kaiser.FastPrototypeMeetsSpec:Kaiser.BalancedPrototypeMeetsSpec:Kaiser.TransparentPrototypeMeetsSpec:Kaiser.EconomyPrototypeMeetsSpec:Kaiser.RationalPhaseCountMeetsSpec:Kaiser.CompensatedSpecsHoldAt16k:MultitoneAnalysis.*")
target_link_libraries(tap_dsp_tests PRIVATE
    tap::dsp
    tap_dsp_warnings)
