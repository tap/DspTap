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

# ==============================================================================
# BEGIN Stage 2a additions (docs/audit-fft-and-code-smells.md, Part 9): the
# independent oracle and the real-time guard join the main executable; the
# Ooura bit-identity gate gets its own targets so it can own its compiler flags.
# Keep this block self-contained and at the end of the file.
# ==============================================================================

target_sources(tap_dsp_tests PRIVATE
    test_fft_oracle.cpp
    test_fft_rt.cpp)

# ------------------------------------------------------------------------------
# Two private reference builds of the vendored C. Neither is the tap_dsp_fft
# the library links: that target disappears at Stage 2c (it survives only under
# TAP_DSP_FFT_CMSIS, where it also carries CMSIS objects), and the gate's whole
# point is to compare against the C independently of how the library builds it.
# The parity executables link one of these and NOT tap::dsp, so their
# rdft/rdft_f resolve here and no backend define (vDSP on Apple) reaches them.
#
# Stage 2c note: fftsg.c moves to tests/reference/ooura/ (D6); this gate needs
# fftsg_float.c to move with it, or its float side becomes port-vs-port.
#
#   ..._nocontract  fp-contraction OFF. Bit identity between the C and its C++
#                   transliteration holds only if both compile to the same
#                   sequence of IEEE operations, and CMake gives neither side an
#                   fp-contract setting: gcc contracts C in gnu17 mode, g++
#                   contracts C++ in every mode, clang contracts everywhere
#                   (Part 4, "fp-contraction policy"). MSVC: nothing is passed;
#                   its default /fp:precise has not contracted since VS 2022
#                   17.0 made /fp:contract opt-in (earlier x64 /arch:AVX2 and
#                   ARM64 builds did contract under /fp:precise; windows-latest
#                   is well past 17.0).
#   ..._default     the same two files at default flags, for the informational
#                   target, so it measures "port at default flags vs C at
#                   default flags" and nothing else.
# ------------------------------------------------------------------------------
set(_tap_dsp_nocontract_c   "$<$<C_COMPILER_ID:GNU,Clang,AppleClang>:-ffp-contract=off>")
set(_tap_dsp_nocontract_cxx "$<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-ffp-contract=off>")
set(_tap_dsp_reference_sources
    ${PROJECT_SOURCE_DIR}/third_party/ooura/fftsg.c
    ${PROJECT_SOURCE_DIR}/third_party/ooura/fftsg_float.c)

add_library(tap_dsp_fft_reference_nocontract STATIC ${_tap_dsp_reference_sources})
target_compile_options(tap_dsp_fft_reference_nocontract PRIVATE ${_tap_dsp_nocontract_c})

add_library(tap_dsp_fft_reference_default STATIC ${_tap_dsp_reference_sources})

# Largest transform the parity suite runs. Hosted, the default covers every
# consumer geometry (2^20: five 8 MB buffers plus tables). Cross-compiled, the
# emulated Cortex-M legs (Part 10) get 4096, the largest size that fits the
# MPS2/MPS3 data regions; a toolchain file may still set its own cache value,
# which this default does not override.
if(CMAKE_CROSSCOMPILING)
    set(_tap_dsp_parity_max_n_default 4096)
else()
    set(_tap_dsp_parity_max_n_default 1048576)
endif()
set(TAP_DSP_PARITY_MAX_N ${_tap_dsp_parity_max_n_default} CACHE STRING
    "Largest FFT size the Ooura parity suite runs (power of two >= 4; 4096 on emulated targets)")

# THE GATE: -ffp-contract=off on both sides, memcmp bit identity.
add_executable(tap_dsp_fft_parity test_fft_parity_ooura.cpp)
target_include_directories(tap_dsp_fft_parity PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_compile_features(tap_dsp_fft_parity PRIVATE cxx_std_20)
target_compile_options(tap_dsp_fft_parity PRIVATE ${_tap_dsp_nocontract_cxx})
target_compile_definitions(tap_dsp_fft_parity PRIVATE TAP_DSP_PARITY_MAX_N=${TAP_DSP_PARITY_MAX_N})
target_link_libraries(tap_dsp_fft_parity PRIVATE
    tap_dsp_fft_reference_nocontract
    tap_dsp_warnings)

# INFORMATIONAL: the same source at default flags on both sides, measuring the
# max-ulp deviation per N and never failing. It is not a pass/fail gate with an
# assumed bound because a different fusion choice per butterfly stage
# accumulates over log2 N stages; the number goes on the record per platform
# instead (Part 6, item N3). Two channels carry it, because ctest hides the
# stdout of a passing test under --output-on-failure: (1) CI runs this label
# with -V as its own step (Part 13), (2) --gtest_output=xml writes the
# RecordProperty values to parity-ulp.xml in the build tree for upload. The
# report is ONE test covering both precisions so a single invocation produces
# the complete XML (each ctest entry is a separate process writing that file).
add_executable(tap_dsp_fft_parity_default_flags test_fft_parity_ooura.cpp)
target_include_directories(tap_dsp_fft_parity_default_flags PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_compile_features(tap_dsp_fft_parity_default_flags PRIVATE cxx_std_20)
target_compile_definitions(tap_dsp_fft_parity_default_flags PRIVATE
    TAP_DSP_PARITY_INFORMATIONAL
    TAP_DSP_PARITY_MAX_N=${TAP_DSP_PARITY_MAX_N})
target_link_libraries(tap_dsp_fft_parity_default_flags PRIVATE
    tap_dsp_fft_reference_default
    tap_dsp_warnings)

include(GoogleTest)
target_link_libraries(tap_dsp_fft_parity PRIVATE GTest::gtest_main)
gtest_discover_tests(tap_dsp_fft_parity
    DISCOVERY_TIMEOUT 120
    PROPERTIES TIMEOUT 900 LABELS parity)
target_link_libraries(tap_dsp_fft_parity_default_flags PRIVATE GTest::gtest_main)
gtest_discover_tests(tap_dsp_fft_parity_default_flags
    EXTRA_ARGS --gtest_output=xml:${CMAKE_BINARY_DIR}/parity-ulp.xml
    DISCOVERY_TIMEOUT 120
    PROPERTIES TIMEOUT 900 LABELS parity)

# ==============================================================================
# END Stage 2a additions
# ==============================================================================
