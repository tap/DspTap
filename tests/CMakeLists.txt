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
# every float-tracks-double oracle check, the FFT engine suites
# (test_fft_backend.cpp, typed over the engines the leg can build: the
# split-radix engine everywhere and CMSIS-DSP Helium beside it on the M55;
# test_fft_engine.cpp, the Stage 4 engine-parameter contract — default
# selection, ABI tag, per-engine size ranges, shareability) and the routing
# proof (test_fft_routing.cpp: basic_real_fft is its engine byte for byte,
# for the split-radix engine named explicitly on every leg and for the
# selected default), fir_kernels over float/Q15/Q31, decimate,
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
        test_fft_engine.cpp
        test_fft_fixed.cpp
        test_fft_routing.cpp
        test_fir_kernels.cpp
        test_kaiser.cpp
        test_log_mel.cpp
        test_nn.cpp
        test_psola.cpp
        test_pvoc.cpp
        test_quantize.cpp
        test_sample_traits.cpp
        test_spectrum.cpp
        test_yin.cpp
    MAIN_FILTER
        "*-psola_test/1.*:pvoc_test/1.*:yin_test/1.*:log_mel_test/1.*:nn_test/1.*:Kaiser.FastPrototypeMeetsSpec:Kaiser.BalancedPrototypeMeetsSpec:Kaiser.TransparentPrototypeMeetsSpec:Kaiser.EconomyPrototypeMeetsSpec:Kaiser.RationalPhaseCountMeetsSpec:Kaiser.CompensatedSpecsHoldAt16k:MultitoneAnalysis.*")
target_link_libraries(tap_dsp_tests PRIVATE
    tap::dsp
    tap_dsp_warnings)

# ==============================================================================
# BEGIN Stage 2a additions (docs/audit-fft-and-code-smells.md, Part 9): the
# independent oracle and the real-time guard join the main executable; the
# Ooura bit-identity gate gets its own targets so it can own its compiler flags
# (and, since Stage 2c, the reference C they compile).
# Registration goes through tap_dsp_add_gtest_executable above, so the two
# targets are hosted (discovered) or bare-metal one-shot exactly like
# tap_dsp_tests. Keep this block self-contained and at the end of the file.
# ==============================================================================

# Largest transform the FFT sweeps run. Hosted, the default covers every
# consumer geometry (2^20: five 8 MB buffers plus tables). Cross-compiled, the
# emulated Cortex-M legs (Part 10) get 4096, what fits the MPS2/MPS3 data
# regions; the toolchain files set the same value in the cache themselves and
# a plain cache set here never overrides theirs or a -D on the command line.
# The parity sweep, its 2^20 test (compiled out below the cap: a GTEST_SKIP is
# rc=1 on the target) and the oracle's closed-form sweep all read it.
if(CMAKE_CROSSCOMPILING)
    set(_tap_dsp_parity_max_n_default 4096)
else()
    set(_tap_dsp_parity_max_n_default 1048576)
endif()
set(TAP_DSP_PARITY_MAX_N ${_tap_dsp_parity_max_n_default} CACHE STRING
    "Largest FFT size the Ooura parity and oracle sweeps run (power of two >= 4; 4096 on emulated targets)")

target_sources(tap_dsp_tests PRIVATE
    test_fft_oracle.cpp
    test_fft_rt.cpp)
target_compile_definitions(tap_dsp_tests PRIVATE TAP_DSP_PARITY_MAX_N=${TAP_DSP_PARITY_MAX_N})

# ------------------------------------------------------------------------------
# Two private reference builds of the reference C, tests/reference/ooura/
# (fftsg.c and its float instantiation fftsg_float.c, moved there from
# third_party/ooura/ at Stage 2c, Decision D6; readme.txt stays at
# third_party/ooura/readme.txt as the license record). Nothing that ships
# compiles these files: since 2c the library is header-only (tap_dsp_fft
# exists only under TAP_DSP_FFT_CMSIS and carries CMSIS objects alone), and
# the gate's whole point is to compare against the C independently of how the
# library builds. The parity executables link one of these and NOT tap::dsp,
# so their rdft/rdft_f resolve here (declared by tests/reference/ooura_rdft.h)
# and no backend define (vDSP on Apple, CMSIS on the M55) reaches them.
# Both files move together because the float half of the gate compares
# against rdft_f; without fftsg_float.c it would be port-vs-port.
#
# EXPORT_COMPILE_COMMANDS OFF: the tidy gate (.github/workflows/style.yml and
# its local mirror scripts/tidy.sh, a taphouse file) runs over every TU in
# the compile database except third_party/ and fetched deps; the reference C
# is vendored code that was excluded by that path until it moved, and is
# excluded by this property now. clang-format is disabled for the directory
# by tests/reference/ooura/.clang-format for the same reason; the declaration
# header the gate includes (tests/reference/ooura_rdft.h) is DspTap's own and
# sits above that directory so the house formatter still covers it.
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
    ${CMAKE_CURRENT_SOURCE_DIR}/reference/ooura/fftsg.c
    ${CMAKE_CURRENT_SOURCE_DIR}/reference/ooura/fftsg_float.c)

add_library(tap_dsp_fft_reference_nocontract STATIC ${_tap_dsp_reference_sources})
target_compile_options(tap_dsp_fft_reference_nocontract PRIVATE ${_tap_dsp_nocontract_c})
set_target_properties(tap_dsp_fft_reference_nocontract PROPERTIES EXPORT_COMPILE_COMMANDS OFF)

add_library(tap_dsp_fft_reference_default STATIC ${_tap_dsp_reference_sources})
set_target_properties(tap_dsp_fft_reference_default PROPERTIES EXPORT_COMPILE_COMMANDS OFF)

# THE GATE: -ffp-contract=off on both sides, memcmp bit identity. On the QEMU
# legs it runs whole (no MAIN_FILTER): at N <= 4096 the sweep is milliseconds.
tap_dsp_add_gtest_executable(tap_dsp_fft_parity
    SOURCES test_fft_parity_ooura.cpp
    LABELS parity)
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
# stdout of a passing test under --output-on-failure: (1) CI runs the parity
# label with -V as its own step (ci.yml, Part 13); (2) the JUnit XML below.
# The report is ONE test covering both precisions so a single invocation
# produces the complete XML.
tap_dsp_add_gtest_executable(tap_dsp_fft_parity_default_flags
    SOURCES test_fft_parity_ooura.cpp
    LABELS parity)
target_include_directories(tap_dsp_fft_parity_default_flags PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_compile_features(tap_dsp_fft_parity_default_flags PRIVATE cxx_std_20)
target_compile_definitions(tap_dsp_fft_parity_default_flags PRIVATE
    TAP_DSP_PARITY_INFORMATIONAL
    TAP_DSP_PARITY_MAX_N=${TAP_DSP_PARITY_MAX_N})
target_link_libraries(tap_dsp_fft_parity_default_flags PRIVATE
    tap_dsp_fft_reference_default
    tap_dsp_warnings)

# JUnit XML for the informational report (hosted only: bare metal has no file
# system). gtest reads GTEST_OUTPUT from the environment at start-up, and the
# helper exposes no per-test properties or extra arguments, so the ENVIRONMENT
# property is attached the way gtest_discover_tests itself attaches its
# registrations: a script on the directory's TEST_INCLUDE_FILES, which ctest
# runs after the discovery script that defines the test (if(TEST ...) is not
# available in that context, so the call is unconditional: a discovery failure
# then errors here as well as in its own _NOT_BUILT entry). Should the helper
# grow an EXTRA_ARGS/PROPERTIES pass-through, this collapses to one line there.
if(NOT TAP_DSP_BARE_METAL)
    set(_tap_dsp_parity_xml_script ${CMAKE_CURRENT_BINARY_DIR}/tap_dsp_fft_parity_default_flags_env.cmake)
    file(WRITE ${_tap_dsp_parity_xml_script}
        "set_tests_properties(fft_parity_ooura_default_flags.ReportsMaxUlpVersusOoura PROPERTIES\n"
        "    ENVIRONMENT \"GTEST_OUTPUT=xml:${CMAKE_BINARY_DIR}/parity-ulp.xml\")\n")
    set_property(DIRECTORY APPEND PROPERTY TEST_INCLUDE_FILES ${_tap_dsp_parity_xml_script})
endif()

# ==============================================================================
# END Stage 2a additions
# ==============================================================================

# ==============================================================================
# BEGIN Stage 3c addition (docs/audit-fft-and-code-smells.md, Part 3): the typed
# analysis instruments' battery. Its own target, like the parity gate above, so
# it can carry -ffp-contract=off: its pins are same-binary bit-identity
# comparisons between the templated instruments and verbatim copies of the
# pre-template ones, which GCC's cross-statement contraction can move by an
# ulp on FMA hardware (test_analysis_typed.cpp, file comment). The double
# joint-fit test is excluded on the QEMU legs for the same budget reason as
# MultitoneAnalysis in the main battery. Keep this block at the end of the file.
# ==============================================================================
tap_dsp_add_gtest_executable(tap_dsp_analysis_typed
    SOURCES test_analysis_typed.cpp
    MAIN_FILTER "*-MultitoneAnalysisTyped.*")
target_compile_options(tap_dsp_analysis_typed PRIVATE ${_tap_dsp_nocontract_cxx})
target_link_libraries(tap_dsp_analysis_typed PRIVATE
    tap::dsp
    tap_dsp_warnings)
# ==============================================================================
# END Stage 3c addition
# ==============================================================================

# ==============================================================================
# BEGIN Stage 6 addition: the tidy gate's view of bench/ and tools/.
# The style job (.github/workflows/style.yml, a taphouse file) configures the
# compile database with -DTAP_DSP_BUILD_TESTS=ON only, so the TUs behind
# TAP_DSP_BUILD_BENCH (bench/) and TAP_DSP_BUILD_CAPI (tools/capi/), and the
# standalone fingerprint tool (tools/fingerprint/), were never in it and never
# saw clang-tidy. Declared here as OBJECT libraries excluded from ALL: they
# land in compile_commands.json (which is all the tidy gate reads) and are
# never compiled by `cmake --build`, so the hosted battery's build time and
# everything the QEMU legs build are unchanged. Only where the real targets are
# absent (a configuration that builds them already lists them) and never on
# bare metal. icount_main.cpp is listed twice, floating (precision 0) and fixed
# point (precision 2), because its body is two preprocessor branches and
# clang-tidy checks every compile command a file has.
# =======================================================================# ==============================================================================

# ==============================================================================
if(NOT TAP_DSP_BARE_METAL)
    # Each entry is name:path-from-the-project-root:definition ("-" for none).
    # Relative paths, because an absolute Windows path's drive colon would split
    # the entry. The definition is the one the TU's real target passes: the
    # icount scenario's precision, or DSPTAP_BUILDING for the capi library.
    #
    # Ordering: this file is add_subdirectory'd before the root CMakeLists.txt
    # declares option(TAP_DSP_BUILD_BENCH) and option(TAP_DSP_BUILD_CAPI), so
    # the tests below read the cache (a -D on the command line) or, when neither
    # was given, an undefined variable, which if(NOT ...) treats as OFF, the
    # options' default. Moving those option() calls above add_subdirectory(tests)
    # would not change the result.
    set(_tap_dsp_tidy_only_tus)
    if(NOT TAP_DSP_BUILD_BENCH)
        list(APPEND _tap_dsp_tidy_only_tus
            bench_fft:bench/bench_fft.cpp:-
            size_probe:bench/size_probe.cpp:-
            icount_float:bench/icount/icount_main.cpp:TAP_DSP_SC_PRECISION=0
            icount_fixed:bench/icount/icount_main.cpp:TAP_DSP_SC_PRECISION=2)
    endif()
    if(NOT TAP_DSP_BUILD_CAPI)
        list(APPEND _tap_dsp_tidy_only_tus
            capi:tools/capi/dsptap_capi.cpp:DSPTAP_BUILDING
            capi_reference:tools/capi/test_capi_reference.cpp:-)
    endif()
    list(APPEND _tap_dsp_tidy_only_tus fingerprint:tools/fingerprint/fingerprint.cpp:-)
    foreach(_tu IN LISTS _tap_dsp_tidy_only_tus)
        string(REPLACE ":" ";" _parts "${_tu}")
        list(GET _parts 0 _name)
        list(GET _parts 1 _src)
        list(GET _parts 2 _definition)
        add_library(tap_dsp_tidy_only_${_name} OBJECT EXCLUDE_FROM_ALL ${PROJECT_SOURCE_DIR}/${_src})
        target_include_directories(tap_dsp_tidy_only_${_name} PRIVATE
            ${PROJECT_SOURCE_DIR}/bench ${PROJECT_SOURCE_DIR}/tools/capi)
        if(NOT _definition STREQUAL "-")
            target_compile_definitions(tap_dsp_tidy_only_${_name} PRIVATE ${_definition})
        endif()
        target_link_libraries(tap_dsp_tidy_only_${_name} PRIVATE tap::dsp tap_dsp_warnings)
    endforeach()
endif()
# ==============================================================================
# END Stage 6 addition
=======
# BEGIN Stage 4 addition (docs/audit-fft-and-code-smells.md, D4; audit F4): the
# two-image test of the ABI tag. One translation unit, tests/abi/abi_tag_image.cpp,
# built into two loadable modules — image A with the configured build's own
# defines (tap::dsp: fft_split_radix on linux, fft_vdsp on macOS), image B
# with TAP_DSP_FFT_CMSIS and tests/abi/stub/ first on the include path, so
# fft.h selects a stub engine with a different object layout (tag fft_cmsis)
# — and a gtest host (test_fft_abi_tag.cpp) that dlopens both with
# RTLD_GLOBAL and asks each image what its embedders see. Image B must NOT
# link tap::dsp: the target's INTERFACE definitions would add the host's own
# backend define (TAP_DSP_FFT_ACCELERATE on macOS) beside TAP_DSP_FFT_CMSIS
# and fft.h rejects the pair; it takes the include directory directly, like
# the parity targets. Hosted, non-Windows only: dlopen and MODULE libraries,
# and the coalescing under test has no path on Windows (nothing is exported
# without __declspec); the QEMU legs have no loader. Keep this block at the
# end of the file.
# ==============================================================================
if(NOT TAP_DSP_BARE_METAL AND NOT WIN32)
    add_library(tap_dsp_abi_image_a MODULE abi/abi_tag_image.cpp)
    target_link_libraries(tap_dsp_abi_image_a PRIVATE tap::dsp tap_dsp_warnings)

    add_library(tap_dsp_abi_image_b MODULE abi/abi_tag_image.cpp)
    target_include_directories(tap_dsp_abi_image_b BEFORE PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/abi/stub)
    target_include_directories(tap_dsp_abi_image_b PRIVATE ${PROJECT_SOURCE_DIR}/include)
    target_compile_features(tap_dsp_abi_image_b PRIVATE cxx_std_20)
    target_compile_definitions(tap_dsp_abi_image_b PRIVATE TAP_DSP_FFT_CMSIS)
    target_link_libraries(tap_dsp_abi_image_b PRIVATE tap_dsp_warnings)

    tap_dsp_add_gtest_executable(tap_dsp_fft_abi_tag
        SOURCES test_fft_abi_tag.cpp)
    add_dependencies(tap_dsp_fft_abi_tag tap_dsp_abi_image_a tap_dsp_abi_image_b)
    target_compile_definitions(tap_dsp_fft_abi_tag PRIVATE
        "TAP_DSP_ABI_IMAGE_A=\"$<TARGET_FILE:tap_dsp_abi_image_a>\""
        "TAP_DSP_ABI_IMAGE_B=\"$<TARGET_FILE:tap_dsp_abi_image_b>\"")
    target_link_libraries(tap_dsp_fft_abi_tag PRIVATE
        tap::dsp
        tap_dsp_warnings
        ${CMAKE_DL_LIBS})
endif()
# ==============================================================================
# END Stage 4 addition
# ==============================================================================
