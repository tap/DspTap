include(FetchContent)

# The hosted test build exercises the default Ooura path and the Apple vDSP
# backend on macOS. The four bare-metal legs (cmake/arm-cortex-m4-mps2.cmake
# in both flavours, arm-cortex-m33-mps2.cmake, arm-cortex-m55-mps3.cmake) run
# the emulation-sized selection in tests/bare_metal_main.cpp under
# qemu-system-arm; the CMSIS-DSP Helium backend and its parity suite run for
# real on the M55 leg.
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

add_executable(tap_dsp_tests
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
    test_yin.cpp)
target_link_libraries(tap_dsp_tests PRIVATE
    tap::dsp
    tap_dsp_warnings)

if(TAP_DSP_BARE_METAL)
    # One-shot mode (Cortex-M under qemu-system-arm): a custom main bakes in
    # the emulation-appropriate test selection (no argv on the target), and
    # the run is judged on gtest's summary text because semihosting does not
    # reliably propagate exit codes through the emulator. CTest prepends
    # CMAKE_CROSSCOMPILING_EMULATOR (set by the toolchain file) to the
    # command automatically.
    target_sources(tap_dsp_tests PRIVATE bare_metal_main.cpp)
    target_link_libraries(tap_dsp_tests PRIVATE GTest::gtest)
    add_test(NAME tap_dsp_tests_emulated COMMAND tap_dsp_tests)
    set_tests_properties(tap_dsp_tests_emulated PROPERTIES
        PASS_REGULAR_EXPRESSION "TAP_DSP_TESTS_COMPLETE rc=0"
        FAIL_REGULAR_EXPRESSION "\\[  FAILED  \\]"
        TIMEOUT 1800)
else()
    target_link_libraries(tap_dsp_tests PRIVATE GTest::gtest_main)
    include(GoogleTest)
    gtest_discover_tests(tap_dsp_tests DISCOVERY_TIMEOUT 120 PROPERTIES TIMEOUT 900)
endif()
