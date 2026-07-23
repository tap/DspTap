include(FetchContent)

# The hosted test build (this file) exercises the default Ooura path and the
# Apple vDSP backend on macOS. The CMSIS-DSP Helium backend is compile-verified
# under the bare-metal M55 toolchain (cmake/arm-cortex-m55-mps3.cmake); running
# the parity suite under qemu is left to the consuming library's harness.
find_package(Threads QUIET)
if(NOT Threads_FOUND)
    set(gtest_disable_pthreads ON CACHE BOOL "" FORCE)
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
    test_fft.cpp
    test_fft_backend.cpp
    test_fir_kernels.cpp
    test_kaiser.cpp
    test_psola.cpp
    test_pvoc.cpp
    test_quantize.cpp
    test_sample_traits.cpp
    test_yin.cpp)
target_link_libraries(tap_dsp_tests PRIVATE
    tap::dsp
    tap_dsp_warnings)

target_link_libraries(tap_dsp_tests PRIVATE GTest::gtest_main)
include(GoogleTest)
gtest_discover_tests(tap_dsp_tests DISCOVERY_TIMEOUT 120 PROPERTIES TIMEOUT 900)
