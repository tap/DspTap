# DspTap benchmarks (root option TAP_DSP_BUILD_BENCH; policy in README.md here).
#
#   tap_dsp_bench_fft   host wall-clock microbenchmark, min-of-N, informational
#                       only — never a gate. Host builds only (steady_clock).
#   icount/             deterministic instruction-count ratchet workloads, one
#                       binary per scenario; built on the host as a smoke test
#                       and on the QEMU cross legs, where
#                       .github/workflows/bench.yml counts them.
#
# TAP_DSP_BENCH_ENGINE selects the engine both measure (bench_common.h):
# reference_c today; split_radix from Stage 2a.
set(TAP_DSP_BENCH_ENGINE reference_c CACHE STRING
    "FFT engine the benchmark scenarios measure (reference_c; split_radix from Stage 2a)")

if(NOT CMAKE_CROSSCOMPILING)
    add_executable(tap_dsp_bench_fft bench_fft.cpp)
    target_compile_definitions(tap_dsp_bench_fft PRIVATE
        TAP_DSP_BENCH_ENGINE=${TAP_DSP_BENCH_ENGINE})
    target_include_directories(tap_dsp_bench_fft PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
    target_link_libraries(tap_dsp_bench_fft PRIVATE
        tap::dsp
        tap_dsp_warnings)
endif()

add_subdirectory(icount)
