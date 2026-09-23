# DspTap benchmarks (root option TAP_DSP_BUILD_BENCH; policy in README.md here).
#
#   tap_dsp_bench_fft        host wall-clock microbenchmark, min-of-N, informational
#                            only — never a gate. Host-class builds only.
#   icount/                  deterministic instruction-count ratchet workloads, one
#                            binary per scenario; built on the host as a smoke test
#                            (ctest runs each twice and compares) and on the QEMU
#                            cross legs, where .github/workflows/bench.yml counts them.
#   tap_dsp_size_probe_*     one profile alone, for the per-leg .text ceilings
#                            (built MinSizeRel by bench.yml, measured with size -A).
#
# "Host-class" is decided by CMAKE_SYSTEM_NAME STREQUAL "Generic" (bare metal)
# throughout, never by CMAKE_CROSSCOMPILING: a Linux-arm64 cross build is
# host-class and gets the double scenario; the Cortex-M legs are not.
#
# Everything here measures tap::dsp::basic_real_fft as built (bench_common.h):
# the split-radix engine, or CMSIS / vDSP behind the class for float where the
# build says so, and the int32 kernel for the fixed-point profiles. From
# Stage 2b until 2c a TAP_DSP_BENCH_ENGINE cache variable could also point the
# bench at the vendored C, and the float scenarios and probe were built a
# second time against it under a `_c` suffix; Stage 2c retired the C from the
# shipping tree and that mechanism with it (README.md).

if(CMAKE_SYSTEM_NAME STREQUAL "Generic")
    set(_tap_dsp_bench_host_class OFF)
else()
    set(_tap_dsp_bench_host_class ON)
endif()

if(_tap_dsp_bench_host_class)
    add_executable(tap_dsp_bench_fft bench_fft.cpp)
    target_include_directories(tap_dsp_bench_fft PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
    target_link_libraries(tap_dsp_bench_fft PRIVATE
        tap::dsp
        tap_dsp_warnings)
endif()

add_subdirectory(icount)

# Size probes: name:precision:n, the same profiles as the icount scenarios at
# N = 512 (the ceilings are stated per profile, not per N); precision codes as
# in icount_main.cpp (0 float, 1 double, 2 Q15, 3 Q31). The float probe is the
# Part 4 `.text` assertion for the float instantiation; the two fixed-point
# probes are Stage 3b's. Each key in bench.yml carries a ceiling per probe
# (bench/README.md, "Sizes").
set(_tap_dsp_size_probes
    rfft_f32_512:0:512
    rfft_q15_512:2:512
    rfft_q31_512:3:512)
if(_tap_dsp_bench_host_class)
    list(APPEND _tap_dsp_size_probes rfft_f64_512:1:512)
endif()
foreach(_probe IN LISTS _tap_dsp_size_probes)
    string(REPLACE ":" ";" _parts "${_probe}")
    list(GET _parts 0 _name)
    list(GET _parts 1 _precision)
    list(GET _parts 2 _n)
    add_executable(tap_dsp_size_probe_${_name} size_probe.cpp)
    target_compile_definitions(tap_dsp_size_probe_${_name} PRIVATE
        TAP_DSP_SC_PRECISION=${_precision}
        TAP_DSP_SC_N=${_n})
    target_include_directories(tap_dsp_size_probe_${_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
    target_link_libraries(tap_dsp_size_probe_${_name} PRIVATE
        tap::dsp
        tap_dsp_warnings)
endforeach()
