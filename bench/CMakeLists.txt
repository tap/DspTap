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
# TAP_DSP_BENCH_ENGINE selects the engine everything here measures
# (bench_common.h): reference_c (the vendored C, the default) or split_radix
# (the Stage 2a port). From Stage 2a until 2c the float icount scenarios and
# the float size probe are also built against the port under a `_port`
# suffix, whatever this variable says, so every bench run reports both.
set(TAP_DSP_BENCH_ENGINE reference_c CACHE STRING
    "FFT engine the benchmark scenarios measure (reference_c or split_radix)")

if(CMAKE_SYSTEM_NAME STREQUAL "Generic")
    set(_tap_dsp_bench_host_class OFF)
else()
    set(_tap_dsp_bench_host_class ON)
endif()

if(_tap_dsp_bench_host_class)
    add_executable(tap_dsp_bench_fft bench_fft.cpp)
    target_compile_definitions(tap_dsp_bench_fft PRIVATE
        TAP_DSP_BENCH_ENGINE=${TAP_DSP_BENCH_ENGINE})
    target_include_directories(tap_dsp_bench_fft PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
    target_link_libraries(tap_dsp_bench_fft PRIVATE
        tap::dsp
        tap_dsp_warnings)
endif()

add_subdirectory(icount)

# Size probes: name:precision:n:engine, the same profiles as the icount
# scenarios at N = 512 (the ceilings are stated per profile, not per N). The
# float probe builds twice, C and port (`_port`), so bench.yml reports the
# port's .text beside the C's; the ceiling applies to the C probe only until
# Stage 2b routes the port (bench/README.md, "Sizes").
set(_tap_dsp_size_probes
    rfft_f32_512:0:512:${TAP_DSP_BENCH_ENGINE}
    rfft_f32_512_port:0:512:split_radix)
if(_tap_dsp_bench_host_class)
    list(APPEND _tap_dsp_size_probes rfft_f64_512:1:512:${TAP_DSP_BENCH_ENGINE})
endif()
foreach(_probe IN LISTS _tap_dsp_size_probes)
    string(REPLACE ":" ";" _parts "${_probe}")
    list(GET _parts 0 _name)
    list(GET _parts 1 _precision)
    list(GET _parts 2 _n)
    list(GET _parts 3 _engine)
    add_executable(tap_dsp_size_probe_${_name} size_probe.cpp)
    target_compile_definitions(tap_dsp_size_probe_${_name} PRIVATE
        TAP_DSP_SC_PRECISION=${_precision}
        TAP_DSP_SC_N=${_n}
        TAP_DSP_BENCH_ENGINE=${_engine})
    target_include_directories(tap_dsp_size_probe_${_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
    target_link_libraries(tap_dsp_size_probe_${_name} PRIVATE
        tap::dsp
        tap_dsp_warnings)
endforeach()
