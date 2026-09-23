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
# (bench_common.h): basic_real_fft (what ships: the split-radix engine, or
# CMSIS / vDSP behind the class for float where the build says so; the
# default) or reference_c (the vendored C called directly). From Stage 2b
# until 2c the float icount scenarios and the float size probe are also built
# against the C under a `_c` suffix, so every bench run at the default reports
# both — unless this variable is already reference_c, in which case the
# bare-keyed binaries are the C and the `_c` pairs are not built (same rule in
# icount/).
set(TAP_DSP_BENCH_ENGINE basic_real_fft CACHE STRING
    "FFT engine the benchmark scenarios measure (basic_real_fft or reference_c)")

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
# float probe builds twice, what ships and the vendored C (`_c`), so bench.yml
# reports the C's .text beside the shipping probe's until Stage 2c retires
# the C; the ceiling applies to the bare probe (bench/README.md, "Sizes").
# The `_c` probe is skipped when the engine is already reference_c, exactly
# as icount/ skips the `_c` scenarios: the probe would be the C compared to
# itself.
set(_tap_dsp_size_probes rfft_f32_512:0:512:${TAP_DSP_BENCH_ENGINE})
if(NOT TAP_DSP_BENCH_ENGINE STREQUAL "reference_c")
    list(APPEND _tap_dsp_size_probes rfft_f32_512_c:0:512:reference_c)
endif()
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
