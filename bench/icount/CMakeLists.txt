# Instruction-count ratchet workloads (bench/README.md). One binary per
# scenario: bare-metal targets have no argv, and per-binary instruction totals
# are what scripts/icount.py compares. The scenario name (after the
# tap_dsp_icount_ prefix) is the key in bench/baselines.json.
#
# Entries are name:precision:n — precision 0 = float, 1 = double (see
# icount_main.cpp). Adapted from MuTap's bench/icount/CMakeLists.txt.
set(_tap_dsp_icount_scenarios
    rfft_f32_512:0:512
    rfft_f32_2048:0:2048)

# The double scenario is host-class only: on the bare-metal Cortex-M legs
# double is soft-float, which is a correctness check elsewhere in the repo but
# not a profile, so it is neither built nor baselined there.
if(NOT CMAKE_SYSTEM_NAME STREQUAL "Generic")
    list(APPEND _tap_dsp_icount_scenarios rfft_f64_512:1:512)
endif()

# TODO(Stage 2a): when TAP_DSP_BENCH_ENGINE gains split_radix, build each float
# scenario twice (the C and the port, e.g. a `_port` suffix on the second
# binary's key) so the job prints both counts and their ratio.
foreach(_sc IN LISTS _tap_dsp_icount_scenarios)
    string(REPLACE ":" ";" _parts "${_sc}")
    list(GET _parts 0 _name)
    list(GET _parts 1 _precision)
    list(GET _parts 2 _n)
    add_executable(tap_dsp_icount_${_name} icount_main.cpp)
    target_compile_definitions(tap_dsp_icount_${_name} PRIVATE
        TAP_DSP_SC_NAME="${_name}"
        TAP_DSP_SC_PRECISION=${_precision}
        TAP_DSP_SC_N=${_n}
        TAP_DSP_BENCH_ENGINE=${TAP_DSP_BENCH_ENGINE})
    target_include_directories(tap_dsp_icount_${_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/..)
    target_link_libraries(tap_dsp_icount_${_name} PRIVATE
        tap::dsp
        tap_dsp_warnings)
    # The QEMU legs link the semihosting startup through the toolchain file's
    # linker flags exactly as the test binary does; nothing extra here.
endforeach()
