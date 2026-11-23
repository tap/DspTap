# Instruction-count ratchet workloads (bench/README.md). One binary per
# scenario: bare-metal targets have no argv, and per-binary instruction totals
# are what scripts/icount.py compares. The scenario name (after the
# tap_dsp_icount_ prefix) is the key in bench/baselines.json.
#
# Entries are name:precision:n — precision 0 = float, 1 = double, 2 = Q15,
# 3 = Q31 (see icount_main.cpp). The float pair is Stage 1b's; the three
# fixed-point scenarios are the ones Stage 3b specified (audit Part 11),
# seeded at Stage 2c, and run on every leg: fixed point is a profile on the
# soft-float M4 too. Adapted from MuTap's bench/icount/CMakeLists.txt.
set(_tap_dsp_icount_scenarios
    rfft_f32_512:0:512
    rfft_f32_2048:0:2048
    rfft_q15_512:2:512
    rfft_q31_512:3:512
    rfft_q31_2048:3:2048)

# The double scenario is host-class only: on the bare-metal Cortex-M legs
# double is soft-float, which is a correctness check elsewhere in the repo but
# not a profile, so it is neither built nor baselined there.
if(NOT CMAKE_SYSTEM_NAME STREQUAL "Generic")
    list(APPEND _tap_dsp_icount_scenarios rfft_f64_512:1:512)
endif()

# One scenario binary: tap_dsp_icount_<key>, measuring basic_real_fft as
# built (bench_common.h) over precision/n.
function(tap_dsp_add_icount_scenario _key _precision _n)
    add_executable(tap_dsp_icount_${_key} icount_main.cpp)
    target_compile_definitions(tap_dsp_icount_${_key} PRIVATE
        TAP_DSP_SC_NAME="${_key}"
        TAP_DSP_SC_PRECISION=${_precision}
        TAP_DSP_SC_N=${_n})
    target_include_directories(tap_dsp_icount_${_key} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/..)
    target_link_libraries(tap_dsp_icount_${_key} PRIVATE
        tap::dsp
        tap_dsp_warnings)
    # The QEMU legs link the semihosting startup through the toolchain file's
    # linker flags exactly as the test binary does; nothing extra here.

    # Host smoke test, the README's rule made executable: the binary runs
    # twice, both runs report ok=1 and print the same DONE line. There is no
    # instruction count on the host; that is the QEMU job's.
    if(NOT CMAKE_SYSTEM_NAME STREQUAL "Generic")
        add_test(NAME tap_dsp_icount_${_key}_smoke
            COMMAND ${CMAKE_COMMAND}
                -DBINARY=$<TARGET_FILE:tap_dsp_icount_${_key}>
                -P ${CMAKE_CURRENT_SOURCE_DIR}/smoke.cmake)
    endif()
endfunction()

# Every scenario is a gated binary keyed by its bare name. (From Stage 2a
# until 2c a second loop here built each float scenario again against the
# vendored Ooura C under an informational `_port` / `_c` suffix, so the job
# could print the C/port ratio and whether the two checksums agreed; Stage 2c
# retired the C from the shipping tree and the pair with it.)
foreach(_sc IN LISTS _tap_dsp_icount_scenarios)
    string(REPLACE ":" ";" _parts "${_sc}")
    list(GET _parts 0 _name)
    list(GET _parts 1 _precision)
    list(GET _parts 2 _n)
    tap_dsp_add_icount_scenario(${_name} ${_precision} ${_n})
endforeach()

# The Butterworth stage's scenarios (stage_main.cpp): name:kind:order:ch:call.
# Their counts are bracket sums, not whole runs (the file's header); the DONE
# line says so and scripts/icount.py reads it.
set(_tap_dsp_stage_scenarios
    bw4_f32_ch2:0:4:2:128
    bw8_f32_ch2:0:8:2:128
    bw8_f32_ch1:0:8:1:128
    bw4_f32_ch2_n4:0:4:2:4)

# set_cutoff() runs the design, which is double. Where the target computes
# double in hardware (the hexagon key, the M55 keys, a host) the two
# scenarios are gated; where double is libgcc's soft float (both M4 keys, the
# M33) they would gate the toolchain's library, the reason rfft_f64 is not
# built there, so they are built under the informational prefix instead and
# bench.yml reports them without a baseline.
include(CheckCXXSourceCompiles)
check_cxx_source_compiles("
#if defined(__arm__) && !(defined(__ARM_FP) && (__ARM_FP & 8))
#error soft double
#endif
int main() { return 0; }" TAP_DSP_BENCH_HARD_DOUBLE)
set(_tap_dsp_stage_design_scenarios
    bw8_f32_set_cutoff:1:8:2:128
    bw8_f32_set_cutoff_deviation:2:8:2:128)

function(tap_dsp_add_stage_scenario _prefix _key _kind _order _ch _call)
    add_executable(${_prefix}${_key} stage_main.cpp)
    target_compile_definitions(${_prefix}${_key} PRIVATE
        TAP_DSP_ST_NAME="${_key}"
        TAP_DSP_ST_KIND=${_kind}
        TAP_DSP_ST_ORDER=${_order}
        TAP_DSP_ST_CH=${_ch}
        TAP_DSP_ST_CALL=${_call})
    target_include_directories(${_prefix}${_key} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/..)
    target_link_libraries(${_prefix}${_key} PRIVATE
        tap::dsp
        tap_dsp_warnings)
    if(NOT CMAKE_SYSTEM_NAME STREQUAL "Generic")
        add_test(NAME ${_prefix}${_key}_smoke
            COMMAND ${CMAKE_COMMAND}
                -DBINARY=$<TARGET_FILE:${_prefix}${_key}>
                -P ${CMAKE_CURRENT_SOURCE_DIR}/smoke.cmake)
    endif()
endfunction()

foreach(_sc IN LISTS _tap_dsp_stage_scenarios _tap_dsp_stage_design_scenarios)
    string(REPLACE ":" ";" _parts "${_sc}")
    list(GET _parts 0 _name)
    list(GET _parts 1 _kind)
    list(GET _parts 2 _order)
    list(GET _parts 3 _ch)
    list(GET _parts 4 _call)
    set(_prefix tap_dsp_icount_)
    if(NOT _kind EQUAL 0 AND NOT TAP_DSP_BENCH_HARD_DOUBLE)
        set(_prefix tap_dsp_info_)
    endif()
    tap_dsp_add_stage_scenario(${_prefix} ${_name} ${_kind} ${_order} ${_ch} ${_call})
endforeach()
