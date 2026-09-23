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

# One scenario binary: tap_dsp_icount_<key>, measuring <engine>
# (bench_common.h: basic_real_fft or reference_c) over precision/n.
function(tap_dsp_add_icount_scenario _key _precision _n _engine)
    add_executable(tap_dsp_icount_${_key} icount_main.cpp)
    target_compile_definitions(tap_dsp_icount_${_key} PRIVATE
        TAP_DSP_SC_NAME="${_key}"
        TAP_DSP_SC_PRECISION=${_precision}
        TAP_DSP_SC_N=${_n}
        TAP_DSP_BENCH_ENGINE=${_engine})
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

foreach(_sc IN LISTS _tap_dsp_icount_scenarios)
    string(REPLACE ":" ";" _parts "${_sc}")
    list(GET _parts 0 _name)
    list(GET _parts 1 _precision)
    list(GET _parts 2 _n)
    # The gated binary: the engine the cache variable selects (basic_real_fft,
    # what ships, unless overridden), keyed by the bare scenario name.
    tap_dsp_add_icount_scenario(${_name} ${_precision} ${_n} ${TAP_DSP_BENCH_ENGINE})
    # Stage 2b until 2c: every float scenario also builds against the vendored
    # C called directly (TAP_DSP_BENCH_ENGINE=reference_c), key suffix `_c`.
    # The suffix is scripts/icount.py's INFORMATIONAL_SUFFIX: the binary is
    # counted and printed with its ratio to the gated sibling and whether the
    # two output checksums agree, never gated and never written to
    # bench/baselines.json. Until 2b the same mechanism ran the other way
    # (`_port`: the port beside the C while the C shipped); the flip
    # re-recorded the bare keys to the port and the C became the informational
    # side. Stage 2c deletes the C and this pair with it.
    if(_precision EQUAL 0 AND NOT TAP_DSP_BENCH_ENGINE STREQUAL "reference_c")
        tap_dsp_add_icount_scenario(${_name}_c ${_precision} ${_n} reference_c)
    endif()
endforeach()
