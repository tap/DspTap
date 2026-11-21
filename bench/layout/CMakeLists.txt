# SPDX-License-Identifier: MIT
# Copyright 2026 Timothy Place and the DspTap contributors.
#
# Interleaved against planar (bench/README.md, "Layout"): informational
# instruction counts, never a gate. The binaries carry their own prefix,
# tap_dsp_layout_, so the ratchet's glob (tap_dsp_icount_*) never sees them;
# bench.yml counts them with scripts/icount.py --prefix tap_dsp_layout_
# --report-only, which prints each count and each planar / interleaved
# ratio and fails only when a binary yields no count.
#
# Scenarios: <kind>_f32_ch<C>_<i|p>, kind in {bw, fir, chain, gain}, C in
# {1, 2, 6, 12}; at one channel the two layouts are the same buffer, so only
# _i is built.

set(_tap_dsp_layout_kinds bw:0 fir:1 chain:2 gain:3)
set(_tap_dsp_layout_channels 1 2 6 12)

foreach(_kind IN LISTS _tap_dsp_layout_kinds)
    string(REPLACE ":" ";" _parts "${_kind}")
    list(GET _parts 0 _kname)
    list(GET _parts 1 _kid)
    foreach(_ch IN LISTS _tap_dsp_layout_channels)
        set(_layouts i:0)
        if(_ch GREATER 1)
            list(APPEND _layouts p:1)
        endif()
        foreach(_layout IN LISTS _layouts)
            string(REPLACE ":" ";" _lparts "${_layout}")
            list(GET _lparts 0 _lname)
            list(GET _lparts 1 _planar)
            set(_key ${_kname}_f32_ch${_ch}_${_lname})
            add_executable(tap_dsp_layout_${_key} layout_main.cpp)
            target_compile_definitions(tap_dsp_layout_${_key} PRIVATE
                TAP_DSP_LY_NAME="${_key}"
                TAP_DSP_LY_KIND=${_kid}
                TAP_DSP_LY_CH=${_ch}
                TAP_DSP_LY_PLANAR=${_planar})
            target_include_directories(tap_dsp_layout_${_key} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/..)
            target_link_libraries(tap_dsp_layout_${_key} PRIVATE
                tap::dsp
                tap_dsp_warnings)
            if(NOT CMAKE_SYSTEM_NAME STREQUAL "Generic")
                add_test(NAME tap_dsp_layout_${_key}_smoke
                    COMMAND ${CMAKE_COMMAND}
                        -DBINARY=$<TARGET_FILE:tap_dsp_layout_${_key}>
                        -P ${CMAKE_CURRENT_SOURCE_DIR}/../icount/smoke.cmake)
            endif()
        endforeach()
        # The two layouts of a pair must print the same checksum: the same
        # output, bit for bit.
        if(_ch GREATER 1 AND NOT CMAKE_SYSTEM_NAME STREQUAL "Generic")
            set(_pair ${_kname}_f32_ch${_ch})
            add_test(NAME tap_dsp_layout_${_pair}_pair_agrees
                COMMAND ${CMAKE_COMMAND}
                    -DFIRST=$<TARGET_FILE:tap_dsp_layout_${_pair}_i>
                    -DSECOND=$<TARGET_FILE:tap_dsp_layout_${_pair}_p>
                    -P ${CMAKE_CURRENT_SOURCE_DIR}/pair.cmake)
        endif()
    endforeach()
endforeach()
