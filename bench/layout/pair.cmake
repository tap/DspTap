# SPDX-License-Identifier: MIT
# Copyright 2026 Timothy Place and the DspTap contributors.
#
# ctest -P script: the interleaved and the planar binary of one layout pair
# print the same checksum (bench/layout/layout_main.cpp folds each channel in
# frame order, so equal checksums mean bit-identical output).
foreach(_var FIRST SECOND)
    if(NOT ${_var})
        message(FATAL_ERROR "pair.cmake: -D${_var}=<binary> is required")
    endif()
endforeach()
set(_sums)
foreach(_bin "${FIRST}" "${SECOND}")
    execute_process(COMMAND "${_bin}" RESULT_VARIABLE _rc OUTPUT_VARIABLE _out)
    if(NOT _rc EQUAL 0)
        message(FATAL_ERROR "${_bin}: exited ${_rc}:\n${_out}")
    endif()
    string(REGEX MATCH "checksum=0x[0-9a-f]+" _sum "${_out}")
    if(NOT _sum)
        message(FATAL_ERROR "${_bin}: no checksum:\n${_out}")
    endif()
    list(APPEND _sums "${_sum}")
endforeach()
list(GET _sums 0 _a)
list(GET _sums 1 _b)
if(NOT _a STREQUAL _b)
    message(FATAL_ERROR "layouts disagree: ${FIRST} ${_a}, ${SECOND} ${_b}")
endif()
message(STATUS "both layouts: ${_a}")
