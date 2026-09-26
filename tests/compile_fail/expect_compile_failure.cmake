# Runs one build target that must FAIL to compile, and checks how it fails.
#
#   cmake -DBUILD_DIR=<dir> -DTARGET=<target> -DCONFIG=<config>
#         -DEXPECT=<regex> -P expect_compile_failure.cmake
#
# Passes only if (1) the build fails, (2) its output matches EXPECT, and
# (3) the output has exactly one compiler "error:" line — the expected one.
# (2) is what keeps the test from passing vacuously (a missing include or a
# typo also fails the build, with a different message); (3) pins that the
# expected diagnostic is the only one (no cascade after it). GCC and Clang
# print one "error:" per error; `make` and `ninja` report the failed step
# as "Error 1" / "FAILED:", which the count does not see.
foreach(var BUILD_DIR TARGET EXPECT)
    if(NOT DEFINED ${var})
        message(FATAL_ERROR "expect_compile_failure.cmake: ${var} is required")
    endif()
endforeach()
set(config_args)
if(CONFIG)
    set(config_args --config ${CONFIG})
endif()
execute_process(
    COMMAND ${CMAKE_COMMAND} --build ${BUILD_DIR} --target ${TARGET} ${config_args}
    RESULT_VARIABLE result
    OUTPUT_VARIABLE out
    ERROR_VARIABLE out)
message("${out}")
if(result EQUAL 0)
    message(FATAL_ERROR "COMPILE-FAIL CHECK: ${TARGET} compiled; it must not")
endif()
if(NOT out MATCHES "${EXPECT}")
    message(FATAL_ERROR "COMPILE-FAIL CHECK: ${TARGET} failed without the expected diagnostic /${EXPECT}/")
endif()
string(REGEX MATCHALL "error:" errors "${out}")
list(LENGTH errors n_errors)
if(NOT n_errors EQUAL 1)
    message(FATAL_ERROR "COMPILE-FAIL CHECK: ${TARGET} produced ${n_errors} errors; exactly 1 (the expected one) is required")
endif()
message("COMPILE-FAIL CHECK PASSED: ${TARGET} rejected with /${EXPECT}/ as its only error")
