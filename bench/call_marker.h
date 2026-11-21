/// @file call_marker.h
/// @brief The per-call markers the counting plugin brackets (bench/layout/).
// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// tools/qemu_insn_plugin/insn_count.c finds these two functions by symbol and
// reports the instructions executed between each call to tap_dsp_call_begin()
// and the next tap_dsp_call_end(): the most in one bracket and their sum
// (TAP_DSP_CALL_MAX). Bench harness code, never shipped; SampleRateTap's
// tools/probe/call_marker.h with the names changed.
//
// extern "C" so the symbols are the plain names the plugin looks for,
// noinline so each is a real call the guest executes, and different (each
// stores its own value) so no identical-code folding can make them one
// address. Each is a store and a return, so a bracket's count includes the
// end marker's few instructions, the same in every call.
#pragma once

#include <cstdint>

/// noinline in every compiler's spelling (MSVC rejects [[gnu::noinline]] under /WX).
#if defined(_MSC_VER) && !defined(__clang__)
#define TAP_DSP_BENCH_NOINLINE __declspec(noinline)
#else
#define TAP_DSP_BENCH_NOINLINE __attribute__((noinline))
#endif

namespace tap::dsp::bench::detail {
    inline volatile std::uint32_t g_call_marker = 0;
} // namespace tap::dsp::bench::detail

extern "C" {
// NOLINTBEGIN(readability-identifier-naming)
TAP_DSP_BENCH_NOINLINE inline void tap_dsp_call_begin() noexcept {
    tap::dsp::bench::detail::g_call_marker = 1;
}
TAP_DSP_BENCH_NOINLINE inline void tap_dsp_call_end() noexcept {
    tap::dsp::bench::detail::g_call_marker = 2;
}
// NOLINTEND(readability-identifier-naming)
}
