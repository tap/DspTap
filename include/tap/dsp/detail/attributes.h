/// @file attributes.h
/// @brief TAP_DSP_NOINLINE — keep a setup-time function out of its caller's body.
// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// TAP_DSP_NOINLINE marks a function the compiler must not inline: GCC and
// Clang's [[gnu::noinline]], MSVC's __declspec(noinline), nothing elsewhere.
// It is for setup-time code whose inlining would change the code AROUND it
// rather than itself — the accelerated FFT engines' constructors
// (fft/backends/cmsis.h, accelerate.h), which allocate scratch, call into a
// vendor library and, on Apple, may throw. Inlined into a function that
// also runs a processing loop, such a constructor's live ranges and cleanup
// paths take part in that function's register allocation. Measured on the
// Cortex-M55 ratchet key (bench/README.md, Stage 4): with the CMSIS
// constructor inlined, the bench's checksum fold spilled the low word of
// its 64-bit hash and rematerialized the FNV prime on every element, +12 %
// instructions over the whole scenario with the transform byte-identical
// and its own code unchanged instruction for instruction; out of line, the
// count is the recorded baseline to within +5 / +45 instructions. A
// consumer's own processing function has the same exposure whenever it
// constructs an engine in the function it transforms in.
//
// Not a performance knob for hot code and not used on any transform; the
// split-radix engine's constructor (the table build) is left to the
// compiler, as the recorded baselines already include whatever it did with
// it.

#pragma once

#if defined(__GNUC__) || defined(__clang__)
#define TAP_DSP_NOINLINE [[gnu::noinline]]
#elif defined(_MSC_VER)
#define TAP_DSP_NOINLINE __declspec(noinline)
#else
#define TAP_DSP_NOINLINE
#endif
