/// @file version.h
/// @brief DspTap's version, for consumers to assert the minimum they were tested with.
// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// One header for the one tree. The version is the CMake project()'s, kept
// equal by tests/test_version.cpp, and a tag vX.Y.Z marks the commit that
// carries it. DspTap is tagged on its own cadence, never in lockstep with
// a consumer: a consumer pins a commit (a submodule gitlink is a SHA), names
// the tag in the commit that moves the pin, and states the oldest DspTap it
// was tested with as a compile-time floor:
//
//   #include "tap/dsp/version.h"
//   static_assert(TAP_DSP_VERSION >= TAP_DSP_VERSION_ENCODE(0, 1, 0),
//                 "this library needs DspTap 0.1.0 or newer");
//
// The packing is the Tap family's C ABI version probe's,
// (major << 16) | (minor << 8) | patch, so versions compare as integers in
// #if and static_assert alike. Before 1.0 a minor bump may break a
// consumer; a floor is a statement of what was tested, not a promise of
// what will keep working.
#pragma once

#define TAP_DSP_VERSION_MAJOR 0
#define TAP_DSP_VERSION_MINOR 1
#define TAP_DSP_VERSION_PATCH 0

/// (major << 16) | (minor << 8) | patch: comparable as an integer.
#define TAP_DSP_VERSION_ENCODE(major, minor, patch) (((major) << 16) | ((minor) << 8) | (patch))

/// This tree's version, encoded.
#define TAP_DSP_VERSION TAP_DSP_VERSION_ENCODE(TAP_DSP_VERSION_MAJOR, TAP_DSP_VERSION_MINOR, TAP_DSP_VERSION_PATCH)
