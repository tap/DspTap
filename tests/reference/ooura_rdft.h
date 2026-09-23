/// @file ooura_rdft.h
/// @brief Declarations of the reference C (Takuya Ooura's rdft) for the parity gate.
// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// The two entry points of the reference copy of Ooura's split-radix real FFT
// kept under tests/reference/ooura/: rdft from fftsg.c (double) and rdft_f
// from fftsg_float.c (the same source instantiated for float). This header is
// DspTap's own (MIT, house-formatted) and sits one level above that directory
// on purpose: the directory's .clang-format disables formatting for the two
// vendored files alone, and the house hook keeps checking this one.
//
// Until Stage 2c these declarations lived in include/tap/dsp/fft.h, where the
// shipping wrapper called them; since the flip (Stage 2b) nothing that ships
// reaches them, and Stage 2c moved the C out of the shipping tree
// (docs/audit-fft-and-code-smells.md, Part 3; Decision D6). The only caller
// is tests/test_fft_parity_ooura.cpp, the bit-identity gate for the C++20
// engine that replaced the C, linked against the private reference builds
// tests/CMakeLists.txt makes of these two files.
//
// Calling convention, as fftsg.c documents it: n is the transform length
// (power of two), isgn +1 forward / -1 inverse (unnormalized), a the n
// samples in place, ip an int workspace of 2 + sqrt(n/2) entries with
// ip[0] = 0 requesting table initialization on the first call, w a sample
// workspace of n/2.
#pragma once

extern "C" {
void rdft(int n, int isgn, double* a, int* ip, double* w);
void rdft_f(int n, int isgn, float* a, int* ip, float* w);
}
