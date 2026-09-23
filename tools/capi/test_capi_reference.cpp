/// @file test_capi_reference.cpp
/// @brief The C smoke test's independent reference: decimate.h called once over a whole stream.
// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// dsptap_decimator_process stages its input in fixed-size blocks, so a capi-vs-capi comparison
// cannot see a block-boundary bug. This helper is the other side: the header's
// basic_decimator<float, M>, constructed fresh and fed the whole stream in ONE process() call,
// with the same double -> float -> double conversions the capi documents. Compiled into
// dsptap_capi_smoke only (never into the shared library); the declaration lives in test_capi.c.

#include <cstddef>
#include <vector>

#include "tap/dsp/decimate.h"

namespace {

    template <std::size_t M>
    int reference(int transparent, const double* in, int n, double* out) {
        const auto p =
            transparent != 0 ? tap::dsp::decimate_profile::transparent() : tap::dsp::decimate_profile::economy();
        tap::dsp::basic_decimator<float, M> dec(p);
        const auto                          len = static_cast<std::size_t>(n);
        std::vector<float>                  fin(len);
        std::vector<float>                  fout(dec.outputs_for(len));
        for (std::size_t i = 0; i < len; ++i) {
            fin[i] = static_cast<float>(in[i]);
        }
        const std::size_t made = dec.process(fin.data(), len, fout.data());
        for (std::size_t i = 0; i < made; ++i) {
            out[i] = static_cast<double>(fout[i]);
        }
        return static_cast<int>(made);
    }

} // namespace

extern "C" int dsptap_test_reference_decimate(int ratio, int transparent, const double* in, int n, double* out) {
    switch (ratio) {
    case 2:
        return reference<2>(transparent, in, n, out);
    case 3:
        return reference<3>(transparent, in, n, out);
    case 6:
        return reference<6>(transparent, in, n, out);
    default:
        return -1;
    }
}
