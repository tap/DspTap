/// @file pitch.h
/// @brief The pitch-shifter batteries' shared measurement helpers.
// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// run_sine, measure_hz, tail_peak and cents were copy-pasted verbatim between
// test_psola.cpp and test_pvoc.cpp (docs/audit-fft-and-code-smells.md,
// Part 2); they live here since Stage 6, with the same arithmetic, so both
// batteries' pinned numbers are unchanged. The pitch oracle is
// tap::dsp::yin (double), certified by its own battery.

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include <gtest/gtest.h>

#include "signals.h"
#include "tap/dsp/yin.h"

namespace tap::dsp::test {

    /// Feed `seconds` of a unit sine at freq_hz (sine_at, sampled at
    /// sample_rate) through a shifter at a fixed ratio, one sample at a time,
    /// and return its output. A shifter whose process() takes the source period
    /// (psola) is given sample_rate / freq_hz, rounded to its Sample.
    template <template <typename> class Shifter, typename Sample>
    std::vector<Sample> run_sine(Shifter<Sample>& shifter, double freq_hz, double ratio, double seconds,
                                 double sample_rate) {
        const int           n      = static_cast<int>(seconds * sample_rate);
        const Sample        period = static_cast<Sample>(sample_rate / freq_hz);
        std::vector<Sample> out(static_cast<std::size_t>(n));
        for (int t = 0; t < n; ++t) {
            const Sample x = static_cast<Sample>(sine_at(freq_hz, static_cast<double>(t), sample_rate));
            if constexpr (requires { shifter.process(x, period, x); }) {
                out[static_cast<std::size_t>(t)] = shifter.process(x, period, static_cast<Sample>(ratio));
            }
            else {
                out[static_cast<std::size_t>(t)] = shifter.process(x, static_cast<Sample>(ratio));
            }
        }
        return out;
    }

    /// The fundamental of the signal's tail in Hz, by tap::dsp::yin (double)
    /// over one frame of the last samples: tau in [sample_rate / 2000,
    /// ceil(sample_rate / 55)], window = tau_max. Expects (gtest) a voiced
    /// verdict; returns 0 when unvoiced.
    template <typename Sample>
    double measure_hz(const std::vector<Sample>& x, double sample_rate) {
        const auto          tau_min = static_cast<std::size_t>(sample_rate / 2000.0);
        const auto          tau_max = static_cast<std::size_t>(std::ceil(sample_rate / 55.0));
        tap::dsp::yin       det(tau_max, tau_min, tau_max);
        std::vector<double> tail(det.frame_size());
        for (std::size_t i = 0; i < tail.size(); ++i) {
            tail[i] = static_cast<double>(x[x.size() - tail.size() + i]);
        }
        const auto r = det.analyze(tail.data());
        EXPECT_TRUE(r.voiced());
        return (r.period > 0.0) ? sample_rate / r.period : 0.0;
    }

    /// Peak |x| over the last `samples` samples.
    template <typename Sample>
    double tail_peak(const std::vector<Sample>& x, std::size_t samples = 4800) {
        double peak = 0.0;
        for (std::size_t i = x.size() - samples; i < x.size(); ++i) {
            peak = std::max(peak, std::abs(static_cast<double>(x[i])));
        }
        return peak;
    }

    /// Interval from ref to f in cents: 1200 log2(f / ref).
    inline double cents(double f, double ref) {
        return 1200.0 * std::log2(f / ref);
    }

} // namespace tap::dsp::test
