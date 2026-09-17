// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// Same-host A/B fingerprint of the primitives whose output must not move when
// their internals do. Runs tap::dsp::basic_pvoc (float and double, ratio 1.0
// and 1.5, formant off and on; N = 1024) and tap::dsp::basic_log_mel (float
// and double, plain log and PCEN; the default geometry) over one fixed corpus
// (xorshift32 noise plus two tones, 48000 samples) and prints the FNV-1a-64
// hash of each case's raw output bytes, one line per case.
//
// HOW TO READ THE NUMBERS. A hash is a function of the source tree AND of the
// host, the compiler, its flags and its libm: pvoc runs atan2/sin/cos in
// double, and floating-point contraction alone changes every pvoc line. So the
// hashes are never golden values, never committed, never asserted in a test
// (a pinned hash would fail on two of the three CI hosts by design). They are
// diffed: build this tool against the tree before a change and the tree after
// it, on the same host with the same compiler and flags, and the twelve lines
// must be identical when the change claims bit identity (the gate for every
// migration of pvoc.h or log_mel.h: the packed-spectrum view, the FFT port's
// routing flip, the engine parameter, the Hann/pi consolidation).
//
//     cmake -B build_fp -S tools/fingerprint -DCMAKE_BUILD_TYPE=Release
//     cmake --build build_fp && build_fp/dsptap_fingerprint > before.txt
//     ... apply the change, rebuild ...
//     build_fp/dsptap_fingerprint > after.txt && diff before.txt after.txt
//
// Every case is a chain of noexcept process() calls, so an assert-enabled
// (Debug) build of the same tree exercises the preconditions as well; its
// hashes are comparable only with another Debug build.

#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <numbers>
#include <vector>

#include "tap/dsp/log_mel.h"
#include "tap/dsp/pvoc.h"

namespace {

    constexpr double      k_sample_rate = 48000.0;
    constexpr std::size_t k_samples     = 48000;
    constexpr std::size_t k_pvoc_fft    = 1024;

    class xorshift32 {
      public:
        explicit xorshift32(std::uint32_t seed)
            : m_s(seed) {}
        double next() noexcept {
            m_s ^= m_s << 13;
            m_s ^= m_s >> 17;
            m_s ^= m_s << 5;
            return (static_cast<double>(m_s % 65536U) - 32768.0) / 32768.0;
        }

      private:
        std::uint32_t m_s;
    };

    std::uint64_t fnv1a64(const void* p, std::size_t n) noexcept {
        const auto*   b = static_cast<const unsigned char*>(p);
        std::uint64_t h = 1469598103934665603ULL;
        for (std::size_t i = 0; i < n; ++i) {
            h ^= b[i];
            h *= 1099511628211ULL;
        }
        return h;
    }

    // The corpus: -20 dB xorshift32 noise under a 220 Hz tone (so pvoc has a
    // peak to lock) and an off-bin 1234.5 Hz tone (so it has one to translate
    // fractionally), built in double and cast per profile.
    std::vector<double> corpus() {
        std::vector<double> x(k_samples);
        xorshift32          rng(0x2545F491U);
        for (std::size_t i = 0; i < k_samples; ++i) {
            // Association order is part of the corpus: 2*pi*f*i/sr, left to right.
            const double i_d = static_cast<double>(i);
            x[i]             = 0.1 * rng.next() + 0.5 * std::sin(2.0 * std::numbers::pi * 220.0 * i_d / k_sample_rate)
                   + 0.3 * std::sin(2.0 * std::numbers::pi * 1234.5 * i_d / k_sample_rate);
        }
        return x;
    }

    void print(const char* primitive, const char* profile, const char* config, std::uint64_t hash) {
        std::printf("%-8s %-6s %-22s %016llx\n", primitive, profile, config, static_cast<unsigned long long>(hash));
    }

    template <typename Sample>
    void fingerprint_pvoc(const std::vector<double>& x, const char* profile, double ratio, bool formant,
                          const char* config) {
        tap::dsp::basic_pvoc<Sample> shifter(k_pvoc_fft);
        shifter.set_formant(formant);
        std::vector<Sample> out(x.size());
        for (std::size_t i = 0; i < x.size(); ++i) {
            out[i] = shifter.process(static_cast<Sample>(x[i]), static_cast<Sample>(ratio));
        }
        print("pvoc", profile, config, fnv1a64(out.data(), out.size() * sizeof(Sample)));
    }

    template <typename Sample>
    void fingerprint_log_mel(const std::vector<double>& x, const char* profile, bool pcen, const char* config) {
        tap::dsp::log_mel_geometry g;
        g.pcen.enabled = pcen;
        tap::dsp::basic_log_mel<Sample> front_end(g);
        std::vector<Sample>             in(x.size());
        for (std::size_t i = 0; i < x.size(); ++i) {
            in[i] = static_cast<Sample>(x[i]);
        }
        const std::size_t   frames = front_end.frames_for(in.size());
        std::vector<Sample> out(frames * g.bands);
        const std::size_t   written = front_end.process(in.data(), in.size(), out.data(), frames);
        print("log_mel", profile, config, fnv1a64(out.data(), written * g.bands * sizeof(Sample)));
    }

    template <typename Sample>
    void fingerprint_profile(const std::vector<double>& x, const char* profile) {
        fingerprint_pvoc<Sample>(x, profile, 1.0, false, "ratio=1.0 formant=off");
        fingerprint_pvoc<Sample>(x, profile, 1.5, false, "ratio=1.5 formant=off");
        fingerprint_pvoc<Sample>(x, profile, 1.0, true, "ratio=1.0 formant=on");
        fingerprint_pvoc<Sample>(x, profile, 1.5, true, "ratio=1.5 formant=on");
        fingerprint_log_mel<Sample>(x, profile, false, "log");
        fingerprint_log_mel<Sample>(x, profile, true, "pcen");
    }

} // namespace

int main() {
    const std::vector<double> x = corpus();
    fingerprint_profile<float>(x, "float");
    fingerprint_profile<double>(x, "double");
    return 0;
}
