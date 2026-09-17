/// @file bench_fft.cpp
/// @brief Host wall-clock microbenchmark for the real FFT: min-of-N ns per transform.
// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// Informational only (bench/README.md): wall clock on a shared runner is
// noise, so nothing gates on these numbers. This is the local tool for the
// desktop and Apple vDSP claims — record the machine, compiler, load average
// and date with any number that goes into docs/fft-design.md. The gate is the
// instruction-count ratchet in icount/.
//
// Same scenarios and engine selector as the ratchet (bench_common.h):
// rfft_f32_512, rfft_f32_2048, rfft_f64_512. Each rep runs a batch of
// out-of-place forward() calls and, separately, a batch of inverse() calls
// (the class's scaled inverse), and the minimum over reps is reported per
// transform. From Stage 2a the port and the C build side by side here via
// TAP_DSP_BENCH_ENGINE; after Stage 2c the recorded C numbers are the
// comparison.
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <vector>

#include "bench_common.h"

namespace {

    constexpr std::size_t k_reps          = 25;
    constexpr std::size_t k_batch_samples = std::size_t{1} << 20; // per batch, per direction
    constexpr std::size_t k_corpus_blocks = 4;

    struct result {
        double forward_ns;
        double inverse_ns;
    };

    using clock = std::chrono::steady_clock;

    double ns_since(clock::time_point t0) noexcept {
        return std::chrono::duration<double, std::nano>(clock::now() - t0).count();
    }

    template <typename Sample>
    result measure(std::size_t n, Sample& sink) {
        tap::dsp::bench::fft_under_test<Sample> fft(n);
        const std::size_t                       batch = k_batch_samples / n;

        std::vector<Sample>         corpus(k_corpus_blocks * n);
        std::vector<Sample>         spectrum(n);
        std::vector<Sample>         out(n);
        tap::dsp::bench::xorshift32 rng(0x9E3779B9u);
        rng.fill(corpus.data(), corpus.size());

        Sample acc = Sample(0);
        result best{1e300, 1e300};
        for (std::size_t rep = 0; rep < k_reps; ++rep) {
            const auto t0 = clock::now();
            for (std::size_t i = 0; i < batch; ++i) {
                fft.forward(corpus.data() + (i % k_corpus_blocks) * n, spectrum.data());
                acc += spectrum[i % n];
            }
            const double fwd = ns_since(t0) / static_cast<double>(batch);

            const auto t1 = clock::now();
            for (std::size_t i = 0; i < batch; ++i) {
                fft.inverse(spectrum.data(), out.data());
                acc += out[i % n];
            }
            const double inv = ns_since(t1) / static_cast<double>(batch);

            best.forward_ns = fwd < best.forward_ns ? fwd : best.forward_ns;
            best.inverse_ns = inv < best.inverse_ns ? inv : best.inverse_ns;
        }
        sink += acc;
        return best;
    }

    void report(const char* name, const result& r) {
        std::printf("%-16s %12.1f %12.1f %14.1f\n", name, r.forward_ns, r.inverse_ns, r.forward_ns + r.inverse_ns);
    }

} // namespace

int main() {
    std::printf("tap_dsp_bench_fft  engine=%s  min of %zu reps, %zu samples per batch per direction\n",
                tap::dsp::bench::k_engine_name.data(), k_reps, k_batch_samples);
    std::printf("%-16s %12s %12s %14s\n", "scenario", "forward ns", "inverse ns", "round-trip ns");

    float  sink32 = 0.0f;
    double sink64 = 0.0;
    report("rfft_f32_512", measure<float>(512, sink32));
    report("rfft_f32_2048", measure<float>(2048, sink32));
    report("rfft_f64_512", measure<double>(512, sink64));

    // Keep the transforms observable; wall clock is not gated on this value.
    const bool ok = sink32 == sink32 && sink64 == sink64;
    std::printf("checksum f32=%.9g f64=%.17g ok=%d\n", static_cast<double>(sink32), sink64, ok ? 1 : 0);
    return ok ? 0 : 1;
}
