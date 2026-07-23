// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// Contract battery for the FIR dot kernels. The load-bearing promise —
// extracted from SampleRateTap's multichannel suite, where it gates the
// frame-major fast path — is bit-exactness: the channel-parallel kernel must
// produce the identical bits as the planar dot for every sample type, because
// consumers switch between the layouts by channel count and target.

#include <cstdint>
#include <vector>

#include <gtest/gtest.h>

#include "tap/dsp/fir_kernels.h"

namespace {

    using tap::dsp::dot_row;
    using tap::dsp::dot_rows_frame_major;
    using tap::dsp::sample_traits;

    // Deterministic pseudo-random generator (xorshift), mapped per sample type
    // to well-inside-full-scale values so no test depends on saturation.
    inline std::uint32_t next(std::uint32_t& s) {
        s ^= s << 13;
        s ^= s >> 17;
        s ^= s << 5;
        return s;
    }

    template <typename S>
    S sample_from(std::uint32_t r);
    template <>
    float sample_from<float>(std::uint32_t r) {
        return (static_cast<float>(r % 65536) - 32768.0f) / 65536.0f;
    }
    template <>
    std::int16_t sample_from<std::int16_t>(std::uint32_t r) {
        return static_cast<std::int16_t>(static_cast<int>(r % 32768) - 16384);
    }
    template <>
    std::int32_t sample_from<std::int32_t>(std::uint32_t r) {
        return static_cast<std::int32_t>(r) / 4;
    }

    template <typename S>
    typename sample_traits<S>::coeff coeff_from(std::uint32_t r) {
        // Coefficients in roughly [-0.5, 0.5) of the format's unity.
        return sample_traits<S>::make_coeff((static_cast<double>(r % 4096) - 2048.0) / 4096.0);
    }

    template <typename S>
    class fir_kernels_test : public ::testing::Test {};
    using sample_types = ::testing::Types<float, std::int16_t, std::int32_t>;
    TYPED_TEST_SUITE(fir_kernels_test, sample_types, );

    // dot_row must equal the reference accumulation: mac per tap in order,
    // one finalize. (On SMLALD targets this also pins the dual-MAC pairing;
    // pairing is bit-free because each 16x16 product is exact in int32.)
    TYPED_TEST(fir_kernels_test, DotRowMatchesReferenceAccumulation) {
        using sample                           = TypeParam;
        using tr                               = sample_traits<sample>;
        constexpr std::size_t           k_taps = 48;
        std::uint32_t                   seed   = 0x12345678u;
        std::vector<sample>             hist(k_taps);
        std::vector<typename tr::coeff> row(k_taps);
        for (std::size_t t = 0; t < k_taps; ++t) {
            hist[t] = sample_from<sample>(next(seed));
            row[t]  = coeff_from<sample>(next(seed));
        }
        typename tr::accum acc{};
        for (std::size_t t = 0; t < k_taps; ++t) {
            acc = tr::mac(acc, hist[t], row[t]);
        }
        EXPECT_EQ(dot_row<sample>(row.data(), hist.data(), k_taps), tr::finalize(acc));
    }

    // The frame-major channel-parallel kernel is bit-exact against the planar
    // dot for every channel count that exercises the 8/4/2/1 tiling, float
    // included (lanes are channels, not taps: each channel's double
    // accumulation order is unchanged).
    TYPED_TEST(fir_kernels_test, ChannelParallelMatchesPlanarBitExact) {
        using sample                 = TypeParam;
        using tr                     = sample_traits<sample>;
        constexpr std::size_t k_taps = 48;
        for (std::size_t channels : {1u, 2u, 3u, 4u, 7u, 8u, 11u, 12u, 16u}) {
            std::uint32_t                   seed = 0x9e3779b9u + static_cast<std::uint32_t>(channels);
            std::vector<sample>             x(k_taps * channels); // frame-major
            std::vector<typename tr::coeff> row(k_taps);
            for (auto& v : x) {
                v = sample_from<sample>(next(seed));
            }
            for (auto& c : row) {
                c = coeff_from<sample>(next(seed));
            }
            std::vector<sample> out(channels);
            dot_rows_frame_major<sample>(row.data(), x.data(), k_taps, channels, out.data());
            for (std::size_t c = 0; c < channels; ++c) {
                std::vector<sample> planar(k_taps);
                for (std::size_t t = 0; t < k_taps; ++t) {
                    planar[t] = x[t * channels + c];
                }
                EXPECT_EQ(out[c], dot_row<sample>(row.data(), planar.data(), k_taps))
                    << "channels=" << channels << " c=" << c;
            }
        }
    }

    // Odd tap counts exercise dot_row's scalar tail on SMLALD targets; on
    // hosts this just pins the same contract at a second geometry.
    TYPED_TEST(fir_kernels_test, OddTapCountMatchesReference) {
        using sample                           = TypeParam;
        using tr                               = sample_traits<sample>;
        constexpr std::size_t           k_taps = 33;
        std::uint32_t                   seed   = 0xdeadbeefu;
        std::vector<sample>             hist(k_taps);
        std::vector<typename tr::coeff> row(k_taps);
        for (std::size_t t = 0; t < k_taps; ++t) {
            hist[t] = sample_from<sample>(next(seed));
            row[t]  = coeff_from<sample>(next(seed));
        }
        typename tr::accum acc{};
        for (std::size_t t = 0; t < k_taps; ++t) {
            acc = tr::mac(acc, hist[t], row[t]);
        }
        EXPECT_EQ(dot_row<sample>(row.data(), hist.data(), k_taps), tr::finalize(acc));
    }

} // namespace
