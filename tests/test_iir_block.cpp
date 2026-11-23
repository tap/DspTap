// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// Pins iir_kernels.h's biquad_cascade_block against the per-sample
// biquad_cascade, bit for bit: the output and the state it leaves, for
// K = 1 ... 8 sections, channel counts that leave every remainder of the
// group (1, 2, 3, 6, 12), any split of the frames into calls (empty, one
// frame, odd, uneven), in place, and the two kernels alternating on one
// state; a section with b2 = a2 = 0 (a first-order row) among them; a NaN
// poisoning the same words; the group the kernel picks per format and
// section count. Built without fp contraction (tap_dsp_chain_stages,
// tests/CMakeLists.txt): that is the condition the kernel states for the
// identity, and with it on a compiler may fuse the two kernels differently.
// Runs on every leg, the QEMU ones too.
#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <span>
#include <vector>

#include <gtest/gtest.h>

#include "support/iir_support.h"
#include "tap/dsp/butterworth.h"
#include "tap/dsp/iir_kernels.h"

namespace {

    using tap::dsp::biquad_cascade;
    using tap::dsp::biquad_cascade_block;
    using tap::dsp::biquad_cascade_table;
    using tap::dsp::iir_section;
    using tap::dsp::test::flush_to_zero_fixture;

    template <typename S>
    class biquad_block_test : public flush_to_zero_fixture {};
    using floating_types = ::testing::Types<float, double>;
    TYPED_TEST_SUITE(biquad_block_test, floating_types, );

    /// K stable sections: a Butterworth's where 2 K is an order it designs,
    /// with the last row made first order (b2 = a2 = 0) when asked.
    template <std::size_t K>
    std::array<iir_section, K> sections(bool first_order_last) {
        std::array<iir_section, K> sec = tap::dsp::design_butterworth_lowpass<2 * K>(0.11);
        if (first_order_last) {
            sec[K - 1] = {0.3, 0.3, 0.0, -0.4, 0.0};
        }
        return sec;
    }

    /// Bits, not values: -0 against +0 and one NaN against another must show.
    template <typename S>
    bool same_bits(const std::vector<S>& a, const std::vector<S>& b) {
        return a.size() == b.size() && (a.empty() || std::memcmp(a.data(), b.data(), a.size() * sizeof(S)) == 0);
    }

    template <typename S>
    std::vector<S> noise(std::size_t n, std::uint32_t seed) {
        std::vector<S> x(n);
        for (S& v : x) {
            seed ^= seed << 13;
            seed ^= seed >> 17;
            seed ^= seed << 5;
            v = static_cast<S>(static_cast<double>(seed) / 4294967296.0 - 0.5);
        }
        return x;
    }

    /// The reference: biquad_cascade per channel per frame, as the stage ran it through 0.6.0.
    template <typename S, std::size_t K, std::size_t Ch>
    void per_sample(const biquad_cascade_table<S>& t, std::vector<S>& state, const S* in, std::size_t n, S* out) {
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t c = 0; c < Ch; ++c) {
                out[i * Ch + c] = biquad_cascade<S>(t, state.data() + c * 2 * K, in[i * Ch + c]);
            }
        }
    }

    /// Frames per call, cycled: empty, single, odd and uneven calls.
    constexpr std::array<std::size_t, 9> k_calls{1, 0, 7, 1, 64, 2, 33, 128, 5};

    template <typename S, std::size_t K, std::size_t Ch>
    void expect_block_is_per_sample(bool first_order_last) {
        SCOPED_TRACE(::testing::Message()
                     << "K = " << K << ", Ch = " << Ch << ", first-order last row " << first_order_last);
        const auto                    sec = sections<K>(first_order_last);
        const biquad_cascade_table<S> t{std::span<const iir_section>(sec)};
        constexpr std::size_t         k_frames = 600;
        const std::vector<S>          x        = noise<S>(k_frames * Ch, 0x2545F491u + static_cast<std::uint32_t>(K));

        std::vector<S> ref_state(Ch * 2 * K, S(0));
        std::vector<S> ref(x.size());
        per_sample<S, K, Ch>(t, ref_state, x.data(), k_frames, ref.data());

        // Uneven calls, out of place.
        std::vector<S> state(Ch * 2 * K, S(0));
        std::vector<S> y(x.size(), S(-1));
        std::size_t    done = 0;
        for (std::size_t call = 0; done < k_frames; ++call) {
            const std::size_t n = std::min(k_calls[call % k_calls.size()], k_frames - done);
            biquad_cascade_block<S, K, Ch>(t, state.data(), x.data() + done * Ch, n, y.data() + done * Ch);
            done += n;
        }
        EXPECT_TRUE(same_bits(y, ref)) << "uneven calls";
        EXPECT_TRUE(same_bits(state, ref_state)) << "the state after uneven calls";

        // One call, in place.
        std::vector<S> in_place = x;
        std::fill(state.begin(), state.end(), S(0));
        biquad_cascade_block<S, K, Ch>(t, state.data(), in_place.data(), k_frames, in_place.data());
        EXPECT_TRUE(same_bits(in_place, ref)) << "in place";
        EXPECT_TRUE(same_bits(state, ref_state)) << "the state after one call in place";

        // The two kernels alternating on one state.
        std::fill(state.begin(), state.end(), S(0));
        std::fill(y.begin(), y.end(), S(-1));
        done = 0;
        for (std::size_t call = 0; done < k_frames; ++call) {
            const std::size_t n = std::min(k_calls[(call + 2) % k_calls.size()], k_frames - done);
            if (call % 2 == 0) {
                biquad_cascade_block<S, K, Ch>(t, state.data(), x.data() + done * Ch, n, y.data() + done * Ch);
            }
            else {
                per_sample<S, K, Ch>(t, state, x.data() + done * Ch, n, y.data() + done * Ch);
            }
            done += n;
        }
        EXPECT_TRUE(same_bits(y, ref)) << "the kernels alternating";
        EXPECT_TRUE(same_bits(state, ref_state)) << "the state after the kernels alternating";
    }

    template <typename S, std::size_t K>
    void expect_every_channel_count(bool first_order_last) {
        expect_block_is_per_sample<S, K, 1>(first_order_last);
        expect_block_is_per_sample<S, K, 2>(first_order_last);
        expect_block_is_per_sample<S, K, 3>(first_order_last);
        expect_block_is_per_sample<S, K, 6>(first_order_last);
        expect_block_is_per_sample<S, K, 12>(first_order_last);
    }

    TYPED_TEST(biquad_block_test, BlockIsThePerSampleKernelBitForBit) {
        expect_every_channel_count<TypeParam, 1>(false);
        expect_every_channel_count<TypeParam, 2>(false);
        expect_every_channel_count<TypeParam, 3>(false);
        expect_every_channel_count<TypeParam, 4>(false);
        expect_every_channel_count<TypeParam, 5>(false);
        expect_every_channel_count<TypeParam, 6>(false);
        expect_every_channel_count<TypeParam, 7>(false);
        expect_every_channel_count<TypeParam, 8>(false);
    }

    TYPED_TEST(biquad_block_test, AFirstOrderRowIsThePerSampleKernelBitForBit) {
        expect_every_channel_count<TypeParam, 1>(true);
        expect_every_channel_count<TypeParam, 2>(true);
        expect_every_channel_count<TypeParam, 4>(true);
        expect_every_channel_count<TypeParam, 8>(true);
    }

    // The honest limit both kernels share: a non-finite input poisons the
    // state. The block kernel leaves the same words as the per-sample one.
    template <typename S>
    void expect_a_nan_poisons_the_same_words() {
        constexpr std::size_t         k   = 2;
        constexpr std::size_t         ch  = 3;
        const auto                    sec = sections<k>(false);
        const biquad_cascade_table<S> t{std::span<const iir_section>(sec)};
        std::vector<S>                x = noise<S>(40 * ch, 7u);
        x[17 * ch + 1]                  = std::numeric_limits<S>::quiet_NaN(); // channel 1 alone
        std::vector<S> ref_state(ch * 2 * k, S(0)), state(ch * 2 * k, S(0));
        std::vector<S> ref(x.size()), y(x.size());
        per_sample<S, k, ch>(t, ref_state, x.data(), 40, ref.data());
        biquad_cascade_block<S, k, ch>(t, state.data(), x.data(), 40, y.data());
        EXPECT_TRUE(same_bits(y, ref));
        EXPECT_TRUE(same_bits(state, ref_state));
        EXPECT_NE(y[39 * ch + 1], y[39 * ch + 1]) << "channel 1 stays poisoned";
        EXPECT_EQ(y[39 * ch + 0], y[39 * ch + 0]) << "channel 0 does not see it";
    }

    TYPED_TEST(biquad_block_test, ANonFiniteInputPoisonsTheSameWords) {
        expect_a_nan_poisons_the_same_words<TypeParam>();
    }

    // The group is a measured choice (bench/README.md, "The block kernel"):
    // a change here is a change to every count the ratchet holds for the stage.
    TEST(biquad_block, TheGroupIsTwoForFloatAtOneOrTwoSectionsAndOneOtherwise) {
        using tap::dsp::detail::k_biquad_block_group;
        EXPECT_EQ((k_biquad_block_group<float, 1>), 2u);
        EXPECT_EQ((k_biquad_block_group<float, 2>), 2u);
        EXPECT_EQ((k_biquad_block_group<float, 3>), 1u);
        EXPECT_EQ((k_biquad_block_group<float, 8>), 1u);
        EXPECT_EQ((k_biquad_block_group<double, 1>), 1u);
        EXPECT_EQ((k_biquad_block_group<double, 2>), 1u);
        EXPECT_EQ((k_biquad_block_group<double, 8>), 1u);
    }

} // namespace
