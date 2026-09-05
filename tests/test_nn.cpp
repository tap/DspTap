// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// Contract battery for tap::dsp::nn — the dense and GRU inference kernels.
// Every promise in nn.h's header comment is pinned here as a number: the
// row-major [out x in] layout, the activation forms, the r/z/n gate order
// and the PyTorch GRU formula (against an independent long-double
// restatement that sums in the opposite order, so agreement proves the
// formula rather than the code), reset semantics, ownership of moved-in
// weights, noexcept processing, and float-tracks-double as a measured
// number at MuTap's suppressor geometry.
#include <cmath>
#include <cstdint>
#include <type_traits>
#include <vector>

#include <gtest/gtest.h>

#include "tap/dsp/nn.h"

namespace {

    using tap::dsp::nn::activation;
    using tap::dsp::nn::basic_dense;
    using tap::dsp::nn::basic_gru;

    class xorshift32 {
      public:
        explicit xorshift32(std::uint32_t seed)
            : m_s(seed) {}
        /// Uniform in [-1, 1).
        float next() {
            m_s ^= m_s << 13;
            m_s ^= m_s >> 17;
            m_s ^= m_s << 5;
            return static_cast<float>(m_s) / 2147483648.0f - 1.0f;
        }

      private:
        std::uint32_t m_s;
    };

    std::vector<float> random_weights(xorshift32& rng, std::size_t n, float scale) {
        std::vector<float> w(n);
        for (auto& v : w) {
            v = rng.next() * scale;
        }
        return w;
    }

    // Independent restatement of the PyTorch GRU cell in long double, summing
    // each dot product from the LAST input to the first (the kernel sums
    // bias-first, ascending), so agreement with the kernel tests the formula
    // and the gate order, not shared code.
    struct gru_reference {
        std::vector<float>       w_ih, w_hh, b_ih, b_hh;
        std::size_t              hidden, in;
        std::vector<long double> h;

        static long double sigmoid(long double v) { return 1.0L / (1.0L + std::exp(-v)); }

        void step(const std::vector<long double>& x) {
            std::vector<long double> gi(3 * hidden), gh(3 * hidden);
            for (std::size_t i = 0; i < 3 * hidden; ++i) {
                long double a = 0.0L;
                for (std::size_t j = in; j-- > 0;) {
                    a += static_cast<long double>(w_ih[i * in + j]) * x[j];
                }
                gi[i] = a + static_cast<long double>(b_ih[i]);
                a     = 0.0L;
                for (std::size_t j = hidden; j-- > 0;) {
                    a += static_cast<long double>(w_hh[i * hidden + j]) * h[j];
                }
                gh[i] = a + static_cast<long double>(b_hh[i]);
            }
            for (std::size_t i = 0; i < hidden; ++i) {
                const long double r = sigmoid(gi[i] + gh[i]);
                const long double z = sigmoid(gi[hidden + i] + gh[hidden + i]);
                const long double n = std::tanh(gi[2 * hidden + i] + r * gh[2 * hidden + i]);
                h[i]                = (1.0L - z) * n + z * h[i];
            }
        }
    };

    template <typename Sample>
    class nn_test : public ::testing::Test {};
    using sample_types = ::testing::Types<float, double>;
    TYPED_TEST_SUITE(nn_test, sample_types);

    // ------------------------------------------------------------- dense

    TYPED_TEST(nn_test, DenseLayoutIsRowMajorOutByIn) {
        // W = [1 2 3; 10 20 30], b = [0.5, -0.5], x = [1, -1, 2]:
        // y0 = 0.5 + 1 - 2 + 6 = 5.5, y1 = -0.5 + 10 - 20 + 60 = 49.5.
        basic_dense<TypeParam> d({1, 2, 3, 10, 20, 30}, {0.5f, -0.5f}, 2, 3);
        const TypeParam        x[3] = {TypeParam(1), TypeParam(-1), TypeParam(2)};
        TypeParam              y[2] = {TypeParam(0), TypeParam(0)};
        d.apply(x, y);
        EXPECT_EQ(y[0], TypeParam(5.5));
        EXPECT_EQ(y[1], TypeParam(49.5));
        EXPECT_EQ(d.out(), 2u);
        EXPECT_EQ(d.in(), 3u);
        EXPECT_EQ(d.act(), activation::linear);
    }

    TYPED_TEST(nn_test, DenseActivationsAreTanhAndLogistic) {
        const TypeParam        x[1] = {TypeParam(0.75)};
        TypeParam              y[1];
        basic_dense<TypeParam> lin({2.0f}, {-0.5f}, 1, 1, activation::linear);
        basic_dense<TypeParam> th({2.0f}, {-0.5f}, 1, 1, activation::tanh);
        basic_dense<TypeParam> sg({2.0f}, {-0.5f}, 1, 1, activation::sigmoid);
        lin.apply(x, y);
        EXPECT_EQ(y[0], TypeParam(1)); // -0.5 + 2 * 0.75
        th.apply(x, y);
        EXPECT_EQ(y[0], std::tanh(TypeParam(1)));
        sg.apply(x, y);
        EXPECT_EQ(y[0], TypeParam(1) / (TypeParam(1) + std::exp(TypeParam(-1))));
    }

    TYPED_TEST(nn_test, SigmoidIsTheLogisticFunction) {
        using tap::dsp::nn::sigmoid;
        EXPECT_EQ(sigmoid(TypeParam(0)), TypeParam(0.5));
        for (const double v : {-6.0, -1.5, 0.25, 3.0}) {
            const TypeParam s = sigmoid(static_cast<TypeParam>(v));
            EXPECT_EQ(s, TypeParam(1) / (TypeParam(1) + std::exp(static_cast<TypeParam>(-v))));
            constexpr double k_tol = std::is_same_v<TypeParam, float> ? 1e-6 : 1e-15;
            EXPECT_NEAR(static_cast<double>(s + sigmoid(static_cast<TypeParam>(-v))), 1.0, k_tol);
        }
    }

    // --------------------------------------------------------------- GRU

    TYPED_TEST(nn_test, GruGatesAreOrderedResetUpdateNew) {
        // hidden 1, in 1, all weights zero: the biases alone drive the gates,
        // so each block of the 3*hidden axis can be identified by its effect.
        const std::vector<float> zero_w = {0.0f, 0.0f, 0.0f};
        const TypeParam          x[1]   = {TypeParam(0)};

        // z block (index 1) saturated high: h' = h, the state never moves.
        basic_gru<TypeParam> hold(zero_w, zero_w, {0.0f, 40.0f, 1.0f}, {0.0f, 0.0f, 0.0f}, 1, 1);
        hold.step(x);
        EXPECT_NEAR(static_cast<double>(hold.state()[0]), 0.0, 1e-15);

        // z block saturated low: h' = n = tanh(b_in + r * b_hn) with
        // r = sigmoid(b_ir + b_hr). Choose b_ir = b_hr = 0 (r = 0.5),
        // b_in = 0.2, b_hn = 0.6: n = tanh(0.5).
        basic_gru<TypeParam> take(zero_w, zero_w, {0.0f, -40.0f, 0.2f}, {0.0f, 0.0f, 0.6f}, 1, 1);
        take.step(x);
        EXPECT_NEAR(static_cast<double>(take.state()[0]), std::tanh(0.5), 1e-6);

        // r block saturated low (r = 0): the hidden new-gate term drops out,
        // n = tanh(b_in) = tanh(0.2).
        basic_gru<TypeParam> no_r(zero_w, zero_w, {-40.0f, -40.0f, 0.2f}, {0.0f, 0.0f, 0.6f}, 1, 1);
        no_r.step(x);
        EXPECT_NEAR(static_cast<double>(no_r.state()[0]), std::tanh(0.2), 1e-6);
    }

    TYPED_TEST(nn_test, GruMatchesAnIndependentReference) {
        constexpr std::size_t k_hidden = 8;
        constexpr std::size_t k_in     = 5;
        constexpr std::size_t k_steps  = 50;
        xorshift32            rng(0x2545F491u);
        gru_reference         ref{random_weights(rng, 3 * k_hidden * k_in, 0.5f),
                          random_weights(rng, 3 * k_hidden * k_hidden, 0.5f),
                          random_weights(rng, 3 * k_hidden, 0.5f),
                          random_weights(rng, 3 * k_hidden, 0.5f),
                          k_hidden,
                          k_in,
                          std::vector<long double>(k_hidden, 0.0L)};
        basic_gru<TypeParam>  g(ref.w_ih, ref.w_hh, ref.b_ih, ref.b_hh, k_hidden, k_in);

        double worst = 0.0;
        for (std::size_t t = 0; t < k_steps; ++t) {
            std::vector<TypeParam>   x(k_in);
            std::vector<long double> xl(k_in);
            for (std::size_t j = 0; j < k_in; ++j) {
                x[j]  = static_cast<TypeParam>(rng.next());
                xl[j] = static_cast<long double>(x[j]);
            }
            g.step(x.data());
            ref.step(xl);
            for (std::size_t i = 0; i < k_hidden; ++i) {
                worst = std::max(worst, std::abs(static_cast<double>(g.state()[i]) - static_cast<double>(ref.h[i])));
            }
        }
        // Measured 2026-09: double 2.2e-16, float 1.0e-7. Pinned at ~10x.
        constexpr double k_tol = std::is_same_v<TypeParam, float> ? 1e-6 : 1e-14;
        EXPECT_LT(worst, k_tol);
    }

    TYPED_TEST(nn_test, GruResetClearsTheState) {
        xorshift32           rng(7u);
        basic_gru<TypeParam> g(random_weights(rng, 3 * 4 * 3, 0.5f), random_weights(rng, 3 * 4 * 4, 0.5f),
                               random_weights(rng, 12, 0.5f), random_weights(rng, 12, 0.5f), 4, 3);
        const TypeParam      x[3] = {TypeParam(0.3), TypeParam(-0.2), TypeParam(0.9)};
        g.step(x);
        bool moved = false;
        for (std::size_t i = 0; i < 4; ++i) {
            moved = moved || g.state()[i] != TypeParam(0);
        }
        EXPECT_TRUE(moved);
        std::vector<TypeParam> first(g.state(), g.state() + 4);
        g.reset();
        for (std::size_t i = 0; i < 4; ++i) {
            EXPECT_EQ(g.state()[i], TypeParam(0));
        }
        g.step(x); // deterministic: the same first step after reset
        for (std::size_t i = 0; i < 4; ++i) {
            EXPECT_EQ(g.state()[i], first[i]);
        }
    }

    TYPED_TEST(nn_test, LayersOwnTheirWeightsAndCopyDeeply) {
        std::vector<float>     w = {1, 2, 3, 4};
        std::vector<float>     b = {0, 0};
        basic_dense<TypeParam> d(std::move(w), std::move(b), 2, 2);
        EXPECT_TRUE(w.empty()); // moved in, not copied  // NOLINT(bugprone-use-after-move)
        basic_dense<TypeParam> copy = d;
        const TypeParam        x[2] = {TypeParam(1), TypeParam(1)};
        TypeParam              y1[2], y2[2];
        d.apply(x, y1);
        copy.apply(x, y2);
        EXPECT_EQ(y1[0], TypeParam(3));
        EXPECT_EQ(y1[1], TypeParam(7));
        EXPECT_EQ(y2[0], y1[0]);
        EXPECT_EQ(y2[1], y1[1]);
        EXPECT_NE(copy.weights().data(), d.weights().data());
    }

    TYPED_TEST(nn_test, ProcessingPathIsNoexcept) {
        static_assert(noexcept(std::declval<const basic_dense<TypeParam>&>().apply(nullptr, nullptr)));
        static_assert(noexcept(std::declval<basic_gru<TypeParam>&>().step(nullptr)));
        static_assert(noexcept(std::declval<basic_gru<TypeParam>&>().reset()));
        SUCCEED();
    }

    // ---------------------------------------------- cross-precision pin

    TEST(NnCrossPrecision, FloatTracksDoubleAtTheSuppressorGeometry) {
        // MuTap's residual suppressor: features -> dense(64, tanh) -> GRU(96)
        // -> dense(bands, sigmoid), 0.3-scaled random weights, 200 steps.
        constexpr std::size_t k_feat  = 56;
        constexpr std::size_t k_dense = 64;
        constexpr std::size_t k_gru   = 96;
        constexpr std::size_t k_bands = 26;
        constexpr std::size_t k_steps = 200;
        xorshift32            rng(0x9E3779B9u);
        const auto            w_in  = random_weights(rng, k_dense * k_feat, 0.3f);
        const auto            b_in  = random_weights(rng, k_dense, 0.3f);
        const auto            w_ih  = random_weights(rng, 3 * k_gru * k_dense, 0.3f);
        const auto            w_hh  = random_weights(rng, 3 * k_gru * k_gru, 0.3f);
        const auto            b_ih  = random_weights(rng, 3 * k_gru, 0.3f);
        const auto            b_hh  = random_weights(rng, 3 * k_gru, 0.3f);
        const auto            w_out = random_weights(rng, k_bands * k_gru, 0.3f);
        const auto            b_out = random_weights(rng, k_bands, 0.3f);

        basic_dense<double> din(w_in, b_in, k_dense, k_feat, activation::tanh);
        basic_dense<float>  din32(w_in, b_in, k_dense, k_feat, activation::tanh);
        basic_gru<double>   g(w_ih, w_hh, b_ih, b_hh, k_gru, k_dense);
        basic_gru<float>    g32(w_ih, w_hh, b_ih, b_hh, k_gru, k_dense);
        basic_dense<double> dout(w_out, b_out, k_bands, k_gru, activation::sigmoid);
        basic_dense<float>  dout32(w_out, b_out, k_bands, k_gru, activation::sigmoid);

        std::vector<double> f(k_feat), h(k_dense), gains(k_bands);
        std::vector<float>  f32(k_feat), h32(k_dense), gains32(k_bands);
        double              worst_gain = 0.0, worst_state = 0.0;
        for (std::size_t t = 0; t < k_steps; ++t) {
            for (std::size_t j = 0; j < k_feat; ++j) {
                f32[j] = rng.next();
                f[j]   = static_cast<double>(f32[j]);
            }
            din.apply(f.data(), h.data());
            din32.apply(f32.data(), h32.data());
            g.step(h.data());
            g32.step(h32.data());
            dout.apply(g.state(), gains.data());
            dout32.apply(g32.state(), gains32.data());
            for (std::size_t i = 0; i < k_gru; ++i) {
                worst_state = std::max(worst_state, std::abs(static_cast<double>(g32.state()[i]) - g.state()[i]));
            }
            for (std::size_t b = 0; b < k_bands; ++b) {
                worst_gain = std::max(worst_gain, std::abs(static_cast<double>(gains32[b]) - gains[b]));
            }
        }
        // Measured 2026-09: state 6.3e-7, gains 2.5e-7. Pinned at ~4-5x.
        EXPECT_LT(worst_state, 3e-6);
        EXPECT_LT(worst_gain, 1e-6);
    }

} // namespace
