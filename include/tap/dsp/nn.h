/// @file nn.h
/// @brief Dense and GRU inference kernels with a fixed numeric contract.
// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// The two layers a small recurrent gain network is made of — a fully
// connected layer and a gated recurrent unit cell — as allocation-free,
// noexcept kernels over a Sample type. Lifted from MuTap's learned residual
// suppressor (its wake-word plan, milestone M3), whose float-tracks-double
// and Python-parity pins are the oracle this header must keep satisfying;
// the keyword spotter is a different head on the same arithmetic. The
// formulas are the textbook ones (the GRU is Cho et al. 2014 in PyTorch's
// nn.GRU convention); nothing here is reverse-engineered from a product.
//
// The contract, as numbers (k_contract_version changes when any of it does):
//   - Weights are stored as float32 whatever the Sample type — the storage
//     precision of a trained model's file — and each weight is converted to
//     Sample at the point of use. The double profile therefore computes
//     with float32-valued weights in double arithmetic, which is exactly
//     what a float64 numpy reference over the same file computes.
//   - Layouts are row-major [out x in]: W[i * in + j] multiplies input j
//     into output i. Bias vectors are [out].
//   - Dense: y[i] = act(b[i] + sum_{j=0}^{in-1} W[i*in + j] * x[j]). The
//     accumulator starts at the bias and takes the terms in ascending j, in
//     Sample; then the activation: linear (identity), tanh (std::tanh), or
//     sigmoid, sigmoid(v) = 1 / (1 + exp(-v)) in Sample.
//   - GRU cell (PyTorch nn.GRU, gates ordered r, z, n along the 3*hidden
//     axis of every weight and bias):
//       g_ih = b_ih + W_ih x        (3*hidden, W_ih is [3*hidden x in])
//       g_hh = b_hh + W_hh h        (3*hidden, W_hh is [3*hidden x hidden])
//       r = sigmoid(g_ih[i]            + g_hh[i])
//       z = sigmoid(g_ih[hidden + i]   + g_hh[hidden + i])
//       n = tanh   (g_ih[2*hidden + i] + r * g_hh[2*hidden + i])
//       h'[i] = (1 - z) * n + z * h[i]
//     evaluated in that order, every sum in Sample with the same
//     bias-first ascending-j accumulation as the dense layer. The state is
//     zero after construction and after reset().
//   - Cost per step: 3*hidden*(in + hidden) multiply-adds for the GRU,
//     out*in for a dense layer, plus one transcendental per gate element.
//
// Profiles: double is the golden model, float the embedded profile. Every
// dot product accumulates in Sample, so the float profile contains no
// double arithmetic and runs natively on a single-precision FPU (the
// Cortex-M33 of the RP2350 among them). Fixed-point profiles are not
// provided here: a Q-format design for these layers is a documented
// per-primitive decision (CLAUDE.md), not a template instantiation.
//
// Real-time safety: geometry and weights are fixed at construction (the
// weight vectors are moved in and owned, so a copied layer is a deep copy
// that stays valid), every buffer is allocated there, and apply()/step()
// are noexcept and allocation-free.
#pragma once

#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

namespace tap::dsp::nn {

    /// Bumped when any formula or layout above changes; a trained model
    /// records the version its runtime was built against.
    inline constexpr std::uint32_t k_contract_version = 1;

    /// Element-wise nonlinearity of a dense layer.
    enum class activation : std::uint8_t {
        linear,  ///< identity
        tanh,    ///< std::tanh
        sigmoid, ///< 1 / (1 + exp(-v))
    };

    /// The logistic function in Sample arithmetic (the GRU gate form).
    template <typename Sample>
    inline Sample sigmoid(Sample v) noexcept {
        return Sample(1) / (Sample(1) + std::exp(-v));
    }

    /// Fully connected layer, y = act(W x + b), W row-major [out x in].
    template <typename Sample>
    class basic_dense {
      public:
        /// @pre w.size() == out * in, b.size() == out, out >= 1, in >= 1.
        basic_dense(std::vector<float> w, std::vector<float> b, std::size_t out, std::size_t in,
                    activation act = activation::linear)
            : m_w(std::move(w))
            , m_b(std::move(b))
            , m_out(out)
            , m_in(in)
            , m_act(act) {
            assert(out >= 1 && in >= 1);
            assert(m_w.size() == out * in && m_b.size() == out);
        }

        /// y[0..out) from x[0..in); the two must not alias.
        void apply(const Sample* x, Sample* y) const noexcept {
            const float* w = m_w.data();
            for (std::size_t i = 0; i < m_out; ++i) {
                Sample acc = static_cast<Sample>(m_b[i]);
                for (std::size_t j = 0; j < m_in; ++j) {
                    acc += static_cast<Sample>(w[i * m_in + j]) * x[j];
                }
                switch (m_act) {
                case activation::tanh:
                    y[i] = std::tanh(acc);
                    break;
                case activation::sigmoid:
                    y[i] = sigmoid(acc);
                    break;
                case activation::linear:
                default:
                    y[i] = acc;
                    break;
                }
            }
        }

        std::size_t               in() const noexcept { return m_in; }
        std::size_t               out() const noexcept { return m_out; }
        activation                act() const noexcept { return m_act; }
        const std::vector<float>& weights() const noexcept { return m_w; }
        const std::vector<float>& bias() const noexcept { return m_b; }

      private:
        std::vector<float> m_w;
        std::vector<float> m_b;
        std::size_t        m_out;
        std::size_t        m_in;
        activation         m_act;
    };

    /// Gated recurrent unit cell, PyTorch nn.GRU convention (gates r, z, n).
    template <typename Sample>
    class basic_gru {
      public:
        /// @pre w_ih.size() == 3*hidden*in, w_hh.size() == 3*hidden*hidden,
        ///      b_ih.size() == b_hh.size() == 3*hidden, hidden >= 1, in >= 1.
        basic_gru(std::vector<float> w_ih, std::vector<float> w_hh, std::vector<float> b_ih, std::vector<float> b_hh,
                  std::size_t hidden, std::size_t in)
            : m_w_ih(std::move(w_ih))
            , m_w_hh(std::move(w_hh))
            , m_b_ih(std::move(b_ih))
            , m_b_hh(std::move(b_hh))
            , m_hidden(hidden)
            , m_in(in)
            , m_g_ih(3 * hidden, Sample(0))
            , m_g_hh(3 * hidden, Sample(0))
            , m_state(hidden, Sample(0)) {
            assert(hidden >= 1 && in >= 1);
            assert(m_w_ih.size() == 3 * hidden * in && m_w_hh.size() == 3 * hidden * hidden);
            assert(m_b_ih.size() == 3 * hidden && m_b_hh.size() == 3 * hidden);
        }

        /// h = 0.
        void reset() noexcept {
            for (auto& h : m_state) {
                h = Sample(0);
            }
        }

        /// One time step on x[0..in); the new state is readable via state().
        void step(const Sample* x) noexcept {
            const std::size_t nh = m_hidden;
            const std::size_t ni = m_in;
            const float*      wi = m_w_ih.data();
            const float*      wh = m_w_hh.data();
            for (std::size_t i = 0; i < 3 * nh; ++i) {
                Sample acc = static_cast<Sample>(m_b_ih[i]);
                for (std::size_t j = 0; j < ni; ++j) {
                    acc += static_cast<Sample>(wi[i * ni + j]) * x[j];
                }
                m_g_ih[i] = acc;
                acc       = static_cast<Sample>(m_b_hh[i]);
                for (std::size_t j = 0; j < nh; ++j) {
                    acc += static_cast<Sample>(wh[i * nh + j]) * m_state[j];
                }
                m_g_hh[i] = acc;
            }
            for (std::size_t i = 0; i < nh; ++i) {
                const Sample r = sigmoid(m_g_ih[i] + m_g_hh[i]);
                const Sample z = sigmoid(m_g_ih[nh + i] + m_g_hh[nh + i]);
                const Sample n = std::tanh(m_g_ih[2 * nh + i] + r * m_g_hh[2 * nh + i]);
                m_state[i]     = (Sample(1) - z) * n + z * m_state[i];
            }
        }

        const Sample* state() const noexcept { return m_state.data(); }
        std::size_t   hidden() const noexcept { return m_hidden; }
        std::size_t   in() const noexcept { return m_in; }

      private:
        std::vector<float>  m_w_ih;
        std::vector<float>  m_w_hh;
        std::vector<float>  m_b_ih;
        std::vector<float>  m_b_hh;
        std::size_t         m_hidden;
        std::size_t         m_in;
        std::vector<Sample> m_g_ih;
        std::vector<Sample> m_g_hh;
        std::vector<Sample> m_state;
    };

    using dense   = basic_dense<double>;
    using dense32 = basic_dense<float>;
    using gru     = basic_gru<double>;
    using gru32   = basic_gru<float>;

} // namespace tap::dsp::nn
