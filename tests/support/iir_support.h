/// @file iir_support.h
/// @brief The recursive batteries' support: the host's flush-to-zero, direct-form golden references, stream drivers.
// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// - The host's job (subnormals.h): set_flush_to_zero() writes the control
//   bits a host sets on its audio thread, per target, and the fixtures set
//   them in SetUp and restore the caller's register in TearDown, so the rest
//   of the battery on the same thread is untouched (IIR_PLAN I9's tests).
//   Writing the register is the tests' business only; tap/dsp never does.
// - The golden legs, independent of the polyphase structure: each transfer
//   function run in direct form at the higher rate, in double, around
//   zero-stuffing or down-sampling (IIR_PLAN section 4's golden leg, with
//   the structures written out here rather than through scipy).
// - Stream drivers: a mono commutator over the step kernels, the minimal
//   stage machine the kernels leave to their caller (iir_kernels.h).

#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

#include <gtest/gtest.h>

#include "tap/dsp/iir_design.h"
#include "tap/dsp/iir_kernels.h"
#include "tap/dsp/subnormals.h"

namespace tap::dsp::test {

    // ------------------------------------------------------------------------
    // The host's flush-to-zero
    // ------------------------------------------------------------------------

    /// Writes the calling thread's floating-point control register (where the
    /// target has one).
    inline void write_floating_point_control(std::uint64_t v) noexcept {
#if defined(TAP_DSP_FP_CONTROL_X86_64)
        _mm_setcsr(static_cast<unsigned>(v));
#elif defined(TAP_DSP_FP_CONTROL_AARCH64)
#if defined(_MSC_VER) && !defined(__clang__)
        _WriteStatusReg(ARM64_FPCR, static_cast<__int64>(v));
#else
        __asm__ volatile("msr fpcr, %0" : : "r"(v));
#endif
#elif defined(TAP_DSP_FP_CONTROL_AARCH32)
        __builtin_arm_set_fpscr(static_cast<unsigned>(v));
#else
        static_cast<void>(v);
#endif
    }

    /// Sets (or clears) the host's flush-to-zero bits; returns the register
    /// as it was.
    inline std::uint64_t set_flush_to_zero(bool on) noexcept {
        const std::uint64_t saved = floating_point_control();
        write_floating_point_control(on ? (saved | k_flush_to_zero_bits) : (saved & ~k_flush_to_zero_bits));
        return saved;
    }

    /// The control register's bits that are not cumulative exception flags
    /// (which an arithmetic call may legitimately raise): what "the kernel
    /// never touches the environment" compares.
    inline constexpr std::uint64_t k_control_mask =
#if defined(TAP_DSP_FP_CONTROL_X86_64)
        ~std::uint64_t{0x3f}; // MXCSR flags IE DE ZE OE UE PE
#elif defined(TAP_DSP_FP_CONTROL_AARCH64)
        ~std::uint64_t{0}; // FPCR holds no flags (they are FPSR's)
#elif defined(TAP_DSP_FP_CONTROL_AARCH32)
        ~std::uint64_t{0xF800009Fu}; // FPSCR N Z C V QC and IDC IXC UFC OFC DZC IOC
#else
        ~std::uint64_t{0};
#endif

    /// The fixture of every recursive test that runs as a host would: FTZ
    /// (and DAZ on x86) set in SetUp, the caller's register restored in
    /// TearDown.
    class flush_to_zero_fixture : public ::testing::Test {
      protected:
        void SetUp() override { m_saved = set_flush_to_zero(true); }
        void TearDown() override { write_floating_point_control(m_saved); }

      private:
        std::uint64_t m_saved = 0;
    };

    /// The same with the bits cleared: the honest limit's host.
    class no_flush_to_zero_fixture : public ::testing::Test {
      protected:
        void SetUp() override { m_saved = set_flush_to_zero(false); }
        void TearDown() override { write_floating_point_control(m_saved); }

      private:
        std::uint64_t m_saved = 0;
    };

    // ------------------------------------------------------------------------
    // Golden legs: direct form at the higher rate, double
    // ------------------------------------------------------------------------

    /// (a + z^-M) / (1 + a z^-M), direct form, in place.
    inline void allpass_in_zm(std::vector<double>& x, double a, std::size_t m) {
        std::vector<double> y(x.size());
        for (std::size_t n = 0; n < x.size(); ++n) {
            const double xm = n >= m ? x[n - m] : 0.0;
            const double ym = n >= m ? y[n - m] : 0.0;
            y[n]            = a * x[n] + xm - a * ym;
        }
        x.swap(y);
    }

    /// (c2 + c1 z^-M + z^-2M) / (1 + c1 z^-M + c2 z^-2M), direct form, in place.
    inline void allpass2_in_zm(std::vector<double>& x, double c1, double c2, std::size_t m) {
        std::vector<double> y(x.size());
        for (std::size_t n = 0; n < x.size(); ++n) {
            const double x1 = n >= m ? x[n - m] : 0.0;
            const double x2 = n >= 2 * m ? x[n - 2 * m] : 0.0;
            const double y1 = n >= m ? y[n - m] : 0.0;
            const double y2 = n >= 2 * m ? y[n - 2 * m] : 0.0;
            y[n]            = c2 * x[n] + c1 * x1 + x2 - c1 * y1 - c2 * y2;
        }
        x.swap(y);
    }

    /// z^-d, in place.
    inline void delay_by(std::vector<double>& x, std::size_t d) {
        std::vector<double> y(x.size(), 0.0);
        for (std::size_t n = d; n < x.size(); ++n) {
            y[n] = x[n - d];
        }
        x.swap(y);
    }

    /// The nonlinear half-band H(z) = 1/2 [A0(z^2) + z^-1 A1(z^2)] at the
    /// higher rate, direct form.
    inline std::vector<double> halfband_direct(const nonlinear_halfband_design& d, std::vector<double> x) {
        std::vector<double> p0 = x;
        std::vector<double> p1 = std::move(x);
        for (std::size_t i = 0; i < d.path_size(0); ++i) {
            allpass_in_zm(p0, d.path_coefficient(0, i), 2);
        }
        delay_by(p1, 1);
        for (std::size_t i = 0; i < d.path_size(1); ++i) {
            allpass_in_zm(p1, d.path_coefficient(1, i), 2);
        }
        for (std::size_t n = 0; n < p0.size(); ++n) {
            p0[n] = 0.5 * (p0[n] + p1[n]);
        }
        return p0;
    }

    /// The near-linear half-band 1/2 [z^-2K + z^-1 A1(z^2)], direct form.
    inline std::vector<double> halfband_direct(const near_linear_halfband_design& d, std::vector<double> x) {
        std::vector<double> p0 = x;
        std::vector<double> p1 = std::move(x);
        delay_by(p0, 2 * static_cast<std::size_t>(d.delay));
        delay_by(p1, 1);
        std::size_t i = 0;
        for (int left = d.order; left > 0;) {
            if (left >= 2) {
                allpass2_in_zm(p1, d.sections[i], d.sections[i + 1], 2);
                i += 2;
                left -= 2;
            }
            else {
                allpass_in_zm(p1, d.sections[i], 2);
                i += 1;
                left -= 1;
            }
        }
        for (std::size_t n = 0; n < p0.size(); ++n) {
            p0[n] = 0.5 * (p0[n] + p1[n]);
        }
        return p0;
    }

    /// The four-path H(z) = 1/4 sum_k z^-k A_k(z^4), direct form.
    inline std::vector<double> fourpath_direct(const fourpath_design& d, const std::vector<double>& x) {
        std::vector<double> sum(x.size(), 0.0);
        for (std::size_t k = 0; k < 4; ++k) {
            std::vector<double> p = x;
            delay_by(p, k);
            for (const double a : d.paths[k]) {
                allpass_in_zm(p, a, 4);
            }
            for (std::size_t n = 0; n < x.size(); ++n) {
                sum[n] += p[n];
            }
        }
        for (double& v : sum) {
            v *= 0.25;
        }
        return sum;
    }

    /// A section cascade in direct form I, double.
    inline std::vector<double> sections_direct(std::span<const iir_section> s, std::vector<double> x) {
        for (const iir_section& c : s) {
            std::vector<double> y(x.size());
            double              x1 = 0.0;
            double              x2 = 0.0;
            double              y1 = 0.0;
            double              y2 = 0.0;
            for (std::size_t n = 0; n < x.size(); ++n) {
                y[n] = c.b0 * x[n] + c.b1 * x1 + c.b2 * x2 - c.a1 * y1 - c.a2 * y2;
                x2   = x1;
                x1   = x[n];
                y2   = y1;
                y1   = y[n];
            }
            x.swap(y);
        }
        return x;
    }

    /// Zero-stuffed by L: u[L n] = x[n], zeros between.
    inline std::vector<double> zero_stuffed(const std::vector<double>& x, std::size_t l) {
        std::vector<double> u(x.size() * l, 0.0);
        for (std::size_t n = 0; n < x.size(); ++n) {
            u[l * n] = x[n];
        }
        return u;
    }

    // ------------------------------------------------------------------------
    // Stream drivers over the step kernels (the caller's commutator)
    // ------------------------------------------------------------------------

    template <typename S>
    double to_unit(S v) {
        if constexpr (std::is_floating_point_v<S>) {
            return static_cast<double>(v);
        }
        else {
            return static_cast<double>(v) / static_cast<double>(std::int64_t{1} << std::numeric_limits<S>::digits);
        }
    }

    /// The nonlinear decimator over a whole even-length stream, mono:
    /// out[k] = y[k] from in[2k] and in[2k - 1].
    template <typename S>
    std::vector<S> run_halfband_decimate(const nonlinear_halfband_table<S>& t, const std::vector<S>& in) {
        std::vector<iir_state<S>> st(t.state_size(), iir_state<S>{});
        std::vector<S>            out(in.size() / 2);
        S                         odd{};
        for (std::size_t k = 0; k < out.size(); ++k) {
            out[k] = halfband_decimate(t, st.data(), in[2 * k], odd);
            odd    = in[2 * k + 1];
        }
        return out;
    }

    template <typename S>
    std::vector<S> run_halfband_interpolate(const nonlinear_halfband_table<S>& t, const std::vector<S>& in) {
        std::vector<iir_state<S>> st(t.state_size(), iir_state<S>{});
        std::vector<S>            out(2 * in.size());
        for (std::size_t n = 0; n < in.size(); ++n) {
            halfband_interpolate(t, st.data(), in[n], out.data() + 2 * n);
        }
        return out;
    }

    template <typename S>
    std::vector<S> run_near_linear_decimate(const near_linear_halfband_table<S>& t, const std::vector<S>& in) {
        std::vector<iir_state<S>> st(t.state_size(), iir_state<S>{});
        std::vector<S>            out(in.size() / 2);
        S                         odd{};
        for (std::size_t k = 0; k < out.size(); ++k) {
            const S delayed = k >= t.delay ? in[2 * (k - t.delay)] : S{};
            out[k]          = near_linear_decimate(t, st.data(), delayed, odd);
            odd             = in[2 * k + 1];
        }
        return out;
    }

    template <typename S>
    std::vector<S> run_near_linear_interpolate(const near_linear_halfband_table<S>& t, const std::vector<S>& in) {
        std::vector<iir_state<S>> st(t.state_size(), iir_state<S>{});
        std::vector<S>            out(2 * in.size());
        for (std::size_t n = 0; n < in.size(); ++n) {
            const S delayed = n >= t.delay ? in[n - t.delay] : S{};
            near_linear_interpolate(t, st.data(), delayed, in[n], out.data() + 2 * n);
        }
        return out;
    }

    template <typename S>
    std::vector<S> run_fourpath_decimate(const fourpath_table<S>& t, const std::vector<S>& in) {
        std::vector<iir_state<S>> st(t.state_size(), iir_state<S>{});
        std::vector<S>            out(in.size() / 4);
        for (std::size_t n = 0; n < out.size(); ++n) {
            std::array<S, 4> x{};
            for (std::size_t k = 0; k < 4; ++k) {
                x[k] = 4 * n >= k ? in[4 * n - k] : S{};
            }
            out[n] = fourpath_decimate(t, st.data(), x.data());
        }
        return out;
    }

    template <typename S>
    std::vector<S> run_fourpath_interpolate(const fourpath_table<S>& t, const std::vector<S>& in) {
        std::vector<iir_state<S>> st(t.state_size(), iir_state<S>{});
        std::vector<S>            out(4 * in.size());
        for (std::size_t n = 0; n < in.size(); ++n) {
            fourpath_interpolate(t, st.data(), in[n], out.data() + 4 * n);
        }
        return out;
    }

    template <typename S>
    std::vector<S> run_biquad_cascade(const biquad_cascade_table<S>& t, const std::vector<S>& in) {
        std::vector<iir_state<S>> st(t.state_size(), iir_state<S>{});
        std::vector<S>            out(in.size());
        for (std::size_t n = 0; n < in.size(); ++n) {
            out[n] = biquad_cascade(t, st.data(), in[n]);
        }
        return out;
    }

    /// FNV-1a-64 over the outputs' bytes (the fixed-point bit pins).
    template <typename S>
    std::uint64_t fnv1a64(const std::vector<S>& v) {
        std::uint64_t h = 0xcbf29ce484222325ull;
        for (const S s : v) {
            auto u = static_cast<std::uint64_t>(static_cast<std::make_unsigned_t<S>>(s));
            for (std::size_t b = 0; b < sizeof(S); ++b) {
                h ^= (u >> (8 * b)) & 0xffu;
                h *= 0x100000001b3ull;
            }
        }
        return h;
    }

} // namespace tap::dsp::test
