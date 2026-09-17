/// @file spectrum.h
/// @brief Non-owning view over a DspTap packed real spectrum.
// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// One home for the bin arithmetic that every consumer of basic_real_fft used
// to re-derive by hand (data[0] is DC, data[1] is Nyquist, data[2k]/[2k+1]
// are bin k). The view adds no state beyond the pointer and the size, no
// virtuals and no allocation; every accessor is constexpr, noexcept and
// inlines to the same index expression the hand-written code used, so a
// migration onto it is bit-identical by construction.

#pragma once

#include <cassert>
#include <complex>
#include <cstddef>
#include <type_traits>

namespace tap::dsp {

    /// Non-owning view over the packed spectrum that basic_real_fft's forward
    /// transform leaves in place and its inverse consumes.
    ///
    /// The packing, as numbers, for a transform of N real samples (N even;
    /// basic_real_fft requires a power of two >= 4). N/2 + 1 bins live in N
    /// values a[0..N):
    ///   - bin[0]   = a[0]                (DC;      its imaginary part is zero and not stored)
    ///   - bin[N/2] = a[1]                (Nyquist; its imaginary part is zero and not stored)
    ///   - bin[k]   = a[2k] + i * a[2k+1]  for 1 <= k < N/2
    ///
    /// Sign convention: bin[k] = sum_j x[j] * W^(jk) with W = exp(+2*pi*i/N),
    /// so the imaginary parts are CONJUGATED relative to the engineering DFT
    /// (exp(-2*pi*i/N)). The inverse transform is unnormalized: an in-place
    /// round trip needs a 2/N scaling.
    ///
    /// The native accessors (dc(), nyquist(), re(k), im(k), power(k)) read
    /// the spectrum exactly as stored and are the primary interface; spectral
    /// products between spectra from this library need no conjugation.
    /// bin_engineering(k) is the one convention-flipping accessor, for code
    /// that applies textbook phase formulas verbatim; it says so in its name.
    ///
    /// `Sample` may be const-qualified: packed_spectrum<const float> is a
    /// read-only view, packed_spectrum<float> also writes through re()/im()/
    /// dc()/nyquist(). A mutable view converts implicitly to the const one.
    ///
    /// @pre data points at N values; N is even and >= 2. Bin indices are
    ///      asserted in debug builds and unchecked in release, as everywhere
    ///      in the Tap libraries.
    template <typename Sample>
    class packed_spectrum {
      public:
        using value_type = std::remove_const_t<Sample>;

        static_assert(std::is_same_v<value_type, float> || std::is_same_v<value_type, double>,
                      "packed_spectrum views the two basic_real_fft profiles: float and double");

        /// View N packed values at data.
        constexpr packed_spectrum(Sample* data, std::size_t n) noexcept
            : m_data(data)
            , m_n(n) {
            assert(data != nullptr && n >= 2 && (n % 2) == 0);
        }

        /// A mutable view converts to the read-only view over the same buffer.
        template <typename Other>
            requires(std::is_const_v<Sample> && std::is_same_v<Other, value_type>)
        constexpr packed_spectrum(const packed_spectrum<Other>& other) noexcept
            : m_data(other.data())
            , m_n(other.size()) {}

        /// Transform size N: the number of packed values.
        [[nodiscard]] constexpr std::size_t size() const noexcept { return m_n; }
        /// N/2 + 1: the number of distinct bins, DC and Nyquist included.
        [[nodiscard]] constexpr std::size_t num_bins() const noexcept { return m_n / 2 + 1; }
        /// The packed buffer itself, for handing to a transform.
        [[nodiscard]] constexpr Sample* data() const noexcept { return m_data; }

        /// bin[0], the DC term (a[0]); real by construction.
        [[nodiscard]] constexpr Sample& dc() const noexcept { return m_data[0]; }
        /// bin[N/2], the Nyquist term (a[1]); real by construction.
        [[nodiscard]] constexpr Sample& nyquist() const noexcept { return m_data[1]; }

        /// Real part of bin k (a[2k]).
        /// @pre 1 <= k < N/2 — DC and Nyquist have their own accessors.
        [[nodiscard]] constexpr Sample& re(std::size_t k) const noexcept {
            assert(k >= 1 && k < m_n / 2);
            return m_data[2 * k];
        }
        /// Imaginary part of bin k (a[2k+1]), in the native exp(+i) convention.
        /// @pre 1 <= k < N/2 — DC and Nyquist have their own accessors.
        [[nodiscard]] constexpr Sample& im(std::size_t k) const noexcept {
            assert(k >= 1 && k < m_n / 2);
            return m_data[2 * k + 1];
        }

        /// |bin[k]|^2 for any bin, DC and Nyquist included, computed in Sample
        /// precision as re*re + im*im (in that order; consumers that pinned
        /// the hand-written product keep their bits).
        /// @pre 0 <= k <= N/2.
        [[nodiscard]] constexpr value_type power(std::size_t k) const noexcept {
            assert(k <= m_n / 2);
            if (k == 0) {
                return m_data[0] * m_data[0];
            }
            if (k == m_n / 2) {
                return m_data[1] * m_data[1];
            }
            const value_type re = m_data[2 * k];
            const value_type im = m_data[2 * k + 1];
            return re * re + im * im;
        }

        /// Bin k as a complex number in the ENGINEERING convention
        /// (exp(-2*pi*i/N)): this CONJUGATES the stored value, returning
        /// a[2k] - i * a[2k+1], so textbook phase-vocoder formulas apply
        /// verbatim. Not the native reading of the spectrum; a value written
        /// back must be conjugated again (see pvoc.h's synthesis pack).
        /// @pre 1 <= k < N/2.
        [[nodiscard]] constexpr std::complex<value_type> bin_engineering(std::size_t k) const noexcept {
            assert(k >= 1 && k < m_n / 2);
            return std::complex<value_type>(m_data[2 * k], -m_data[2 * k + 1]);
        }

      private:
        Sample*     m_data;
        std::size_t m_n;
    };

    template <typename Sample>
    packed_spectrum(Sample*, std::size_t) -> packed_spectrum<Sample>;

} // namespace tap::dsp
