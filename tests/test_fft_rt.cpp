// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// THE REAL-TIME GUARD for tap::dsp::basic_real_fft (Part 9 of
// docs/audit-fft-and-code-smells.md names this file).
//
// fft.h promises that the transforms are noexcept and allocation-free once
// the object is constructed, so they are safe on an audio thread. This file
// makes both halves of that promise checkable:
//
//   - noexcept is a static_assert on the four transform entry points, on
//     all six instantiations -- float, double, and the Q15 / Q31 fixed-point
//     profiles under both scaling policies (the pattern test_nn.cpp uses);
//   - "allocation-free" is checked by REPLACING THE GLOBAL operator new /
//     operator delete in this translation unit with counting versions and
//     asserting the count does not move across a transform. The replacement
//     is program-wide (that is how replaceable allocation functions work),
//     so every window below is kept tight: nothing between the two reads of
//     the counter but the call under test, and every buffer allocated before
//     the first read.
//
// WHAT THE GUARD SEES AND WHAT IT DOES NOT. It counts C++ allocation
// functions only. An engine that allocates through malloc/calloc directly, or
// inside a vendor library (Apple's vDSP on the macOS leg, CMSIS on the M55),
// is invisible to it. For the split-radix engine (the double profile, and
// float without a backend, since Stage 2b) the claim is complete: the
// engine's tables are two std::vectors sized in the constructor, and nothing
// in a transform touches an allocator. For the backends the guard
// covers the wrapper's own code and nothing more; a malloc interposer
// (glibc's __libc_malloc, or DYLD_INTERPOSE) is the tool for the vendor
// layer and is out of scope here. The self-test CountsAVectorAllocation
// keeps the counter honest: if the replacement ever stopped being picked up,
// that test fails first.
//
// The FIRST call after construction is covered as well as a steady-state one:
// the vendored C built its trig and bit-reversal tables lazily on the first
// transform (Part 1, item F6), and that initialization had to be
// allocation-free too, because the first transform a consumer runs is very
// often on the audio thread already. The engine builds its tables in the
// constructor (Part 4; routed at Stage 2b); this test is indifferent to where
// it happens, only to what it allocates, so it keeps covering the first call.
//
// NOT covered, deliberately: the float-I/O convenience overloads on the
// double engine (forward(const float*, float*) / inverse(const float*,
// float*)). fft.h documents them as a setup-time path that allocates a
// staging buffer per call (Part 1, item F7) and Decision D5 deprecates them
// for one consumer cycle; guarding them would pin a behaviour the plan is
// removing. They are also not noexcept, consistently with that.
//
// Copy and copy-assignment are pinned to produce bit-identical output to the
// source in both directions and all six instantiations (for the fixed-point
// profiles that includes the returned exponent, and the Q15 profile's int32
// work buffer travelling with the copy). test_fft_backend.cpp's
// CopiesAgreeWithTheirSource covers the float FORWARD at the certified
// geometries while walking the heap (its purpose is vDSP's alignment
// dispatch); this one is the general value-semantics check that the engine's
// tables travel with the object, including through operator=.
//
// Stage 4 (Part 9): k_is_shareable is what the header says — asserted below
// per instantiation, together with the compiler-checked half (a shareable
// engine's transforms are callable on a const engine, a non-shareable one's
// are not; test_fft_engine.cpp holds the per-engine numbers). The std::span
// overloads Part 9 also listed do not exist: Stage 4 added no overload to
// the transform surface (the pointer forms are the contract every consumer
// holds), so there is nothing to equate and nothing is asserted about them.

#include <atomic>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <new>
#include <type_traits>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include "support/signals.h"
#include "tap/dsp/fft.h"
#include "tap/dsp/fft/split_radix.h"

// ----------------------------------------------------------------------------
// Counting replacements for the replaceable global allocation functions
// ([new.delete.single] and [new.delete.array]). All eight allocating forms
// (throwing / nothrow x plain / aligned x single / array) go through one
// counter; every deallocating form is matched so nothing falls back to the
// library default. Aligned forms are served by over-allocating and stashing
// the raw pointer just below the aligned block, which needs no platform
// aligned-malloc.
// ----------------------------------------------------------------------------
namespace {

    std::atomic<std::size_t> allocation_count{0};

    void* counted_malloc(std::size_t n) {
        allocation_count.fetch_add(1, std::memory_order_relaxed);
        return std::malloc(n == 0 ? 1 : n);
    }

    void* counted_aligned_malloc(std::size_t n, std::size_t alignment) {
        allocation_count.fetch_add(1, std::memory_order_relaxed);
        const std::size_t slack = alignment + sizeof(void*);
        void* const       raw   = std::malloc(n + slack);
        if (raw == nullptr) {
            return nullptr;
        }
        const std::uintptr_t base    = reinterpret_cast<std::uintptr_t>(raw) + sizeof(void*);
        const std::uintptr_t aligned = (base + alignment - 1) & ~(static_cast<std::uintptr_t>(alignment) - 1);
        void* const          p       = reinterpret_cast<void*>(aligned);
        std::memcpy(reinterpret_cast<void**>(aligned) - 1, &raw, sizeof(void*));
        return p;
    }

    void counted_aligned_free(void* p) noexcept {
        if (p == nullptr) {
            return;
        }
        void* raw = nullptr;
        std::memcpy(&raw, static_cast<void**>(p) - 1, sizeof(void*));
        std::free(raw);
    }

} // namespace

// GCC pairs the `new` expressions it inlines in this TU with the replaced
// operator delete below, sees std::free on the result of an "allocation
// function" it treats as opaque, and reports -Wmismatched-new-delete once per
// inlining site (GCC bug 101480: replaced allocation functions that forward to
// malloc/free). The pairing is correct by construction here — every allocating
// form returns std::malloc's result and every deallocating form frees it — so
// the diagnostic is silenced around the replacement functions only.
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmismatched-new-delete"
#endif

void* operator new(std::size_t n) {
    void* p = counted_malloc(n);
    if (p == nullptr) {
        throw std::bad_alloc();
    }
    return p;
}
void* operator new[](std::size_t n) {
    void* p = counted_malloc(n);
    if (p == nullptr) {
        throw std::bad_alloc();
    }
    return p;
}
void* operator new(std::size_t n, const std::nothrow_t&) noexcept {
    return counted_malloc(n);
}
void* operator new[](std::size_t n, const std::nothrow_t&) noexcept {
    return counted_malloc(n);
}
void* operator new(std::size_t n, std::align_val_t al) {
    void* p = counted_aligned_malloc(n, static_cast<std::size_t>(al));
    if (p == nullptr) {
        throw std::bad_alloc();
    }
    return p;
}
void* operator new[](std::size_t n, std::align_val_t al) {
    void* p = counted_aligned_malloc(n, static_cast<std::size_t>(al));
    if (p == nullptr) {
        throw std::bad_alloc();
    }
    return p;
}
void* operator new(std::size_t n, std::align_val_t al, const std::nothrow_t&) noexcept {
    return counted_aligned_malloc(n, static_cast<std::size_t>(al));
}
void* operator new[](std::size_t n, std::align_val_t al, const std::nothrow_t&) noexcept {
    return counted_aligned_malloc(n, static_cast<std::size_t>(al));
}

void operator delete(void* p) noexcept {
    std::free(p);
}
void operator delete[](void* p) noexcept {
    std::free(p);
}
void operator delete(void* p, std::size_t) noexcept {
    std::free(p);
}
void operator delete[](void* p, std::size_t) noexcept {
    std::free(p);
}
void operator delete(void* p, const std::nothrow_t&) noexcept {
    std::free(p);
}
void operator delete[](void* p, const std::nothrow_t&) noexcept {
    std::free(p);
}
void operator delete(void* p, std::align_val_t) noexcept {
    counted_aligned_free(p);
}
void operator delete[](void* p, std::align_val_t) noexcept {
    counted_aligned_free(p);
}
void operator delete(void* p, std::size_t, std::align_val_t) noexcept {
    counted_aligned_free(p);
}
void operator delete[](void* p, std::size_t, std::align_val_t) noexcept {
    counted_aligned_free(p);
}
void operator delete(void* p, std::align_val_t, const std::nothrow_t&) noexcept {
    counted_aligned_free(p);
}
void operator delete[](void* p, std::align_val_t, const std::nothrow_t&) noexcept {
    counted_aligned_free(p);
}

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic pop
#endif

namespace {

    // The suite is typed over the FFT class itself, so the two scaling
    // policies of each fixed-point profile are distinct rows (the second
    // argument is the engine for the floating rows and the policy for the
    // fixed-point ones; fft.h, "Two parameter lists in one").
    template <typename Fft>
    struct sample_of;
    template <typename Sample, typename Policy>
    struct sample_of<tap::dsp::basic_real_fft<Sample, Policy>> {
        using type = Sample;
    };
    template <typename Fft>
    using sample_of_t = typename sample_of<Fft>::type;

    // ------------------------------------------------------------------------
    // noexcept, as a compile-time fact on every instantiation.
    // ------------------------------------------------------------------------
    template <typename Fft>
    constexpr bool transforms_are_noexcept() {
        using fft    = Fft;
        using sample = sample_of_t<Fft>;
        static_assert(noexcept(std::declval<fft&>().forward_inplace(std::declval<sample*>())));
        static_assert(noexcept(std::declval<fft&>().inverse_inplace(std::declval<sample*>())));
        static_assert(noexcept(std::declval<fft&>().forward(std::declval<const sample*>(), std::declval<sample*>())));
        static_assert(noexcept(std::declval<fft&>().inverse(std::declval<const sample*>(), std::declval<sample*>())));
        static_assert(noexcept(std::declval<const fft&>().size()));
        static_assert(noexcept(std::declval<const fft&>().num_bins()));
        return true;
    }
    static_assert(transforms_are_noexcept<tap::dsp::real_fft32>());
    static_assert(transforms_are_noexcept<tap::dsp::real_fft>());
    static_assert(transforms_are_noexcept<tap::dsp::real_fft_q15>());
    static_assert(transforms_are_noexcept<tap::dsp::real_fft_q31>());
    static_assert(transforms_are_noexcept<tap::dsp::real_fft_q15_bfp>());
    static_assert(transforms_are_noexcept<tap::dsp::real_fft_q31_bfp>());
    // The fixed-point profiles' constant exponent is a noexcept constexpr too.
    static_assert(noexcept(tap::dsp::real_fft_q15::fixed_scaling_exponent(4)));
    static_assert(noexcept(tap::dsp::real_fft_q31_bfp::fixed_scaling_exponent(4)));

    // ------------------------------------------------------------------------
    // k_is_shareable, per instantiation, as the header states it (Stage 4):
    // the split-radix engine and Q31 true, the two scratch-carrying
    // accelerated engines and Q15 false. The float default's value follows
    // the engine the build selected, so it is derived, not spelled.
    // ------------------------------------------------------------------------
    template <typename Fft>
    constexpr bool expected_shareable() {
        using sample = sample_of_t<Fft>;
        if constexpr (std::is_same_v<sample, double>) {
            return true;
        }
        else if constexpr (std::is_same_v<sample, float>) {
            return std::is_same_v<typename Fft::engine, tap::dsp::detail::split_radix_rdft<float>>;
        }
        else {
            return std::is_same_v<sample, std::int32_t>; // Q31 true, Q15 false
        }
    }

    /// The compiler-checked half of the trait: the engine's transforms are
    /// callable on a const engine exactly when the trait says shareable.
    template <typename Fft>
    constexpr bool k_engine_transforms_are_const = requires(const typename Fft::engine& e, sample_of_t<Fft>* a) {
        e.forward_inplace(a);
        e.inverse_inplace(a);
    };

    template <typename Fft>
    constexpr bool shareability_is_the_headers_number() {
        static_assert(std::is_same_v<decltype(Fft::k_is_shareable), const bool>);
        static_assert(Fft::k_is_shareable == expected_shareable<Fft>());
        static_assert(k_engine_transforms_are_const<Fft> == Fft::k_is_shareable);
        return true;
    }
    static_assert(shareability_is_the_headers_number<tap::dsp::real_fft32>());
    static_assert(shareability_is_the_headers_number<tap::dsp::real_fft>());
    static_assert(shareability_is_the_headers_number<tap::dsp::real_fft_q15>());
    static_assert(shareability_is_the_headers_number<tap::dsp::real_fft_q31>());
    static_assert(shareability_is_the_headers_number<tap::dsp::real_fft_q15_bfp>());
    static_assert(shareability_is_the_headers_number<tap::dsp::real_fft_q31_bfp>());

    // ------------------------------------------------------------------------
    // The allocation guard.
    // ------------------------------------------------------------------------
    class allocation_guard {
      public:
        allocation_guard() noexcept
            : m_start(allocation_count.load(std::memory_order_relaxed)) {}
        std::size_t allocations_since() const noexcept {
            return allocation_count.load(std::memory_order_relaxed) - m_start;
        }

      private:
        std::size_t m_start;
    };

    constexpr std::size_t k_guarded_sizes[] = {512, 4096};

    template <typename Fft>
    class fft_rt_test : public ::testing::Test {};

    using fft_types = ::testing::Types<tap::dsp::real_fft32, tap::dsp::real_fft, tap::dsp::real_fft_q15,
                                       tap::dsp::real_fft_q31, tap::dsp::real_fft_q15_bfp, tap::dsp::real_fft_q31_bfp>;
    TYPED_TEST_SUITE(fft_rt_test, fft_types);

    // Already proved by the namespace-scope static_asserts above; these exist
    // so each promise has a row in the test listing per profile.
    TYPED_TEST(fft_rt_test, TransformsAreNoexcept) {
        EXPECT_TRUE(transforms_are_noexcept<TypeParam>());
    }

    TYPED_TEST(fft_rt_test, ShareabilityIsTheHeadersNumber) {
        EXPECT_TRUE(shareability_is_the_headers_number<TypeParam>());
        EXPECT_EQ(TypeParam::k_is_shareable, expected_shareable<TypeParam>());
    }

    // The guard itself must see allocations, or every test below passes for
    // the wrong reason.
    TEST(fft_rt_guard, CountsAVectorAllocation) {
        allocation_guard guard;
        {
            std::vector<double> v(64);
            EXPECT_NE(v.data(), nullptr);
        }
        EXPECT_GE(guard.allocations_since(), 1u);
    }

    template <typename Fft, typename Call>
    void expect_no_allocation(std::size_t n, const char* what, Call&& call) {
        using sample = sample_of_t<Fft>;
        Fft                 fft(n);
        std::vector<sample> in  = tap::dsp::test::random_signal<sample>(n, 0x9E3779B9u);
        std::vector<sample> out = in;

        // First call after construction (lazy table init today, F6).
        {
            allocation_guard guard;
            call(fft, in.data(), out.data());
            const std::size_t count = guard.allocations_since();
            EXPECT_EQ(count, 0u) << what << " N=" << n << " allocated " << count
                                 << " time(s) on the FIRST call after construction";
        }
        // Steady state.
        {
            allocation_guard guard;
            call(fft, in.data(), out.data());
            const std::size_t count = guard.allocations_since();
            EXPECT_EQ(count, 0u) << what << " N=" << n << " allocated " << count << " time(s) in steady state";
        }
    }

    TYPED_TEST(fft_rt_test, ForwardInplaceAllocatesNothing) {
        using sample = sample_of_t<TypeParam>;
        for (const std::size_t n : k_guarded_sizes) {
            expect_no_allocation<TypeParam>(
                n, "forward_inplace", [](TypeParam& fft, sample*, sample* out) { (void)fft.forward_inplace(out); });
        }
    }

    TYPED_TEST(fft_rt_test, InverseInplaceAllocatesNothing) {
        using sample = sample_of_t<TypeParam>;
        for (const std::size_t n : k_guarded_sizes) {
            expect_no_allocation<TypeParam>(
                n, "inverse_inplace", [](TypeParam& fft, sample*, sample* out) { (void)fft.inverse_inplace(out); });
        }
    }

    TYPED_TEST(fft_rt_test, ForwardOutOfPlaceAllocatesNothing) {
        using sample = sample_of_t<TypeParam>;
        for (const std::size_t n : k_guarded_sizes) {
            expect_no_allocation<TypeParam>(
                n, "forward", [](TypeParam& fft, const sample* in, sample* out) { (void)fft.forward(in, out); });
        }
    }

    TYPED_TEST(fft_rt_test, InverseOutOfPlaceAllocatesNothing) {
        using sample = sample_of_t<TypeParam>;
        for (const std::size_t n : k_guarded_sizes) {
            expect_no_allocation<TypeParam>(
                n, "inverse", [](TypeParam& fft, const sample* in, sample* out) { (void)fft.inverse(in, out); });
        }
    }

    // ------------------------------------------------------------------------
    // Value semantics: a copy computes exactly what its source computes.
    // ------------------------------------------------------------------------
    template <typename Sample>
    void expect_bit_identical(const std::vector<Sample>& a, const std::vector<Sample>& b, const char* what) {
        ASSERT_EQ(a.size(), b.size());
        for (std::size_t i = 0; i < a.size(); ++i) {
            // memcmp, not ==: the comparison is over the exact bits.
            ASSERT_EQ(std::memcmp(&a[i], &b[i], sizeof(Sample)), 0)
                << what << " differs from its source at index " << i;
        }
    }

    template <typename Sample>
    struct both_directions {
        std::vector<Sample> spectrum;
        std::vector<Sample> time;
        int                 forward_exponent = 0; ///< fixed-point profiles; 0 for float/double
        int                 inverse_exponent = 0;
    };

    /// The returned exponent, or 0 where the transform returns void.
    template <typename Call>
    int exponent_of(Call&& call) {
        if constexpr (std::is_void_v<decltype(call())>) {
            call();
            return 0;
        }
        else {
            return call();
        }
    }

    template <typename Fft>
    both_directions<sample_of_t<Fft>> run_both(Fft& fft, const std::vector<sample_of_t<Fft>>& x) {
        both_directions<sample_of_t<Fft>> r;
        r.spectrum         = x;
        r.forward_exponent = exponent_of([&] { return fft.forward_inplace(r.spectrum.data()); });
        r.time             = r.spectrum;
        r.inverse_exponent = exponent_of([&] { return fft.inverse_inplace(r.time.data()); });
        return r;
    }

    template <typename Sample>
    void expect_same_result(const both_directions<Sample>& a, const both_directions<Sample>& b, const char* what) {
        expect_bit_identical(a.spectrum, b.spectrum, what);
        expect_bit_identical(a.time, b.time, what);
        EXPECT_EQ(a.forward_exponent, b.forward_exponent) << what << ": forward exponent differs";
        EXPECT_EQ(a.inverse_exponent, b.inverse_exponent) << what << ": inverse exponent differs";
    }

    TYPED_TEST(fft_rt_test, CopyProducesBitIdenticalOutput) {
        using sample            = sample_of_t<TypeParam>;
        constexpr std::size_t n = 512;
        const auto            x = tap::dsp::test::random_signal<sample>(n, 0x2545F491u);

        TypeParam original(n);
        // Source already warmed (tables built) so the copy carries built tables.
        const auto from_original = run_both(original, x);

        TypeParam  copy(original);
        const auto from_copy = run_both(copy, x);
        expect_same_result(from_copy, from_original, "copy");

        // And the source is unaffected by having been copied.
        const auto again = run_both(original, x);
        expect_same_result(again, from_original, "source after copy");
    }

    TYPED_TEST(fft_rt_test, CopyOfUnwarmedSourceProducesBitIdenticalOutput) {
        using sample            = sample_of_t<TypeParam>;
        constexpr std::size_t n = 512;
        const auto            x = tap::dsp::test::random_signal<sample>(n, 0x2545F491u);

        // Copied BEFORE any transform: with lazy tables (F6) both objects
        // build their own; with constructor-built tables both carry the same.
        TypeParam  original(n);
        TypeParam  copy(original);
        const auto from_copy     = run_both(copy, x);
        const auto from_original = run_both(original, x);
        expect_same_result(from_copy, from_original, "unwarmed copy");
    }

    TYPED_TEST(fft_rt_test, CopyAssignmentProducesBitIdenticalOutput) {
        using sample            = sample_of_t<TypeParam>;
        constexpr std::size_t n = 512;
        const auto            x = tap::dsp::test::random_signal<sample>(n, 0x2545F491u);

        TypeParam  original(n);
        const auto from_original = run_both(original, x);

        // Assigned over an engine of a DIFFERENT geometry, so the assignment
        // has to replace every table (and the Q15 work buffer), not just
        // refresh one of the same size.
        TypeParam target(4096);
        target = original;
        EXPECT_EQ(target.size(), n);
        EXPECT_EQ(target.num_bins(), n / 2 + 1);
        const auto from_target = run_both(target, x);
        expect_same_result(from_target, from_original, "assigned");
    }

} // namespace
