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
//     both instantiations (the pattern test_nn.cpp uses);
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
// is invisible to it. For today's Ooura path the claim is complete: makewt
// and makect write into the caller's preallocated ip/w tables (fftsg.c
// 655-756) and the C never calls an allocator. For the backends the guard
// covers the wrapper's own code and nothing more; a malloc interposer
// (glibc's __libc_malloc, or DYLD_INTERPOSE) is the tool for the vendor
// layer and is out of scope here. The self-test CountsAVectorAllocation
// keeps the counter honest: if the replacement ever stopped being picked up,
// that test fails first.
//
// The FIRST call after construction is covered as well as a steady-state one:
// today Ooura builds its trig and bit-reversal tables lazily on the first
// transform (Part 1, item F6), and that initialization must be allocation-free
// too, because the first transform a consumer runs is very often on the audio
// thread already. The port moves table construction into the constructor
// (Part 4); this test is indifferent to where it happens, only to what it
// allocates.
//
// NOT covered, deliberately: the float-I/O convenience overloads on the
// double engine (forward(const float*, float*) / inverse(const float*,
// float*)). fft.h documents them as a setup-time path that allocates a
// staging buffer per call (Part 1, item F7) and Decision D5 deprecates them
// for one consumer cycle; guarding them would pin a behaviour the plan is
// removing. They are also not noexcept, consistently with that.
//
// Copy and copy-assignment are pinned to produce bit-identical output to the
// source in both directions and both precisions. test_fft_backend.cpp's
// CopiesAgreeWithTheirSource covers the float FORWARD at the certified
// geometries while walking the heap (its purpose is vDSP's alignment
// dispatch); this one is the general value-semantics check that the engine's
// tables travel with the object, including through operator=.
//
// Stage 4 extension points, per Part 9: the std::span overloads equal the
// pointer overloads, and Engine::is_shareable is what the header says. Neither
// exists yet, so neither is asserted here.

#include <atomic>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <new>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include "support/signals.h"
#include "tap/dsp/fft.h"

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

    template <typename Sample>
    using fft_t = tap::dsp::basic_real_fft<Sample>;

    // ------------------------------------------------------------------------
    // noexcept, as a compile-time fact on both instantiations.
    // ------------------------------------------------------------------------
    template <typename Sample>
    constexpr bool transforms_are_noexcept() {
        using fft = fft_t<Sample>;
        static_assert(noexcept(std::declval<fft&>().forward_inplace(std::declval<Sample*>())));
        static_assert(noexcept(std::declval<fft&>().inverse_inplace(std::declval<Sample*>())));
        static_assert(noexcept(std::declval<fft&>().forward(std::declval<const Sample*>(), std::declval<Sample*>())));
        static_assert(noexcept(std::declval<fft&>().inverse(std::declval<const Sample*>(), std::declval<Sample*>())));
        static_assert(noexcept(std::declval<const fft&>().size()));
        static_assert(noexcept(std::declval<const fft&>().num_bins()));
        return true;
    }
    static_assert(transforms_are_noexcept<float>());
    static_assert(transforms_are_noexcept<double>());

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

    template <typename Sample>
    class fft_rt_test : public ::testing::Test {};

    using sample_types = ::testing::Types<float, double>;
    TYPED_TEST_SUITE(fft_rt_test, sample_types);

    // Already proved by the namespace-scope static_asserts above; this exists
    // so the promise has a row in the test listing per profile.
    TYPED_TEST(fft_rt_test, TransformsAreNoexcept) {
        EXPECT_TRUE(transforms_are_noexcept<TypeParam>());
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

    template <typename Sample, typename Call>
    void expect_no_allocation(std::size_t n, const char* what, Call&& call) {
        fft_t<Sample>       fft(n);
        std::vector<Sample> in  = tap::dsp::test::random_signal<Sample>(n, 0x9E3779B9u);
        std::vector<Sample> out = in;

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
        for (const std::size_t n : k_guarded_sizes) {
            expect_no_allocation<TypeParam>(
                n, "forward_inplace",
                [](fft_t<TypeParam>& fft, TypeParam*, TypeParam* out) { fft.forward_inplace(out); });
        }
    }

    TYPED_TEST(fft_rt_test, InverseInplaceAllocatesNothing) {
        for (const std::size_t n : k_guarded_sizes) {
            expect_no_allocation<TypeParam>(
                n, "inverse_inplace",
                [](fft_t<TypeParam>& fft, TypeParam*, TypeParam* out) { fft.inverse_inplace(out); });
        }
    }

    TYPED_TEST(fft_rt_test, ForwardOutOfPlaceAllocatesNothing) {
        for (const std::size_t n : k_guarded_sizes) {
            expect_no_allocation<TypeParam>(
                n, "forward", [](fft_t<TypeParam>& fft, const TypeParam* in, TypeParam* out) { fft.forward(in, out); });
        }
    }

    TYPED_TEST(fft_rt_test, InverseOutOfPlaceAllocatesNothing) {
        for (const std::size_t n : k_guarded_sizes) {
            expect_no_allocation<TypeParam>(
                n, "inverse", [](fft_t<TypeParam>& fft, const TypeParam* in, TypeParam* out) { fft.inverse(in, out); });
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
    };

    template <typename Sample>
    both_directions<Sample> run_both(fft_t<Sample>& fft, const std::vector<Sample>& x) {
        both_directions<Sample> r;
        r.spectrum = x;
        fft.forward_inplace(r.spectrum.data());
        r.time = r.spectrum;
        fft.inverse_inplace(r.time.data());
        return r;
    }

    TYPED_TEST(fft_rt_test, CopyProducesBitIdenticalOutput) {
        constexpr std::size_t n = 512;
        const auto            x = tap::dsp::test::random_signal<TypeParam>(n, 0x2545F491u);

        fft_t<TypeParam> original(n);
        // Source already warmed (tables built) so the copy carries built tables.
        const auto from_original = run_both(original, x);

        fft_t<TypeParam> copy(original);
        const auto       from_copy = run_both(copy, x);
        expect_bit_identical(from_copy.spectrum, from_original.spectrum, "copy forward");
        expect_bit_identical(from_copy.time, from_original.time, "copy inverse");

        // And the source is unaffected by having been copied.
        const auto again = run_both(original, x);
        expect_bit_identical(again.spectrum, from_original.spectrum, "source forward after copy");
        expect_bit_identical(again.time, from_original.time, "source inverse after copy");
    }

    TYPED_TEST(fft_rt_test, CopyOfUnwarmedSourceProducesBitIdenticalOutput) {
        constexpr std::size_t n = 512;
        const auto            x = tap::dsp::test::random_signal<TypeParam>(n, 0x2545F491u);

        // Copied BEFORE any transform: with lazy tables (F6) both objects
        // build their own; with constructor-built tables both carry the same.
        fft_t<TypeParam> original(n);
        fft_t<TypeParam> copy(original);
        const auto       from_copy     = run_both(copy, x);
        const auto       from_original = run_both(original, x);
        expect_bit_identical(from_copy.spectrum, from_original.spectrum, "unwarmed copy forward");
        expect_bit_identical(from_copy.time, from_original.time, "unwarmed copy inverse");
    }

    TYPED_TEST(fft_rt_test, CopyAssignmentProducesBitIdenticalOutput) {
        constexpr std::size_t n = 512;
        const auto            x = tap::dsp::test::random_signal<TypeParam>(n, 0x2545F491u);

        fft_t<TypeParam> original(n);
        const auto       from_original = run_both(original, x);

        // Assigned over an engine of a DIFFERENT geometry, so the assignment
        // has to replace every table, not just refresh one of the same size.
        fft_t<TypeParam> target(4096);
        target = original;
        EXPECT_EQ(target.size(), n);
        EXPECT_EQ(target.num_bins(), n / 2 + 1);
        const auto from_target = run_both(target, x);
        expect_bit_identical(from_target.spectrum, from_original.spectrum, "assigned forward");
        expect_bit_identical(from_target.time, from_original.time, "assigned inverse");
    }

} // namespace
