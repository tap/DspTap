// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// Pins tap/dsp/analysis/heap_meter.h: the requested bytes of every form of the
// global operator new are counted while a meter is open, live and at peak,
// freed exactly, over-aligned forms included, and an allocation made before
// the meter opened is never counted either way. This executable carries the
// replacement operators (tests/CMakeLists.txt gives it a target of its own),
// so it is also the proof that a gtest binary runs whole on them.

#define TAP_DSP_HEAP_METER_REPLACE_GLOBAL_NEW
#include <cstdint>
#include <memory>
#include <new>
#include <vector>

#include <gtest/gtest.h>

#include "tap/dsp/analysis/heap_meter.h"

namespace {

    using tap::dsp::analysis::heap_meter;

    // A new-expression whose result goes nowhere may be elided ([expr.new];
    // clang does, under optimization), so every pointer a test allocates
    // escapes through this.
    void* volatile g_escape = nullptr;
    template <typename T>
    T* escape(T* p) {
        g_escape = p;
        return p;
    }

    TEST(HeapMeter, IsInstalledInThisExecutable) {
        EXPECT_TRUE(heap_meter::installed());
    }

    TEST(HeapMeter, CountsTheRequestedBytesLiveAndAtPeak) {
        const heap_meter meter;
        auto             a = std::make_unique<std::vector<std::int32_t>>(100); // the vector object, then its 400 bytes
        {
            const std::vector<double> b(50); // 400 bytes, freed at the brace
            EXPECT_EQ(meter.stats().live_bytes, sizeof(std::vector<std::int32_t>) + 400 + 400);
        }
        const auto s = meter.stats();
        EXPECT_EQ(s.live_bytes, sizeof(std::vector<std::int32_t>) + 400);
        EXPECT_EQ(s.peak_bytes, sizeof(std::vector<std::int32_t>) + 800);
        EXPECT_EQ(s.total_bytes, sizeof(std::vector<std::int32_t>) + 800);
        EXPECT_EQ(s.allocations, 3U);
        a.reset();
        EXPECT_EQ(meter.stats().live_bytes, 0U);
    }

    TEST(HeapMeter, OverAlignedAllocationsAreAlignedCountedAndFreed) {
        struct alignas(256) wide { // exactly two alignments: no padding (MSVC C4324)
            unsigned char bytes[512];
        };
        const heap_meter meter;
        auto             p = std::unique_ptr<wide>(escape(new wide));
        EXPECT_EQ(reinterpret_cast<std::uintptr_t>(p.get()) % 256U, 0U);
        EXPECT_EQ(meter.stats().live_bytes, sizeof(wide));
        auto* q = escape(new wide[3]);
        EXPECT_EQ(reinterpret_cast<std::uintptr_t>(q) % 256U, 0U);
        EXPECT_GE(meter.stats().live_bytes, sizeof(wide) * 4); // an array may carry a cookie
        delete[] q;
        p.reset();
        EXPECT_EQ(meter.stats().live_bytes, 0U);
        EXPECT_EQ(meter.stats().allocations, 2U);
    }

    TEST(HeapMeter, NothrowAndArrayFormsAreCounted) {
        const heap_meter meter;
        auto*            a = escape(new (std::nothrow) std::int64_t(7));
        auto*            b = escape(new std::int16_t[10]);
        EXPECT_EQ(meter.stats().live_bytes, sizeof(std::int64_t) + 10 * sizeof(std::int16_t));
        delete a;
        delete[] b;
        EXPECT_EQ(meter.stats().live_bytes, 0U);
    }

    TEST(HeapMeter, AnAllocationFromBeforeTheMeterIsNeverCounted) {
        auto             early = std::make_unique<std::vector<float>>(1000);
        const heap_meter meter;
        early.reset(); // its bytes were never counted, so they are not subtracted
        const auto s = meter.stats();
        EXPECT_EQ(s.live_bytes, 0U);
        EXPECT_EQ(s.total_bytes, 0U);
        EXPECT_EQ(s.allocations, 0U);
    }

    TEST(HeapMeter, AZeroByteRequestIsOneAllocationOfNoBytes) {
        const heap_meter meter;
        void*            p = ::operator new(0);
        EXPECT_NE(p, nullptr);
        EXPECT_EQ(meter.stats().allocations, 1U);
        EXPECT_EQ(meter.stats().live_bytes, 0U);
        ::operator delete(p);
    }

    TEST(HeapMeter, EachMeterStartsFromZero) {
        {
            const heap_meter        first;
            const std::vector<char> v(64);
            EXPECT_EQ(first.stats().live_bytes, 64U);
        }
        const heap_meter second;
        EXPECT_EQ(second.stats().peak_bytes, 0U);
        EXPECT_EQ(second.stats().allocations, 0U);
    }

} // namespace
