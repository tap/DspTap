/// @file heap_meter.h
/// @brief Heap accounting for footprint tests: the bytes a scope allocates through
/// the global operator new, live and at peak.
// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// A footprint claim ("one instance holds this many bytes") is only as good as
// the instrument that checks it. The family's converters allocate through
// std::vector with the default allocator, so the one place every byte passes is
// the global operator new; this header measures there.
//
// Two parts. tap::dsp::analysis::heap_meter is a scope: while one is open,
// every allocation through the global operator new is counted, and stats()
// reports the bytes still live, the peak of the live bytes, the total and the
// number of allocations. The replacements of the global operator new and
// operator delete that do the counting are emitted only in the translation
// unit that defines TAP_DSP_HEAP_METER_REPLACE_GLOBAL_NEW before including this
// header: exactly one per executable (a second definition is a link error, as
// for any replacement), and in an executable of its own, because the
// replacement serves every allocation the program makes. heap_meter::installed()
// says whether the executable carries it; without it a meter counts nothing.
//
// The numbers are the bytes the caller asked for, not the allocator's: each
// allocation carries a header in front of it (its size, its base pointer and
// whether it was made inside a meter), so the replacement can account sized and
// unsized deletes alike, over-aligned forms included. An allocation made before
// the meter opened and freed inside it is not counted either way. One meter at
// a time; the counters are atomic, so allocations from other threads while a
// meter is open are counted too.
//
// What a container asks operator new for is not always the bytes it stores:
// MSVC's standard library on x86 and x64 pads every std::allocator request of
// 4096 bytes or more (for a type of default alignment) to align the block to 32
// bytes, sizeof(void*) + 31 bytes in front of it (one pointer more under
// _DEBUG). heap_meter::allocator_request(n) is that request for n stored bytes,
// n on every other library, so a footprint formula written per allocation
// through it holds on every host.
//
// Not under AddressSanitizer, ThreadSanitizer or MemorySanitizer: their runtimes
// define the global operator new themselves, so the replacement cannot link
// beside them (an executable that defines the macro is left out of those builds).

#pragma once

#include <atomic>
#include <cstddef>

namespace tap::dsp::analysis {

    /// What a heap_meter has seen since it opened.
    struct heap_stats {
        std::size_t live_bytes  = 0; ///< requested bytes allocated while open and not yet freed
        std::size_t peak_bytes  = 0; ///< the largest live_bytes reached while open
        std::size_t total_bytes = 0; ///< requested bytes of every allocation while open
        std::size_t allocations = 0; ///< number of allocations while open
    };

    namespace detail::heap_meter_state {
        inline std::atomic<bool>        installed{false};
        inline std::atomic<bool>        open{false};
        inline std::atomic<std::size_t> live{0};
        inline std::atomic<std::size_t> peak{0};
        inline std::atomic<std::size_t> total{0};
        inline std::atomic<std::size_t> count{0};

        inline void on_allocate(std::size_t bytes) noexcept {
            const std::size_t now = live.fetch_add(bytes, std::memory_order_relaxed) + bytes;
            total.fetch_add(bytes, std::memory_order_relaxed);
            count.fetch_add(1, std::memory_order_relaxed);
            std::size_t seen = peak.load(std::memory_order_relaxed);
            while (now > seen && !peak.compare_exchange_weak(seen, now, std::memory_order_relaxed)) {
            }
        }

        inline void on_free(std::size_t bytes) noexcept {
            live.fetch_sub(bytes, std::memory_order_relaxed);
        }
    } // namespace detail::heap_meter_state

    /// Counts the global heap while it is open (see the file comment). One at a time.
    class heap_meter {
      public:
        heap_meter() noexcept {
            namespace s = detail::heap_meter_state;
            s::live.store(0, std::memory_order_relaxed);
            s::peak.store(0, std::memory_order_relaxed);
            s::total.store(0, std::memory_order_relaxed);
            s::count.store(0, std::memory_order_relaxed);
            s::open.store(true, std::memory_order_release);
        }
        ~heap_meter() { detail::heap_meter_state::open.store(false, std::memory_order_release); }
        heap_meter(const heap_meter&)            = delete;
        heap_meter& operator=(const heap_meter&) = delete;
        heap_meter(heap_meter&&)                 = delete;
        heap_meter& operator=(heap_meter&&)      = delete;

        /// The counts so far.
        heap_stats stats() const noexcept {
            namespace s = detail::heap_meter_state;
            return {s::live.load(std::memory_order_relaxed), s::peak.load(std::memory_order_relaxed),
                    s::total.load(std::memory_order_relaxed), s::count.load(std::memory_order_relaxed)};
        }

        /// True when this executable carries the replacement operators.
        static bool installed() noexcept { return detail::heap_meter_state::installed.load(); }

        /// The bytes std::allocator asks operator new for to store n bytes of a
        /// type of default alignment (see the file comment): n plus
        /// k_big_allocation_padding from k_big_allocation_bytes up.
        static constexpr std::size_t allocator_request(std::size_t n) noexcept {
            return n >= k_big_allocation_bytes ? n + k_big_allocation_padding : n;
        }

        /// The request size from which MSVC's std::allocator pads, and by how
        /// much (0 on every other standard library).
        static constexpr std::size_t k_big_allocation_bytes = 4096;
#if defined(_MSVC_STL_VERSION) && (defined(_M_IX86) || (defined(_M_X64) && !defined(_M_ARM64EC)))
#if defined(_DEBUG)
        static constexpr std::size_t k_big_allocation_padding = 2 * sizeof(void*) + 31;
#else
        static constexpr std::size_t k_big_allocation_padding = sizeof(void*) + 31;
#endif
#else
        static constexpr std::size_t k_big_allocation_padding = 0;
#endif
    };

} // namespace tap::dsp::analysis

#if defined(TAP_DSP_HEAP_METER_REPLACE_GLOBAL_NEW)

#include <cstdint>
#include <cstdlib>
#include <new>

#include "tap/dsp/detail/expects.h"

namespace tap::dsp::analysis::detail::heap_meter_state {

    // In front of every block: what delete needs. The space it takes is
    // rounded up to the default new alignment, so the block behind it keeps
    // that alignment (a rounded size, not alignas: MSVC's /W4 warns when an
    // alignment specifier pads a structure).
    struct block_header {
        std::size_t bytes;
        void*       base;
        bool        counted;
    };
    inline constexpr std::size_t k_header_bytes = (sizeof(block_header) + __STDCPP_DEFAULT_NEW_ALIGNMENT__ - 1)
                                                  / __STDCPP_DEFAULT_NEW_ALIGNMENT__ * __STDCPP_DEFAULT_NEW_ALIGNMENT__;

    inline void* allocate(std::size_t bytes, std::size_t align) noexcept {
        if (align < __STDCPP_DEFAULT_NEW_ALIGNMENT__) {
            align = __STDCPP_DEFAULT_NEW_ALIGNMENT__;
        }
        // Room for the header and for aligning the block behind it.
        void* const base = std::malloc(k_header_bytes + align + bytes);
        if (base == nullptr) {
            return nullptr;
        }
        const auto first = reinterpret_cast<std::uintptr_t>(base) + k_header_bytes;
        const auto user  = (first + align - 1) & ~(static_cast<std::uintptr_t>(align) - 1);
        auto*      head  = reinterpret_cast<block_header*>(user - k_header_bytes);
        head->bytes      = bytes;
        head->base       = base;
        head->counted    = open.load(std::memory_order_acquire);
        if (head->counted) {
            on_allocate(bytes);
        }
        return reinterpret_cast<void*>(user);
    }

    inline void release(void* p) noexcept {
        if (p == nullptr) {
            return;
        }
        const auto* head = reinterpret_cast<const block_header*>(static_cast<const char*>(p) - k_header_bytes);
        if (head->counted) {
            on_free(head->bytes);
        }
        std::free(head->base);
    }

    inline void* allocate_or_raise(std::size_t bytes, std::size_t align) {
        void* p = allocate(bytes, align);
        if (p == nullptr) {
            tap::dsp::raise<std::bad_alloc>();
        }
        return p;
    }

    // Marks the executable as carrying the replacement before main runs.
    inline const bool k_installed_at_startup = [] {
        installed.store(true);
        return true;
    }();

} // namespace tap::dsp::analysis::detail::heap_meter_state

// The replaceable global allocation functions ([new.delete]), every form.
// NOLINTBEGIN(readability-identifier-naming)
void* operator new(std::size_t n) {
    return tap::dsp::analysis::detail::heap_meter_state::allocate_or_raise(n, 0);
}
void* operator new[](std::size_t n) {
    return tap::dsp::analysis::detail::heap_meter_state::allocate_or_raise(n, 0);
}
void* operator new(std::size_t n, std::align_val_t a) {
    return tap::dsp::analysis::detail::heap_meter_state::allocate_or_raise(n, static_cast<std::size_t>(a));
}
void* operator new[](std::size_t n, std::align_val_t a) {
    return tap::dsp::analysis::detail::heap_meter_state::allocate_or_raise(n, static_cast<std::size_t>(a));
}
void* operator new(std::size_t n, const std::nothrow_t&) noexcept {
    return tap::dsp::analysis::detail::heap_meter_state::allocate(n, 0);
}
void* operator new[](std::size_t n, const std::nothrow_t&) noexcept {
    return tap::dsp::analysis::detail::heap_meter_state::allocate(n, 0);
}
void* operator new(std::size_t n, std::align_val_t a, const std::nothrow_t&) noexcept {
    return tap::dsp::analysis::detail::heap_meter_state::allocate(n, static_cast<std::size_t>(a));
}
void* operator new[](std::size_t n, std::align_val_t a, const std::nothrow_t&) noexcept {
    return tap::dsp::analysis::detail::heap_meter_state::allocate(n, static_cast<std::size_t>(a));
}
void operator delete(void* p) noexcept {
    tap::dsp::analysis::detail::heap_meter_state::release(p);
}
void operator delete[](void* p) noexcept {
    tap::dsp::analysis::detail::heap_meter_state::release(p);
}
void operator delete(void* p, std::size_t) noexcept {
    tap::dsp::analysis::detail::heap_meter_state::release(p);
}
void operator delete[](void* p, std::size_t) noexcept {
    tap::dsp::analysis::detail::heap_meter_state::release(p);
}
void operator delete(void* p, std::align_val_t) noexcept {
    tap::dsp::analysis::detail::heap_meter_state::release(p);
}
void operator delete[](void* p, std::align_val_t) noexcept {
    tap::dsp::analysis::detail::heap_meter_state::release(p);
}
void operator delete(void* p, std::size_t, std::align_val_t) noexcept {
    tap::dsp::analysis::detail::heap_meter_state::release(p);
}
void operator delete[](void* p, std::size_t, std::align_val_t) noexcept {
    tap::dsp::analysis::detail::heap_meter_state::release(p);
}
void operator delete(void* p, const std::nothrow_t&) noexcept {
    tap::dsp::analysis::detail::heap_meter_state::release(p);
}
void operator delete[](void* p, const std::nothrow_t&) noexcept {
    tap::dsp::analysis::detail::heap_meter_state::release(p);
}
void operator delete(void* p, std::align_val_t, const std::nothrow_t&) noexcept {
    tap::dsp::analysis::detail::heap_meter_state::release(p);
}
void operator delete[](void* p, std::align_val_t, const std::nothrow_t&) noexcept {
    tap::dsp::analysis::detail::heap_meter_state::release(p);
}
// NOLINTEND(readability-identifier-naming)

#endif // TAP_DSP_HEAP_METER_REPLACE_GLOBAL_NEW
