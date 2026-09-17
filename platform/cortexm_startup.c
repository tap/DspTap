/* Minimal bare-metal startup for Cortex-M targets under QEMU (Cortex-M4 on
 * mps2-an386, Cortex-M33 on mps2-an505, Cortex-M55 on mps3-an547),
 * replacing the toolchain crt0 (linked with -nostartfiles):
 *   - vector table with initial stack pointer and Reset_Handler
 *   - MSPLIM stack-limit guard on Armv8-M Mainline (M33/M55); the register
 *     does not exist on Armv7E-M (M4), so the write is compiled out there
 *   - FPU/MVE coprocessor enable before any FP instruction executes
 *   - .bss zeroing (QEMU's ELF loader already placed .data)
 *   - semihosting stdio (librdimon), C++ static constructors, main, exit
 *   - fault handlers that print TAP_DSP_TESTS_FAULT over semihosting and
 *     exit, so a fault fails the CTest run in seconds with a diagnostic
 *     instead of parking until the timeout
 *   - deterministic _sbrk over the linker-defined heap region (overrides
 *     librdimon's weak version, whose limit depends on the semihosting
 *     SYS_HEAPINFO call returning sensible values for this board)
 *   - 64-bit __atomic_* helpers: M-profile has no 64-bit exclusives and the
 *     bare-metal toolchain ships no libatomic; PRIMASK critical sections
 *     are sufficient on a single-core part
 *
 * The toolchain file passes this to the link line with `-x c`: under C the
 * vector table's address-constant initializers are guaranteed link-time
 * constants (a C++ compile could legally lower them to dynamic
 * initialization, leaving the table zeroed at reset). The extern "C"
 * guards keep the file safe if it is ever compiled as C++ anyway.
 *
 * Only the __atomic_* helpers the library and runtime currently need are
 * provided; any future use of others (e.g. compare-exchange) fails loudly
 * at link time.
 */
// SPDX-License-Identifier: MIT
// Copyright 2025-2026 Timothy Place and the DspTap contributors.
// Copied from MuTap's platform/armv8m_startup.c at 142361b (2026-09-17), itself
// ported from SampleRateTap's (same license); generalized here to Armv7E-M and
// given the fault-exit path. See platform/README.md for which copy is canonical.
#include <errno.h>
#include <stddef.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>

#ifdef __cplusplus
extern "C" {
#endif

extern uint32_t __bss_start__, __bss_end__;
extern uint32_t __stack_top;
extern uint32_t __stack_limit;
extern char     __heap_start__, __heap_end__;

extern void __libc_init_array(void);
extern void initialise_monitor_handles(void);
extern int  main(int argc, char** argv);
extern void exit(int) __attribute__((noreturn));

void* __dso_handle;

/* Referenced by newlib's fini machinery; nothing to do without crti/crtn. */
void _init(void) {}
void _fini(void) {}

void* _sbrk(ptrdiff_t increment) {
    static char* brk = &__heap_start__;
    if (brk + increment > &__heap_end__) {
        errno = ENOMEM;
        return (void*)-1;
    }
    char* prev = brk;
    brk += increment;
    return prev;
}

static inline uint32_t irq_lock(void) {
    uint32_t primask;
    __asm volatile("mrs %0, PRIMASK\n cpsid i" : "=r"(primask)::"memory");
    return primask;
}

static inline void irq_restore(uint32_t primask) {
    __asm volatile("msr PRIMASK, %0" ::"r"(primask) : "memory");
}

uint64_t __atomic_load_8(const volatile void* ptr, int memorder) {
    (void)memorder;
    const uint32_t m = irq_lock();
    const uint64_t v = *(const volatile uint64_t*)ptr;
    irq_restore(m);
    return v;
}

void __atomic_store_8(volatile void* ptr, uint64_t value, int memorder) {
    (void)memorder;
    const uint32_t m         = irq_lock();
    *(volatile uint64_t*)ptr = value;
    irq_restore(m);
}

uint64_t __atomic_fetch_add_8(volatile void* ptr, uint64_t value, int memorder) {
    (void)memorder;
    const uint32_t m         = irq_lock();
    const uint64_t prev      = *(volatile uint64_t*)ptr;
    *(volatile uint64_t*)ptr = prev + value;
    irq_restore(m);
    return prev;
}

uint64_t __atomic_exchange_8(volatile void* ptr, uint64_t value, int memorder) {
    (void)memorder;
    const uint32_t m         = irq_lock();
    const uint64_t prev      = *(volatile uint64_t*)ptr;
    *(volatile uint64_t*)ptr = value;
    irq_restore(m);
    return prev;
}

void Reset_Handler(void) {
#if defined(__ARM_ARCH_8M_MAIN__) || defined(__ARM_ARCH_8_1M_MAIN__)
    /* MSPLIM exists on Armv8-M Mainline only (M33/M55): a main-stack
     * overflow past __stack_limit raises a fault instead of silently
     * corrupting whatever sits below the stack. GCC defines
     * __ARM_ARCH_8M_MAIN__ for both cores; clang spells the M55's Armv8.1-M
     * as __ARM_ARCH_8_1M_MAIN__. Armv7E-M (M4) has no such register — its
     * linker script still defines __stack_limit so the symbol layout is
     * identical, but nothing enforces it. */
    __asm volatile("msr msplim, %0" ::"r"(&__stack_limit));
#else
    (void)&__stack_limit;
#endif

    /* Grant full access to CP10/CP11 (scalar FPU + MVE) first: code below
     * may legitimately use FP registers once newlib is involved. On a core
     * without an FPU (or a soft-float build) the write is harmless. */
    volatile uint32_t* const cpacr = (volatile uint32_t*)0xE000ED88u;
    *cpacr |= 0xFu << 20;
    __asm volatile("dsb\n isb" ::: "memory");

    memset(&__bss_start__, 0, (size_t)((char*)&__bss_end__ - (char*)&__bss_start__));

    initialise_monitor_handles(); /* semihosting stdin/stdout/stderr */
    __libc_init_array();          /* C++ static constructors */
    exit(main(0, (char**)0));
}

/* Fault path. Print a marker CTest treats as failure
 * (FAIL_REGULAR_EXPRESSION in tests/CMakeLists.txt) together with the active
 * exception number from IPSR, then leave through semihosting SYS_EXIT
 * (_exit), so a fault costs seconds rather than the test TIMEOUT and leaves a
 * diagnostic in the log. write(2) is unbuffered and needs no stdio state; the
 * buffer is static so the handler itself touches almost no stack (an MSPLIM
 * overflow has none to give). Define TAP_DSP_FAULT_BKPT to stop in a
 * debugger first. */
static void fault_exit(const char* what) {
    static char buf[64];
    uint32_t    ipsr;
    __asm volatile("mrs %0, ipsr" : "=r"(ipsr));
    size_t n = 0;
    for (const char* p = "TAP_DSP_TESTS_FAULT "; *p != '\0'; ++p) {
        buf[n++] = *p;
    }
    for (const char* p = what; *p != '\0'; ++p) {
        buf[n++] = *p;
    }
    for (const char* p = " ipsr="; *p != '\0'; ++p) {
        buf[n++] = *p;
    }
    char     digits[10];
    size_t   nd = 0;
    uint32_t v  = ipsr & 0x1FFu;
    do {
        digits[nd++] = (char)('0' + (v % 10u));
        v /= 10u;
    } while (v != 0u);
    while (nd > 0) {
        buf[n++] = digits[--nd];
    }
    buf[n++] = '\n';
#ifdef TAP_DSP_FAULT_BKPT
    __asm volatile("bkpt #0");
#endif
    (void)write(2, buf, n);
    _exit(2);
}

void Default_Handler(void) {
    fault_exit("unexpected exception");
}

void HardFault_Handler(void) {
    fault_exit("HardFault");
}

__attribute__((section(".vectors"), used)) static const uintptr_t vectors[16] = {
    (uintptr_t)&__stack_top,
    (uintptr_t)&Reset_Handler,
    (uintptr_t)&Default_Handler,   /* NMI */
    (uintptr_t)&HardFault_Handler, /* HardFault */
    (uintptr_t)&Default_Handler,   /* MemManage */
    (uintptr_t)&Default_Handler,   /* BusFault */
    (uintptr_t)&Default_Handler,   /* UsageFault */
    (uintptr_t)&Default_Handler,   /* SecureFault */
    0,
    0,
    0,
    (uintptr_t)&Default_Handler, /* SVCall */
    (uintptr_t)&Default_Handler, /* DebugMonitor */
    0,
    (uintptr_t)&Default_Handler, /* PendSV */
    (uintptr_t)&Default_Handler, /* SysTick */
};

#ifdef __cplusplus
} /* extern "C" */
#endif
