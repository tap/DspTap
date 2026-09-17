/* Linker script for the Arm MPS2 AN386 FPGA image (Cortex-M4) as modeled by
 * QEMU's mps2-an386 machine. Armv7E-M: VTOR resets to 0, so the vector table
 * goes at the base of SSRAM1. QEMU's -kernel loader places the ELF directly
 * into RAM (VMA == LMA, no load-time copy).
 *
 * Memory map (QEMU model, hw/arm/mps2.c, mps2_common_init(), FPGA_AN386):
 *   ZBT SSRAM1     4 MB @ 0x00000000  - vector table + code + rodata
 *                        (mirrored at 0x00400000)
 *   ZBT SSRAM2/3   4 MB @ 0x20000000  - data + bss + heap + stack
 *                        (mirrored at 0x20400000)
 *   PSRAM         16 MB @ 0x21000000  - unused here; the machine's default
 *                        RAM, available if a test ever needs more heap
 *   block RAM     16 KB @ 0x01000000  - unused
 */
MEMORY
{
    CODE (rx) : ORIGIN = 0x00000000, LENGTH = 4M
    DATA (rw) : ORIGIN = 0x20000000, LENGTH = 4M
}

__stack_top = ORIGIN(DATA) + LENGTH(DATA);

ENTRY(Reset_Handler)

SECTIONS
{
    .vectors : {
        KEEP(*(.vectors))
    } > CODE

    .text : {
        *(.text*)
        *(.rodata*)
        KEEP(*(.init))
        KEEP(*(.fini))
    } > CODE

    .ARM.extab : {
        *(.ARM.extab* .gnu.linkonce.armextab.*)
    } > CODE

    .ARM.exidx : {
        __exidx_start = .;
        *(.ARM.exidx* .gnu.linkonce.armexidx.*)
        __exidx_end = .;
    } > CODE

    .preinit_array : {
        PROVIDE_HIDDEN(__preinit_array_start = .);
        KEEP(*(.preinit_array*))
        PROVIDE_HIDDEN(__preinit_array_end = .);
    } > CODE

    .init_array : {
        PROVIDE_HIDDEN(__init_array_start = .);
        KEEP(*(SORT(.init_array.*)))
        KEEP(*(.init_array*))
        PROVIDE_HIDDEN(__init_array_end = .);
    } > CODE

    .fini_array : {
        PROVIDE_HIDDEN(__fini_array_start = .);
        KEEP(*(SORT(.fini_array.*)))
        KEEP(*(.fini_array*))
        PROVIDE_HIDDEN(__fini_array_end = .);
    } > CODE

    .data : {
        *(.data*)
    } > DATA

    .bss (NOLOAD) : {
        __bss_start__ = .;
        *(.bss*)
        *(COMMON)
        __bss_end__ = .;
    } > DATA

    /* Stack lives at the top of DATA; cap the heap 64 KB below it. */
    .heap (NOLOAD) : ALIGN(8) {
        __heap_start__ = .;
        . = ORIGIN(DATA) + LENGTH(DATA) - 64K;
        __heap_end__ = .;
    } > DATA

    /* Same symbol as the Armv8-M scripts so the shared startup links; the
     * M4 has no MSPLIM, so this is documentation, not enforcement. */
    __stack_limit = __heap_end__;

    /* librdimon's (unused, weak) _sbrk references `end`; satisfy it. */
    PROVIDE(end = __heap_start__);
}
