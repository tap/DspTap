# platform/ — bare-metal rig for the QEMU Cortex-M legs

Startup file and linker scripts for running the test battery bare-metal under
`qemu-system-arm` (see `cmake/arm-cortex-*.cmake` and `tests/CMakeLists.txt`).

| file | board / core | origin |
|---|---|---|
| `cortexm_startup.c` | all four legs | MuTap `platform/armv8m_startup.c` @ `142361b`, generalized to Armv7E-M (MSPLIM write guarded) and given the `TAP_DSP_TESTS_FAULT` exit path |
| `mps2_an386.ld` | `mps2-an386`, Cortex-M4 / M4F | new here; memory map from QEMU `hw/arm/mps2.c` (`FPGA_AN386`) |
| `mps2_an505.ld` | `mps2-an505`, Cortex-M33 | MuTap `platform/mps2_an505.ld` @ `142361b`, byte-identical apart from the header |
| `mps3_an547.ld` | `mps3-an547`, Cortex-M55 | MuTap `platform/mps3_an547.ld` @ `142361b`, byte-identical apart from the header |

All MIT (SPDX lines in each file). Provenance chain: SampleRateTap → RatioTap
(the MPS2 script) → MuTap → here, copied on 2026-09-17.

**Which copy is canonical.** This one, from the day a consumer pins a tree
containing it: DspTap is the shared substrate, the M4 support exists only
here, and the fault-exit path is an improvement MuTap does not have yet. The
duplication with MuTap's `platform/` is deliberate and temporary — consumers
retire their copies when they bump their submodule pin, and the family-wide
home is the taphouse consolidation the plan names
(`docs/audit-fft-and-code-smells.md`, Part 11 "Sharing"). Until then, fix
here first and port back, never the other way round.
