# Vendored CMSIS-DSP subset (float32 real FFT, Helium backend)

This is a **minimal subset** of Arm's CMSIS-DSP and CMSIS-Core, vendored so
the optional Cortex-M55 (Helium/MVE) FFT backend can be built without a
submodule or network fetch. It is compiled **only** when the `MUTAP_FFT_CMSIS`
CMake option is ON (ARM cross builds); nothing here is touched by the default
desktop/Hexagon builds, which stay on the Ooura FFT.

## Provenance

| Component  | Upstream                                   | Pinned commit  | Date       |
|------------|--------------------------------------------|----------------|------------|
| CMSIS-DSP  | github.com/ARM-software/CMSIS-DSP          | `918014f0ba96` | 2026-07-17 |
| CMSIS-Core | github.com/ARM-software/CMSIS_6 (Core)     | `7f62ddc8ab8e` | 2026-06-30 |

License: Apache-2.0 (see `LICENSE`; SPDX headers retained in every file).

## Contents (why each file is here)

The closure was computed with `gcc -M` over the eight compiled sources for
`-mcpu=cortex-m55`; only files that closure actually opens are vendored.

- `Source/TransformFunctions/` — `arm_rfft_fast_f32`, `arm_rfft_fast_init_f32`,
  `arm_cfft_f32`, `arm_cfft_init_f32`, `arm_bitreversal2` (the real-FFT entry
  point and its complex-FFT + bit-reversal dependencies).
- `Source/CommonTables/` — `arm_common_tables`, `arm_const_structs`,
  `arm_mve_tables` (twiddle factors, radix structs, MVE-specific tables).
- `Include/` — `arm_math.h` umbrella and the headers it pulls (types, memory,
  helium utils, the table declarations, and the `dsp/*.h` group headers the
  umbrella includes even though only the transform group is used).
- `PrivateInclude/` — `arm_compiler_specific.h`, `arm_vec_fft.h`.
- `cmsis-core/Include/` — `cmsis_compiler.h`, `cmsis_gcc.h`,
  `m-profile/cmsis_gcc_m.h` (the only CMSIS-Core headers the closure needs;
  the GCC-only compiler shims, no device/core header required).

## Refreshing

To bump versions: re-clone both repos at the new commits, re-run the `gcc -M`
closure (see `bench/README.md` / the FFT-backend section of
`docs/optimization.md`), copy the resulting file list here, update the table
above, and re-run the parity harness (`tests/test_fft_backend.cpp`) plus the
full float32 battery on the M55 leg. Do **not** hand-edit vendored sources.
