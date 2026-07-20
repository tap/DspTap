# Third-party notices

DspTap's own wrapper code (`include/dsptap/`, `tests/`, the build files) is
licensed under the MIT License — see [`LICENSE`](LICENSE).

It bundles the following third-party code under `third_party/`, each retaining
its own license and SPDX headers:

## Ooura General Purpose FFT Package

- Path: `third_party/ooura/` (`fftsg.c`, and `fftsg_float.c`, a
  single-precision instantiation of the same source).
- Author: Takuya Ooura.
- License: permissive/public — the package's terms permit free use provided
  the copyright notice is retained. See `third_party/ooura/readme.txt`.

## CMSIS-DSP and CMSIS-Core (subset)

- Path: `third_party/cmsis-dsp/`.
- Origin: Arm Limited — a minimal subset (float32 real FFT + Helium/MVE
  transform closure) vendored for the optional Cortex-M55 backend.
- License: Apache-2.0 — see `third_party/cmsis-dsp/LICENSE`; SPDX headers
  retained in every vendored file. Provenance and refresh procedure in
  `third_party/cmsis-dsp/VENDOR.md`.
- Compiled only when `DSPTAP_FFT_CMSIS` is ON (Arm cross builds); untouched by
  the default desktop / Apple / Hexagon builds.
