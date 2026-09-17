# Third-party notices

DspTap's own code (`include/tap/dsp/`, `tests/`, `tools/`, `notebooks/`, the
build files) is licensed under the MIT License — see [`LICENSE`](LICENSE).

It bundles the following third-party code under `third_party/`, each retaining
its own license text and headers:

## Ooura General Purpose FFT Package

- Path: `third_party/ooura/` (`fftsg.c`, the split-radix "Fast Version III",
  and `fftsg_float.c`, a single-precision instantiation of the same source).
- Author: Takuya Ooura. Copyright (C) 1996-2001 Takuya OOURA.
- License: the package's own terms, stated in full in `readme.txt` and in the
  `fftsg.c` banner. There is no other license text; in particular the package
  is neither public domain nor under a named open-source license. The full
  grant reads:

  > You may use, copy, modify this code for any purpose and without fee. You
  > may distribute this ORIGINAL package.

  What that grants: use, copying and modification of the code, for any purpose
  and without fee; and distribution of the *original* package. Distribution of
  a modified derivative is not expressly granted. Today DspTap distributes the
  original `fftsg.c` unmodified (`fftsg_float.c` is a 62-line DspTap wrapper
  that `#include`s it under a type remap), which the distribution grant covers.

- What DspTap relies on going forward: the planned C++20 port of `rdft`
  (`include/tap/dsp/fft/split_radix.h`, a statement-for-statement
  transliteration) is a **derivative work, not the ORIGINAL package**, and its
  redistribution relies on the **modification grant** ("modify this code for
  any purpose") together with the widespread practice of shipping modified
  Ooura code with the notice retained — practice plus the author's informally
  reported lack of objection, not license text. Before the vendored C leaves
  the main tree, the maintainer will attempt to reach the address in Ooura's
  notice for an explicit statement on derivative distribution and record the
  outcome here either way (draft in `docs/fft-design.md`).
- SPDX plan: the port header carries `SPDX-License-Identifier:
  LicenseRef-Ooura AND MIT` — Ooura's notice verbatim as the governing terms
  for the derived portion, MIT for the wrapper and DspTap's additions — plus a
  line stating that it is a derivative work with the modification copyright.
  `fft.h` (the wrapper) and everything else DspTap wrote stay plain MIT.
- `third_party/ooura/readme.txt` stays in-tree permanently, as the license
  record for the derived code, after `fftsg.c` itself has moved to the test
  reference tree and been retired.

## CMSIS-DSP and CMSIS-Core (subset)

- Path: `third_party/cmsis-dsp/`.
- Origin: Arm Limited — a minimal subset (float32 real FFT + Helium/MVE
  transform closure) vendored for the optional Cortex-M55 backend.
- License: Apache-2.0 — see `third_party/cmsis-dsp/LICENSE`; SPDX headers
  retained in every vendored file. Provenance and refresh procedure in
  `third_party/cmsis-dsp/VENDOR.md` (which still names the option by its
  pre-rename MuTap spelling; the option below is the one this repo defines).
- Compiled only when `TAP_DSP_FFT_CMSIS` is ON (Arm cross builds); untouched by
  the default desktop / Apple / Hexagon builds.
