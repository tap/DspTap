# Third-party notices

DspTap's own code — everything outside `third_party/`: `include/tap/dsp/`,
`tests/`, `tools/`, `notebooks/`, `bench/`, `scripts/`, `cmake/`, `docs/`
and the build files — is licensed under the MIT License — see
[`LICENSE`](LICENSE).

It bundles the following third-party code under `third_party/`, each retaining
its own license text:

## Ooura General Purpose FFT Package

- Path: `third_party/ooura/` (`fftsg.c`, the split-radix "Fast Version III",
  and `fftsg_float.c`, a 62-line DspTap wrapper that `#include`s it under a
  type remap).
- Author: Takuya Ooura. Copyright(C) 1996-2001 Takuya OOURA.
- License: the package's own terms, stated in `readme.txt`, the only upstream
  license text (the banner at the top of the vendored `fftsg.c` is Tap's copy
  of it; the upstream file carries no notice of its own). There is no other
  license text; in particular the package is neither public domain nor under
  a named open-source license. The full grant reads:

  > You may use, copy, modify this code for any purpose and without fee. You
  > may distribute this ORIGINAL package.

  What that grants: use, copying and modification of the code, for any purpose
  and without fee; and distribution of the *original* package. Distribution of
  a modified derivative is not expressly granted.

- What DspTap ships today: one source file of the package, `fftsg.c`, plus its
  `readme.txt`. The vendored `fftsg.c` is textually identical to the
  2006-12-28 `fft.tgz` except for a provenance banner Tap added at the top
  (the upstream file carries no notice of its own; the notice lives in
  `readme.txt`) and stripped trailing whitespace. That is a partial copy of
  the original package with the notice attached, not a modified transform;
  the maintainer's reading is that it is within the intent of the
  distribution grant, and the draft email in `docs/fft-design.md` puts the
  question to the author.
- What DspTap relies on going forward: the planned C++20 port of `rdft`
  (`include/tap/dsp/fft/split_radix.h`, a statement-for-statement
  transliteration) is a **derivative work, not the ORIGINAL package**, and its
  redistribution relies on the **modification grant** ("modify this code for
  any purpose"). Precedent exists but is not relied on: WebRTC/Chromium ship a
  modified `fft4g.c` under `common_audio/third_party/ooura/`, and their
  `LICENSE` quotes a broader notice — "You may use, copy, modify and
  distribute this code for any purpose (include commercial use) and without
  fee. Please refer to this package when you modify this code." — that is not
  in the `fft.tgz` readme. DspTap has not traced the origin of that text and
  does not rely on it. The maintainer decides whether to ask the author for an
  explicit statement on derivative distribution (draft in
  `docs/fft-design.md`); the decision and any outcome are recorded here at
  Stage 2c.
- SPDX plan: the port header carries `SPDX-License-Identifier:
  LicenseRef-Ooura AND MIT` — Ooura's notice verbatim as the governing terms
  for the derived portion, MIT for the wrapper and DspTap's additions — plus a
  line stating that it is a derivative work with the modification copyright.
  `LicenseRef-Ooura` denotes the `Copyright:` block of
  `third_party/ooura/readme.txt` (lines 140–145), reproduced verbatim as
  [`LICENSES/LicenseRef-Ooura.txt`](LICENSES/LicenseRef-Ooura.txt) so the
  reference resolves. `fft.h` (the wrapper) and everything else DspTap wrote
  stay plain MIT.
- `readme.txt` stays at `third_party/ooura/readme.txt` permanently, as the
  license record for the derived code; only `fftsg.c` moves to the test
  reference tree at Stage 2c and is retired afterwards.
- This is a maintainer judgement call, not legal advice. This file is the
  canonical statement; `README.md` and `docs/fft-design.md` point here.

## CMSIS-DSP and CMSIS-Core (subset)

- Path: `third_party/cmsis-dsp/`.
- Origin: Arm Limited — a minimal subset (float32 real FFT + Helium/MVE
  transform closure) vendored for the optional Cortex-M55 backend.
- License: Apache-2.0 — see `third_party/cmsis-dsp/LICENSE`; SPDX headers
  retained in every vendored file. Provenance and refresh procedure in
  `third_party/cmsis-dsp/VENDOR.md`.
- Compiled only when `TAP_DSP_FFT_CMSIS` is ON (Arm cross builds); untouched by
  the default desktop / Apple / Hexagon builds.
