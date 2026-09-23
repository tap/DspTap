# Third-party notices

DspTap's own code — everything outside `third_party/` and the two Ooura
files under `tests/reference/ooura/` (`fftsg.c`, `fftsg_float.c`; the
directory-local `.clang-format` beside them and the declaration header
`tests/reference/ooura_rdft.h` above them are DspTap's): `include/tap/dsp/`,
`tests/`, `tools/`, `notebooks/`, `bench/`, `scripts/`, `cmake/`, `platform/`
(the bare-metal startup file and linker scripts the QEMU legs link, carried
from MuTap's MIT copies; `platform/README.md` records their origin), `docs/`
and the build files — is licensed under the MIT License — see
[`LICENSE`](LICENSE).

It bundles the following third-party code, each retaining its own license
text:

## Ooura General Purpose FFT Package

- Paths: `third_party/ooura/readme.txt` (the package's readme: the only
  upstream license text, kept at this path permanently) and, **outside the
  shipping tree**, `tests/reference/ooura/` (`fftsg.c`, the split-radix "Fast
  Version III", and `fftsg_float.c`, a 65-line DspTap wrapper that
  `#include`s it under a type remap): the reference copy that the parity gate
  `tests/test_fft_parity_ooura.cpp` compiles and compares the shipping C++20
  port against. Nothing that ships compiles these two files; since Stage 2c
  of `docs/audit-fft-and-code-smells.md` (Decision D6) the library is
  header-only and no vendored C is part of it (the one compiled library,
  `tap_dsp_fft`, exists only under `TAP_DSP_FFT_CMSIS` and carries the
  CMSIS-DSP objects below). D6 retires the reference copy after both MuTap
  and MuTap-Max pin a tree containing 2c; `readme.txt` stays.
- Author: Takuya Ooura. Copyright(C) 1996-2001 Takuya OOURA.
- License: the package's own terms, stated in `readme.txt`, the only upstream
  license text (the banner at the top of the reference `fftsg.c` is Tap's
  copy of it; the upstream file carries no notice of its own). There is no
  other license text; in particular the package is neither public domain nor
  under a named open-source license. The full grant reads:

  > You may use, copy, modify this code for any purpose and without fee. You
  > may distribute this ORIGINAL package.

  What that grants: use, copying and modification of the code, for any purpose
  and without fee; and distribution of the *original* package. Distribution of
  a modified derivative is not expressly granted.

- What DspTap ships and relies on: the C++20 port of `rdft`
  (`include/tap/dsp/fft/split_radix.h`, a statement-for-statement
  transliteration). It landed at Stage 2a (tap/DspTap#28) beside the vendored
  C, bit-identical to it for both precisions under the parity gate
  (`tests/test_fft_parity_ooura.cpp`); Stage 2b (tap/DspTap#31) routed
  `basic_real_fft` at it, so every floating transform a consumer runs is the
  port; Stage 2c moved the C out of the shipping tree. The port is a
  **derivative work, not the ORIGINAL package**, and its redistribution
  relies on the **modification grant** ("modify this code for any purpose").
  Precedent exists but is not relied on: WebRTC/Chromium ship a modified
  `fft4g.c` under `common_audio/third_party/ooura/`, and their `LICENSE`
  quotes a broader notice — "You may use, copy, modify and distribute this
  code for any purpose (include commercial use) and without fee. Please refer
  to this package when you modify this code." — that is not in the `fft.tgz`
  readme. DspTap has not traced the origin of that text and does not rely on
  it.
- What DspTap still carries from the package itself: `readme.txt`, and the
  reference `fftsg.c` under `tests/reference/ooura/` for the gate. The
  reference `fftsg.c` is textually identical to the 2006-12-28 `fft.tgz`
  except for a provenance banner Tap added at the top (the upstream file
  carries no notice of its own; the notice lives in `readme.txt`; the banner's
  pointer to the readme was updated when the file moved at 2c) and stripped
  trailing whitespace. That is a partial copy of the original package with
  the notice attached, not a modified transform; the maintainer's reading is
  that it is within the intent of the distribution grant.
- Contact with the author (audit Part 5: "before Stage 2c merges, attempt the
  address in the notice for an explicit statement on derivative distribution
  and record the outcome either way"): **not yet attempted as of
  2026-09-23.** The draft is at `docs/fft-design.md`, "Draft email to the
  author — for the maintainer to send" (addresses
  `ooura@kurims.kyoto-u.ac.jp`, then `ooura@mmm.t.u-tokyo.ac.jp`). Sending it
  is the maintainer's decision and has not been delegated; nothing has been
  sent by anyone. When it is sent, or the decision is taken not to send it,
  this bullet records the date and the outcome (reply / no reply by date).
- SPDX: the port header carries `SPDX-License-Identifier:
  LicenseRef-Ooura AND MIT` — Ooura's notice verbatim as the governing terms
  for the derived portion, MIT for the wrapper and DspTap's additions — plus a
  line stating that it is a derivative work with the modification copyright.
  `LicenseRef-Ooura` denotes the `Copyright:` block of
  `third_party/ooura/readme.txt` (lines 140–145), reproduced verbatim as
  [`LICENSES/LicenseRef-Ooura.txt`](LICENSES/LicenseRef-Ooura.txt) so the
  reference resolves. `fft.h` (the wrapper) and everything else DspTap wrote
  stay plain MIT.
- `readme.txt` stays at `third_party/ooura/readme.txt` permanently, as the
  license record for the derived code and the one fixed path the port
  header's banner cites.
- This is a maintainer judgement call, not legal advice. This file is the
  canonical statement; `README.md` and `docs/fft-design.md` point here.

## CMSIS-DSP and CMSIS-Core (subset)

- Path: `third_party/cmsis-dsp/`.
- Origin: Arm Limited — a minimal subset (float32 real FFT + Helium/MVE
  transform closure) vendored for the optional Cortex-M55 backend.
- License: Apache-2.0 — see `third_party/cmsis-dsp/LICENSE`; every vendored
  file keeps the header it has upstream, which is an SPDX Apache-2.0 header
  in all of them except `PrivateInclude/arm_compiler_specific.h`, which
  carries no license header upstream either (checked against CMSIS-DSP at
  the pinned commit `918014f0ba96`) and is covered by the package's
  `LICENSE`. Provenance and refresh procedure in
  `third_party/cmsis-dsp/VENDOR.md`.
- Compiled only when `TAP_DSP_FFT_CMSIS` is ON (Arm cross builds), into the
  `tap_dsp_fft` library that exists only then; untouched by the default
  desktop / Apple / Hexagon builds, which link no compiled library at all.
