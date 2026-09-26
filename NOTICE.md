# Third-party notices

DspTap's own code — everything outside `third_party/`, except the portion
of `include/tap/dsp/fft/split_radix.h` derived from Ooura's package (below):
`include/tap/dsp/`, `tests/`, `tools/`, `notebooks/`, `bench/`, `scripts/`,
`cmake/`, `platform/`
(the bare-metal startup file and linker scripts the QEMU legs link, carried
from, or written beside, MuTap's MIT copies; `platform/README.md` records
their origin), `docs/` and the build files — is licensed under the MIT
License — see [`LICENSE`](LICENSE).

It bundles the following third-party code, each retaining its own license
text:

## Ooura General Purpose FFT Package

- Paths: `third_party/ooura/readme.txt` (the package's readme: the only
  upstream license text, kept at this path permanently as the license record
  for the derived code) and `include/tap/dsp/fft/split_radix.h` (the derived
  C++20 port, below). No source file of the package is carried any more:
  since Stage 2c of `docs/audit-fft-and-code-smells.md` the library is
  header-only and no vendored C is part of it (the one compiled library,
  `tap_dsp_fft`, exists only under `TAP_DSP_FFT_CMSIS` and carries the
  CMSIS-DSP objects below), and the test-only reference copy of `fftsg.c`
  (with DspTap's `fftsg_float.c` wrapper) that the parity gate compiled under
  `tests/reference/ooura/` from 2c on was deleted at Decision D6, once both
  MuTap and MuTap-Max pinned a tree containing 2c.
- Author: Takuya Ooura. Copyright(C) 1996-2001 Takuya OOURA.
- License: the package's own terms, stated in `readme.txt`, the only upstream
  license text (the banner that was at the top of the reference `fftsg.c`
  was Tap's copy of it; the upstream file carries no notice of its own).
  There is no other license text; in particular the package is neither
  public domain nor under a named open-source license. The full grant reads:

  > You may use, copy, modify this code for any purpose and without fee. You
  > may distribute this ORIGINAL package.

  What that grants: use, copying and modification of the code, for any purpose
  and without fee; and distribution of the *original* package. Distribution of
  a modified derivative is not expressly granted.

- What DspTap ships and relies on: the C++20 port of `rdft`
  (`include/tap/dsp/fft/split_radix.h`, a statement-for-statement
  transliteration). It landed at Stage 2a (tap/DspTap#28) beside the vendored
  C, bit-identical to it for both precisions under the parity gate that ran
  until D6 (pinned since as output fingerprints,
  `tests/test_fft_split_radix_fingerprint.cpp`); Stage 2b (tap/DspTap#31)
  routed `basic_real_fft` at it, so every floating transform a consumer runs
  is the port; Stage 2c moved the C out of the shipping tree and D6 deleted
  the test-only reference copy. The port is a **derivative work, not the
  ORIGINAL package**, and its redistribution
  relies on the **modification grant** ("modify this code for any purpose").
  Precedent exists but is not relied on: WebRTC/Chromium ship a modified
  `fft4g.c` under `common_audio/third_party/ooura/`, and their `LICENSE`
  quotes a broader notice — "You may use, copy, modify and distribute this
  code for any purpose (include commercial use) and without fee. Please refer
  to this package when you modify this code." — that is not in the `fft.tgz`
  readme. DspTap has not traced the origin of that text and does not rely on
  it.
- What DspTap still carries from the package itself: `readme.txt` alone.
  Until D6 it also carried the reference `fftsg.c` for the gate, textually
  identical to the 2006-12-28 `fft.tgz` except for a provenance banner Tap
  added at the top and stripped trailing whitespace — a partial copy of the
  original package with the notice attached, which the maintainer read as
  within the intent of the distribution grant; with that copy gone, the
  derivative (the port) is the only question the grant has to answer.
  `docs/fft-design.md` ("The bit-identity record after D6") records the
  upstream file's hash and how to re-verify the port against it.
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
- The fixed-point engine (`include/tap/dsp/fft/fixed_point.h`,
  `include/tap/dsp/fft/tables.h`) is DspTap's own, MIT, and carries nothing
  of the package. History, recorded because consumers pinned the trees
  concerned: from Stage 3b (tap/DspTap#27, `b08f6c6`) until tap/DspTap#39,
  its real post-pass, the inverse's pre-pass, their coefficient table and
  the DC/Nyquist handling were a transcription of the package's
  `rftfsub` / `rftbsub` / `makect` formulas into the fixed-point
  arithmetic, and those trees fall under the same modification-grant
  analysis as the port. tap/DspTap#39 removed that code and re-derived the
  pass from the published literature (Cooley, Lewis & Welch 1970; Sorensen,
  Jones, Heideman & Burrus 1987) under a recorded clean-room procedure: a
  hand-off commit removed the transcription and every passage stating its
  formulas; the implementer worked from a brief holding the published
  half-length method in DspTap's convention and the numeric contract, was
  barred from the package, the port, the removed text and its history, and
  every other FFT library's source, and reported what it accessed. The
  derivation converged on the same arithmetic: the arrangement of the
  published method with the fewest roundings within the pass's one-bit
  growth budget is one complex product per bin pair by (1 + i W_N^k)/2, and
  with that coefficient rounded once from its double it reproduces the
  previous output bit for bit. The code, the table generator and the
  derivation were written independently; the arithmetic is the published
  method's. `docs/fft-design.md` ("The fixed-point post-pass, re-derived")
  records the procedure and the access statement. The maintainer's
  judgement is that the fixed-point engine is not a derivative of the
  package.
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
