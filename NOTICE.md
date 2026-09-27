# Third-party notices

DspTap's own code — everything outside `third_party/`:
`include/tap/dsp/`, `tests/`, `tools/`, `notebooks/`, `bench/`, `scripts/`,
`cmake/`, `platform/`
(the bare-metal startup file and linker scripts the QEMU legs link, carried
from, or written beside, MuTap's MIT copies; `platform/README.md` records
their origin), `docs/` and the build files — is licensed under the MIT
License — see [`LICENSE`](LICENSE). As of tap/DspTap#42, in the
maintainer's judgement (below), none of it is derived from third-party code;
the trees before #42 carried the C of Ooura's package, code derived from it,
or both, as recorded below ("Which trees carry what").

It bundles the following third-party code, each retaining its own license
text:

## Ooura General Purpose FFT Package (license record only since tap/DspTap#42)

- **What DspTap carries: the license record, no code.**
  `third_party/ooura/readme.txt` (the package's readme: the only upstream
  license text) and `LICENSES/LicenseRef-Ooura.txt` (the readme's
  `Copyright:` block, lines 140–145, verbatim, in the REUSE layout). As of
  tap/DspTap#42 DspTap ships, in the maintainer's judgement (the #42 and
  #39 bullets below), **no code derived from the package**: no source file
  of it, no port, no transcription. Both files stay at these paths
  permanently as the license record for the historical trees consumers
  pinned, which do carry derived code (below): in those trees the SPDX
  reference `LicenseRef-Ooura` resolves to the second file, and a
  consumer's notices can cite one fixed path for either.
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

- **History: what DspTap carried from the package, and until when.** Recorded
  because consumers pinned the trees concerned.
  - *The vendored C.* `fftsg.c` (with DspTap's `fftsg_float.c` wrapper,
    `#define double float` plus renames) shipped in the compiled
    `tap_dsp_fft` library until Stage 2c of
    `docs/audit-fft-and-code-smells.md` (tap/DspTap#32, `8350f13`), which
    made the library header-only; from 2c a test-only reference copy under
    `tests/reference/ooura/` served the parity gate, until Decision D6
    (tap/DspTap#36) deleted it once MuTap and MuTap-Max pinned a tree
    containing 2c. That copy was textually identical to the 2006-12-28
    `fft.tgz` except for a provenance banner Tap added at the top and
    stripped trailing whitespace — a partial copy of the original package
    with the notice attached, which the maintainer read as within the intent
    of the distribution grant. `docs/fft-design.md` ("The bit-identity
    record after D6") records the upstream file's hash.
  - *The port.* `include/tap/dsp/fft/split_radix.h`, a statement-for-statement
    C++20 transliteration of `rdft`, landed at Stage 2a (tap/DspTap#28)
    bit-identical to the vendored C for both precisions, and from Stage 2b
    (tap/DspTap#31, `bbfa48d`) until tap/DspTap#42 every floating transform
    a consumer ran was the port. It was a **derivative work, not the
    ORIGINAL package**, and its redistribution relied on the **modification
    grant** ("modify this code for any purpose"). Its header carried
    `SPDX-License-Identifier: LicenseRef-Ooura AND MIT` — Ooura's notice
    verbatim as the governing terms for the derived portion, MIT for the
    wrapper and DspTap's additions — plus a line stating that it was a
    derivative work with the modification copyright; `fft.h` (the wrapper)
    and everything else DspTap wrote stayed plain MIT. Precedent existed but
    was not relied on: WebRTC/Chromium ship a modified `fft4g.c` under
    `common_audio/third_party/ooura/`, and their `LICENSE` quotes a broader
    notice — "You may use, copy, modify and distribute this code for any
    purpose (include commercial use) and without fee. Please refer to this
    package when you modify this code." — that is not in the `fft.tgz`
    readme; DspTap never traced the origin of that text. tap/DspTap#42
    deleted the port, together with its output-fingerprint pins.
  - *The fixed-point transcription.* From Stage 3b (tap/DspTap#27,
    `b08f6c6`) until tap/DspTap#39 the fixed-point engine's real post-pass,
    the inverse's pre-pass, their coefficient table and the DC/Nyquist
    handling were a transcription of the package's `rftfsub` / `rftbsub` /
    `makect` formulas into the fixed-point arithmetic, under `SPDX: MIT`,
    and those trees fall under the same modification-grant analysis as the
    port. tap/DspTap#39 removed that code and re-derived the pass from the
    published literature (Cooley, Lewis & Welch 1970; Sorensen, Jones,
    Heideman & Burrus 1987) under a recorded clean-room procedure: a
    hand-off commit removed the transcription and every passage stating its
    formulas; the implementer worked from a brief holding the published
    half-length method in DspTap's convention and the numeric contract, was
    barred from the package, the port, the removed text and its history, and
    every other FFT library's source, and reported what it accessed. The
    hand-off left the old table checksums and output fingerprints pinned in
    the tests, so the result was checked against the removed code's bits,
    not obtained blind. The derivation converged on the same arithmetic: an
    arrangement of the published method with the fewest roundings within the
    pass's one-bit growth budget is one complex product per bin pair (by
    (1 + i W_N^k)/2, or equally (1 − i W_N^k)/2), every such arrangement over
    a once-rounded table yields the same integers but for exact half-LSB
    ties, and the result reproduces the previous output bit for bit. The
    text — code, table generator, derivation — was written without access to
    the package or the removed code; the dataflow, statement order, signs and
    table layout are nonetheless the same as `rftfsub` / `rftbsub` /
    `makect`'s, the published method in the arrangement that code also takes.
    `docs/fft-design.md` ("The fixed-point post-pass, re-derived") records the
    procedure, its shortfalls and the access statement. The maintainer's
    judgement is that the fixed-point engine
    (`include/tap/dsp/fft/fixed_point.h`, `include/tap/dsp/fft/tables.h`) is
    DspTap's own, MIT, and not a derivative of the package.
  - *The replacement of the port.* tap/DspTap#42 replaced it with the srdif
    engine (`include/tap/dsp/fft/srdif.h`), written under a recorded
    clean-room procedure: hand-off commits (`17db855`, `05c81f1`) removed the
    port, its per-N bit-exact output pins and every comment naming its
    internal routines or table recurrences; unlike at #39 no per-N oracle
    remained, but one whole-scenario checksum of the port did
    (`bench/README.md`, `0x662dd085b5b88325`, the x86-64 `rfft_f32_512`
    value), in a paragraph the implementer read and rewrote, and the
    replacement does not match it; the implementer worked from the published
    literature, DspTap's own fixed-point engine and a targets sheet of
    numbers, was barred from every copy of the port and the package, the
    history before the hand-off, the provenance records and every other FFT
    library's source, and reported what it accessed. The brief, written by
    the party that had read the port, named the split-radix literature and
    the half-length DIF / bit-reversal / post-pass pipeline, which is also
    the port's; independence is claimed for the arrangement and the text,
    not the family. Every floating output bit changed. What the hand-off did
    not remove is stated rather than implied: history prose naming the port
    and structural facts at that level (a split-radix `rdft`, tables from
    libm); the targets sheet's heap formula (which implies the port's
    two-table layout) and its deterministic cells (three-digit fingerprints
    of the port's output); the port's routine names and verbatim statements
    in files the brief forbade but the worktree still held, fenced by the
    brief alone; and the fixed-point post-pass, which the srdif engine
    generalizes to floating point and which #39 records as arithmetically
    the package's `rftfsub` / `rftbsub`.
    `docs/fft-design.md` ("The floating engine, replaced clean-room") records
    the procedure, the access statement and its near misses, and the
    structural comparison against `fftsg.c` made by #42's provenance review;
    "Comparison with other FFT libraries (2026-09-27)" there records a
    further audit against FFTW, KissFFT, pocketfft, PFFFT, CMSIS-DSP and
    Numerical Recipes (no copying or close following found; not legal
    advice).
    The maintainer's judgement is that the srdif engine is DspTap's own,
    MIT, and not a derivative of the package.
  - *Which trees carry what* (for a consumer's notices): trees before
    tap/DspTap#32 ship `fftsg.c` itself; trees from #32 (`8350f13`) to
    before #36 (D6, `0db95b6`) carry it as a test-only reference copy under
    `tests/reference/ooura/` (`fftsg.c`, `fftsg_float.c`: in the checkout,
    not shipped); trees from #28 to the base of #42 (`7a58ebe`) carry the
    port, shipped by routing from #31; trees from #27 to before #39 also
    carry the fixed-point transcription; trees from #42 on carry only the
    license record above.
- Contact with the author (audit Part 5: "before Stage 2c merges, attempt the
  address in the notice for an explicit statement on derivative distribution
  and record the outcome either way"): **not sent; superseded by the
  replacement** (maintainer decision, 2026-09-26). It had not been attempted
  as of 2026-09-23 (Stage 2c); rather than ask, the maintainer removed the
  derived code (#39, #42), so no derivative remains in what DspTap ships and
  the question has no subject. Nothing was sent, by the maintainer or on the
  maintainer's behalf. The draft stays in `docs/fft-design.md`, "Draft email
  to the author (not sent; superseded by the replacement, 2026-09-26)", as
  the record of what would have been asked.
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
