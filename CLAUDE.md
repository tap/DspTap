# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

**DspTap** — shared DSP primitives for the Tap family of audio libraries. Header-only, plain
portable C++20 (standard library only, no frameworks; `tap::dsp` is a pure INTERFACE target, and
the one compiled library, `tap_dsp_fft`, exists only under `TAP_DSP_FFT_CMSIS` to carry the
vendored CMSIS-DSP objects), consumed as a git submodule by the
individual libraries (TapTools pins it as `submodules/dsptap`; the AmbiTap/MuTap lineage is where
the FFT came from). Seven primitives today: the real FFT (`fft.h`), the YIN pitch detector
(`yin.h`), two pitch shifters (`psola.h`, `pvoc.h`), the wake-word front end —
the log-mel/PCEN feature extractor (`log_mel.h`) and the fixed-ratio decimators to
16 kHz (`decimate.h`) — and the dense/GRU inference kernels (`nn.h`) — plus the FIR substrate carried from
SampleRateTap for the two rate converters (SampleRateTap, RatioTap): Kaiser prototype design
(`kaiser.h`), the sample-format traits (`sample_traits.h`: double/float/Q15/Q31), the FIR dot kernels
(`fir_kernels.h`), row-sum-preserving quantization (`quantize.h`), and the measurement
instruments (`analysis/`). See `README.md` for each asset's contract summary.

## The design discipline (load-bearing — every primitive follows it)

- **Fixed numeric contracts.** Each header documents its packing, conventions, normalization, and
  latency as *numbers*, and the tests pin them. Changing a documented contract point is a breaking
  change for every consumer.
- **Double is the golden model; float32 is the embedded profile; Q15/Q31 are format-limited
  embedded profiles.** `basic_*<Sample>` templates with `using x = basic_x<double>` /
  `x32 = basic_x<float>` aliases. The double path never changes for speed; accelerated float
  backends (vDSP, CMSIS-Helium for the FFT) must re-present the *exact* golden contract, so the
  double test battery stays a valid oracle. Cross-precision agreement is pinned by tests. The
  fixed-point profiles (`sample_traits.h`) carry their contracts as numbers the same way — Q
  formats, the single rounding point, saturation — and exist for M33/M55-class targets
  (Bluetooth-adjacent converters, eurorack/pedal deployments) where double or any float is
  unaffordable. Per-primitive fixed-point adoption is opt-in and is a documented Q-format design
  each time, via traits over raw sample types, never wrapper classes. The real FFT is the
  four-profile ladder in full: `double` / `float` / Q15 / Q31 over one contract, the fixed
  profiles as `basic_real_fft<int16_t|int32_t, Scaling>` returning an exponent, with a
  fixed-scaling and a block-floating policy (`fft.h`, `fft/fixed_point.h`,
  `docs/fft-design.md`), each floor a measured number against the double profile. The floating
  profiles run the C++20 split-radix engine (`fft/split_radix.h`), bit-identical to the vendored
  Ooura C it replaced; since Stage 4 the engine is a template parameter
  (`basic_real_fft<Sample, Policy = detail::default_real_fft_policy_t<Sample>>`: the engine for
  `float`/`double`, defaulting to `default_real_fft_engine_t<Sample>`, the scaling policy for
  Q15/Q31; the accelerated engines under `fft/backends/`, each stating its size range and
  shareability as contract numbers), and
  the build's selection opens an inline-namespace ABI tag on `tap::dsp` (`fft_split_radix` /
  `fft_cmsis` / `fft_vdsp`) that every class embedding the FFT by value lives in; since Stage 2c
  the C is not part of what ships — the reference copy under `tests/reference/ooura/` exists for
  the bit-identity gate (`tests/test_fft_parity_ooura.cpp`) alone and goes once both consumers
  pin a tree past 2c (Decision D6); `third_party/ooura/readme.txt` stays as the license record.
- **Real-time safe by construction.** Geometry fixed at construction, every buffer allocated
  there; processing is `noexcept` and allocation-free. Numerically fragile recursions (e.g. the
  order-48 Levinson–Durbin inside `pvoc`) run in double even in the float profile — documented
  where it happens.
- **Hot loops stay plain.** The expensive inner loops (YIN's difference function above all) are
  written as contiguous, branch-light scalar code on purpose: they are the designated seams for
  future MVE/HVX backends behind the same contract.
- **Implement from published literature only.** This is an IP policy, not a citation habit: cite
  the paper (YIN: de Cheveigné & Kawahara 2002; peak locking: Laroche–Dolson; LPC: standard
  source-filter literature) and never reverse-engineer a shipping product's behavior.
- **Document honest limits.** When an algorithm has a known failure mode (PSOLA thinning pure
  tones, YIN's first-dip rule on missing-fundamental synthetics), the header says so and a test
  *pins* the behavior rather than papering over it.

## Style

`STYLE.md` is the shared Tap house style; `.clang-format` and `.clang-tidy` enforce it and CI runs
both (plus a drift check that the config files match the canonical taphouse copies — never edit
them locally). Run `pre-commit install` once per clone so commits are formatted automatically; on
Claude Code web the checked-in SessionStart hook (`.claude/hooks/session-start.sh`) does this for
you at session start. clang-tidy compiles with a *clang* front end — code that GCC accepts can
still fail there, so treat the tidy job as a second compiler.

## Build & test

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
```

Tests are GoogleTest (FetchContent), typed over `float`/`double` (`TYPED_TEST_SUITE`). Patterns to
preserve: contract tests named for the promise they pin (`RoundTripReproducesInput`,
`PureToneOctaveUpThinsOut`); the certified-geometry parity tests for accelerated FFT backends; and
`tap::dsp::yin` (certified by its own battery) used as the pitch oracle for the shifters.

CI runs the battery on three hosts (linux, windows, macOS with vDSP) and on four bare-metal
Cortex-M legs under `qemu-system-arm` (`ci.yml`, `embedded`): M4 soft-float and M4F
(`cmake/arm-cortex-m4-mps2.cmake`, `TAP_DSP_M4_FPU`), M33 (`cmake/arm-cortex-m33-mps2.cmake`) and
M55 (`cmake/arm-cortex-m55-mps3.cmake`), the one leg where the CMSIS-Helium backend is ON and its
parity suite runs for real. Each leg builds `tap_dsp_tests` one-shot (`tests/bare_metal_main.cpp`)
with the `MAIN_FILTER` in `tests/CMakeLists.txt` naming what the emulation budget excludes, plus
the FFT parity gate at `TAP_DSP_PARITY_MAX_N=4096`. The instruction-count ratchet (`bench.yml`,
`bench/README.md`) runs the same four cores under five baseline keys at ±3 %, with `.text`
ceilings per profile; a red ratchet is a failing check. The hosted build and battery need no Arm
toolchain; where `gcc-arm-none-eabi` and `qemu-system-arm` are installed,
`cmake -S . -B build-m33 -DCMAKE_TOOLCHAIN_FILE=cmake/arm-cortex-m33-mps2.cmake` reproduces a leg.

## Verification layer (C ABI + notebooks)

`tools/capi/` exposes the primitives through a minimal C ABI (`cmake -B build_capi -S tools/capi`);
`notebooks/dsptap_py.py` is the ctypes bridge (builds `build_capi/` on first import), and
`notebooks/pitchshift.ipynb` and `notebooks/fft.ipynb` are committed *executed* — they measure the
shipping C++, not a Python re-implementation (the one deliberate exception is a labeled
naive-phase-vocoder strawman). If you
change a primitive's behavior, re-execute the notebook; if you add a primitive that notebooks
should measure, extend the capi + bridge alongside it.

## Adding a primitive (checklist)

1. `include/tap/dsp/<name>.h` — `basic_<name><Sample>` + aliases, full contract docstring
   (geometry, conventions, latency, limits), MIT SPDX banner.
2. `tests/test_<name>.cpp` — typed battery pinning every contract point, plus float/double
   cross-precision agreement; add to `tests/CMakeLists.txt`. State each profile's numbers
   (scale, floor, saturation-free level, exponent where there is one), and name the on-target
   leg that runs its battery (the M4/M4F/M33/M55 QEMU legs run `tap_dsp_tests` whole unless the
   `MAIN_FILTER` excludes a suite by name).
3. `README.md` section (and bump the primitive count in the intro line).
4. capi + `dsptap_py` exposure if the notebooks need it.
5. Vendored code goes under `third_party/` with license retained and a `NOTICE.md` entry — and is
   excluded from formatting (a reference copy kept for a test alone, like `tests/reference/ooura/`,
   carries its own directory-local `.clang-format` with `DisableFormat: true` and is kept out of the
   compile database the tidy gate reads; the canonical style files are never edited for it).

## Consumers & release flow

Changes land here first, then consumers bump their submodule pin (TapTools pins this repo;
TapTools-Max picks it up transitively through TapTools). After a PR merges via rebase/squash,
repoint any open consumer pins at the identical tree on `main` so they stay reachable after branch
cleanup.
