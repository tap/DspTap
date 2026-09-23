# DspTap benchmarks and the performance ratchet

The FFT is the primitive under MuTap's, SampleRateTap's and RatioTap's own
instruction-count ratchets, so a regression here shows up downstream as a
mysterious chain-level delta someone "updates" past. The primitive needs its
own gate, per profile, per target, seeded before anything changes — which is
why this scaffold lands **before** the C++20 port (Stage 1b of
[`docs/audit-fft-and-code-smells.md`](../docs/audit-fft-and-code-smells.md),
design in its Part 11).

Two instruments, one policy:

| instrument | what it measures | where it runs | role |
|---|---|---|---|
| `icount/` | executed guest instructions per scenario binary, under QEMU's TCG plugin | the QEMU legs (`.github/workflows/bench.yml`) | **the gate**: ±3 %, two-sided |
| `bench_fft.cpp` | wall-clock ns per transform, min-of-N | the host, locally | informational only, never a gate |

## Build and run (host)

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DTAP_DSP_BUILD_BENCH=ON
cmake --build build
./build/bench/tap_dsp_bench_fft                       # wall clock
./build/bench/icount/tap_dsp_icount_rfft_f32_512      # prints its checksum
```

`bench_fft` prints ns per forward and per inverse transform for each scenario
(minimum of 25 reps of 2^20 samples per direction). Only same-machine,
same-load comparisons mean anything; record the machine, compiler, load
average and date with any number that goes into `docs/fft-design.md`.

On the host the `tap_dsp_icount_*` binaries are a smoke test, registered with
ctest (`tap_dsp_icount_<scenario>_smoke`, `bench/icount/smoke.cmake`): each
binary runs twice, both runs must report `ok=1`, and the two `DONE` lines
must be byte-identical. There is no instruction count on the host; the count
comes from QEMU on the legs.

## Scenarios

One binary per scenario, selected by compile definitions (bare metal has no
argv), all built from `icount/icount_main.cpp`:

| scenario | profile | N | what is measured | targets | gated |
|---|---|---|---|---|---|
| `rfft_f32_512` | float | 512 | `basic_real_fft<float>` as built: the split-radix engine since Stage 2b; CMSIS-DSP Helium behind it on `m55` | all legs | yes |
| `rfft_f32_2048` | float | 2048 | as above | all legs | yes |
| `rfft_f64_512` | double | 512 | `basic_real_fft<double>`: the split-radix engine | host-class only: soft-float double is not a profile, so the bare-metal legs neither build nor baseline it | no key counts it |
| `rfft_q15_512` | Q15 (`std::int16_t`, `scaling::fixed`) | 512 | `basic_real_fft<std::int16_t>`: the int32 fixed-point kernel, Q15 I/O | all legs (fixed point is a profile on the soft-float M4 too) | yes, since Stage 2c |
| `rfft_q31_512` | Q31 (`std::int32_t`, `scaling::fixed`) | 512 | `basic_real_fft<std::int32_t>` | all legs | yes, since Stage 2c |
| `rfft_q31_2048` | Q31 | 2048 | as above | all legs | yes, since Stage 2c |

The three fixed-point scenarios are the ones Stage 3b specified (Part 11),
seeded at Stage 2c, so every later `SMMULR`, Helium or table-layout change
to the Q15/Q31 kernel has a number to beat. They run under `scaling::fixed`
only: its count is data-independent (no headroom scan), which is what a
ratchet wants; block floating point is measured by the battery's pins, not
here.

The bench measures only what ships. From Stage 2a until Stage 2c a second
binary per float scenario (`_port`, then `_c` after the flip; the
`INFORMATIONAL_SUFFIX` of `scripts/icount.py`) counted the vendored Ooura C
beside the shipping engine and the job log printed the ratio and whether
the two `DONE` checksums agreed, never gating on it; that is how Stage 2b
was judged (the tables below keep the record). Stage 2c moved the C out of
the shipping tree — the reference copy under `tests/reference/ooura/`
serves the bit-identity gate alone — and retired the selector
(`TAP_DSP_BENCH_ENGINE`), the adapter and the `_c` pairs. The `DONE` line
still prints `engine=basic_real_fft`, the constant the shipping class always
printed, so the recorded counts (which include the print) did not move.

Each scenario constructs one transform (`tap::dsp::basic_real_fft<Sample>`),
then runs a forward + inverse loop over a four-block xorshift corpus for 2^20 samples per
direction (2048 iterations at N = 512, 512 at N = 2048). Every output word
of every iteration — the spectrum and the scaled inverse, or for the
fixed-point scenarios the spectrum, the unnormalized inverse and the
exponent each transform returned — goes through an
integer FNV-1a-64 fold over its bit pattern (`bench_common.h`), printed at
the end:

```
TAP_DSP_ICOUNT_DONE ok=1 engine=basic_real_fft backend=split_radix scenario=rfft_f32_512 checksum=0x662dd085b5b88325
```

The fold is exact and order-sensitive: two runs of the same binary print the
same value, and a 1-ulp change in any single output changes it (verified by
nudging one spectrum bin at one iteration with `nextafterf`: `0x662dd085b5b88325`
becomes `0x259510dc8721cba8`). A floating running sum cannot promise that — it
absorbs differences below the accumulator's ulp — which is why the checksum
is an integer hash: it was the fingerprint the job compared between the C
and the port on the QEMU legs, where the host parity TU does not run (the
checksum above is the same value the C printed before the flip: the port is
bit-identical), and it is what a same-key checksum change means today: an
output bit moved. `ok` is a sanity check that the last iteration
round-trips its input — within 1e-3 in the floating scenarios' own
precision; within eight output LSB referred to the input, i.e.
|out · 2^s − x| ≤ 8 · 2^s with s = e_fwd + e_inv + 1 − log₂N, in the
fixed-point ones (the contract's round-trip identity, `fft.h`; eight is the
header's pinned Q31 single-transform worst case rounded to a power of two,
against 0.52 / 3.6 LSB measured in this workload for Q15 / Q31); `backend=`
names what the build routed `basic_real_fft` through: `split_radix`,
`cmsis` on the `m55` key, `accelerate` on Apple, `fixed_point` for Q15 /
Q31.

Why these numbers: on the host (x86-64, GCC 13 `-O2`, callgrind) the three
floating scenarios execute 109 M, 125 M and 108 M instructions (the three
fixed-point ones 434 M, 407 M and 495 M at `-O3`, of which the kernel's two
radix-4 stage functions are 71-80 %; construction is under 0.5 %).
Construction plus the one-time table build plus the print is under
0.2 M — well under the 1 % the design asks for. The share of the count that
is not the transform itself — the class's out-of-place copies, the 2/N
scaling loop and the fold — is 16.1 % / 14.0 % / 16.6 %
(`main` inclusive minus the transform's inclusive count, measured on the C
before the flip); it is constant per
scenario, so the ratchet works, but a 3 % change in the FFT alone reads as
roughly 3 % × (1 − that share) at the gate, and the share is larger on a
scalar Cortex-M. The Arm counts are larger than the host's; the ratios
hold. No `<random>` (a toolchain's libstdc++ would move the count), no
allocation in the loop, and no double anywhere in the float scenarios (the
M4 soft-float leg would otherwise measure libgcc).

A build directory configured before Stage 2c may still carry the retired
`TAP_DSP_BENCH_ENGINE` cache variable; it is read by nothing now and a fresh
configure ignores it. The port/C ratio the retired pairs printed was what
Stage 2b was judged on (the port within ±3 % of the C on every QEMU leg,
audit Part 3); the outcome is in the table below.

## Baselines and targets

[`baselines.json`](baselines.json) is keyed by target, one dict per key:

| key | core | QEMU machine | float32 FFT |
|---|---|---|---|
| `m4-softfp` | Cortex-M4, `-mfloat-abi=soft` | `mps2-an386` | split-radix engine (every float op a library call) |
| `m4f` | Cortex-M4F, `fpv4-sp-d16` | `mps2-an386` | split-radix engine |
| `m33` | Cortex-M33, no MVE | `mps2-an505` | split-radix engine |
| `m55` | Cortex-M55, Helium | `mps3-an547` | CMSIS-DSP — the deployed profile |
| `m55-ooura` | Cortex-M55, `-DTAP_DSP_FFT_CMSIS=OFF` | `mps3-an547` | split-radix engine — the fallback (the key keeps its pre-flip name; the engine is bit-identical to the Ooura C it names) |

JSON carries no comments, so the provenance of every recorded set lives here,
in the table below: the `main` run that recorded it, and the GCC and QEMU
versions, because the counts are only comparable within a toolchain/QEMU
pair. The plugin header is pinned to the QEMU that Ubuntu 24.04 ships
(8.2.2; that header defines `QEMU_PLUGIN_VERSION 1`, plugin API v1) and
digest-verified on download, as MuTap does.

### Recorded baselines and updates

Every seeding or `--update` commit adds one row per key it touched. The
shape is fixed so the record stays greppable:

| key | scenario | before | after | delta | reason | main run (URL) | main SHA | arm-none-eabi-gcc | qemu-system-arm |
|---|---|---|---|---|---|---|---|---|---|
| m33 | rfft_f32_2048 | — | 116,385,409 | — | **seeded** from the first push-to-main run after #21 merged (vendored Ooura C, `TAP_DSP_BENCH_ENGINE=reference_c`; m55 = CMSIS-DSP Helium, m55-ooura = the Ooura fallback) | [run 35281280300](https://github.com/tap/DspTap/actions/runs/35281280300) | df482d1 | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m33 | rfft_f32_512 | — | 102,248,169 | — | **seeded** from the first push-to-main run after #21 merged (vendored Ooura C, `TAP_DSP_BENCH_ENGINE=reference_c`; m55 = CMSIS-DSP Helium, m55-ooura = the Ooura fallback) | [run 35281280300](https://github.com/tap/DspTap/actions/runs/35281280300) | df482d1 | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m4-softfp | rfft_f32_2048 | — | 2,296,984,479 | — | **seeded** from the first push-to-main run after #21 merged (vendored Ooura C, `TAP_DSP_BENCH_ENGINE=reference_c`; m55 = CMSIS-DSP Helium, m55-ooura = the Ooura fallback) | [run 35281280300](https://github.com/tap/DspTap/actions/runs/35281280300) | df482d1 | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m4-softfp | rfft_f32_512 | — | 1,868,441,244 | — | **seeded** from the first push-to-main run after #21 merged (vendored Ooura C, `TAP_DSP_BENCH_ENGINE=reference_c`; m55 = CMSIS-DSP Helium, m55-ooura = the Ooura fallback) | [run 35281280300](https://github.com/tap/DspTap/actions/runs/35281280300) | df482d1 | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m4f | rfft_f32_2048 | — | 111,859,257 | — | **seeded** from the first push-to-main run after #21 merged (vendored Ooura C, `TAP_DSP_BENCH_ENGINE=reference_c`; m55 = CMSIS-DSP Helium, m55-ooura = the Ooura fallback) | [run 35281280300](https://github.com/tap/DspTap/actions/runs/35281280300) | df482d1 | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m4f | rfft_f32_512 | — | 98,090,666 | — | **seeded** from the first push-to-main run after #21 merged (vendored Ooura C, `TAP_DSP_BENCH_ENGINE=reference_c`; m55 = CMSIS-DSP Helium, m55-ooura = the Ooura fallback) | [run 35281280300](https://github.com/tap/DspTap/actions/runs/35281280300) | df482d1 | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m55 | rfft_f32_2048 | — | 54,858,120 | — | **seeded** from the first push-to-main run after #21 merged (vendored Ooura C, `TAP_DSP_BENCH_ENGINE=reference_c`; m55 = CMSIS-DSP Helium, m55-ooura = the Ooura fallback) | [run 35281280300](https://github.com/tap/DspTap/actions/runs/35281280300) | df482d1 | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m55 | rfft_f32_512 | — | 52,382,331 | — | **seeded** from the first push-to-main run after #21 merged (vendored Ooura C, `TAP_DSP_BENCH_ENGINE=reference_c`; m55 = CMSIS-DSP Helium, m55-ooura = the Ooura fallback) | [run 35281280300](https://github.com/tap/DspTap/actions/runs/35281280300) | df482d1 | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m55-ooura | rfft_f32_2048 | — | 107,806,480 | — | **seeded** from the first push-to-main run after #21 merged (vendored Ooura C, `TAP_DSP_BENCH_ENGINE=reference_c`; m55 = CMSIS-DSP Helium, m55-ooura = the Ooura fallback) | [run 35281280300](https://github.com/tap/DspTap/actions/runs/35281280300) | df482d1 | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m55-ooura | rfft_f32_512 | — | 94,561,954 | — | **seeded** from the first push-to-main run after #21 merged (vendored Ooura C, `TAP_DSP_BENCH_ENGINE=reference_c`; m55 = CMSIS-DSP Helium, m55-ooura = the Ooura fallback) | [run 35281280300](https://github.com/tap/DspTap/actions/runs/35281280300) | df482d1 | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |

| m33 | rfft_f32_2048 | 116,385,409 | 115,465,626 | -0.79 % | **re-recorded at the Stage 2b flip** (tap/DspTap#31): the bare key now measures the split-radix engine that ships through `basic_real_fft`; the C it replaced is the `_c` informational sibling of the same run (output checksums identical on every key). Per Part 11 / D11 the exact new count is recorded, never absorbed: m55-ooura sat outside the +-3% band (IMPROVED), the M4 and M33 keys inside it, and all four move together | [run 35844483811](https://github.com/tap/DspTap/actions/runs/35844483811) (workflow_dispatch on the PR branch at `b078a80`; the PR's compare-mode run confirms +0.00 %) | `bbfa48d` (the #31 squash on `main`, filled in by tap/DspTap Stage 2c per #31's after-merge list) | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m33 | rfft_f32_512 | 102,248,169 | 100,945,841 | -1.27 % | **re-recorded at the Stage 2b flip** (tap/DspTap#31): the bare key now measures the split-radix engine that ships through `basic_real_fft`; the C it replaced is the `_c` informational sibling of the same run (output checksums identical on every key). Per Part 11 / D11 the exact new count is recorded, never absorbed: m55-ooura sat outside the +-3% band (IMPROVED), the M4 and M33 keys inside it, and all four move together | [run 35844483811](https://github.com/tap/DspTap/actions/runs/35844483811) (workflow_dispatch on the PR branch at `b078a80`; the PR's compare-mode run confirms +0.00 %) | `bbfa48d` (the #31 squash on `main`, filled in by tap/DspTap Stage 2c per #31's after-merge list) | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m4-softfp | rfft_f32_2048 | 2,296,984,479 | 2,294,360,259 | -0.11 % | **re-recorded at the Stage 2b flip** (tap/DspTap#31): the bare key now measures the split-radix engine that ships through `basic_real_fft`; the C it replaced is the `_c` informational sibling of the same run (output checksums identical on every key). Per Part 11 / D11 the exact new count is recorded, never absorbed: m55-ooura sat outside the +-3% band (IMPROVED), the M4 and M33 keys inside it, and all four move together | [run 35844483811](https://github.com/tap/DspTap/actions/runs/35844483811) (workflow_dispatch on the PR branch at `b078a80`; the PR's compare-mode run confirms +0.00 %) | `bbfa48d` (the #31 squash on `main`, filled in by tap/DspTap Stage 2c per #31's after-merge list) | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m4-softfp | rfft_f32_512 | 1,868,441,244 | 1,864,929,141 | -0.19 % | **re-recorded at the Stage 2b flip** (tap/DspTap#31): the bare key now measures the split-radix engine that ships through `basic_real_fft`; the C it replaced is the `_c` informational sibling of the same run (output checksums identical on every key). Per Part 11 / D11 the exact new count is recorded, never absorbed: m55-ooura sat outside the +-3% band (IMPROVED), the M4 and M33 keys inside it, and all four move together | [run 35844483811](https://github.com/tap/DspTap/actions/runs/35844483811) (workflow_dispatch on the PR branch at `b078a80`; the PR's compare-mode run confirms +0.00 %) | `bbfa48d` (the #31 squash on `main`, filled in by tap/DspTap Stage 2c per #31's after-merge list) | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m4f | rfft_f32_2048 | 111,859,257 | 111,278,416 | -0.52 % | **re-recorded at the Stage 2b flip** (tap/DspTap#31): the bare key now measures the split-radix engine that ships through `basic_real_fft`; the C it replaced is the `_c` informational sibling of the same run (output checksums identical on every key). Per Part 11 / D11 the exact new count is recorded, never absorbed: m55-ooura sat outside the +-3% band (IMPROVED), the M4 and M33 keys inside it, and all four move together | [run 35844483811](https://github.com/tap/DspTap/actions/runs/35844483811) (workflow_dispatch on the PR branch at `b078a80`; the PR's compare-mode run confirms +0.00 %) | `bbfa48d` (the #31 squash on `main`, filled in by tap/DspTap Stage 2c per #31's after-merge list) | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m4f | rfft_f32_512 | 98,090,666 | 97,138,544 | -0.97 % | **re-recorded at the Stage 2b flip** (tap/DspTap#31): the bare key now measures the split-radix engine that ships through `basic_real_fft`; the C it replaced is the `_c` informational sibling of the same run (output checksums identical on every key). Per Part 11 / D11 the exact new count is recorded, never absorbed: m55-ooura sat outside the +-3% band (IMPROVED), the M4 and M33 keys inside it, and all four move together | [run 35844483811](https://github.com/tap/DspTap/actions/runs/35844483811) (workflow_dispatch on the PR branch at `b078a80`; the PR's compare-mode run confirms +0.00 %) | `bbfa48d` (the #31 squash on `main`, filled in by tap/DspTap Stage 2c per #31's after-merge list) | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m55-ooura | rfft_f32_2048 | 107,806,480 | 102,784,096 | -4.66 % | **re-recorded at the Stage 2b flip** (tap/DspTap#31): the bare key now measures the split-radix engine that ships through `basic_real_fft`; the C it replaced is the `_c` informational sibling of the same run (output checksums identical on every key). Per Part 11 / D11 the exact new count is recorded, never absorbed: m55-ooura sat outside the +-3% band (IMPROVED), the M4 and M33 keys inside it, and all four move together | [run 35844483811](https://github.com/tap/DspTap/actions/runs/35844483811) (workflow_dispatch on the PR branch at `b078a80`; the PR's compare-mode run confirms +0.00 %) | `bbfa48d` (the #31 squash on `main`, filled in by tap/DspTap Stage 2c per #31's after-merge list) | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m55-ooura | rfft_f32_512 | 94,561,954 | 89,276,321 | -5.59 % | **re-recorded at the Stage 2b flip** (tap/DspTap#31): the bare key now measures the split-radix engine that ships through `basic_real_fft`; the C it replaced is the `_c` informational sibling of the same run (output checksums identical on every key). Per Part 11 / D11 the exact new count is recorded, never absorbed: m55-ooura sat outside the +-3% band (IMPROVED), the M4 and M33 keys inside it, and all four move together | [run 35844483811](https://github.com/tap/DspTap/actions/runs/35844483811) (workflow_dispatch on the PR branch at `b078a80`; the PR's compare-mode run confirms +0.00 %) | `bbfa48d` (the #31 squash on `main`, filled in by tap/DspTap Stage 2c per #31's after-merge list) | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m55 | rfft_f32_2048 | 54,858,120 | 54,858,120 | +0.00 % | **not re-recorded**: CMSIS-DSP Helium is untouched by the flip; the same run measured 54,858,195 (+75 instructions, +0.00 %: the wrapper lost its two unused Ooura workspace vectors and the DONE line's `engine=` string grew by three characters) | run 35844483811, compare mode | — | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m55 | rfft_f32_512 | 52,382,331 | 52,382,331 | +0.00 % | as above; measured 52,382,366 (+35 instructions, +0.00 %) | run 35844483811, compare mode | — | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |

| m4-softfp | rfft_q15_512 | — | 746,311,273 | — | **seeded at Stage 2c** (the Q15/Q31 scenarios Stage 3b specified, Part 11; `basic_real_fft<std::int16_t \| std::int32_t>` under `scaling::fixed`, the int32 kernel on every key: fixed point has no backend, so the two M55 keys measure the same binary) | [run 35861317022](https://github.com/tap/DspTap/actions/runs/35861317022) (workflow_dispatch on the PR branch at `21b3488`; the PR's compare-mode run confirms +0.00 %) | **pending** — the Stage 2c squash SHA on `main`, filled in by the follow-up that pins it (the MuTap bump) | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m4-softfp | rfft_q31_512 | — | 709,161,574 | — | **seeded at Stage 2c** (the Q15/Q31 scenarios Stage 3b specified, Part 11; `basic_real_fft<std::int16_t \| std::int32_t>` under `scaling::fixed`, the int32 kernel on every key: fixed point has no backend, so the two M55 keys measure the same binary) | [run 35861317022](https://github.com/tap/DspTap/actions/runs/35861317022) (workflow_dispatch on the PR branch at `21b3488`; the PR's compare-mode run confirms +0.00 %) | **pending** — the Stage 2c squash SHA on `main`, filled in by the follow-up that pins it (the MuTap bump) | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m4-softfp | rfft_q31_2048 | — | 869,187,566 | — | **seeded at Stage 2c** (the Q15/Q31 scenarios Stage 3b specified, Part 11; `basic_real_fft<std::int16_t \| std::int32_t>` under `scaling::fixed`, the int32 kernel on every key: fixed point has no backend, so the two M55 keys measure the same binary) | [run 35861317022](https://github.com/tap/DspTap/actions/runs/35861317022) (workflow_dispatch on the PR branch at `21b3488`; the PR's compare-mode run confirms +0.00 %) | **pending** — the Stage 2c squash SHA on `main`, filled in by the follow-up that pins it (the MuTap bump) | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m4f | rfft_q15_512 | — | 753,215,637 | — | **seeded at Stage 2c** (the Q15/Q31 scenarios Stage 3b specified, Part 11; `basic_real_fft<std::int16_t \| std::int32_t>` under `scaling::fixed`, the int32 kernel on every key: fixed point has no backend, so the two M55 keys measure the same binary) | [run 35861317022](https://github.com/tap/DspTap/actions/runs/35861317022) (workflow_dispatch on the PR branch at `21b3488`; the PR's compare-mode run confirms +0.00 %) | **pending** — the Stage 2c squash SHA on `main`, filled in by the follow-up that pins it (the MuTap bump) | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m4f | rfft_q31_512 | — | 715,019,447 | — | **seeded at Stage 2c** (the Q15/Q31 scenarios Stage 3b specified, Part 11; `basic_real_fft<std::int16_t \| std::int32_t>` under `scaling::fixed`, the int32 kernel on every key: fixed point has no backend, so the two M55 keys measure the same binary) | [run 35861317022](https://github.com/tap/DspTap/actions/runs/35861317022) (workflow_dispatch on the PR branch at `21b3488`; the PR's compare-mode run confirms +0.00 %) | **pending** — the Stage 2c squash SHA on `main`, filled in by the follow-up that pins it (the MuTap bump) | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m4f | rfft_q31_2048 | — | 876,486,909 | — | **seeded at Stage 2c** (the Q15/Q31 scenarios Stage 3b specified, Part 11; `basic_real_fft<std::int16_t \| std::int32_t>` under `scaling::fixed`, the int32 kernel on every key: fixed point has no backend, so the two M55 keys measure the same binary) | [run 35861317022](https://github.com/tap/DspTap/actions/runs/35861317022) (workflow_dispatch on the PR branch at `21b3488`; the PR's compare-mode run confirms +0.00 %) | **pending** — the Stage 2c squash SHA on `main`, filled in by the follow-up that pins it (the MuTap bump) | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m33 | rfft_q15_512 | — | 752,189,619 | — | **seeded at Stage 2c** (the Q15/Q31 scenarios Stage 3b specified, Part 11; `basic_real_fft<std::int16_t \| std::int32_t>` under `scaling::fixed`, the int32 kernel on every key: fixed point has no backend, so the two M55 keys measure the same binary) | [run 35861317022](https://github.com/tap/DspTap/actions/runs/35861317022) (workflow_dispatch on the PR branch at `21b3488`; the PR's compare-mode run confirms +0.00 %) | **pending** — the Stage 2c squash SHA on `main`, filled in by the follow-up that pins it (the MuTap bump) | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m33 | rfft_q31_512 | — | 714,626,257 | — | **seeded at Stage 2c** (the Q15/Q31 scenarios Stage 3b specified, Part 11; `basic_real_fft<std::int16_t \| std::int32_t>` under `scaling::fixed`, the int32 kernel on every key: fixed point has no backend, so the two M55 keys measure the same binary) | [run 35861317022](https://github.com/tap/DspTap/actions/runs/35861317022) (workflow_dispatch on the PR branch at `21b3488`; the PR's compare-mode run confirms +0.00 %) | **pending** — the Stage 2c squash SHA on `main`, filled in by the follow-up that pins it (the MuTap bump) | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m33 | rfft_q31_2048 | — | 875,801,377 | — | **seeded at Stage 2c** (the Q15/Q31 scenarios Stage 3b specified, Part 11; `basic_real_fft<std::int16_t \| std::int32_t>` under `scaling::fixed`, the int32 kernel on every key: fixed point has no backend, so the two M55 keys measure the same binary) | [run 35861317022](https://github.com/tap/DspTap/actions/runs/35861317022) (workflow_dispatch on the PR branch at `21b3488`; the PR's compare-mode run confirms +0.00 %) | **pending** — the Stage 2c squash SHA on `main`, filled in by the follow-up that pins it (the MuTap bump) | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m55 | rfft_q15_512 | — | 684,157,511 | — | **seeded at Stage 2c** (the Q15/Q31 scenarios Stage 3b specified, Part 11; `basic_real_fft<std::int16_t \| std::int32_t>` under `scaling::fixed`, the int32 kernel on every key: fixed point has no backend, so the two M55 keys measure the same binary) | [run 35861317022](https://github.com/tap/DspTap/actions/runs/35861317022) (workflow_dispatch on the PR branch at `21b3488`; the PR's compare-mode run confirms +0.00 %) | **pending** — the Stage 2c squash SHA on `main`, filled in by the follow-up that pins it (the MuTap bump) | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m55 | rfft_q31_512 | — | 657,595,969 | — | **seeded at Stage 2c** (the Q15/Q31 scenarios Stage 3b specified, Part 11; `basic_real_fft<std::int16_t \| std::int32_t>` under `scaling::fixed`, the int32 kernel on every key: fixed point has no backend, so the two M55 keys measure the same binary) | [run 35861317022](https://github.com/tap/DspTap/actions/runs/35861317022) (workflow_dispatch on the PR branch at `21b3488`; the PR's compare-mode run confirms +0.00 %) | **pending** — the Stage 2c squash SHA on `main`, filled in by the follow-up that pins it (the MuTap bump) | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m55 | rfft_q31_2048 | — | 806,142,635 | — | **seeded at Stage 2c** (the Q15/Q31 scenarios Stage 3b specified, Part 11; `basic_real_fft<std::int16_t \| std::int32_t>` under `scaling::fixed`, the int32 kernel on every key: fixed point has no backend, so the two M55 keys measure the same binary) | [run 35861317022](https://github.com/tap/DspTap/actions/runs/35861317022) (workflow_dispatch on the PR branch at `21b3488`; the PR's compare-mode run confirms +0.00 %) | **pending** — the Stage 2c squash SHA on `main`, filled in by the follow-up that pins it (the MuTap bump) | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m55-ooura | rfft_q15_512 | — | 684,157,511 | — | **seeded at Stage 2c** (the Q15/Q31 scenarios Stage 3b specified, Part 11; `basic_real_fft<std::int16_t \| std::int32_t>` under `scaling::fixed`, the int32 kernel on every key: fixed point has no backend, so the two M55 keys measure the same binary) | [run 35861317022](https://github.com/tap/DspTap/actions/runs/35861317022) (workflow_dispatch on the PR branch at `21b3488`; the PR's compare-mode run confirms +0.00 %) | **pending** — the Stage 2c squash SHA on `main`, filled in by the follow-up that pins it (the MuTap bump) | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m55-ooura | rfft_q31_512 | — | 657,595,969 | — | **seeded at Stage 2c** (the Q15/Q31 scenarios Stage 3b specified, Part 11; `basic_real_fft<std::int16_t \| std::int32_t>` under `scaling::fixed`, the int32 kernel on every key: fixed point has no backend, so the two M55 keys measure the same binary) | [run 35861317022](https://github.com/tap/DspTap/actions/runs/35861317022) (workflow_dispatch on the PR branch at `21b3488`; the PR's compare-mode run confirms +0.00 %) | **pending** — the Stage 2c squash SHA on `main`, filled in by the follow-up that pins it (the MuTap bump) | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m55-ooura | rfft_q31_2048 | — | 806,142,635 | — | **seeded at Stage 2c** (the Q15/Q31 scenarios Stage 3b specified, Part 11; `basic_real_fft<std::int16_t \| std::int32_t>` under `scaling::fixed`, the int32 kernel on every key: fixed point has no backend, so the two M55 keys measure the same binary) | [run 35861317022](https://github.com/tap/DspTap/actions/runs/35861317022) (workflow_dispatch on the PR branch at `21b3488`; the PR's compare-mode run confirms +0.00 %) | **pending** — the Stage 2c squash SHA on `main`, filled in by the follow-up that pins it (the MuTap bump) | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |

`before` is `—` for a seed. `main SHA` is the commit on `main` whose push
run measured the numbers: a pull-request head SHA stops resolving after
this repo's rebase/squash + branch-delete flow, so a seed or update is never
taken from a PR run. The one exception is a re-record forced by a routing
change, where the numbers cannot be measured on `main` before the change is
on `main`: the Stage 2b rows above were seeded by a `workflow_dispatch` on
the PR branch (the branch SHA is named; it stops resolving after the squash)
and the PR's own compare-mode run is what proves the committed numbers are
the branch's, at +0.00 %. The `main SHA` column then names the squash —
which cannot be known while the PR is open, so the cell reads **pending**
until the follow-up that pins the squash (the MuTap bump) fills it in.

The Stage 2b re-record, read against the ±3 % band the C baselines
defined (this run's `_c` siblings reproduced the seeded C counts to within
+28 … +34 instructions, i.e. the adapter's own prologue; the C column below
is the seeded baseline):

| key | `rfft_f32_512`: C → port | ratio | `rfft_f32_2048`: C → port | ratio | checksums |
|---|---|---|---|---|---|
| m4-softfp | 1,868,441,244 → 1,864,929,141 | 0.9981 | 2,296,984,479 → 2,294,360,259 | 0.9989 | identical |
| m4f | 98,090,666 → 97,138,544 | 0.9903 | 111,859,257 → 111,278,416 | 0.9948 | identical |
| m33 | 102,248,169 → 100,945,841 | 0.9873 | 116,385,409 → 115,465,626 | 0.9921 | identical |
| m55-ooura | 94,561,954 → 89,276,321 | 0.9441 | 107,806,480 → 102,784,096 | 0.9534 | identical |
| m55 (CMSIS, not re-recorded) | 52,382,331 → 52,382,366 | 1.0000 | 54,858,120 → 54,858,195 | 1.0000 | n/a (CMSIS vs CMSIS; the `_c` sibling is the Ooura C and differs, as expected) |

### Sizes: recorded `.text` and ceilings

One row per (key, probe) each time a ceiling is recorded or re-recorded; the
policy is under "Size runs in the same job" below. `.text` is the row of
`arm-none-eabi-size -A` on the MinSizeRel probe; the ceiling is what
`bench.yml` carries for that key and profile.

| key | probe | `.text` (bytes) | ceiling | reason | run (URL) | arm-none-eabi-gcc |
|---|---|---|---|---|---|---|
| m4-softfp | `rfft_f32_512` | 53,289 | 54,912 | **recorded at Stage 2c**: measured + 3 %, rounded up to 64 bytes | [run 35861317022](https://github.com/tap/DspTap/actions/runs/35861317022) (`21b3488`; MinSizeRel, `size -A` `.text`) | 13.2.1 (15:13.2.rel1-2) |
| m4-softfp | `rfft_q15_512` | 31,473 | 32,448 | **recorded at Stage 2c**: measured + 3 %, rounded up to 64 bytes | [run 35861317022](https://github.com/tap/DspTap/actions/runs/35861317022) (`21b3488`; MinSizeRel, `size -A` `.text`) | 13.2.1 (15:13.2.rel1-2) |
| m4-softfp | `rfft_q31_512` | 31,001 | 31,936 | **recorded at Stage 2c**: measured + 3 %, rounded up to 64 bytes | [run 35861317022](https://github.com/tap/DspTap/actions/runs/35861317022) (`21b3488`; MinSizeRel, `size -A` `.text`) | 13.2.1 (15:13.2.rel1-2) |
| m4f | `rfft_f32_512` | 44,601 | 45,952 | **recorded at Stage 2c**: measured + 3 %, rounded up to 64 bytes | [run 35861317022](https://github.com/tap/DspTap/actions/runs/35861317022) (`21b3488`; MinSizeRel, `size -A` `.text`) | 13.2.1 (15:13.2.rel1-2) |
| m4f | `rfft_q15_512` | 31,681 | 32,640 | **recorded at Stage 2c**: measured + 3 %, rounded up to 64 bytes | [run 35861317022](https://github.com/tap/DspTap/actions/runs/35861317022) (`21b3488`; MinSizeRel, `size -A` `.text`) | 13.2.1 (15:13.2.rel1-2) |
| m4f | `rfft_q31_512` | 31,209 | 32,192 | **recorded at Stage 2c**: measured + 3 %, rounded up to 64 bytes | [run 35861317022](https://github.com/tap/DspTap/actions/runs/35861317022) (`21b3488`; MinSizeRel, `size -A` `.text`) | 13.2.1 (15:13.2.rel1-2) |
| m33 | `rfft_f32_512` | 44,001 | 45,376 | **recorded at Stage 2c**: measured + 3 %, rounded up to 64 bytes | [run 35861317022](https://github.com/tap/DspTap/actions/runs/35861317022) (`21b3488`; MinSizeRel, `size -A` `.text`) | 13.2.1 (15:13.2.rel1-2) |
| m33 | `rfft_q15_512` | 31,105 | 32,064 | **recorded at Stage 2c**: measured + 3 %, rounded up to 64 bytes | [run 35861317022](https://github.com/tap/DspTap/actions/runs/35861317022) (`21b3488`; MinSizeRel, `size -A` `.text`) | 13.2.1 (15:13.2.rel1-2) |
| m33 | `rfft_q31_512` | 30,633 | 31,552 | **recorded at Stage 2c**: measured + 3 %, rounded up to 64 bytes | [run 35861317022](https://github.com/tap/DspTap/actions/runs/35861317022) (`21b3488`; MinSizeRel, `size -A` `.text`) | 13.2.1 (15:13.2.rel1-2) |
| m55 | `rfft_f32_512` | 107,657 | 110,912 | **recorded at Stage 2c**: measured + 3 %, rounded up to 64 bytes | [run 35861317022](https://github.com/tap/DspTap/actions/runs/35861317022) (`21b3488`; MinSizeRel, `size -A` `.text`) | 13.2.1 (15:13.2.rel1-2) |
| m55 | `rfft_q15_512` | 27,097 | 27,968 | **recorded at Stage 2c**: measured + 3 %, rounded up to 64 bytes | [run 35861317022](https://github.com/tap/DspTap/actions/runs/35861317022) (`21b3488`; MinSizeRel, `size -A` `.text`) | 13.2.1 (15:13.2.rel1-2) |
| m55 | `rfft_q31_512` | 26,593 | 27,392 | **recorded at Stage 2c**: measured + 3 %, rounded up to 64 bytes | [run 35861317022](https://github.com/tap/DspTap/actions/runs/35861317022) (`21b3488`; MinSizeRel, `size -A` `.text`) | 13.2.1 (15:13.2.rel1-2) |
| m55-ooura | `rfft_f32_512` | 39,281 | 40,512 | **recorded at Stage 2c**: measured + 3 %, rounded up to 64 bytes | [run 35861317022](https://github.com/tap/DspTap/actions/runs/35861317022) (`21b3488`; MinSizeRel, `size -A` `.text`) | 13.2.1 (15:13.2.rel1-2) |
| m55-ooura | `rfft_q15_512` | 27,097 | 27,968 | **recorded at Stage 2c**: measured + 3 %, rounded up to 64 bytes | [run 35861317022](https://github.com/tap/DspTap/actions/runs/35861317022) (`21b3488`; MinSizeRel, `size -A` `.text`) | 13.2.1 (15:13.2.rel1-2) |
| m55-ooura | `rfft_q31_512` | 26,593 | 27,392 | **recorded at Stage 2c**: measured + 3 %, rounded up to 64 bytes | [run 35861317022](https://github.com/tap/DspTap/actions/runs/35861317022) (`21b3488`; MinSizeRel, `size -A` `.text`) | 13.2.1 (15:13.2.rel1-2) |

### Seeding, and how the job decides what to do

The workflow decides per key from `baselines.json`:

- **compare** — the key has baselines. The normal mode: ±3 %, two-sided,
  and a recorded scenario with no binary fails too (`STALE BASELINE`).
- **seed** — the key is empty **and** the event is a push to `main` or a
  `workflow_dispatch`. The job runs `scripts/icount.py --update`, prints
  every scenario's count into the step summary and uploads
  `baselines-<key>`; the `seed-summary` job (which fails if any key failed,
  and is never a required check) merges the per-key files into one
  `baselines-merged` artifact.
- **refuse** — the key is empty on a `pull_request`. The job **fails** with
  the seeding instructions. Seeding from a PR would be a gate bypass (empty
  the key, go green), so it is never allowed; likewise a PR whose
  `baselines.json` empties a key that its base branch had seeded fails in
  the "Baselines state" step, before anything is counted.

Every run, in every mode, uploads `measured-<key>` — the counts this run
measured, as JSON — so a scenario added later (Stage 3b's `rfft_q15_512`,
`rfft_q31_*`) gets its number by "run CI, copy the artifact, commit", never
by typing one in.

The seeding commit copies `baselines-merged` to `bench/baselines.json`,
adds the rows above with the `main` run's URL and SHA and the versions the
job printed in its "Toolchain versions" step, sets the `text_ceiling_*`
numbers in `bench.yml` from the same run's size step (policy under "Sizes"
above), and nothing else rides in it. A scenario added to a key that is
already seeded (the Stage 2c fixed-point scenarios were the first) never
sees seed mode: its count is read from a `workflow_dispatch` run's log or
`measured-<key>` artifact — compare mode prints every new scenario as
`NO BASELINE` and records it — and committed the same way, with its row.
From that commit on, the key's job compares and a red ratchet is a
failing check.

The first push to `main` after the QEMU legs (#17) and this scaffold have
merged is therefore the seeding run for all five keys (wave 2 of Part 12);
re-recording later is a `workflow_dispatch` on `main`, or the hand
procedure below.

Re-recording by hand, in the target's environment (cross toolchain, QEMU and
the plugin on PATH):

```sh
gcc -shared -fPIC $(pkg-config --cflags glib-2.0) -I<dir with qemu-plugin.h> \
    -o /tmp/libinsncount.so tools/qemu_insn_plugin/insn_count.c
cmake -S . -B build-m33 -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_TOOLCHAIN_FILE=cmake/arm-cortex-m33-mps2.cmake \
    -DTAP_DSP_BUILD_TESTS=OFF -DTAP_DSP_BUILD_BENCH=ON
cmake --build build-m33 -j
python3 scripts/icount.py --target m33 --build-dir build-m33 \
    --plugin /tmp/libinsncount.so [--update] [--record measured-m33.json]
python3 scripts/icount.py --merge a.json b.json    # fold per-key files into one
```

## Policy

- **Tolerance ±3 %**, the family's number. Tighter is not better: link-order
  and libgcc changes move counts by fractions of a percent legitimately.
- **Two-sided.** A regression beyond tolerance fails; an *improvement* beyond
  tolerance also fails, so a stale, too-high baseline can never let a later
  regression hide in the slack — the winning commit re-records.
- **The ratchet runs on every pull request and on every push to `main`**
  on every QEMU leg (one run per ref at a time). A red ratchet is a failing
  check, not a warning. Once seeded, the five `icount <key>` jobs are the
  required checks; the artifact-merge job never is.
- **`--update` is a written commit on its own**: the measured before/after per
  key and the reason go into the table above, in its fixed shape. An expected
  regression (a correctness fix that costs instructions) is written down;
  there is no silent absorb. A key is never emptied to "reset" it — the
  workflow fails a pull request that does.
- **Size runs in the same job, MinSizeRel.** The job builds the three
  shipping size probes — `tap_dsp_size_probe_rfft_f32_512`,
  `tap_dsp_size_probe_rfft_q15_512`, `tap_dsp_size_probe_rfft_q31_512`
  (`bench/size_probe.cpp`: one profile, no stdio) — with
  `-DCMAKE_BUILD_TYPE=MinSizeRel` — the build type the test legs and Part 10
  item 4 use, so a ceiling here and a size in `ci.yml` describe the same
  object — and reads the `.text` row of `arm-none-eabi-size -A` (Berkeley
  format folds `.rodata` and NOLOAD sections into "text"). Each key carries
  one ceiling per probe in `bench.yml` (`text_ceiling_f32`, the Part 4
  assertion for the float instantiation — CMSIS-DSP on the `m55` key;
  `text_ceiling_q15`; `text_ceiling_q31`); `0` means not yet recorded and
  the step only prints. **Ceiling policy:** a ceiling is the `.text` figure
  measured on the run that recorded it, plus 3 % (the family's number,
  one-sided here: only exceeding it fails), rounded up to the next multiple
  of 64 bytes; it is re-recorded the way a count is — a written commit with
  a row in the sizes table above (before, after, reason, run) — when a
  change legitimately moves the figure, and never widened to absorb one. The
  step runs with `pipefail` and fails on a `size` that fails or prints no
  numeric `.text` row, for any probe: a measurement that cannot be taken is
  an error, not an empty figure.
- **Wall clock is never a gate.** `bench_fft` is the local tool for the
  desktop and Apple vDSP claims; its numbers go into `docs/fft-design.md`
  with machine and date.

## What the ratchet catches, by stage

| stage | what a red ratchet means |
|---|---|
| 2a (port beside the C, nothing routed) | nothing: the `_port` scenarios were informational, and the bare keys measured the unchanged C at 0 % |
| 2b (routing flip to the port; done) | a port that vectorizes worse than the C; the C-vs-port ratio printed from 2a on was the gate, and the bare keys were re-recorded to the port's exact counts (table above) |
| 2b → 2c (the C as the `_c` informational sibling; done) | nothing gated; the job log showed C/shipping and whether the checksums still agreed, until 2c retired the C |
| 2c / 3b (fixed point; the Q15/Q31 scenarios seeded at 2c) | the seed itself gates nothing; every later `SMMULR`, Helium or table-layout change to the Q15/Q31 kernel has a number to beat |
| 4 (engine as a parameter) | a backend regression on the deployed `m55` profile, or the `m55-ooura` fallback quietly getting slower |
| 6 (hygiene) | the hygiene pass slowing the hot path unnoticed |

## Provenance

`icount/icount_main.cpp`, `scripts/icount.py` and
`tools/qemu_insn_plugin/insn_count.c` are adapted from MuTap's (MIT, MuTap
contributors; the same pattern lives in SampleRateTap and RatioTap). Copied
rather than shared because they are small; a taphouse-style consolidation is
the eventual home (Part 11, "Sharing"). Everything under `bench/` and
`scripts/` is DspTap-authored MIT.

The plugin is compiled against QEMU's `qemu-plugin.h`, which is
`SPDX-License-Identifier: GPL-2.0-or-later`. That header is fetched at CI
time (digest-verified), never vendored into this repo, and used only to
build a test tool that runs in CI and is not shipped; nothing in the DspTap
tree or in what consumers link is GPL, so `NOTICE.md` carries no entry for
it.
