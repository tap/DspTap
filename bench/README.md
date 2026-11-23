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
| `rfft_f32_512` | float | 512 | `basic_real_fft<float>` as built: the srdif engine (`fft/srdif.h`) since tap/DspTap#42; the port of Ooura's `rdft` from Stage 2b until #42, the vendored Ooura C before that; CMSIS-DSP Helium behind it on `m55` | all legs | yes |
| `rfft_f32_2048` | float | 2048 | as above | all legs | yes |
| `rfft_f64_512` | double | 512 | `basic_real_fft<double>`: the srdif engine | host-class only: soft-float double is not a profile, so the bare-metal legs neither build nor baseline it | no key counts it |
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
the shipping tree — a reference copy under `tests/reference/ooura/` served
the bit-identity gate alone until Decision D6 deleted it — and retired the
selector (`TAP_DSP_BENCH_ENGINE`), the adapter and the `_c` pairs. The `DONE` line
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
TAP_DSP_ICOUNT_DONE ok=1 engine=basic_real_fft backend=srdif scenario=rfft_f32_512 checksum=0x322a64b478d14325
```

(the `m4f` key's line on the srdif engine since its Hexagon tuning;
`0xd26ebf9b45534325` from #42 until then: the same outputs at
`-ffp-contract=off`, fused differently by GCC). The fold is exact and
order-sensitive: two runs of the same binary print the same value, and a
1-ulp change in any single output changes it (verified on the port, the
engine of the time, by nudging one spectrum bin at one iteration with
`nextafterf`: its x86-64 value `0x662dd085b5b88325` became
`0x259510dc8721cba8`). A floating running sum cannot promise that — it
absorbs differences below the accumulator's ulp — which is why the checksum
is an integer hash: until Stage 2c it was the fingerprint the job compared between
the C and the port on each QEMU leg at the bench's own build flags (each key
its own value; the port's `m4f` value was `0xa686b767a9d68725`), beside the
parity TU's comparison at `-ffp-contract=off` (retired at D6; the port was
bit-identical to the C, so `0x662dd085b5b88325` was also the value the C
printed on x86-64 before the flip), and it is what a same-key checksum
change means today: an output bit moved. The srdif engine prints
`0x5caf505901a97f25` on x86-64 and the `m4f` line above; neither matches
the port's. (The port's x86-64 value stayed in this paragraph through #42's
clean-room hand-off: `docs/fft-design.md`, "The floating engine, replaced
clean-room", item 1.) `ok` is a sanity check that the last iteration
round-trips its input — within 1e-3 in the floating scenarios' own
precision; within eight output LSB referred to the input, i.e.
|out · 2^s − x| ≤ 8 · 2^s with s = e_fwd + e_inv + 1 − log₂N, in the
fixed-point ones (the contract's round-trip identity, `fft.h`; eight is the
header's pinned Q31 single-transform worst case rounded to a power of two,
against 0.52 / 3.6 LSB measured in this workload for Q15 / Q31); `backend=`
names what the build routed `basic_real_fft` through: `srdif`,
`cmsis` on the `m55` key, `accelerate` on Apple, `fixed_point` for Q15 /
Q31.

Why these numbers: on the host (x86-64, GCC 13 `-O2`, callgrind) the three
floating scenarios execute 109 M, 125 M and 108 M instructions (the three
fixed-point ones 434 M, 407 M and 495 M at `-O3`, of which the kernel's two
radix-4 stage functions are 71-80 %; construction is under 0.5 %).
Construction plus the one-time table build plus the print is under
0.2 M — well under the 1 % the design asks for (for the srdif engine, whose
tables are integer arithmetic, construction measures 0.10 M / 0.40 M
instructions at N = 512 / 2048 on the Cortex-M4F: 0.11 % / 0.37 %). The share of the count that
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

### The Butterworth stage's scenarios (bracket sums)

Built from `icount/stage_main.cpp`, each measuring
`basic_butterworth_lowpass<float, order, channels>` (`butterworth.h`) at
fc = 0.1 over 16,384 interleaved frames:

| scenario | order | channels | frames per `process()` | what it is there for | targets | gated |
|---|---|---|---|---|---|---|
| `bw4_f32_ch2` | 4 | 2 | 128 | the stage as an oversampling kernel runs it | all legs | yes |
| `bw8_f32_ch2` | 8 | 2 | 128 | the order-8 cascade | all legs | yes |
| `bw8_f32_ch1` | 8 | 1 | 128 | a mono stage (what a planar host runs per channel) | all legs | yes |
| `bw4_f32_ch2_n4` | 4 | 2 | 4 | short calls: what a stage after a decimator in a `chain` is fed | all legs | yes |
| `bw8_f32_set_cutoff` | 8 | 2 | n/a | four `set_cutoff()` calls at 0.45 / os, os = 2, 4, 8, 16 | all legs | where double is hardware: `hexagon`, `m55`, `m55-ooura` |
| `bw8_f32_set_cutoff_deviation` | 8 | 2 | n/a | the same four calls, each followed by one read of `latency_deviation_output_frames()` | all legs | as above |

**These counts are bracket sums, not whole runs.** Each step (one
`process()` call, or one `set_cutoff()`) runs between `bench/call_marker.h`'s
two markers, and the plugin sums the instructions between them
(`TAP_DSP_CALL_MAX ... total=`). The workload's DONE line states
`brackets=<n>`; `scripts/icount.py` then takes the bracket sum as the
scenario's count and prints `brackets=<n>` beside it, and it aborts the run
when the plugin's bracket line is missing or counts another number of
brackets, so a bracketed scenario can never fall back to its whole-run count.
The whole run would not do as the gate: the stage's design at construction is
double, soft float on the M4 and the M33, where it is about 17 M instructions
against 1.8 M of streaming for `bw4_f32_ch2`, so a quarter of the kernel's
cost is 2.5 % of the whole run, inside the tolerance. On the `hexagon` key
the bracket sum also leaves out qemu-hexagon's per-process offset (below): it
is the same number for any environment and path.

The two `set_cutoff()` scenarios run the design, which is double. On the keys
where double is libgcc's soft float (`m4-softfp`, `m4f`, `m33`) they would
gate the toolchain's library, the reason `rfft_f64_512` is not built there;
there they are built as `tap_dsp_info_*`, outside the ratchet's glob, and the
job's informational step prints their counts. `bench/icount/CMakeLists.txt`
decides by compiling a check of `__ARM_FP`.

## Baselines and targets

[`baselines.json`](baselines.json) is keyed by target, one dict per key:

| key | core | QEMU machine | float32 FFT |
|---|---|---|---|
| `m4-softfp` | Cortex-M4, `-mfloat-abi=soft` | `mps2-an386` | srdif engine (every float op a library call) |
| `m4f` | Cortex-M4F, `fpv4-sp-d16` | `mps2-an386` | srdif engine |
| `m33` | Cortex-M33, no MVE | `mps2-an505` | srdif engine |
| `m55` | Cortex-M55, Helium | `mps3-an547` | CMSIS-DSP — the deployed profile |
| `m55-ooura` | Cortex-M55, `-DTAP_DSP_FFT_CMSIS=OFF` | `mps3-an547` | srdif engine — the fallback (the key keeps the historical name of the vendored C the ratchet was seeded on) |
| `hexagon` | Hexagon v68 + HVX-128, CodeLinaro clang 19.1.5, static musl | `qemu-hexagon` user mode (8.2.2, built with plugins) | srdif engine; hosted, so `rfft_f64_512` is built and gated too |

The `hexagon` key counts **packets**, not instructions: qemu-hexagon reports
one guest instruction per VLIW packet (up to four instructions), as the
per-address profile of the Hexagon tuning pass showed (`docs/fft-design.md`,
"Hexagon (clang) tuning"). So on that key the ratchet rewards what the
compiler packs side by side as well as what it executes. The key was added
when MuTap's Hexagon ratchet found the srdif engine 12–13 % above the port on
its chain workloads, which no key here had measured. Its counts are
comparable with MuTap's Hexagon ratchet: same toolchain, flags and QEMU.
Absolute counts carry a per-process offset: user-mode qemu hands the guest
its environment and the binary's path, and the count moves with them. A
44-character longer path read +731 packets, and `env -i` read −1,584. The
offset is constant across the scenarios of one setup. Local setups have read
+473 to +1,723 against the runner, and CI runs agree with each other to the
packet. Compare within a run, and take baselines from CI only (as for every
key). The predecessor cells in the `hexagon` rows below are
local measurements converted to the runner's offset, ±28 packets: a rebuild
of `7a58ebe` at another path read its float scenarios 28 below them and its
fixed-point scenarios equal to #42's.

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
| m4-softfp | rfft_q15_512 | — | 746,311,273 | — | **seeded at Stage 2c** (the Q15/Q31 scenarios Stage 3b specified, Part 11; `basic_real_fft<std::int16_t \| std::int32_t>` under `scaling::fixed`, the int32 kernel on every key: fixed point has no backend, so the two M55 keys measure the same binary) | [run 35861317022](https://github.com/tap/DspTap/actions/runs/35861317022) (workflow_dispatch on the PR branch at `21b3488`; the PR's compare-mode run confirms +0.00 %) | `8350f13` (the #32 squash on `main`, filled in after tap/MuTap#55 pinned it, per #32's reviewer notes) | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m4-softfp | rfft_q31_512 | — | 709,161,574 | — | **seeded at Stage 2c** (the Q15/Q31 scenarios Stage 3b specified, Part 11; `basic_real_fft<std::int16_t \| std::int32_t>` under `scaling::fixed`, the int32 kernel on every key: fixed point has no backend, so the two M55 keys measure the same binary) | [run 35861317022](https://github.com/tap/DspTap/actions/runs/35861317022) (workflow_dispatch on the PR branch at `21b3488`; the PR's compare-mode run confirms +0.00 %) | `8350f13` (the #32 squash on `main`, filled in after tap/MuTap#55 pinned it, per #32's reviewer notes) | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m4-softfp | rfft_q31_2048 | — | 869,187,566 | — | **seeded at Stage 2c** (the Q15/Q31 scenarios Stage 3b specified, Part 11; `basic_real_fft<std::int16_t \| std::int32_t>` under `scaling::fixed`, the int32 kernel on every key: fixed point has no backend, so the two M55 keys measure the same binary) | [run 35861317022](https://github.com/tap/DspTap/actions/runs/35861317022) (workflow_dispatch on the PR branch at `21b3488`; the PR's compare-mode run confirms +0.00 %) | `8350f13` (the #32 squash on `main`, filled in after tap/MuTap#55 pinned it, per #32's reviewer notes) | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m4f | rfft_q15_512 | — | 753,215,637 | — | **seeded at Stage 2c** (the Q15/Q31 scenarios Stage 3b specified, Part 11; `basic_real_fft<std::int16_t \| std::int32_t>` under `scaling::fixed`, the int32 kernel on every key: fixed point has no backend, so the two M55 keys measure the same binary) | [run 35861317022](https://github.com/tap/DspTap/actions/runs/35861317022) (workflow_dispatch on the PR branch at `21b3488`; the PR's compare-mode run confirms +0.00 %) | `8350f13` (the #32 squash on `main`, filled in after tap/MuTap#55 pinned it, per #32's reviewer notes) | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m4f | rfft_q31_512 | — | 715,019,447 | — | **seeded at Stage 2c** (the Q15/Q31 scenarios Stage 3b specified, Part 11; `basic_real_fft<std::int16_t \| std::int32_t>` under `scaling::fixed`, the int32 kernel on every key: fixed point has no backend, so the two M55 keys measure the same binary) | [run 35861317022](https://github.com/tap/DspTap/actions/runs/35861317022) (workflow_dispatch on the PR branch at `21b3488`; the PR's compare-mode run confirms +0.00 %) | `8350f13` (the #32 squash on `main`, filled in after tap/MuTap#55 pinned it, per #32's reviewer notes) | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m4f | rfft_q31_2048 | — | 876,486,909 | — | **seeded at Stage 2c** (the Q15/Q31 scenarios Stage 3b specified, Part 11; `basic_real_fft<std::int16_t \| std::int32_t>` under `scaling::fixed`, the int32 kernel on every key: fixed point has no backend, so the two M55 keys measure the same binary) | [run 35861317022](https://github.com/tap/DspTap/actions/runs/35861317022) (workflow_dispatch on the PR branch at `21b3488`; the PR's compare-mode run confirms +0.00 %) | `8350f13` (the #32 squash on `main`, filled in after tap/MuTap#55 pinned it, per #32's reviewer notes) | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m33 | rfft_q15_512 | — | 752,189,619 | — | **seeded at Stage 2c** (the Q15/Q31 scenarios Stage 3b specified, Part 11; `basic_real_fft<std::int16_t \| std::int32_t>` under `scaling::fixed`, the int32 kernel on every key: fixed point has no backend, so the two M55 keys measure the same binary) | [run 35861317022](https://github.com/tap/DspTap/actions/runs/35861317022) (workflow_dispatch on the PR branch at `21b3488`; the PR's compare-mode run confirms +0.00 %) | `8350f13` (the #32 squash on `main`, filled in after tap/MuTap#55 pinned it, per #32's reviewer notes) | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m33 | rfft_q31_512 | — | 714,626,257 | — | **seeded at Stage 2c** (the Q15/Q31 scenarios Stage 3b specified, Part 11; `basic_real_fft<std::int16_t \| std::int32_t>` under `scaling::fixed`, the int32 kernel on every key: fixed point has no backend, so the two M55 keys measure the same binary) | [run 35861317022](https://github.com/tap/DspTap/actions/runs/35861317022) (workflow_dispatch on the PR branch at `21b3488`; the PR's compare-mode run confirms +0.00 %) | `8350f13` (the #32 squash on `main`, filled in after tap/MuTap#55 pinned it, per #32's reviewer notes) | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m33 | rfft_q31_2048 | — | 875,801,377 | — | **seeded at Stage 2c** (the Q15/Q31 scenarios Stage 3b specified, Part 11; `basic_real_fft<std::int16_t \| std::int32_t>` under `scaling::fixed`, the int32 kernel on every key: fixed point has no backend, so the two M55 keys measure the same binary) | [run 35861317022](https://github.com/tap/DspTap/actions/runs/35861317022) (workflow_dispatch on the PR branch at `21b3488`; the PR's compare-mode run confirms +0.00 %) | `8350f13` (the #32 squash on `main`, filled in after tap/MuTap#55 pinned it, per #32's reviewer notes) | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m55 | rfft_q15_512 | — | 684,157,511 | — | **seeded at Stage 2c** (the Q15/Q31 scenarios Stage 3b specified, Part 11; `basic_real_fft<std::int16_t \| std::int32_t>` under `scaling::fixed`, the int32 kernel on every key: fixed point has no backend, so the two M55 keys measure the same binary) | [run 35861317022](https://github.com/tap/DspTap/actions/runs/35861317022) (workflow_dispatch on the PR branch at `21b3488`; the PR's compare-mode run confirms +0.00 %) | `8350f13` (the #32 squash on `main`, filled in after tap/MuTap#55 pinned it, per #32's reviewer notes) | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m55 | rfft_q31_512 | — | 657,595,969 | — | **seeded at Stage 2c** (the Q15/Q31 scenarios Stage 3b specified, Part 11; `basic_real_fft<std::int16_t \| std::int32_t>` under `scaling::fixed`, the int32 kernel on every key: fixed point has no backend, so the two M55 keys measure the same binary) | [run 35861317022](https://github.com/tap/DspTap/actions/runs/35861317022) (workflow_dispatch on the PR branch at `21b3488`; the PR's compare-mode run confirms +0.00 %) | `8350f13` (the #32 squash on `main`, filled in after tap/MuTap#55 pinned it, per #32's reviewer notes) | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m55 | rfft_q31_2048 | — | 806,142,635 | — | **seeded at Stage 2c** (the Q15/Q31 scenarios Stage 3b specified, Part 11; `basic_real_fft<std::int16_t \| std::int32_t>` under `scaling::fixed`, the int32 kernel on every key: fixed point has no backend, so the two M55 keys measure the same binary) | [run 35861317022](https://github.com/tap/DspTap/actions/runs/35861317022) (workflow_dispatch on the PR branch at `21b3488`; the PR's compare-mode run confirms +0.00 %) | `8350f13` (the #32 squash on `main`, filled in after tap/MuTap#55 pinned it, per #32's reviewer notes) | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m55-ooura | rfft_q15_512 | — | 684,157,511 | — | **seeded at Stage 2c** (the Q15/Q31 scenarios Stage 3b specified, Part 11; `basic_real_fft<std::int16_t \| std::int32_t>` under `scaling::fixed`, the int32 kernel on every key: fixed point has no backend, so the two M55 keys measure the same binary) | [run 35861317022](https://github.com/tap/DspTap/actions/runs/35861317022) (workflow_dispatch on the PR branch at `21b3488`; the PR's compare-mode run confirms +0.00 %) | `8350f13` (the #32 squash on `main`, filled in after tap/MuTap#55 pinned it, per #32's reviewer notes) | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m55-ooura | rfft_q31_512 | — | 657,595,969 | — | **seeded at Stage 2c** (the Q15/Q31 scenarios Stage 3b specified, Part 11; `basic_real_fft<std::int16_t \| std::int32_t>` under `scaling::fixed`, the int32 kernel on every key: fixed point has no backend, so the two M55 keys measure the same binary) | [run 35861317022](https://github.com/tap/DspTap/actions/runs/35861317022) (workflow_dispatch on the PR branch at `21b3488`; the PR's compare-mode run confirms +0.00 %) | `8350f13` (the #32 squash on `main`, filled in after tap/MuTap#55 pinned it, per #32's reviewer notes) | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m55-ooura | rfft_q31_2048 | — | 806,142,635 | — | **seeded at Stage 2c** (the Q15/Q31 scenarios Stage 3b specified, Part 11; `basic_real_fft<std::int16_t \| std::int32_t>` under `scaling::fixed`, the int32 kernel on every key: fixed point has no backend, so the two M55 keys measure the same binary) | [run 35861317022](https://github.com/tap/DspTap/actions/runs/35861317022) (workflow_dispatch on the PR branch at `21b3488`; the PR's compare-mode run confirms +0.00 %) | `8350f13` (the #32 squash on `main`, filled in after tap/MuTap#55 pinned it, per #32's reviewer notes) | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m4-softfp | rfft_f32_2048 | 2,294,360,259 | 2,243,212,021 | −2.23 % | **re-recorded at the srdif replacement** (tap/DspTap#42): the floating engine behind `basic_real_fft<float>` on these keys is now `fft/srdif.h` (a split-radix DIF kernel of length N/2 written from the literature), which replaced the port of Ooura's `rdft` (`fft/split_radix.h`, deleted at #42); every float output bit changes (new checksums) and the counts fall on every key, beyond the −3 % band on every key but m4-softfp, so the exact new counts are recorded (Part 11 / D11). The counts are those of the fix pass for review A of #42 (the bit-reversal permutation without an index table, GCC's SLP vectorizer off for the engine; output bits unchanged), 0.04 … 2.62 % below the first srdif record (1,815,063,718 / 2,244,317,934 m4-softfp, 94,003,130 / 108,030,219 m4f, 94,382,323 / 108,470,414 m33, 86,188,216 / 100,091,803 m55-ooura at 512 / 2048). Fixed-point scenarios and the `m55` (CMSIS) key are not re-recorded: measured inside the band (fixed point +0.49 … +1.15 %, CMSIS +0.00 %) with no change to that code. `docs/fft-design.md`, "The floating engine (srdif)" | measured locally with `scripts/icount.py` exactly as `bench.yml` runs it (fresh Release build per key, TAP_DSP_BUILD_TESTS=OFF, the pinned plugin header), 2026-09-27; confirmed to the instruction by the push-to-`main` [run 36348787347](https://github.com/tap/DspTap/actions/runs/36348787347) (compare mode, +0.00 %) | `72977aa` (the #42 squash) | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m4-softfp | rfft_f32_512 | 1,864,929,141 | 1,814,303,702 | −2.71 % | **re-recorded at the srdif replacement** (tap/DspTap#42): the floating engine behind `basic_real_fft<float>` on these keys is now `fft/srdif.h` (a split-radix DIF kernel of length N/2 written from the literature), which replaced the port of Ooura's `rdft` (`fft/split_radix.h`, deleted at #42); every float output bit changes (new checksums) and the counts fall on every key, beyond the −3 % band on every key but m4-softfp, so the exact new counts are recorded (Part 11 / D11). The counts are those of the fix pass for review A of #42 (the bit-reversal permutation without an index table, GCC's SLP vectorizer off for the engine; output bits unchanged), 0.04 … 2.62 % below the first srdif record (1,815,063,718 / 2,244,317,934 m4-softfp, 94,003,130 / 108,030,219 m4f, 94,382,323 / 108,470,414 m33, 86,188,216 / 100,091,803 m55-ooura at 512 / 2048). Fixed-point scenarios and the `m55` (CMSIS) key are not re-recorded: measured inside the band (fixed point +0.49 … +1.15 %, CMSIS +0.00 %) with no change to that code. `docs/fft-design.md`, "The floating engine (srdif)" | measured locally with `scripts/icount.py` exactly as `bench.yml` runs it (fresh Release build per key, TAP_DSP_BUILD_TESTS=OFF, the pinned plugin header), 2026-09-27; confirmed to the instruction by the push-to-`main` [run 36348787347](https://github.com/tap/DspTap/actions/runs/36348787347) (compare mode, +0.00 %) | `72977aa` (the #42 squash) | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m4f | rfft_f32_2048 | 111,278,416 | 106,282,752 | −4.49 % | **re-recorded at the srdif replacement** (tap/DspTap#42): the floating engine behind `basic_real_fft<float>` on these keys is now `fft/srdif.h` (a split-radix DIF kernel of length N/2 written from the literature), which replaced the port of Ooura's `rdft` (`fft/split_radix.h`, deleted at #42); every float output bit changes (new checksums) and the counts fall on every key, beyond the −3 % band on every key but m4-softfp, so the exact new counts are recorded (Part 11 / D11). The counts are those of the fix pass for review A of #42 (the bit-reversal permutation without an index table, GCC's SLP vectorizer off for the engine; output bits unchanged), 0.04 … 2.62 % below the first srdif record (1,815,063,718 / 2,244,317,934 m4-softfp, 94,003,130 / 108,030,219 m4f, 94,382,323 / 108,470,414 m33, 86,188,216 / 100,091,803 m55-ooura at 512 / 2048). Fixed-point scenarios and the `m55` (CMSIS) key are not re-recorded: measured inside the band (fixed point +0.49 … +1.15 %, CMSIS +0.00 %) with no change to that code. `docs/fft-design.md`, "The floating engine (srdif)" | measured locally with `scripts/icount.py` exactly as `bench.yml` runs it (fresh Release build per key, TAP_DSP_BUILD_TESTS=OFF, the pinned plugin header), 2026-09-27; confirmed to the instruction by the push-to-`main` [run 36348787347](https://github.com/tap/DspTap/actions/runs/36348787347) (compare mode, +0.00 %) | `72977aa` (the #42 squash) | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m4f | rfft_f32_512 | 97,138,544 | 92,575,609 | −4.70 % | **re-recorded at the srdif replacement** (tap/DspTap#42): the floating engine behind `basic_real_fft<float>` on these keys is now `fft/srdif.h` (a split-radix DIF kernel of length N/2 written from the literature), which replaced the port of Ooura's `rdft` (`fft/split_radix.h`, deleted at #42); every float output bit changes (new checksums) and the counts fall on every key, beyond the −3 % band on every key but m4-softfp, so the exact new counts are recorded (Part 11 / D11). The counts are those of the fix pass for review A of #42 (the bit-reversal permutation without an index table, GCC's SLP vectorizer off for the engine; output bits unchanged), 0.04 … 2.62 % below the first srdif record (1,815,063,718 / 2,244,317,934 m4-softfp, 94,003,130 / 108,030,219 m4f, 94,382,323 / 108,470,414 m33, 86,188,216 / 100,091,803 m55-ooura at 512 / 2048). Fixed-point scenarios and the `m55` (CMSIS) key are not re-recorded: measured inside the band (fixed point +0.49 … +1.15 %, CMSIS +0.00 %) with no change to that code. `docs/fft-design.md`, "The floating engine (srdif)" | measured locally with `scripts/icount.py` exactly as `bench.yml` runs it (fresh Release build per key, TAP_DSP_BUILD_TESTS=OFF, the pinned plugin header), 2026-09-27; confirmed to the instruction by the push-to-`main` [run 36348787347](https://github.com/tap/DspTap/actions/runs/36348787347) (compare mode, +0.00 %) | `72977aa` (the #42 squash) | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m33 | rfft_f32_2048 | 115,465,626 | 106,773,786 | −7.53 % | **re-recorded at the srdif replacement** (tap/DspTap#42): the floating engine behind `basic_real_fft<float>` on these keys is now `fft/srdif.h` (a split-radix DIF kernel of length N/2 written from the literature), which replaced the port of Ooura's `rdft` (`fft/split_radix.h`, deleted at #42); every float output bit changes (new checksums) and the counts fall on every key, beyond the −3 % band on every key but m4-softfp, so the exact new counts are recorded (Part 11 / D11). The counts are those of the fix pass for review A of #42 (the bit-reversal permutation without an index table, GCC's SLP vectorizer off for the engine; output bits unchanged), 0.04 … 2.62 % below the first srdif record (1,815,063,718 / 2,244,317,934 m4-softfp, 94,003,130 / 108,030,219 m4f, 94,382,323 / 108,470,414 m33, 86,188,216 / 100,091,803 m55-ooura at 512 / 2048). Fixed-point scenarios and the `m55` (CMSIS) key are not re-recorded: measured inside the band (fixed point +0.49 … +1.15 %, CMSIS +0.00 %) with no change to that code. `docs/fft-design.md`, "The floating engine (srdif)" | measured locally with `scripts/icount.py` exactly as `bench.yml` runs it (fresh Release build per key, TAP_DSP_BUILD_TESTS=OFF, the pinned plugin header), 2026-09-27; confirmed to the instruction by the push-to-`main` [run 36348787347](https://github.com/tap/DspTap/actions/runs/36348787347) (compare mode, +0.00 %) | `72977aa` (the #42 squash) | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m33 | rfft_f32_512 | 100,945,841 | 92,979,212 | −7.89 % | **re-recorded at the srdif replacement** (tap/DspTap#42): the floating engine behind `basic_real_fft<float>` on these keys is now `fft/srdif.h` (a split-radix DIF kernel of length N/2 written from the literature), which replaced the port of Ooura's `rdft` (`fft/split_radix.h`, deleted at #42); every float output bit changes (new checksums) and the counts fall on every key, beyond the −3 % band on every key but m4-softfp, so the exact new counts are recorded (Part 11 / D11). The counts are those of the fix pass for review A of #42 (the bit-reversal permutation without an index table, GCC's SLP vectorizer off for the engine; output bits unchanged), 0.04 … 2.62 % below the first srdif record (1,815,063,718 / 2,244,317,934 m4-softfp, 94,003,130 / 108,030,219 m4f, 94,382,323 / 108,470,414 m33, 86,188,216 / 100,091,803 m55-ooura at 512 / 2048). Fixed-point scenarios and the `m55` (CMSIS) key are not re-recorded: measured inside the band (fixed point +0.49 … +1.15 %, CMSIS +0.00 %) with no change to that code. `docs/fft-design.md`, "The floating engine (srdif)" | measured locally with `scripts/icount.py` exactly as `bench.yml` runs it (fresh Release build per key, TAP_DSP_BUILD_TESTS=OFF, the pinned plugin header), 2026-09-27; confirmed to the instruction by the push-to-`main` [run 36348787347](https://github.com/tap/DspTap/actions/runs/36348787347) (compare mode, +0.00 %) | `72977aa` (the #42 squash) | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m55-ooura | rfft_f32_2048 | 102,784,096 | 97,602,563 | −5.04 % | **re-recorded at the srdif replacement** (tap/DspTap#42): the floating engine behind `basic_real_fft<float>` on these keys is now `fft/srdif.h` (a split-radix DIF kernel of length N/2 written from the literature), which replaced the port of Ooura's `rdft` (`fft/split_radix.h`, deleted at #42); every float output bit changes (new checksums) and the counts fall on every key, beyond the −3 % band on every key but m4-softfp, so the exact new counts are recorded (Part 11 / D11). The counts are those of the fix pass for review A of #42 (the bit-reversal permutation without an index table, GCC's SLP vectorizer off for the engine; output bits unchanged), 0.04 … 2.62 % below the first srdif record (1,815,063,718 / 2,244,317,934 m4-softfp, 94,003,130 / 108,030,219 m4f, 94,382,323 / 108,470,414 m33, 86,188,216 / 100,091,803 m55-ooura at 512 / 2048). Fixed-point scenarios and the `m55` (CMSIS) key are not re-recorded: measured inside the band (fixed point +0.49 … +1.15 %, CMSIS +0.00 %) with no change to that code. `docs/fft-design.md`, "The floating engine (srdif)" | measured locally with `scripts/icount.py` exactly as `bench.yml` runs it (fresh Release build per key, TAP_DSP_BUILD_TESTS=OFF, the pinned plugin header), 2026-09-27; confirmed to the instruction by the push-to-`main` [run 36348787347](https://github.com/tap/DspTap/actions/runs/36348787347) (compare mode, +0.00 %) | `72977aa` (the #42 squash) | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m55-ooura | rfft_f32_512 | 89,276,321 | 83,933,381 | −5.98 % | **re-recorded at the srdif replacement** (tap/DspTap#42): the floating engine behind `basic_real_fft<float>` on these keys is now `fft/srdif.h` (a split-radix DIF kernel of length N/2 written from the literature), which replaced the port of Ooura's `rdft` (`fft/split_radix.h`, deleted at #42); every float output bit changes (new checksums) and the counts fall on every key, beyond the −3 % band on every key but m4-softfp, so the exact new counts are recorded (Part 11 / D11). The counts are those of the fix pass for review A of #42 (the bit-reversal permutation without an index table, GCC's SLP vectorizer off for the engine; output bits unchanged), 0.04 … 2.62 % below the first srdif record (1,815,063,718 / 2,244,317,934 m4-softfp, 94,003,130 / 108,030,219 m4f, 94,382,323 / 108,470,414 m33, 86,188,216 / 100,091,803 m55-ooura at 512 / 2048). Fixed-point scenarios and the `m55` (CMSIS) key are not re-recorded: measured inside the band (fixed point +0.49 … +1.15 %, CMSIS +0.00 %) with no change to that code. `docs/fft-design.md`, "The floating engine (srdif)" | measured locally with `scripts/icount.py` exactly as `bench.yml` runs it (fresh Release build per key, TAP_DSP_BUILD_TESTS=OFF, the pinned plugin header), 2026-09-27; confirmed to the instruction by the push-to-`main` [run 36348787347](https://github.com/tap/DspTap/actions/runs/36348787347) (compare mode, +0.00 %) | `72977aa` (the #42 squash) | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m4-softfp | rfft_f32_2048 | 2,243,212,021 | 2,241,935,605 | −0.06 % | **re-recorded at the Hexagon tuning of the srdif engine**: `fft/srdif.h` was tuned for the `hexagon` key (clang 19; forced inlining under clang, the post-pass's two bin pairs loaded before either is stored, the small blocks' split-radix levels in memory down to pairs instead of register leaves, one run-time group per loop step in double, a multiply-subtract spelling for clang's contraction), and the float count falls on every key that runs the engine, −0.06 … −1.11 %, inside the band; the exact new counts are recorded (bench/README.md policy: the winning commit re-records). No output bit moves at `-ffp-contract=off` (the srdif fingerprints pass unchanged on the host, the four QEMU legs and Hexagon); at the bench's default flags the FPU keys' checksums move (GCC contracts across the former leaves differently), `m4-softfp`'s do not. Fixed-point scenarios and the `m55` (CMSIS) key are not re-recorded: unchanged to the instruction. `docs/fft-design.md`, "Hexagon (clang) tuning" | measured locally with `scripts/icount.py` exactly as `bench.yml` runs it (fresh Release build per key, TAP_DSP_BUILD_TESTS=OFF, the pinned plugin header), 2026-09-28; confirmed to the instruction by the `workflow_dispatch` [run 36363710480](https://github.com/tap/DspTap/actions/runs/36363710480) on the branch (`2a5fc69`, compare mode, +0.00 %) ; confirmed again by the push-to-`main` [run 36367784465](https://github.com/tap/DspTap/actions/runs/36367784465) (compare mode, +0.00 %) | `d9c1e33` (the #44 squash) | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m4-softfp | rfft_f32_512 | 1,814,303,702 | 1,813,126,102 | −0.06 % | **re-recorded at the Hexagon tuning of the srdif engine**: `fft/srdif.h` was tuned for the `hexagon` key (clang 19; forced inlining under clang, the post-pass's two bin pairs loaded before either is stored, the small blocks' split-radix levels in memory down to pairs instead of register leaves, one run-time group per loop step in double, a multiply-subtract spelling for clang's contraction), and the float count falls on every key that runs the engine, −0.06 … −1.11 %, inside the band; the exact new counts are recorded (bench/README.md policy: the winning commit re-records). No output bit moves at `-ffp-contract=off` (the srdif fingerprints pass unchanged on the host, the four QEMU legs and Hexagon); at the bench's default flags the FPU keys' checksums move (GCC contracts across the former leaves differently), `m4-softfp`'s do not. Fixed-point scenarios and the `m55` (CMSIS) key are not re-recorded: unchanged to the instruction. `docs/fft-design.md`, "Hexagon (clang) tuning" | measured locally with `scripts/icount.py` exactly as `bench.yml` runs it (fresh Release build per key, TAP_DSP_BUILD_TESTS=OFF, the pinned plugin header), 2026-09-28; confirmed to the instruction by the `workflow_dispatch` [run 36363710480](https://github.com/tap/DspTap/actions/runs/36363710480) on the branch (`2a5fc69`, compare mode, +0.00 %) ; confirmed again by the push-to-`main` [run 36367784465](https://github.com/tap/DspTap/actions/runs/36367784465) (compare mode, +0.00 %) | `d9c1e33` (the #44 squash) | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m4f | rfft_f32_2048 | 106,282,752 | 105,151,232 | −1.06 % | **re-recorded at the Hexagon tuning of the srdif engine**: `fft/srdif.h` was tuned for the `hexagon` key (clang 19; forced inlining under clang, the post-pass's two bin pairs loaded before either is stored, the small blocks' split-radix levels in memory down to pairs instead of register leaves, one run-time group per loop step in double, a multiply-subtract spelling for clang's contraction), and the float count falls on every key that runs the engine, −0.06 … −1.11 %, inside the band; the exact new counts are recorded (bench/README.md policy: the winning commit re-records). No output bit moves at `-ffp-contract=off` (the srdif fingerprints pass unchanged on the host, the four QEMU legs and Hexagon); at the bench's default flags the FPU keys' checksums move (GCC contracts across the former leaves differently), `m4-softfp`'s do not. Fixed-point scenarios and the `m55` (CMSIS) key are not re-recorded: unchanged to the instruction. `docs/fft-design.md`, "Hexagon (clang) tuning" | measured locally with `scripts/icount.py` exactly as `bench.yml` runs it (fresh Release build per key, TAP_DSP_BUILD_TESTS=OFF, the pinned plugin header), 2026-09-28; confirmed to the instruction by the `workflow_dispatch` [run 36363710480](https://github.com/tap/DspTap/actions/runs/36363710480) on the branch (`2a5fc69`, compare mode, +0.00 %) ; confirmed again by the push-to-`main` [run 36367784465](https://github.com/tap/DspTap/actions/runs/36367784465) (compare mode, +0.00 %) | `d9c1e33` (the #44 squash) | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m4f | rfft_f32_512 | 92,575,609 | 91,545,465 | −1.11 % | **re-recorded at the Hexagon tuning of the srdif engine**: `fft/srdif.h` was tuned for the `hexagon` key (clang 19; forced inlining under clang, the post-pass's two bin pairs loaded before either is stored, the small blocks' split-radix levels in memory down to pairs instead of register leaves, one run-time group per loop step in double, a multiply-subtract spelling for clang's contraction), and the float count falls on every key that runs the engine, −0.06 … −1.11 %, inside the band; the exact new counts are recorded (bench/README.md policy: the winning commit re-records). No output bit moves at `-ffp-contract=off` (the srdif fingerprints pass unchanged on the host, the four QEMU legs and Hexagon); at the bench's default flags the FPU keys' checksums move (GCC contracts across the former leaves differently), `m4-softfp`'s do not. Fixed-point scenarios and the `m55` (CMSIS) key are not re-recorded: unchanged to the instruction. `docs/fft-design.md`, "Hexagon (clang) tuning" | measured locally with `scripts/icount.py` exactly as `bench.yml` runs it (fresh Release build per key, TAP_DSP_BUILD_TESTS=OFF, the pinned plugin header), 2026-09-28; confirmed to the instruction by the `workflow_dispatch` [run 36363710480](https://github.com/tap/DspTap/actions/runs/36363710480) on the branch (`2a5fc69`, compare mode, +0.00 %) ; confirmed again by the push-to-`main` [run 36367784465](https://github.com/tap/DspTap/actions/runs/36367784465) (compare mode, +0.00 %) | `d9c1e33` (the #44 squash) | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m33 | rfft_f32_2048 | 106,773,786 | 106,248,986 | −0.49 % | **re-recorded at the Hexagon tuning of the srdif engine**: `fft/srdif.h` was tuned for the `hexagon` key (clang 19; forced inlining under clang, the post-pass's two bin pairs loaded before either is stored, the small blocks' split-radix levels in memory down to pairs instead of register leaves, one run-time group per loop step in double, a multiply-subtract spelling for clang's contraction), and the float count falls on every key that runs the engine, −0.06 … −1.11 %, inside the band; the exact new counts are recorded (bench/README.md policy: the winning commit re-records). No output bit moves at `-ffp-contract=off` (the srdif fingerprints pass unchanged on the host, the four QEMU legs and Hexagon); at the bench's default flags the FPU keys' checksums move (GCC contracts across the former leaves differently), `m4-softfp`'s do not. Fixed-point scenarios and the `m55` (CMSIS) key are not re-recorded: unchanged to the instruction. `docs/fft-design.md`, "Hexagon (clang) tuning" | measured locally with `scripts/icount.py` exactly as `bench.yml` runs it (fresh Release build per key, TAP_DSP_BUILD_TESTS=OFF, the pinned plugin header), 2026-09-28; confirmed to the instruction by the `workflow_dispatch` [run 36363710480](https://github.com/tap/DspTap/actions/runs/36363710480) on the branch (`2a5fc69`, compare mode, +0.00 %) ; confirmed again by the push-to-`main` [run 36367784465](https://github.com/tap/DspTap/actions/runs/36367784465) (compare mode, +0.00 %) | `d9c1e33` (the #44 squash) | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m33 | rfft_f32_512 | 92,979,212 | 92,624,908 | −0.38 % | **re-recorded at the Hexagon tuning of the srdif engine**: `fft/srdif.h` was tuned for the `hexagon` key (clang 19; forced inlining under clang, the post-pass's two bin pairs loaded before either is stored, the small blocks' split-radix levels in memory down to pairs instead of register leaves, one run-time group per loop step in double, a multiply-subtract spelling for clang's contraction), and the float count falls on every key that runs the engine, −0.06 … −1.11 %, inside the band; the exact new counts are recorded (bench/README.md policy: the winning commit re-records). No output bit moves at `-ffp-contract=off` (the srdif fingerprints pass unchanged on the host, the four QEMU legs and Hexagon); at the bench's default flags the FPU keys' checksums move (GCC contracts across the former leaves differently), `m4-softfp`'s do not. Fixed-point scenarios and the `m55` (CMSIS) key are not re-recorded: unchanged to the instruction. `docs/fft-design.md`, "Hexagon (clang) tuning" | measured locally with `scripts/icount.py` exactly as `bench.yml` runs it (fresh Release build per key, TAP_DSP_BUILD_TESTS=OFF, the pinned plugin header), 2026-09-28; confirmed to the instruction by the `workflow_dispatch` [run 36363710480](https://github.com/tap/DspTap/actions/runs/36363710480) on the branch (`2a5fc69`, compare mode, +0.00 %) ; confirmed again by the push-to-`main` [run 36367784465](https://github.com/tap/DspTap/actions/runs/36367784465) (compare mode, +0.00 %) | `d9c1e33` (the #44 squash) | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m55-ooura | rfft_f32_2048 | 97,602,563 | 97,109,507 | −0.51 % | **re-recorded at the Hexagon tuning of the srdif engine**: `fft/srdif.h` was tuned for the `hexagon` key (clang 19; forced inlining under clang, the post-pass's two bin pairs loaded before either is stored, the small blocks' split-radix levels in memory down to pairs instead of register leaves, one run-time group per loop step in double, a multiply-subtract spelling for clang's contraction), and the float count falls on every key that runs the engine, −0.06 … −1.11 %, inside the band; the exact new counts are recorded (bench/README.md policy: the winning commit re-records). No output bit moves at `-ffp-contract=off` (the srdif fingerprints pass unchanged on the host, the four QEMU legs and Hexagon); at the bench's default flags the FPU keys' checksums move (GCC contracts across the former leaves differently), `m4-softfp`'s do not. Fixed-point scenarios and the `m55` (CMSIS) key are not re-recorded: unchanged to the instruction. `docs/fft-design.md`, "Hexagon (clang) tuning" | measured locally with `scripts/icount.py` exactly as `bench.yml` runs it (fresh Release build per key, TAP_DSP_BUILD_TESTS=OFF, the pinned plugin header), 2026-09-28; confirmed to the instruction by the `workflow_dispatch` [run 36363710480](https://github.com/tap/DspTap/actions/runs/36363710480) on the branch (`2a5fc69`, compare mode, +0.00 %) ; confirmed again by the push-to-`main` [run 36367784465](https://github.com/tap/DspTap/actions/runs/36367784465) (compare mode, +0.00 %) | `d9c1e33` (the #44 squash) | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m55-ooura | rfft_f32_512 | 83,933,381 | 83,628,229 | −0.36 % | **re-recorded at the Hexagon tuning of the srdif engine**: `fft/srdif.h` was tuned for the `hexagon` key (clang 19; forced inlining under clang, the post-pass's two bin pairs loaded before either is stored, the small blocks' split-radix levels in memory down to pairs instead of register leaves, one run-time group per loop step in double, a multiply-subtract spelling for clang's contraction), and the float count falls on every key that runs the engine, −0.06 … −1.11 %, inside the band; the exact new counts are recorded (bench/README.md policy: the winning commit re-records). No output bit moves at `-ffp-contract=off` (the srdif fingerprints pass unchanged on the host, the four QEMU legs and Hexagon); at the bench's default flags the FPU keys' checksums move (GCC contracts across the former leaves differently), `m4-softfp`'s do not. Fixed-point scenarios and the `m55` (CMSIS) key are not re-recorded: unchanged to the instruction. `docs/fft-design.md`, "Hexagon (clang) tuning" | measured locally with `scripts/icount.py` exactly as `bench.yml` runs it (fresh Release build per key, TAP_DSP_BUILD_TESTS=OFF, the pinned plugin header), 2026-09-28; confirmed to the instruction by the `workflow_dispatch` [run 36363710480](https://github.com/tap/DspTap/actions/runs/36363710480) on the branch (`2a5fc69`, compare mode, +0.00 %) ; confirmed again by the push-to-`main` [run 36367784465](https://github.com/tap/DspTap/actions/runs/36367784465) (compare mode, +0.00 %) | `d9c1e33` (the #44 squash) | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| hexagon | rfft_f32_2048 | — | 51,608,903 | — | **seeded** (the `hexagon` key is new): Hexagon v68 + HVX-128, CodeLinaro clang 19.1.5, static musl, under qemu-hexagon user mode; the key counts packets, not instructions (above). Seeded on the tuned srdif engine; the engine as #42 shipped it read 65,684,295 on the first seeding run ([run 36357132506](https://github.com/tap/DspTap/actions/runs/36357132506), `28befbd`), and the engine it replaced 54,834,714 (local, the same harness, converted to the runner's offset; float ±28). `docs/fft-design.md`, "Hexagon (clang) tuning" | `workflow_dispatch` [run 36363710480](https://github.com/tap/DspTap/actions/runs/36363710480) on the branch (`2a5fc69`, seed mode); the push-to-`main` [run 36367784465](https://github.com/tap/DspTap/actions/runs/36367784465) reads +22 packets on every scenario, fixed point included (the process-environment offset, `bench/README.md` above; inside the band, not re-recorded) | `d9c1e33` (the #44 squash) | clang 19.1.5 (CodeLinaro, hexagon-unknown-linux-musl) | qemu-hexagon 8.2.2 (built from the release with `--enable-plugins`) |
| hexagon | rfft_f32_512 | — | 45,934,423 | — | **seeded** (the `hexagon` key is new): Hexagon v68 + HVX-128, CodeLinaro clang 19.1.5, static musl, under qemu-hexagon user mode; the key counts packets, not instructions (above). Seeded on the tuned srdif engine; the engine as #42 shipped it read 57,659,223 on the first seeding run ([run 36357132506](https://github.com/tap/DspTap/actions/runs/36357132506), `28befbd`), and the engine it replaced 49,732,931 (local, the same harness, converted to the runner's offset; float ±28). `docs/fft-design.md`, "Hexagon (clang) tuning" | `workflow_dispatch` [run 36363710480](https://github.com/tap/DspTap/actions/runs/36363710480) on the branch (`2a5fc69`, seed mode); the push-to-`main` [run 36367784465](https://github.com/tap/DspTap/actions/runs/36367784465) reads +22 packets on every scenario, fixed point included (the process-environment offset, `bench/README.md` above; inside the band, not re-recorded) | `d9c1e33` (the #44 squash) | clang 19.1.5 (CodeLinaro, hexagon-unknown-linux-musl) | qemu-hexagon 8.2.2 (built from the release with `--enable-plugins`) |
| hexagon | rfft_f64_512 | — | 99,939,789 | — | **seeded** (the `hexagon` key is new): Hexagon v68 + HVX-128, CodeLinaro clang 19.1.5, static musl, under qemu-hexagon user mode; the key counts packets, not instructions (above). Seeded on the tuned srdif engine; the engine as #42 shipped it read 110,683,597 on the first seeding run ([run 36357132506](https://github.com/tap/DspTap/actions/runs/36357132506), `28befbd`), and the engine it replaced 103,239,537 (local, the same harness, converted to the runner's offset; float ±28). `docs/fft-design.md`, "Hexagon (clang) tuning" | `workflow_dispatch` [run 36363710480](https://github.com/tap/DspTap/actions/runs/36363710480) on the branch (`2a5fc69`, seed mode); the push-to-`main` [run 36367784465](https://github.com/tap/DspTap/actions/runs/36367784465) reads +22 packets on every scenario, fixed point included (the process-environment offset, `bench/README.md` above; inside the band, not re-recorded) | `d9c1e33` (the #44 squash) | clang 19.1.5 (CodeLinaro, hexagon-unknown-linux-musl) | qemu-hexagon 8.2.2 (built from the release with `--enable-plugins`) |
| hexagon | rfft_q15_512 | — | 159,976,886 | — | **seeded** (the `hexagon` key is new): Hexagon v68 + HVX-128, CodeLinaro clang 19.1.5, static musl, under qemu-hexagon user mode; the key counts packets, not instructions (above). Seeded on the tuned srdif engine; the engine as #42 shipped it read 159,976,886 on the first seeding run ([run 36357132506](https://github.com/tap/DspTap/actions/runs/36357132506), `28befbd`), and the engine it replaced 159,976,886, the same code (local, the same harness, converted to the runner's offset; float ±28). `docs/fft-design.md`, "Hexagon (clang) tuning" | `workflow_dispatch` [run 36363710480](https://github.com/tap/DspTap/actions/runs/36363710480) on the branch (`2a5fc69`, seed mode); the push-to-`main` [run 36367784465](https://github.com/tap/DspTap/actions/runs/36367784465) reads +22 packets on every scenario, fixed point included (the process-environment offset, `bench/README.md` above; inside the band, not re-recorded) | `d9c1e33` (the #44 squash) | clang 19.1.5 (CodeLinaro, hexagon-unknown-linux-musl) | qemu-hexagon 8.2.2 (built from the release with `--enable-plugins`) |
| hexagon | rfft_q31_2048 | — | 184,890,581 | — | **seeded** (the `hexagon` key is new): Hexagon v68 + HVX-128, CodeLinaro clang 19.1.5, static musl, under qemu-hexagon user mode; the key counts packets, not instructions (above). Seeded on the tuned srdif engine; the engine as #42 shipped it read 184,890,581 on the first seeding run ([run 36357132506](https://github.com/tap/DspTap/actions/runs/36357132506), `28befbd`), and the engine it replaced 184,890,581, the same code (local, the same harness, converted to the runner's offset; float ±28). `docs/fft-design.md`, "Hexagon (clang) tuning" | `workflow_dispatch` [run 36363710480](https://github.com/tap/DspTap/actions/runs/36363710480) on the branch (`2a5fc69`, seed mode); the push-to-`main` [run 36367784465](https://github.com/tap/DspTap/actions/runs/36367784465) reads +22 packets on every scenario, fixed point included (the process-environment offset, `bench/README.md` above; inside the band, not re-recorded) | `d9c1e33` (the #44 squash) | clang 19.1.5 (CodeLinaro, hexagon-unknown-linux-musl) | qemu-hexagon 8.2.2 (built from the release with `--enable-plugins`) |
| hexagon | rfft_q31_512 | — | 157,303,237 | — | **seeded** (the `hexagon` key is new): Hexagon v68 + HVX-128, CodeLinaro clang 19.1.5, static musl, under qemu-hexagon user mode; the key counts packets, not instructions (above). Seeded on the tuned srdif engine; the engine as #42 shipped it read 157,303,237 on the first seeding run ([run 36357132506](https://github.com/tap/DspTap/actions/runs/36357132506), `28befbd`), and the engine it replaced 157,303,237, the same code (local, the same harness, converted to the runner's offset; float ±28). `docs/fft-design.md`, "Hexagon (clang) tuning" | `workflow_dispatch` [run 36363710480](https://github.com/tap/DspTap/actions/runs/36363710480) on the branch (`2a5fc69`, seed mode); the push-to-`main` [run 36367784465](https://github.com/tap/DspTap/actions/runs/36367784465) reads +22 packets on every scenario, fixed point included (the process-environment offset, `bench/README.md` above; inside the band, not re-recorded) | `d9c1e33` (the #44 squash) | clang 19.1.5 (CodeLinaro, hexagon-unknown-linux-musl) | qemu-hexagon 8.2.2 (built from the release with `--enable-plugins`) |
| m4-softfp | bw4_f32_ch2 | — | 27,434,928 | — | **recorded** (the Butterworth stage's scenarios are new; `butterworth.h` as 0.6.0 shipped it, the per-sample `biquad_cascade`). A bracket sum, not a whole run ("The Butterworth stage's scenarios", above) | `workflow_dispatch` [run 45251667308](https://github.com/tap/DspTap/actions/runs/45251667308) on the branch (`f46a688`); the same numbers measured locally before it, to the instruction | `b5f8214` | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m4-softfp | bw4_f32_ch2_n4 | — | 27,436,208 | — | **recorded** (the Butterworth stage's scenarios are new; `butterworth.h` as 0.6.0 shipped it, the per-sample `biquad_cascade`). A bracket sum, not a whole run ("The Butterworth stage's scenarios", above) | `workflow_dispatch` [run 45251667308](https://github.com/tap/DspTap/actions/runs/45251667308) on the branch (`f46a688`); the same numbers measured locally before it, to the instruction | `b5f8214` | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m4-softfp | bw8_f32_ch1 | — | 26,252,841 | — | **recorded** (the Butterworth stage's scenarios are new; `butterworth.h` as 0.6.0 shipped it, the per-sample `biquad_cascade`). A bracket sum, not a whole run ("The Butterworth stage's scenarios", above) | `workflow_dispatch` [run 45251667308](https://github.com/tap/DspTap/actions/runs/45251667308) on the branch (`f46a688`); the same numbers measured locally before it, to the instruction | `b5f8214` | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m4-softfp | bw8_f32_ch2 | — | 53,936,305 | — | **recorded** (the Butterworth stage's scenarios are new; `butterworth.h` as 0.6.0 shipped it, the per-sample `biquad_cascade`). A bracket sum, not a whole run ("The Butterworth stage's scenarios", above) | `workflow_dispatch` [run 45251667308](https://github.com/tap/DspTap/actions/runs/45251667308) on the branch (`f46a688`); the same numbers measured locally before it, to the instruction | `b5f8214` | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m4f | bw4_f32_ch2 | — | 1,837,440 | — | **recorded** (the Butterworth stage's scenarios are new; `butterworth.h` as 0.6.0 shipped it, the per-sample `biquad_cascade`). A bracket sum, not a whole run ("The Butterworth stage's scenarios", above) | `workflow_dispatch` [run 45251667308](https://github.com/tap/DspTap/actions/runs/45251667308) on the branch (`f46a688`); the same numbers measured locally before it, to the instruction | `b5f8214` | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m4f | bw4_f32_ch2_n4 | — | 1,855,488 | — | **recorded** (the Butterworth stage's scenarios are new; `butterworth.h` as 0.6.0 shipped it, the per-sample `biquad_cascade`). A bracket sum, not a whole run ("The Butterworth stage's scenarios", above) | `workflow_dispatch` [run 45251667308](https://github.com/tap/DspTap/actions/runs/45251667308) on the branch (`f46a688`); the same numbers measured locally before it, to the instruction | `b5f8214` | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m4f | bw8_f32_ch1 | — | 1,444,224 | — | **recorded** (the Butterworth stage's scenarios are new; `butterworth.h` as 0.6.0 shipped it, the per-sample `biquad_cascade`). A bracket sum, not a whole run ("The Butterworth stage's scenarios", above) | `workflow_dispatch` [run 45251667308](https://github.com/tap/DspTap/actions/runs/45251667308) on the branch (`f46a688`); the same numbers measured locally before it, to the instruction | `b5f8214` | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m4f | bw8_f32_ch2 | — | 3,148,160 | — | **recorded** (the Butterworth stage's scenarios are new; `butterworth.h` as 0.6.0 shipped it, the per-sample `biquad_cascade`). A bracket sum, not a whole run ("The Butterworth stage's scenarios", above) | `workflow_dispatch` [run 45251667308](https://github.com/tap/DspTap/actions/runs/45251667308) on the branch (`f46a688`); the same numbers measured locally before it, to the instruction | `b5f8214` | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m33 | bw4_f32_ch2 | — | 1,837,440 | — | **recorded** (the Butterworth stage's scenarios are new; `butterworth.h` as 0.6.0 shipped it, the per-sample `biquad_cascade`). A bracket sum, not a whole run ("The Butterworth stage's scenarios", above) | `workflow_dispatch` [run 45251667308](https://github.com/tap/DspTap/actions/runs/45251667308) on the branch (`f46a688`); the same numbers measured locally before it, to the instruction | `b5f8214` | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m33 | bw4_f32_ch2_n4 | — | 1,855,488 | — | **recorded** (the Butterworth stage's scenarios are new; `butterworth.h` as 0.6.0 shipped it, the per-sample `biquad_cascade`). A bracket sum, not a whole run ("The Butterworth stage's scenarios", above) | `workflow_dispatch` [run 45251667308](https://github.com/tap/DspTap/actions/runs/45251667308) on the branch (`f46a688`); the same numbers measured locally before it, to the instruction | `b5f8214` | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m33 | bw8_f32_ch1 | — | 1,444,224 | — | **recorded** (the Butterworth stage's scenarios are new; `butterworth.h` as 0.6.0 shipped it, the per-sample `biquad_cascade`). A bracket sum, not a whole run ("The Butterworth stage's scenarios", above) | `workflow_dispatch` [run 45251667308](https://github.com/tap/DspTap/actions/runs/45251667308) on the branch (`f46a688`); the same numbers measured locally before it, to the instruction | `b5f8214` | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m33 | bw8_f32_ch2 | — | 3,148,160 | — | **recorded** (the Butterworth stage's scenarios are new; `butterworth.h` as 0.6.0 shipped it, the per-sample `biquad_cascade`). A bracket sum, not a whole run ("The Butterworth stage's scenarios", above) | `workflow_dispatch` [run 45251667308](https://github.com/tap/DspTap/actions/runs/45251667308) on the branch (`f46a688`); the same numbers measured locally before it, to the instruction | `b5f8214` | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m55 | bw4_f32_ch2 | — | 1,837,184 | — | **recorded** (the Butterworth stage's scenarios are new; `butterworth.h` as 0.6.0 shipped it, the per-sample `biquad_cascade`). A bracket sum, not a whole run ("The Butterworth stage's scenarios", above) | `workflow_dispatch` [run 45251667308](https://github.com/tap/DspTap/actions/runs/45251667308) on the branch (`f46a688`); the same numbers measured locally before it, to the instruction | `b5f8214` | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m55 | bw4_f32_ch2_n4 | — | 1,908,736 | — | **recorded** (the Butterworth stage's scenarios are new; `butterworth.h` as 0.6.0 shipped it, the per-sample `biquad_cascade`). A bracket sum, not a whole run ("The Butterworth stage's scenarios", above) | `workflow_dispatch` [run 45251667308](https://github.com/tap/DspTap/actions/runs/45251667308) on the branch (`f46a688`); the same numbers measured locally before it, to the instruction | `b5f8214` | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m55 | bw8_f32_ch1 | — | 1,476,992 | — | **recorded** (the Butterworth stage's scenarios are new; `butterworth.h` as 0.6.0 shipped it, the per-sample `biquad_cascade`). A bracket sum, not a whole run ("The Butterworth stage's scenarios", above) | `workflow_dispatch` [run 45251667308](https://github.com/tap/DspTap/actions/runs/45251667308) on the branch (`f46a688`); the same numbers measured locally before it, to the instruction | `b5f8214` | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m55 | bw8_f32_ch2 | — | 3,082,368 | — | **recorded** (the Butterworth stage's scenarios are new; `butterworth.h` as 0.6.0 shipped it, the per-sample `biquad_cascade`). A bracket sum, not a whole run ("The Butterworth stage's scenarios", above) | `workflow_dispatch` [run 45251667308](https://github.com/tap/DspTap/actions/runs/45251667308) on the branch (`f46a688`); the same numbers measured locally before it, to the instruction | `b5f8214` | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m55 | bw8_f32_set_cutoff | — | 3,534,019 | — | **recorded** (the Butterworth stage's scenarios are new; `butterworth.h` as 0.6.0 shipped it, the per-sample `biquad_cascade`). A bracket sum, not a whole run ("The Butterworth stage's scenarios", above) | `workflow_dispatch` [run 45251667308](https://github.com/tap/DspTap/actions/runs/45251667308) on the branch (`f46a688`); the same numbers measured locally before it, to the instruction | `b5f8214` | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m55 | bw8_f32_set_cutoff_deviation | — | 3,534,035 | — | **recorded** (the Butterworth stage's scenarios are new; `butterworth.h` as 0.6.0 shipped it, the per-sample `biquad_cascade`). A bracket sum, not a whole run ("The Butterworth stage's scenarios", above) | `workflow_dispatch` [run 45251667308](https://github.com/tap/DspTap/actions/runs/45251667308) on the branch (`f46a688`); the same numbers measured locally before it, to the instruction | `b5f8214` | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m55-ooura | bw4_f32_ch2 | — | 1,837,184 | — | **recorded** (the Butterworth stage's scenarios are new; `butterworth.h` as 0.6.0 shipped it, the per-sample `biquad_cascade`). A bracket sum, not a whole run ("The Butterworth stage's scenarios", above) | `workflow_dispatch` [run 45251667308](https://github.com/tap/DspTap/actions/runs/45251667308) on the branch (`f46a688`); the same numbers measured locally before it, to the instruction | `b5f8214` | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m55-ooura | bw4_f32_ch2_n4 | — | 1,908,736 | — | **recorded** (the Butterworth stage's scenarios are new; `butterworth.h` as 0.6.0 shipped it, the per-sample `biquad_cascade`). A bracket sum, not a whole run ("The Butterworth stage's scenarios", above) | `workflow_dispatch` [run 45251667308](https://github.com/tap/DspTap/actions/runs/45251667308) on the branch (`f46a688`); the same numbers measured locally before it, to the instruction | `b5f8214` | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m55-ooura | bw8_f32_ch1 | — | 1,476,992 | — | **recorded** (the Butterworth stage's scenarios are new; `butterworth.h` as 0.6.0 shipped it, the per-sample `biquad_cascade`). A bracket sum, not a whole run ("The Butterworth stage's scenarios", above) | `workflow_dispatch` [run 45251667308](https://github.com/tap/DspTap/actions/runs/45251667308) on the branch (`f46a688`); the same numbers measured locally before it, to the instruction | `b5f8214` | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m55-ooura | bw8_f32_ch2 | — | 3,082,368 | — | **recorded** (the Butterworth stage's scenarios are new; `butterworth.h` as 0.6.0 shipped it, the per-sample `biquad_cascade`). A bracket sum, not a whole run ("The Butterworth stage's scenarios", above) | `workflow_dispatch` [run 45251667308](https://github.com/tap/DspTap/actions/runs/45251667308) on the branch (`f46a688`); the same numbers measured locally before it, to the instruction | `b5f8214` | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m55-ooura | bw8_f32_set_cutoff | — | 3,534,019 | — | **recorded** (the Butterworth stage's scenarios are new; `butterworth.h` as 0.6.0 shipped it, the per-sample `biquad_cascade`). A bracket sum, not a whole run ("The Butterworth stage's scenarios", above) | `workflow_dispatch` [run 45251667308](https://github.com/tap/DspTap/actions/runs/45251667308) on the branch (`f46a688`); the same numbers measured locally before it, to the instruction | `b5f8214` | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m55-ooura | bw8_f32_set_cutoff_deviation | — | 3,534,035 | — | **recorded** (the Butterworth stage's scenarios are new; `butterworth.h` as 0.6.0 shipped it, the per-sample `biquad_cascade`). A bracket sum, not a whole run ("The Butterworth stage's scenarios", above) | `workflow_dispatch` [run 45251667308](https://github.com/tap/DspTap/actions/runs/45251667308) on the branch (`f46a688`); the same numbers measured locally before it, to the instruction | `b5f8214` | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| hexagon | bw4_f32_ch2 | — | 739,968 | — | **recorded** (the Butterworth stage's scenarios are new; `butterworth.h` as 0.6.0 shipped it, the per-sample `biquad_cascade`). A bracket sum, not a whole run ("The Butterworth stage's scenarios", above); packets, as every `hexagon` row | `workflow_dispatch` [run 45251667308](https://github.com/tap/DspTap/actions/runs/45251667308) on the branch (`f46a688`); the same numbers measured locally before it, to the instruction | `b5f8214` | clang 19.1.5 (CodeLinaro, hexagon-unknown-linux-musl) | qemu-hexagon 8.2.2 (built from the release with `--enable-plugins`) |
| hexagon | bw4_f32_ch2_n4 | — | 823,296 | — | **recorded** (the Butterworth stage's scenarios are new; `butterworth.h` as 0.6.0 shipped it, the per-sample `biquad_cascade`). A bracket sum, not a whole run ("The Butterworth stage's scenarios", above); packets, as every `hexagon` row | `workflow_dispatch` [run 45251667308](https://github.com/tap/DspTap/actions/runs/45251667308) on the branch (`f46a688`); the same numbers measured locally before it, to the instruction | `b5f8214` | clang 19.1.5 (CodeLinaro, hexagon-unknown-linux-musl) | qemu-hexagon 8.2.2 (built from the release with `--enable-plugins`) |
| hexagon | bw8_f32_ch1 | — | 526,080 | — | **recorded** (the Butterworth stage's scenarios are new; `butterworth.h` as 0.6.0 shipped it, the per-sample `biquad_cascade`). A bracket sum, not a whole run ("The Butterworth stage's scenarios", above); packets, as every `hexagon` row | `workflow_dispatch` [run 45251667308](https://github.com/tap/DspTap/actions/runs/45251667308) on the branch (`f46a688`); the same numbers measured locally before it, to the instruction | `b5f8214` | clang 19.1.5 (CodeLinaro, hexagon-unknown-linux-musl) | qemu-hexagon 8.2.2 (built from the release with `--enable-plugins`) |
| hexagon | bw8_f32_ch2 | — | 1,034,752 | — | **recorded** (the Butterworth stage's scenarios are new; `butterworth.h` as 0.6.0 shipped it, the per-sample `biquad_cascade`). A bracket sum, not a whole run ("The Butterworth stage's scenarios", above); packets, as every `hexagon` row | `workflow_dispatch` [run 45251667308](https://github.com/tap/DspTap/actions/runs/45251667308) on the branch (`f46a688`); the same numbers measured locally before it, to the instruction | `b5f8214` | clang 19.1.5 (CodeLinaro, hexagon-unknown-linux-musl) | qemu-hexagon 8.2.2 (built from the release with `--enable-plugins`) |
| hexagon | bw8_f32_set_cutoff | — | 13,726,637 | — | **recorded** (the Butterworth stage's scenarios are new; `butterworth.h` as 0.6.0 shipped it, the per-sample `biquad_cascade`). A bracket sum, not a whole run ("The Butterworth stage's scenarios", above); packets, as every `hexagon` row | `workflow_dispatch` [run 45251667308](https://github.com/tap/DspTap/actions/runs/45251667308) on the branch (`f46a688`); the same numbers measured locally before it, to the instruction | `b5f8214` | clang 19.1.5 (CodeLinaro, hexagon-unknown-linux-musl) | qemu-hexagon 8.2.2 (built from the release with `--enable-plugins`) |
| hexagon | bw8_f32_set_cutoff_deviation | — | 13,726,653 | — | **recorded** (the Butterworth stage's scenarios are new; `butterworth.h` as 0.6.0 shipped it, the per-sample `biquad_cascade`). A bracket sum, not a whole run ("The Butterworth stage's scenarios", above); packets, as every `hexagon` row | `workflow_dispatch` [run 45251667308](https://github.com/tap/DspTap/actions/runs/45251667308) on the branch (`f46a688`); the same numbers measured locally before it, to the instruction | `b5f8214` | clang 19.1.5 (CodeLinaro, hexagon-unknown-linux-musl) | qemu-hexagon 8.2.2 (built from the release with `--enable-plugins`) |
| m4-softfp | bw4_f32_ch2 | 27,434,928 | 26,049,456 | -5.05 % | **re-recorded at the block kernel**: `butterworth.h`'s `process()` runs `iir_kernels.h`'s `biquad_cascade_block` (the coefficients and state in locals for the call, a compile-time section count, pairs of channels in float at orders 2 and 4), where it called the per-sample `biquad_cascade` per channel per frame. No output bit moves with contraction off (`tests/test_iir_block.cpp`); every layout pair bit-identical on every key. "The block kernel", below | `workflow_dispatch` [run 45283542111](https://github.com/tap/DspTap/actions/runs/45283542111) on the branch (`5634a73`); the same numbers measured locally before it, to the instruction | `e183a79` | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m4-softfp | bw4_f32_ch2_n4 | 27,436,208 | 26,240,176 | -4.36 % | **re-recorded at the block kernel**: `butterworth.h`'s `process()` runs `iir_kernels.h`'s `biquad_cascade_block` (the coefficients and state in locals for the call, a compile-time section count, pairs of channels in float at orders 2 and 4), where it called the per-sample `biquad_cascade` per channel per frame. No output bit moves with contraction off (`tests/test_iir_block.cpp`); every layout pair bit-identical on every key. "The block kernel", below | `workflow_dispatch` [run 45283542111](https://github.com/tap/DspTap/actions/runs/45283542111) on the branch (`5634a73`); the same numbers measured locally before it, to the instruction | `e183a79` | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m4-softfp | bw8_f32_ch1 | 26,252,841 | 25,968,681 | -1.08 % | **re-recorded at the block kernel**: `butterworth.h`'s `process()` runs `iir_kernels.h`'s `biquad_cascade_block` (the coefficients and state in locals for the call, a compile-time section count, pairs of channels in float at orders 2 and 4), where it called the per-sample `biquad_cascade` per channel per frame. No output bit moves with contraction off (`tests/test_iir_block.cpp`); every layout pair bit-identical on every key. "The block kernel", below | `workflow_dispatch` [run 45283542111](https://github.com/tap/DspTap/actions/runs/45283542111) on the branch (`5634a73`); the same numbers measured locally before it, to the instruction | `e183a79` | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m4-softfp | bw8_f32_ch2 | 53,936,305 | 51,994,801 | -3.60 % | **re-recorded at the block kernel**: `butterworth.h`'s `process()` runs `iir_kernels.h`'s `biquad_cascade_block` (the coefficients and state in locals for the call, a compile-time section count, pairs of channels in float at orders 2 and 4), where it called the per-sample `biquad_cascade` per channel per frame. No output bit moves with contraction off (`tests/test_iir_block.cpp`); every layout pair bit-identical on every key. "The block kernel", below | `workflow_dispatch` [run 45283542111](https://github.com/tap/DspTap/actions/runs/45283542111) on the branch (`5634a73`); the same numbers measured locally before it, to the instruction | `e183a79` | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m4f | bw4_f32_ch2 | 1,837,440 | 628,608 | -65.79 % | **re-recorded at the block kernel**: `butterworth.h`'s `process()` runs `iir_kernels.h`'s `biquad_cascade_block` (the coefficients and state in locals for the call, a compile-time section count, pairs of channels in float at orders 2 and 4), where it called the per-sample `biquad_cascade` per channel per frame. No output bit moves with contraction off (`tests/test_iir_block.cpp`); every layout pair bit-identical on every key. "The block kernel", below | `workflow_dispatch` [run 45283542111](https://github.com/tap/DspTap/actions/runs/45283542111) on the branch (`5634a73`); the same numbers measured locally before it, to the instruction | `e183a79` | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m4f | bw4_f32_ch2_n4 | 1,855,488 | 667,648 | -64.02 % | **re-recorded at the block kernel**: `butterworth.h`'s `process()` runs `iir_kernels.h`'s `biquad_cascade_block` (the coefficients and state in locals for the call, a compile-time section count, pairs of channels in float at orders 2 and 4), where it called the per-sample `biquad_cascade` per channel per frame. No output bit moves with contraction off (`tests/test_iir_block.cpp`); every layout pair bit-identical on every key. "The block kernel", below | `workflow_dispatch` [run 45283542111](https://github.com/tap/DspTap/actions/runs/45283542111) on the branch (`5634a73`); the same numbers measured locally before it, to the instruction | `e183a79` | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m4f | bw8_f32_ch1 | 1,444,224 | 600,064 | -58.45 % | **re-recorded at the block kernel**: `butterworth.h`'s `process()` runs `iir_kernels.h`'s `biquad_cascade_block` (the coefficients and state in locals for the call, a compile-time section count, pairs of channels in float at orders 2 and 4), where it called the per-sample `biquad_cascade` per channel per frame. No output bit moves with contraction off (`tests/test_iir_block.cpp`); every layout pair bit-identical on every key. "The block kernel", below | `workflow_dispatch` [run 45283542111](https://github.com/tap/DspTap/actions/runs/45283542111) on the branch (`5634a73`); the same numbers measured locally before it, to the instruction | `e183a79` | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m4f | bw8_f32_ch2 | 3,148,160 | 1,264,896 | -59.82 % | **re-recorded at the block kernel**: `butterworth.h`'s `process()` runs `iir_kernels.h`'s `biquad_cascade_block` (the coefficients and state in locals for the call, a compile-time section count, pairs of channels in float at orders 2 and 4), where it called the per-sample `biquad_cascade` per channel per frame. No output bit moves with contraction off (`tests/test_iir_block.cpp`); every layout pair bit-identical on every key. "The block kernel", below | `workflow_dispatch` [run 45283542111](https://github.com/tap/DspTap/actions/runs/45283542111) on the branch (`5634a73`); the same numbers measured locally before it, to the instruction | `e183a79` | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m33 | bw4_f32_ch2 | 1,837,440 | 661,376 | -64.01 % | **re-recorded at the block kernel**: `butterworth.h`'s `process()` runs `iir_kernels.h`'s `biquad_cascade_block` (the coefficients and state in locals for the call, a compile-time section count, pairs of channels in float at orders 2 and 4), where it called the per-sample `biquad_cascade` per channel per frame. No output bit moves with contraction off (`tests/test_iir_block.cpp`); every layout pair bit-identical on every key. "The block kernel", below | `workflow_dispatch` [run 45283542111](https://github.com/tap/DspTap/actions/runs/45283542111) on the branch (`5634a73`); the same numbers measured locally before it, to the instruction | `e183a79` | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m33 | bw4_f32_ch2_n4 | 1,855,488 | 659,456 | -64.46 % | **re-recorded at the block kernel**: `butterworth.h`'s `process()` runs `iir_kernels.h`'s `biquad_cascade_block` (the coefficients and state in locals for the call, a compile-time section count, pairs of channels in float at orders 2 and 4), where it called the per-sample `biquad_cascade` per channel per frame. No output bit moves with contraction off (`tests/test_iir_block.cpp`); every layout pair bit-identical on every key. "The block kernel", below | `workflow_dispatch` [run 45283542111](https://github.com/tap/DspTap/actions/runs/45283542111) on the branch (`5634a73`); the same numbers measured locally before it, to the instruction | `e183a79` | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m33 | bw8_f32_ch1 | 1,444,224 | 649,344 | -55.04 % | **re-recorded at the block kernel**: `butterworth.h`'s `process()` runs `iir_kernels.h`'s `biquad_cascade_block` (the coefficients and state in locals for the call, a compile-time section count, pairs of channels in float at orders 2 and 4), where it called the per-sample `biquad_cascade` per channel per frame. No output bit moves with contraction off (`tests/test_iir_block.cpp`); every layout pair bit-identical on every key. "The block kernel", below | `workflow_dispatch` [run 45283542111](https://github.com/tap/DspTap/actions/runs/45283542111) on the branch (`5634a73`); the same numbers measured locally before it, to the instruction | `e183a79` | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m33 | bw8_f32_ch2 | 3,148,160 | 1,363,456 | -56.69 % | **re-recorded at the block kernel**: `butterworth.h`'s `process()` runs `iir_kernels.h`'s `biquad_cascade_block` (the coefficients and state in locals for the call, a compile-time section count, pairs of channels in float at orders 2 and 4), where it called the per-sample `biquad_cascade` per channel per frame. No output bit moves with contraction off (`tests/test_iir_block.cpp`); every layout pair bit-identical on every key. "The block kernel", below | `workflow_dispatch` [run 45283542111](https://github.com/tap/DspTap/actions/runs/45283542111) on the branch (`5634a73`); the same numbers measured locally before it, to the instruction | `e183a79` | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m55 | bw4_f32_ch2 | 1,837,184 | 629,632 | -65.73 % | **re-recorded at the block kernel**: `butterworth.h`'s `process()` runs `iir_kernels.h`'s `biquad_cascade_block` (the coefficients and state in locals for the call, a compile-time section count, pairs of channels in float at orders 2 and 4), where it called the per-sample `biquad_cascade` per channel per frame. No output bit moves with contraction off (`tests/test_iir_block.cpp`); every layout pair bit-identical on every key. "The block kernel", below | `workflow_dispatch` [run 45283542111](https://github.com/tap/DspTap/actions/runs/45283542111) on the branch (`5634a73`); the same numbers measured locally before it, to the instruction | `e183a79` | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m55 | bw4_f32_ch2_n4 | 1,908,736 | 684,032 | -64.16 % | **re-recorded at the block kernel**: `butterworth.h`'s `process()` runs `iir_kernels.h`'s `biquad_cascade_block` (the coefficients and state in locals for the call, a compile-time section count, pairs of channels in float at orders 2 and 4), where it called the per-sample `biquad_cascade` per channel per frame. No output bit moves with contraction off (`tests/test_iir_block.cpp`); every layout pair bit-identical on every key. "The block kernel", below | `workflow_dispatch` [run 45283542111](https://github.com/tap/DspTap/actions/runs/45283542111) on the branch (`5634a73`); the same numbers measured locally before it, to the instruction | `e183a79` | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m55 | bw8_f32_ch1 | 1,476,992 | 633,216 | -57.13 % | **re-recorded at the block kernel**: `butterworth.h`'s `process()` runs `iir_kernels.h`'s `biquad_cascade_block` (the coefficients and state in locals for the call, a compile-time section count, pairs of channels in float at orders 2 and 4), where it called the per-sample `biquad_cascade` per channel per frame. No output bit moves with contraction off (`tests/test_iir_block.cpp`); every layout pair bit-identical on every key. "The block kernel", below | `workflow_dispatch` [run 45283542111](https://github.com/tap/DspTap/actions/runs/45283542111) on the branch (`5634a73`); the same numbers measured locally before it, to the instruction | `e183a79` | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m55 | bw8_f32_ch2 | 3,082,368 | 1,331,456 | -56.80 % | **re-recorded at the block kernel**: `butterworth.h`'s `process()` runs `iir_kernels.h`'s `biquad_cascade_block` (the coefficients and state in locals for the call, a compile-time section count, pairs of channels in float at orders 2 and 4), where it called the per-sample `biquad_cascade` per channel per frame. No output bit moves with contraction off (`tests/test_iir_block.cpp`); every layout pair bit-identical on every key. "The block kernel", below | `workflow_dispatch` [run 45283542111](https://github.com/tap/DspTap/actions/runs/45283542111) on the branch (`5634a73`); the same numbers measured locally before it, to the instruction | `e183a79` | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m55-ooura | bw4_f32_ch2 | 1,837,184 | 629,632 | -65.73 % | **re-recorded at the block kernel**: `butterworth.h`'s `process()` runs `iir_kernels.h`'s `biquad_cascade_block` (the coefficients and state in locals for the call, a compile-time section count, pairs of channels in float at orders 2 and 4), where it called the per-sample `biquad_cascade` per channel per frame. No output bit moves with contraction off (`tests/test_iir_block.cpp`); every layout pair bit-identical on every key. "The block kernel", below | `workflow_dispatch` [run 45283542111](https://github.com/tap/DspTap/actions/runs/45283542111) on the branch (`5634a73`); the same numbers measured locally before it, to the instruction | `e183a79` | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m55-ooura | bw4_f32_ch2_n4 | 1,908,736 | 684,032 | -64.16 % | **re-recorded at the block kernel**: `butterworth.h`'s `process()` runs `iir_kernels.h`'s `biquad_cascade_block` (the coefficients and state in locals for the call, a compile-time section count, pairs of channels in float at orders 2 and 4), where it called the per-sample `biquad_cascade` per channel per frame. No output bit moves with contraction off (`tests/test_iir_block.cpp`); every layout pair bit-identical on every key. "The block kernel", below | `workflow_dispatch` [run 45283542111](https://github.com/tap/DspTap/actions/runs/45283542111) on the branch (`5634a73`); the same numbers measured locally before it, to the instruction | `e183a79` | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m55-ooura | bw8_f32_ch1 | 1,476,992 | 633,216 | -57.13 % | **re-recorded at the block kernel**: `butterworth.h`'s `process()` runs `iir_kernels.h`'s `biquad_cascade_block` (the coefficients and state in locals for the call, a compile-time section count, pairs of channels in float at orders 2 and 4), where it called the per-sample `biquad_cascade` per channel per frame. No output bit moves with contraction off (`tests/test_iir_block.cpp`); every layout pair bit-identical on every key. "The block kernel", below | `workflow_dispatch` [run 45283542111](https://github.com/tap/DspTap/actions/runs/45283542111) on the branch (`5634a73`); the same numbers measured locally before it, to the instruction | `e183a79` | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m55-ooura | bw8_f32_ch2 | 3,082,368 | 1,331,456 | -56.80 % | **re-recorded at the block kernel**: `butterworth.h`'s `process()` runs `iir_kernels.h`'s `biquad_cascade_block` (the coefficients and state in locals for the call, a compile-time section count, pairs of channels in float at orders 2 and 4), where it called the per-sample `biquad_cascade` per channel per frame. No output bit moves with contraction off (`tests/test_iir_block.cpp`); every layout pair bit-identical on every key. "The block kernel", below | `workflow_dispatch` [run 45283542111](https://github.com/tap/DspTap/actions/runs/45283542111) on the branch (`5634a73`); the same numbers measured locally before it, to the instruction | `e183a79` | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| hexagon | bw4_f32_ch2 | 739,968 | 209,024 | -71.75 % | **re-recorded at the block kernel**: `butterworth.h`'s `process()` runs `iir_kernels.h`'s `biquad_cascade_block` (the coefficients and state in locals for the call, a compile-time section count, pairs of channels in float at orders 2 and 4), where it called the per-sample `biquad_cascade` per channel per frame. No output bit moves with contraction off (`tests/test_iir_block.cpp`); every layout pair bit-identical on every key. "The block kernel", below | `workflow_dispatch` [run 45283542111](https://github.com/tap/DspTap/actions/runs/45283542111) on the branch (`5634a73`); the same numbers measured locally before it, to the instruction | `e183a79` | clang 19.1.5 (CodeLinaro, hexagon-unknown-linux-musl) | qemu-hexagon 8.2.2 (built from the release with `--enable-plugins`) |
| hexagon | bw4_f32_ch2_n4 | 823,296 | 339,968 | -58.71 % | **re-recorded at the block kernel**: `butterworth.h`'s `process()` runs `iir_kernels.h`'s `biquad_cascade_block` (the coefficients and state in locals for the call, a compile-time section count, pairs of channels in float at orders 2 and 4), where it called the per-sample `biquad_cascade` per channel per frame. No output bit moves with contraction off (`tests/test_iir_block.cpp`); every layout pair bit-identical on every key. "The block kernel", below | `workflow_dispatch` [run 45283542111](https://github.com/tap/DspTap/actions/runs/45283542111) on the branch (`5634a73`); the same numbers measured locally before it, to the instruction | `e183a79` | clang 19.1.5 (CodeLinaro, hexagon-unknown-linux-musl) | qemu-hexagon 8.2.2 (built from the release with `--enable-plugins`) |
| hexagon | bw8_f32_ch1 | 526,080 | 268,928 | -48.88 % | **re-recorded at the block kernel**: `butterworth.h`'s `process()` runs `iir_kernels.h`'s `biquad_cascade_block` (the coefficients and state in locals for the call, a compile-time section count, pairs of channels in float at orders 2 and 4), where it called the per-sample `biquad_cascade` per channel per frame. No output bit moves with contraction off (`tests/test_iir_block.cpp`); every layout pair bit-identical on every key. "The block kernel", below | `workflow_dispatch` [run 45283542111](https://github.com/tap/DspTap/actions/runs/45283542111) on the branch (`5634a73`); the same numbers measured locally before it, to the instruction | `e183a79` | clang 19.1.5 (CodeLinaro, hexagon-unknown-linux-musl) | qemu-hexagon 8.2.2 (built from the release with `--enable-plugins`) |
| hexagon | bw8_f32_ch2 | 1,034,752 | 537,472 | -48.06 % | **re-recorded at the block kernel**: `butterworth.h`'s `process()` runs `iir_kernels.h`'s `biquad_cascade_block` (the coefficients and state in locals for the call, a compile-time section count, pairs of channels in float at orders 2 and 4), where it called the per-sample `biquad_cascade` per channel per frame. No output bit moves with contraction off (`tests/test_iir_block.cpp`); every layout pair bit-identical on every key. "The block kernel", below | `workflow_dispatch` [run 45283542111](https://github.com/tap/DspTap/actions/runs/45283542111) on the branch (`5634a73`); the same numbers measured locally before it, to the instruction | `e183a79` | clang 19.1.5 (CodeLinaro, hexagon-unknown-linux-musl) | qemu-hexagon 8.2.2 (built from the release with `--enable-plugins`) |

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
`arm-none-eabi-size -A` (the toolchain's `llvm-size -A` on the `hexagon` key) on the MinSizeRel probe; the ceiling is what
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
| m4-softfp | `rfft_f32_512` | 31,081 | 32,064 | **re-recorded at the srdif replacement** (tap/DspTap#42) (ceiling 54,912 before): measured + 3 %, rounded up to 64 bytes; the srdif engine builds its tables from integer arithmetic, so the float probe no longer links libm's double `sin` / `cos`; the figure is the fix pass's for review A of #42 (no swap-list builder, one table allocation), 1,112 … 1,136 B below the first srdif record | measured locally as `bench.yml` builds the probe (MinSizeRel, `size -A` `.text`), 2026-09-27; confirmed to the byte by the push-to-`main` [run 36348787347](https://github.com/tap/DspTap/actions/runs/36348787347) on `72977aa` | 13.2.1 (15:13.2.rel1-2) |
| m4f | `rfft_f32_512` | 28,921 | 29,824 | **re-recorded at the srdif replacement** (tap/DspTap#42) (ceiling 45,952 before): measured + 3 %, rounded up to 64 bytes; the srdif engine builds its tables from integer arithmetic, so the float probe no longer links libm's double `sin` / `cos`; the figure is the fix pass's for review A of #42 (no swap-list builder, one table allocation), 1,112 … 1,136 B below the first srdif record | measured locally as `bench.yml` builds the probe (MinSizeRel, `size -A` `.text`), 2026-09-27; confirmed to the byte by the push-to-`main` [run 36348787347](https://github.com/tap/DspTap/actions/runs/36348787347) on `72977aa` | 13.2.1 (15:13.2.rel1-2) |
| m33 | `rfft_f32_512` | 28,305 | 29,184 | **re-recorded at the srdif replacement** (tap/DspTap#42) (ceiling 45,376 before): measured + 3 %, rounded up to 64 bytes; the srdif engine builds its tables from integer arithmetic, so the float probe no longer links libm's double `sin` / `cos`; the figure is the fix pass's for review A of #42 (no swap-list builder, one table allocation), 1,112 … 1,136 B below the first srdif record | measured locally as `bench.yml` builds the probe (MinSizeRel, `size -A` `.text`), 2026-09-27; confirmed to the byte by the push-to-`main` [run 36348787347](https://github.com/tap/DspTap/actions/runs/36348787347) on `72977aa` | 13.2.1 (15:13.2.rel1-2) |
| m55-ooura | `rfft_f32_512` | 28,161 | 29,056 | **re-recorded at the srdif replacement** (tap/DspTap#42) (ceiling 40,512 before): measured + 3 %, rounded up to 64 bytes; the srdif engine builds its tables from integer arithmetic, so the float probe no longer links libm's double `sin` / `cos`; the figure is the fix pass's for review A of #42 (no swap-list builder, one table allocation), 1,112 … 1,136 B below the first srdif record | measured locally as `bench.yml` builds the probe (MinSizeRel, `size -A` `.text`), 2026-09-27; confirmed to the byte by the push-to-`main` [run 36348787347](https://github.com/tap/DspTap/actions/runs/36348787347) on `72977aa` | 13.2.1 (15:13.2.rel1-2) |
| hexagon | `rfft_f32_512` | 238,628 | 245,824 | **recorded for the new key**: measured + 3 %, rounded up to 64 bytes; `llvm-size -A` of the statically linked probe, so the figure includes musl libc and libc++ (most of it), not DspTap alone | `workflow_dispatch` [run 36363710480](https://github.com/tap/DspTap/actions/runs/36363710480) (`2a5fc69`; MinSizeRel) | clang 19.1.5 (CodeLinaro) |
| hexagon | `rfft_q15_512` | 244,580 | 251,968 | **recorded for the new key**: measured + 3 %, rounded up to 64 bytes; `llvm-size -A` of the statically linked probe, so the figure includes musl libc and libc++ (most of it), not DspTap alone | `workflow_dispatch` [run 36363710480](https://github.com/tap/DspTap/actions/runs/36363710480) (`2a5fc69`; MinSizeRel) | clang 19.1.5 (CodeLinaro) |
| hexagon | `rfft_q31_512` | 244,132 | 251,456 | **recorded for the new key**: measured + 3 %, rounded up to 64 bytes; `llvm-size -A` of the statically linked probe, so the figure includes musl libc and libc++ (most of it), not DspTap alone | `workflow_dispatch` [run 36363710480](https://github.com/tap/DspTap/actions/runs/36363710480) (`2a5fc69`; MinSizeRel) | clang 19.1.5 (CodeLinaro) |

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
  check, not a warning. Once seeded, the six `icount <key>` jobs are the
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
| srdif (the floating engine replaced, tap/DspTap#42) | the replacement costing more than the engine it replaced on any key: the gate was "no float scenario above its baseline", met on all four portable-engine keys (−2.2 … −6.5 % at the first srdif record, −2.2 … −7.9 % after the fix pass for review A of #42) before the counts were re-recorded (table above) |

**Stage 4 (engine parameter + ABI tag), tap/DspTap#35: what the ratchet
caught, and a harness defect fixed without a re-record.** Measured locally
with the CI toolchain (arm-none-eabi-gcc 13.2.1, qemu 8.2.2, the same
plugin), and reproduced to the instruction by the hostile review. After the
backends moved out of `fft.h` and every engine took its size in its
constructor, the `m55` key read `rfft_f32_512` 58,667,632 vs 52,382,331
(+12.00 %) and `rfft_f32_2048` 61,148,069 vs 54,858,120 (+11.47 %), with
every fixed-point scenario at +0.00 % and both float checksums identical to
the baseline binaries'. The CMSIS archive was byte-identical and the
wrapper's loops identical instruction for instruction. **The defect was the
harness.** The count includes the bench's own FNV-1a fold, and the fold's
codegen depended on what was inlined before it in `run()`: the recorded
baselines were taken when the CMSIS engine's construction reached the
workload through an out-of-line `make_floating_engine<float>`; with the
constructor inlined into `run()` (vector allocation, memset,
`arm_rfft_fast_init_f32`, EH cleanup) GCC's register allocation changed and
the fold went from 7 to 10 instructions per element (the 64-bit hash's low
word spilled, the FNV prime rematerialized), i.e. +6 per sample per
iteration, exactly the delta. A first version of the fix put a `noinline`
attribute on the two accelerated engines' constructors in shipping code
(m55 +5 / +45); the review measured that the harness fix alone lands
closer and costs nothing in the library, and that the attribute added 860
bytes of MinSizeRel `.text` to the one-construction m55 f32 probe. So:
`bench/icount/icount_main.cpp` constructs the transform under test through
a non-inlined `make_fft()` in both `run()` bodies — the baseline's own
shape — and shipping code carries no benchmark-shaped attribute. With that
and nothing else (all five keys, `rfft_f32_512` / `rfft_f32_2048` /
`rfft_q15_512` / `rfft_q31_512` / `rfft_q31_2048`, delta vs baseline):
m55 −2 / +38 / +5 / +5 / +5; m33 −10,236 (−0.01 %) / −5,116 (−0.00 %) /
−2,043 / +5 / +5; m55-ooura −8,193 (−0.01 %) / −1,025 (−0.00 %) / +5 / +5 /
+5; m4f +6 / +6 / +5 / +5 / +5; m4-softfp +4 / −1,020 / +5 / +5 / +5. The
+5 on every fixed-point scenario and on the M4 float scenarios is the
factory call itself. No baseline moves. `.text` probes (MinSizeRel, `size
-A`), all under their ceilings and none re-recorded: m55 f32 **106,725**
(−932 vs the recorded 107,657), Q15 27,105 (+8), Q31 26,553 (−40); m33
44,009 / 31,113 / 30,593 (+8 / +8 / −40); m55-ooura 39,281 / 27,105 /
26,553 (0 / +8 / −40); m4f 44,601 / 31,745 / 31,209 (0 / +64 / 0);
m4-softfp identical to the record. Lesson, recorded in
`docs/audit-fft-and-code-smells.md` Part 13: a workload that constructs its
engine in the function it transforms in measures its own register
allocation as much as the transform; this harness no longer does.

## The block kernel (`biquad_cascade_block`)

Through 0.6.0 the Butterworth stage called `iir_kernels.h`'s per-sample
`biquad_cascade` for every channel of every frame. Its section count is a
run-time number and a store to the output may alias the coefficient table
and the state, so the compiled loop reloaded five coefficients and two
state words per section per sample (read in the x86-64, Cortex-M33 and
Cortex-M4F disassembly). `biquad_cascade_block` runs the same expressions
with a compile-time section count, the coefficients and a group's states in
locals for the call, and the group's channels through the whole block.

Its count as a fraction of the per-sample stage's, measured locally with
`scripts/icount.py` as `bench.yml` runs it (bracket sums, 16,384 frames,
calls of 128 frames unless stated; GCC 13.2.1 and QEMU 8.2.2 for the Arm
keys, CodeLinaro clang 19.1.5 and qemu-hexagon 8.2.2; 2026-11-23; the
`m55` column with the srdif FFT, which these binaries do not run):

| stage | `m4f` | `m33` | `m55` | `hexagon` |
|---|---:|---:|---:|---:|
| order 2, 1 ch | 0.40 | 0.43 | 0.34 | 0.34 |
| order 2, 2 ch | 0.31 | 0.34 | 0.31 | 0.29 |
| order 2, 6 ch | 0.35 | 0.38 | 0.35 | 0.32 |
| order 2, 12 ch | 0.36 | 0.39 | 0.36 | 0.28 |
| order 4, 1 ch | 0.40 | 0.42 | 0.37 | 0.35 |
| order 4, 2 ch | 0.34 | 0.36 | 0.34 | 0.28 |
| order 4, 6 ch | 0.37 | 0.39 | 0.37 | 0.30 |
| order 4, 12 ch | 0.38 | 0.40 | 0.37 | 0.26 |
| order 6, 1 ch | 0.40 | 0.42 | 0.39 | 0.38 |
| order 6, 2 ch | 0.39 | 0.40 | 0.39 | 0.39 |
| order 6, 6 ch | 0.41 | 0.43 | 0.41 | 0.40 |
| order 6, 12 ch | 0.42 | 0.43 | 0.42 | 0.32 |
| order 8, 1 ch | 0.42 | 0.45 | 0.43 | 0.51 |
| order 8, 2 ch | 0.40 | 0.43 | 0.43 | 0.52 |
| order 8, 6 ch | 0.42 | 0.45 | 0.45 | 0.52 |
| order 8, 12 ch | 0.42 | 0.46 | 0.45 | 0.48 |
| order 16, 1 ch | 0.60 | 0.59 | 0.59 | 0.67 |
| order 16, 2 ch | 0.58 | 0.58 | 0.59 | 0.67 |
| order 16, 6 ch | 0.60 | 0.60 | 0.60 | 0.67 |
| order 16, 12 ch | 0.60 | 0.60 | 0.60 | 0.64 |
| order 4, 2 ch, calls of 16 | 0.36 | 0.38 | 0.37 | 0.31 |
| order 8, 2 ch, calls of 16 | 0.44 | 0.47 | 0.48 | 0.60 |
| order 4, 2 ch, calls of 4 | 0.36 | 0.36 | 0.36 | 0.41 |
| order 8, 2 ch, calls of 4 | 0.46 | 0.47 | 0.47 | 0.85 |
| order 4, 1 ch, calls of 1 | 0.79 | 0.79 | 0.84 | 1.03 |
| order 4, 2 ch, calls of 1 | 0.63 | 0.63 | 0.68 | 0.68 |
| order 8, 1 ch, calls of 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| order 8, 2 ch, calls of 1 | 1.01 | 1.01 | 1.03 | 0.88 |
| order 4, 1 ch, double | — | — | 0.37 | 0.72 |
| order 4, 2 ch, double | — | — | 0.38 | 0.72 |
| order 4, 6 ch, double | — | — | 0.41 | 0.70 |
| order 8, 1 ch, double | — | — | 0.53 | 0.88 |
| order 8, 2 ch, double | — | — | 0.53 | 0.85 |
| order 8, 6 ch, double | — | — | 0.55 | 0.83 |

Double is measured where the target computes it in hardware. The gated
scenarios' before and after are in the recorded table above.

**The group** (`detail::k_biquad_block_group`): two channels together for
float at orders 2 and 4, one at a time for float from order 6 and for
double at every order. Measured over groups of 1, 2 and 4 on the same
matrix, with the group loop still inline: pairs are 3 to 6 points cheaper
than singles at orders 2 and 4 (one order-4 channel is 14 values, a pair
18, inside the 32 single registers of an M-class FPU); at order 6 the two
are within 4 points on the M-class cores and singles 5 to 9 points cheaper
on Hexagon; from order 8 singles win on the M-class cores by 4 to 15
points and on Hexagon by 6 to 8 at order 8 (within 2 at order 16); fours
never win by more than 3 points and lose up to 12. In double (16 registers
on the M55) singles win on the M55 by 4 to 9 points and tie on Hexagon.

**Tried and not kept:**

- *Channels in vector lanes* (every channel through one section, then the
  next): at 12 channels 3.4 times under planar on x86-64 and 30 % under on
  the M55, but GCC vectorized it at 12 channels only, it is worse on the
  M4F and Hexagon, and its last bits differ from the scalar path's under
  contraction.
- *The coefficients through a `restrict` pointer instead of a copy:* on
  Hexagon (clang) the same or a few points better; on the M33 (GCC) 15 to
  20 points worse at every block size.
- *The state in place through a `restrict` pointer* (no copies): erratic on
  Hexagon, 0.50 to 1.27 of the per-sample count.
- *The group loop inlined into `process()`:* the same counts within 5
  points where the stage is called once, but inlined into a caller's loop
  over mono stages GCC on the M55 vectorized across sections: 1.48 times
  the count of the same stages interleaved, and other last bits under
  contraction (the layout bench's order-8 pairs read `DIFFERENT`). Out of
  line, every layout pair is bit-identical on every key.

**Short calls.** A single frame through four or more sections goes to the
per-sample kernel (the block's copies cost 1.49 times the per-sample count
there on Hexagon, one channel, before this rule). The rule keys on the
frame count alone, so under contraction a stage's last bits may depend on
how its input is chunked; with contraction off they cannot
(`tests/test_iir_block.cpp`).

## Layout: interleaved against planar (informational, never a gate)

`bench/layout/` counts 1:1 chain stages in both channel layouts on every key
the ratchet runs, so the question "does a planar host lose by calling an
interleaved API?" has a number per core. It is informational by decision
(the maintainer's, 2026-11-21): no baselines, no tolerance, no failure.

- **Binaries:** `tap_dsp_layout_<kind>_f32_ch<C>_<i|p>`, kind `bw4` or `bw8`
  (the Butterworth at order 4 or 8, fc = 0.1) or `gain` (a per-channel gain,
  `bench_gain` in `layout_main.cpp`: the stateless case, bench-local, not a
  primitive); C = 1, 2, 6, 12 (mono, stereo, 5.1, 7.1.4); `i` = one stage of
  C channels over interleaved frames, `p` = C mono stages over contiguous
  channel buffers. One channel builds `_i` only. 16,384 frames per channel
  in calls of 128, float.
- **What a count includes:** the bracket sum, as for the Butterworth
  scenarios above: each 128-frame step (one `process()` interleaved, C calls
  planar) between the markers. The worst single step is printed beside it.
- **Agreement:** both layouts of a pair print the same checksum where the
  compiler contracts the two instantiations alike; the job prints
  `bit-identical` or `DIFFERENT` per pair, and on the host ctest pins it
  (`tap_dsp_layout_<pair>_pair_agrees`, beside a smoke test per binary).
- **The job:** the `Layout counts (informational, never a gate)` step of
  `bench.yml`, after the ratchet's count, whatever its verdict:
  `icount.py --prefix tap_dsp_layout_ --report-only` prints every count and
  the p / i table to the step summary and uploads `layout-<key>.json`. The
  binaries' prefix keeps them out of the ratchet's glob, and the step is
  `continue-on-error`.

Per core, planar / interleaved bracketed instructions (`bench.yml` run 173,
`workflow_dispatch` at `5634a73`, 2026-11-23: the Butterworth on
`biquad_cascade_block`; every pair bit-identical on every key; the counts
themselves are in each key's `layout-<key>.json`; the `hexagon` key counts
packets):

| key | `bw4` 2 | `bw4` 6 | `bw4` 12 | `bw8` 2 | `bw8` 6 | `bw8` 12 | `gain` 2 | `gain` 6 | `gain` 12 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `m4-softfp` | 1.003 | 1.003 | 1.003 | 0.999 | 0.999 | 0.999 | 1.018 | 1.041 | 1.050 |
| `m4f` | 1.004 | 1.001 | 1.001 | 0.948 | 0.949 | 0.948 | 1.012 | 1.301 | 1.403 |
| `m33` | 1.004 | 1.001 | 1.001 | 0.952 | 0.952 | 0.952 | 1.012 | 1.301 | 1.403 |
| `m55` | 1.005 | 1.004 | 1.005 | 0.950 | 0.951 | 0.951 | 1.248 | 0.276 | 0.291 |
| `m55-ooura` | 1.005 | 1.004 | 1.005 | 0.950 | 0.951 | 0.951 | 1.248 | 0.276 | 0.291 |
| `hexagon` | 1.284 | 1.287 | 1.288 | 1.000 | 1.002 | 1.002 | 1.008 | 1.013 | 1.014 |

Below 1 the mono stages are cheaper. On the block kernel the Butterworth's
two layouts are within 1 % of each other at order 4 on the M-class cores
(the interleaved stage runs its channels in pairs, the mono stages one at a
time) and the mono stages 5 % cheaper at order 8, where the interleaved
stage runs one channel at a time too. On
Hexagon interleaved is 1.28 times cheaper at order 4 (pairs fill its
packets) and the two are even at order 8. With the per-sample kernel (run
170, `f46a688`) the mono stages were 1 to 14 % cheaper on the M-class cores
with an FPU, and both layouts cost 1.9 to 3.5 times what they do now. The
stateless gain stage goes the other way where the vector unit is wide: on
the M55 planar is 3.5 times cheaper at 6 and 12 channels, on the scalar M4F
and M33 interleaved is, by 1.3 to 1.4 times.

The `set_cutoff()` scenarios on the keys that do not gate them (the same run;
four calls, order 8): `bw8_f32_set_cutoff` 136,672,436 on `m4-softfp`, 136,651,481 on `m4f`, 137,417,360 on `m33`;
`bw8_f32_set_cutoff_deviation` reads 16 more on each (today the deviation is
computed by `set_cutoff()` and the read returns it).

The bench came from the branch `ccr-ef01ff18-nxz8n4` (`0fd75b2`), where it
also measured an FFT FIR stage and its chain with the Butterworth;
`fft_fir.h` is not on `main`, so those two kinds are not built here.

## Provenance

`bench/call_marker.h`, the plugin's bracket feature and `bench/layout/` are
this repository's, the bracket mechanism ported from SampleRateTap's
(`tools/qemu_insn_plugin`, its `tools/probe/call_marker.h`; MIT, the
SampleRateTap contributors) with the marker names and the output line
renamed.


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
