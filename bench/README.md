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

| scenario | profile | N | engine | targets | gated |
|---|---|---|---|---|---|
| `rfft_f32_512` | float | 512 | `TAP_DSP_BENCH_ENGINE` (`basic_real_fft`, what ships: the split-radix engine since Stage 2b; CMSIS behind it on `m55`) | all legs | yes |
| `rfft_f32_2048` | float | 2048 | as above | all legs | yes |
| `rfft_f64_512` | double | 512 | as above | host-class only: soft-float double is not a profile, so the bare-metal legs neither build nor baseline it | no key counts it |
| `rfft_f32_512_c` | float | 512 | `reference_c`, the vendored Ooura C called directly | all legs | **no** — informational until Stage 2c retires the C |
| `rfft_f32_2048_c` | float | 2048 | `reference_c` | all legs | **no** — informational until Stage 2c |

Stage 3b's `rfft_q15_512`, `rfft_q31_512` and `rfft_q31_2048` (Part 11) are
a follow-up after the Stage 2b flip merges.

The `_c` scenarios exist from Stage 2b until Stage 2c retires the C: each
float scenario builds twice (`bench/icount/CMakeLists.txt`), and
`scripts/icount.py` treats the `_c` suffix (`INFORMATIONAL_SUFFIX`) as
informational — the binary is counted and printed after the gated
scenarios, beside its shipping sibling, with the ratio C/shipping and whether
the two `DONE` checksums agree — but it never enters the verdict: a `_c`
binary that times out, faults or does not print `ok=1` is reported as
`informational binary failed: <reason>` and the run continues to the gated
scenarios' verdict; `--update` never writes it to `baselines.json`; a
baseline that names it is reported and ignored; and `--record` files it
under a separate top-level `informational` key that `--merge` skips
(`--merge` also drops a `_c` key found under a real target, so a hand-edited
file cannot seed one). That is what keeps the C-vs-port comparison honest in
every job log until the C is gone: the gated bare key measures what ships,
and the `_c` sibling is the C it replaced, with the checksum line saying
whether the two are still bit-identical at the bench's default flags. On the
`m55` key the shipping sibling is the CMSIS-DSP Helium backend, so the ratio
there is C-vs-CMSIS; the C-vs-port ratio on the M55 is the `m55-ooura` key's.
(From Stage 2a until 2b the same mechanism ran the other way round, as
`_port`: the port beside the shipping C. The flip re-recorded the bare keys to
the port's counts — see the table below — and swapped the sides.)

Each scenario constructs one transform (`tap::dsp::basic_real_fft<Sample>`,
or for the `_c` scenarios the vendored C behind an adapter presenting the same
out-of-place surface with the same copy loop and 2/N arithmetic,
`bench_common.h`), then runs a forward + inverse loop over a four-block xorshift corpus for 2^20 samples per
direction (2048 iterations at N = 512, 512 at N = 2048). Every output word
of every iteration — the spectrum and the scaled inverse — goes through an
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
is an integer hash: it is the fingerprint the job compares between the C
and the port on the QEMU legs, where the host parity TU does not run (the
checksum above is the same value the C printed before the flip: the port is
bit-identical). `ok` is a sanity check that the last iteration round-trips
its input within 1e-3 in the scenario's own precision; `backend=` names what
the build routed `basic_real_fft<float>` through (`split_radix`, `cmsis` on
the `m55` key, `accelerate` on Apple), or `ooura_c` for the `_c` scenarios,
independently of the engine value.

Why these numbers: on the host (x86-64, GCC 13 `-O2`, callgrind) the three
scenarios execute 109 M, 125 M and 108 M instructions.
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

`TAP_DSP_BENCH_ENGINE` (`bench_common.h`, a CMake cache variable of the same
name) selects the engine the bare-keyed binaries measure: `basic_real_fft`
(the default) — the class as built, which is the split-radix engine
(`include/tap/dsp/fft/split_radix.h`) or, on the `m55` key, the CMSIS-DSP
Helium backend behind the same class — or `reference_c`, the vendored Ooura C
called directly (`rdft` / `rdft_f`), printed as `engine=reference_c
backend=ooura_c`. At the default the `_c` binaries above (and the `_c` size
probe) are also built, against the C, so every bench run from Stage 2b until
2c reports both counts and their ratio; when the variable is already
`reference_c` the bare-keyed binaries *are* the C and the `_c` pairs are not
built (`bench/icount/CMakeLists.txt`, `bench/CMakeLists.txt`: the C compared
to itself measures nothing). A build directory configured before Stage 2b
still carries `reference_c` in its cache (the old default, which then meant
the class as built): reconfigure with `-UTAP_DSP_BENCH_ENGINE` or a fresh
directory, or the bare keys silently measure the C. Before the flip the
port/C ratio these pairs printed was what Stage 2b was judged on (the port
within ±3 % of the C on every QEMU leg, audit Part 3); the outcome is in the
table below.

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

| m33 | rfft_f32_2048 | 116,385,409 | 115,465,626 | -0.79 % | **re-recorded at the Stage 2b flip** (tap/DspTap#31): the bare key now measures the split-radix engine that ships through `basic_real_fft`; the C it replaced is the `_c` informational sibling of the same run (output checksums identical on every key). Per Part 11 / D11 the exact new count is recorded, never absorbed: m55-ooura sat outside the +-3% band (IMPROVED), the M4 and M33 keys inside it, and all four move together | [run 35844483811](https://github.com/tap/DspTap/actions/runs/35844483811) (workflow_dispatch on the PR branch at `b078a80`; the PR's compare-mode run confirms +0.00 %) | the #31 squash SHA on `main` | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m33 | rfft_f32_512 | 102,248,169 | 100,945,841 | -1.27 % | **re-recorded at the Stage 2b flip** (tap/DspTap#31): the bare key now measures the split-radix engine that ships through `basic_real_fft`; the C it replaced is the `_c` informational sibling of the same run (output checksums identical on every key). Per Part 11 / D11 the exact new count is recorded, never absorbed: m55-ooura sat outside the +-3% band (IMPROVED), the M4 and M33 keys inside it, and all four move together | [run 35844483811](https://github.com/tap/DspTap/actions/runs/35844483811) (workflow_dispatch on the PR branch at `b078a80`; the PR's compare-mode run confirms +0.00 %) | the #31 squash SHA on `main` | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m4-softfp | rfft_f32_2048 | 2,296,984,479 | 2,294,360,259 | -0.11 % | **re-recorded at the Stage 2b flip** (tap/DspTap#31): the bare key now measures the split-radix engine that ships through `basic_real_fft`; the C it replaced is the `_c` informational sibling of the same run (output checksums identical on every key). Per Part 11 / D11 the exact new count is recorded, never absorbed: m55-ooura sat outside the +-3% band (IMPROVED), the M4 and M33 keys inside it, and all four move together | [run 35844483811](https://github.com/tap/DspTap/actions/runs/35844483811) (workflow_dispatch on the PR branch at `b078a80`; the PR's compare-mode run confirms +0.00 %) | the #31 squash SHA on `main` | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m4-softfp | rfft_f32_512 | 1,868,441,244 | 1,864,929,141 | -0.19 % | **re-recorded at the Stage 2b flip** (tap/DspTap#31): the bare key now measures the split-radix engine that ships through `basic_real_fft`; the C it replaced is the `_c` informational sibling of the same run (output checksums identical on every key). Per Part 11 / D11 the exact new count is recorded, never absorbed: m55-ooura sat outside the +-3% band (IMPROVED), the M4 and M33 keys inside it, and all four move together | [run 35844483811](https://github.com/tap/DspTap/actions/runs/35844483811) (workflow_dispatch on the PR branch at `b078a80`; the PR's compare-mode run confirms +0.00 %) | the #31 squash SHA on `main` | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m4f | rfft_f32_2048 | 111,859,257 | 111,278,416 | -0.52 % | **re-recorded at the Stage 2b flip** (tap/DspTap#31): the bare key now measures the split-radix engine that ships through `basic_real_fft`; the C it replaced is the `_c` informational sibling of the same run (output checksums identical on every key). Per Part 11 / D11 the exact new count is recorded, never absorbed: m55-ooura sat outside the +-3% band (IMPROVED), the M4 and M33 keys inside it, and all four move together | [run 35844483811](https://github.com/tap/DspTap/actions/runs/35844483811) (workflow_dispatch on the PR branch at `b078a80`; the PR's compare-mode run confirms +0.00 %) | the #31 squash SHA on `main` | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m4f | rfft_f32_512 | 98,090,666 | 97,138,544 | -0.97 % | **re-recorded at the Stage 2b flip** (tap/DspTap#31): the bare key now measures the split-radix engine that ships through `basic_real_fft`; the C it replaced is the `_c` informational sibling of the same run (output checksums identical on every key). Per Part 11 / D11 the exact new count is recorded, never absorbed: m55-ooura sat outside the +-3% band (IMPROVED), the M4 and M33 keys inside it, and all four move together | [run 35844483811](https://github.com/tap/DspTap/actions/runs/35844483811) (workflow_dispatch on the PR branch at `b078a80`; the PR's compare-mode run confirms +0.00 %) | the #31 squash SHA on `main` | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m55-ooura | rfft_f32_2048 | 107,806,480 | 102,784,096 | -4.66 % | **re-recorded at the Stage 2b flip** (tap/DspTap#31): the bare key now measures the split-radix engine that ships through `basic_real_fft`; the C it replaced is the `_c` informational sibling of the same run (output checksums identical on every key). Per Part 11 / D11 the exact new count is recorded, never absorbed: m55-ooura sat outside the +-3% band (IMPROVED), the M4 and M33 keys inside it, and all four move together | [run 35844483811](https://github.com/tap/DspTap/actions/runs/35844483811) (workflow_dispatch on the PR branch at `b078a80`; the PR's compare-mode run confirms +0.00 %) | the #31 squash SHA on `main` | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m55-ooura | rfft_f32_512 | 94,561,954 | 89,276,321 | -5.59 % | **re-recorded at the Stage 2b flip** (tap/DspTap#31): the bare key now measures the split-radix engine that ships through `basic_real_fft`; the C it replaced is the `_c` informational sibling of the same run (output checksums identical on every key). Per Part 11 / D11 the exact new count is recorded, never absorbed: m55-ooura sat outside the +-3% band (IMPROVED), the M4 and M33 keys inside it, and all four move together | [run 35844483811](https://github.com/tap/DspTap/actions/runs/35844483811) (workflow_dispatch on the PR branch at `b078a80`; the PR's compare-mode run confirms +0.00 %) | the #31 squash SHA on `main` | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m55 | rfft_f32_2048 | 54,858,120 | 54,858,120 | +0.00 % | **not re-recorded**: CMSIS-DSP Helium is untouched by the flip; the same run measured 54,858,195 (+75 instructions, +0.00 %: the wrapper lost its two unused Ooura workspace vectors and the DONE line's `engine=` string grew by three characters) | run 35844483811, compare mode | — | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |
| m55 | rfft_f32_512 | 52,382,331 | 52,382,331 | +0.00 % | as above; measured 52,382,366 (+35 instructions, +0.00 %) | run 35844483811, compare mode | — | 13.2.1 (15:13.2.rel1-2) | 8.2.2 (1:8.2.2+ds-0ubuntu1.18) |

`before` is `—` for a seed. `main SHA` is the commit on `main` whose push
run measured the numbers: a pull-request head SHA stops resolving after
this repo's rebase/squash + branch-delete flow, so a seed or update is never
taken from a PR run. The one exception is a re-record forced by a routing
change, where the numbers cannot be measured on `main` before the change is
on `main`: the Stage 2b rows above were seeded by a `workflow_dispatch` on
the PR branch (the branch SHA is named; it stops resolving after the squash)
and the PR's own compare-mode run is what proves the committed numbers are
the branch's, at +0.00 %. The `main SHA` column then names the squash.

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
job printed in its "Toolchain versions" step, sets the `text_ceiling`
numbers in `bench.yml` from the same run's size step, and nothing else rides
in it. From that commit on, the key's job compares and a red ratchet is a
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
- **Size runs in the same job, MinSizeRel.** The job builds
  `tap_dsp_size_probe_rfft_f32_512` (`bench/size_probe.cpp`: one profile,
  no stdio) a second time with `-DCMAKE_BUILD_TYPE=MinSizeRel` — the build
  type the test legs and Part 10 item 4 use, so a ceiling here and a size in
  `ci.yml` describe the same object — and reads the `.text` row of
  `arm-none-eabi-size -A` (Berkeley format folds `.rodata` and NOLOAD
  sections into "text"). Each key carries a `text_ceiling` in `bench.yml`;
  `0` means not yet recorded and the step only prints. The ceilings are a
  **promised item for wave 2**: the seeding commit sets them from the same
  `main` run that seeds the counts, and they are updated the same way as the
  counts, in the table above. From Stage 2b until 2c the step also builds
  `tap_dsp_size_probe_rfft_f32_512_c` (the vendored C) and prints its `.text`
  and the shipping/C ratio beside the shipping probe's, informational: the
  ceiling applies to the shipping probe, and the C probe's numbers go into
  `docs/fft-design.md` from the job log. The step runs with `pipefail`
  and fails on a `size` that fails or prints no numeric `.text` row, for
  either probe: a measurement that cannot be taken is an error, not a
  `0.0000` ratio.
- **Wall clock is never a gate.** `bench_fft` is the local tool for the
  desktop and Apple vDSP claims; its numbers go into `docs/fft-design.md`
  with machine and date.

## What the ratchet catches, by stage

| stage | what a red ratchet means |
|---|---|
| 2a (port beside the C, nothing routed) | nothing: the `_port` scenarios were informational, and the bare keys measured the unchanged C at 0 % |
| 2b (routing flip to the port; done) | a port that vectorizes worse than the C; the C-vs-port ratio printed from 2a on was the gate, and the bare keys were re-recorded to the port's exact counts (table above) |
| 2b → 2c (the C as the `_c` informational sibling) | nothing gated; the job log keeps showing C/shipping and whether the checksums still agree |
| 3b (fixed point) | nothing yet — it seeds; every later `SMMULR`, Helium or table-layout change to the Q15/Q31 kernel then has a number to beat |
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
