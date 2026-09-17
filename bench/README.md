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

On the host the `tap_dsp_icount_*` binaries are a smoke test: they must run,
and two runs of the same binary must print the same checksum. There is no
instruction count on the host; the count comes from QEMU on the legs.

## Scenarios

One binary per scenario, selected by compile definitions (bare metal has no
argv), all built from `icount/icount_main.cpp`:

| scenario | profile | N | targets |
|---|---|---|---|
| `rfft_f32_512` | float | 512 | all legs |
| `rfft_f32_2048` | float | 2048 | all legs |
| `rfft_f64_512` | double | 512 | host-class only: soft-float double is not a profile, so the bare-metal legs neither build nor baseline it |

Stage 3b adds `rfft_q15_512`, `rfft_q31_512` and `rfft_q31_2048` (Part 11).

Each scenario constructs one `tap::dsp::basic_real_fft<Sample>`, then runs a
forward + inverse loop over a four-block xorshift corpus for 2^20 samples per
direction (2048 iterations at N = 512, 512 at N = 2048), folding every output
into a checksum of the scenario's own sample type, printed at the end as its
bit pattern:

```
TAP_DSP_ICOUNT_DONE ok=1 engine=reference_c scenario=rfft_f32_512 checksum=0x446b9819
```

Why these numbers: on the host (x86-64, GCC 13 `-O2`, callgrind) the three
scenarios execute 98 M, 114 M and 102 M instructions, and construction plus
the one-time Ooura table build plus the print is under 0.2 M — well under the
1 % the design asks for; the Arm counts are larger but the ratio holds. No
`<random>` (a toolchain's libstdc++ would move the count), no allocation in
the loop, and no double anywhere in the float scenarios (the M4 soft-float
leg would otherwise measure libgcc).

`TAP_DSP_BENCH_ENGINE` (`bench_common.h`, a CMake cache variable of the same
name) selects the engine: `reference_c` today — the class as built, which is
the vendored Ooura C, or on the `m55` key the CMSIS-DSP Helium backend behind
the same class. `split_radix` is the Stage 2a port; from then until Stage 2c
each float scenario builds twice, C and port, and the job prints both counts
and their ratio. That ratio is the Stage 2b gate.

## Baselines and targets

[`baselines.json`](baselines.json) is keyed by target, one dict per key:

| key | core | QEMU machine | float32 FFT |
|---|---|---|---|
| `m4-softfp` | Cortex-M4, `-mfloat-abi=soft` | `mps2-an386` | Ooura (every float op a library call) |
| `m4f` | Cortex-M4F, `fpv4-sp-d16` | `mps2-an386` | Ooura |
| `m33` | Cortex-M33, no MVE | `mps2-an505` | Ooura |
| `m55` | Cortex-M55, Helium | `mps3-an547` | CMSIS-DSP — the deployed profile |
| `m55-ooura` | Cortex-M55, `-DTAP_DSP_FFT_CMSIS=OFF` | `mps3-an547` | Ooura — the fallback |

JSON carries no comments, so the provenance of every recorded set lives here,
in the table below: the commit that recorded it, and the GCC and QEMU
versions, because the counts are only comparable within a toolchain/QEMU
pair. The plugin header is pinned to the QEMU that Ubuntu 24.04 ships
(8.2.2, plugin API v2) and digest-verified on download, as MuTap does.

### Recorded baselines

| key | recorded at (commit) | arm-none-eabi-gcc | qemu-system-arm | note |
|---|---|---|---|---|
| all | — | — | — | **unseeded**: the skeleton is empty by design; numbers are measured on CI, never typed in |

### Seeding, and how the job decides what to do

The workflow reads `baselines.json` per key. A key whose dict is **empty** is
unseeded: its job runs `scripts/icount.py --update`, prints every scenario's
count, writes it to the job summary and uploads the file as the artifact
`baselines-<key>`; the `seed-summary` job then merges the per-key artifacts
into one `baselines-merged` artifact. The seeding commit copies that file to
`bench/baselines.json`, fills the table above with the run's commit and the
versions the job printed in its "Toolchain versions" step, and nothing else
rides in it. From that commit on, the key's job **compares** and a red
ratchet is a failing check.

The first run of `bench.yml` after the QEMU legs are green is therefore a
seeding run for all five keys (wave 2 of Part 12).

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
    --plugin /tmp/libinsncount.so [--update]
python3 scripts/icount.py --merge a.json b.json    # fold per-key files into one
```

## Policy

- **Tolerance ±3 %**, the family's number. Tighter is not better: link-order
  and libgcc changes move counts by fractions of a percent legitimately.
- **Two-sided.** A regression beyond tolerance fails; an *improvement* beyond
  tolerance also fails, so a stale, too-high baseline can never let a later
  regression hide in the slack — the winning commit re-records.
- **The ratchet runs on push and pull request** on every QEMU leg. A red
  ratchet is a failing check, not a warning.
- **`--update` is a written commit on its own**: the measured before/after per
  key and the reason go into the table above. An expected regression (a
  correctness fix that costs instructions) is written down; there is no
  silent absorb. A key is never emptied to "reset" it — that is an update
  without its reason.
- **Size runs in the same job.** Today the job prints `arm-none-eabi-size`
  for every scenario binary (report only: a whole binary carries newlib and
  semihosting). The `.text` ceilings of Part 10 item 4 — asserted per key on
  a size-probe object per profile — become numbers in `bench.yml` once the
  legs are seeded, and are recorded the same way as the counts.
- **Wall clock is never a gate.** `bench_fft` is the local tool for the
  desktop and Apple vDSP claims; its numbers go into `docs/fft-design.md`
  with machine and date.

## What the ratchet catches, by stage

| stage | what a red ratchet means |
|---|---|
| 2b (routing flip to the port) | a port that vectorizes worse than the C; the C-vs-port ratio printed from 2a on is the gate |
| 3b (fixed point) | nothing yet — it seeds; every later `SMMULR`, Helium or table-layout change to the Q15/Q31 kernel then has a number to beat |
| 4 (engine as a parameter) | a backend regression on the deployed `m55` profile, or the `m55-ooura` fallback quietly getting slower |
| 6 (hygiene) | the hygiene pass slowing the hot path unnoticed |

## Provenance

`icount/icount_main.cpp`, `scripts/icount.py` and
`tools/qemu_insn_plugin/insn_count.c` are adapted from MuTap's (MIT, MuTap
contributors; the same pattern lives in SampleRateTap and RatioTap). Copied
rather than shared because they are small; a taphouse-style consolidation is
the eventual home (Part 11, "Sharing").
