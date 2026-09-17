#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright 2026 Timothy Place and the DspTap contributors.
# Adapted from MuTap's scripts/icount.py (MIT, MuTap contributors).
"""Deterministic instruction-count ratchet (see bench/README.md).

Runs every tap_dsp_icount_* binary in a build directory under QEMU with the
instruction-counting plugin (tools/qemu_insn_plugin), then compares against
bench/baselines.json.

  icount.py --target {m4-softfp,m4f,m33,m55,m55-ooura} --build-dir DIR
            --plugin LIB [--update] [--record FILE]
            [--baselines bench/baselines.json] [--tolerance 0.03]
  icount.py --merge FILE [FILE ...] [--baselines bench/baselines.json]

The gate is two-sided: exit nonzero if any scenario regresses beyond
tolerance, improves beyond tolerance (the baseline must be re-recorded so
the gate stays tight), has no recorded baseline, or has a recorded baseline
but no binary (a renamed or removed workload must not linger as a dead gate
entry). --update rewrites the target's entry to exactly the measured
scenarios instead. --record FILE always writes {target: measured} as JSON,
whatever the verdict, so a run that fails on NO BASELINE still hands back
the numbers to commit (bench.yml uploads it as measured-<key> on every run).
Seeding runs on pushes to main, never from a pull request (bench.yml); --merge
folds several per-target files (each with one target filled in) into one,
which is how the seeding commit is assembled from the job's artifacts.

The QEMU machine per target, the binary prefix and the output markers are
DspTap's; the gate logic is MuTap's.
"""
import argparse
import glob
import json
import os
import pathlib
import re
import subprocess
import sys

# Both M4 flavours (soft-float and M4F) run on the same board model; the
# difference is the toolchain's -mfloat-abi. m55 and m55-ooura likewise share
# the AN547: the key records which float32 backend the binaries were built
# with (CMSIS-DSP Helium, the deployed profile, or the Ooura fallback).
MACHINES = {
    "m4-softfp": "mps2-an386",
    "m4f": "mps2-an386",
    "m33": "mps2-an505",
    "m55": "mps3-an547",
    "m55-ooura": "mps3-an547",
}
PREFIX = "tap_dsp_icount_"
DONE_MARKER = "TAP_DSP_ICOUNT_DONE ok=1"
COUNT_RE = re.compile(r"TAP_DSP_INSN_COUNT (\d+)")


def qemu_cmd(target: str, plugin: str, binary: str) -> list[str]:
    # "-d plugin" routes qemu_plugin_outs() to stderr; without it the count
    # line is silently dropped.
    machine = MACHINES.get(target)
    if machine is None:
        raise SystemExit(f"unknown target {target}")
    return ["qemu-system-arm", "-M", machine, "-nographic", "-semihosting",
            "-d", "plugin", "-plugin", plugin, "-kernel", binary]


def measure(target: str, plugin: str, binary: str) -> int:
    try:
        proc = subprocess.run(qemu_cmd(target, plugin, binary), timeout=1200,
                              capture_output=True, text=True)
    except subprocess.TimeoutExpired:
        raise SystemExit(f"{binary}: timed out after 1200 s under QEMU")
    out = proc.stdout + proc.stderr
    if DONE_MARKER not in out:
        print(out, file=sys.stderr)
        raise SystemExit(f"{binary}: workload did not complete cleanly")
    m = COUNT_RE.search(out)
    if not m:
        print(out, file=sys.stderr)
        raise SystemExit(f"{binary}: no TAP_DSP_INSN_COUNT (plugin not loaded?)")
    return int(m.group(1))


def write(path: pathlib.Path, baselines: dict) -> None:
    path.write_text(json.dumps(baselines, indent=2, sort_keys=True) + "\n")


def merge(path: pathlib.Path, files: list[str]) -> int:
    baselines = json.loads(path.read_text()) if path.exists() else {}
    for f in files:
        for target, scenarios in json.loads(pathlib.Path(f).read_text()).items():
            if scenarios:
                baselines[target] = scenarios
                print(f"{target}: {len(scenarios)} scenario(s) from {f}")
    write(path, baselines)
    print(f"merged into {path}")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", choices=sorted(MACHINES))
    ap.add_argument("--build-dir")
    ap.add_argument("--plugin")
    ap.add_argument("--baselines", default="bench/baselines.json")
    ap.add_argument("--tolerance", type=float, default=0.03)
    ap.add_argument("--update", action="store_true")
    ap.add_argument("--record", metavar="FILE",
                    help="always write {target: measured} here, whatever the verdict")
    ap.add_argument("--merge", nargs="+", metavar="FILE")
    args = ap.parse_args()

    path = pathlib.Path(args.baselines)
    if args.merge:
        return merge(path, args.merge)
    if not (args.target and args.build_dir and args.plugin):
        ap.error("--target, --build-dir and --plugin are required (or --merge)")

    binaries = sorted(glob.glob(os.path.join(args.build_dir, "**", PREFIX + "*"),
                                recursive=True))
    binaries = [b for b in binaries if os.access(b, os.X_OK) and os.path.isfile(b)]
    if not binaries:
        raise SystemExit(f"no {PREFIX}* binaries under {args.build_dir}")

    baselines = json.loads(path.read_text()) if path.exists() else {}
    base = baselines.get(args.target, {})

    failures = []
    measured = {}
    for binary in binaries:
        scenario = os.path.basename(binary).removeprefix(PREFIX)
        count = measure(args.target, args.plugin, binary)
        measured[scenario] = count
        recorded = base.get(scenario)
        if recorded is None:
            print(f"{scenario}: {count} insns (NO BASELINE — commit this value)")
            if not args.update:
                failures.append(scenario)
        elif recorded == 0:
            print(f"{scenario}: {count} insns vs baseline 0 (INVALID BASELINE)")
            failures.append(scenario)
        else:
            delta = (count - recorded) / recorded
            verdict = "ok"
            if delta > args.tolerance:
                verdict = "REGRESSION"
                failures.append(scenario)
            elif delta < -args.tolerance:
                # Two-sided: a stale (too-high) baseline would let future
                # regressions hide inside the slack, so improvements must be
                # committed too.
                verdict = ("IMPROVED beyond tolerance — run icount.py --update "
                           "and commit bench/baselines.json")
                failures.append(scenario)
            print(f"{scenario}: {count} insns vs baseline {recorded} "
                  f"({delta:+.2%}) {verdict}")

    # A recorded scenario with no binary is a dead gate entry (renamed or
    # removed workload); compare mode fails on it, --update drops it.
    for scenario in sorted(set(base) - set(measured)):
        print(f"{scenario}: baseline {base[scenario]} but no binary "
              "(STALE BASELINE — run icount.py --update and commit)")
        if not args.update:
            failures.append(scenario)

    if args.record:
        pathlib.Path(args.record).write_text(
            json.dumps({args.target: measured}, indent=2, sort_keys=True) + "\n")
        print(f"recorded {args.record}")

    if args.update:
        # Exactly the measured scenarios: stale keys for renamed/removed
        # workloads must not linger as dead gate entries.
        baselines[args.target] = measured
        write(path, baselines)
        print(f"updated {path}")
        return 0
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
