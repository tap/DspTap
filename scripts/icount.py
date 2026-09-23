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

Informational scenarios. A scenario whose key ends in INFORMATIONAL_SUFFIX
("_c": the vendored Ooura C called directly, built beside the shipping
split-radix engine from Stage 2b until Stage 2c retires the C) is counted and
printed with its ratio to the sibling scenario (the key without the suffix)
and whether the two output checksums agree, but is never a gate entry: it
never enters the verdict (a `_c` binary that times out, faults or does not
print ok=1 is reported as "informational binary failed: <reason>" and the run
continues to the gated scenarios' verdict, where only a gated binary's
failure aborts the run), --update never writes it to the baselines, a
baseline that names it is reported and ignored, and --record files it under
a separate top-level "informational" key that --merge skips (--merge also
drops a `_c` key found under a real target, so a hand-edited record file
cannot seed one). This is what keeps the C-vs-port comparison honest in every
job log until 2c deletes the C: the gated bare key measures what ships, and
the `_c` sibling is the C it replaced. (Before 2b the same mechanism ran the
other way round, as `_port`: the port beside the shipping C.)

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
DONE_RE = re.compile(r"TAP_DSP_ICOUNT_DONE ok=1 (.*)")
# Two characters, tested with endswith(): a future GATED key that happened to
# end in "_c" would be classified informational. None does (the gated keys
# are rfft_<precision>_<n>); the suffix goes with the C at Stage 2c.
INFORMATIONAL_SUFFIX = "_c"
INFORMATIONAL_KEY = "informational"


def qemu_cmd(target: str, plugin: str, binary: str) -> list[str]:
    # "-d plugin" routes qemu_plugin_outs() to stderr; without it the count
    # line is silently dropped.
    machine = MACHINES.get(target)
    if machine is None:
        raise SystemExit(f"unknown target {target}")
    return ["qemu-system-arm", "-M", machine, "-nographic", "-semihosting",
            "-d", "plugin", "-plugin", plugin, "-kernel", binary]


class MeasurementError(Exception):
    """A binary that did not yield a count: QEMU timed out, the workload did
    not print the DONE marker (fault, hang caught by the timeout, ok=0), or
    the plugin's count line is missing. Fatal for a gated scenario; reported
    and skipped for an informational one (see the module docstring)."""


def measure(target: str, plugin: str, binary: str) -> tuple[int, dict[str, str]]:
    """Returns the instruction count and the DONE line's key=value fields
    (engine, backend, scenario, checksum). Raises MeasurementError when the
    binary produced no count; the QEMU output is echoed to stderr first."""
    try:
        proc = subprocess.run(qemu_cmd(target, plugin, binary), timeout=1200,
                              capture_output=True, text=True)
    except subprocess.TimeoutExpired:
        raise MeasurementError(f"{binary}: timed out after 1200 s under QEMU")
    out = proc.stdout + proc.stderr
    if DONE_MARKER not in out:
        print(out, file=sys.stderr)
        raise MeasurementError(f"{binary}: workload did not complete cleanly")
    m = COUNT_RE.search(out)
    if not m:
        print(out, file=sys.stderr)
        raise MeasurementError(f"{binary}: no TAP_DSP_INSN_COUNT (plugin not loaded?)")
    done = DONE_RE.search(out)
    fields = dict(kv.split("=", 1) for kv in done.group(1).split() if "=" in kv) if done else {}
    return int(m.group(1)), fields


def write(path: pathlib.Path, baselines: dict) -> None:
    path.write_text(json.dumps(baselines, indent=2, sort_keys=True) + "\n")


def merge(path: pathlib.Path, files: list[str]) -> int:
    baselines = json.loads(path.read_text()) if path.exists() else {}
    for f in files:
        for target, scenarios in json.loads(pathlib.Path(f).read_text()).items():
            if target == INFORMATIONAL_KEY:
                continue  # never a gate entry (see the module docstring)
            # Nor is a `_c` key nested under a real target: --update and
            # --record never write one there, so it can only come from a
            # hand-edited file, and it is dropped rather than seeded.
            gated = {k: v for k, v in scenarios.items() if not k.endswith(INFORMATIONAL_SUFFIX)}
            for dropped in sorted(set(scenarios) - set(gated)):
                print(f"{target}: {dropped} from {f} DROPPED (informational scenario; never a gate entry)")
            if gated:
                baselines[target] = gated
                print(f"{target}: {len(gated)} scenario(s) from {f}")
    write(path, baselines)
    print(f"merged into {path}")
    return 0


def describe(fields: dict[str, str]) -> str:
    return f"[engine={fields.get('engine', '?')} backend={fields.get('backend', '?')}]"


def report_informational(informational: dict, failed: dict, measured: dict, fields: dict) -> None:
    if not informational and not failed:
        return
    print(f"--- informational: '{INFORMATIONAL_SUFFIX}' scenarios (the vendored C beside what ships; "
          "counted, never gated, never baselined; Stage 2b until 2c) ---")
    for scenario, reason in sorted(failed.items()):
        print(f"{scenario}: informational binary failed: {reason}")
    for scenario, count in sorted(informational.items()):
        sibling = scenario[: -len(INFORMATIONAL_SUFFIX)]
        line = f"{scenario}: {count} insns {describe(fields[scenario])}"
        if sibling in measured:
            base = measured[sibling]
            same = fields[scenario].get("checksum") == fields[sibling].get("checksum")
            line += (f"; {sibling}: {base} insns {describe(fields[sibling])}"
                     f"; ratio {scenario}/{sibling} = {count / base:.4f}"
                     f"; output checksums {'identical' if same else 'DIFFER'}")
        else:
            line += f"; no sibling {sibling} to compare against"
        print(line)


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
    informational = {}
    informational_failed = {}
    fields = {}
    for binary in binaries:
        scenario = os.path.basename(binary).removeprefix(PREFIX)
        try:
            count, fields[scenario] = measure(args.target, args.plugin, binary)
        except MeasurementError as e:
            if not scenario.endswith(INFORMATIONAL_SUFFIX):
                raise SystemExit(str(e))  # a gated binary that yields no count aborts the run
            # An informational binary that yields no count is reported (in
            # the informational block below, so it lands in the job log and
            # the step summary) and never touches the verdict.
            informational_failed[scenario] = str(e)
            continue
        if scenario.endswith(INFORMATIONAL_SUFFIX):
            informational[scenario] = count
            if scenario in base:
                print(f"{scenario}: baseline {base[scenario]} IGNORED (informational scenario; "
                      "remove it from bench/baselines.json)")
            continue
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
    # removed workload); compare mode fails on it, --update drops it. An
    # informational key in the baselines is not a gate entry either way.
    for scenario in sorted(set(base) - set(measured)):
        if scenario.endswith(INFORMATIONAL_SUFFIX):
            continue
        print(f"{scenario}: baseline {base[scenario]} but no binary "
              "(STALE BASELINE — run icount.py --update and commit)")
        if not args.update:
            failures.append(scenario)

    report_informational(informational, informational_failed, measured, fields)

    if args.record:
        record = {args.target: measured}
        if informational:
            record[INFORMATIONAL_KEY] = {args.target: informational}
        pathlib.Path(args.record).write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
        print(f"recorded {args.record}")

    if args.update:
        # Exactly the measured (gated) scenarios: stale keys for renamed or
        # removed workloads must not linger as dead gate entries, and the
        # informational scenarios never become gate entries.
        baselines[args.target] = measured
        write(path, baselines)
        print(f"updated {path}")
        return 0
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
