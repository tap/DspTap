#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright 2026 Timothy Place and the DspTap contributors.
"""Generates include/tap/dsp/detail/iir_tables_data.h: the committed recursive
designs iir_tables.h declares (SampleRateTap's docs/IIR_PLAN.md, I4 and I13).

    python3 tools/iir/gen_iir_tables.py     # from the repository root, in the
                                            # venv tools/iir/requirements.txt pins

Versions the committed header was generated with: Python 3.11.15, numpy 2.4.6,
scipy 1.17.1, mpmath 1.4.1 (tools/iir/requirements.txt). The optimizers'
results can move in their last digits under another scipy; the tests pin what
each row promises (the spec it meets), not the digits.

The rows, at the four spec points (iir_design.SPEC_POINTS):
- the near-linear two-path half-band at the divisors 1, 2, 4 and 8 (I13), by
  IIR-0's method: the lowest allpass order N meeting -(A + 1 dB) on the
  16384-point stopband grid, its delay K searched from N - 1 to N + 2, the
  iterated linear program then the SLSQP polish (iir_design.nearlinear_*);
- the four-path stage at design divisor 2 for the three 70 dB points, the
  placements IIR-0 found to pass I12's bar on instructions: six sections over
  the distinct permutations of the layout [2, 2, 1, 1], 24 SLSQP starts each
  (seed 1), the best on the folding bands, which must meet -(A + 1 dB);
- the Chebyshev II's L-inf gains at the odd orders 3 to 11, computed at 40
  digits from the sections the constexpr designer builds (unity DC gain per
  section, ordered by increasing pole radius): each partial cascade's peak on
  a 2048-point grid over [0, 1/2], refined by golden-section search;
- the type II FIR guard band filter matching each Chebyshev order's stop
  edge, by scipy.signal.remez (the shortest even length flat within 0.01 dB
  peak to peak to p and A dB down from the edge), where it converges.

Run time: a few minutes for the four-path searches, about twenty for the
near-linear ones (N19 at 120 dB, 5/12 dominates).
"""
from __future__ import annotations

import pathlib
import sys

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
import iir_design as d  # noqa: E402

OUT = ROOT / "include" / "tap" / "dsp" / "detail" / "iir_tables_data.h"


def fmt(x):
    s = f"{float(x):.17g}"
    if "." not in s and "e" not in s and "inf" not in s and "nan" not in s:
        s += ".0"
    return s


def arr_lines(xs, indent, per_line=4):
    xs = [fmt(x) for x in xs]
    lines = []
    for i in range(0, len(xs), per_line):
        lines.append(indent + ", ".join(xs[i:i + per_line]) + ("," if i + per_line < len(xs) else ""))
    return lines


def tag(A, pn, pd):
    return f"{A}_{pn}_{pd}"


def main():
    body = []
    nl_rows, fp_rows, g_rows, fir_rows = [], [], [], []
    print("near-linear half-bands", flush=True)
    for A, pn, pd in d.SPEC_POINTS:
        p = d.spec_p(pn, pd)
        for div in d.DIVISORS:
            N, K, den = d.nearlinear_design(A, p, div)
            sos = d.nearlinear_sos(den)
            fp_hi = p / (2.0 * div)
            sb = float(d.db(d.nearlinear_response_sos(sos, K, np.linspace(0.5 - fp_hi, 0.5, d.GRID))).max())
            tau = d.nearlinear_group_delay(den, K, np.linspace(0.0, fp_hi, d.GRID))
            dev = float(np.abs(tau - 2 * K).max())
            radius = float(np.max(np.abs(np.roots(den))))
            flat = [c for s in sos for c in s]
            name = f"k_near_linear_{tag(A, pn, pd)}_{div}"
            body.append(f"    inline constexpr double {name}[] = {{")
            body += arr_lines(flat, "        ")
            body.append("    };")
            nl_rows.append(f"        {{{A}, {pn}, {pd}, {div}, {N}, {K}, {name}, {fmt(sb)}, {fmt(dev)}, {fmt(radius)}}},")
            print(f"  {A} dB {pn}/{pd} d{div}: N{N} K{K} stopband {sb:.2f} dB deviation {dev:.4f}", flush=True)

    print("four-path stages", flush=True)
    for A, pn, pd in d.SPEC_POINTS:
        if A != 70:
            continue
        p = d.spec_p(pn, pd)
        D = 2
        sb, layout, paths = d.fourpath_design(A, p, D, d.distinct_permutations((2, 2, 1, 1)))
        if sb > -(A + d.MARGIN_DB):
            raise RuntimeError(f"four-path {A} dB {pn}/{pd} D{D}: {sb:.2f} dB misses the spec")
        names = []
        for k, pth in enumerate(paths):
            name = f"k_fourpath_{tag(A, pn, pd)}_{D}_p{k}"
            body.append(f"    inline constexpr double {name}[] = {{")
            body += arr_lines(pth, "        ")
            body.append("    };")
            names.append(name)
        fp_rows.append(f"        {{{A}, {pn}, {pd}, {D}, {{{', '.join(names)}}}, {fmt(sb)}}},")
        print(f"  {A} dB {pn}/{pd} D{D}: layout {layout} stopband {sb:.2f} dB", flush=True)

    print("Chebyshev II gains and type II FIRs", flush=True)
    for A, pn, pd in d.SPEC_POINTS:
        p = d.spec_p(pn, pd)
        for order in d.CHEBY2_ORDERS:
            secs, edge = d.cheby2_mp(order, p, A)
            gains = d.cheby2_linf_gains(secs)
            padded = [float(g) for g in gains] + [0.0] * (7 - len(gains))
            g_rows.append(f"        {{{A}, {pn}, {pd}, {order}, {{{', '.join(fmt(g) for g in padded)}}}}},")
            fir = d.typeii_min_length(p, float(edge), A)
            if fir is None:
                print(f"  {A} dB {pn}/{pd} order {order}: gains {[f'{float(g):.6f}' for g in gains]}; "
                      f"no type II FIR (remez did not converge)", flush=True)
                continue
            n, h, ripple, sbf = fir
            name = f"k_guard_fir_{tag(A, pn, pd)}_{order}"
            body.append(f"    inline constexpr double {name}[] = {{")
            body += arr_lines(h, "        ")
            body.append("    };")
            fir_rows.append(f"        {{{A}, {pn}, {pd}, {order}, {fmt(float(edge))}, {name}, {fmt(ripple)}, "
                            f"{fmt(sbf)}}},")
            print(f"  {A} dB {pn}/{pd} order {order}: gains {[f'{float(g):.6f}' for g in gains]}; "
                  f"type II {n} taps, ripple {ripple:.4f} dB, stopband {sbf:.2f} dB", flush=True)

    lines = [
        "/// @file iir_tables_data.h",
        "/// @brief Generated by tools/iir/gen_iir_tables.py: the committed recursive designs. Do not edit.",
        "// SPDX-License-Identifier: MIT",
        "// Copyright 2026 Timothy Place and the DspTap contributors.",
        "// NOLINTBEGIN(readability-identifier-naming)",
        "// Generated data, four values to a line; clang-format would re-pack the",
        "// tables, so the file opts out rather than changing what the generator",
        "// writes. The row types and what each field promises: tap/dsp/iir_tables.h.",
        "// clang-format off",
        "#pragma once",
        "",
        "#include \"tap/dsp/iir_tables.h\"",
        "",
        "namespace tap::dsp {",
        "",
    ]
    lines += body
    lines += ["", "    /// The near-linear half-bands: {A, p num, p den, divisor, N, K, sections, stopband dB, "
              "delay deviation, pole radius}."]
    lines += ["    inline constexpr near_linear_halfband_design k_near_linear_halfband_tables[] = {"] + nl_rows + ["    };"]
    lines += ["", "    /// The four-path stages: {A, p num, p den, divisor, paths, stopband dB on the folding bands}."]
    lines += ["    inline constexpr fourpath_design k_fourpath_tables[] = {"] + fp_rows + ["    };"]
    lines += ["", "    /// The Chebyshev II L-inf gains: {A, p num, p den, order, gains (sections 0 .. count - 2)}."]
    lines += ["    inline constexpr cheby2_gains k_cheby2_gain_tables[] = {"] + g_rows + ["    };"]
    lines += ["", "    /// The type II guard band FIRs: {A, p num, p den, Chebyshev order, stop edge, taps, "
              "ripple dB, stopband dB}."]
    lines += ["    inline constexpr guard_band_fir_design k_guard_band_fir_tables[] = {"] + fir_rows + ["    };"]
    lines += ["", "} // namespace tap::dsp", "", "// clang-format on", "// NOLINTEND(readability-identifier-naming)", ""]
    OUT.write_text("\n".join(lines))
    print(f"wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
