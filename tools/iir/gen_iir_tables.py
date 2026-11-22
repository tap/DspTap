#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright 2026 Timothy Place and the DspTap contributors.
"""Generates include/tap/dsp/detail/iir_tables_data.h: the committed recursive
designs iir_tables.h declares (SampleRateTap's docs/IIR_PLAN.md, I4 and I13).

    python3 tools/iir/gen_iir_tables.py     # from the repository root, in the
                                            # venv tools/iir/requirements.txt pins

Versions the committed header was generated with: Python 3.11.18, numpy 2.4.6,
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
  peak to peak to p and A dB down from the edge), where it converges;
- its Q15 rows (DspTap 0.4.0) at the three 70 dB points and the same orders
  and edges: the shortest even length (from 4) at which a remez design,
  targeted at A to A + 30 dB in 0.01 dB steps (one weight ratio per target),
  keeps that spec on the 16384-point grids once its taps are quantized to
  Q1.14 by quantize_symmetric_row_preserving_sum<std::int16_t> (emulated bit
  for bit, iir_design.quantize_symmetric_q15); at that length, the design
  whose Q1.14 taps attain the deepest stopband. No row at 120 dB: Q1.14
  cannot carry it (tap/dsp/iir_tables.h). The Q1.14 taps the emulation made
  go to tests/reference/guard_band_q15_taps.h, where the C++ battery checks
  DspTap's own quantizer against them.

Run time: a few minutes for the four-path searches, about twenty for the
near-linear ones (N19 at 120 dB, 5/12 dominates), and a few for the Q15 rows.
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
OUT_Q15_REF = ROOT / "tests" / "reference" / "guard_band_q15_taps.h"


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


def int_lines(xs, indent, per_line=16):
    lines = []
    for i in range(0, len(xs), per_line):
        lines.append(indent + ", ".join(str(x) for x in xs[i:i + per_line]) + ("," if i + per_line < len(xs) else ""))
    return lines


def main():
    body = []
    nl_rows, fp_rows, g_rows, fir_rows = [], [], [], []
    q15_body, q15_rows, ref_body, ref_rows = [], [], [], []
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
            if A != 70:
                continue
            q15 = d.typeii_q15_min_length(p, float(edge), A)
            if q15 is None:
                raise RuntimeError(f"no Q15 type II row at {A} dB {pn}/{pd} order {order}")
            nq, target, hq, taps, rq, sbq, passing, tried = q15
            if any(float(fmt(x)) != float(x) for x in hq):
                raise RuntimeError("a design tap does not survive its decimal form")
            name = f"k_guard_fir_q15_{tag(A, pn, pd)}_{order}"
            q15_body.append(f"    inline constexpr double {name}[] = {{")
            q15_body += arr_lines(hq, "        ")
            q15_body.append("    };")
            q15_rows.append(f"        {{{A}, {pn}, {pd}, {order}, {fmt(float(edge))}, {name}, {target:.2f}, "
                            f"{fmt(rq)}, {fmt(sbq)}}},")
            ref = f"k_{tag(A, pn, pd)}_{order}"
            ref_body.append(f"    inline constexpr std::int16_t {ref}[] = {{")
            ref_body += int_lines(taps, "        ")
            ref_body.append("    };")
            ref_rows.append(f"        {{{A}, {pn}, {pd}, {order}, {ref}}},")
            print(f"  {A} dB {pn}/{pd} order {order}: Q15 type II {nq} taps (double {n}), designed for "
                  f"{target:.2f} dB ({passing} of {tried} designs pass at {nq}), Q1.14 ripple {rq:.4f} dB, "
                  f"stopband {sbq:.2f} dB", flush=True)

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
    lines += body + q15_body
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
    lines += ["", "    /// The type II guard band FIRs for Q15: {A, p num, p den, Chebyshev order, stop edge, taps, "
              "design dB, ripple dB, stopband dB}, the last two attained by the Q1.14 taps."]
    lines += ["    inline constexpr guard_band_fir_q15_design k_guard_band_fir_q15_tables[] = {"] + q15_rows + ["    };"]
    lines += ["", "} // namespace tap::dsp", "", "// clang-format on", "// NOLINTEND(readability-identifier-naming)", ""]
    OUT.write_text("\n".join(lines))
    print(f"wrote {OUT.relative_to(ROOT)}")

    ref = [
        "// Generated by tools/iir/gen_iir_tables.py: the Q1.14 taps its emulation of",
        "// tap::dsp::quantize_symmetric_row_preserving_sum<std::int16_t> made of each",
        "// committed Q15 guard band row (tap/dsp/iir_tables.h, k_guard_band_fir_q15_tables).",
        "// Do not edit; regenerate. tests/test_iir_design.cpp quantizes the same rows",
        "// with DspTap's own quantizer and checks it against these, bit for bit.",
        "// SPDX-License-Identifier: MIT",
        "// Copyright 2026 Timothy Place and the DspTap contributors.",
        "// NOLINTBEGIN(readability-identifier-naming)",
        "// Generated data; clang-format would re-pack the tables, so the file opts out.",
        "// clang-format off",
        "#pragma once",
        "",
        "#include <cstdint>",
        "#include <span>",
        "",
        "namespace tap::dsp::test::guard_band_q15_ref {",
        "",
        "    struct row {",
        "        int                           stop_db;",
        "        int                           pass_num;",
        "        int                           pass_den;",
        "        int                           order;",
        "        std::span<const std::int16_t> taps;",
        "    };",
        "",
    ]
    ref += ref_body
    ref += ["", "    inline constexpr row k_rows[] = {"] + ref_rows + ["    };"]
    ref += ["", "} // namespace tap::dsp::test::guard_band_q15_ref", "", "// clang-format on",
            "// NOLINTEND(readability-identifier-naming)", ""]
    OUT_Q15_REF.write_text("\n".join(ref))
    print(f"wrote {OUT_Q15_REF.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
