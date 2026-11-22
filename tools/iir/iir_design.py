# SPDX-License-Identifier: MIT
# Copyright 2026 Timothy Place and the DspTap contributors.
"""The recursive designs behind include/tap/dsp/iir_design.h and iir_tables.h.

Imported by gen_iir_tables.py (the committed tables) and gen_cx_reference.py
(the mpmath leg of the constexpr designers). Adapted from SampleRateTap's
IIR-0 spike (rational/tools/spikes/iir/python/iir_design.py, MIT, the same
copyright holder), whose README carries the derivations and the numbers this
module reproduces; NOTICE.md keeps the clean-room record and the sources.
Versions: tools/iir/requirements.txt (Python 3.11, numpy 2.4.6, scipy 1.17.1,
mpmath 1.4.1).

Conventions (iir_design.h's file header):
- A spec point is (A dB, p = pass_num / pass_den of the chain's lowest rate)
  and a design divisor d; a by-2 stage at divisor d has its passband edge at
  f_pass = p / (2 d) of its higher rate and its stopband from 1/2 - f_pass
  to 1/2. Every design meets -(A + MARGIN_DB) on GRID points of it.
- Frequencies f are cycles per sample of the rate named; theta is the path
  (low) rate's angle.
- H(z) = 1/2 [A0(z^2) + z^-1 A1(z^2)] (Valenzuela and Constantinides 1983;
  Vaidyanathan 1993, the allpass decomposition of a power-symmetric filter).
"""
from __future__ import annotations

import itertools
import math

import mpmath
import numpy as np
from scipy import optimize, signal

# The four spec points (the family's profiles' (A, p); the names are
# rational's and appear here only in comments), the divisors a power-of-two
# cascade uses, the margin and the grid (rational's design.h:
# k_design_margin_db, k_design_grid_points).
SPEC_POINTS = ((70, 1, 3), (70, 3, 8), (70, 19, 48), (120, 5, 12))
DIVISORS = (1, 2, 4, 8)
CHEBY2_ORDERS = (3, 5, 7, 9, 11)
RIPPLE_DB = 0.01
MARGIN_DB = 1.0
GRID = 16384


def db(x):
    return 20.0 * np.log10(np.maximum(np.abs(x), 1e-300))


def spec_p(pass_num, pass_den):
    """p as the C++ computes it: double(pass_num) / double(pass_den)."""
    return float(pass_num) / float(pass_den)


# ---------------------------------------------------------------------------
# The nonlinear (elliptic) half-band: Gazsi 1985's explicit formulas in the
# bireciprocal case, at high precision for the mpmath leg. Evaluated in the
# published form (e = (1 - k'^(1/2)) / (2 (1 + k'^(1/2))), a = (1 - x) /
# (1 + x)): at 60 digits the cancellations the C++ removes by algebra cost
# nothing, so the leg also checks the rewrites.
# ---------------------------------------------------------------------------
def halfband_mp(order, f_pass, dps=60):
    """The (order - 1) / 2 allpass coefficients, ascending, as mpf."""
    with mpmath.workdps(dps):
        f = mpmath.mpf(f_pass)
        k = mpmath.tan(mpmath.pi * f) ** 2
        kk = mpmath.sqrt(mpmath.sqrt(1 - k * k))
        e = (1 - kk) / (2 * (1 + kk))
        q = e + 2 * e ** 5 + 15 * e ** 9 + 150 * e ** 13
        n = order
        out = []
        for c in range(1, (order - 1) // 2 + 1):
            num = mpmath.fsum((-1) ** i * q ** (i * (i + 1)) * mpmath.sin((2 * i + 1) * c * mpmath.pi / n)
                              for i in range(20))
            den = 1 + 2 * mpmath.fsum((-1) ** i * q ** (i * i) * mpmath.cos(2 * i * c * mpmath.pi / n)
                                      for i in range(1, 20))
            w = 2 * q ** mpmath.mpf(0.25) * num / den
            w2 = w * w
            x = mpmath.sqrt((1 - w2 * k) * (1 - w2 / k)) / (1 + w2)
            out.append((1 - x) / (1 + x))
        return sorted(out)


def halfband_order(stop_db, f_pass):
    """The smallest odd N >= 3 with 10 log10(1 + q^(-N/2) / 4) >= stop_db."""
    with mpmath.workdps(40):
        f = mpmath.mpf(f_pass)
        k = mpmath.tan(mpmath.pi * f) ** 2
        kk = mpmath.sqrt(mpmath.sqrt(1 - k * k))
        e = (1 - kk) / (2 * (1 + kk))
        q = e + 2 * e ** 5 + 15 * e ** 9 + 150 * e ** 13
        n = 3
        while 10 * mpmath.log10(1 + q ** (-mpmath.mpf(n) / 2) / 4) < stop_db:
            n += 2
        return n


def nonlinear_design(A, p, d):
    """(order, coefficients as float) at a spec point."""
    f_pass = p / (2.0 * d)
    order = halfband_order(A + MARGIN_DB, f_pass)
    return order, [float(a) for a in halfband_mp(order, f_pass)]


def first_order_cascade(coefs, w):
    out = np.ones_like(w)
    wi = 1.0 / w
    for a in coefs:
        out = out * (a + wi) / (1 + a * wi)
    return out


def nonlinear_response(coefs, f):
    """1/2 [A0(z^2) + z^-1 A1(z^2)], coefficients alternating between the paths."""
    z = np.exp(2j * np.pi * np.asarray(f, dtype=float))
    w = z * z
    a = sorted(coefs)
    return 0.5 * (first_order_cascade(a[0::2], w) + first_order_cascade(a[1::2], w) / z)


# ---------------------------------------------------------------------------
# The near-linear half-band (IIR-0's method): A0 = w^-K, A1 a real allpass of
# order N in w, the minimax of the stopband, which by the structure's identity
# is the passband phase error |sin(eps / 2)| at the mirror frequency.
# ---------------------------------------------------------------------------
def allpass_phase(den, theta):
    den = np.asarray(den, dtype=float)
    N = len(den) - 1
    E = np.exp(-1j * np.outer(theta, np.arange(N + 1)))
    D = E @ den
    return -N * theta - 2 * np.unwrap(np.angle(D))


def allpass_response(den, theta):
    den = np.asarray(den, dtype=float)
    N = len(den) - 1
    E = np.exp(-1j * np.outer(theta, np.arange(N + 1)))
    D = E @ den
    return np.exp(-1j * N * theta) * np.conj(D) / D


def allpass_group_delay(den, theta):
    den = np.asarray(den, dtype=float)
    N = len(den) - 1
    k = np.arange(N + 1)
    E = np.exp(-1j * np.outer(theta, k))
    D = E @ den
    dD = E @ (-1j * k * den)
    return N + 2 * np.imag(dD / D)


def nearlinear_phase_error(den, K, theta):
    """eps(theta) = arg A1 + (K - 1/2) theta."""
    return allpass_phase(den, theta) + (K - 0.5) * theta


def nearlinear_lp(N, K, theta_p, n_grid=None, iters=12):
    """The iterated linear program in the denominator coefficients d_k:
    sum_k d_k sin(k theta + beta) = |D| sin(eps / 2), beta = (K - 1/2 - N)
    theta / 2, is linear in d_k, so with the weight 1 / |D| frozen from the
    previous iterate the minimax of |sin(eps / 2)| is a linear program
    (scipy's HiGHS), iterated until it stops improving."""
    if n_grid is None:
        n_grid = max(400, 40 * N)
    theta = np.linspace(0.0, theta_p, n_grid)
    beta = (K - 0.5 - N) * theta / 2.0
    S = np.sin(np.outer(theta, np.arange(N + 1)) + beta[:, None])
    wgt = np.ones(n_grid)
    best = None
    for _ in range(iters):
        A = S[:, 1:] * wgt[:, None]
        b0 = S[:, 0] * wgt
        A_ub = np.vstack([np.hstack([A, -np.ones((n_grid, 1))]), np.hstack([-A, -np.ones((n_grid, 1))])])
        b_ub = np.concatenate([-b0, b0])
        c = np.zeros(N + 1)
        c[-1] = 1.0
        res = optimize.linprog(c, A_ub=A_ub, b_ub=b_ub, bounds=[(None, None)] * N + [(0, None)], method="highs")
        if res.status != 0:
            break
        den = np.concatenate([[1.0], res.x[:N]])
        if np.max(np.abs(np.roots(den))) >= 1.0:
            break
        err = np.abs(np.sin(nearlinear_phase_error(den, K, theta) / 2)).max()
        if best is None or err < best[1] * (1 - 1e-9):
            best = (den, err)
        else:
            break
        E = np.exp(-1j * np.outer(theta, np.arange(N + 1)))
        wgt = 1.0 / np.abs(E @ den)
    return best


def nearlinear_polish(den, K, theta_p, n=1500):
    """SLSQP on the true objective's epigraph, analytic Jacobian; kept only
    when stable and better on the dense grid."""
    N = len(den) - 1
    th = np.linspace(0, theta_p, n)
    k = np.arange(N + 1)
    E = np.exp(-1j * np.outer(th, k))

    def eps(x):
        return nearlinear_phase_error(np.r_[1, x], K, th)

    def cons(v):
        e = eps(v[:-1])
        return np.r_[v[-1] - e, v[-1] + e]

    def cjac(v):
        D = E @ np.r_[1, v[:-1]]
        J = -2 * np.imag(E[:, 1:] / D[:, None])
        one = np.ones((n, 1))
        return np.vstack([np.hstack([-J, one]), np.hstack([J, one])])

    v0 = np.r_[den[1:], np.abs(eps(den[1:])).max()]
    try:
        r = optimize.minimize(lambda v: v[-1], v0, jac=lambda v: np.r_[np.zeros(N), 1.0], method="SLSQP",
                              constraints=[{"type": "ineq", "fun": cons, "jac": cjac}],
                              options={"maxiter": 500, "ftol": 1e-20})
    except Exception:  # noqa: BLE001 - the polish is optional
        return den
    d2 = np.r_[1.0, r.x[:-1]]
    if not np.all(np.isfinite(d2)) or np.max(np.abs(np.roots(d2))) >= 1.0:
        return den
    dense = np.linspace(0.0, theta_p, GRID)
    e1 = np.abs(np.sin(nearlinear_phase_error(den, K, dense) / 2)).max()
    e2 = np.abs(np.sin(nearlinear_phase_error(d2, K, dense) / 2)).max()
    return d2 if e2 < e1 else den


def nearlinear_sos(den):
    """A1 as sections: second-order (c1, c2) by increasing pole radius (real
    poles paired, the smallest radii together), then one first-order c1."""
    poles = np.roots(den)
    cplx = sorted([p for p in poles if p.imag > 1e-12], key=abs)
    real = sorted([p.real for p in poles if abs(p.imag) <= 1e-12], key=abs)
    secs = []
    for p in cplx:
        secs.append((-2 * p.real, abs(p) ** 2, abs(p)))
    while len(real) >= 2:
        r1, r2 = real.pop(0), real.pop(0)
        secs.append((-(r1 + r2), r1 * r2, max(abs(r1), abs(r2))))
    secs.sort(key=lambda s: s[2])
    out = [(float(s[0]), float(s[1])) for s in secs]
    if real:
        out.append((float(-real[0]),))
    return out


def sos_response(sos, theta):
    """A1 from its sections at the path angle theta."""
    wi = np.exp(-1j * np.asarray(theta, dtype=float))
    out = np.ones_like(wi)
    for s in sos:
        if len(s) == 2:
            c1, c2 = s
            out = out * (c2 + c1 * wi + wi * wi) / (1 + c1 * wi + c2 * wi * wi)
        else:
            (c1,) = s
            out = out * (c1 + wi) / (1 + c1 * wi)
    return out


def nearlinear_response_sos(sos, K, f):
    f = np.asarray(f, dtype=float)
    z = np.exp(2j * np.pi * f)
    theta = 4 * np.pi * f
    return 0.5 * (np.exp(-1j * K * theta) + sos_response(sos, theta) / z)


def nearlinear_group_delay(den, K, f):
    """High-rate samples (-Im(H'/H) from the paths' delays)."""
    f = np.asarray(f, dtype=float)
    theta = 4 * np.pi * f
    z = np.exp(2j * np.pi * f)
    A0 = np.exp(-1j * K * theta)
    A1 = allpass_response(den, theta)
    t0 = np.full_like(f, 2.0 * K)
    t1 = 2 * allpass_group_delay(den, theta)
    H = 0.5 * (A0 + A1 / z)
    dH = 0.5 * (-1j * t0 * A0 - 1j * (1 + t1) * A1 / z)
    return -np.imag(dH / H)


def nearlinear_design(A, p, d, max_order=80, K_offsets=(-1, 0, 1, 2)):
    """The lowest N (and its best K) meeting -(A + MARGIN_DB): (N, K, den)."""
    fp_hi = p / (2.0 * d)
    theta_p = 4.0 * np.pi * fp_hi
    target = 10 ** (-(A + MARGIN_DB) / 20.0)
    dense = np.linspace(0.0, theta_p, GRID)
    for N in range(1, max_order + 1):
        cands = []
        for off in K_offsets:
            K = N + off
            if K < 1:
                continue
            r = nearlinear_lp(N, K, theta_p)
            if r is None:
                continue
            den = nearlinear_polish(r[0], K, theta_p)
            err = np.abs(np.sin(nearlinear_phase_error(den, K, dense) / 2)).max()
            cands.append((err, K, den))
        if not cands:
            continue
        err, K, den = min(cands, key=lambda c: c[0])
        if err <= target:
            return N, K, den
    raise RuntimeError("no near-linear design")


# ---------------------------------------------------------------------------
# The four-path stage (harris 2004; Renfors and Saramaki 1987): H(z) = 1/4
# sum_k z^-k A_k(z^4). IIR-0's design: SLSQP minimax over the bands that fold
# onto the passband, from random starts.
# ---------------------------------------------------------------------------
def fourpath_response(paths, f):
    z = np.exp(2j * np.pi * np.asarray(f, dtype=float))
    zi4 = z ** -4
    tot = np.zeros_like(z)
    for k, cs in enumerate(paths):
        A = np.ones_like(z)
        for a in cs:
            A = A * (a + zi4) / (1 + a * zi4)
        tot = tot + A * z ** (-k)
    return tot / 4


def fourpath_bands(p, D):
    """[1/4 - f, 1/4 + f] and [1/2 - f, 1/2] of the input rate 4 D r_min, f = p / (4 D)."""
    fp = p / (4.0 * D)
    return [(0.25 - fp, 0.25 + fp), (0.5 - fp, 0.5)]


def fourpath_design(A, p, D, layouts, starts=24, seed=1):
    """The best minimax design over the layouts: (stopband dB, layout, paths)."""
    bands = fourpath_bands(p, D)
    fst = np.concatenate([np.linspace(lo, hi, 240) for lo, hi in bands])
    fst_dense = np.concatenate([np.linspace(lo, hi, GRID // 2) for lo, hi in bands])
    rng = np.random.default_rng(seed)
    best = None
    for layout in layouts:
        n = sum(layout)

        def unpack(x, layout=layout):
            out, i = [], 0
            for m in layout:
                out.append(list(x[i:i + m]))
                i += m
            return out

        def cons(v, unpack=unpack):
            return v[-1] - np.abs(fourpath_response(unpack(v[:-1]), fst)) ** 2

        for _ in range(starts):
            x0 = rng.uniform(0.02, 0.95, n)
            v0 = np.concatenate([x0, [np.abs(fourpath_response(unpack(x0), fst)).max() ** 2]])
            try:
                r = optimize.minimize(lambda v: v[-1], v0, method="SLSQP",
                                      constraints=[{"type": "ineq", "fun": cons}],
                                      bounds=[(-0.999, 0.999)] * n + [(0, 1)],
                                      options={"maxiter": 400, "ftol": 1e-18})
            except Exception:  # noqa: BLE001 - a failed start is skipped
                continue
            paths = unpack(r.x[:-1])
            sb = db(fourpath_response(paths, fst_dense)).max()
            if best is None or sb < best[0]:
                best = (float(sb), list(layout), [sorted(float(a) for a in pth) for pth in paths])
    return best


def distinct_permutations(layout):
    return sorted(set(itertools.permutations(layout)))


# ---------------------------------------------------------------------------
# The odd-order Chebyshev II (Antoniou 2006; the inverse-Chebyshev case of
# Gazsi 1985), at high precision: sections at unity DC gain, ordered by
# increasing pole radius, exactly as design_cheby2_lowpass builds them, and
# the L-inf gains (Jackson 1970).
# ---------------------------------------------------------------------------
def cheby2_mp(order, p, A, rp=RIPPLE_DB, dps=60):
    """[(b0, b1, b2, a1, a2, radius)] as mpf, and the stop edge."""
    with mpmath.workdps(dps):
        n = order
        m = (n - 1) // 2
        a10 = mpmath.mpf(10) ** (mpmath.mpf(A) / 10) - 1
        r10 = mpmath.mpf(10) ** (mpmath.mpf(rp) / 10) - 1
        op = mpmath.tan(mpmath.pi * mpmath.mpf(p))
        D = mpmath.sqrt(a10 / r10)
        os_ = op * mpmath.cosh(mpmath.acosh(D) / n)
        mu = mpmath.asinh(mpmath.sqrt(a10)) / n
        sh, ch = mpmath.sinh(mu), mpmath.cosh(mu)
        secs = []
        for k in range(1, m + 1):
            th = mpmath.pi * (2 * k - 1) / (2 * n)
            cheb = mpmath.mpc(-sh * mpmath.sin(th), ch * mpmath.cos(th))
            pz = os_ / cheb
            zp = (1 + pz) / (1 - pz)
            oz = os_ / mpmath.cos(th)
            cz = (1 - oz * oz) / (1 + oz * oz)
            g = abs(1 - zp) ** 2 / (2 * (1 - cz))
            secs.append((g, -2 * cz * g, g, -2 * zp.real, abs(zp) ** 2, abs(zp)))
        pr = -os_ / sh
        zr = (1 + pr) / (1 - pr)
        g = (1 - zr) / 2
        secs.append((g, g, mpmath.mpf(0), -zr, mpmath.mpf(0), abs(zr)))
        secs.sort(key=lambda s: s[5])
        edge = mpmath.atan(os_) / mpmath.pi
        return secs, edge


def sections_mag_mp(secs, f):
    zi = mpmath.expj(-2 * mpmath.pi * f)
    zi2 = zi * zi
    r = mpmath.mpc(1)
    for b0, b1, b2, a1, a2, *_ in secs:
        r *= (b0 + b1 * zi + b2 * zi2) / (1 + a1 * zi + a2 * zi2)
    return abs(r)


def cheby2_linf_gains(secs, grid=2048, dps=40):
    """gains[i] = 1 / peak of the partial cascade 0 .. i (earlier sections
    scaled), for every section but the last; the peak on a grid over [0, 1/2]
    refined by golden-section search around the largest."""
    with mpmath.workdps(dps):
        scaled = [list(s) for s in secs]
        gains = []
        inv_phi = (mpmath.sqrt(5) - 1) / 2
        for i in range(len(secs) - 1):
            part = scaled[:i + 1]
            vals = [sections_mag_mp(part, mpmath.mpf(g) / (2 * grid)) for g in range(grid + 1)]
            best = max(range(grid + 1), key=lambda g: vals[g])
            lo = mpmath.mpf(max(best - 1, 0)) / (2 * grid)
            hi = mpmath.mpf(min(best + 1, grid)) / (2 * grid)
            for _ in range(90):
                m1 = hi - inv_phi * (hi - lo)
                m2 = lo + inv_phi * (hi - lo)
                if sections_mag_mp(part, m1) < sections_mag_mp(part, m2):
                    lo = m1
                else:
                    hi = m2
            peak = max(vals[best], sections_mag_mp(part, (lo + hi) / 2))
            sc = 1 / peak
            gains.append(sc)
            for j in range(3):
                scaled[i][j] *= sc
        return gains


# ---------------------------------------------------------------------------
# The linear guard band form: the shortest even-length (type II) equiripple
# FIR (Parks-McClellan, scipy.signal.remez) flat within 0.01 dB peak to peak
# to p and A dB down from f_stop.
# ---------------------------------------------------------------------------
def typeii_min_length(p_rate, f_stop, A, rp_pp=RIPPLE_DB, max_len=400):
    dp = (10 ** (rp_pp / 20) - 1) / (10 ** (rp_pp / 20) + 1)
    ds = 10 ** (-A / 20)
    for N in range(4, max_len, 2):
        try:
            h = signal.remez(N, [0, p_rate, f_stop, 0.5], [1, 0], weight=[1 / dp, 1 / ds], fs=1.0, maxiter=400)
        except Exception:  # noqa: BLE001 - remez fails to converge at some lengths
            continue
        w, H = signal.freqz(h, worN=2 ** 15, fs=1.0)
        pb = db(H[w <= p_rate])
        sb = db(H[w >= f_stop]).max()
        if pb.max() - pb.min() <= rp_pp + 1e-9 and sb <= -A:
            return N, h, float(pb.max() - pb.min()), float(sb)
    return None


# ---------------------------------------------------------------------------
# The linear guard band form in Q15 (DspTap 0.4.0): the design whose Q1.14
# taps, quantized by include/tap/dsp/quantize.h's
# quantize_symmetric_row_preserving_sum<std::int16_t>, keep the spec. The
# quantizer is emulated here bit for bit: every step below is one IEEE double
# operation the C++ performs in the same order (the product by 2^14 is exact,
# so a contracting compiler cannot move it either), and
# tests/test_iir_design.cpp checks DspTap's own quantizer against the taps
# this emulation produced (tests/reference/guard_band_q15_taps.h).
# ---------------------------------------------------------------------------
Q15_SCALE = 16384.0  # sample_traits<std::int16_t>::k_coeff_scale: Q1.14 units per 1.0
Q15_MIN, Q15_MAX = -32768, 32767
# The design targets the search tries at every length: A to A + 30 dB in
# 0.01 dB steps. In remez only the ratio of the two band weights matters, so
# a higher target and more stopband weight are one lever, and this one
# parameter spans both.
Q15_TARGET_STEPS = 3000
Q15_TARGET_STEP_DB = 0.01
# A pass must clear the spec by this much on the search's evaluation, so a
# re-measurement in another summation order (the C++ test's) cannot flip it.
Q15_MEASURE_GUARD_DB = 1e-6


def round_sat_q15(v):
    """sample_traits' detail::round_sat<std::int16_t>: half away from zero
    (v -/+ 0.5 in double, then truncation), saturating."""
    r = v - 0.5 if v < 0.0 else v + 0.5
    if r <= float(Q15_MIN):
        return Q15_MIN
    if r >= float(Q15_MAX):
        return Q15_MAX
    return int(r)  # truncation toward zero, as static_cast


def llround(x):
    """std::llround: to nearest, half away from zero (a - floor(a) is exact)."""
    a = abs(x)
    f = math.floor(a)
    if a - f >= 0.5:
        f += 1
    return -int(f) if x < 0.0 else int(f)


def quantize_row_preserving_sum_q15(src):
    """quantize_row_preserving_sum<std::int16_t>, step for step."""
    n = len(src)
    dst = [0] * n
    remainder = [0.0] * n
    exact_sum = 0.0
    quant_sum = 0
    for t in range(n):
        scaled = float(src[t]) * Q15_SCALE
        q = round_sat_q15(float(src[t]) * Q15_SCALE)
        dst[t] = q
        remainder[t] = scaled - float(q)
        exact_sum += scaled
        quant_sum += q
    residual = llround(exact_sum) - quant_sum
    while residual != 0:
        sgn = 1.0 if residual > 0 else -1.0
        rail = Q15_MAX if residual > 0 else Q15_MIN
        best = n
        for u in range(n):
            if dst[u] != rail and (best == n or sgn * remainder[u] > sgn * remainder[best]):
                best = u
        if best == n:
            break
        step = 1 if residual > 0 else -1
        dst[best] += step
        remainder[best] -= sgn
        residual -= step
    return dst


def quantize_symmetric_q15(h):
    """quantize_symmetric_row_preserving_sum<std::int16_t>: the first half,
    then its mirror."""
    half = quantize_row_preserving_sum_q15(list(h[: len(h) // 2]))
    return half + half[::-1]


def typeii_cos_matrix(f, n):
    """|H(f)| of an even-length symmetric row is |2 sum_{i < n/2} h[i]
    cos(2 pi f ((n - 1) / 2 - i))|: the matrix over the frequencies f."""
    i = np.arange(n // 2)
    return 2.0 * np.cos(2.0 * np.pi * np.outer(f, (n - 1) / 2.0 - i))


def band_grids(p_rate, f_stop, grid=GRID):
    """measure_fir_bands' grids (iir_design.h): grid + 1 points over
    [f_stop, 1/2] and over [0, p]."""
    t = np.arange(grid + 1) / grid
    return f_stop + (0.5 - f_stop) * t, p_rate * t


def fir_bands(h, p_rate, f_stop, grid=GRID):
    """(worst stopband dB over [f_stop, 1/2], passband ripple dB peak to peak
    over [0, p]) on measure_fir_bands' grids, by the complex sum."""
    fs, fp = band_grids(p_rate, f_stop, grid)
    h = np.asarray(h, dtype=float)
    i = np.arange(len(h))
    hs = np.abs(np.exp(-2j * np.pi * np.outer(fs, i)) @ h)
    hp = np.abs(np.exp(-2j * np.pi * np.outer(fp, i)) @ h)
    return float(db(hs.max())), float(db(hp.max()) - db(hp.min()))


def typeii_q15_min_length(p_rate, f_stop, A, rp_pp=RIPPLE_DB, max_len=400):
    """The shortest even length N (from 4) at which a remez design of the
    target grid, its taps quantized to Q1.14 by quantize_symmetric_q15, is
    flat within rp_pp peak to peak to p and A dB down from f_stop on the
    16384-point grids, each by Q15_MEASURE_GUARD_DB; at that length, the
    design whose Q1.14 taps attain the deepest stopband (the lowest target on
    a tie). Returns (N, target dB, design taps, Q1.14 taps, ripple dB and
    stopband dB of the Q1.14 taps, designs passing at N, designs tried at N),
    or None."""
    dp = (10 ** (rp_pp / 20) - 1) / (10 ** (rp_pp / 20) + 1)
    targets = [(100 * A + k) / 100 for k in range(0, Q15_TARGET_STEPS + 1)]
    assert Q15_TARGET_STEP_DB == 0.01
    fs, fp = band_grids(p_rate, f_stop)
    for N in range(4, max_len, 2):
        cs, cp = typeii_cos_matrix(fs, N), typeii_cos_matrix(fp, N)
        tried = []
        for target in targets:
            ds = 10 ** (-target / 20)
            try:
                h = signal.remez(N, [0, p_rate, f_stop, 0.5], [1, 0], weight=[1 / dp, 1 / ds], fs=1.0, maxiter=400)
            except Exception:  # noqa: BLE001 - remez fails to converge at some lengths
                continue
            tried.append((target, h, quantize_symmetric_q15(h)))
        if not tried:
            continue
        halves = np.array([q[: N // 2] for _, _, q in tried], dtype=float).T / Q15_SCALE
        sb = db(np.abs(cs @ halves).max(axis=0))
        mag = np.abs(cp @ halves)
        rp = db(mag.max(axis=0)) - db(mag.min(axis=0))
        g = Q15_MEASURE_GUARD_DB
        passing = [k for k in range(len(tried)) if sb[k] <= -A - g and rp[k] <= rp_pp - g]
        if passing:
            best = min(passing, key=lambda k: (sb[k], tried[k][0]))
            target, h, q = tried[best]
            sbq, rpq = fir_bands([x / Q15_SCALE for x in q], p_rate, f_stop)
            if sbq > -A - g / 2 or rpq > rp_pp - g / 2:
                raise RuntimeError(f"Q15 row N{N}: the complex sum disagrees with the search ({sbq}, {rpq})")
            return N, target, h, q, rpq, sbq, len(passing), len(tried)
    return None


def ulps(x, ref):
    """|x - ref| in units of the last place of x (a double)."""
    return float(abs(mpmath.mpf(x) - ref) / math.ulp(x))
