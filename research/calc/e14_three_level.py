"""E14: 211・322・544 の三準位の運動論(322 から育つ BH と、211 の臨界より下への降下)。

1. 固定点の構造: 211 の臨界より上は (211, 322)、下は (322, 544)。どちらも安定でない幅と、降下の時間の E12 との比。
2. 三準位の ODE(fp3.evolve3、1 個の量子から)と chi_sr を比べる: 211 から育つ BH(固定点、自由場の側)と 322 から育つ BH。
3. 窓の数値が変わらないこと。
4. f の走査(図のデータ)。
5. 1 世代目の自然スピン Beta(0.69, 2.28) での印の割合と、計算時間。
出力: out/e14.json、../figs/data/e14_scan.dat、../figs/data/e14_ode.dat
実行: cd research && ~/miniforge3/bin/python3 calc/e14_three_level.py
"""
import json
import os
import subprocess
import sys
import time
import warnings

import numpy as np
from scipy.integrate import trapezoid

warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(__file__))
import chisr  # noqa: E402
import fp3  # noqa: E402
import sr  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
FIG = os.path.join(HERE, "..", "figs", "data")
MU = 0.71e-12

ODE_CASES = [  # (M, chi0, f, [T])。211 から育つ BH と 322 から育つ BH
    (19.6, 0.78, 1e15, [1e8, 1e9, 1e10]), (19.6, 0.78, 1e16, [1e7, 1e8, 1e9, 1e10]),
    (19.6, 0.78, 3e16, [1e7, 1e8, 1e9, 1e10]), (19.6, 0.78, 1e17, [1e7, 1e8, 1e9]), (19.6, 0.78, 1e18, [1e7, 1e8, 1e9]),
    (25.0, 0.8, 1e16, [1e7, 1e8, 3e8, 1e9, 3e9, 1e10]), (30.0, 0.8, 1e15, [1e8, 1e9, 1e10]),
    (45.0, 0.7, 1e13, [1e8, 1e9, 1e10]), (45.0, 0.7, 1e14, [1e8, 1e9, 1e10]), (45.0, 0.7, 1e15, [1e7, 1e8, 1e9]),
    (20.0, 0.3, 1e14, [1e9, 1e10]), (20.0, 0.3, 1e15, [1e8, 1e9, 1e10]), (20.0, 0.3, 1e16, [1e8, 1e9, 1e10]),
    (30.0, 0.5, 1e14, [1e8, 1e9, 1e10]), (30.0, 0.5, 1e15, [1e8, 1e9, 1e10]), (15.0, 0.25, 1e16, [1e9, 1e10]),
]
SCAN_ODE_F = [1e13, 1e14, 1e15, 1e16]


def ode_worker(M, c0, f, Ts):
    s = fp3.evolve3(M, c0, MU, f, Ts[-1], Ts)
    print(json.dumps({"status": int(s.status), "chi": [float(v) for v in s.y[3]],
                      "t_stop": [float(t) / sr.YR for t in s.t_events[0]] if len(s.t_events[0]) else None,
                      "log10eps_end": [float(v) / np.log(10) for v in s.y[:3, -1]] if len(s.t) else None}))


def run_ode(M, c0, f, Ts, timeout=600):
    try:
        r = subprocess.run([sys.executable, os.path.abspath(__file__), "--ode", str(M), str(c0), str(f),
                            ",".join(str(t) for t in Ts)], capture_output=True, text=True, timeout=timeout)
        return json.loads(r.stdout.strip().splitlines()[-1])
    except (subprocess.TimeoutExpired, IndexError, json.JSONDecodeError):
        return {"status": -9, "chi": []}


def structure(M):
    """固定点の型が変わる場所と、322 の臨界 + 0.01 から 211 の臨界までの時間の積分(E12 の二準位との比)。"""
    a = sr.alpha_of(M, MU)
    c1, c2 = chisr._zero(a, 2, 1, 1), chisr._zero(a, 3, 2, 2)
    x = np.concatenate([np.linspace(c2 + 0.01, c1 - 1e-2, 300, endpoint=False), c1 - np.geomspace(1e-2, 1e-10, 140)])
    st2, st3, h2, h3 = [], [], [], []
    for chi in x:
        g, k = fp3.coeffs(a, chi)
        p2, p3 = fp3.fp_P2(g, k), fp3.fp_P3(g, k)
        st2.append(fp3.percap(p2, g, k)[0] < 0)
        st3.append(fp3.percap(p3, g, k)[2] < 0)
        h2.append(fp3.H_of(p2, g))
        h3.append(fp3.H_of(p3, g))
    st2, st3, h2, h3 = map(np.array, (st2, st3, h2, h3))
    gap = ~st2 & ~st3
    Hs = np.where(st3, h3, np.where(st2, h2, np.maximum(h2, h3)))
    Hf = np.where(st3, h3, np.where(st2, h2, np.minimum(h2, h3)))
    I_s, I_f, I_2 = trapezoid(1 / np.abs(Hs), x), trapezoid(1 / np.abs(Hf), x), trapezoid(1 / np.abs(h3), x)
    return {"M": M, "alpha": round(a, 4), "chi_c211": round(c1, 5), "chi_c322": round(c2, 5),
            "P3_stable_within_below_c211": float(c1 - x[st3].min()) if st3.any() else 0.0,
            "gap_width_below_c211": float(c1 - x[gap].min()) if gap.any() else 0.0,
            "gap_share_of_time": float(I_s / I_f - 1), "slower_than_E12": float(I_s / I_2)}


def main():
    res = {}
    tab = chisr.build()

    # ---------- 1. 固定点の構造 ----------
    res["structure"] = [structure(M) for M in (19.6, 25.0, 30.0, 45.0)]

    # ---------- 2. ODE と ----------
    rows, worst = [], {"fp": 0.0, "free": 0.0}
    for M, c0, f, Ts in ODE_CASES:
        o = run_ode(M, c0, f, Ts)
        c, r = chisr.chi_sr(tab, M, c0, np.array(Ts), MU, f)
        for k, T in enumerate(Ts):
            ode = o["chi"][k] if k < len(o["chi"]) else None
            rows.append({"M": M, "chi0": c0, "f": f, "T": T, "ode": None if ode is None else round(ode, 5),
                         "chi_sr": round(float(c[k]), 5), "regime": int(r[k]), "ode_status": o["status"],
                         "ode_stop_yr": o.get("t_stop")})
            if ode is not None:
                key = "free" if r[k] == 2 else "fp"
                worst[key] = max(worst[key], abs(float(c[k]) - ode))
    res["vs_ode3"] = rows
    res["vs_ode3_max_abs_diff"] = worst

    # ---------- 3. 窓 ----------
    res["window_check"] = {"GW241011_T5e7_f1e12": float(chisr.chi_sr(tab, 19.6, 0.78, 5e7, MU, 1e12)[0]),
                           "M15_chi0.7_T6e9_f1.4e12": float(chisr.chi_sr(tab, 15.0, 0.7, 6e9, MU, 1.38e12)[0])}

    # ---------- 4. f の走査 ----------
    fs = np.logspace(11, 19, 321)
    cols, hdr = [fs], ["f"]
    for (M, c0, T, lab) in [(19.6, 0.78, 1e9, "GW241011"), (45.0, 0.7, 1e9, "M45"), (20.0, 0.3, 6e9, "M20")]:
        c, r = chisr.chi_sr(tab, M, c0, T, MU, fs)
        cols += [c, r]
        hdr += [f"chi_{lab}", f"reg_{lab}"]
        res.setdefault("scan_max_step", {})[lab] = float(np.abs(np.diff(c)).max())
    np.savetxt(os.path.join(FIG, "e14_scan.dat"), np.column_stack(cols), header=" ".join(hdr), comments="", fmt="%.6e")
    pts = []
    for (M, c0, T, lab) in [(19.6, 0.78, 1e9, "GW241011"), (45.0, 0.7, 1e9, "M45"), (20.0, 0.3, 6e9, "M20")]:
        for f in SCAN_ODE_F:
            o = run_ode(M, c0, f, [T])
            pts.append((f, o["chi"][0] if o["chi"] else np.nan, {"GW241011": 0, "M45": 1, "M20": 2}[lab]))
    np.savetxt(os.path.join(FIG, "e14_ode.dat"), np.array(pts), header="f chi which", comments="", fmt="%.6e")

    # ---------- 5. 印の割合と計算時間 ----------
    rng = np.random.default_rng(12)
    n = 200000
    Mpop = rng.uniform(8, 50, n)
    c0pop = np.minimum(rng.beta(0.69, 2.28, n), 0.999)
    res["natal_beta_regimes_T6Gyr"] = {}
    for f in [1e12, 1e14, 1e16, 1e18]:
        _, r = chisr.chi_sr(tab, Mpop, c0pop, 6e9, MU, f)
        res["natal_beta_regimes_T6Gyr"][f"{f:.0e}"] = {str(k): round(float(np.mean(r == k)), 3) for k in range(5)}
    nev = 1_000_000
    M = rng.uniform(5, 80, nev)
    c0 = rng.uniform(0.0, 0.99, nev)
    T = 10 ** rng.uniform(6, 10, nev)
    mus = rng.uniform(0.5e-12, 1.0e-12, nev)
    fs = 10 ** rng.uniform(11, 18, nev)
    t0 = time.time()
    chisr.chi_sr(tab, M, c0, T, mus, fs)
    res["timing_1e6_s"] = time.time() - t0

    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "e14.json"), "w") as fo:
        json.dump(res, fo, indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--ode":
        ode_worker(float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), [float(x) for x in sys.argv[5].split(",")])
    else:
        main()
