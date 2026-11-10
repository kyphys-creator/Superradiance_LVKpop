"""E12: 二度目の査読で見つかった問題の検証と修正。

1. 211 の臨界より下への降下: 新しい chi_sr(322 の超放射を含む二準位の固定点)を E5 の ODE と比べる。
   E11 の関数は 211 の臨界で止まっていたので、その値(新しい値と 211 の臨界の大きいほう)も並べる。
2. f の走査(図のデータ)。
3. 窓の年齢依存: GW241011 の主星の年齢 T2 = 5e7, 1e8, 3e8 yr での f の上限と、回転を失う 1 世代目の最小の質量。
4. delta 依存: 15 Msun、6 Gyr で、211 の臨界 + delta まで落とすのに要る f。
5. 1 世代目の自然スピン Beta(0.69, 2.28)(2607.01317 の事後分布の中央値)で、211 が超放射する・322 から育つ・
   超放射しない割合と、chi_sr の印の割合。
6. gamma_BH の相対論的な補正(E7 の C_BH)で、f の上限と 1 世代目の最小の質量がどう動くか。
出力: out/e12.json、../figs/data/e12_scan.dat、../figs/data/e12_ode.dat
実行: cd research && ~/miniforge3/bin/python3 calc/e12_review.py
"""
import json
import os
import subprocess
import sys
import warnings

import numpy as np
from scipy.optimize import brentq
from scipy.stats import beta

warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(__file__))
import chisr  # noqa: E402
import si  # noqa: E402
import sr  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
FIG = os.path.join(HERE, "..", "figs", "data")
MU = 0.71e-12

ODE_CASES = [  # (M, chi0, f, [T])
    (19.6, 0.78, 1e13, [1e9]), (19.6, 0.78, 1e14, [1e9]), (19.6, 0.78, 1e15, [1e8, 1e9, 1e10]),
    (19.6, 0.78, 1e16, [1e7, 1e8, 3e8, 1e9]), (19.6, 0.78, 3e16, [1e7, 3e7, 1e8, 3e8]),
    (25.0, 0.8, 1e14, [1e9, 1e10]), (25.0, 0.8, 1e15, [1e8, 3e8, 1e9, 3e9]),
    (30.0, 0.8, 1e13, [1e9, 1e10]), (30.0, 0.8, 1e14, [1e8, 1e9]), (15.0, 0.7, 1e13, [1e9, 1e10]),
    (19.6, 0.78, 1e17, [1e7, 1e8]),
]


def ode_worker(M, c0, f, Ts):
    """E5 の ODE を 1 本解いて、各時刻のスピンを JSON で出す(別のプロセスで呼ぶ: 解法が止まっても本体は待たない)。"""
    import e05_kinetics as e5
    s = e5.evolve(M, c0, MU, f, Ts[-1], Ts)
    ys = [float(v) for v in s.y[2]] if s.status == 0 else []
    print(json.dumps({"status": int(s.status), "chi": ys}))


def run_ode(M, c0, f, Ts, timeout=300):
    try:
        r = subprocess.run([sys.executable, os.path.abspath(__file__), "--ode", str(M), str(c0), str(f),
                            ",".join(str(t) for t in Ts)], capture_output=True, text=True, timeout=timeout)
        return json.loads(r.stdout.strip().splitlines()[-1])
    except (subprocess.TimeoutExpired, IndexError, json.JSONDecodeError):
        return {"status": -9, "chi": []}


# ---------- 窓(E8 の f_limit、E10 の m_min を初期スピンについて一般にしたもの) ----------
def f_limit(M, chi_i, chi_f, T, mu=MU):
    from e08_window import f_limit as fl
    return fl(M, chi_i, chi_f, T, mu)


def f_1g(M1, c0, dl, T1=6e9):
    thr = si.chi_threshold_211(M1, MU)
    return f_limit(M1, c0, thr + dl, T1) if c0 > thr + dl else np.inf


def m_min(f_up, c0, dl, lo=5.0):
    hi = brentq(lambda M: si.chi_threshold_211(M, MU) + dl - c0, 1.0, 200.0, xtol=1e-3) - 0.01
    g = lambda M: np.log(f_1g(M, c0, dl)) - np.log(f_up)   # noqa: E731
    if g(hi) > 0:
        return np.nan
    if g(lo) < 0:
        return lo
    return brentq(g, lo, hi, xtol=0.02)


def main():
    res = {}
    tab = chisr.build()

    # ---------- 1. ODE と ----------
    rows, ode_pts = [], []
    for M, c0, f, Ts in ODE_CASES:
        o = run_ode(M, c0, f, Ts)
        c_new, reg = chisr.chi_sr(tab, M, c0, np.array(Ts), MU, f)
        c1 = float(chisr.chic_vec(sr.alpha_of(M, MU), 1, 2, tab))
        c2 = float(chisr.chic_vec(sr.alpha_of(M, MU), 2, 3, tab))
        for k, T in enumerate(Ts):
            ode = o["chi"][k] if o["chi"] else None
            rows.append({"M": M, "chi0": c0, "f": f, "T": T, "chi_c211": round(c1, 5), "chi_c322": round(c2, 5),
                         "ode": None if ode is None else round(ode, 5), "new": round(float(c_new[k]), 5),
                         "E11": round(max(float(c_new[k]), c1), 5), "regime": int(reg[k]), "ode_status": o["status"]})
            if ode is not None and M == 19.6 and T == 1e9:
                ode_pts.append((f, ode))
    res["vs_ode"] = rows
    ok = [r for r in rows if r["ode"] is not None]
    res["vs_ode_max_abs_diff_new"] = max(abs(r["new"] - r["ode"]) for r in ok)
    res["vs_ode_max_abs_diff_E11"] = max(abs(r["E11"] - r["ode"]) for r in ok)
    np.savetxt(os.path.join(FIG, "e12_ode.dat"), np.array(sorted(ode_pts)), header="f chi_T9", comments="", fmt="%.6e")

    # ---------- 2. f の走査 ----------
    fs = np.logspace(10, 19, 361)
    cols, hdr = [fs], ["f"]
    for (M, c0, lab) in [(19.6, 0.78, "GW241011"), (15.0, 0.7, "M15"), (30.0, 0.8, "M30")]:
        for T in [1e7, 1e9]:
            c, r = chisr.chi_sr(tab, M, c0, T, MU, fs)
            cols += [c, r]
            hdr += [f"chi_{lab}_T{int(np.log10(T))}", f"reg_{lab}_T{int(np.log10(T))}"]
            res.setdefault("scan_max_step", {})[f"{lab}_T{T:.0e}"] = float(np.abs(np.diff(c)).max())
    np.savetxt(os.path.join(FIG, "e12_scan.dat"), np.column_stack(cols), header=" ".join(hdr), comments="", fmt="%.6e")

    # ---------- 3. 窓の年齢依存 ----------
    win = []
    for T2 in [5e7, 1e8, 3e8]:
        for lab, c0g in [("opt", 0.78), ("real", 0.72)]:
            fup = f_limit(19.6, c0g, 0.69, T2)
            row = {"T2": T2, "set": lab, "chi0_GW241011": c0g, "f_max": fup}
            for c1g in (0.7, 0.5):
                for dl in (0.02, 0.05):
                    row[f"Mmin_chi0{c1g}_d{dl}"] = m_min(fup, c1g, dl)
            win.append(row)
    res["window_vs_age"] = win

    # ---------- 4. delta 依存 ----------
    res["f_low_vs_delta_15Msun_6Gyr"] = {f"chi0={c0}": {f"delta={dl}": f_1g(15.0, c0, dl) for dl in (0.18, 0.05, 0.02, 0.01, 0.001)}
                                         for c0 in (0.5, 0.7, 0.9)}

    # ---------- 5. 1 世代目の自然スピンでの割合 ----------
    b = beta(0.69, 2.28)
    zone = []
    for M in [8, 10, 15, 20, 25, 30, 40, 50]:
        a = sr.alpha_of(M, MU)
        c1 = float(chisr.chic_vec(a, 1, 2, tab))
        c2 = float(chisr.chic_vec(a, 2, 3, tab))
        zone.append({"M": M, "alpha": round(a, 4), "chi_c211": round(c1, 4), "chi_c322": round(c2, 4),
                     "P_211": round(1 - b.cdf(c1), 3), "P_322first": round(b.cdf(c1) - b.cdf(c2), 3),
                     "P_none": round(b.cdf(c2), 3)})
    res["natal_beta_zones"] = zone
    rng = np.random.default_rng(12)
    n = 200000
    Mpop = rng.uniform(8, 50, n)
    c0pop = np.minimum(rng.beta(0.69, 2.28, n), 0.999)
    res["natal_beta_regimes_T6Gyr"] = {}
    for f in [1e12, 1e14, 1e16, 1e18]:
        _, r = chisr.chi_sr(tab, Mpop, c0pop, 6e9, MU, f)
        res["natal_beta_regimes_T6Gyr"][f"{f:.0e}"] = {str(k): round(float(np.mean(r == k)), 3) for k in range(5)}

    # ---------- 6. C_BH ----------
    cb = {}
    for on in [False, True]:
        si.CBH_ON = on
        fup = f_limit(19.6, 0.72, 0.69, 5e7)
        cb["cbh" if on else "lo"] = {"f_max_real_5e7": fup, "Mmin_chi0.7_d0.02": m_min(fup, 0.7, 0.02),
                                     "Mmin_chi0.7_d0.05": m_min(fup, 0.7, 0.05)}
    si.CBH_ON = False
    res["cbh_effect"] = cb

    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "e12.json"), "w") as fo:
        json.dump(res, fo, indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--ode":
        ode_worker(float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), [float(x) for x in sys.argv[5].split(",")])
    else:
        main()
