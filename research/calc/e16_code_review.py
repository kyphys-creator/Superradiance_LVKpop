"""E16: 計算コードの査読で見つかった問題の修正の前後比較。

修正前の chisr.py は git の 72eac8a のもの(git show で読み込み、chisr_old として使う)。
1. tau -> 0 で chi_after が chi0 を返すか(表の引き方)。
2. 322 から育つ BH の自由場の飽和(BH の質量の減少)。
3. alpha > 0.30 の臨界の帯(211 が超放射する BH の経路)。
4. 臨界の表の補間の誤差。
5. 三準位の ODE(fp3.evolve3)との一致(E14 の固定点の例)と、窓の値。
6. 1 世代目の自然スピン Beta(0.69, 2.28) での印の割合と、計算時間。
出力: out/e16.json
実行: cd research && ~/miniforge3/bin/python3 calc/e16_code_review.py
"""
import importlib.util
import json
import os
import subprocess
import sys
import time
import warnings

import numpy as np

warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(__file__))
import chisr  # noqa: E402
import fp3  # noqa: E402
import sr  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
MU = 0.71e-12
OLD_COMMIT = "72eac8a"


def load_old():
    """修正前の chisr(72eac8a)を別の名前で読み込む。"""
    src = subprocess.run(["git", "show", f"{OLD_COMMIT}:research/calc/chisr.py"], cwd=HERE, capture_output=True,
                         text=True, check=True).stdout
    path = os.path.join(OUT, "_chisr_old.py")
    with open(path, "w") as fo:
        fo.write(src)
    spec = importlib.util.spec_from_file_location("chisr_old", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    os.remove(path)
    return mod


def roundtrip(mod, tab, alphas):
    worst = []
    for a in alphas:
        c2 = float(mod.chic_vec(a, 2, 3, tab))
        xs = np.linspace(c2 + 1e-6, 0.99, 300)
        d = np.abs(mod.chi_after(tab, np.full_like(xs, a), xs, 1e-300) - xs)
        worst.append(float(d.max()))
    return dict(zip([str(a) for a in alphas], worst))


def main():
    res = {}
    old = load_old()
    t_old = old.build()
    t_new = chisr.build()

    # 1. tau -> 0
    al = [0.02, 0.05, 0.1, 0.15, 0.2, 0.25, 0.29]
    res["roundtrip_old"] = roundtrip(old, t_old, al)
    res["roundtrip_new"] = roundtrip(chisr, t_new, al)

    # 2. 322 から育つ BH の自由場
    Ms, cs = np.array([20.0, 30.0, 45.0, 50.0]), np.array([0.3, 0.5, 0.7, 0.6])
    co, _ = old.chi_sr(t_old, Ms, cs, 1e10, MU, 1e19)
    cn, rn = chisr.chi_sr(t_new, Ms, cs, 1e10, MU, 1e19)
    res["free322"] = [{"M": float(M), "chi0": float(c), "old": float(x), "new": float(y), "regime_new": int(r)}
                      for M, c, x, y, r in zip(Ms, cs, co, cn, rn)]

    # 3. alpha > 0.30 の臨界の帯
    band = []
    for a in [0.35, 0.42]:
        M = a / sr.alpha_of(1.0, 1e-12)
        z = chisr._zero(a, 2, 1, 1)
        c_an = float(old.chic_vec(a, 1, 2, t_old))
        x = 0.5 * (z + c_an)
        for f in [1e17, 3e17, 1e19]:
            xo, ro = old.chi_sr(t_old, M, x, 1e10, 1e-12, f)
            xn, rnn = chisr.chi_sr(t_new, M, x, 1e10, 1e-12, f)
            band.append({"alpha": a, "chi0": x, "g211": float(chisr.g_vec(a, x, 2, 1, 1)), "f": f,
                         "old": [float(xo), int(ro)], "new": [float(xn), int(rnn)]})
    res["alpha_gt_0.3_band"] = band

    # 4. 臨界の表の補間の誤差
    err = 0.0
    for a in np.geomspace(0.005, 0.5, 120):
        for m, (n, l) in ((1, (2, 1)), (2, (3, 2))):
            z = chisr._zero(a, n, l, m)
            if z < 0.999:
                err = max(err, abs(float(chisr.chic_vec(a, m, tab=t_new)) - z))
    res["threshold_interp_max_err"] = err

    # 5. ODE と窓
    rows, worst = [], 0.0
    for M, c0, f, Ts in [(19.6, 0.78, 1e16, [1e8, 1e9, 1e10]), (25.0, 0.8, 1e16, [1e7, 1e8, 1e9, 1e10]),
                         (45.0, 0.7, 1e14, [1e8, 1e9, 1e10]), (20.0, 0.3, 1e16, [1e9, 1e10]), (30.0, 0.5, 1e15, [1e8, 1e9, 1e10])]:
        s = fp3.evolve3(M, c0, MU, f, Ts[-1], Ts)
        cnw, _ = chisr.chi_sr(t_new, M, c0, np.array(Ts), MU, f)
        cod, _ = old.chi_sr(t_old, M, c0, np.array(Ts), MU, f)
        for k, T in enumerate(Ts):
            if k < len(s.y[3]):
                rows.append({"M": M, "chi0": c0, "f": f, "T": T, "ode": float(s.y[3][k]), "new": float(cnw[k]), "old": float(cod[k])})
                worst = max(worst, abs(float(cnw[k]) - float(s.y[3][k])))
    res["vs_ode3"] = rows
    res["vs_ode3_max_new"] = worst
    res["vs_ode3_max_old"] = max(abs(r["old"] - r["ode"]) for r in rows)
    res["window"] = {"GW241011_T5e7_f1e12_new": float(chisr.chi_sr(t_new, 19.6, 0.78, 5e7, MU, 1e12)[0]),
                     "GW241011_T5e7_f1e12_old": float(old.chi_sr(t_old, 19.6, 0.78, 5e7, MU, 1e12)[0])}

    # 6. 印の割合と計算時間
    rng = np.random.default_rng(12)
    n = 200000
    Mp = rng.uniform(8, 50, n)
    cp = np.minimum(rng.beta(0.69, 2.28, n), 0.999)
    res["natal_beta_regimes_T6Gyr"] = {}
    for f in [1e12, 1e14, 1e16]:
        _, r = chisr.chi_sr(t_new, Mp, cp, 6e9, MU, f)
        res["natal_beta_regimes_T6Gyr"][f"{f:.0e}"] = {str(k): round(float(np.mean(r == k)), 3) for k in range(5)}
    N = 1_000_000
    M = rng.uniform(5, 80, N)
    c0 = rng.uniform(0.0, 0.99, N)
    T = 10 ** rng.uniform(6, 10, N)
    mus = rng.uniform(0.5e-12, 1.0e-12, N)
    fs = 10 ** rng.uniform(11, 18, N)
    t0 = time.time()
    chisr.chi_sr(t_new, M, c0, T, mus, fs)
    res["timing_1e6_s"] = time.time() - t0

    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "e16.json"), "w") as fo:
        json.dump(res, fo, indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
