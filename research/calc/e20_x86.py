"""E20: 単体版 share/chi_sr を、x86-64 の命令で実際に走らせて確かめる(E19 は x86 の挙動を模しただけだった)。

いまの版 = E19 の版に、f の下限(F_MIN_GEV = 1e-40 GeV。入力の検査)だけを加えたもの。値を計算する所は E19 と同じ。

Apple Silicon の Mac の Rosetta 2(x86-64 の命令の変換)で、Command Line Tools の Python 3.9 の x86-64 の部分と、PyPI の
x86-64 版の numpy・scipy を使う。wheel は展開するだけで、conda の環境にもシステムにも入れない:
    numpy-2.0.2-cp39-cp39-macosx_10_9_x86_64.whl   21,169,942 バイト
    scipy-1.13.1-cp39-cp39-macosx_10_9_x86_64.whl  39,352,927 バイト
    (https://pypi.org/project/numpy/2.0.2/#files、https://pypi.org/project/scipy/1.13.1/#files。sha256 は out/e20.json)
    mkdir site && cd site && unzip -q ../numpy-2.0.2-...whl && unzip -q ../scipy-1.13.1-...whl

  run SITE [WHEELS] : 入力を作り、修正の前(git の OLD_COMMIT = E18)と後(いまの share/chi_sr)の単体版を、この計算機の Python
                      (arm64)と x86-64 の Python で子プロセスとして呼ぶ。x86-64 は二通り: そのまま(numpy が見つける SIMD は
                      SSE4.2 まで)と、環境変数 ROSETTA_ADVERTISE_AVX を立てたもの(AVX2 と FMA3 まで。numpy は別の SIMD の
                      ループを選ぶ)。数えるもの:
                      (a) 例外で止まるか(E18 の版は x86-64 で IndexError になるはず)、
                      (b) いまの版の値と印が arm64 と x86-64 で同じか(警告は例外にする)、
                      (c) share/chi_sr/selftest.py が x86-64 で通るか。
                      WHEELS(wheel を置いたフォルダ)を渡すと、ファイル名・大きさ・sha256 も記録する。
  extras            : run の最後の部分(x86-64 の二通りどうしの比較と、f の小さい極限の例: E19 の版といまの版)だけをやり直す。
  lenient SITE      : run の最後の部分(警告を例外にすると止まった組を、警告を無視して呼び直す)だけをやり直す。
  versions PY ...   : ほかの版の numpy・scipy の入った Python(この計算機にある環境)で share/chi_sr/selftest.py を走らせ、
                      out/e20.json の selftest_other_versions に書き足す(run のあとで)。
  worker ...        : 子プロセス。Python 3.9 と numpy 2.0 でも動くように、単体版と numpy だけを読む。

入力の組:
  e19a_*、e19b_* : E19 の usability(x86 の挙動を模した確かめ)と同じ乱数の列。a は E19 で修正の前の版に、b は後の版に渡した
                   20 万個ずつ(M 5--60 Msun、chi0 0--0.99、年齢 1e5--1e10 yr、mu 0.6--1.0e-12 eV、f は名前の範囲)。
  wide           : 値の表の全体とその外(M 1--316 Msun、chi0 0--0.999、年齢 10^4.5--10^10.3 yr、mu 10^-13.5--10^-10.5 eV、
                   f 1e8--10^20.5 GeV)、30 万個。
  val_<seed>     : E19 の照合(e19_validate.py eval)の入力、4,800 点 × 200 年齢。保存してある arm64 の出力
                   (out/e19_validate_<seed>.npz、値は float32)とも比べる。いまの版だけ。
  edge_*、absurd_* : 端の値。e19 の組と同じ範囲の 2,000 個の入力のうち、一つの引数だけを EDGE(年齢 0、chi0 = 0 や 1 の
                   すぐ下、表のずっと外の M・mu・f。f は下限 F_MIN_GEV = 1e-40 GeV まで)か ABSURD(1e-300、1e300 など)の値に
                   置き換える。組ごとに 1 回呼ぶ。いまの版だけ。absurd の f = 1e-300 は入力の検査で ValueError になるはず。
  box、boxlow    : 引数を同時に広い範囲に振った 30 万個ずつ(対数で一様)。M 1e-3--1e6 Msun、mu 1e-20--1e-5 eV、
                   年齢 1e-3--1e15 yr、chi0 0--1。f は box が 1e-10--1e30 GeV、boxlow が 1e-40--1e-10 GeV。いまの版だけ。
出力: out/e20.json
実行: cd research && ~/miniforge3/envs/python313/bin/python calc/e20_x86.py run /path/to/site [/path/to/wheels]
      cd research && ~/miniforge3/envs/python313/bin/python calc/e20_x86.py versions ~/miniforge3/envs/qft/bin/python ...
"""
import json
import os
import subprocess
import sys
import time
import warnings

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
PKG = os.path.abspath(os.path.join(HERE, "..", "..", "share", "chi_sr"))
TABLES = os.path.join(PKG, "chi_sr_tables.npz")
OLD_COMMIT = "7846eae"
E19_COMMIT = "fd650ff"
X86 = ["arch", "-x86_64", "/usr/bin/python3"]
CHUNK = 200000
SEEDS = (1901, 1902, 1903)
ARGS = ("M", "chi0", "T", "mu", "f")
EDGE = {"T": (0.0, 1e-30, 1e-5, 1.0, 1e12, 1e20, 1e100),
        "chi0": (0.0, 1e-12, 0.995, 0.9951, 0.999, 0.999999, 1 - 1e-12, 1 - 2.0 ** -53),
        "M": (1e-6, 1e-2, 0.5, 1e3, 1e6, 1e12),
        "mu": (1e-30, 1e-20, 1e-15, 1e-9, 1e-5, 1.0),
        "f": (1e-40, 1e-30, 1e-10, 1.0, 1e5, 1e25, 1e40, 1e100)}
ABSURD = {"T": (5e-324, 1e-300, 1e300), "chi0": (5e-324, 1e-300), "M": (1e-300, 1e300), "mu": (1e-300, 1e300),
          "f": (1e-300, 1e300)}


# ---------------------------------------------------------------- 子プロセス
def worker(pkgdir, tables, inp, outp, strict, site):
    if site:
        sys.path.insert(0, site)
    import platform

    import numpy as np
    import scipy
    sys.path.insert(0, pkgdir)
    import chi_sr

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cast = int(np.floor(np.array([np.nan])).astype(int)[0])
    tr = subprocess.run(["sysctl", "-n", "sysctl.proc_translated"], capture_output=True, text=True).stdout.strip()
    env = {"python": sys.version.split()[0], "machine": platform.machine(), "numpy": np.__version__,
           "scipy": scipy.__version__, "rosetta_translated": tr == "1", "nan_to_int64": cast,
           "ROSETTA_ADVERTISE_AVX": os.environ.get("ROSETTA_ADVERTISE_AVX"),
           "simd": [k for k, v in np._core._multiarray_umath.__cpu_features__.items() if v]}
    t0 = time.time()
    model = chi_sr.ChiSR(tables)
    env["load_s"] = round(time.time() - t0, 2)
    z = np.load(inp)
    names = sorted({k.split("/")[0] for k in z.files})
    res, arrs = {}, {}
    for nm in names:
        a = [z[nm + "/" + k] for k in ARGS]
        n = len(a[0])
        chi = np.full(n, np.nan)
        reg = np.full(n, -1, np.int8)
        t0 = time.time()
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("error" if strict else "ignore")
                for i in range(0, n, CHUNK):
                    c, r = model(*[x[i:i + CHUNK] for x in a])
                    chi[i:i + CHUNK] = c
                    reg[i:i + CHUNK] = r
            res[nm] = {"result": "ok", "n": int(n), "call_s": round(time.time() - t0, 2)}
            arrs[nm + "/chi"] = chi
            arrs[nm + "/regime"] = reg
        except Exception as e:                          # E18 の版が x86-64 で止まる所
            res[nm] = {"result": type(e).__name__ + ": " + str(e)[:80], "n": int(n)}
    np.savez(outp, **arrs)
    print(json.dumps({"env": env, "sets": res}))


# ---------------------------------------------------------------- 親
def make_inputs(path):
    import numpy as np
    d = {}
    rng = np.random.default_rng(1)                      # E19 の usability と同じ列(a: 前の版に、b: 後の版に渡した標本)
    n = 200000
    for grp in ("e19a", "e19b"):
        for name, (flo, fhi) in {"f1e10-1e19": (10, 19), "f1e11-1e13": (11, 13), "f1e17-1e20": (17, 20)}.items():
            M = rng.uniform(5, 60, n)
            c0 = rng.uniform(0, 0.99, n)
            T = 10 ** rng.uniform(5, 10, n)
            mu = rng.uniform(0.6e-12, 1.0e-12, n)
            f = 10 ** rng.uniform(flo, fhi, n)
            d[f"{grp}_{name}"] = (M, c0, T, mu, f)
    rng = np.random.default_rng(2001)
    n = 300000
    d["wide"] = (10 ** rng.uniform(0, 2.5, n), rng.uniform(0, 0.999, n), 10 ** rng.uniform(4.5, 10.3, n),
                 10 ** rng.uniform(-13.5, -10.5, n), 10 ** rng.uniform(8, 20.5, n))
    rng = np.random.default_rng(2002)
    n = 2000
    for grp, tab in (("edge", EDGE), ("absurd", ABSURD)):
        for k, key in enumerate(ARGS):
            for i, val in enumerate(tab[key]):
                a = [rng.uniform(5, 60, n), rng.uniform(0, 0.99, n), 10 ** rng.uniform(5, 10, n),
                     rng.uniform(0.6e-12, 1.0e-12, n), 10 ** rng.uniform(10, 19, n)]
                a[k] = np.full(n, val)
                d[f"{grp}_{key}_{i}"] = tuple(a)
    rng = np.random.default_rng(2003)
    n = 300000
    for name, (flo, fhi) in (("box", (-10, 30)), ("boxlow", (-40, -10))):
        d[name] = (10 ** rng.uniform(-3, 6, n), rng.uniform(0, 1, n), 10 ** rng.uniform(-3, 15, n),
                   10 ** rng.uniform(-20, -5, n), 10 ** rng.uniform(flo, fhi, n))
    stored = {}
    for seed in SEEDS:
        z = np.load(os.path.join(OUT, f"e19_validate_{seed}.npz"))
        nT = len(z["T"])
        rep = lambda x: np.repeat(np.asarray(x, float), nT)                      # noqa: E731
        d[f"val_{seed}"] = (rep(z["M"]), rep(z["chi0"]), np.tile(z["T"], len(z["M"])), rep(z["mu"]), rep(z["f"]))
        stored[f"val_{seed}"] = (z["chisr"].ravel(), z["regime"].ravel())
    np.savez(path, **{f"{nm}/{k}": v for nm, a in d.items() for k, v in zip(ARGS, a)})
    return {nm: len(a[0]) for nm, a in d.items()}, stored


def subset(path, sub, keep):
    """keep(名前) が真の組だけを入れた入力のファイルを作る。"""
    import numpy as np
    z = np.load(path)
    np.savez(sub, **{k: z[k] for k in z.files if keep(k.split("/")[0])})


def call(py, pkgdir, inp, outp, strict, site, env=None):
    cmd = list(py) + ["-s", "-B", os.path.abspath(__file__), "worker", pkgdir, TABLES, inp, outp, str(int(strict)), site]
    r = subprocess.run(cmd, capture_output=True, text=True, env=env)
    if r.returncode != 0:
        raise RuntimeError(f"worker failed: {' '.join(cmd)}\n{r.stderr[-2000:]}")
    return json.loads(r.stdout.strip().splitlines()[-1])


def compare(fa, fb, names):
    import numpy as np
    za, zb = np.load(fa), np.load(fb)
    out = {}
    for nm in names:
        if nm + "/chi" not in za.files or nm + "/chi" not in zb.files:
            continue
        ca, cb, ra, rb = za[nm + "/chi"], zb[nm + "/chi"], za[nm + "/regime"], zb[nm + "/regime"]
        fin = np.isfinite(ca) & np.isfinite(cb)
        out[nm] = {"n": int(len(ca)), "max_abs_dchi": float(np.max(np.abs(ca[fin] - cb[fin]))) if fin.any() else 0.0,
                   "n_chi_not_bit_identical": int(np.sum(~((ca == cb) | (np.isnan(ca) & np.isnan(cb))))),
                   "nan_pattern_same": bool(np.array_equal(np.isnan(ca), np.isnan(cb))),
                   "n_nan": int(np.isnan(ca).sum()), "regime_differs": int(np.sum(ra != rb)),
                   "regime_share": [round(float(np.mean(ra == k)), 4) for k in range(6)]}
    return out


def total(cmp, sets, pick):
    """compare の結果を、名前が pick を満たす組について足し合わせる(sets: その版・その計算機の組ごとの結果)。"""
    nm = [k for k in cmp if pick(k)]
    return {"n_sets": len(nm), "n": int(sum(cmp[k]["n"] for k in nm)),
            "max_abs_dchi": max([cmp[k]["max_abs_dchi"] for k in nm], default=0.0),
            "n_chi_not_bit_identical": int(sum(cmp[k]["n_chi_not_bit_identical"] for k in nm)),
            "regime_differs": int(sum(cmp[k]["regime_differs"] for k in nm)),
            "nan_pattern_same": all(cmp[k]["nan_pattern_same"] for k in nm),
            "n_nan": int(sum(cmp[k]["n_nan"] for k in nm)),
            "n_sets_raised": int(sum(1 for k, v in sets.items() if pick(k) and v["result"] != "ok"))}


def arch_list(site):
    env0 = {k: v for k, v in os.environ.items() if k != "ROSETTA_ADVERTISE_AVX"}
    return (("arm64", [sys.executable], "", env0), ("x86_64", X86, site, env0),
            ("x86_64_avx2", X86, site, dict(env0, ROSETTA_ADVERTISE_AVX="1")))


def lenient(site, res):
    """警告を例外にすると止まる組(極端な値)を、警告を無視して呼び直す: x86-64 で例外にならないか、値と印が同じか。"""
    work = os.path.join(OUT, "work", "e20")
    f_ = lambda arch: os.path.join(work, f"out_lenient_{arch}.npz")              # noqa: E731
    archs = arch_list(os.path.abspath(site))
    raised = sorted({k for a, *_ in archs for k, v in res["new"][a].items() if v["result"] != "ok"})
    out = {"sets": raised, "strict_result": {k: res["new"]["arm64"][k]["result"] for k in raised}}
    if raised:
        sub = os.path.join(work, "inputs_lenient.npz")
        subset(os.path.join(work, "inputs.npz"), sub, lambda nm: nm in raised)
        runs = {arch: call(py, PKG, sub, f_(arch), False, st, env) for arch, py, st, env in archs}
        out["result"] = {arch: {k: v["result"] for k, v in runs[arch]["sets"].items()} for arch in runs}
        import numpy as np
        z, zi = np.load(f_("arm64")), np.load(sub)
        out["arm64_output"] = {k: {"n": int(len(z[k + "/chi"])), "n_nan": int(np.isnan(z[k + "/chi"]).sum()),
                                   "regime_share": [round(float(np.mean(z[k + "/regime"] == j)), 4) for j in range(6)],
                                   "max_abs_chi_minus_chi0": float(np.max(np.abs(z[k + "/chi"] - zi[k + "/chi0"])))}
                               for k in raised if k + "/chi" in z.files}
        for arch, *_ in archs[1:]:
            out[f"arm64_vs_{arch}"] = compare(f_("arm64"), f_(arch), raised)
    res["new_lenient"] = out
    return out


def extras(res):
    """x86-64 の二通りどうしの比較と、f の小さい極限の例(E19 の版 = git の E19_COMMIT と、いまの版)。"""
    import importlib.util
    import numpy as np
    work = os.path.join(OUT, "work", "e20")
    names = sorted(res["n_inputs"])
    pair = compare(os.path.join(work, "out_new_x86_64.npz"), os.path.join(work, "out_new_x86_64_avx2.npz"), names)
    if res.get("new_lenient", {}).get("sets"):
        pair.update(compare(os.path.join(work, "out_lenient_x86_64.npz"), os.path.join(work, "out_lenient_x86_64_avx2.npz"),
                            res["new_lenient"]["sets"]))
    res["new_x86_64_vs_x86_64_avx2"] = {"n_sets": len(pair), "n": int(sum(v["n"] for v in pair.values())),
                                        "max_abs_dchi": max(v["max_abs_dchi"] for v in pair.values()),
                                        "n_chi_not_bit_identical": int(sum(v["n_chi_not_bit_identical"] for v in pair.values())),
                                        "regime_differs": int(sum(v["regime_differs"] for v in pair.values()))}
    e19pkg = os.path.join(work, "e19_pkg")
    os.makedirs(e19pkg, exist_ok=True)
    src = subprocess.run(["git", "show", f"{E19_COMMIT}:share/chi_sr/chi_sr.py"], cwd=HERE, capture_output=True, text=True,
                         check=True).stdout
    with open(os.path.join(e19pkg, "chi_sr.py"), "w") as fo:
        fo.write(src)
    mods = {}
    for tag, d in (("e19", e19pkg), ("now", PKG)):
        spec = importlib.util.spec_from_file_location(f"chi_sr_{tag}", os.path.join(d, "chi_sr.py"))
        mods[tag] = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mods[tag])
    models = {tag: m.ChiSR(TABLES) for tag, m in mods.items()}
    M, c0, T, mu = 30.0, 0.9, 1e9, 0.8e-12
    out = {"M": M, "chi0": c0, "T_yr": T, "mu_eV": mu, "alpha": float(models["now"].alpha(np.array([M]), np.array([mu]))[0]),
           "f_overflow_GeV": float(models["now"].MPL_GEV / np.finfo(float).max ** 0.25), "F_MIN_GEV": mods["now"].F_MIN_GEV,
           "rows": []}
    for f in (1e20, 1e12, 1e-10, 1e-40, 1e-57, 1e-58, 1e-100, 1e-300):
        row = {"f_GeV": f}
        for tag, model in models.items():
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                try:
                    chi, reg = model(M, c0, T, mu, f)
                    row[tag] = {"chi": float(chi), "regime": int(reg), "warnings": sorted({str(x.message) for x in w})}
                except ValueError as e:
                    row[tag] = {"error": "ValueError: " + str(e)}
        out["rows"].append(row)
    # 下限 F_MIN_GEV のところで、alpha <= ALPHA_MAX の点の率があふれるのに要る mu(r_+ <= 2 GM)。余裕の目安
    m = models["now"]
    r4 = (m.MPL_GEV / mods["now"].F_MIN_GEV) ** 4
    big = np.finfo(float).max
    am = float(m.ALPHA_MAX)
    out["mu_eV_overflow_at_F_MIN"] = float(min(big / (m.KAPPA_BH * am ** 11 * r4 * 2.0), big / (m.KAPPA_INF * am ** 8 * r4),
                                               big / (m.K_C * am ** 11 * r4 * 2.0), big / (m.K_D * am ** 8 * r4)) * m.HBAR_EVS)
    res["small_f"] = out
    # 警告(あふれ)が出る広い箱(boxlow): 自己相互作用の率があふれる点はどこか、その点の印は何か
    ov = {}
    zi = np.load(os.path.join(work, "inputs.npz"))
    lz = os.path.join(work, "out_lenient_arm64.npz")
    if os.path.exists(lz) and "boxlow/regime" in np.load(lz).files:
        m = models["now"]
        Mb, cb, mub, fb = (zi["boxlow/" + k] for k in ("M", "chi0", "mu", "f"))
        a = m.alpha(Mb, mub)
        with np.errstate(all="ignore"):
            r4 = (m.MPL_GEV / fb) ** 4
            gbh = m.KAPPA_BH * a ** 11 * r4 * m.rplus(np.minimum(cb, 0.999)) * (mub / m.HBAR_EVS)
            ginf = m.KAPPA_INF * a ** 8 * r4 * (mub / m.HBAR_EVS)
        bad = ~(np.isfinite(gbh) & np.isfinite(ginf))
        reg = np.load(lz)["boxlow/regime"]
        ov = {"n": int(len(a)), "n_rate_overflow": int(bad.sum()), "alpha_min_overflow": float(a[bad].min()) if bad.any() else None,
              "alpha_max_no_overflow": float(a[~bad].max()), "regimes_at_overflow": sorted(int(x) for x in set(reg[bad])),
              "ALPHA_MAX": float(m.ALPHA_MAX)}
    res["boxlow_overflow"] = ov
    return {"new_x86_64_vs_x86_64_avx2": res["new_x86_64_vs_x86_64_avx2"], "small_f": out, "boxlow_overflow": ov}


def run(site, wheels=None):
    import hashlib

    import numpy as np
    site = os.path.abspath(site)
    work = os.path.join(OUT, "work", "e20")
    oldpkg = os.path.join(work, "old_pkg")
    os.makedirs(oldpkg, exist_ok=True)
    src = subprocess.run(["git", "show", f"{OLD_COMMIT}:share/chi_sr/chi_sr.py"], cwd=HERE, capture_output=True, text=True,
                         check=True).stdout
    with open(os.path.join(oldpkg, "chi_sr.py"), "w") as fo:
        fo.write(src)
    inp = os.path.join(work, "inputs.npz")
    inp_old = os.path.join(work, "inputs_old.npz")
    sizes, stored = make_inputs(inp)
    main = lambda nm: nm.startswith(("e19a_", "e19b_", "wide"))                 # noqa: E731
    subset(inp, inp_old, main)
    res = {"old_commit": OLD_COMMIT, "e19_commit": E19_COMMIT,          # いまの版は、ファイルの sha256 で特定する
           "new_chi_sr_py_sha256": hashlib.sha256(open(os.path.join(PKG, "chi_sr.py"), "rb").read()).hexdigest(),
           "tables_sha256": hashlib.sha256(open(TABLES, "rb").read()).hexdigest(),
           "n_inputs": sizes, "edge_values": {"edge": EDGE, "absurd": ABSURD}}
    if wheels:
        res["wheels"] = []
        for fn in sorted(os.listdir(wheels)):
            if fn.endswith(".whl"):
                with open(os.path.join(wheels, fn), "rb") as fo:
                    b = fo.read()
                res["wheels"].append({"file": fn, "bytes": len(b), "sha256": hashlib.sha256(b).hexdigest()})
    archs = arch_list(site)
    xs = [a[0] for a in archs[1:]]
    f_ = lambda ver, arch: os.path.join(work, f"out_{ver}_{arch}.npz")           # noqa: E731
    runs = {}
    for ver, pkgdir, strict, ifile in (("old", oldpkg, False, inp_old), ("new", PKG, True, inp)):
        for arch, py, st, env in archs:
            t0 = time.time()
            runs[ver, arch] = call(py, pkgdir, ifile, f_(ver, arch), strict, st, env)
            print(f"--- {ver} {arch} ({time.time() - t0:.0f} s)", flush=True)
    res["env"] = {arch: runs["new", arch]["env"] for arch, *_ in archs}
    for ver in ("old", "new"):
        res[ver] = {arch: runs[ver, arch]["sets"] for arch, *_ in archs}
    names = sorted(sizes)
    val = lambda nm: nm.startswith("val_")                                      # noqa: E731
    edge = lambda nm: nm.startswith(("edge_", "absurd_"))                       # noqa: E731
    box = lambda nm: nm.startswith("box")                                       # noqa: E731
    res["summary"] = {}
    for x in xs:
        cn = compare(f_("new", "arm64"), f_("new", x), names)
        co = compare(f_("old", "arm64"), f_("old", x), names)
        res[f"new_arm64_vs_{x}"] = cn
        res[f"old_arm64_vs_{x}"] = co
        res["summary"][x] = {"new_main": total(cn, res["new"][x], main), "new_val": total(cn, res["new"][x], val),
                             "new_edge": total(cn, res["new"][x], edge), "new_box": total(cn, res["new"][x], box),
                             "old_main": total(co, res["old"][x], main)}
        # E19 の照合で保存した出力(arm64、値は float32)と、x86-64 の出力
        zx = np.load(f_("new", x))
        d = {}
        for nm, (c_st, r_st) in stored.items():
            cx, rx = zx[nm + "/chi"], zx[nm + "/regime"]
            fin = np.isfinite(c_st) & np.isfinite(cx)
            d[nm] = {"n": int(len(cx)), "regime_differs": int(np.sum(rx != r_st)),
                     "n_float32_not_identical": int(np.sum(~((cx.astype(np.float32) == c_st) | (np.isnan(cx) & np.isnan(c_st))))),
                     "max_abs_dchi": float(np.max(np.abs(cx[fin] - c_st[fin])))}
        res[f"{x}_vs_stored_validation"] = d
    res["new_vs_old_arm64"] = compare(f_("new", "arm64"), f_("old", "arm64"), names)
    res["summary"]["arm64"] = {"new_sets_raised": {k: v["result"] for k, v in res["new"]["arm64"].items() if v["result"] != "ok"},
                               "old_sets_raised": {k: v["result"] for k, v in res["old"]["arm64"].items() if v["result"] != "ok"}}
    for x in xs:
        res["summary"][x]["new_sets_raised"] = {k: v["result"] for k, v in res["new"][x].items() if v["result"] != "ok"}
        res["summary"][x]["old_sets_raised"] = {k: v["result"] for k, v in res["old"][x].items() if v["result"] != "ok"}
    # 自己診断
    res["selftest"] = {}
    for arch, py, st, env in archs:
        e = dict(env, PYTHONPATH=st) if st else env
        r = subprocess.run(list(py) + ["-s", "-B", os.path.join(PKG, "selftest.py")], capture_output=True, text=True, env=e)
        res["selftest"][arch] = {"returncode": r.returncode, "stdout": r.stdout.strip().splitlines(),
                                 "stderr_tail": r.stderr.strip().splitlines()[-3:]}
    lenient(site, res)
    extras(res)
    with open(os.path.join(OUT, "e20.json"), "w") as fo:
        json.dump(res, fo, indent=1)
    print(json.dumps({k: res[k] for k in ("env", "summary", "selftest", "new_lenient", "small_f")}, indent=1))


def versions(pythons):
    path = os.path.join(OUT, "e20.json")
    res = json.load(open(path))
    out = []
    for py in pythons:
        v = subprocess.run([py, "-s", "-c", "import sys, platform, numpy, scipy; print(sys.version.split()[0], "
                            "platform.machine(), numpy.__version__, scipy.__version__)"], capture_output=True, text=True,
                           check=True).stdout.split()
        r = subprocess.run([py, "-s", "-B", os.path.join(PKG, "selftest.py")], capture_output=True, text=True)
        out.append({"python": v[0], "machine": v[1], "numpy": v[2], "scipy": v[3], "returncode": r.returncode,
                    "stdout": r.stdout.strip().splitlines()})
        print(json.dumps(out[-1], indent=1))
    res["selftest_other_versions"] = out
    with open(path, "w") as fo:
        json.dump(res, fo, indent=1)


if __name__ == "__main__":
    if sys.argv[1] == "worker":
        worker(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6] == "1", sys.argv[7] if len(sys.argv) > 7 else "")
    elif sys.argv[1] == "run":
        run(sys.argv[2], sys.argv[3] if len(sys.argv) > 3 else None)
    elif sys.argv[1] == "versions":
        versions(sys.argv[2:])
    elif sys.argv[1] == "extras":                       # run のあと、この部分だけをやり直す
        _p = os.path.join(OUT, "e20.json")
        _res = json.load(open(_p))
        print(json.dumps(extras(_res), indent=1, ensure_ascii=False))
        with open(_p, "w") as _fo:
            json.dump(_res, _fo, indent=1)
    elif sys.argv[1] == "lenient":                      # run のあと、この部分だけをやり直す
        _p = os.path.join(OUT, "e20.json")
        _res = json.load(open(_p))
        print(json.dumps(lenient(sys.argv[2], _res), indent=1))
        with open(_p, "w") as _fo:
            json.dump(_res, _fo, indent=1)
    else:
        sys.exit(__doc__)
