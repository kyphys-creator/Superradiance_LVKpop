"""share/chi_sr/chi_sr_tables.npz の表をヒートマップにして、1 つの PDF にまとめる(見るための図。ノートの図ではない)。

ページ:
  1. 目次と、各配列の意味
  2. 固定点の表: log10 I、log10 J、log10 J422 を (chi, alpha) の面で
  3. 相対論的な率の表: F(211)、F(322) と、臨界スピン thr1--thr3
  4. 値の表、211 から育つ節(s = 1.25): 211 の臨界・322 の臨界に着く時刻、届いた進み具合
  5. 値の表、322 から育つ節(s = 0.25): 322 の臨界まで半分・臨界に着く時刻、届いた進み具合
  6. 値の表、f = 1e12 GeV: 初期スピンの位置 s と alpha の面で、最初の段階を終える時刻と届いた進み具合
  7. 値の表の診断: 解けなかった節の割合、曲がりの指標 KI の最大、戻りの最大
表の値の読み方(log10 tau の符号の戻し方)は share/chi_sr/chi_sr.py の ChiSR をそのまま使う。
出力: out/chi_sr_tables_heatmaps.pdf
実行: cd research && ~/miniforge3/envs/python313/bin/python calc/plot_chisr_tables.py   (約 20 秒)
"""
import os
import subprocess
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402
from matplotlib.colors import FuncNorm, LinearSegmentedColormap, Normalize  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..")
TAB = os.path.join(ROOT, "share", "chi_sr", "chi_sr_tables.npz")
OUT = os.path.join(HERE, "out", "chi_sr_tables_heatmaps.pdf")
STYLE = os.path.expanduser("~/Dropbox/GraduateResearch/physrev_mplstyle-main/physrev.mplstyle")
sys.path.insert(0, os.path.join(ROOT, "share", "chi_sr"))
from chi_sr import ChiSR  # noqa: E402

plt.style.use(STYLE)
plt.rcParams["legend.frameon"] = False
plt.rcParams["font.family"] = ["DejaVu Sans", "Arial Unicode MS"]   # 日本語は TrueType の Arial Unicode MS で補う(Hiragino は CFF で、Type 42 に埋めると PDF の字体の型が食い違う)
plt.rcParams["axes.titlesize"] = 9
plt.rcParams["pdf.compression"] = 9
plt.rcParams["pdf.fonttype"] = 42                               # 文字を PDF の中で検索できるように

# 一つの色相の連続(薄い = 小さい、濃い = 大きい)と、灰色を中心にした二色の発散
SEQ = LinearSegmentedColormap.from_list("seq", ["#f1f4f9", "#a9bedb", "#5d81b4", "#2c4a78", "#14233d"])
DIV = LinearSegmentedColormap.from_list("div", ["#1f4e85", "#7d9fc9", "#e6e6e6", "#e0a17f", "#a3401a"])
INK = "#222222"


def hatch_bg(ax, label=True):
    """NaN(届かない・表の外)の所は斜線で見せる: 軸の背景に斜線の四角を敷き、NaN のセルは透明にする。"""
    ax.add_patch(Rectangle((0, 0), 1, 1, transform=ax.transAxes, facecolor="white", edgecolor="#b0b0b0",
                           hatch="////", lw=0, zorder=-10))


def seq_cmap():
    c = SEQ.copy()
    c.set_bad((0, 0, 0, 0))
    return c


def page_title(fig, text, sub=None):
    fig.suptitle(text, fontsize=11, x=0.02, ha="left")
    if sub:
        fig.text(0.02, 0.005, sub, fontsize=7, color="#555555", ha="left", va="bottom", wrap=True)


m = ChiSR(TAB)
z = np.load(TAB)
rev = subprocess.run(["git", "log", "-1", "--format=%h %ad", "--date=short", "--", TAB], cwd=ROOT,
                     capture_output=True, text=True).stdout.strip()

# ---------- 固定点の表の座標: w を chi に戻す(行ごとに臨界が違うので、2 次元の網) ----------
A = m.ALPHAS
c1 = m.chic(A, 1)[:, None]
c2 = m.chic(A, 2)[:, None]
CHI = m.chi_of_w(m.W[None, :], c1, c2)
AL = np.broadcast_to(A[:, None], CHI.shape)

# ---------- 値の表 ----------
rt = m.rt
AN, LF, SA, SB, MN = rt["A_NODES"], rt["LF_NODES"], rt["SA"], rt["SB"], rt["M_NODES"]
PL = m.P_LEV2
nreach = z["rt_nreach"].astype(int)
PMAX = np.where(nreach > 0, PL[np.clip(nreach - 1, 0, len(PL) - 1)], np.nan)
LT_RANGE = (-14.0, 0.2)
CONT = [-8, -6, -4, -2]                                  # log10 tau の等高線


def lt_at(piece, js, il, p):
    """log10 tau(tau = mu T [eV yr])が進み具合 p を初めて越える時刻、(alpha の節, log10 f の節) の面。"""
    k = int(round(p / (PL[1] - PL[0])))
    v = np.asarray(rt["LT"][piece][:, js, :, il, k], float)
    return np.where(np.isfinite(v), v, np.nan)


def tau_panel(ax, X, Y, V, title, xlog=True, contours=True):
    hatch_bg(ax)
    pc = ax.pcolormesh(X, Y, V.T, cmap=seq_cmap(), vmin=LT_RANGE[0], vmax=LT_RANGE[1], shading="nearest",
                       rasterized=True)
    if contours and np.isfinite(V).sum() > 4:
        cs = ax.contour(X, Y, V.T, levels=CONT, colors="white", linewidths=0.5)
        ax.clabel(cs, fmt=lambda v: f"{v:.0f}", fontsize=6, inline_spacing=2)
    if xlog:
        ax.set_xscale("log")
    ax.set_title(title, loc="left")
    return pc


def pmax_panel(ax, X, Y, V, title, xlog=True):
    hatch_bg(ax)
    pc = ax.pcolormesh(X, Y, V.T, cmap=seq_cmap(), vmin=0, vmax=4, shading="nearest", rasterized=True)
    if xlog:
        ax.set_xscale("log")
    ax.set_title(title, loc="left")
    return pc


def alpha_f_axes(ax, row, col, ncol):
    if row == 1:
        ax.set_xlabel(r"$\alpha$")
    else:
        ax.tick_params(labelbottom=False)
    if col == 0:
        ax.set_ylabel(r"$\log_{10}(f/\mathrm{GeV})$")
    else:
        ax.tick_params(labelleft=False)


with PdfPages(OUT) as pdf:
    # ---------- 1. 目次 ----------
    fig = plt.figure(figsize=(7.0, 8.0))
    fig.text(0.06, 0.95, "chi_sr_tables.npz のヒートマップ", fontsize=14, va="top")
    fig.text(0.06, 0.915, f"表: share/chi_sr/chi_sr_tables.npz(最後に変えたコミット {rev})。"
             f"図: research/calc/plot_chisr_tables.py", fontsize=8, va="top", color="#555555")
    lines = [
        ("2", "固定点の表(ALPHAS 161 × W 1031)",
         ["w は臨界に対するスピンの位置(322 の臨界 = 0、211 の臨界 = 1、χ = 0.995 = 2)。図では χ に戻して描く。",
          "I(χ): 自己相互作用の固定点で降りる時間。I(χ_T) − I(χ_0) = (μT/ħ)(f/M_Pl)²(無次元)。",
          "J(χ): 三準位の外の準位(411、433、211 の臨界より上の 544)の e-fold の累計 × (f/M_Pl)²。J422 は 422 の同じ量。"]),
        ("3", "相対論的な率の表(rel211、rel322)と臨界スピン(thr1–thr3)",
         ["F: 超放射の率を Γ/μ = 2F α^(4ℓ+4) と書いたときの係数(連分数法)。F > 0 で超放射、F < 0 で BH が吸収。",
          "臨界スピン χ_c(α): 率が 0 になるスピン(m = 1, 2, 3)。"]),
        ("4–6", "値の表(rt_*: 211 か 322 から育つ × α 47 × s 25 × log10 f 51 × M 2 × 進み具合 401)",
         ["基準の運動論を 1 個の量子から解き、スピンの「進み具合」p が各準位を初めて越える時刻 log10(μT/eV yr) を持つ。",
          "211 から育つ節: p = 1 で 211 の臨界、p = 2 で 322 の臨界に着く。322 から育つ節: p = 1 で 322 の臨界に着く。",
          "s は初期スピンの位置(w と同じ座標)。211 から育つ節は s = 1 + SB、322 から育つ節は s = SA。",
          "斜線: 届かない(τ ≤ 1.3 eV yr、解けなかった節は解けた所までのうちに、その p を越えない)。",
          "白い等高線は log10 τ = −8, −6, −4, −2。",
          f"目安: μ = 0.71×10⁻¹² eV で T = 10⁷ 年は log10 τ = {np.log10(0.71e-12 * 1e7):.1f}、"
          f"10¹⁰ 年は {np.log10(0.71e-12 * 1e10):.1f}。".replace("-", "−"),
          "M は 1 と 300 M☉ の 2 節(M は ln(N0 M²) を通してだけ効き、間は時刻を線形に混ぜる)。"]),
        ("7", "値の表の診断",
         ["解けなかった節(status ≠ 0: ODE が止まった、または飛ばした)の割合、曲がりの指標 KI の最大 [dex]、",
          "戻りの最大(スピンが包絡から戻る量。b の単位)。いずれも s と M について最大か割合をとる。"]),
    ]
    y = 0.86
    for pg, head, body in lines:
        fig.text(0.06, y, f"p.{pg}", fontsize=9, va="top", color="#5d81b4")
        fig.text(0.14, y, head, fontsize=9, va="top")
        y -= 0.032
        for b in body:
            fig.text(0.14, y, b, fontsize=7.2, va="top", color="#333333")
            y -= 0.026
        y -= 0.02
    fig.text(0.06, 0.05, "値はすべて表のまま(補間していない)。値の表の時刻は節での値で、chi_sr.py が節のあいだで混ぜる前のもの。",
             fontsize=7.2, color="#555555")
    pdf.savefig(fig)
    plt.close(fig)

    # ---------- 2. 固定点の表 ----------
    fig, axs = plt.subplots(1, 3, figsize=(7.0, 3.6), layout="constrained", sharey=True)
    log10I = m.lnI / np.log(10)
    with np.errstate(divide="ignore"):
        lJ = np.where(m.J > 0, np.log10(m.J), np.nan)
        lJ4 = np.where(m.J422 > 0, np.log10(m.J422), np.nan)
    # 色の範囲は値の全体(I)か、0 でない値の 1 パーセント点から最大まで(J、J422)
    rng = lambda V, q: (np.floor(np.nanpercentile(V, q) * 2) / 2, np.ceil(np.nanmax(V) * 10) / 10)
    panels = [(log10I, r"$\log_{10} I$", "(a) 固定点で降りる時間 I", rng(log10I, 0)),
              (lJ, r"$\log_{10} J$", "(b) 外の準位の e-fold の累計 J", rng(lJ, 1)),
              (lJ4, r"$\log_{10} J_{422}$", "(c) 422 の e-fold の累計", rng(lJ4, 1))]
    LINE = "#e09b24"                                      # 臨界の線(青の上で見える physrev の C1)
    for ax, (V, lab, title, (lo, hi)) in zip(axs, panels):
        hatch_bg(ax)
        pc = ax.pcolormesh(CHI, AL, V, cmap=seq_cmap(), vmin=lo, vmax=hi, shading="gouraud", rasterized=True)
        ax.plot(c1[:, 0], A, color=LINE, lw=0.8, ls="--")
        ax.plot(c2[:, 0], A, color=LINE, lw=0.8, ls="-")
        ax.set_yscale("log")
        ax.set_xlim(0, 1)
        ax.set_xlabel(r"$\chi$")
        ax.set_title(title, loc="left")
        cb = fig.colorbar(pc, ax=ax, location="top", shrink=0.9, extend="min" if lo > np.nanmin(V) else "neither",
                          pad=0.02)
        cb.set_label(lab, fontsize=8)
        cb.ax.tick_params(labelsize=7)
    axs[0].set_ylabel(r"$\alpha$")
    axs[0].text(0.04, 0.97, "実線: 322 の臨界\n破線: 211 の臨界", transform=axs[0].transAxes, fontsize=6.5,
                va="top", color=INK, bbox=dict(facecolor="white", edgecolor="none", pad=1.5))
    fig.suptitle("固定点の表(ALPHAS × W)", fontsize=10, x=0.01, ha="left")
    fig.text(0.01, -0.02, "斜線: 表の外(322 の臨界より下)か、J = 0(外の準位がまだ育たない)。"
             "I は 322 の臨界に近づくと発散する。", fontsize=6.5, color="#555555")
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)

    # ---------- 3. 相対論的な率の表と臨界スピン ----------
    fig, axs = plt.subplots(1, 3, figsize=(7.0, 3.6), layout="constrained")
    LT0 = 1e-4                                            # これより小さい |F| は線形

    def make_norm(Fmin, Fmax):
        sn, sp = np.log10(1 + abs(Fmin) / LT0), np.log10(1 + Fmax / LT0)

        def fwd(F):
            F = np.asarray(F, float)
            return np.where(F < 0, -np.log10(1 + np.abs(F) / LT0) / sn, np.log10(1 + np.abs(F) / LT0) / sp)

        def inv(u):
            u = np.asarray(u, float)
            return np.where(u < 0, -LT0 * (10 ** (np.abs(u) * sn) - 1), LT0 * (10 ** (np.abs(u) * sp) - 1))
        return FuncNorm((fwd, inv), vmin=Fmin, vmax=Fmax)

    for ax, key, mm, title in [(axs[0], "211", 1, "(a) 211 の率の係数 F"), (axs[1], "322", 2, "(b) 322 の率の係数 F")]:
        a_ = z[f"rel{key}_alpha"]
        c_ = z[f"rel{key}_chi"]
        F = z[f"rel{key}_F"]
        Fmin, Fmax = float(F.min()), float(F.max())
        pc = ax.pcolormesh(c_, a_, F, cmap=DIV, norm=make_norm(Fmin, Fmax), shading="gouraud", rasterized=True)
        ax.contour(c_, a_, F, levels=[0], colors=INK, linewidths=0.6)
        ax.plot(z[f"thr{mm}"], z["thr_alpha"], color="white", lw=0.5, ls=":")
        lo, hi = z[f"rel{key}_range"]
        for v in (lo, hi):
            ax.axhline(v, color="#555555", lw=0.4, ls="--")
        ax.set_yscale("log")
        ax.set_ylim(a_[0], a_[-1])
        ax.set_xlabel(r"$\chi$")
        ax.set_ylabel(r"$\alpha$")
        ax.set_title(title, loc="left")
        # 目盛り: 0 と、各側の大きい桁から順に、色の帯の上で 0.25 以上離れるものだけ
        norm = pc.norm
        cand = [0.0] + [-10.0 ** k for k in range(0, -7, -1)] + [10.0 ** k for k in range(0, -7, -1)]
        ticks = []
        for t in cand:
            if Fmin <= t <= Fmax and all(abs(float(norm(t)) - float(norm(u))) >= 0.25 for u in ticks):
                ticks.append(t)
        ticks.sort()
        cb = fig.colorbar(pc, ax=ax, location="top", shrink=0.9, pad=0.02, ticks=ticks)
        cb.ax.set_xticklabels([("0" if t == 0 else (f"$-10^{{{int(np.log10(-t))}}}$" if t < 0 else f"$10^{{{int(np.log10(t))}}}$"))
                               for t in ticks])
        cb.set_label(r"$F=\Gamma/(2\mu\alpha^{4\ell+4})$", fontsize=8)
        cb.ax.tick_params(labelsize=7)
    ax = axs[2]
    for mm, ls in [(1, "-"), (2, "--"), (3, "-.")]:
        ax.plot(z["thr_alpha"], z[f"thr{mm}"], color=INK, ls=ls, lw=0.8, label=f"$m={mm}$")
    ax.legend(loc="upper left", fontsize=7)
    ax.set_xscale("log")
    ax.set_ylim(0, 1.02)
    ax.set_xlabel(r"$\alpha$")
    ax.set_ylabel(r"$\chi_c$")
    ax.set_title("(c) 臨界スピン thr1–thr3", loc="left")
    fig.suptitle("相対論的な率の表と臨界スピン", fontsize=10, x=0.01, ha="left")
    fig.text(0.01, -0.02, "黒い線: F = 0。白い点線: 臨界スピンの表 thr(F = 0 と重なるはず)。灰色の破線: 表を使う α の範囲"
             "(外は解析式)。青: BH が吸収、橙: 超放射。", fontsize=6.5, color="#555555")
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)

    # ---------- 4, 5. 値の表、(alpha, log10 f) の面 ----------
    for piece, js, cols, head in [
            (1, 12, [(1.0, "211 の臨界に着く"), (2.0, "322 の臨界に着く")], f"211 から育つ節(s = 1 + {SB[12]:.2f})"),
            (0, 12, [(0.5, "322 の臨界まで半分"), (1.0, "322 の臨界に着く")], f"322 から育つ節(s = {SA[12]:.2f})")]:
        fig, axs = plt.subplots(2, 3, figsize=(7.0, 5.6), layout="constrained", sharex=True, sharey=True)
        for r, il in enumerate((0, 1)):
            for c, (p, lab) in enumerate(cols):
                pc = tau_panel(axs[r, c], AN, LF, lt_at(piece, js, il, p),
                               f"({'abcdef'[3 * r + c]}) p = {p:g}: {lab}\n$M={MN[il]:g}\\,M_\\odot$")
                alpha_f_axes(axs[r, c], r, c, 3)
            pcm = pmax_panel(axs[r, 2], AN, LF, PMAX[piece, :, js, :, il],
                             f"({'abcdef'[3 * r + 2]}) 届いた進み具合\n$M={MN[il]:g}\\,M_\\odot$")
            alpha_f_axes(axs[r, 2], r, 2, 3)
        cb = fig.colorbar(pc, ax=axs[:, :2], location="bottom", shrink=0.7, pad=0.02, extend="min")
        # 日本語と数式は別の行に(同じ行に $ があると、行全体が数式の字体になり日本語が出ない)
        cb.set_label("その進み具合を初めて越える時刻\n" r"$\log_{10}(\mu T/\mathrm{eV\,yr})$", fontsize=8)
        cb2 = fig.colorbar(pcm, ax=axs[:, 2], location="bottom", shrink=0.9, pad=0.02)
        cb2.set_label("τ ≤ 1.3 eV yr のうちに届いた進み具合\n" r"$p_\mathrm{max}$", fontsize=8)
        fig.suptitle(f"値の表: {head}", fontsize=10, x=0.01, ha="left")
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

    # ---------- 6. 値の表、f = 1e12 GeV、(alpha, s) の面 ----------
    kf = int(np.argmin(np.abs(LF - 12.0)))
    S_ALL = np.concatenate([SA, 1 + SB])
    fig, axs = plt.subplots(2, 2, figsize=(7.0, 5.6), layout="constrained", sharex=True, sharey=True)
    for r, il in enumerate((0, 1)):
        V = np.concatenate([np.stack([lt_at(0, j, il, 1.0)[:, kf] for j in range(len(SA))], 1),
                            np.stack([lt_at(1, j, il, 1.0)[:, kf] for j in range(len(SB))], 1)], 1)
        P = np.concatenate([PMAX[0, :, :, kf, il], PMAX[1, :, :, kf, il]], 1)
        pc = tau_panel(axs[r, 0], AN, S_ALL, V, f"({'ac'[r]}) p = 1: 最初の段階を終える(近いほうの臨界に着く)\n"
                       f"$M={MN[il]:g}\\,M_\\odot$")
        pcm = pmax_panel(axs[r, 1], AN, S_ALL, P, f"({'bd'[r]}) 届いた進み具合\n$M={MN[il]:g}\\,M_\\odot$")
        for ax in axs[r]:
            ax.axhline(1.0, color=INK, lw=0.5, ls="--")
        axs[r, 0].set_ylabel("初期スピンの位置 s")
        if r == 1:
            for ax in axs[r]:
                ax.set_xlabel(r"$\alpha$")
    box = dict(facecolor="white", edgecolor="none", pad=1.0)
    axs[0, 0].text(0.03, 0.53, "211 から育つ", transform=axs[0, 0].transAxes, fontsize=6.5, va="bottom", bbox=box)
    axs[0, 0].text(0.03, 0.47, "322 から育つ", transform=axs[0, 0].transAxes, fontsize=6.5, va="top", bbox=box)
    cb = fig.colorbar(pc, ax=axs[:, 0], location="bottom", shrink=0.9, pad=0.02, extend="min")
    cb.set_label(r"$\log_{10}(\mu T/\mathrm{eV\,yr})$", fontsize=8)
    cb2 = fig.colorbar(pcm, ax=axs[:, 1], location="bottom", shrink=0.9, pad=0.02)
    cb2.set_label(r"$p_\mathrm{max}$", fontsize=8)
    sup = str.maketrans("0123456789-.", "⁰¹²³⁴⁵⁶⁷⁸⁹⁻·")
    fig.suptitle(f"値の表: f = 10{f'{LF[kf]:g}'.translate(sup)} GeV、初期スピンの位置 s と α の面", fontsize=10, x=0.01,
                 ha="left")
    fig.text(0.01, -0.03, "s = 0: 322 の臨界、s = 1: 211 の臨界(破線)、s = 2: χ = 0.995。", fontsize=6.5, color="#555555")
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)

    # ---------- 7. 診断 ----------
    st = rt["status"]
    KI = z["rt_KI"].astype(float) * 0.01                       # [dex](255 = 2.55 dex 以上)
    RG = z["rt_RG"].astype(float)
    fig, axs = plt.subplots(2, 3, figsize=(7.0, 5.6), layout="constrained", sharex=True, sharey=True)
    for r, (piece, name) in enumerate([(1, "211 から育つ"), (0, "322 から育つ")]):
        bad = (st[piece] != 0).mean(axis=(1, 3))               # (alpha, f): s と M についての割合
        ki = KI[piece].max(axis=(1, 3, 4))
        dip = np.where(np.isfinite(RG[piece]), -RG[piece], -np.inf).max(axis=(1, 3, 4))
        dip = np.where(np.isfinite(dip), dip, np.nan)          # どの節も値がない所は斜線
        for c, (V, lab, vmax, title) in enumerate([(bad, "割合", 1.0, "解けなかった節の割合"),
                                                   (ki, "dex", 2.55, "KI の最大"),
                                                   (dip, "b の単位", 0.2, "戻りの最大")]):
            ax = axs[r, c]
            hatch_bg(ax)
            pc = ax.pcolormesh(AN, LF, V.T, cmap=seq_cmap(), vmin=0, vmax=vmax, shading="nearest", rasterized=True)
            ax.set_xscale("log")
            ax.set_title(f"({'abcdef'[3 * r + c]}) {name}\n{title}", loc="left")
            alpha_f_axes(ax, r, c, 3)
            if r == 1:
                cb = fig.colorbar(pc, ax=axs[:, c], location="bottom", shrink=0.9, pad=0.02)
                cb.set_label(lab, fontsize=8)
    fig.suptitle("値の表の診断(s 25 点と M 2 点について割合・最大)", fontsize=10, x=0.01, ha="left")
    fig.text(0.01, -0.03, f"status: 0 = 解けた {np.sum(st == 0)}、1 = 止まった {np.sum(st == 1)}、2 = 飛ばした "
             f"{np.sum(st == 2)} 節。KI は 255(2.55 dex)で頭打ち。", fontsize=6.5, color="#555555")
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)

print(f"wrote {os.path.relpath(OUT, ROOT)} ({os.path.getsize(OUT) / 1e6:.2f} MB)")
