"""chi_SR の単体版の自己診断: この計算機で、作った側と同じ値・同じ印を返すかを確かめる。

    python selftest.py

selftest.npz に、いろいろな領域の入力 5000 個と、作った側(Apple Silicon、numpy 2.3)での出力が入っている。警告を例外にして
呼び、値の差の最大と、印の違う個数を表示する。CPU や numpy の版が違うと、丸めの違いで閾値のきわの印が変わることはあるが、
値の差は 1e-9 以下のはずである。それより大きい差や例外(IndexError など)が出たら、表(chi_sr_tables.npz)と chi_sr.py の
組が合っているかを確かめ、作った側に知らせてほしい。
"""
import os
import pickle
import sys
import time
import warnings

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from chi_sr import ChiSR  # noqa: E402


def main():
    z = np.load(os.path.join(HERE, "selftest.npz"))
    t0 = time.time()
    model = ChiSR(os.path.join(HERE, "chi_sr_tables.npz"))
    t_load = time.time() - t0
    with warnings.catch_warnings():
        warnings.simplefilter("error")                   # 警告が出たら例外で止まる
        t0 = time.time()
        chi, reg = model(z["M"], z["chi0"], z["T"], z["mu"], z["f"])
        t_call = time.time() - t0
        chi1, reg1 = model(z["M"][:50], z["chi0"][:50], z["T"][:50], z["mu"][:50], z["f"][:50])   # 分け方によらない
        m2 = pickle.loads(pickle.dumps(model))           # pickle(中身は表の場所だけ)
        chi2, reg2 = m2(z["M"][:50], z["chi0"][:50], z["T"][:50], z["mu"][:50], z["f"][:50])
    dmax = float(np.max(np.abs(chi - z["chi"])))
    nreg = int(np.sum(reg != z["regime"]))
    same = bool(np.array_equal(chi1, chi[:50]) and np.array_equal(reg1, reg[:50]) and np.array_equal(chi2, chi1)
                and np.array_equal(reg2, reg1))
    print(f"inputs: {len(chi)}   load {t_load:.1f} s   call {t_call:.2f} s")
    print(f"max |chi - expected| = {dmax:.3e}")
    print(f"regime differs in {nreg} of {len(chi)} inputs")
    print(f"regime share (0..5): {[round(float(np.mean(reg == k)), 3) for k in range(6)]}")
    print(f"same result in chunks and after pickling: {same}")
    ok = dmax <= 1e-9 and nreg <= 5 and same
    print("OK" if ok else "MISMATCH: この計算機では作った側と違う結果になる")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
