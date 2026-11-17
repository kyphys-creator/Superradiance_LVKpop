"""Mac と Hoffman2 の照合用の小さな試し(1 コア、15 分ほど)。 python h2/compare.py <out.json>
kinn.py check の 3 項目と、kinn.corr_grid のうち軽い 4 例(kinn._corr_job)を解いて JSON に書く。Mac で同じものを解いて比べる。"""
import json
import os
import platform
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np  # noqa: E402
import scipy  # noqa: E402
import kinn  # noqa: E402
import sr  # noqa: E402

PICK = [(0.03, 1e13), (0.12, 1e15), (0.30, 1e13), (0.45, 1e17)]


def main(out):
    t0 = time.time()
    rep = dict(host=platform.node(), machine=platform.machine(), python=platform.python_version(),
               numpy=np.__version__, scipy=scipy.__version__)
    rep["check"] = {"evolve3m": kinn.check_evolve3m(), "fixed_mass": kinn.check_fixed_mass(),
                    "conservation": kinn.check_conservation()}
    mu1 = sr.alpha_of(30.0, 1e-12) / 1e-12
    cases = [c for c in kinn.corr_grid() if any(abs(c[2] * mu1 - a) < 1e-9 and c[3] == f for a, f in PICK)]
    rep["cases"] = [kinn._corr_job(c) for c in cases]
    rep["secs"] = time.time() - t0
    json.dump(rep, open(out, "w"), default=float)
    print("done", len(rep["cases"]), "cases", round(rep["secs"]), "s")


if __name__ == "__main__":
    main(sys.argv[1])
