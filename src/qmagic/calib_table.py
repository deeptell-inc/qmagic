"""Measured data-block logical error table p_L(p, d) = p_X + p_Z per d rounds, from the shipped memory calibrations
(data/calib_gidney.json plus the deep points of data/calib_deep.json; same circuits, noise and decoder, independent seeds).
A point enters the table when it has at least MIN_ERR logical errors (X+Z). Beyond the largest measured distance of a given p
the model extrapolates with a two-parameter fit A x^{(d+1)/2} to the three largest measured distances of that p.
The global three-parameter fit A (p/p_th)^{s (d+1)/2} is kept for error rates outside the table (hardware mapping) and is
reported with residuals and bootstrap uncertainties.
    python -m qmagic.calib_table -> results/pL_table.json
"""
from qmagic.paths import RESULTS, find
import json, math, numpy as np
MIN_ERR = 30
def load():
    agg = {}
    for name in ("calib_gidney.json", "calib_deep.json"):
        try: rows = json.load(open(find(name)))["rows"]
        except FileNotFoundError: continue
        for r in rows:
            a = agg.setdefault((r["p"], r["d"]), {}).setdefault(r["basis"], [0, 0]); a[0] += r["shots"]; a[1] += r["errors"]
    pts = {}
    for (p, d), v in agg.items():
        if "X" not in v or "Z" not in v: continue
        (sx, ex), (sz, ez) = v["X"], v["Z"]; pl = ex / sx + ez / sz; n = ex + ez
        pts[(p, d)] = dict(p=p, d=d, shots_X=sx, shots_Z=sz, errors=n, pL=pl, rel_err=(math.sqrt(ex / sx ** 2 + ez / sz ** 2) / pl) if pl > 0 else None)
    return pts
def global_fit(pts, errs=None):
    used = [(p, d, v["pL"] if errs is None else errs[(p, d)]) for (p, d), v in pts.items() if v["errors"] >= MIN_ERR]
    used = [(p, d, pl) for p, d, pl in used if pl > 0]
    X = np.array([[1.0, (d + 1) / 2, (d + 1) / 2 * math.log(p)] for p, d, _ in used]); y = np.array([math.log(pl) for *_, pl in used])
    logA, b, c = np.linalg.lstsq(X, y, rcond=None)[0]; return math.exp(logA), math.exp(-b / c), c
def main():
    pts = load(); out = {}; rng = np.random.default_rng(0)
    for p in sorted({p for p, _ in pts}):
        rows = sorted((v for (pp, d), v in pts.items() if pp == p), key=lambda v: v["d"])
        meas = [v for v in rows if v["errors"] >= MIN_ERR]; top = meas[-3:]
        X = np.array([[1.0, (v["d"] + 1) / 2] for v in top]); y = np.array([math.log(v["pL"]) for v in top]); w = np.array([v["errors"] for v in top], float)
        lnA, lnx = np.linalg.lstsq(X * np.sqrt(w)[:, None], y * np.sqrt(w), rcond=None)[0]
        out[f"{p:g}"] = dict(measured={str(v["d"]): v["pL"] for v in meas}, rel_err={str(v["d"]): v["rel_err"] for v in meas},
                             errors={str(v["d"]): v["errors"] for v in rows}, shots={str(v["d"]): v["shots_X"] + v["shots_Z"] for v in rows},
                             unmeasured={str(v["d"]): dict(errors=v["errors"], pL=v["pL"]) for v in rows if v["errors"] < MIN_ERR},
                             A=math.exp(lnA), x=math.exp(lnx), Lambda=math.exp(-lnx), fit_distances=[v["d"] for v in top])
    A, pth, s = global_fit(pts)
    boot = []
    for _ in range(300):                                   # Poisson bootstrap of the error counts
        e = {k: (rng.poisson(v["errors"]) / max(v["errors"], 1)) * v["pL"] for k, v in pts.items()}; boot.append(global_fit(pts, e))
    boot = np.array(boot); cov = np.cov(np.log(boot[:, :2]).T.tolist() + [boot[:, 2].tolist()])
    resid = [dict(p=p, d=d, errors=v["errors"], observed=v["pL"], global_fit=A * (p / pth) ** (s * (d + 1) / 2),
                  ratio=(v["pL"] / (A * (p / pth) ** (s * (d + 1) / 2))) if v["pL"] > 0 else None, in_table=v["errors"] >= MIN_ERR)
             for (p, d), v in sorted(pts.items())]
    res = dict(min_err=MIN_ERR, table=out, global_fit=dict(A=A, p_th=pth, s=s, sd=dict(A=float(boot[:, 0].std()), p_th=float(boot[:, 1].std()), s=float(boot[:, 2].std())),
                                                          cov_lnA_lnpth_s=cov.tolist(), n_points=sum(1 for v in pts.values() if v["errors"] >= MIN_ERR),
                                                          max_ratio=max(max(r["ratio"], 1 / r["ratio"]) for r in resid if r["in_table"])), residuals=resid)
    json.dump(res, open(RESULTS / "pL_table.json", "w"), indent=1)
    print(f"global fit: A={A:.3f}±{boot[:,0].std():.3f} p_th={pth:.5f}±{boot[:,1].std():.5f} s={s:.3f}±{boot[:,2].std():.3f} over {res['global_fit']['n_points']} points, max ratio {res['global_fit']['max_ratio']:.2f}")
    for p, r in out.items():
        print(f"p={p}: measured " + ", ".join(f"d={d}: {v:.3e} (n={r['errors'][d]}, ±{100*r['rel_err'][d]:.0f}%)" for d, v in r["measured"].items())
              + f" | extrapolation from d={r['fit_distances']}: Lambda={r['Lambda']:.2f}" + "".join(f" | d={d} below threshold: {u['errors']} errors" for d, u in r["unmeasured"].items()))
if __name__ == "__main__": main()
