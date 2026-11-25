"""Re-derive the surviving claims (b) and (c) under the unified clock, separating measured from scaled cultivation rows.
    python -m qmagic.claims_v3  -> results/claims_v3.json + stdout tables
Claim (b): serial-ish architectures (u = pi/n_L small): cultivation adds little to Q, no N_max gain.
Claim (c): parallel (u ~ 1): hybrid saves most of the factory area vs QEM alone; N_max gain modest and conditional.
"""
from qmagic.paths import DATA, RESULTS, find
import json, math, copy
from qmagic import costmodel
prm0 = json.load(open(find("params_v3.json"))); G = 100.0
nmax = lambda prm, kind, p, d_max, f_T: costmodel.nmax(prm, kind, p, d_max, f_T, G)
# label, proxy mode, model options.  "stage" is the main result; "scale" the pessimistic proxy correction; "proxy" the sampled S-proxy values.
VARIANTS = (("stage", "stage", {}), ("scale", "scale", {}), ("proxy", "none", {}),
            ("stage_route", "stage", {"kappa_mode": "route"}),            # routing tiles charged in proportion to duty cycle u
            ("stage_growth1", "stage", {"growth_rounds": 1}),            # 1 idle round at d2 before growth (d2 < d rows)
            ("stage_growthd", "stage", {"growth_rounds": "d"}),          # whole consumption at d2: upper variant
            ("stage_route_growthd", "stage", {"kappa_mode": "route", "growth_rounds": "d"}),
            ("scale_route_growthd", "scale", {"kappa_mode": "route", "growth_rounds": "d"}),
            ("stage_rT0.3", "stage", {"bias_rT": 0.3}),                   # T-state error known to 30%: Sigma_T * 0.3 <= N_e on mitigated arms
            ("stage_rT0.3_Ne1e-2", "stage", {"bias_rT": 0.3, "N_e": 1e-2, "bias_r": 1e-2}),
            ("stage_rT0.03", "stage", {"bias_rT": 0.03}))
from scipy.stats import chi2
def stat(prm, side):
    """Sampled cultivation errors moved to the end of their 95% Poisson interval (per row, from its error count); stage correction kept."""
    q = copy.deepcopy(prm)
    for s in q["cultivate"]:
        n = s.get("errors")
        if not n: continue
        f = chi2.ppf(0.975, 2 * n + 2) / 2 / n if side == "hi" else chi2.ppf(0.025, 2 * n) / 2 / n
        s["eps_T"] = s["eps_T_proxy"] * f + (s["eps_T_stage"] - s["eps_T_proxy"])
    return q
def no_d5_low_p(prm):
    """Pessimistic correction when the d1=5 multiplier at p=5e-4 is unbounded (no proxy error was observed there): d1=5 rows removed at that p."""
    q = copy.deepcopy(prm); q["cultivate"] = [s for s in q["cultivate"] if not (s["p"] == 5e-4 and s.get("d1") == 5)]; return q
TRANSFORM = {"stage_statHi": lambda q: stat(q, "hi"), "stage_statLo": lambda q: stat(q, "lo"), "scale_noD5": no_d5_low_p}
VARIANTS += (("stage_statHi", "stage", {}), ("stage_statLo", "stage", {}), ("scale_noD5", "scale", {}))
out = {}
for label, mode, opts in VARIANTS:
    prm = costmodel.with_proxy_mode(prm0, mode); prm.update(opts)
    if label in TRANSFORM: prm = TRANSFORM[label](prm)
    res = {}
    for u in (0.0333, 0.1, 0.3, 1.0):
        prm["parallelism"] = u * prm["n_L"]; prm["p"] = [1e-3, 5e-4]; prm["d_max"] = [11, 13]; prm["N"] = [1e4, 1e5, 1e6]; prm["f_T"] = [0.3]; prm["C_max"] = [G]
        rows = costmodel.sweep(prm)
        q = {}
        for r in rows:
            if r["arm"] in ("A", "C"): q[(r["arm"], r["p"], r["d_max"], r["N"])] = r
        sav, dq, redA = [], [], []
        for p in (1e-3, 5e-4):
            for N in (1e4, 1e5, 1e6):
                a, c = q[("A", p, 11, N)], q[("C", p, 11, N)]
                if a["feasible"] and c["feasible"]: sav.append(c["saving_C_vs_A"]); dq.append(c["deltaQ_C_minus_A"]); redA.append(a["reduction_vs_ref"])
        nm = {f"{p:g}_{dm}": (nmax(prm, "distill", p, dm, 0.3), nmax(prm, "cultivate", p, dm, 0.3)) for p in (1e-3, 5e-4) for dm in (11, 13)}
        res[u] = dict(saving_range=(min(sav), max(sav)) if sav else None, dQ_range=(min(dq), max(dq)) if dq else None, redA_range=(min(redA), max(redA)) if redA else None,
                      n_feasible=len(sav), nmax={k: dict(A=v[0], C=v[1], ratio=v[1] / v[0]) for k, v in nm.items()})
        print(f"[{label}] u={u:g} (rho={u*prm['n_L']*0.3:.0f}): saving={None if not sav else (round(min(sav),3), round(max(sav),3))} feasible={len(sav)}/6 | N_max C/A: " +
              ", ".join(f"{k}: {v['ratio']:.2f} (A {v['A']:.1e}, C {v['C']:.1e})" for k, v in res[u]["nmax"].items()))
    out[label] = res
json.dump(out, open(RESULTS / "claims_v3.json", "w"), indent=1, default=float)
