"""Panel-corrected parameter set (v3; v4 = full factory grid + bias constraint): Gidney-noise X+Z calibration, factory grid over ALL odd d_X (5..25) with exit-only p_out,
physical d_req for every cultivation row, symmetric consumption 0.75*p_L(d) on both arms, unified clock (rho = pi*f_T, kappa=1).
    python -m qmagic.build_params_v3  -> results/params_v3.json (+ prints factory floors)
"""
from qmagic.paths import DATA, RESULTS, EXTERNAL, find
import json, glob
from qmagic import factory as F
tab = json.load(open(find("pL_table.json"))); gf = tab["global_fit"]          # measured table + per-p extrapolation (calib_table)
base = json.load(open(find("params_corr_d2.json")))
prm = dict(base); prm["pL_table"] = tab["table"]
prm["pL_fit"] = dict(A=gf["A"], p_th=gf["p_th"], s=gf["s"], source="pL_table.json global fit (Gidney memory circuit, uniform_depolarizing incl. idle, X+Z, correlated matching); used only outside the measured table")
from qmagic import costmodel as _cm
pL_cycle = lambda p, d: _cm.pL(prm, p, d) / d                                 # per cycle (X+Z): measured where available
prm.update(unified_clock=True, kappa=1.0, consumption=True, parallelism=None,
           bias_r=1e-3, kappa_mode="fixed", growth_rounds=0)   # common precision target N_e for the data-block part (Sigma*bias_r <= N_e); kappa(u), growth, bias_rT are sensitivities
prm["_note"] = "v5: measured p_L table; exit-only factory p_out on the full (dX,dZ,dm) grid; symmetric consumption; physical d_req; unified clock; cultivation errors = sampled S-proxy value corrected to actual T gates with the exact stage rates of Wan & Zapirain (eps_T = stage correction, eps_T_scale = pessimistic, eps_T_proxy = as sampled)."
# distillation: exit-only p_out, all odd dX, Litinski's half-weights use pL(X+Z) -> pass pL_cycle/2 per type? Litinski's rule: X and Z each 0.5*pL where pL is the TOTAL per-cycle logical rate -> pass the total.
dist = []
for p in (1e-3, 5e-4):
    for dX in range(5, 27, 2):
        rs = [F.simulate(dX, dZ, dm, p, pL=pL_cycle, consume=False) for dZ in (3, 5, 7, 9) for dm in (3, 5, 7, 9) if dZ <= dX and dm <= dX]
        for r in rs:                                   # v4: keep the FULL (dZ, dm) grid; the sweep picks the Q-optimal factory under the error constraint
            dist.append(dict(name=f"15to1_{r['dX']},{r['dZ']},{r['dm']}", p=p, d_req=dX, qubits=r["qubits"], cycles_per_T=round(r["cycles"], 1), eps_T=max(r["p_out"], 1e-16),
                             source="factory.simulate(consume=False), Gidney-noise X+Z calibration"))
    print(f"p={p:g}: exit-only distillation floor per dX: " + ", ".join(f"{dX}: {min(s['eps_T'] for s in dist if s['p']==p and s['d_req']==dX):.1e}" for dX in range(5, 27, 2)))
prm["distill"] = dist
# cultivation: physical d_req = d2 for every row; rows that transfer published d2=15 tails to smaller d2 are dropped
cult = []
for s in base["cultivate"]:
    if s.get("scaled"): continue
    s = dict(s); s["d_req"] = int(s.get("d2") or 15); cult.append(s)
# deep-tail d1=5, d2=11 reruns (5e9 / 2e9 shots): measured rows, physical d_req=11
def thin(curve, minerr=4):
    rs = sorted([x for x in curve if x["errors_raw"] >= minerr], key=lambda x: (x["err"], x["v"])); out, prev = [], None
    for x in rs:
        if prev is not None and x["v"] > prev * 0.9 and x["err"] != 0: continue
        prev = x["v"]; out.append(x)
    return out
for f in ("cultiv_d5deep_d2_11_p1e3_5e9.json", "cultiv_d5deep_d2_11_p5e4_2e9.json"):
    try: runs = json.load(open(find(f)))
    except FileNotFoundError: continue
    for r in runs:
        m = r["meta"]
        for x in thin(r["curve"]):
            cult.append(dict(name=f"cult_d1=5_d2=11_deep_att{x['attempts']:.1f}", p=m["p"], d_req=11, eps_T=x["err"], volume_per_T=x["v"], estimate=False, scaled=False,
                             attempts=round(x["attempts"], 2), errors=x["errors_raw"], d1=5, d2=11, q=r["q"], source=f"deep-tail rerun {f} (5e9/2e9 shots), errors x2 (S->T)"))
# S-gate proxy -> actual T gates. The sampled error (x2 convention) of a row is e_proxy = 2*(stage proxy) + escape part.
# Wan & Zapirain give the exact T-gate error T_stage of the injection+cultivation stage, so
#   stage  : eps_T = e_proxy + (T_stage - 2*proxy_stage)        (only the non-Clifford stage differs; default)
#   scale  : eps_T = e_proxy * T_stage / (2*proxy_stage)         (the whole sampled error scales like the stage; pessimistic)
from fractions import Fraction
wz = json.load(open(find("wan_zapirain_2026.json")))
def T_stage(d1, p):
    x = Fraction(p).limit_denominator(10 ** 12); x = x / (1 - x)
    return float(sum(Fraction(c) * x ** int(k) for k, c in wz[f"d{d1}"]["series"].items()))
corr = {}
for d1 in (3, 5):
    for p in (1e-3, 5e-4):
        T = T_stage(d1, p); px = wz["proxy_stage"][f"d{d1}"].get(f"{p:g}")
        if px is None: px_add, px_ratio = 0.0, wz["proxy_stage"][f"d{d1}_upper95"][f"{p:g}"]      # no proxy error observed: no subtraction; ratio from the 95% limit
        else: px_add = px_ratio = px
        corr[(d1, p)] = dict(T_stage=T, proxy_stage=px, add=T - 2 * px_add, ratio=T / (2 * px_ratio))
        print(f"d1={d1} p={p:g}: exact T stage {T:.3e}, proxy stage {px}, stage correction +{corr[(d1,p)]['add']:.2e}, scale factor x{corr[(d1,p)]['ratio']:.2f}")
for s in cult:
    c = corr.get((s.get("d1"), s["p"]))
    s["eps_T_proxy"] = s["eps_T"]
    if c: s["eps_T_stage"] = s["eps_T"] + c["add"]; s["eps_T_scale"] = s["eps_T"] * c["ratio"]; s["eps_T"] = s["eps_T_stage"]
prm["proxy_correction"] = {f"d{d1}_{p:g}": v for (d1, p), v in corr.items()}
prm["cultivate"] = cult
json.dump(prm, open(RESULTS / "params_v3.json", "w"), indent=1)
print(len(dist), "distillation rows;", len(cult), "cultivation rows")
