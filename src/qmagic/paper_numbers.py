"""Every number printed in the manuscript, generated from the result files: LaTeX macros (numbers.tex) and table bodies (tab_*.tex).
    python -m qmagic.paper_numbers -> results/paper/   (run after claims_v3, sensitivity_v3, q3_lambda, noise_costmodel, gap, pec_exact)
"""
from qmagic.paths import RESULTS, find
import json, csv, math
OUT = RESULTS / "paper"; OUT.mkdir(exist_ok=True)
J = lambda n: json.load(open(find(n)))
prm, cl, se, tab = J("params_v3.json"), J("claims_v3.json"), J("sensitivity_v3.json"), J("pL_table.json")
M = {}
def sci(x, nd=1):
    if x is None or (isinstance(x, float) and math.isnan(x)): return "--"
    e = math.floor(math.log10(abs(x))); m = round(x / 10 ** e, nd)
    if m >= 10: m /= 10; e += 1
    return f"{m:.{nd}f}\\times10^{{{e}}}"
def num(x, nd=0): return f"{x:,.{nd}f}".replace(",", "{,}")
def rng(lo, hi, nd=0, scale=100): 
    a, b = f"{lo*scale:.{nd}f}", f"{hi*scale:.{nd}f}"; return a if a == b else f"{a}$--${b}"
def fx(x, nd=1): return f"{x:.{nd}f}"
US = (("0.0333", "Ser", "1/30"), ("0.1", "Mid", "0.1"), ("0.3", "Par", "0.3"), ("1.0", "Full", "1"))
# ---- savings (stage correction; d_max = 11)
for u, tag, _ in US:
    s = cl["stage"][u]["saving_range"]; M[f"sav{tag}"] = rng(*s); M[f"savScale{tag}"] = rng(*cl["scale"][u]["saving_range"])
    M[f"savRoute{tag}"] = rng(*cl["stage_route"][u]["saving_range"]); M[f"feas{tag}"] = str(cl["stage"][u]["n_feasible"])
    M[f"dQ{tag}"] = rng(*cl["stage"][u]["dQ_range"]); M[f"savNoFive{tag}"] = rng(*cl["scale_noD5"][u]["saving_range"])
    M[f"savStat{tag}"] = rng(min(cl[v][u]["saving_range"][0] for v in ("stage_statHi", "stage_statLo")), max(cl[v][u]["saving_range"][1] for v in ("stage_statHi", "stage_statLo")))
allred = [x for u, _, _ in US for x in (cl["stage"][u]["redA_range"] or [])]; M["redA"] = rng(min(allred), max(allred))
# ---- circuit size ratios at u = 1
SENS = ("stage", "stage_route", "stage_growth1", "stage_growthd", "stage_route_growthd")
ALLV = SENS + ("scale", "scale_route_growthd", "scale_noD5")
KEYS = (("0.001_11", "AEleven"), ("0.001_13", "AThirteen"), ("0.0005_11", "BEleven"), ("0.0005_13", "BThirteen"))
for key, tag in KEYS:
    n = cl["stage"]["1.0"]["nmax"][key]; M[f"rat{tag}"] = fx(n["ratio"], 1) if n["ratio"] >= 1 else fx(n["ratio"], 2)
    M[f"ratScale{tag}"] = fx(cl["scale"]["1.0"]["nmax"][key]["ratio"], 2)
    lo = min(cl[v]["1.0"]["nmax"][key]["ratio"] for v in ALLV); hi = max(cl[v]["1.0"]["nmax"][key]["ratio"] for v in ALLV)
    M[f"ratLo{tag}"] = fx(lo, 2 if lo < 1 else 1); M[f"ratHi{tag}"] = fx(hi, 1); M[f"nmaxA{tag}"] = sci(n["A"]); M[f"nmaxC{tag}"] = sci(n["C"])
    los = min(cl[v]["1.0"]["nmax"][key]["ratio"] for v in SENS); M[f"ratStageLo{tag}"] = fx(los, 2 if los < 1 else 1)
    a, b = cl["stage_statHi"]["1.0"]["nmax"][key]["ratio"], cl["stage_statLo"]["1.0"]["nmax"][key]["ratio"]   # sampled errors at the ends of their 95% Poisson intervals
    M[f"ratStatLo{tag}"] = fx(a, 2 if a < 1 else 1); M[f"ratStatHi{tag}"] = fx(b, 2 if b < 1 else 1); M[f"ratNoFive{tag}"] = fx(cl["scale_noD5"]["1.0"]["nmax"][key]["ratio"], 2)
for u, tag, _ in US[:2]:
    r = [cl["stage"][u]["nmax"][k]["ratio"] for k in ("0.001_11", "0.0005_11")]; M[f"ratEleven{tag}"] = rng(min(r), max(r), 2, 1)
M["nmaxBest"] = sci(max(cl[v][u]["nmax"][k][a] for v in ("stage",) for u, _, _ in US for k, _ in KEYS for a in ("A", "C")))
# ---- T-state calibration (u = 1)
for lab, tag in (("stage_rT0.3", "RTa"), ("stage_rT0.3_Ne1e-2", "RTb"), ("stage_rT0.03", "RTc")):
    for key, kt in (("0.001_11", "A"), ("0.0005_11", "B")):
        n = cl[lab]["1.0"]["nmax"][key]; M[f"nA{tag}{kt}"] = sci(n["A"]); M[f"nC{tag}{kt}"] = sci(n["C"]); M[f"rat{tag}{kt}"] = fx(n["ratio"], 1)
# ---- floors
for p, pt in ((1e-3, "A"), (5e-4, "B")):
    fl = lambda dmax: min(s["eps_T"] for s in prm["distill"] if s["p"] == p and s["d_req"] <= dmax)
    M[f"distFloor{pt}"] = sci(fl(11)); M[f"distFloorThirteen{pt}"] = sci(fl(13)); M[f"distSat{pt}"] = sci(fl(25))
    for d1, dt in ((3, "Three"), (5, "Five")):
        c = prm["proxy_correction"][f"d{d1}_{p:g}"]; M[f"stage{dt}{pt}"] = sci(c["T_stage"]); M[f"scaleFac{dt}{pt}"] = fx(c["ratio"], 1)
        M[f"proxyStage{dt}{pt}"] = sci(c["proxy_stage"]) if c["proxy_stage"] else "--"
        best = {}
        for s in prm["cultivate"]:
            if s["p"] != p or s.get("d1") != d1 or s.get("errors", 0) < 4 or "_d2=" not in s["name"] or (d1 == 5 and "deep" not in s["name"]): continue
            if s["d2"] not in best or s["eps_T"] < best[s["d2"]]["eps_T"]: best[s["d2"]] = s
        v = [best[d]["eps_T"] for d in sorted(best) if d <= 11]
        if v:
            M[f"cultLo{dt}{pt}"] = sci(min(v)); M[f"cultHi{dt}{pt}"] = sci(max(v)); M[f"cultOverStage{dt}{pt}"] = (rng(min(v) / c["T_stage"], max(v) / c["T_stage"], 1, 1) if max(v) / c["T_stage"] < 10 else rng(min(v) / c["T_stage"], max(v) / c["T_stage"], 0, 1))
            M[f"distOverCult{dt}{pt}"] = fx(fl(11) / min(v), 1); M[f"distOverStage{dt}{pt}"] = (fx(fl(11) / c["T_stage"], 1) if fl(11) / c["T_stage"] < 100 else num(round(fl(11) / c["T_stage"], -1)))
# ---- p_L table
T = tab["table"]; gf = tab["global_fit"]
from qmagic import costmodel
for p, pt in (("0.001", "A"), ("0.0005", "B")):
    M[f"pLEleven{pt}"] = sci(costmodel.pL(prm, float(p), 11)); M[f"pLNine{pt}"] = sci(costmodel.pL(prm, float(p), 9)); M[f"lam{pt}"] = fx(T[p]["Lambda"], 1)
    M[f"maxMeasured{pt}"] = max(T[p]["measured"], key=int)
M.update(fitA=fx(gf["A"], 3), fitPth=sci(gf["p_th"], 2), fitS=fx(gf["s"], 3), fitSdA=fx(gf["sd"]["A"], 3), fitSdS=fx(gf["sd"]["s"], 3), fitSdPth=sci(gf["sd"]["p_th"], 0),
         fitN=str(gf["n_points"]), fitMaxRatio=fx(gf["max_ratio"], 2))
rows = []
for p in ("0.001", "0.0005"):
    for d in ("3", "5", "7", "9", "11", "13"):
        if d in T[p]["measured"]: rows.append(f"${float(p):g}$ & {d} & ${sci(T[p]['measured'][d], 2)}$ & {num(T[p]['errors'][d])} & ${sci(T[p]['shots'][d], 1)}$ & meas. \\\\")
        elif int(d) <= 13: rows.append(f"${float(p):g}$ & {d} & ${sci(costmodel.pL(prm, float(p), int(d)), 2)}$ & {T[p]['errors'].get(d, '--')} & " + (f"${sci(T[p]['shots'][d], 1)}$" if d in T[p]["shots"] else "--") + " & extrap. \\\\")
open(OUT / "tab_pl.tex", "w").write("\n".join(rows) + "\n")
# ---- u table
rows = [f"${ul}$ & {round(float(u) * 100 * 0.3):d} & ${M['sav'+tag]}\\%$ & ${M['savScale'+tag]}\\%$ & ${M['savRoute'+tag]}\\%$ & {M['feas'+tag]}/6 \\\\" for u, tag, ul in US]
open(OUT / "tab_u.tex", "w").write("\n".join(rows) + "\n")
# ---- worked example
rows = []
for u, ul in (("0.0333", "1/30"), ("0.1", "0.1"), ("1", "1")):
    ex = se[f"example_u={u}"]
    if not ex.get("A") and not ex.get("C"): rows.append(f"${ul}$ & A, C & -- & \\multicolumn{{4}}{{c}}{{infeasible (idle error)}} \\\\"); continue
    for a in ("A", "C"):
        r = ex.get(a)
        rows.append(f"${ul}$ & {a} & " + (f"{r['d']} & {num(r['Q_data'])} & {num(r['Q_fac'])} & {fx(r['Gamma'], 1)} & ${sci(r['T_cycles'])}$ \\\\" if r else "-- & \\multicolumn{4}{c}{infeasible} \\\\"))
        if r:
            t = {"0.1": "Mid", "1": "Full", "0.0333": "Ser"}[u]
            M[f"ex{a}{t}Q"] = num(round(r["Q"], -2)); M[f"ex{a}{t}D"] = str(r["d"]); M[f"ex{a}{t}Fac"] = num(round(r["Q_fac"], -1 if r["Q_fac"] < 1000 else -2)); M[f"ex{a}{t}K"] = str(r["k"]); M[f"ex{a}{t}G"] = fx(r["Gamma"], 0)
            M[f"ex{a}{t}Src"] = r["src"].replace("_", "\\_")
    if ex.get("A") and ex.get("C"): M[f"exSav{ {'0.1': 'Mid', '1': 'Full', '0.0333': 'Ser'}[u] }"] = fx(100 * (1 - ex["C"]["Q"] / ex["A"]["Q"]), 0)
open(OUT / "tab_example.tex", "w").write("\n".join(rows) + "\n")
q = se["qtg"]; M["qtgN"] = str(sum(v["n"] for v in q.values())); M["qtgBetter"] = str(sum(v["n_C_better"] for v in q.values()))
M["qtgLo"] = fx(min(v["min"] for v in q.values()), 2); M["qtgHi"] = fx(max(v["max"] for v in q.values()), 2); M["qtgFull"] = rng(q["u=1"]["min"], q["u=1"]["max"], 2, 1)
# ---- N_max table (u = 1)
rows = []
for key, tag in KEYS:
    p, dm = key.split("_"); ext = "" if str(dm) in T[p]["measured"] else "$^{\\dagger}$"
    rows.append(f"${float(p):g}$ & {dm}{ext} & ${M['nmaxA'+tag]}$ & ${M['nmaxC'+tag]}$ & {M['rat'+tag]} & {M['ratStatLo'+tag]}--{M['ratStatHi'+tag]} & {M['ratStageLo'+tag]} & {M['ratScale'+tag]} \\\\")
open(OUT / "tab_nmax.tex", "w").write("\n".join(rows) + "\n")
# ---- T-state calibration table
rows = []
for lab, name in (("stage", "$r_T=0$"), ("stage_rT0.03", "$r_T=0.03$"), ("stage_rT0.3", "$r_T=0.3$"), ("stage_rT0.3_Ne1e-2", "$r_T=0.3$, $\\delta=10^{-2}$")):
    a, b = cl[lab]["1.0"]["nmax"]["0.001_11"], cl[lab]["1.0"]["nmax"]["0.0005_11"]
    rows.append(f"{name} & ${sci(a['A'])}$ & ${sci(a['C'])}$ & {fx(a['ratio'], 1)} & ${sci(b['A'])}$ & ${sci(b['C'])}$ & {fx(b['ratio'], 1)} \\\\")
open(OUT / "tab_rt.tex", "w").write("\n".join(rows) + "\n")
# ---- calibration cost of the data block (1 microsecond per cycle)
def dur(t):
    for lim, u, n in ((60, 1, "s"), (3600, 60, "min"), (86400, 3600, "h"), (86400 * 365, 86400, "d"), (1e30, 86400 * 365, "y")):
        if t < lim: v = t / u; return (f"{v:.0f}" if v >= 10 else f"{v:.1f}") + f"\\,{n}"
rows = []
for p in (1e-3, 5e-4):
    for d in (9, 11):
        e = costmodel.pL(prm, p, d); t2, t3 = d * 1e-6 / (e * 1e-4), d * 1e-6 / (e * 1e-6) / prm["n_L"]
        assert abs(t2 / t3 - 1) < 1e-9          # per patch to 1e-2 and block to 1e-3 over n_L=100 parallel patches need the same time
        rows.append(f"${p:g}$ & {d} & ${sci(e)}$ & ${sci(1e4 / e)}$ & {dur(t2)} & {dur(t2 * 1e3)} \\\\")
open(OUT / "tab_calib.tex", "w").write("\n".join(rows) + "\n")
# ---- sensitivity table
lab = [("baseline", "nL=100"), ("$n_L=20$", "nL=20"), ("$n_L=300$", "nL=300"), ("$f_T=0.1$", "fT=0.1"), ("$f_T=1$", "fT=1.0"), ("$r_{\\rm route}=1.5$", "r=1.5"), ("$r_{\\rm route}=3$", "r=3.0"),
       ("$\\Gamma_{\\max}=10$", "G=10"), ("$\\Gamma_{\\max}=10^{4}$", "G=10000"), ("$r=10^{-2}$ ($\\Sigma\\le0.1$)", "r=1e-2"), ("routing tiles charged, $\\kappa(u)$", "kappa=route"),
       ("growth, $d$ rounds at $d_2$", "growth=d"), ("pessimistic proxy correction", "proxy=scale"), ("sampled proxy values, no correction", "proxy=none"), ("$r_T=0.3$", "rT=0.3")]
rows = []
for name, k in lab:
    cells = []
    for u in ("0.0333", "0.1", "1"):
        v = se.get(f"u={u},{k}", [None])[0]; cells.append(f"${rng(v[0], v[1])}\\%$ ({v[2]})" if v else "none feasible")
    rows.append(f"{name} & " + " & ".join(cells) + " \\\\")
rows.append("\\hline\nspacetime ratio C/A & " + " & ".join(f"${rng(q[f'u={u}']['min'], q[f'u={u}']['max'], 2, 1)}$" for u in ("0.0333", "0.1", "1")) + " \\\\")
open(OUT / "tab_sens.tex", "w").write("\n".join(rows) + "\n")
# ---- hardware mapping
h = {}
for r in csv.DictReader(open(find("q3_nmax.csv"))):
    if r["d_max"] == "11": h[r["source"]] = r
for src, tag in (("Lacroix2025_SC_colour_code", "Lac"), ("Bluvstein2026_neutral_atoms", "Blu"), ("GoogleWillow2025_surface_code", "Wil"), ("Lavasani2026_ion_colour_code", "Lav"), ("reference p=2e-3", "Ref")):
    r = h[src]; M[f"hw{tag}Lam"] = fx(float(r["Lambda"]), 2 if float(r["Lambda"]) < 3 else 1); M[f"hw{tag}P"] = sci(float(r["p"])); n = float(r["N_data"])
    M[f"hw{tag}N"] = (f"{n:.0f}" if n < 1000 else sci(n)); M[f"hw{tag}NA"] = sci(float(r["Nmax_A"])) if float(r["Nmax_A"]) > 10 else fx(float(r["Nmax_A"]), 0)
# ---- noise models, gap, PEC, two-level factory
nz = J("noise_costmodel.json")
rows = []
for m, name in (("uniform", "uniform"), ("biased_10", "bias $\\eta=10$"), ("biased_100", "bias $\\eta=100$"), ("erasure_0.5", "erasure 50\\%"), ("erasure_0.9", "erasure 90\\%")):
    r = nz[m]; f = r["fit"]; res = m.startswith(("uniform", "biased"))
    rows.append(f"{name} & {fx(f['Lambda_1e3'], 1 if f['Lambda_1e3'] < 10 else 0)} & ${sci(f['A'] * (1e-3 / f['p_th']) ** (f['s'] * 6))}$ & ${sci(r['dist_floor'])}$ & " + (f"${sci(r['cult_floor'])}$ & {fx(r['NmaxC'] / r['NmaxA'], 2)}" if res else "-- & --") + " \\\\")
    M["nz" + m.replace("_", "").replace(".", "").replace("0", "o").replace("1", "i").replace("5", "v").replace("9", "n") + "Lam"] = fx(f["Lambda_1e3"], 1 if f["Lambda_1e3"] < 10 else 0)
open(OUT / "tab_noise.tex", "w").write("\n".join(rows) + "\n")
M["nzUniformOverMeas"] = fx(nz["uniform"]["fit"]["A"] * (1e-3 / nz["uniform"]["fit"]["p_th"]) ** (nz["uniform"]["fit"]["s"] * 6) / T["0.001"]["measured"]["11"], 1)
M["nzFallLo"] = fx(min(nz["uniform"]["NmaxA"] / nz[m]["NmaxA"] for m in ("biased_10", "biased_100")), 0); M["nzFallHi"] = fx(max(nz["uniform"]["NmaxA"] / nz[m]["NmaxA"] for m in ("biased_10", "biased_100")), 0)
M["nzRatUniform"] = fx(nz["uniform"]["NmaxC"] / nz["uniform"]["NmaxA"], 2); M["nzRatBias"] = rng(min(nz[m]["NmaxC"] / nz[m]["NmaxA"] for m in ("biased_10", "biased_100")), max(nz[m]["NmaxC"] / nz[m]["NmaxA"] for m in ("biased_10", "biased_100")), 1, 1)
g = J("gap_qem.json"); sl = [v["fit"][0] for k, v in g["maps"].items() if k.startswith("d1=3")]; M["gapSlope"] = fx(sum(sl) / len(sl), 3); M["gapSlopeLo"] = fx(max(sl), 3); M["gapSlopeHi"] = fx(min(sl), 3)
rr = [v for k, v in g["cross"].items() if k.endswith("err>=4") and k.startswith("d1=3") and v]; M["gapPredLo"] = fx(min(v[0] for v in rr), 2); M["gapPredHi"] = fx(max(v[1] for v in rr), 1)
pe = J("pec_exact.json"); b = [x for x in pe["b"] if x["r_mis"] != 0]; M["pecBiasLo"] = fx(min(x["bias"] / x["rule_sigma_r"] for x in b), 1); M["pecBiasHi"] = fx(max(x["bias"] / x["rule_sigma_r"] for x in b), 1)
M["pecGammaDev"] = fx(100 * max(abs(x["gamma2"] / x["exp4sum"] - 1) for x in pe["a"]), 1)
tl = J("factories_2level_d11.json"); one = min(s["eps_T"] for s in prm["distill"] if s["p"] == 1e-3 and s["d_req"] <= 11)
two = min((r for r in tl if r["p"] == 1e-3), key=lambda r: r["p_out"]); M["twoLevelGain"] = fx(100 * (1 - two["p_out"] / one), 0); M["twoLevelQubits"] = fx(two["qubits"] / 2066, 0)
# ---- extra traced numbers used in the prose
u11 = T["0.0005"]["unmeasured"].get("11")
if u11: M["pLElevenBObs"] = sci(u11["pL"]); M["pLElevenBErr"] = str(u11["errors"])
M["pLThirteenA"] = sci(T["0.001"]["measured"]["13"]) if "13" in T["0.001"]["measured"] else "--"
for u, t in (("0.0333", "Ser"), ("0.1", "Mid"), ("1", "Full")):
    v = se.get(f"u={u},rT=0.3", [None])[0]
    if v: M[f"savRT{t}"] = rng(v[0], v[1]); M[f"feasRT{t}"] = str(v[2])
    ex = se[f"example_u={u}"]
    if ex.get("A") and ex.get("C"): M[f"exRun{t}"] = fx(ex["C"]["T_cycles"] / ex["A"]["T_cycles"], 1 if ex["C"]["T_cycles"] / ex["A"]["T_cycles"] < 10 else 0)
bv = {}
for f in ("cultiv_d2scan_d3.json", "cultiv_validate_d15.json"):
    for r in J(f):
        if r["meta"]["p"] == 1e-3 and r["meta"]["d1"] == 3: bv[r["meta"]["d2"]] = r["baseline_volume"]
for d2, t in ((7, "Seven"), (9, "Nine"), (11, "Eleven")): M[f"volFrac{t}"] = fx(100 * bv[d2] / bv[15], 0)
cnt = {d2: min((s for s in prm["cultivate"] if s["p"] == 1e-3 and s.get("d1") == 3 and s.get("d2") == d2 and s.get("errors", 0) >= 4 and "_d2=" in s["name"]), key=lambda s: s["eps_T"])["errors"] for d2 in (7, 9, 11)}
M["cultErrLo"] = str(min(cnt.values())); M["cultErrHi"] = str(max(cnt.values()))
M["fitOverThirteenA"] = fx(next(1 / r["ratio"] for r in tab["residuals"] if r["p"] == 1e-3 and r["d"] == 13), 2)
# ---- write macros
import re
bad = [k for k in M if not re.fullmatch(r"[A-Za-z]+", k)]; assert not bad, bad
open(OUT / "numbers.tex", "w").write("% generated by qmagic.paper_numbers; do not edit\n" + "".join(f"\\newcommand{{\\{k}}}{{{v}}}\n" for k, v in sorted(M.items())))
json.dump(M, open(OUT / "numbers.json", "w"), indent=1); print(len(M), "macros;", "tables:", sorted(p.name for p in OUT.glob("tab_*.tex")))
