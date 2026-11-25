"""Fig. 3 of the paper. Left: qubit saving of mitigated cultivation (C) over mitigated distillation (A) versus u (rho = n_L u f_T),
band over p and N. Right: ratio of the largest feasible circuit sizes at d_max = 11; solid: stage-only correction of the S-gate
proxy; dashed: pessimistic scaling of the whole sampled error; shaded: routing-tile error and growth charge."""
from qmagic.paths import find, RESULTS
import json, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
c = json.load(open(find("claims_v3.json")))
us = sorted(float(u) for u in c["stage"]); g = lambda lab, u: c[lab][str(u) if str(u) in c[lab] else f"{u:g}"]
fig, ax = plt.subplots(1, 2, figsize=(10.5, 4.2))
lo = [g("stage", u)["saving_range"][0] * 100 for u in us]; hi = [g("stage", u)["saving_range"][1] * 100 for u in us]
ax[0].fill_between(us, lo, hi, alpha=.25, label="range over $p$ and $N$"); ax[0].plot(us, hi, "o-", color="C0"); ax[0].plot(us, lo, "o-", color="C0")
ax[0].plot(us, [g("scale", u)["saving_range"][0] * 100 for u in us], "k--", lw=1, label="pessimistic proxy correction"); ax[0].plot(us, [g("scale", u)["saving_range"][1] * 100 for u in us], "k--", lw=1)
ax[0].set_xscale("log"); ax[0].set_xlabel("$u=\\pi/n_L$  ($\\rho=n_L u f_T$ $T$ states per step; $f_T=0.3$)"); ax[0].set_ylabel("qubit saving $1-Q_{\\mathrm{C}}/Q_{\\mathrm{A}}$ (%)")
ax[0].set_title("qubits saved by cultivation ($d_{\\max}=11$)"); ax[0].legend(fontsize=7, loc="upper left"); ax[0].grid(alpha=.3); ax[0].set_ylim(0, 100)
SENS = ("stage", "stage_route", "stage_growth1", "stage_growthd", "stage_route_growthd")
for key, col, lab in (("0.001_11", "C0", "$p=10^{-3}$"), ("0.0005_11", "C2", "$p=5\\times10^{-4}$")):
    r = lambda v, u: g(v, u)["nmax"][key]["ratio"]
    ax[1].plot(us, [r("stage", u) for u in us], "o-", color=col, label=f"{lab}, stage correction")
    ax[1].plot(us, [r("scale", u) for u in us], "s--", color=col, mfc="none", label=f"{lab}, pessimistic correction")
    ax[1].fill_between(us, [min(r(v, u) for v in SENS) for u in us], [max(r(v, u) for v in SENS) for u in us], alpha=.15, color=col, label=f"{lab}, routing and growth charged")
ax[1].axhline(1, color="k", ls=":", lw=.8); ax[1].set_xscale("log"); ax[1].set_yscale("log"); ax[1].set_xlabel("$u=\\pi/n_L$"); ax[1].set_ylabel("$N_{\\max}(\\mathrm{C})/N_{\\max}(\\mathrm{A})$")
ax[1].set_title("circuit-size ceiling ratio ($d_{\\max}=11$)"); ax[1].legend(fontsize=6.5, loc="upper left"); ax[1].grid(alpha=.3, which="both")
ax[1].set_yticks([0.2, 0.5, 1, 2, 3]); ax[1].set_yticklabels(["0.2", "0.5", "1", "2", "3"]); ax[1].set_ylim(0.13, 6)
fig.tight_layout(); fig.savefig(RESULTS / "fig10_unified_clock.png", dpi=300); print("fig10 written")
