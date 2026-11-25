"""Fig. 1 of the paper: exit error of 15-to-1 distillation versus its largest patch distance d_X (Litinski model on the measured
p_L table), the exact T-gate error of the cultivation stage (Wan & Zapirain 2026, horizontal lines), and the lowest sampled
end-to-end cultivation errors at escape distance d2 (S-gate proxy, corrected to T gates by the stage difference), with 95%
Poisson intervals. Filled: this work; open: published d2=15 data of Gidney et al."""
from qmagic.paths import find, RESULTS
import json, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.stats import chi2
prm = json.load(open(find("params_v3.json")))
fig, ax = plt.subplots(1, 2, figsize=(10, 4.2), sharey=True)
for i, p in enumerate((1e-3, 5e-4)):
    dx = sorted({s["d_req"] for s in prm["distill"] if s["p"] == p}); fl = [min(s["eps_T"] for s in prm["distill"] if s["p"] == p and s["d_req"] == d) for d in dx]
    ax[i].semilogy(dx, fl, "k-o", label="15-to-1 distillation, best $(d_X,d_Z,d_m)$, exit error")
    for d1, mk, col in ((3, "s", "C0"), (5, "^", "C3")):
        c = prm["proxy_correction"][f"d{d1}_{p:g}"]
        ax[i].axhline(c["T_stage"], color=col, ls=":", lw=1.5, label=f"$d_1={d1}$ cultivation stage, exact $T$-gate error [Wan & Zapirain]")
        best = {}
        for s in prm["cultivate"]:
            if s["p"] != p or s.get("d1") != d1 or s.get("errors", 0) < 4: continue
            ours = "_d2=" in s["name"]
            if d1 == 5 and ours and "deep" not in s["name"]: continue            # shallow d1=5 runs do not resolve a floor
            k = (s["d2"], ours)
            if k not in best or s["eps_T"] < best[k]["eps_T"]: best[k] = s
        for (d2, ours), s in sorted(best.items()):
            n = s["errors"]; lo, hi = chi2.ppf(0.025, 2 * n) / (2 * n), chi2.ppf(0.975, 2 * n + 2) / (2 * n)
            e = s["eps_T_proxy"]; y = s["eps_T"]; err = [[e * (1 - lo)], [e * (hi - 1)]]
            ax[i].errorbar([d2], [y], yerr=err, fmt=mk, color=col, ms=8, capsize=3, mfc=col if ours else "none",
                           label=(f"$d_1={d1}$ cultivation, lowest sampled error at escape distance $d_2$" if ours and d2 == min(k[0] for k in best if k[1]) else None))
    ax[i].axvline(11, color="gray", ls=":", lw=1); ax[i].text(11.3, 0.60, "$d_{\\max}=11$", fontsize=8, color="gray", transform=ax[i].get_xaxis_transform()); ax[i].set_xticks(dx)
    ax[i].set_title(f"$p={p:g}$"); ax[i].set_xlabel("largest patch distance ($d_X$ for distillation, $d_2$ for cultivation)"); ax[i].grid(alpha=.3, which="both")
ax[0].set_ylabel("error per output $|T\\rangle$ state (exit, before consumption)"); ax[0].legend(fontsize=6.5, loc="upper right")
fig.tight_layout(); fig.savefig(RESULTS / "fig11_floors.png", dpi=300); print("fig11 written")
