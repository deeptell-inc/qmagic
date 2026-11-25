"""Low-weight enumeration of undetected logical error sets of the inject+cultivate stage in the S-GATE PROXY circuit
(Gidney's ErrorEnumerationReport, external repo). The result is the proxy error of the stage, truncated at `max_weight`:
converged for d1=3 at weight 5; leading order only for d1=5. It is NOT the error of the circuit with actual T gates, which has
lower fault distance (2 and 3) and is given exactly by Wan & Zapirain, arXiv:2609.18922 (data/wan_zapirain_2026.json).
Kept as a cross-check of the published proxy statistics (d1=3: 3.42e-7 at p=1e-3). "floor_T" in the output is 2x the proxy value.
    python -m qmagic.cultiv_enum [max_weight=5] [d1=3] -> results/cultiv_enum_d<d1>.json   (d1=5, weight 5: 2797 mechanisms, ~45 min)
"""
from qmagic.paths import RESULTS, find
import json, sys, time
from qmagic import noise as NM
cultiv = NM.cultiv; gen = NM.gen
from cultiv import ErrorEnumerationReport
def main(max_weight=5, d1=3):
    out = {}; cache = {}
    meas = {s["p"]: {} for s in json.load(open(find("params_v3.json")))["cultivate"] if s.get("d1") == d1}
    for s in json.load(open(find("params_v3.json")))["cultivate"]:
        if s.get("d1") == d1 and not s.get("scaled") and s.get("errors", 0) >= 4:
            meas[s["p"]][s["d2"]] = min(meas[s["p"]].get(s["d2"], 1.0), s["eps_T"])
    for p in (1e-3, 5e-4):
        c = cultiv.make_inject_and_cultivate_circuit(inject_style="unitary", dcolor=d1, basis="Y")
        t0 = time.time(); rep = ErrorEnumerationReport.from_circuit(c, max_weight=max_weight, noise=p, cache=cache); sec = time.time() - t0
        byw = {int(w): float(v) for w, v in rep.distance_to_heralded_error_rate.items()}
        rec = dict(p=p, max_weight=max_weight, keep_rate=rep.keep_rate, heralded_error_rate_proxy=rep.heralded_error_rate,
                   floor_T=2 * rep.heralded_error_rate, by_weight_proxy=byw, n_logical_sets={int(w): sum(1 for e in rep.logical_errs if len(e.src_errors) == w) for w in byw},
                   seconds=sec, measured_floor_T_by_d2=meas.get(p, {}))
        out[f"{p:g}"] = rec
        print(f"p={p:g}: keep={rep.keep_rate:.3f}  exact proxy error (w<={max_weight}) = {rep.heralded_error_rate:.3e}  x2 -> {rec['floor_T']:.3e}  "
              f"| by weight: " + ", ".join(f"w{w}: {v:.2e} ({rec['n_logical_sets'][w]} sets)" for w, v in sorted(byw.items()) if v > 0)
              + f" | measured floors (x2, >=4 errors) by d2: " + ", ".join(f"d2={d}: {e:.2e}" for d, e in sorted(meas.get(p, {}).items())) + f" | {sec:.1f} s")
    json.dump(out, open(RESULTS / f"cultiv_enum_d{d1}.json", "w"), indent=1)
if __name__ == "__main__": main(int(sys.argv[1]) if len(sys.argv) > 1 else 5, int(sys.argv[2]) if len(sys.argv) > 2 else 3)
