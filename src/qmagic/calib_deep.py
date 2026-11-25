"""Deep data-block calibration points (direct measurement instead of extrapolation): Gidney memory circuit, uniform
depolarizing noise incl. idle, X and Z bases, PyMatching correlated matching; sequential stopping at `target` errors per basis.
    python -m qmagic.calib_deep <tag> <target> <maxshots_per_basis> <p,d> [<p,d> ...]   -> results/calib_deep_<tag>.json (checkpointed)
"""
from qmagic.paths import RESULTS
import json, sys
from qmagic import seqcalib, calib_gidney
def main():
    tag, target, maxshots = sys.argv[1], int(sys.argv[2]), int(float(sys.argv[3]))
    pts = [(float(a.split(",")[0]), int(a.split(",")[1]), b) for a in sys.argv[4:] for b in ("X", "Z")]
    out = RESULTS / f"calib_deep_{tag}.json"
    def save(rows): json.dump(dict(target=target, maxshots=maxshots, rows=calib_gidney.to_rows(rows)), open(out, "w"), indent=1)
    rows, _, wall = seqcalib.run_points(pts, [], calib_gidney.build, lambda c, k: None, target=target, maxshots=maxshots, procs=16,
                                        seed_base=11_000_000, chunk=5_000_000, checkpoint=save)
    for r in sorted(calib_gidney.to_rows(rows), key=lambda r: (r["p"], r["d"], r["basis"])):
        print(f"p={r['p']:g} d={r['d']:2} {r['basis']}: shots={r['shots']:.3e} err={r['errors']:4d} pL={r['pL']:.3e} {r['sec']} core-s", flush=True)
    print(f"wall {wall:.0f} s", flush=True)
if __name__ == "__main__": main()
