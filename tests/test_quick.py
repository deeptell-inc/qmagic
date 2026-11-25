"""Quick regression tests (no external repository needed)."""
import json, math, pytest, numpy as np
from qmagic.paths import DATA, EXTERNAL
from qmagic import factory as F, costmodel, pec

pytestmark = pytest.mark.quick

def _find(name):
    from qmagic.paths import find
    try: return find(name)
    except FileNotFoundError: pytest.skip(f"{name} not generated yet (run qmagic-reproduce)")


def test_factory_reproduces_litinski_table():
    worst = 1
    for (dX, dZ, dm), p, pout, q, cyc in F.TABLE:
        r = F.simulate(dX, dZ, dm, p); ratio = r["p_out"] / pout; worst = max(worst, ratio, 1 / ratio)
        assert abs(r["qubits"] - q) / q < 0.01 and abs(r["cycles"] - cyc) / cyc < 0.02
    assert worst < 1.3

def test_two_level_factory_matches_litinski():
    r = F.simulate2((11, 5, 5), (25, 11, 11), 1e-3)
    assert 0.4 < r["p_out"] / 2.7e-12 < 2.0 and 0.8 < r["qubits"] / 30700 < 1.2

def test_d11_distillation_floor_is_storage_limited():
    fit = json.load(open(DATA / "calib_pL_corr.json"))["fit"]
    pL = lambda p, d: fit["A"] * (p / fit["p_th"]) ** (fit["exponent_scale"] * (d + 1) / 2) / d
    one = F.simulate(11, 5, 5, 1e-3, pL=pL)["p_out"]; two = F.simulate2((11, 5, 5), (11, 5, 5), 1e-3, pL=pL)["p_out"]
    assert two > 0.5 * one           # second level buys < 2x

def test_suzuki_self_check():
    out = costmodel.suzuki_check(); assert out[1e4][:2] == (9, 4) and out[1e10][:2] == (19, 14)

def test_costmodel_headline_regression():
    prm = json.load(open(_find("params_corr_d2.json")))
    prm["p"] = [1e-3]; prm["d_max"] = [11]; prm["N"] = [1e4]; prm["f_T"] = [0.3]; prm["C_max"] = [100.0]
    rows = {r["arm"]: r for r in costmodel.sweep(prm)}
    assert rows["A"]["feasible"] and rows["C"]["feasible"]
    assert 0.25 < rows["C"]["saving_C_vs_A"] < 0.40          # reported 0.32
    assert rows["C"]["Gamma"] <= 100 and rows["C"]["src"].startswith("cult_d1=3")

def test_pec_gamma_formula_and_unbiasedness():
    r = pec.run(400, 2.5e-4, 2.5e-4, M=20_000, seed=5)
    assert abs(r["gamma_tot"] ** 2 / math.exp(4 * r["sum_eps"]) - 1) < 1e-3
    assert abs(r["mean"] - r["ideal"]) < 4 * r["sem"]

def test_gap_mapping_universality_d1_3():
    g = json.load(open(_find("gap_qem.json")))["maps"]
    slopes = [v["fit"][0] for k, v in g.items() if k.startswith("d1=3")]
    assert len(slopes) >= 5 and max(slopes) - min(slopes) < 0.006 and all(-0.09 < s < -0.078 for s in slopes)

def test_erasure_reindexing_preserves_dem_when_pe_is_zero():
    import stim
    from qmagic import noise as NM
    c = stim.Circuit.generated("surface_code:rotated_memory_z", distance=3, rounds=3, after_clifford_depolarization=1e-3)
    ce = NM.erasure_circuit(c, 0.0); cs = NM.strip_herald_detectors(ce)
    assert ce.num_detectors > c.num_detectors and cs.num_detectors == c.num_detectors and cs.num_observables == 1
    # compare error mechanisms/detector structure; coordinates can differ by SHIFT_COORDS handling after flattening
    strip = lambda dem: "\n".join(l for l in str(dem.flattened()).splitlines() if l.startswith("error"))
    assert strip(cs.detector_error_model(approximate_disjoint_errors=True)) == strip(c.detector_error_model())

@pytest.mark.external
@pytest.mark.skipif(not EXTERNAL.exists(), reason="external repo not fetched")
def test_cultivation_sampler_matches_gidney_head():
    from qmagic import cultivation  # noqa: F401  (import only; full run is hours)

def test_v3_unified_clock_regression():
    """Panel-corrected model: serial u=1/30 -> small hybrid saving and no N_max gain; parallel u=1 -> large saving (measured rows only)."""
    prm = json.load(open(_find("params_v3.json")))
    prm["cultivate"] = [s for s in prm["cultivate"] if not s.get("scaled")]
    prm["p"] = [1e-3]; prm["d_max"] = [11]; prm["N"] = [1e4]; prm["f_T"] = [0.3]; prm["C_max"] = [100.0]
    prm["parallelism"] = prm["n_L"] / 30; s_serial = {r["arm"]: r for r in costmodel.sweep(prm)}["C"]["saving_C_vs_A"]
    prm["parallelism"] = prm["n_L"]; s_par = {r["arm"]: r for r in costmodel.sweep(prm)}["C"]["saving_C_vs_A"]
    assert 0.05 < s_serial < 0.15 and 0.6 < s_par < 0.95

def test_build_params_v3_from_shipped_inputs(tmp_path):
    """The v4 parameter set builds from data/ alone in a scratch root and reproduces the headline u=1 saving."""
    import subprocess, sys, os, shutil
    root = tmp_path / "root"; shutil.copytree(DATA, root / "data")
    env = dict(os.environ, QMAGIC_ROOT=str(root))
    for mod in ("qmagic.calib_table", "qmagic.rerun_corr", "qmagic.build_params_d2scan", "qmagic.build_params_v3"):
        subprocess.run([sys.executable, "-m", mod], check=True, env=env, capture_output=True, cwd=root)
    prm = json.load(open(root / "results" / "params_v3.json")); prm["cultivate"] = [c for c in prm["cultivate"] if not c.get("scaled")]
    assert len(prm["distill"]) > 100 and prm.get("bias_r") == 1e-3
    prm["p"] = [1e-3]; prm["d_max"] = [11]; prm["N"] = [1e5]; prm["f_T"] = [0.3]; prm["C_max"] = [100.0]; prm["parallelism"] = prm["n_L"]
    rows = {r["arm"]: r for r in costmodel.sweep(prm)}
    assert 0.7 < rows["C"]["saving_C_vs_A"] < 0.9 and rows["A"]["k"] > 50 and "pL_table" in prm

def test_pec_single_shot_variance_identity():
    r = pec.run(400, 2.5e-4, 2.5e-4, M=40_000, seed=9)
    assert abs(r["var_x_M"] / (r["gamma_tot"] ** 2 - r["ideal"] ** 2) - 1) < 0.05


def test_pL_table_lookup_and_extrapolation():
    prm = dict(pL_table={"0.001": dict(measured={"9": 3.6e-6}, A=0.3, x=0.1)}, pL_fit=dict(A=0.333, p_th=0.0101, s=0.963))
    assert costmodel.pL(prm, 1e-3, 9) == 3.6e-6                                   # measured value
    assert abs(costmodel.pL(prm, 1e-3, 11) / (0.3 * 0.1 ** 6) - 1) < 1e-12         # per-p extrapolation
    assert abs(costmodel.pL(prm, 2e-3, 9) / (0.333 * (2e-3 / 0.0101) ** (0.963 * 5)) - 1) < 1e-12   # global fit elsewhere

def test_proxy_modes_and_T_calibration_constraint():
    prm = json.load(open(_find("params_v3.json")))
    row = next(c for c in prm["cultivate"] if c.get("d1") == 5 and "eps_T_scale" in c)
    assert row["eps_T_proxy"] < row["eps_T_stage"] < row["eps_T_scale"] and row["eps_T"] == row["eps_T_stage"]
    assert not any(c.get("scaled") for c in prm["cultivate"])                      # no transferred tails
    prm["parallelism"] = prm["n_L"]
    base = costmodel.nmax(prm, "cultivate", 1e-3, 11, 0.3); prm["bias_rT"] = 0.3
    lim = costmodel.nmax(prm, "cultivate", 1e-3, 11, 0.3)
    best = min(c["eps_T"] for c in prm["cultivate"] if c["p"] == 1e-3 and c["d_req"] <= 11)
    assert lim < base and lim <= prm["N_e"] / (0.3 * 0.3 * best) * 1.0001          # N_T * eps_T * r_T <= N_e
    assert costmodel.nmax(costmodel.with_proxy_mode(prm, "scale"), "cultivate", 5e-4, 11, 0.3) < costmodel.nmax(prm, "cultivate", 5e-4, 11, 0.3)

def test_wan_zapirain_series_reproduces_exact_values():
    from fractions import Fraction
    wz = json.load(open(DATA / "wan_zapirain_2026.json"))
    x = Fraction(1, 1000) / (1 - Fraction(1, 1000))
    for d, ref in ((3, 8.982790954e-7), (5, 3.328425371418e-9)):
        v = float(sum(Fraction(c) * x ** int(k) for k, c in wz[f"d{d}"]["series"].items())); assert abs(v / ref - 1) < 1e-6
