# qmagic — reproduction package

Code, data and reproduction scripts for

> H. Wakaura and T. Tanimae, *Qubit cost of magic-state cultivation with logical error mitigation in early
> fault-tolerant surface-code machines* (2026).

The package evaluates four early fault-tolerant architectures (15-to-1 distillation or magic-state cultivation, each
with and without logical probabilistic error cancellation) on one clock, one noise calibration and one precision target.
It reports physical qubits, sampling overhead and the largest feasible circuit size. Scope: rotated surface codes with
lattice surgery, uniform depolarizing circuit noise, estimation of expectation values.

## Quick start

```bash
python3 -m venv .venv && source .venv/bin/activate      # Python >= 3.11
pip install -e ".[dev]"                                 # numpy, scipy, matplotlib, stim 1.16.0, pymatching 2.4.0, pytest
python -m pytest -m quick                               # 15 regression tests, < 1 min
qmagic-reproduce --verify                               # ~10 s: regenerates every table and figure into results/
                                                        # and checks the JSON/CSV outputs against reference/
```

`qmagic-reproduce` runs, in order: the measured logical-error table, the factory reconstruction, the parameter-set
builders, the claims and sensitivity sweeps, the acceptance-cutoff and hardware analyses, the exact error-cancellation
calculation, the gap analysis, the noise-model comparison, the figures, and the manuscript numbers. All of it is
deterministic. `--verify` compares the regenerated files with `reference/` to a relative tolerance of 1e-6 (1e-2 for the
optional Monte-Carlo demo, `--mc`).

## What is shipped and what is regenerated

| directory | content |
|---|---|
| `data/` | inputs that the fast path does not regenerate: surface-code memory calibrations (`calib_gidney.json`, `calib_deep.json`, `calib_noise.json`; Stim + PyMatching), cultivation reruns of Gidney's circuits (`cultiv_*.json`), Gidney 2024 published statistics (`gidney2024_*`), exact T-gate cultivation-stage error rates of Wan and Zapirain, arXiv:2609.18922 (`wan_zapirain_2026.json`), base parameter set `params.json` |
| `reference/` | the outputs of `qmagic-reproduce` as used in the paper (JSON/CSV), for `--verify` |
| `results/` | created by `qmagic-reproduce` (ignored by git); `results/paper/` holds every number and table body printed in the paper |
| `src/qmagic/` | the package; `tests/` pytest; `scripts/fetch_external.sh`, `patches/` for the slow paths |

Where the paper's numbers come from:

| paper | file(s) |
|---|---|
| every number in the text | `numbers.json` (written by `qmagic.paper_numbers` from the files below) |
| Table 2 (logical error), Appendix C | `pL_table.json` from `data/calib_gidney.json` + `data/calib_deep.json` |
| Fig. 1 (error floors) | `params_v3.json`, `factories_2level_d11.json`, `data/wan_zapirain_2026.json` |
| Fig. 2 (cultivation at small escape distance) | `data/cultiv_*.json`, `data/gidney2024_fig1_points.json` |
| Tables 3–6, Fig. 3 | `claims_v3.json`, `sensitivity_v3.json` (from `params_v3.json`) |
| Appendix D (sensitivity table) | `sensitivity_v3.json` |
| Hardware mapping, Appendix E (noise-model table) | `q3_nmax.csv`, `noise_costmodel.json` |
| Appendix F (decoder gap), Appendix G (bias of error cancellation) | `gap_qem.json`, `pec_exact.json` |

## Two inputs that need care

**Data-block error.** `pL_table.json` holds the measured logical error per patch per `d` cycles. A point is used as
measured if it has at least 30 logical errors. Larger distances are extrapolated from the three largest measured
distances of the same physical error rate. A global three-parameter fit is kept only for the hardware mapping; its
residuals are in the same file.

**Cultivation error.** Circuits with T gates cannot be sampled by stabilizer simulation. The sampled values in
`data/cultiv_*.json` replace each T gate by an S gate and double the error. `build_params_v3` corrects them with the
exact T-gate error of the cultivation stage from Wan and Zapirain: `eps_T_stage` adds the difference for that stage
(main results), `eps_T_scale` multiplies the whole sampled error (pessimistic), `eps_T_proxy` is the value as sampled.
`data/cultiv_enum_d{3,5}.json` are enumerations of the S-gate proxy circuit. They are not T-gate error rates and are
not used by the model.

## Slow paths (optional)

The calibrations and cultivation reruns in `data/` were produced with Gidney's open-source cultivation repository
(Apache-2.0), which is not vendored. To re-run them:

```bash
pip install -e ".[dev,cultivation]"     # + sinter 1.16.0, chromobius 1.1.1
scripts/fetch_external.sh               # clones Strilanc/magic-state-cultivation @ 871e68f and applies patches/
qmagic-reproduce --calib                # shallow memory calibrations with sequential stopping (~6 min on 16 cores);
                                        # writes results/calib_*_seq.json next to the shipped data/calib_*.json
qmagic-reproduce --cultiv               # cultivation under biased noise (hours)
python -m qmagic.cultivation <circuit_dir> <out.json> <shots> [procs]   # cultivation sampling runner (d2 scans, deep runs)
python -m pytest -m external            # import check against the external repository
```

`--calib` does not re-run the deep calibration points of `data/calib_deep.json`. Those took about 90 core-hours:

```bash
python -m qmagic.calib_deep a 100 1.5e9 0.001,11 0.0005,9   # 16 core-hours
python -m qmagic.calib_deep b 100 2e9 0.001,13              # 45 core-hours
python -m qmagic.calib_deep c 100 4e9 0.0005,11             # 32 core-hours
```

Arguments are a tag, the target number of logical errors per basis, the shot limit per basis, and the `p,d` points.
Each writes `results/calib_deep_<tag>.json`. Statistical reruns produce new samples; they match `data/` only within
Poisson error.

## Model in one paragraph

Logical error per patch and per $d$ cycles $p_L(p,d)$ from the measured table; data block $2d^2 r_{\rm route} n_L$
qubits; $\pi$ logical operations per step of $d$ cycles, $u=\pi/n_L$; data-block error
$\varepsilon_{\rm data}=(N/\pi)\,\kappa n_L p_L$; $T$ demand $\rho=\pi f_T$ per step; consumed $T$ error
$\varepsilon_T=\varepsilon_{\rm exit}+0.75p_L(d)$; sources limited to $d_{\rm req}\le d_{\max}$. No mitigation:
$\Sigma\le\delta=10^{-3}$. Mitigation: $\Gamma=e^{4\Sigma}\le100$ and bias $\Sigma r\le\delta$ with $r=10^{-3}$; option
`bias_rT` adds $N f_T\varepsilon_T r_T\le\delta$ when the $T$-state error is known only to $r_T$. Options
`kappa_mode="route"` and `growth_rounds` charge routing-patch error and the growth of $d_2<d$ cultivated states. See
`src/qmagic/costmodel.py`.

## License

MIT (code and data of this package). The external cultivation repository is Apache-2.0 and is fetched, not redistributed.
