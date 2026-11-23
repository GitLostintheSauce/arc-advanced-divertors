# DIII-D and ITER: how realistic are the FF′ and P′ profiles in OFT's examples?

## Question
OFT's DIII-D shape-helper example (`DIIID_shapehelper_ex.ipynb`) sets FF′ and P′ from simple power laws.
OFT's ITER H-mode example (`ITER_Hmode_ex.ipynb`) builds them from kinetic profiles and a bootstrap current
calculation, and it has been described as giving realistic FF′ and P′. This study asks:

- what each example's solved FF′ and P′ actually look like, in physical units;
- how far the shape-helper profiles are from a real DIII-D reconstruction;
- what happens to the shape-helper profiles under realistic operating conditions (more current, more pressure);
- which features of the ITER profiles come from the solver and which come from the chosen input profiles.

This covers steps 1–3 of the plan: metrics, reproduction of both examples, and a one-at-a-time DIII-D scan.
Swapping H-mode profiles into DIII-D, and sensitivity tests on the ITER input parameters, are not done yet.

## Metrics
All profiles are TokaMaker's solved output (`get_profiles`, `get_q`). ψ̂ = 0 is the magnetic axis and 1 the LCFS.
P′ and FF′ are per unit of TokaMaker flux (Wb/rad). Integrating P′ over the flux range recovers the axis pressure
to better than 0.01 % (`p_axis_from_pp` against `P_ax` in each case JSON), which confirms the units.

| Metric | Meaning |
|---|---|
| q₀, q₉₅, q_min | safety factor on axis, at ψ̂ = 0.95, and its minimum |
| l_i | normalized internal inductance (`get_stats`, standard normalization) |
| β_N | 100·β_t·a·B₀/I_p[MA], with β_t from the volume-averaged pressure and the vacuum field at R_geo |
| FF′(axis), min FF′, FF′ zero crossings | the value of FF′ on axis, its minimum, and where it changes sign |
| edge/core P′ peak | largest P′ at ψ̂ > 0.85 divided by the largest P′ inside it. Above 1 means an edge pedestal dominates P′. For a monotone P′ it is just P′(0.85)/P′(0) |

Reference values: the EFIT01 reconstruction of DIII-D shot 192185 (file `g192185.02440`, shipped with OFT) and
the same shot solved by TokaMaker from that file's inputs (OFT's `DIIID_baseline_ex.ipynb`). EFIT01 uses only
magnetic measurements with few profile basis functions, so it cannot resolve a pedestal. It is a fair reference
for the core and for global quantities, not for the edge.

## Inputs
From the public OpenFUSIONToolkit repository at commit `8fc4afe` (LGPL-3.0), downloaded on demand into
`work/diiid_iter_flux_function_realism/inputs/` and checked by SHA-256 (`code/common.py`): `DIIID_mesh.h5`,
`g192185.02440`, `ITER_mesh.h5`, and `bootstrap.py`.

`bootstrap.py` at `8fc4afe` includes OFT PR #316, a fix for derivative accuracy at the axis and separatrix in
`solve_with_bootstrap`. The pinned OpenFUSIONToolkit 26.9 wheel does not have it. The ITER example is run with
both versions, the second loaded as a standalone module. No pins change.

The code copies the physics cells of the notebooks at `8fc4afe` (the API the 26.9 wheel supports) and leaves out
plot cells. Upstream `main` has since moved the ITER example to `solve_bootstrap()`, which 26.9 does not have.
TokaMaker runs with 1 thread, so results are bitwise reproducible.

## Cases
**DIII-D** (`code/run_diiid_shapehelper.py`):
- `shapehelper_as_published`: the notebook unchanged. Ip 0.5 MA, axis pressure 7 kPa, FF′ = (1 − ψ̂^1.5)², P′ = 1 − ψ̂⁴.
- `efit_profiles_same_shape`: the same plasma shape and coils, with the FF′ and P′ shapes of the EFIT reconstruction.
- `efit_192185_baseline`: the baseline notebook's inverse solve with the shot's boundary, profiles, Ip, pressure and coil targets.
- One-at-a-time changes from the published case:
  - Ip from 0.75 to 1.5 MA
  - axis pressure from 15 to 120 kPa
  - each FF′ and P′ exponent
- A pressure scan at 1.2 MA, a typical DIII-D H-mode current. These cases change two inputs.
- Cases that do not converge in the notebook's default 40 nonlinear iterations are rerun with 300 (`*_maxits300`).

**ITER** (`code/run_iter_hmode.py`): the notebook as published, with `isolate_edge_jBS` False and then True. Each
is run with the wheel's `bootstrap.py` and with the `8fc4afe` version.

## How to rerun
From the repo root, in the pinned `.venv`:

```bash
.venv/bin/python studies/diiid_iter_flux_function_realism/code/run_diiid_shapehelper.py
```

```bash
.venv/bin/python studies/diiid_iter_flux_function_realism/code/run_iter_hmode.py
```

```bash
.venv/bin/python studies/diiid_iter_flux_function_realism/code/summarise.py
```

Rerun a single DIII-D case with `--case NAME`, and add `--maxits 300 --save-as NAME_maxits300` for the longer-iteration rerun.

## Outputs
- `results/cases/*.json`: one record per case (scalars and profile metrics)
- `results/summary.md`: tables
- `figures/`:
  - `diiid_published_vs_efit.png`
  - `diiid_scan_profiles.png`
  - `diiid_pressure_scan_scalars.png`
  - `iter_hmode_profiles.png`
- Profile arrays (`.npz`) and logs: `work/diiid_iter_flux_function_realism/`

Results and interpretation are in [REPORT.md](REPORT.md).
