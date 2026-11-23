# SPARC: how equilibria respond to current, pressure and shape

Static (fixed-time, free-boundary) TokaMaker equilibria of the public SPARC Primary Reference
Discharge, changing one input at a time: plasma current, pressure, X-point position or elongation.
A side exploration, separate from the ARC V3A work, using only public data.

- **[REPORT.md](REPORT.md)**: findings and analysis.
- **[SECONDARY_NULLS_REPORT.md](SECONDARY_NULLS_REPORT.md)**: follow-up on the extra field nulls that
  appear at lower current or higher pressure.

## Cases

| Case | Change from baseline |
|---|---|
| `baseline` | none (PRD: Ip 8.7 MA, p_axis 2.6 MPa) |
| `current_minus20pct`, `current_minus10pct`, `current_plus10pct`, `current_plus20pct` | Ip target × 0.8, 0.9, 1.1, 1.2 |
| `pressure_minus20pct`, `pressure_plus20pct` | on-axis pressure target × 0.8, 1.2 |
| `xpoints_inward_5cm`, `xpoints_outward_5cm` | both X-point targets moved 5 cm in R |
| `xpoints_toward_midplane_5cm`, `xpoints_away_from_midplane_5cm` | both X-point targets moved 5 cm in Z |
| `upper_xpoint_up_5cm` | upper X-point target only, 5 cm up |
| `elongation_minus5pct`, `elongation_plus5pct` | boundary and X-point targets squashed / stretched 5 % in Z |

## Starting point and data provenance (all public)

- OFT example `src/examples/TokaMaker/SPARC/SPARC_public_baseline_ex.ipynb` and `sparc_public_mesh.h5`,
  added upstream in OpenFUSIONToolkit PR #317 (commit 85fb2f1). Identical to `upstream/main`
  and already present in the pinned fork commit `8fc4afe`.
- Target equilibrium: CFS SPARCPublic, `PrimaryReferenceDischarge/2 - SPARC_DN_PRD_freegs_20221013`
  (double-null PRD, Ip = 8.7 MA, B0 = 12.2 T, p_axis = 2.6 MPa). MIT-style licence; CFS asks for the
  acknowledgement in their README if used in publications.
- **No coil limits are public.** The SPARCPublic device description lists power supplies by name only,
  and its README says the coil turn split is "toy data". So no case can be checked against coil limits.
  `COIL_LIMITS` in the script is left empty. Fill it only from a public source.

## How to run

```
.venv/bin/python studies/sparc_current_pressure_and_shape/code/run_equilibria.py             # all cases, ~2.5 min
.venv/bin/python studies/sparc_current_pressure_and_shape/code/run_equilibria.py --case current_plus10pct
.venv/bin/python studies/sparc_current_pressure_and_shape/code/run_equilibria.py --summary   # table + grid plot only
```

On first run the two inputs (the OFT SPARC mesh and the SPARCPublic PRD gEQDSK) are downloaded
into `work/sparc_current_pressure_and_shape/inputs/` (gitignored) and checked against fixed SHA-256 sums. They are not
committed here.

Each case runs in its own process (one TokaMaker per process). `nthreads=1`: the baseline JSON is
bit-identical on rerun. OFT 26.9, Python 3.12.13 (the ARC repo's `.venv`).

Outputs: `results/summary_by_case.md` (main table), `results/coil_currents_by_case.md`,
`results/<case>.json`, `figures/<case>_equilibrium.png` (each equilibrium vs the baseline boundary and
X-points) and `figures/all_cases_boundary_vs_baseline.png` (all boundaries).
Equilibrium arrays (`*_eq.npz`, ~4 MB each) and solver logs go to `work/sparc_current_pressure_and_shape/`.

`fit_scalings.py` reads the results and writes `results/scaling_fits.json` (the fitted scalings quoted in
the report) and `figures/q95_and_beta_p_scalings.png`:

```
.venv/bin/python studies/sparc_current_pressure_and_shape/code/fit_scalings.py
```

`find_secondary_nulls.py` re-solves each case and finds every field null on the mesh. It uses the divertor
metrics from OFT PR #403, fetched as a single file at a fixed commit (pins unchanged). It writes
`results/secondary_nulls.{json,md}` and `figures/secondary_nulls_upper_divertor.png` (~3.5 min). Findings: [SECONDARY_NULLS_REPORT.md](SECONDARY_NULLS_REPORT.md).

```
.venv/bin/python studies/sparc_current_pressure_and_shape/code/find_secondary_nulls.py
```

## Setup kept fixed in every case

Same mesh, the same F0, and P' and FF' profile *shapes* from the gEQDSK. Ip and p_axis are hard targets.
The gEQDSK LCFS is used as isoflux targets and two mirrored X-points as saddle targets. VS1U/VS1L act as
a virtual VSC. Coils are regularised toward zero current (VS coils weighted more heavily). No coil bounds
are set. Only the one listed input changes.

## Baseline vs the public FreeGS equilibrium

| | TokaMaker | gEQDSK |
|---|---|---|
| Ip [MA] | 8.70 | 8.70 |
| q95 | 3.17 | 3.01 |
| q0 | 0.95 | 0.94 |
| kappa | 1.96 | 1.95 |
| R_axis [m] | 1.891 | 1.890 |

Other baseline values: delta 0.571, l_i 1.28, beta_p 0.361, beta_N 0.81. X-points at
(1.524, +1.124) m and (1.523, −1.123) m. Converged in 20 iterations, diverted double null.
The ~5 % q95 difference is in the upstream example itself, not caused by this scan.

## Results

See `results/summary_by_case.md`. All 14 cases converged. None hit a coil limit, but that is only because there
are no public limits to check against.

Things to know when reading the table:

- **OFT's `beta_pol` is in percent.** The scan stores beta_p as a fraction.
- **delta and kappa from `get_stats` have about ±0.02 noise in this double null.** In the baseline
  the X-points mirror each other to within 1 mm, yet deltaU = 0.590 and deltaL = 0.552.
  So changes in delta smaller than about 0.02 (all the Ip and pressure cases) are not meaningful.
  The X-point coordinates are the reliable shape measure.
- **"Active" X-point upper or lower** flips on sub-millimetre differences in a balanced double
  null, so it is not meaningful except in `upper_xpoint_up_5cm`, where the bias is deliberate.
- **Extra field nulls.** At lower Ip (−10 %, −20 %) or higher pressure (+20 %), TokaMaker reports an
  extra field null at about (1.94, ±1.49) m. It is a saddle of the coil field about 10 cm behind the
  first wall, next to the DV1 coil, and 51 cm from the LCFS. It exists in every case and is reported
  only when it moves into the vessel conductor. See [SECONDARY_NULLS_REPORT.md](SECONDARY_NULLS_REPORT.md). The table
  reports the X-point nearest each target.
- **Coil currents are a static snapshot.** In a fixed-time solve the CS state is only what the
  regularisation (pull toward 0 A) prefers. CS changes here say nothing about flux swing.

## Acknowledgement

The information, data, or work presented herein builds on the SPARC primary reference discharge and
X-point target discharge data provided by Commonwealth Fusion Systems.
