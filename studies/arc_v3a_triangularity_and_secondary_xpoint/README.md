# ARC V3A: can the coils change triangularity or move a secondary X-point?

How do the PF/CS coils respond when you change the plasma's shape targets? Can they change
triangularity, or move a secondary X-point, while keeping the other divertor targets? These are static
free-boundary TokaMaker equilibria with fixed plasma profiles; there is no TORAX coupling yet.

- **[picture_guide.html](picture_guide.html)**: an interactive picture guide to what we did and found, written for
  anyone (age 10 and up). Download it and open it in a browser; GitHub shows HTML files as source code,
  not as a page.
- **[REPORT.md](REPORT.md)**: the full write-up: setup, observations, interpretation, open questions.
- **[HOW_THE_NOTEBOOK_WORKS.md](HOW_THE_NOTEBOOK_WORKS.md)**: how the original notebook sets up and solves
  the equilibrium.

## Where the inputs come from

- **Equilibrium setup:** F. Sheehan's (Columbia) ARC V3A shape-testing notebook
  `2026-06-08_ARC_diverted_shape_testing.ipynb`, shared for this study and published here with
  permission. That covers the X-point and strike-point targets, constraint weights, coil bounds,
  regularization and profiles. The notebook itself is not committed; the script reproduces its
  equilibrium exactly.
- **Mesh:** the public OFT ARC V3A mesh (`src/examples/TokaMaker/ARC/arc_v3a_mesh.h5`, OFT PR #298),
  read from the OFT fork and checked by SHA-256.

## Cases

| Case | Change from baseline |
|---|---|
| `A_baseline_delta065` | none (requested δ = 0.65) |
| `A_baseline_delta065_4threads` | same, run on 4 threads (reproducibility check only) |
| `B_delta060` | requested δ = 0.60 |
| `C_secondary_xpoint_out_5cm` | secondary X-point targets R = 5.00 → 5.05 m |

## Files

| Path | What it is |
|---|---|
| `results/summary_by_case.md` | Requested vs achieved shape, X-points, ψ̂ at each null, strike points, constraint residuals, convergence |
| `results/coil_currents_by_case.md` | Coil currents and differences, kA per turn and MA-turns |
| `results/<case>.json` | Every recorded number for each case |
| `results/shape_measured_from_boundary.json` | δ and κ measured from the traced boundary |
| `results/coil_cancellation_check.json` | Vacuum-field check of how "visible" the coil-current changes are at the plasma |
| `figures/A_vs_B_triangularity_separatrix.png` | Experiment 1: separatrix for δ 0.65 vs 0.60, full view and divertor zoom |
| `figures/A_vs_C_secondary_xpoint_separatrix.png` | Experiment 2: separatrix with the secondary X-points at R = 5.00 vs 5.05 m |
| `figures/coil_currents_by_case.png` | Coil currents by coil, all three cases |
| `figures/<case>_equilibrium.png` | Each equilibrium in the original notebook's plot style |
| `code/` | `run_equilibria.py` (runs the cases), `make_figures_and_tables.py`, `check_coil_cancellation.py`, `build_picture_guide.py` (+ `picture_guide_template.html`) |

`A_baseline_delta065_4threads.json` changes in its last digits each time it is rerun, because multi-threaded OFT is not
bit-reproducible. All other results repeat exactly.

## How to rerun

From the repo root, with the repo's `.venv` (OFT 26.9) and the OFT fork at
`~/Documents/code/OpenFUSIONToolkit` (or pass `--oft-root PATH`):

```
.venv/bin/python studies/arc_v3a_triangularity_and_secondary_xpoint/code/run_equilibria.py       # 4 cases + figures + tables, ~30 s
.venv/bin/python studies/arc_v3a_triangularity_and_secondary_xpoint/code/check_coil_cancellation.py   # coil diagnostic, ~10 s
.venv/bin/python studies/arc_v3a_triangularity_and_secondary_xpoint/code/build_picture_guide.py        # rebuild picture_guide.html from the results
```

Flux arrays and solver logs go to `work/arc_v3a_triangularity_and_secondary_xpoint/` (gitignored). The guide and figures need them,
so run `run_equilibria.py` first.
