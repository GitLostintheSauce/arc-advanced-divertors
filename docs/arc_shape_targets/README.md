# ARC V3A shape-target study

How do the PF/CS coils respond when you change the plasma's shape targets? Can they change
triangularity, or move a secondary X-point, while keeping the other divertor targets? These are static
free-boundary TokaMaker equilibria with fixed plasma profiles; there is no TORAX coupling yet.

- **[guide.html](guide.html)**: an interactive picture guide to what we did and found, written for
  anyone (age 10 and up). Download it and open it in a browser. GitHub shows HTML files as source code,
  not as a page.
- **[REPORT.md](REPORT.md)**: the full write-up: setup, observations, interpretation, open questions.
- **[notebook_walkthrough.md](notebook_walkthrough.md)**: how the original notebook sets up and solves
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
| `A_baseline` | none (requested δ = 0.65) |
| `A_baseline_4thr` | same, run on 4 threads (reproducibility check only) |
| `B_delta060` | requested δ = 0.60 |
| `C_secR505` | secondary X-point targets R = 5.00 → 5.05 m |

## Files

| Path | What it is |
|---|---|
| `results/summary_table.md` | Requested vs achieved shape, X-points, ψ̂ at each null, strike points, constraint residuals, convergence |
| `results/coil_table.md` | Coil currents and differences, kA per turn and MA-turns |
| `results/<case>.json` | Every recorded number for each case |
| `results/measured_shape.json` | δ and κ measured from the traced boundary |
| `results/coil_null_space.json` | Vacuum-field check of how "visible" the coil-current changes are at the plasma |
| `../../figures/arc_shape_targets/` | Equilibrium comparisons, coil bar chart, one plot per case |
| `../../scripts/arc_shape_targets/` | `shape_study.py`, `analyse.py`, `coil_null_space.py`, `make_guide.py` |

`A_baseline_4thr.json` changes in its last digits each time it is rerun, because multi-threaded OFT is not
bit-reproducible. All other results repeat exactly.

## How to rerun

From the repo root, with the repo's `.venv` (OFT 26.9) and the OFT fork at
`~/Documents/code/OpenFUSIONToolkit` (or pass `--oft-root PATH`):

```
.venv/bin/python scripts/arc_shape_targets/shape_study.py       # 4 cases + figures + tables, ~30 s
.venv/bin/python scripts/arc_shape_targets/coil_null_space.py   # coil diagnostic, ~10 s
.venv/bin/python scripts/arc_shape_targets/make_guide.py        # rebuild guide.html from the results
```

Flux arrays and solver logs go to `work/arc_shape_targets/` (gitignored). The guide and figures need them,
so run `shape_study.py` first.
