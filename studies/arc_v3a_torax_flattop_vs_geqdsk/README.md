# ARC V3A: TokaMaker+TORAX flattop vs the public gEQDSK

How does the ARC V3A flattop from the coupled TokaMaker+TORAX workflow compare with the public ARC V3A
reference equilibrium (gEQDSK)? Same machine, same plasma current and field. Different ways of getting
the profiles and the shape.

- **[REPORT.md](REPORT.md)**: results, what they mean, open questions.
- **[results/comparison_table.md](results/comparison_table.md)**: every number side by side.

## The two equilibria

| | Public gEQDSK | TokaMaker+TORAX flattop |
|---|---|---|
| Source | Zenodo record [19498373](https://zenodo.org/records/19498373), file `geqdsk-ARCv3a` (FreeGS, 08/03/2025; Hillesheim et al., *J. Plasma Phys.* 2026) | F. Sheehan's public notebook `ARC_V3A/TMTX_project_parts_1_2_ARC_V3A_flattop_v1.ipynb`, [fs2819/TokaMaker_TORAX_project](https://github.com/fs2819/TokaMaker_TORAX_project) at commit `bf5e8d1` |
| Profiles | Fixed in the file | Predicted by TORAX (QLKNN transport, 4.5 keV pedestal, 50 MW ICRH, gas puff, no sawteeth), iterated with TokaMaker to steady state |
| Shape | FreeGS's own coil solution | Targets from a = 1.18 m, κ = 1.8, δ = 0.65 (`xpoints_from_moments`, `create_isoflux_xpts`) |
| Mesh / machine | ARC V3A | Public OFT ARC V3A mesh (same file as the notebook downloads) |

## How it is run

1. `code/run_flattop.py` downloads the notebook at the pinned commit (SHA-256 checked) and runs
   Section 1, the single flattop run, unchanged except for four small patches listed in the script:
   - skip the pip-install cell;
   - take `Hmode_profiles` from OFT instead of `omfit_classes`;
   - switch the Z_eff sweep off;
   - add one cell that saves the final equilibrium as a gEQDSK on the same 129 × 129 grid as the
     Zenodo file, plus the TokaMaker and TORAX results.

   It takes about 4 minutes.
2. `code/compare_to_geqdsk.py` measures both gEQDSK files with the same functions: shape, X-points,
   a search for an outer-divertor null, profiles, and volume integrals for β and l_i. As a check, it
   compares its TM+TX numbers against TokaMaker's own `get_stats()`. They agree within about 1 %.

### TORAX 1.4.3 for this study only

The notebook's pedestal settings (`pedestal_profile_form`, `explicit_pedestal`) do not exist in
TORAX 1.4.0, which is the repo-wide pin; the notebook was run with TORAX 1.4.3. This study therefore
uses its own environment, [requirements-lock-torax143.txt](requirements-lock-torax143.txt). It is the
repo lock file with exactly two changes: `torax==1.4.3` and `eqdsk==0.9.0`.

The repo stays on 1.4.0 for now. With 1.4.3, the ITER TokaMaker+TORAX example still runs
(25/25 solves), but its L-H transition moves from 136 / 169 s to about 81 s, and the CI check fails.

```
python3.12 -m venv work/venv_torax143
work/venv_torax143/bin/pip install -r studies/arc_v3a_torax_flattop_vs_geqdsk/requirements-lock-torax143.txt
work/venv_torax143/bin/python studies/arc_v3a_torax_flattop_vs_geqdsk/code/run_flattop.py      # ~4 min
work/venv_torax143/bin/python studies/arc_v3a_torax_flattop_vs_geqdsk/code/compare_to_geqdsk.py
```

`run_flattop.py` stops with a clear message if TORAX is not 1.4.x. The executed notebook, the saved
TM+TX gEQDSK and the raw results go to `work/arc_v3a_torax_flattop_vs_geqdsk/` (gitignored).
TORAX's results vary slightly from run to run, so rerun numbers can differ in the last digits.

## Files

| Path | What it is |
|---|---|
| `results/comparison_table.md` | Side-by-side numbers, the measuring check, and TM+TX run context |
| `results/comparison.json` | All measured values |
| `figures/separatrix_tmtx_vs_geqdsk.png` | ψ̂ = 1 contours, X-points and axis of both, full view and upper-divertor zoom |
| `figures/profiles_tmtx_vs_geqdsk.png` | p, q, FF′ and p′ against normalized poloidal flux |
| `figures/torax_final_profiles.png` | TORAX's density, temperature and current profiles (the gEQDSK has none of these) |
| `code/run_flattop.py`, `code/compare_to_geqdsk.py` | The two steps above |
| `requirements-lock-torax143.txt` | The environment this study needs |
