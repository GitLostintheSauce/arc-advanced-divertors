# Advanced divertors on ARC V3A

How hard is it to **form and hold** advanced divertor configurations (X-point target, super-X, X-divertor, snowflake) on ARC V3A with its real coils, shielding and pulse? This repo holds the study-specific work: scans, configurations, figures and the paper. It uses [TokaMaker](https://github.com/OpenFUSIONToolkit/OpenFUSIONToolkit) and the TokaMaker+TORAX workflow.

The plan, timeline and code-flow rules are in the [roadmap](docs/roadmap.md).

## Layout

| Folder | What goes here |
|---|---|
| `studies/` | One folder per study, named for what it studies (for example `arc_v3a_triangularity_and_secondary_xpoint`). Each holds its own `README.md`, `REPORT.md`, `code/`, `results/` and `figures/`. See [studies/README.md](studies/README.md) |
| `notebooks/` | Exploratory and figure-making notebooks for this study |
| `scripts/` | Repo-wide scripts: install check, OFT example runner, CI checks |
| `configs/` | Scan definitions, coil/shape targets and other inputs. Public ARC V3A data only (see below) |
| `figures/` | Repo-wide plots (the install-check equilibrium). Study figures live in each study's `figures/` |
| `docs/` | Roadmap, notes and paper drafts |

Reusable tools (for example divertor metrics or coil-limit checks that another tokamak could use) do **not** live here. They go to the [OFT fork](https://github.com/GitLostintheSauce/OpenFUSIONToolkit) on their own branch and then go upstream to OpenFUSIONToolkit as small PRs.

## Results so far

| Study | What it asks |
|---|---|
| [ARC V3A: triangularity and secondary X-point](studies/arc_v3a_triangularity_and_secondary_xpoint/README.md) | Can the coils change triangularity, or move a secondary X-point, while holding the other divertor targets? Includes a [report](studies/arc_v3a_triangularity_and_secondary_xpoint/REPORT.md) and an interactive [picture guide](studies/arc_v3a_triangularity_and_secondary_xpoint/picture_guide.html) |
| [ARC V3A: TokaMaker+TORAX flattop vs gEQDSK](studies/arc_v3a_torax_flattop_vs_geqdsk/README.md) | How close is the coupled TokaMaker+TORAX flattop to the public ARC V3A reference equilibrium? Includes a [report](studies/arc_v3a_torax_flattop_vs_geqdsk/REPORT.md). Needs TORAX 1.4.3 (study-only environment) |
| [SPARC: current, pressure and shape](studies/sparc_current_pressure_and_shape/README.md) | How does a SPARC equilibrium respond when plasma current, pressure, X-point position or elongation changes one at a time? Includes a [report](studies/sparc_current_pressure_and_shape/REPORT.md) and a [secondary-null follow-up](studies/sparc_current_pressure_and_shape/SECONDARY_NULLS_REPORT.md) |

## Data

ARC V3A inputs come from the public Zenodo record [19498373](https://zenodo.org/records/19498373) (Hillesheim et al., *J. Plasma Phys.* 2026, [doi:10.1017/S0022377826101706](https://doi.org/10.1017/S0022377826101706)). Scripts download them on demand, so they are not committed.

## Setup

Python 3.11+ is needed (TORAX 1.4.0). The environment was built with Python 3.12.13 on macOS arm64.

```bash
python3.12 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt          # top-level pins
# or, to reproduce exactly:
pip install -r requirements-lock.txt
```

| Package | Version |
|---|---|
| OpenFUSIONToolkit | 26.9 (PyPI wheel) |
| torax | 1.4.0 |
| jax / jaxlib | 0.10.0 |
| fusion_surrogates | 0.4.6 |
| OFT fork (examples, mesh) | commit `8fc4afe` |

Do not upgrade TORAX casually: the TokaMaker+TORAX workflow imports TORAX internals. OFT 26.9 needs TORAX 1.4.0 or later (`TokaMaker.pulse_design` imports `ConfinementMode`), even though OFT's `requirements_TokaMaker_TORAX.txt` still lists 1.3.1.

The lock file installs on Linux and macOS. CI (`.github/workflows/ci.yml`) installs it on Linux x86_64 (the default platform) and on macOS arm64, and runs the ARC V3A baseline, the [SPARC current, pressure and shape study](studies/sparc_current_pressure_and_shape/README.md) (checked against the committed results), the TokaMaker learning examples and the ITER TokaMaker+TORAX pulse example (checked against the [reference run](notebooks/iter_tokamaker_torax_reference/README.md)).

## Check the install

```bash
python scripts/run_arcv3a_baseline.py
```

This runs the unmodified OFT example `ARCV3A_baseline_ex.ipynb` from the fork in `work/` and writes [`figures/arcv3a_baseline_equilibrium.png`](figures/arcv3a_baseline_equilibrium.png). It should finish in under a minute and converge in about 24 iterations, giving a diverted double-null equilibrium with Ip = 12 MA and q95 ≈ 3.38.

To run the TokaMaker learning examples from the roadmap (fixed boundary, ITER baseline, DIII-D shape helper, ARC V3A vertical control):

```bash
python scripts/run_oft_examples.py
```

Executed notebooks and their plots go to `work/examples/` (gitignored).

## License

LGPL-3.0, matching OpenFUSIONToolkit. See [LICENSE](LICENSE).
