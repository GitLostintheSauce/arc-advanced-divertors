# Advanced divertors on ARC V3A

How hard is it to **form and hold** advanced divertor configurations (X-point target, super-X, X-divertor, snowflake) on ARC V3A with its real coils, shielding and pulse? This repo holds the study-specific work: scans, configurations, figures and the paper. It uses [TokaMaker](https://github.com/OpenFUSIONToolkit/OpenFUSIONToolkit) and the TokaMaker+TORAX workflow.

The plan, timeline and code-flow rules are in the [roadmap](docs/roadmap.md).

## Layout

| Folder | What goes here |
|---|---|
| `notebooks/` | Exploratory and figure-making notebooks for this study |
| `scripts/` | Runnable scripts: scans, batch runs, figure regeneration |
| `configs/` | Scan definitions, coil/shape targets and other inputs. Public ARC V3A data only (see below) |
| `figures/` | Generated plots. Each one should be reproducible from a script or notebook here |
| `docs/` | Roadmap, notes and paper drafts |

Reusable tools (for example divertor metrics or coil-limit checks that another tokamak could use) do **not** live here. They go to the [OFT fork](https://github.com/GitLostintheSauce/OpenFUSIONToolkit) on their own branch and then go upstream to OpenFUSIONToolkit as small PRs.

## Data

ARC V3A inputs come from the public Zenodo record [19498373](https://zenodo.org/records/19498373) (Hillesheim et al., *J. Plasma Phys.* 2026, [doi:10.1017/S0022377826101706](https://doi.org/10.1017/S0022377826101706)). Scripts download them on demand, so they are not committed.

## Setup

Python 3.11+ is needed (TORAX 1.3.1). The environment was built with Python 3.12.13 on macOS arm64.

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
| torax | 1.3.1 |
| jax / jaxlib | 0.9.2 |
| OFT fork (examples, mesh) | commit `8fc4afe` |

Do not upgrade TORAX casually: the TokaMaker+TORAX workflow imports TORAX internals.

## Check the install

```bash
python scripts/run_arcv3a_baseline.py
```

This runs the unmodified OFT example `ARCV3A_baseline_ex.ipynb` from the fork in `work/` and writes [`figures/arcv3a_baseline_equilibrium.png`](figures/arcv3a_baseline_equilibrium.png). It should finish in under a minute and converge in about 24 iterations, giving a diverted double-null equilibrium with Ip = 12 MA and q95 ≈ 3.38.

## License

LGPL-3.0, matching OpenFUSIONToolkit. See [LICENSE](LICENSE).
