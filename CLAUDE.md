# Working in this repo

Study: advanced divertor configurations (XPT, super-X, X-divertor, snowflake) on ARC V3A with TokaMaker and TokaMaker+TORAX. Plan: [docs/roadmap.md](docs/roadmap.md).

## Hard rules

1. **Never push without asking.** This applies to every repo and every remote, including this one and the OFT fork. Commit locally only when asked. Never push to upstream OpenFUSIONToolkit at all; changes reach it only through PRs from the fork.
2. **This repo is public. Only commit public data.** ARC V3A inputs must come from the public Zenodo record 19498373. Never commit data, drawings, coil/shield numbers or anything else from collaborators unless Eleanor has confirmed in the conversation that it is public. If in doubt, keep it out and ask. This applies to files, notebook outputs, figure annotations, configs and commit messages.
3. **No planning source documents in the repo.** Never add the original roadmap HTML, the project tracker PDF or similar planning files. `docs/roadmap.md` is a custom public summary; update it by hand, and keep it free of collaborator names and non-public details.
4. **Don't upgrade pinned versions** without asking. The TokaMaker+TORAX workflow imports TORAX internals, and TORAX main has already moved its divertor code.

## What goes where

| Kind of work | Location |
|---|---|
| ARC-specific scans, configs, figures, paper | This repo: `scripts/`, `configs/`, `figures/`, `notebooks/`, `docs/` |
| Reusable tools another tokamak could use (divertor metrics, coil limits, pulse example) | OFT fork at `~/Documents/code/OpenFUSIONToolkit`. One branch per tool, each a small PR to upstream OFT |
| TORAX | Used as is, pinned. No changes |

The OFT fork has two remotes: `origin` is `GitLostintheSauce/OpenFUSIONToolkit` (the fork) and `upstream` is `OpenFUSIONToolkit/OpenFUSIONToolkit`. Follow the fork's `AGENTS.md` for OFT coding conventions and PR titles (`TokaMaker: ...`).

## Environment

- Python venv at `.venv/` (Python 3.12.13). Run things with `.venv/bin/python` or activate it first.
- Pins are in `requirements.txt` (top level) and `requirements-lock.txt` (full `pip freeze`). Regenerate the lock file whenever the environment changes.
- `OpenFUSIONToolkit==26.9` (PyPI wheel), `torax==1.4.0`, `jax==0.10.0`, `jaxlib==0.10.0`, `fusion_surrogates==0.4.6`, `pydantic==2.12.5`, `pydantic_core==2.41.5`, `eqdsk==0.7.0`, `tqdm==4.67.3`, `xarray==2026.2.0`. OFT 26.9's `TokaMaker.pulse_design` needs TORAX >= 1.4.0 (it imports `ConfinementMode`), so these pins intentionally differ from `requirements_TokaMaker_TORAX.txt` in OFT, which still lists `torax==1.3.1`.
- The lock file must stay installable on Linux and macOS: mark platform-only packages with environment markers (e.g. `appnope ... ; sys_platform == "darwin"`).
- CI (`.github/workflows/ci.yml`) runs on Linux x86_64 by default and on macOS arm64. It installs the lock file and runs the ARC V3A baseline, the SPARC static scan, the TokaMaker learning examples and the ITER TokaMaker+TORAX pulse example. Use it to test environment changes.
- OFT examples and the ARC mesh are read from the fork at commit `8fc4afe`.
- Only one `OFT_env`/`TokaMaker` instance can exist per Python process. For scans, run each case in its own process.

## Conventions

- Every figure in `figures/` must be reproducible from a script or notebook in this repo.
- Scripts write scratch output to `work/` (gitignored). Downloaded inputs such as `geqdsk-*` are gitignored and fetched from Zenodo on demand.
- Sanity check: `python scripts/run_arcv3a_baseline.py` should converge and regenerate `figures/arcv3a_baseline_equilibrium.png`.
