"""Run the ARC V3A TokaMaker+TORAX flattop and save its final equilibrium and profiles.

The run is F. Sheehan's public notebook `ARC_V3A/TMTX_project_parts_1_2_ARC_V3A_flattop_v1.ipynb`
from github.com/fs2819/TokaMaker_TORAX_project, pinned to commit bf5e8d1 and SHA-256 checked.
It is executed (Section 1, the single flattop run) with four small changes, all listed in
PATCHES below:

  1. the pip-install cell is replaced by a no-op, so nothing is installed into this repo's .venv;
  2. Hmode_profiles is imported from OpenFUSIONToolkit.TokaMaker.bootstrap, OFT's verbatim copy
     of the omfit_classes function, because omfit_classes is not in the pinned environment;
  3. the Z_eff sweep (Section 2) is switched off;
  4. one cell is added after Section 1 that saves the results (see EXTRACT_CELL).

The mesh is the public OFT ARC V3A mesh, copied from the OFT fork. It is byte-identical to the
ARC_V3A_mesh.h5 the notebook would download.

Cases (--case):
  baseline            the notebook as published: targets only for the two main X-points
  secondary_xpoints   the same, plus saddle (null) targets for secondary X-points at
                      SECONDARY_X_POINTS = (4.93, +/-3.64) m, where the public ARC V3A gEQDSK has its
                      outer-divertor nulls. Added to the seed equilibrium and to the coupled loop
                      (set_diverted_shape_targets(secondary_x_point_targets=...), weight = the main
                      X-point weight). Nothing else changes.

Writes to work/arc_v3a_torax_flattop_vs_geqdsk/<run dir>/ (gitignored; run/ for baseline,
run_secondary_xpoints/ for the second case):
    TMTX_..._executed.ipynb     the executed notebook
    extract/tmtx_flattop_final.geqdsk   final TokaMaker equilibrium, same grid as the Zenodo file
    extract/tmtx_flattop.json           TokaMaker stats, X-points, coil currents, TM+TX state,
                                        TORAX final scalars and profiles

Needs an environment installed from requirements-lock.txt (OFT 26.9 with TORAX 1.4.0): OFT 26.9's
TokaMaker+TORAX workflow does not import under TORAX 1.3.x. The script checks this first.

Usage (from the repo root; takes several minutes):
    PYTHON=.venv/bin/python   # any venv installed from requirements-lock.txt
    $PYTHON studies/arc_v3a_torax_flattop_vs_geqdsk/code/run_flattop.py [--case secondary_xpoints] [--oft-root PATH]
"""
import argparse
import hashlib
import os
import shutil
import urllib.request

import nbformat
from nbclient import NotebookClient

HERE = os.path.dirname(os.path.abspath(__file__))
STUDY = os.path.dirname(HERE)
REPO = os.path.dirname(os.path.dirname(STUDY))
SCRATCH = os.path.join(REPO, 'work', 'arc_v3a_torax_flattop_vs_geqdsk')
RUN_DIRS = {'baseline': os.path.join(SCRATCH, 'run'),
            'secondary_xpoints': os.path.join(SCRATCH, 'run_secondary_xpoints')}
# Outer-divertor nulls of the public ARC V3A gEQDSK (found by compare_to_geqdsk.py)
SECONDARY_X_POINTS = [[4.93, 3.64], [4.93, -3.64]]

FS_COMMIT = 'bf5e8d11a0d0146e4755ad53f37d74667fe0c177'
NB_NAME = 'TMTX_project_parts_1_2_ARC_V3A_flattop_v1.ipynb'
NB_URL = f'https://raw.githubusercontent.com/fs2819/TokaMaker_TORAX_project/{FS_COMMIT}/ARC_V3A/{NB_NAME}'
NB_SHA256 = '7c36ac9670a75016a450a3e40c13e3880172212a9cfb5f5dcc96bd27a6483763'
MESH_SHA256 = 'df0eada3ea0bc2503921c736705194d90ed83c3a0d1c25deddc26fcf36c99967'

# Extra patches for the secondary_xpoints case (same format as PATCHES).
SECONDARY_PATCHES = [
    ("mygs.set_saddle_constraints(X_POINTS, weights=np.full(X_POINTS.shape[0], SADDLE_WEIGHT))",
     "mygs.set_saddle_constraints(X_POINTS, weights=np.full(X_POINTS.shape[0], SADDLE_WEIGHT))",
     "SECONDARY_X_POINTS = np.array(" + repr(SECONDARY_X_POINTS) + ")  # [patched] outer-divertor nulls\n"
     "mygs.set_saddle_constraints(np.vstack([X_POINTS, SECONDARY_X_POINTS]),\n"
     "                            weights=np.full(X_POINTS.shape[0] + SECONDARY_X_POINTS.shape[0], SADDLE_WEIGHT))"),
    ("    x_point_weight=100,", "    x_point_weight=100,",
     "    x_point_weight=100,\n    secondary_x_point_targets=SECONDARY_X_POINTS,  # [patched]"),
]

# (marker that identifies the cell, old text, new text). old=None replaces the whole cell.
PATCHES = [
    ("OFT_VERSION = \"26.9\"", None,
     "# [patched] install cell skipped: the repo's pinned .venv already has OFT 26.9 and TORAX 1.4.0\n"
     "print('Install cell skipped (pinned environment).')"),
    ("from omfit_classes.utils_fusion import Hmode_profiles",
     "from omfit_classes.utils_fusion import Hmode_profiles",
     "from OpenFUSIONToolkit.TokaMaker.bootstrap import Hmode_profiles  # [patched] OFT's copy of the omfit function"),
    ("RUN_SWEEPS = True    # True ->", "RUN_SWEEPS = True    # True ->",
     "RUN_SWEEPS = False   # [patched] single flattop run only. True ->"),
]

EXTRACT_CELL = r'''
# [added] Save the final equilibrium, TM+TX state and TORAX results for the gEQDSK comparison.
import json as _json
_out = Path(os.environ['TMTX_EXTRACT_DIR'])
_out.mkdir(parents=True, exist_ok=True)

def _clean(x):
    """Make numpy/xarray values JSON-safe."""
    if isinstance(x, dict):
        return {str(k): _clean(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_clean(v) for v in x]
    if hasattr(x, 'tolist'):
        return _clean(x.tolist())
    if isinstance(x, float) and not np.isfinite(x):
        return None
    return x

_gs = tmtx._tm
# Same grid as the Zenodo gEQDSK (129 x 129, R 0.1-7.0 m, Z -5.5-5.5 m) for a like-for-like comparison.
# lcfs_pad=1e-4 puts the file's boundary at psi_n = 0.9999 (the surface TM+TX itself uses;
# the default 0.01 would label the 99 % surface as the boundary), and truncate_eq=False keeps
# the real flux outside the plasma instead of cutting it off there.
_gs.save_eqdsk(str(_out / 'tmtx_flattop_final.geqdsk'), nr=129, nz=129,
               rbounds=np.array([0.1, 7.0]), zbounds=np.array([-5.5, 5.5]), lcfs_pad=1.E-4, truncate_eq=False,
               run_info='TM+TX ARC V3A flattop final')
_stats = _gs.get_stats()
_xp, _div = _gs.get_xpoints()
_cur, _ = _gs.get_coil_currents()
_psi_q, _q, _, _, _, _ = _gs.get_q(npsi=101, psi_pad=0.001)
_psi_p, _f, _fp, _p, _pp = _gs.get_profiles(npsi=101)

_state = {}
for _k, _v in tmtx._state.items():
    if isinstance(_v, dict):
        if _v:
            _last = sorted(_v.keys())[-1]
            _state[_k] = {'key': _last, 'value': _clean(_v[_last])}
    elif hasattr(_v, 'shape') and getattr(_v, 'ndim', 0) == 1:
        _state[_k] = _clean(_v)

_dt = tmtx._data_tree
_tx_scalars = {}
for _name in _dt.scalars.data_vars:
    try:
        _val = _dt.scalars[_name].isel(time=-1).values
        if np.ndim(_val) == 0:
            _tx_scalars[_name] = float(_val)
    except Exception:
        pass
_tx_profiles = {}
for _name in ['n_e', 'n_i', 'T_e', 'T_i', 'q', 'magnetic_shear', 'pressure_thermal_total', 'pressure_total',
              'psi', 'j_total', 'j_ohmic', 'j_bootstrap', 'j_external', 'p_fusion_ion', 'Z_eff']:
    try:
        _da = _dt.profiles[_name].isel(time=-1)
        _tx_profiles[_name] = {'coord': _da.dims[0], 'x': _clean(_da.coords[_da.dims[0]].values),
                               'y': _clean(_da.values)}
    except Exception as _exc:
        _tx_profiles[_name] = {'error': str(_exc)}

_res = {
    'notebook_commit': os.environ.get('FS_COMMIT'),
    'case': os.environ.get('TMTX_CASE'),
    'tm_stats': _clean(_stats),
    'tm_diverted': bool(_div),
    'tm_x_points': _clean(_xp) if _xp is not None else [],
    'tm_o_point': _clean(_gs.o_point),
    'tm_coil_I_per_turn': _clean({k: v for k, v in _cur.items() if not k.startswith('#')}),
    'tm_q_profile': {'psi_n': _clean(_psi_q), 'q': _clean(_q)},
    'tm_profiles': {'psi_n': _clean(_psi_p), 'F': _clean(_f), 'Fp': _clean(_fp), 'p': _clean(_p), 'pp': _clean(_pp)},
    'tm_lcfs': _clean(_gs.trace_surf(0.9999)),
    'tmtx_state': _state,
    'tx_scalars': _tx_scalars,
    'tx_profiles': _tx_profiles,
}
with open(_out / 'tmtx_flattop.json', 'w') as _fh:
    # Objects that are not plain data (e.g. stored TokaMaker equilibria) are recorded by type name only.
    _json.dump(_res, _fh, indent=1, default=lambda _o: f'<{type(_o).__name__}>')
print('Saved results to', _out)
'''


def sha256(path):
    with open(path, 'rb') as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--oft-root', default=os.path.expanduser('~/Documents/code/OpenFUSIONToolkit'))
    ap.add_argument('--case', choices=list(RUN_DIRS), default='baseline')
    args = ap.parse_args()
    RUN_DIR = RUN_DIRS[args.case]
    EXTRACT_DIR = os.path.join(RUN_DIR, 'extract')
    import torax
    if not torax.__version__.startswith('1.4'):
        raise SystemExit(f'TORAX {torax.__version__} found; this run needs TORAX 1.4.x (install requirements-lock.txt)')
    os.makedirs(RUN_DIR, exist_ok=True)

    nb_path = os.path.join(RUN_DIR, NB_NAME)
    if not os.path.exists(nb_path):
        urllib.request.urlretrieve(NB_URL, nb_path)
    if sha256(nb_path) != NB_SHA256:
        raise RuntimeError(f'{nb_path} does not match the pinned notebook')
    mesh = os.path.join(RUN_DIR, 'ARC_V3A_mesh.h5')
    shutil.copy(os.path.join(args.oft_root, 'src', 'examples', 'TokaMaker', 'ARC', 'arc_v3a_mesh.h5'), mesh)
    if sha256(mesh) != MESH_SHA256:
        raise RuntimeError(f'{mesh} does not match the public ARC V3A mesh')

    nb = nbformat.read(nb_path, as_version=4)
    for marker, old, new in PATCHES + (SECONDARY_PATCHES if args.case == 'secondary_xpoints' else []):
        cells = [c for c in nb.cells if c.cell_type == 'code' and marker in c.source]
        if len(cells) != 1:
            raise RuntimeError(f'patch marker {marker!r} matched {len(cells)} cells')
        cells[0].source = new if old is None else cells[0].source.replace(old, new)
    # Insert the extraction cell right after Section 1.10 (before the Section 2 sweeps).
    i_sweeps = next(i for i, c in enumerate(nb.cells) if c.cell_type == 'markdown' and c.source.startswith('# 2 Sweeps'))
    nb.cells.insert(i_sweeps, nbformat.v4.new_code_cell(EXTRACT_CELL))

    os.environ['TMTX_EXTRACT_DIR'] = EXTRACT_DIR
    os.environ['FS_COMMIT'] = FS_COMMIT
    os.environ['TMTX_CASE'] = args.case
    client = NotebookClient(nb, timeout=-1, kernel_name='python3', resources={'metadata': {'path': RUN_DIR}})
    try:
        client.execute()
    finally:
        nbformat.write(nb, os.path.join(RUN_DIR, NB_NAME.replace('.ipynb', '_executed.ipynb')))
    print('Done:', EXTRACT_DIR)


if __name__ == '__main__':
    main()
