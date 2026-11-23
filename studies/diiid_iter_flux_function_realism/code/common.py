"""Shared paths, inputs and equilibrium metrics for the FF'/P' realism study.

Inputs come from the public OpenFUSIONToolkit repository at commit 8fc4afe (the commit this repo
reads OFT examples from). They are downloaded on demand into work/diiid_iter_flux_function_realism/inputs/
and checked against SHA-256 sums.
"""
import hashlib
import json
import os
import subprocess
import sys
import urllib.request

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
STUDY = os.path.dirname(HERE)
REPO = os.path.dirname(os.path.dirname(STUDY))
RESULTS = os.path.join(STUDY, 'results')
FIGURES = os.path.join(STUDY, 'figures')
SCRATCH = os.path.join(REPO, 'work', 'diiid_iter_flux_function_realism')
INPUTS = os.path.join(SCRATCH, 'inputs')
ARRAYS = os.path.join(SCRATCH, 'arrays')
LOGS = os.path.join(SCRATCH, 'logs')

OFT_COMMIT = '8fc4afea835160e1ef0f6d3e2426ab586996d2f8'
_RAW = f'https://raw.githubusercontent.com/OpenFUSIONToolkit/OpenFUSIONToolkit/{OFT_COMMIT}/'
DIIID_MESH = os.path.join(INPUTS, 'DIIID_mesh.h5')
DIIID_GEQDSK = os.path.join(INPUTS, 'g192185.02440')
ITER_MESH = os.path.join(INPUTS, 'ITER_mesh.h5')
# bootstrap.py at 8fc4afe includes OFT PR #316 (profile-derivative accuracy fix), which the
# pinned OpenFUSIONToolkit 26.9 wheel does not have.
BOOTSTRAP_8FC4AFE = os.path.join(INPUTS, 'bootstrap_8fc4afe.py')
SOURCES = {
    DIIID_MESH: ('src/examples/TokaMaker/DIIID/DIIID_mesh.h5',
                 '28b31f826075393b14053609782bf6bd1afb85000976b5b680ece17529dca55f'),
    DIIID_GEQDSK: ('src/examples/TokaMaker/DIIID/g192185.02440',
                   '6f33a01935847f25aea6edc45e939268ab910570ea6e8adb87528056a77520d5'),
    ITER_MESH: ('src/examples/TokaMaker/ITER/ITER_mesh.h5',
                '7f67f6b83770b6f92fb0afd94290870def0254163aaaa5f2e16b8ef79ac95c09'),
    BOOTSTRAP_8FC4AFE: ('src/python/OpenFUSIONToolkit/TokaMaker/bootstrap.py',
                        'f68d17fc4e4cde46b6c919d1bc947562da130be219538778a76a88a5c3a8f922'),
}
NTHREADS = 1      # 1 thread: bitwise-reproducible TokaMaker output
NPSI = 257        # profile sampling points, as in the ITER notebook
PSI_PAD = 1.E-3
EDGE = 0.85       # psi_hat above which P' counts as "edge" for the pedestal ratio
MU0 = 4.E-7*np.pi


def fetch_inputs():
    """Download any missing input and check every SHA-256."""
    os.makedirs(INPUTS, exist_ok=True)
    for path, (rel, sha) in SOURCES.items():
        if not os.path.exists(path):
            print(f'Downloading {rel}', flush=True)
            urllib.request.urlretrieve(_RAW + rel, path)
        with open(path, 'rb') as fh:
            got = hashlib.sha256(fh.read()).hexdigest()
        if got != sha:
            raise RuntimeError(f'{path}: SHA-256 {got} does not match expected {sha}')


def zero_crossings(x, y):
    """psi_hat locations where y changes sign (linear interpolation)."""
    s = np.sign(y)
    idx = np.where(s[:-1]*s[1:] < 0)[0]
    return [float(x[i] - y[i]*(x[i+1] - x[i])/(y[i+1] - y[i])) for i in idx]


def profile_metrics(psi, ffp, pp):
    """Shape metrics for FF' and P' sampled on psi_hat (0 = axis, 1 = LCFS)."""
    core, edge = psi <= EDGE, psi > EDGE
    return {
        'ffp_axis': float(ffp[0]), 'ffp_edge': float(ffp[-1]), 'ffp_min': float(ffp.min()),
        'ffp_zero_crossings': zero_crossings(psi, ffp),
        'pp_axis': float(pp[0]), 'pp_edge': float(pp[-1]),
        'pp_argmax_psi': float(psi[np.argmax(pp)]),
        'pp_edge_peak_over_core_peak': float(pp[edge].max()/pp[core].max()) if pp[core].max() > 0 else None,
        'pp_zero_crossings': zero_crossings(psi, pp),
    }


def equilibrium_record(mygs):
    """Scalars, profiles and shape metrics of the current TokaMaker equilibrium.

    Profiles are TokaMaker's solved output (get_profiles), in its convention: psi_hat = 0 on axis,
    P in Pa, P' and FF' per unit of TokaMaker flux (Wb/rad). The check 'p_axis_from_pp' integrates
    P' over the flux range and should equal P_ax, which confirms the units used for comparisons.
    """
    stats = mygs.get_stats()
    psi, f, fp, p, pp = mygs.get_profiles(npsi=NPSI, psi_pad=PSI_PAD)
    psi_q, q, _, _, _, _ = mygs.get_q(npsi=NPSI, psi_pad=PSI_PAD)
    psi_bounds = np.array(mygs.psi_bounds, dtype=float)    # [LCFS, axis] in TokaMaker flux
    dpsi = abs(psi_bounds[1] - psi_bounds[0])
    p_axis_from_pp = float(np.trapezoid(pp, psi)*dpsi)
    keep = ('Ip', 'kappa', 'delta', 'deltaU', 'deltaL', 'R_geo', 'a_geo', 'vol', 'q_0', 'q_95',
            'P_ax', 'W_MHD', 'beta_pol', 'l_i')
    rec = {k: float(stats[k]) for k in keep}
    # beta_t and beta_N from volume-averaged pressure and the vacuum field at R_geo
    B0 = mygs._F0/stats['R_geo']
    beta_t = 100.0*2.0*MU0*(stats['W_MHD']/1.5/stats['vol'])/B0**2
    rec.update({'B0_at_Rgeo': float(B0), 'beta_t_pct': float(beta_t),
                'beta_N': float(beta_t*stats['a_geo']*B0/(stats['Ip']/1.E6)),
                'psi_lcfs': float(psi_bounds[0]), 'psi_axis': float(psi_bounds[1]),
                'p_axis_from_pp': p_axis_from_pp, 'q_min': float(q.min()),
                'q_min_psi': float(psi_q[np.argmin(q)])})
    rec.update(profile_metrics(psi, f*fp, pp))
    arrays = {'psi': psi, 'F': f, 'FFp': f*fp, 'P': p, 'Pp': pp, 'psi_q': psi_q, 'q': q}
    return rec, arrays


def save_case(name, rec, arrays):
    os.makedirs(os.path.join(RESULTS, 'cases'), exist_ok=True)
    os.makedirs(ARRAYS, exist_ok=True)
    with open(os.path.join(RESULTS, 'cases', f'{name}.json'), 'w') as fh:
        json.dump(rec, fh, indent=1)
    if arrays is not None:
        np.savez(os.path.join(ARRAYS, f'{name}.npz'), **arrays)


def load_case(name):
    path = os.path.join(RESULTS, 'cases', f'{name}.json')
    if not os.path.exists(path):
        return None
    with open(path) as fh:
        rec = json.load(fh)
    arr = os.path.join(ARRAYS, f'{name}.npz')
    return rec, (dict(np.load(arr)) if os.path.exists(arr) else None)


def run_in_subprocesses(script, names, extra=()):
    """Run each case in its own process (only one TokaMaker instance is allowed per process)."""
    os.makedirs(LOGS, exist_ok=True)
    failed = []
    for name in names:
        log = os.path.join(LOGS, f'{os.path.splitext(os.path.basename(script))[0]}__{name}.log')
        print(f'  {name} ...', end=' ', flush=True)
        with open(log, 'w') as fh:
            proc = subprocess.run([sys.executable, script, '--case', name, *extra],
                                  stdout=fh, stderr=subprocess.STDOUT, cwd=REPO)
        print('ok' if proc.returncode == 0 else f'FAILED (exit {proc.returncode}, see {log})', flush=True)
        if proc.returncode != 0:
            failed.append(name)
    return failed
