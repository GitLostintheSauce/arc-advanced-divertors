"""ARC V3A shape-target sensitivity study (fixed profiles, static free-boundary TokaMaker).

Reproduces the equilibrium in F. Sheehan's ARC V3A shape-testing notebook
(`2026-06-08_ARC_diverted_shape_testing.ipynb`, shared with permission to publish; cells 1-11,
unchanged logic) as a script, then changes ONE input per case:

    A_baseline      delta = 0.65 (notebook as-is), nthreads = 1
    A_baseline_4thr delta = 0.65, nthreads = 4 (as the notebook ran it; reproducibility check)
    B_delta060      delta = 0.60
    C_secR505       secondary X-point targets moved from R = 5.00 to 5.05 m

The notebook is not modified. Its unused imports (omfit_classes, ipywidgets) and the
interactive strike-point picker are left out; they do not affect the equilibrium.

Usage (from the repo root):
    .venv/bin/python scripts/arc_shape_targets/shape_study.py            # all cases + analysis
    .venv/bin/python scripts/arc_shape_targets/shape_study.py --case B_delta060
    .venv/bin/python scripts/arc_shape_targets/shape_study.py --analyse  # figures/tables only
    (add --oft-root PATH if the OFT fork is not at ~/Documents/code/OpenFUSIONToolkit)

Each case runs in a fresh process (one TokaMaker per process, solver re-initialised).
The mesh is the public OFT ARC V3A mesh (src/examples/TokaMaker/ARC/arc_v3a_mesh.h5,
OFT PR #298), read from the OFT fork and checked by SHA-256.

Outputs: docs/arc_shape_targets/results/ (per-case JSON, tables), figures/arc_shape_targets/
(PNGs), work/arc_shape_targets/ (flux arrays, logs; gitignored).
"""
import argparse
import hashlib
import json
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
RES = os.path.join(REPO, 'docs', 'arc_shape_targets', 'results')
FIG = os.path.join(REPO, 'figures', 'arc_shape_targets')
ARR = os.path.join(REPO, 'work', 'arc_shape_targets', 'arrays')
LOGS = os.path.join(REPO, 'work', 'arc_shape_targets', 'logs')
OFT_ROOT = os.environ.get('OFT_ROOT', os.path.expanduser('~/Documents/code/OpenFUSIONToolkit'))
MESH = os.path.join(OFT_ROOT, 'src', 'examples', 'TokaMaker', 'ARC', 'arc_v3a_mesh.h5')
# Public OFT ARC V3A mesh (PR #298); identical to the copy in the shape-testing notebook
MESH_SHA256 = 'df0eada3ea0bc2503921c736705194d90ed83c3a0d1c25deddc26fcf36c99967'

# ---------------------------------------------------------------------------
# Baseline inputs, copied from the notebook
# ---------------------------------------------------------------------------
BASE = dict(
    R0=4.62, B0=11.4, Z0=0.0,
    Ip=12e6, pax=5e5,
    a=1.165, kappa=1.8, delta=0.65, zeta_outer=0.0, zeta_inner=0.0,
    primary_x_points=[[3.9, 2.2], [3.9, -2.2]],
    secondary_x_points=[[5.0, 3.8], [5.0, -3.8]],
    w_primary=100.0, w_secondary=1000.0,
    strike_point_targets=[[3.51, 2.36], [3.51, -2.36], [4.61, 3.87], [4.61, -3.87]],
    lcfs_x_points='all', x_point_isoflux_weight=10.0, secondary_isoflux_weight=3.0,
    strike_point_weight=1.0, trim_lcfs=True, trim_lcfs_perc_limit=0.85,
    coil_bound=100.E4,          # A per turn, every coil (notebook value)
    maxits=500, nthreads=1,
)

CASES = {
    'A_baseline':      {},
    'A_baseline_4thr': {'nthreads': 4},
    'B_delta060':      {'delta': 0.60},
    'C_secR505':       {'secondary_x_points': [[5.05, 3.8], [5.05, -3.8]]},
}
# Cases compared in the figures/tables (A_baseline_4thr is only a reproducibility check)
MAIN = ['A_baseline', 'B_delta060', 'C_secR505']
LABEL = {'A_baseline': 'A: baseline (δ=0.65)', 'B_delta060': 'B: δ=0.60',
         'C_secR505': 'C: secondary X R+5 cm', 'A_baseline_4thr': 'A (4 threads)'}


def build_min_norm_coil_reg(mygs):
    """Notebook cell 5: mild current regularization with U/L symmetry penalties."""
    reg_terms = []
    for name in mygs.coil_sets:
        if name.startswith('CS'):
            weight = 2.E-2 if name.startswith('CS1') else 1.E-2
            reg_terms.append(mygs.coil_reg_term({name: 1.0}, target=0.0, weight=weight))
        elif name.startswith('PF'):
            reg_terms.append(mygs.coil_reg_term({name: 1.0}, target=0.0, weight=1.E-2))
        elif name.startswith('VS'):
            reg_terms.append(mygs.coil_reg_term({name: 1.0}, target=0.0, weight=1.E-2))
    seen_pairs = set()
    for upper_name in mygs.coil_sets:
        if not upper_name.endswith('U'):
            continue
        lower_name = f"{upper_name[:-1]}L"
        if lower_name not in mygs.coil_sets:
            continue
        pair_key = tuple(sorted((upper_name, lower_name)))
        if pair_key in seen_pairs:
            continue
        reg_terms.append(mygs.coil_reg_term({upper_name: 1.0, lower_name: -1.0}, target=0.0, weight=1.E2))
        seen_pairs.add(pair_key)
    reg_terms.append(mygs.coil_reg_term({'#VSC': 1.0}, target=0.0, weight=1.E2))
    mygs.set_coil_reg(reg_terms=reg_terms)


def find_strike_points(mygs, tol=1.0e-12):
    """Notebook cell 9: intersections of the psi_hat = 1 contour with the limiter."""
    import matplotlib.pyplot as plt
    psi = mygs.get_psi(normalized=True)
    tmp_fig, tmp_ax = plt.subplots()
    cs = tmp_ax.tricontour(mygs.r[:, 0], mygs.r[:, 1], mygs.lc, psi, levels=[1.0])
    plt.close(tmp_fig)
    segs = []
    for path in cs.get_paths():
        v, codes = path.vertices, path.codes
        if codes is None:
            segs.append(v)
            continue
        starts = list(np.where(codes == 1)[0]) + [len(v)]
        segs += [v[i:j] for i, j in zip(starts[:-1], starts[1:]) if j - i > 1]
    lim = mygs.lim_contour
    hits = []
    for seg in segs:
        p, r = seg[:-1], seg[1:] - seg[:-1]
        for q, s in zip(lim[:-1], lim[1:] - lim[:-1]):
            denom = r[:, 0] * s[1] - r[:, 1] * s[0]
            ok = np.abs(denom) > tol
            if not ok.any():
                continue
            qp = q - p[ok]
            t = (qp[:, 0] * s[1] - qp[:, 1] * s[0]) / denom[ok]
            u = (qp[:, 0] * r[ok, 1] - qp[:, 1] * r[ok, 0]) / denom[ok]
            on = (t >= 0) & (t <= 1) & (u >= 0) & (u <= 1)
            if on.any():
                hits.append(p[ok][on] + t[on, None] * r[ok][on])
    if not hits:
        return np.zeros((0, 2)), segs
    hits = np.vstack(hits)
    keep = []
    for h in hits:
        if not any(np.hypot(*(h - k)) < 1.0e-3 for k in keep):
            keep.append(h)
    return np.array(keep), segs


def build_shape_targets(c):
    """Notebook cells 7-8: shape points, trimming, strike points and x-point isoflux points."""
    from OpenFUSIONToolkit.TokaMaker.util import create_isoflux_xpts, xpoints_from_moments
    prim = np.array(c['primary_x_points'])
    sec = np.array(c['secondary_x_points'])
    x_points = np.vstack([prim, sec])
    x_w = np.array([c['w_primary']]*len(prim) + [c['w_secondary']]*len(sec))
    k, d = c['kappa'], c['delta']
    x_moments = xpoints_from_moments(c['R0'], 0, c['a'], kappa_upper=k, delta_upper=d,
                                     kappa_lower=k, delta_lower=d)
    shape_all = create_isoflux_xpts(50, c['R0'], 0, c['a'], kappa_upper=k, delta_upper=d,
                                    kappa_lower=k, delta_lower=d,
                                    zeta_outer_upper=c['zeta_outer'], zeta_inner_upper=c['zeta_inner'],
                                    zeta_outer_lower=c['zeta_outer'], zeta_inner_lower=c['zeta_inner'])
    shape = shape_all
    if c['trim_lcfs']:
        z_lim = c['trim_lcfs_perc_limit']*np.max(np.abs(shape[:, 1]))
        shape = shape[np.abs(shape[:, 1]) <= z_lim]       # double null branch of the notebook
    strike = np.atleast_2d(np.array(c['strike_point_targets']))
    xpt_iso, xpt_w = [], []
    if c['lcfs_x_points'] in ('primary', 'all'):
        xpt_iso.append(prim); xpt_w.append(np.full(len(prim), c['x_point_isoflux_weight']))
    if c['lcfs_x_points'] in ('secondary', 'all'):
        xpt_iso.append(sec); xpt_w.append(np.full(len(sec), c['secondary_isoflux_weight']))
    xpt_iso = np.vstack(xpt_iso) if xpt_iso else np.zeros((0, 2))
    xpt_w = np.concatenate(xpt_w) if xpt_w else np.zeros(0)
    iso = np.vstack([shape, strike, xpt_iso])
    iso_w = np.concatenate([np.ones(len(shape)), np.full(len(strike), c['strike_point_weight']), xpt_w])
    groups = (['shape']*len(shape) + ['strike']*len(strike)
              + ['primary_xpt']*len(prim) + ['secondary_xpt']*len(sec))
    return dict(prim=prim, sec=sec, x_points=x_points, x_w=x_w, x_moments=x_moments,
                shape_all=shape_all, shape=shape, strike=strike, iso=iso, iso_w=iso_w, groups=groups)


def run_case(name):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.tri import Triangulation, LinearTriInterpolator
    from OpenFUSIONToolkit import OFT_env
    from OpenFUSIONToolkit.TokaMaker import TokaMaker
    from OpenFUSIONToolkit.TokaMaker.meshing import load_gs_mesh
    from OpenFUSIONToolkit.TokaMaker.util import create_power_flux_fun

    c = dict(BASE, **CASES[name])
    for d in (RES, ARR, FIG):
        os.makedirs(d, exist_ok=True)

    myOFT = OFT_env(nthreads=c['nthreads'])
    mygs = TokaMaker(myOFT)
    mesh_pts, mesh_lc, mesh_reg, coil_dict, cond_dict = load_gs_mesh(MESH)
    mygs.setup_mesh(mesh_pts, mesh_lc, mesh_reg)
    mygs.setup_regions(cond_dict=cond_dict, coil_dict=coil_dict)
    mygs.settings.maxits = c['maxits']
    mygs.setup(order=2, F0=c['R0']*c['B0'])
    mygs.set_coil_vsc({'PF5U': 1.0, 'PF5L': -1.0})

    mygs.set_coil_bounds({key: [-c['coil_bound'], c['coil_bound']] for key in mygs.coil_sets})
    mygs.set_profiles(ffp_prof=create_power_flux_fun(40, 1.5, 2.0),
                      pp_prof=create_power_flux_fun(40, 4.0, 1.0))

    t = build_shape_targets(c)
    mygs.set_saddle_constraints(t['x_points'], weights=t['x_w'])
    mygs.set_isoflux_constraints(t['iso'], weights=t['iso_w'])
    build_min_norm_coil_reg(mygs)
    mygs.set_targets(Ip=c['Ip'], pax=c['pax'])
    mygs.init_psi(c['R0'], 0, c['a'], c['kappa'], c['delta'])

    rec = {'case': name, 'change': CASES[name], 'nthreads': c['nthreads'],
           'delta_requested': c['delta'], 'kappa_requested': c['kappa'],
           'x_point_targets': t['x_points'].tolist(), 'x_point_weights': t['x_w'].tolist(),
           'x_points_from_moments': np.asarray(t['x_moments']).tolist(),
           'strike_targets': t['strike'].tolist(),
           'n_isoflux': {g: t['groups'].count(g) for g in dict.fromkeys(t['groups'])}}
    try:
        _, nits = mygs.solve(return_its=True)
        rec['converged'], rec['nl_its'] = True, int(nits)
    except ValueError as err:
        rec['converged'], rec['error'] = False, str(err)

    if rec['converged']:
        stats = mygs.get_stats()
        psi_n = mygs.get_psi(normalized=True)
        interp = LinearTriInterpolator(Triangulation(mygs.r[:, 0], mygs.r[:, 1], mygs.lc), psi_n)
        psi_at = lambda p: float(interp(p[0], p[1]))
        xpts, diverted = mygs.get_xpoints()
        xpts = np.zeros((0, 2)) if xpts is None else np.asarray(xpts)
        # Match each target to the nearest solved null
        matches = []
        for tgt, w in zip(t['x_points'], t['x_w']):
            if len(xpts) == 0:
                matches.append(None); continue
            i = int(np.argmin(np.hypot(*(xpts - tgt).T)))
            matches.append({'target': tgt.tolist(), 'weight': float(w), 'solved': xpts[i].tolist(),
                            'miss_mm': float(1e3*np.hypot(*(xpts[i] - tgt))),
                            'psi_hat': psi_at(xpts[i])})
        strike, segs = find_strike_points(mygs)
        # Isoflux residuals: psi_hat at every isoflux target (1.0 = constraint met exactly)
        iso_psi = np.array([psi_at(p) for p in t['iso']])
        iso_res = {}
        for g in dict.fromkeys(t['groups']):
            m = np.array(t['groups']) == g
            iso_res[g] = {'max_abs_dev': float(np.max(np.abs(iso_psi[m] - 1))),
                          'mean_abs_dev': float(np.mean(np.abs(iso_psi[m] - 1))),
                          'values': iso_psi[m].tolist()}
        cur, _ = mygs.get_coil_currents()
        cur = {k: float(v) for k, v in cur.items() if not k.startswith('#')}
        lcfs = mygs.trace_surf(0.9999)
        rec.update({
            'diverted': bool(diverted),
            'Ip': stats['Ip'], 'kappa': stats['kappa'], 'kappaU': stats['kappaU'], 'kappaL': stats['kappaL'],
            'delta': stats['delta'], 'deltaU': stats['deltaU'], 'deltaL': stats['deltaL'],
            'R_geo': stats['R_geo'], 'a_geo': stats['a_geo'], 'vol': stats['vol'],
            'q0': abs(stats['q_0']), 'q95': abs(stats['q_95']), 'l_i': stats['l_i'],
            'beta_p': stats['beta_pol']/100.0, 'beta_n': stats.get('beta_n'), 'P_ax': stats['P_ax'],
            'o_point': [float(v) for v in mygs.o_point],
            'all_nulls': [{'R': float(p[0]), 'Z': float(p[1]), 'psi_hat': psi_at(p)} for p in xpts],
            'x_point_matches': matches,
            'strike_points': strike.tolist(),
            'isoflux_residuals': iso_res,
            'coil_I_per_turn': cur,
            'coil_nturns': {k: int(coil_dict[k]['nturns']) for k in cur},
            'coil_Aturns': {k: cur[k]*coil_dict[k]['nturns'] for k in cur},
            'coil_bound_A_per_turn': c['coil_bound'],
            'coil_bound_fraction': {k: abs(cur[k])/c['coil_bound'] for k in cur},
        })
        np.savez(os.path.join(ARR, f'{name}.npz'), r=mygs.r, lc=mygs.lc, psi=psi_n,
                 lcfs=lcfs, lim=mygs.lim_contour, iso=t['iso'], shape_all=t['shape_all'],
                 sep=np.array(segs, dtype=object), allow_pickle=True)

    with open(os.path.join(RES, f'{name}.json'), 'w') as fh:
        json.dump(rec, fh, indent=1)

    # Per-case equilibrium figure (same layout as the notebook's plot)
    fig, ax = plt.subplots(figsize=(7, 8))
    mygs.plot_machine(fig, ax, coil_colormap='seismic', coil_symmap=True, coil_scale=1.E-6,
                      coil_clabel=r'$I_C$ [MA-turns]')
    if rec['converged']:
        mygs.plot_psi(fig, ax, plasma_nlevels=5, vacuum_nlevels=5, xpoint_color=None)
    mygs.plot_constraints(fig, ax, saddle_color=None)
    ax.plot([], [], '+', color='tab:red', ls='none', label='isoflux targets')
    ax.plot(t['x_points'][:, 0], t['x_points'][:, 1], 'x', color='tab:green', mew=2, ms=11, ls='none',
            label='X-point targets')
    if rec['converged']:
        if rec['all_nulls']:
            nn = np.array([[n['R'], n['Z']] for n in rec['all_nulls']])
            ax.plot(nn[:, 0], nn[:, 1], 's', mfc='none', mec='tab:green', mew=1.8, ms=9, ls='none',
                    label='nulls achieved')
        st = np.array(rec['strike_points'])
        if st.size:
            ax.plot(st[:, 0], st[:, 1], 'o', mfc='none', mec='tab:orange', mew=1.8, ms=10, ls='none',
                    label='strike points achieved')
    ax.plot(t['strike'][:, 0], t['strike'][:, 1], 'D', mfc='none', mec='tab:purple', mew=1.5, ms=7,
            ls='none', label='strike targets')
    ax.legend(fontsize=8, loc='upper right')
    ax.set_title(f"{LABEL[name]} — {'converged' if rec['converged'] else 'NOT CONVERGED'}")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, f'{name}_equilibrium.png'), dpi=130)
    plt.close(fig)
    print(f"[{name}] {'converged' if rec['converged'] else 'NOT converged'}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--case', choices=list(CASES))
    ap.add_argument('--analyse', action='store_true')
    ap.add_argument('--oft-root', help='OFT fork checkout (default ~/Documents/code/OpenFUSIONToolkit)')
    args = ap.parse_args()
    if args.oft_root:
        os.environ['OFT_ROOT'] = os.path.abspath(os.path.expanduser(args.oft_root))
        global MESH
        MESH = os.path.join(os.environ['OFT_ROOT'], 'src', 'examples', 'TokaMaker', 'ARC', 'arc_v3a_mesh.h5')
    if args.case:
        run_case(args.case)
        return
    if not args.analyse:
        with open(MESH, 'rb') as fh:
            if hashlib.sha256(fh.read()).hexdigest() != MESH_SHA256:
                raise RuntimeError(f'{MESH} does not match the expected public ARC V3A mesh')
        os.makedirs(LOGS, exist_ok=True)
        for name in CASES:
            with open(os.path.join(LOGS, f'{name}.log'), 'w') as fh:
                p = subprocess.run([sys.executable, __file__, '--case', name], stdout=fh, stderr=subprocess.STDOUT)
            print(f'{name}: exit {p.returncode}', flush=True)
    import analyse
    analyse.main()


if __name__ == '__main__':
    sys.path.insert(0, HERE)
    main()
