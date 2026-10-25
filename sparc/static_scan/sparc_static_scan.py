"""One-input-at-a-time scan of static free-boundary TokaMaker equilibria for SPARC.

Starting point: the upstream OFT example `src/examples/TokaMaker/SPARC/SPARC_public_baseline_ex.ipynb`
(OFT PR #317), which rebuilds the public SPARC Primary Reference Discharge (double null,
Ip = 8.7 MA, B0 = 12.2 T) from CFS's SPARCPublic repository. Every case uses the
notebook's setup unchanged except for ONE input.

Usage (from the repo root):
    .venv/bin/python sparc/static_scan/sparc_static_scan.py              # run all cases + summary
    .venv/bin/python sparc/static_scan/sparc_static_scan.py --case ip_p10  # one case only
    .venv/bin/python sparc/static_scan/sparc_static_scan.py --summary    # rebuild table/plots only

Each case runs in its own subprocess because only one TokaMaker instance may exist per
Python process. nthreads=1 so results are bit-reproducible.

Inputs are public and downloaded on first run into work/sparc_static/inputs/ (gitignored),
then checked against fixed SHA-256 sums:
    sparc_public_mesh.h5                 upstream OFT commit 85fb2f1, src/examples/TokaMaker/SPARC/
    SPARC_DN_PRD_freegs_20221013.geqdsk  SPARCPublic/PrimaryReferenceDischarge/2 - SPARC_DN_PRD_freegs_20221013
Tracked outputs (JSON, PNG, tables) go to sparc/static_scan/results/; bulky equilibrium
arrays and logs go to work/sparc_static/.
"""
import argparse
import hashlib
import json
import os
import subprocess
import sys
import urllib.request

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
OUT = os.path.join(HERE, 'results')                    # tracked: JSON, PNG, tables
SCRATCH = os.path.join(REPO, 'work', 'sparc_static')   # gitignored: inputs, arrays, logs
INPUTS = os.path.join(SCRATCH, 'inputs')
ARRAYS = os.path.join(SCRATCH, 'arrays')
MESH = os.path.join(INPUTS, 'sparc_public_mesh.h5')
GEQDSK = os.path.join(INPUTS, 'SPARC_DN_PRD_freegs_20221013.geqdsk')

# Public sources and their expected SHA-256 sums
SOURCES = {
    MESH: ('https://raw.githubusercontent.com/OpenFUSIONToolkit/OpenFUSIONToolkit/'
           '85fb2f11fb01de69a9d1b1ace011dd3365c666c0/src/examples/TokaMaker/SPARC/sparc_public_mesh.h5',
           'bacd38bd1c64ae8373ad1dd6be8bdb9d625267a12a91a8e1a2b0928f3a793083'),
    GEQDSK: ('https://raw.githubusercontent.com/cfs-energy/SPARCPublic/main/'
             'PrimaryReferenceDischarge/2%20-%20SPARC_DN_PRD_freegs_20221013',
             '6299e73a1ba1f2b8a1bb6d5e0dfa35f424ec09eebe153c96ea7779420eee85c6'),
}
NTHREADS = 1
R0_INIT = 1.85

# Coil current limits [A per turn], as {coil: (min, max)}. The public SPARC device
# description has NO coil limits (its supplies carry names only and its turn counts are
# toy data), so this is empty and the limit check reports "no limits available".
# Fill it only from a public source; it is used for checking, never as a solver bound.
COIL_LIMITS = {}

# ---------------------------------------------------------------------------
# Scan cases: each changes exactly one input relative to the baseline.
#   ip_scale    multiply the Ip target
#   pax_scale   multiply the on-axis pressure target (changes beta)
#   xpt_dR      shift both X-point (saddle) targets in R [m]
#   xpt_dZout   move both X-point targets away from the midplane by this much [m]
#   xpt_upper_dZ  shift only the upper X-point target in Z [m]
#   kappa_scale stretch the boundary (isoflux) targets and X-point targets in Z
# ---------------------------------------------------------------------------
CASES = {
    'baseline':      {},
    'ip_m20':        {'ip_scale': 0.80},
    'ip_m10':        {'ip_scale': 0.90},
    'ip_p10':        {'ip_scale': 1.10},
    'ip_p20':        {'ip_scale': 1.20},
    'pax_m20':       {'pax_scale': 0.80},
    'pax_p20':       {'pax_scale': 1.20},
    'xpt_Rin5':      {'xpt_dR': -0.05},
    'xpt_Rout5':     {'xpt_dR': +0.05},
    'xpt_Zin5':      {'xpt_dZout': -0.05},
    'xpt_Zout5':     {'xpt_dZout': +0.05},
    'xpt_upper_up5': {'xpt_upper_dZ': +0.05},
    'kappa_m5':      {'kappa_scale': 0.95},
    'kappa_p5':      {'kappa_scale': 1.05},
}


def fetch_inputs():
    """Download the public inputs if missing and check their checksums."""
    os.makedirs(INPUTS, exist_ok=True)
    for path, (url, sha) in SOURCES.items():
        if not os.path.exists(path):
            print(f'Downloading {os.path.basename(path)}', flush=True)
            urllib.request.urlretrieve(url, path)
        with open(path, 'rb') as fh:
            got = hashlib.sha256(fh.read()).hexdigest()
        if got != sha:
            raise RuntimeError(f'{path}: SHA-256 {got} does not match expected {sha}')


def run_case(name):
    """Build the baseline setup, apply one change, solve, and save outputs."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from OpenFUSIONToolkit import OFT_env
    from OpenFUSIONToolkit.TokaMaker import TokaMaker
    from OpenFUSIONToolkit.TokaMaker.meshing import load_gs_mesh
    from OpenFUSIONToolkit.TokaMaker.util import read_eqdsk

    mod = CASES[name]
    os.makedirs(OUT, exist_ok=True)
    eqdsk = read_eqdsk(GEQDSK)

    myOFT = OFT_env(nthreads=NTHREADS)
    mygs = TokaMaker(myOFT)
    mygs.settings.maxits = 200

    mesh_pts, mesh_lc, mesh_reg, coil_dict, cond_dict = load_gs_mesh(MESH)
    mygs.setup_mesh(mesh_pts, mesh_lc, mesh_reg)
    mygs.setup_regions(cond_dict=cond_dict, coil_dict=coil_dict)
    mygs.setup(order=2, F0=eqdsk['rcentr']*eqdsk['bcentr'])
    mygs.set_coil_vsc({'VS1U': 1.0, 'VS1L': -1.0})

    # Global targets (hard constraints)
    Ip_target = abs(eqdsk['ip'])*mod.get('ip_scale', 1.0)
    pax_target = eqdsk['pres'][0]*mod.get('pax_scale', 1.0)
    mygs.set_targets(Ip=Ip_target, pax=pax_target)

    # Shape targets (soft constraints): gEQDSK LCFS as isoflux, symmetric X-points as saddles
    isoflux_pts = eqdsk['rzout'].copy()
    x_up = isoflux_pts[np.argmax(isoflux_pts[:, 1])].copy()
    saddle_pts = np.array([x_up, [x_up[0], -x_up[1]]])
    if 'kappa_scale' in mod:
        isoflux_pts[:, 1] *= mod['kappa_scale']
        saddle_pts[:, 1] *= mod['kappa_scale']
    if 'xpt_dR' in mod:
        saddle_pts[:, 0] += mod['xpt_dR']
    if 'xpt_dZout' in mod:
        saddle_pts[:, 1] += np.sign(saddle_pts[:, 1])*mod['xpt_dZout']
    if 'xpt_upper_dZ' in mod:
        saddle_pts[0, 1] += mod['xpt_upper_dZ']
    mygs.set_isoflux_constraints(isoflux_pts)
    mygs.set_saddle_constraints(saddle_pts)

    # Coil regularisation, as in the notebook
    reg = []
    for cname in mygs.coil_sets:
        if cname.startswith('#'):
            continue
        weight = 1.E2 if cname.startswith('VS') else 1.E0
        reg.append(mygs.coil_reg_term({cname: 1.0}, target=0.0, weight=weight))
    reg.append(mygs.coil_reg_term({'#VSC': 1.0}, target=0.0, weight=1.E-2))
    mygs.set_coil_reg(reg_terms=reg)

    # Profiles from the gEQDSK (shape only; amplitudes are set by the Ip/pax targets)
    psi_eq = np.linspace(0.0, 1.0, np.size(eqdsk['ffprim']))
    psi_s = np.linspace(0.0, 1.0, 50)
    mygs.set_profiles(ffp_prof={'type': 'linterp', 'y': np.interp(psi_s, psi_eq, eqdsk['ffprim']), 'x': psi_s},
                      pp_prof={'type': 'linterp', 'y': np.interp(psi_s, psi_eq, eqdsk['pprime']), 'x': psi_s})

    # Initial guess, measured from the (possibly stretched) target boundary
    a = 0.5*(isoflux_pts[:, 0].max() - isoflux_pts[:, 0].min())
    kappa0 = (isoflux_pts[:, 1].max() - isoflux_pts[:, 1].min())/(2.0*a)
    delta0 = (R0_INIT - x_up[0])/a
    mygs.init_psi(R0_INIT, 0.0, a, kappa0, delta0)

    rec = {'case': name, 'mod': mod, 'nthreads': NTHREADS,
           'Ip_target': Ip_target, 'pax_target': pax_target,
           'saddle_targets': saddle_pts.tolist()}
    try:
        _, nits = mygs.solve(return_its=True)
        rec['converged'] = True
        rec['nl_its'] = int(nits)
    except ValueError as err:
        rec['converged'] = False
        rec['error'] = str(err)

    if rec['converged']:
        stats = mygs.get_stats()
        xpts, diverted = mygs.get_xpoints()
        currents, _ = mygs.get_coil_currents()
        rec.update({
            'diverted': bool(diverted),
            'Ip': stats['Ip'], 'q95': abs(stats['q_95']), 'q0': abs(stats['q_0']),
            'kappa': stats['kappa'], 'kappaU': stats['kappaU'], 'kappaL': stats['kappaL'],
            'delta': stats['delta'], 'deltaU': stats['deltaU'], 'deltaL': stats['deltaL'],
            # OFT returns beta_pol and beta_tor in percent; store beta_p as a fraction.
            # beta_n is already in the usual units (%.m.T/MA).
            'l_i': stats['l_i'], 'beta_p': stats['beta_pol']/100.0, 'beta_n': stats.get('beta_n'),
            'R_geo': stats['R_geo'], 'a_geo': stats['a_geo'], 'vol': stats['vol'],
            'P_ax': stats['P_ax'], 'W_MHD': stats['W_MHD'],
            'o_point': [float(v) for v in mygs.o_point],
            'x_points': [] if xpts is None else xpts.tolist(),
        })
        # Per-turn current and total ampere-turns per coil (VSC is virtual, skip it)
        rec['coil_I_per_turn'] = {k: float(v) for k, v in currents.items() if not k.startswith('#')}
        rec['coil_Aturns'] = {k: float(v)*coil_dict[k]['nturns'] for k, v in rec['coil_I_per_turn'].items()}
        # Coil-limit check (only possible if COIL_LIMITS is filled from a public source)
        if COIL_LIMITS:
            rec['coil_violations'] = [k for k, (lo, hi) in COIL_LIMITS.items()
                                      if not (lo <= rec['coil_I_per_turn'][k] <= hi)]
        else:
            rec['coil_violations'] = None
        lcfs = mygs.trace_surf(0.99 if diverted else 1.0)
        psi = mygs.get_psi()
        os.makedirs(ARRAYS, exist_ok=True)
        np.savez(os.path.join(ARRAYS, f'{name}_eq.npz'), lcfs=lcfs, psi=psi,
                 r=mygs.r, lc=mygs.lc)

    with open(os.path.join(OUT, f'{name}.json'), 'w') as fh:
        json.dump(rec, fh, indent=1)

    # Plot this case against the baseline boundary
    fig, ax = plt.subplots(1, 1, figsize=(6, 8))
    mygs.plot_machine(fig, ax, coil_colormap='seismic', coil_symmap=True, coil_scale=1.E-6,
                      coil_clabel=r'$I_C$ [MA-turns]')
    if rec['converged']:
        mygs.plot_psi(fig, ax, vacuum_nlevels=4)
    mygs.plot_constraints(fig, ax)
    base_npz = os.path.join(ARRAYS, 'baseline_eq.npz')
    if name != 'baseline' and os.path.exists(base_npz):
        bl = np.load(base_npz)['lcfs']
        ax.plot(bl[:, 0], bl[:, 1], 'k--', lw=1.5, label='baseline LCFS')
        with open(os.path.join(OUT, 'baseline.json')) as fh:
            bx = np.array(json.load(fh)['x_points'])
        if bx.size:
            ax.plot(bx[:, 0], bx[:, 1], 'ko', mfc='none', ms=8, label='baseline X-pts')
    if rec['converged'] and rec['x_points']:
        xp = np.array(rec['x_points'])
        ax.plot(xp[:, 0], xp[:, 1], 'mx', ms=10, mew=2.5, label='case X-pts')
    status = 'converged' if rec['converged'] else 'NOT CONVERGED'
    ax.set_title(f'{name}: {mod or "baseline"}\n{status}', fontsize=10)
    ax.legend(loc='lower right', fontsize=8)
    ax.set_xlabel('R (m)')
    ax.set_ylabel('Z (m)')
    fig.savefig(os.path.join(OUT, f'{name}.png'), dpi=130, bbox_inches='tight')
    plt.close(fig)
    print(f'[{name}] {status}')


def summarise():
    """Build a comparison table (CSV + Markdown) and a grid plot of all case boundaries."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    recs = {}
    for name in CASES:
        path = os.path.join(OUT, f'{name}.json')
        if os.path.exists(path):
            with open(path) as fh:
                recs[name] = json.load(fh)
    base = recs['baseline']

    def xpt_pair(rec):
        """Return the (upper, lower) divertor X-points: the null closest to each saddle target.

        TokaMaker can also find extra field nulls away from the plasma (e.g. near the
        outer divertor coils), so "highest/lowest null" is not the divertor X-point.
        """
        xp = np.array(rec.get('x_points') or [])
        if xp.size == 0:
            return None, None
        out = []
        for tgt in np.array(rec['saddle_targets']):
            d = np.hypot(*(xp - tgt).T)
            out.append(xp[np.argmin(d)] if d.min() < 0.15 else None)
        return out[0], out[1]

    def active_side(rec):
        """TokaMaker lists the active (LCFS-defining) X-point last."""
        xp = rec.get('x_points') or []
        if not xp:
            return '-'
        return 'upper' if xp[-1][1] > 0 else 'lower'

    cols = ['Ip', 'q95', 'kappa', 'delta', 'deltaU', 'deltaL', 'l_i', 'beta_p', 'beta_n']
    rows = []
    for name, rec in recs.items():
        row = {'case': name, 'change': json.dumps(rec['mod']) if rec['mod'] else '-',
               'converged': rec['converged'], 'its': rec.get('nl_its', '')}
        if rec['converged']:
            for c in cols:
                row[c] = rec[c]
                row['d_' + c] = rec[c] - base[c]
            up, lo = xpt_pair(rec)
            bup, blo = xpt_pair(base)
            row['n_null'] = len(rec['x_points'])
            row['active'] = active_side(rec)
            row['xU'] = None if up is None else tuple(up)
            row['xL'] = None if lo is None else tuple(lo)
            row['dxU_cm'] = None if (up is None or bup is None) else tuple(100*(up - bup))
            row['dxL_cm'] = None if (lo is None or blo is None) else tuple(100*(lo - blo))
            dI = {k: rec['coil_Aturns'][k] - base['coil_Aturns'][k] for k in rec['coil_Aturns']}
            kmax = max(dI, key=lambda k: abs(dI[k]))
            row['max_dcoil'] = f'{kmax} {dI[kmax]/1e6:+.2f} MA-t'
            row['coil_limits'] = ('no public limits' if rec['coil_violations'] is None
                                  else (', '.join(rec['coil_violations']) or 'ok'))
        else:
            row['error'] = rec.get('error', '')
        rows.append(row)

    # Markdown table
    def f(v, fmt):
        return '-' if v is None else fmt.format(v)

    def fx(p):
        return '-' if p is None else f'({p[0]:.3f}, {p[1]:+.3f})'

    def fdx(p):
        return '-' if p is None else f'({p[0]:+.1f}, {p[1]:+.1f})'

    lines = ['| case | change | conv (its) | Ip [MA] | q95 (Δ) | κ (Δ) | δ (Δ) | l_i (Δ) | β_p (Δ) | upper X-pt [m] (Δ cm) | lower X-pt [m] (Δ cm) | nulls / active | biggest coil change | coil limits |',
             '|---|---|---|---|---|---|---|---|---|---|---|---|---|---|']
    for r in rows:
        if not r['converged']:
            lines.append(f"| {r['case']} | {r['change']} | **NO** | | | | | | | | | | | {r['error'][:80]} |")
            continue
        lines.append(
            f"| {r['case']} | {r['change']} | yes ({r['its']}) | {r['Ip']/1e6:.2f} | "
            f"{r['q95']:.2f} ({r['d_q95']:+.2f}) | {r['kappa']:.3f} ({r['d_kappa']:+.3f}) | "
            f"{r['delta']:.3f} ({r['d_delta']:+.3f}) | {r['l_i']:.3f} ({r['d_l_i']:+.3f}) | "
            f"{r['beta_p']:.3f} ({r['d_beta_p']:+.3f}) | {fx(r['xU'])} {fdx(r['dxU_cm'])} | "
            f"{fx(r['xL'])} {fdx(r['dxL_cm'])} | {r['n_null']} / {r['active']} | {r['max_dcoil'] if r['case'] != 'baseline' else '-'} | {r['coil_limits']} |")
    md = '\n'.join(lines)
    with open(os.path.join(OUT, 'summary.md'), 'w') as fh:
        fh.write(md + '\n')

    # Coil current table (MA-turns), baseline and change per case
    coil_names = list(base['coil_Aturns'])
    clines = ['| coil | baseline [MA-t] | ' + ' | '.join(n for n in recs if n != 'baseline') + ' |',
              '|---|---|' + '---|'*(len(recs) - 1)]
    for c in coil_names:
        cells = []
        for n, rec in recs.items():
            if n == 'baseline':
                continue
            cells.append(f"{(rec['coil_Aturns'][c] - base['coil_Aturns'][c])/1e6:+.2f}" if rec['converged'] else 'n/c')
        clines.append(f"| {c} | {base['coil_Aturns'][c]/1e6:+.2f} | " + ' | '.join(cells) + ' |')
    with open(os.path.join(OUT, 'coil_currents.md'), 'w') as fh:
        fh.write('Coil currents in MA-turns: baseline, then change from baseline per case.\n\n')
        fh.write('\n'.join(clines) + '\n')

    # Grid figure: each case's LCFS against the baseline
    names = [n for n in recs if n != 'baseline']
    ncol = 5
    nrow = int(np.ceil(len(names)/ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.0*ncol, 4.6*nrow), sharex=True, sharey=True)
    bl = np.load(os.path.join(ARRAYS, 'baseline_eq.npz'))['lcfs']
    for ax, n in zip(axes.flat, names):
        rec = recs[n]
        ax.plot(bl[:, 0], bl[:, 1], 'k--', lw=1)
        bx = np.array(base['x_points'])
        ax.plot(bx[:, 0], bx[:, 1], 'ko', mfc='none', ms=5)
        if rec['converged']:
            lc = np.load(os.path.join(ARRAYS, f'{n}_eq.npz'))['lcfs']
            ax.plot(lc[:, 0], lc[:, 1], 'tab:red', lw=1.5)
            if rec['x_points']:
                xp = np.array(rec['x_points'])
                ax.plot(xp[:, 0], xp[:, 1], 'x', color='tab:red', ms=6, mew=2)
        else:
            ax.text(0.5, 0.5, 'NOT\nCONVERGED', transform=ax.transAxes, ha='center', color='r')
        st = np.array(rec['saddle_targets'])
        ax.plot(st[:, 0], st[:, 1], '+', color='tab:green', ms=8, mew=2)
        ax.set_title(n, fontsize=9)
        ax.set_aspect('equal')
    for ax in list(axes.flat)[len(names):]:
        ax.axis('off')
    fig.suptitle('LCFS per case (red) vs baseline (black dashed); green + = X-point targets', fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, 'all_cases_lcfs.png'), dpi=130)
    plt.close(fig)
    print(md)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--case', choices=list(CASES))
    parser.add_argument('--summary', action='store_true')
    args = parser.parse_args()
    if args.case:
        fetch_inputs()
        run_case(args.case)
        return
    if not args.summary:
        fetch_inputs()
        os.makedirs(os.path.join(SCRATCH, 'logs'), exist_ok=True)
        for name in CASES:   # baseline first, so other plots can overlay it
            log = os.path.join(SCRATCH, 'logs', f'{name}.log')
            with open(log, 'w') as fh:
                proc = subprocess.run([sys.executable, __file__, '--case', name],
                                      stdout=fh, stderr=subprocess.STDOUT)
            print(f'{name}: exit {proc.returncode} (log: {log})', flush=True)
    summarise()


if __name__ == '__main__':
    main()
