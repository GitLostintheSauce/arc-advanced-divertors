"""Compare the TokaMaker+TORAX ARC V3A flattop runs with the public ARC V3A gEQDSK.

Cases: the public gEQDSK ('ref'), the TM+TX flattop as published ('tmtx') and, if it has been run,
the TM+TX flattop with secondary X-point targets ('tmtx_sec'; run_flattop.py --case
secondary_xpoints). All equilibria are read as gEQDSK files (the TM+TX one is written by run_flattop.py on the same
129 x 129 grid as the Zenodo file) and measured by the SAME functions below, so differences come
from the equilibria, not from the measuring:

  * shape (R_geo, a, kappa, delta upper/lower) from the psi_n = 0.999 contour enclosing the axis;
  * X-points: minima of |grad psi| near the top and bottom of that contour, and psi_n there;
  * profiles q, p, FF', p' against normalized poloidal flux, straight from the files;
  * volume integrals over the plasma: <p>, beta_t, beta_N, beta_p, l_i (standard normalization,
    as in TokaMaker's get_stats), and Ip from the line integral of B_p around the boundary.

The TM+TX values computed here are also checked against TokaMaker's own get_stats() numbers saved
by run_flattop.py, to show the measuring functions are right.

Writes results/comparison.json, results/comparison_table.md and figures/*.png.

Usage (from the repo root, after run_flattop.py):
    .venv/bin/python studies/arc_v3a_torax_flattop_vs_geqdsk/code/compare_to_geqdsk.py
"""
import hashlib
import json
import os
import urllib.request

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.path import Path
from scipy.interpolate import RectBivariateSpline
from scipy.optimize import minimize

HERE = os.path.dirname(os.path.abspath(__file__))
STUDY = os.path.dirname(HERE)
REPO = os.path.dirname(os.path.dirname(STUDY))
RES = os.path.join(STUDY, 'results')
FIG = os.path.join(STUDY, 'figures')
SCRATCH = os.path.join(REPO, 'work', 'arc_v3a_torax_flattop_vs_geqdsk')
RUNS = {'tmtx': os.path.join(SCRATCH, 'run', 'extract'),
        'tmtx_sec': os.path.join(SCRATCH, 'run_secondary_xpoints', 'extract')}
REF_EQDSK = os.path.join(SCRATCH, 'inputs', 'geqdsk-ARCv3a')
REF_URL = 'https://zenodo.org/records/19498373/files/geqdsk-ARCv3a'
REF_SHA256 = None   # filled on first download; see fetch_reference()

MU0 = 4e-7*np.pi
R0_DEVICE = 4.62
INK, INK2, GRID = '#0b0b0b', '#52514e', '#e4e3df'
COLOR = {'ref': '#2a78d6', 'tmtx': '#eb6834', 'tmtx_sec': '#1baf7a'}
STYLE = {'ref': '-', 'tmtx': '--', 'tmtx_sec': ':'}
LABEL = {'ref': 'gEQDSK (Zenodo, FreeGS)', 'tmtx': 'TM+TX flattop (as published)',
         'tmtx_sec': 'TM+TX flattop + secondary X-point targets'}
SHORT = {'ref': 'gEQDSK', 'tmtx': 'TM+TX', 'tmtx_sec': 'TM+TX + secondary X-pts'}
COIL_LIMIT_MAT = 50.0    # MA-turns: the coil bounds set in the TM+TX notebook (+/- 50 MA-turns per coil)


def fetch_reference():
    """Download the public ARC V3A gEQDSK from Zenodo record 19498373 if needed."""
    os.makedirs(os.path.dirname(REF_EQDSK), exist_ok=True)
    if not os.path.exists(REF_EQDSK):
        urllib.request.urlretrieve(REF_URL, REF_EQDSK)
    with open(REF_EQDSK, 'rb') as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def load(path):
    from OpenFUSIONToolkit.TokaMaker.util import read_eqdsk
    e = read_eqdsk(path)
    R = np.linspace(e['rleft'], e['rleft'] + e['rdim'], e['nr'])
    Z = np.linspace(e['zmid'] - e['zdim']/2, e['zmid'] + e['zdim']/2, e['nz'])
    raw = np.asarray(e['psirz'], dtype=float)
    # Pick the array orientation [Z, R] for which psi_n is ~0 at the stated magnetic axis
    # (both candidates are tried because the grid is square).
    best = None
    for psi in (raw, raw.T):
        pn = (psi - e['psimag'])/(e['psibry'] - e['psimag'])
        at_axis = abs(RectBivariateSpline(Z, R, pn)(e['zaxis'], e['raxis'])[0, 0])
        if best is None or at_axis < best[0]:
            best = (at_axis, psi, pn)
    at_axis, psi, pn = best
    assert at_axis < 0.02, f'psi_n at the magnetic axis is {at_axis:.3f}, expected ~0'
    e.update(R=R, Z=Z, psi=psi, pn=pn, psi_1d=np.linspace(0, 1, e['nr']))
    return e


def lcfs_contour(e, level=0.999):
    """Closed psi_n contour that encloses the magnetic axis."""
    fig, ax = plt.subplots()
    cs = ax.contour(e['R'], e['Z'], e['pn'], levels=[level])
    plt.close(fig)
    best = None
    for path in cs.get_paths():
        codes = path.codes
        starts = [0] if codes is None else list(np.where(codes == 1)[0])
        starts.append(len(path.vertices))
        for a, b in zip(starts[:-1], starts[1:]):
            seg = path.vertices[a:b]
            if len(seg) > 20 and Path(seg).contains_point((e['raxis'], e['zaxis'])):
                if best is None or len(seg) > len(best):
                    best = seg
    return best


def separatrix_pieces(e, level=1.0):
    fig, ax = plt.subplots()
    cs = ax.contour(e['R'], e['Z'], e['pn'], levels=[level])
    plt.close(fig)
    out = []
    for path in cs.get_paths():
        codes = path.codes
        starts = [0] if codes is None else list(np.where(codes == 1)[0])
        starts.append(len(path.vertices))
        out += [path.vertices[a:b] for a, b in zip(starts[:-1], starts[1:]) if b - a > 2]
    return out


def shape(lcfs):
    rmin, rmax = lcfs[:, 0].min(), lcfs[:, 0].max()
    R_geo, a = 0.5*(rmin + rmax), 0.5*(rmax - rmin)
    up, lo = lcfs[np.argmax(lcfs[:, 1])], lcfs[np.argmin(lcfs[:, 1])]
    return {'R_geo': R_geo, 'a': a, 'kappa': (up[1] - lo[1])/(2*a),
            'deltaU': (R_geo - up[0])/a, 'deltaL': (R_geo - lo[0])/a,
            'delta': 0.5*((R_geo - up[0]) + (R_geo - lo[0]))/a, 'top': up.tolist(), 'bottom': lo.tolist()}


def x_points(e, lcfs):
    """Minimize |grad psi|^2 starting from the top and bottom of the boundary."""
    spl = RectBivariateSpline(e['Z'], e['R'], e['psi'])
    spn = RectBivariateSpline(e['Z'], e['R'], e['pn'])

    def g2(x):
        R, Z = x
        return spl(Z, R, dy=1)[0, 0]**2 + spl(Z, R, dx=1)[0, 0]**2

    out = []
    for start in (lcfs[np.argmax(lcfs[:, 1])], lcfs[np.argmin(lcfs[:, 1])]):
        res = minimize(g2, start, method='Nelder-Mead', options={'xatol': 1e-5, 'fatol': 1e-12})
        R, Z = res.x
        out.append({'R': float(R), 'Z': float(Z), 'psi_n': float(spn(Z, R)[0, 0]),
                    'Bp_T': float(np.sqrt(g2(res.x))/R)})
    return out


def secondary_nulls(e, starts=((5.0, 3.7), (5.0, -3.7)), null_tol=5e-3):
    """Look for a field null in the outer divertor, starting near (5.0, +/-3.7) m.

    Reports where |grad psi| is smallest nearby, the poloidal field there and psi_n. A point counts
    as a null if |B_p| < null_tol [T] and it stays within 0.4 m of the start.
    """
    spl = RectBivariateSpline(e['Z'], e['R'], e['psi'])
    spn = RectBivariateSpline(e['Z'], e['R'], e['pn'])

    def g2(x):
        R, Z = x
        return spl(Z, R, dy=1)[0, 0]**2 + spl(Z, R, dx=1)[0, 0]**2

    out = []
    for start in starts:
        res = minimize(g2, np.array(start), method='Nelder-Mead', options={'xatol': 1e-5, 'fatol': 1e-14})
        R, Z = res.x
        bp = float(np.sqrt(g2(res.x))/R)
        near = float(np.hypot(R - start[0], Z - start[1]))
        out.append({'R': float(R), 'Z': float(Z), 'Bp_T': bp, 'psi_n': float(spn(Z, R)[0, 0]),
                    'is_null': bool(bp < null_tol and near < 0.4)})
    return out


def integrals(e, lcfs, F0):
    """Volume integrals over the plasma (cells whose centre lies inside the boundary)."""
    R, Z = e['R'], e['Z']
    dR, dZ = R[1] - R[0], Z[1] - Z[0]
    RR, ZZ = np.meshgrid(R, Z)
    inside = Path(lcfs).contains_points(np.column_stack([RR.ravel(), ZZ.ravel()])).reshape(RR.shape)
    dV = 2*np.pi*RR*dR*dZ
    p2d = np.interp(np.clip(e['pn'], 0, 1), e['psi_1d'], np.abs(e['pres']))
    dpsi_dZ, dpsi_dR = np.gradient(e['psi'], dZ, dR)
    Bp2 = (dpsi_dR**2 + dpsi_dZ**2)/RR**2
    V = np.sum(dV[inside])
    p_avg = np.sum((p2d*dV)[inside])/V
    # Ip from Ampere's law around the boundary; boundary length for the beta_p / l_i normalization
    spl = RectBivariateSpline(Z, R, e['psi'])
    seg = np.diff(np.vstack([lcfs, lcfs[:1]]), axis=0)
    mid = 0.5*(np.vstack([lcfs, lcfs[:1]])[1:] + lcfs)
    L = float(np.sum(np.hypot(*seg.T)))
    bp_mid = np.hypot(spl(mid[:, 1], mid[:, 0], dx=1, grid=False), spl(mid[:, 1], mid[:, 0], dy=1, grid=False))/mid[:, 0]
    Ip_loop = float(np.sum(bp_mid*np.hypot(*seg.T))/MU0)
    Ip = abs(e['ip'])
    sh = shape(lcfs)
    B0 = abs(F0)/sh['R_geo']
    beta_t = 2*MU0*p_avg/B0**2
    Bp_edge = MU0*Ip/L
    return {'volume_m3': float(V), 'p_vol_avg_Pa': float(p_avg),
            'beta_t_pct': float(100*beta_t), 'beta_N': float(100*beta_t*sh['a']*B0/(Ip/1e6)),
            'beta_p': float(2*MU0*p_avg/Bp_edge**2), 'l_i': float((np.sum((Bp2*dV)[inside])/V)/Bp_edge**2),
            'Ip_from_boundary_MA': Ip_loop/1e6, 'boundary_length_m': L, 'B0_at_Rgeo_T': float(B0)}


def measure(e):
    lc = lcfs_contour(e)
    F0 = e['rcentr']*e['bcentr']
    m = {'Ip_MA': abs(e['ip'])/1e6, 'F0_Tm': abs(F0), 'B0_at_4.62m_T': abs(F0)/R0_DEVICE,
         'R_axis': float(e['raxis']), 'Z_axis': float(e['zaxis']),
         'q0': float(abs(e['qpsi'][0])), 'q95': float(abs(np.interp(0.95, e['psi_1d'], e['qpsi']))),
         'q_edge_0.995': float(abs(np.interp(0.995, e['psi_1d'], e['qpsi']))),
         'p_axis_MPa': float(abs(e['pres'][0])/1e6), 'psi_axis_to_boundary_Wb_per_rad': float(abs(e['psibry'] - e['psimag']))}
    m.update(shape(lc))
    m['x_points'] = x_points(e, lc)
    m['delta_from_xpoints'] = float((m['R_geo'] - np.mean([x['R'] for x in m['x_points']]))/m['a'])
    m['secondary_nulls'] = secondary_nulls(e)
    m.update(integrals(e, lc, F0))
    return m, lc


def fmt(v, f='{:.3f}'):
    return '—' if v is None else f.format(v)


def coil_turns():
    """Turn count per coil from the public ARC V3A mesh the TM+TX run used."""
    from OpenFUSIONToolkit.TokaMaker.meshing import load_gs_mesh
    mesh = os.path.join(SCRATCH, 'run', 'ARC_V3A_mesh.h5')
    _, _, _, coil_dict, _ = load_gs_mesh(mesh)
    return {k: v['nturns'] for k, v in coil_dict.items()}


def main():
    os.makedirs(RES, exist_ok=True)
    os.makedirs(FIG, exist_ok=True)
    ref_sha = fetch_reference()
    E, J, M, LC = {}, {}, {}, {}
    E['ref'] = load(REF_EQDSK)
    for key, d in RUNS.items():
        if os.path.exists(os.path.join(d, 'tmtx_flattop.json')):
            E[key] = load(os.path.join(d, 'tmtx_flattop_final.geqdsk'))
            with open(os.path.join(d, 'tmtx_flattop.json')) as fh:
                J[key] = json.load(fh)
    cases = [k for k in ('ref', 'tmtx', 'tmtx_sec') if k in E]
    for k in cases:
        M[k], LC[k] = measure(E[k])

    # Check the measuring functions against TokaMaker's own numbers, for each TM+TX run
    checks = {}
    for k in cases[1:]:
        st = J[k]['tm_stats']
        checks[k] = {
            'kappa': (M[k]['kappa'], st['kappa']), 'q95': (M[k]['q95'], abs(st['q_95'])),
            'l_i': (M[k]['l_i'], st['l_i']), 'beta_p': (M[k]['beta_p'], st['beta_pol']/100),
            'beta_N': (M[k]['beta_N'], st.get('beta_n')), 'volume_m3': (M[k]['volume_m3'], st['vol']),
            'Ip_MA (boundary integral vs header)': (M[k]['Ip_from_boundary_MA'], M[k]['Ip_MA'])}

    # TM+TX run context and coil currents
    turns = coil_turns()
    context, coils = {}, {}
    for k in cases[1:]:
        tx, state = J[k]['tx_scalars'], J[k]['tmtx_state']
        last = lambda name: (state.get(name) or [None])[-1] if isinstance(state.get(name), list) else None
        context[k] = {'Q_fusion': tx.get('Q_fusion'),
                      'P_fusion_MW': tx['P_fusion']/1e6 if tx.get('P_fusion') is not None else None,
                      'H98': tx.get('H98'), 'beta_N_torax': tx.get('beta_N'), 'q95_torax': tx.get('q95'),
                      'li3_torax': tx.get('li3'),
                      'l_i_tm (last time)': last('l_i_tm'), 'l_i_tx (last time)': last('l_i_tx'),
                      'beta_N_tm (last time)': last('beta_N_tm'), 'beta_N_tx (last time)': last('beta_N_tx'),
                      'Ip_tm_MA (last time)': (last('Ip_tm') or 0)/1e6, 'Ip_tx_MA (last time)': (last('Ip_tx') or 0)/1e6}
        coils[k] = {c: J[k]['tm_coil_I_per_turn'][c]*turns[c]/1e6 for c in J[k]['tm_coil_I_per_turn']}

    out = {'reference_geqdsk': {'source': REF_URL, 'sha256': ref_sha, 'header': str(E['ref'].get('case', ''))},
           'tmtx_notebook_commit': J['tmtx'].get('notebook_commit') if 'tmtx' in J else None,
           'cases': {k: LABEL[k] for k in cases}, 'measured': M,
           'measuring_check_vs_tokamaker_get_stats': checks, 'tmtx_context': context,
           'coil_currents_MAturns': coils, 'coil_limit_MAturns': COIL_LIMIT_MAT}
    with open(os.path.join(RES, 'comparison.json'), 'w') as fh:
        json.dump(out, fh, indent=1, default=float)

    # ---------------- table ----------------
    rows = [('Plasma current Ip [MA]', 'Ip_MA', '{:.2f}'), ('Vacuum field at R = 4.62 m [T]', 'B0_at_4.62m_T', '{:.2f}'),
            ('Magnetic axis R [m]', 'R_axis', '{:.3f}'), ('Geometric centre R_geo [m]', 'R_geo', '{:.3f}'),
            ('Minor radius a [m]', 'a', '{:.3f}'), ('Elongation κ', 'kappa', '{:.3f}'),
            ('Triangularity δ (mean, from boundary)', 'delta', '{:.3f}'),
            ('Triangularity δ (from X-points: (R_geo − R_x)/a)', 'delta_from_xpoints', '{:.3f}'),
            ('Plasma volume [m³]', 'volume_m3', '{:.1f}'), ('q0', 'q0', '{:.2f}'), ('q95', 'q95', '{:.2f}'),
            ('Axis pressure [MPa]', 'p_axis_MPa', '{:.2f}'), ('Volume-average pressure [MPa]', 'p_vol_avg_MPa', '{:.3f}'),
            ('β_N', 'beta_N', '{:.2f}'), ('β_p', 'beta_p', '{:.3f}'), ('l_i', 'l_i', '{:.3f}'),
            ('Poloidal flux axis→boundary [Wb/rad]', 'psi_axis_to_boundary_Wb_per_rad', '{:.2f}')]
    for k in cases:
        M[k]['p_vol_avg_MPa'] = M[k]['p_vol_avg_Pa']/1e6
    head = '| Quantity | ' + ' | '.join(SHORT[k] for k in cases) + ' |'
    L = ['Columns: ' + '; '.join(f'**{SHORT[k]}** = {LABEL[k]}' for k in cases) + '. '
         'Percentages are differences from the gEQDSK.', '', head, '|---|' + '---|'*len(cases)]
    for name, key, f in rows:
        cells = []
        for k in cases:
            v = M[k][key]
            if k == 'ref' or not M['ref'][key]:
                cells.append(f.format(v))
            else:
                cells.append(f"{f.format(v)} ({100*(v/M['ref'][key]-1):+.0f} %)")
        L.append(f'| {name} | ' + ' | '.join(cells) + ' |')
    L.append('| δ upper / lower (from boundary) | ' + ' | '.join(
        f"{M[k]['deltaU']:.3f} / {M[k]['deltaL']:.3f}" for k in cases) + ' |')
    for i, nm in enumerate(('Upper', 'Lower')):
        L.append(f'| {nm} main X-point (R, Z) [m], ψ̂ there | ' + ' | '.join(
            f"({M[k]['x_points'][i]['R']:.3f}, {M[k]['x_points'][i]['Z']:+.3f}), {M[k]['x_points'][i]['psi_n']:.3f}"
            for k in cases) + ' |')

    def sn(x):
        if x['is_null']:
            return f"({x['R']:.2f}, {x['Z']:+.2f}), ψ̂ = {x['psi_n']:.3f}"
        return f"none (min B_p {1e3*x['Bp_T']:.0f} mT)"
    for i, nm in enumerate(('Upper', 'Lower')):
        L.append(f'| {nm} outer-divertor null near (5.0, ±3.7) m | ' + ' | '.join(
            sn(M[k]['secondary_nulls'][i]) for k in cases) + ' |')

    L += ['', "Check of the measuring code against TokaMaker's own `get_stats()` for each TM+TX run:", '',
          '| Quantity | ' + ' | '.join(f'{SHORT[k]}: here / TokaMaker' for k in cases[1:]) + ' |',
          '|---|' + '---|'*len(cases[1:])]
    for q in (checks[cases[1]] if len(cases) > 1 else {}):
        L.append(f'| {q} | ' + ' | '.join(
            f"{fmt(checks[k][q][0])} / {fmt(checks[k][q][1])}" for k in cases[1:]) + ' |')
    L += ['', 'TM+TX run results (end of the coupled run):', '',
          '| Quantity | ' + ' | '.join(SHORT[k] for k in cases[1:]) + ' |', '|---|' + '---|'*len(cases[1:])]
    for q in (context[cases[1]] if len(cases) > 1 else {}):
        L.append(f'| {q} | ' + ' | '.join(fmt(context[k][q], '{:.4g}') for k in cases[1:]) + ' |')
    if len(cases) > 1:
        L += ['', f'Coil currents [MA-turns]. The TM+TX notebook bounds every coil at ±{COIL_LIMIT_MAT:.0f} MA-turns.', '',
              '| Coil | ' + ' | '.join(SHORT[k] for k in cases[1:]) + (' | Change |' if 'tmtx_sec' in coils else ' |'),
              '|---|' + '---|'*len(cases[1:]) + ('---|' if 'tmtx_sec' in coils else '')]
        for c in coils[cases[1]]:
            vals = [coils[k][c] for k in cases[1:]]
            row = f'| {c} | ' + ' | '.join(f'{v:+.2f}' + (' **at limit**' if abs(v) > 0.99*COIL_LIMIT_MAT else '')
                                         for v in vals)
            if 'tmtx_sec' in coils:
                row += f" | {coils['tmtx_sec'][c] - coils['tmtx'][c]:+.2f}"
            L.append(row + ' |')
        L.append('| Σ\\|I\\| | ' + ' | '.join(f'{sum(abs(v) for v in coils[k].values()):.1f}' for k in cases[1:]) +
                 (' | |' if 'tmtx_sec' in coils else ' |'))
    with open(os.path.join(RES, 'comparison_table.md'), 'w') as fh:
        fh.write('\n'.join(L) + '\n')

    # ---------------- figures ----------------
    def style(ax):
        for side in ('top', 'right'):
            ax.spines[side].set_visible(False)
        ax.grid(color=GRID, lw=0.6)
        ax.set_axisbelow(True)

    fig, axes = plt.subplots(1, 2, figsize=(12, 8), gridspec_kw={'width_ratios': [1, 1.1]})
    for ax, (xl, yl, title) in zip(axes, [((2.4, 7.0), (-5.0, 5.0), 'Full view'),
                                          ((3.2, 5.6), (1.4, 4.4), 'Upper divertor (zoom)')]):
        lim = E['ref']['rzlim']
        ax.plot(lim[:, 0], lim[:, 1], color=INK2, lw=1.2, label='limiter (from gEQDSK)')
        for k in cases:
            for i, seg in enumerate(separatrix_pieces(E[k])):
                ax.plot(seg[:, 0], seg[:, 1], color=COLOR[k], ls=STYLE[k], lw=1.8,
                        label=f'{LABEL[k]}, ψ̂ = 1' if i == 0 else None)
            for xp in M[k]['x_points']:
                ax.plot(xp['R'], xp['Z'], 'x', color=COLOR[k], ms=10, mew=2.5)
            for xp in M[k]['secondary_nulls']:
                if xp['is_null']:
                    ax.plot(xp['R'], xp['Z'], 'o', mfc='none', mec=COLOR[k], ms=10, mew=2)
        if 'tmtx_sec' in cases:
            from run_flattop import SECONDARY_X_POINTS
            sx = np.array(SECONDARY_X_POINTS)
            ax.plot(sx[:, 0], sx[:, 1], '+', color=INK, ms=14, mew=1.5,
                    label='secondary X-point targets' if ax is axes[0] else None)
        ax.set_xlim(*xl); ax.set_ylim(*yl); ax.set_aspect('equal')
        ax.set_xlabel('R [m]'); ax.set_ylabel('Z [m]'); ax.set_title(title, fontsize=10)
        style(ax)
    axes[0].legend(fontsize=8, frameon=False, loc='lower left')
    fig.suptitle('ARC V3A separatrix: × main X-points, ○ outer-divertor nulls found', fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, 'separatrix_tmtx_vs_geqdsk.png'), dpi=140)
    plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True)
    panels = [('pres', 1e-6, 'Pressure p [MPa]'), ('qpsi', 1.0, 'Safety factor q'),
              ('ffprim', 1.0, "FF′ [T² m² / (Wb/rad)]"), ('pprime', 1e-6, "p′ [MPa / (Wb/rad)]")]
    for ax, (k2, sc, lab) in zip(axes.flat, panels):
        for k in cases:
            ax.plot(E[k]['psi_1d'], np.abs(np.asarray(E[k][k2]))*sc, color=COLOR[k], ls=STYLE[k], lw=2, label=LABEL[k])
        ax.set_ylabel(lab); style(ax)
    for ax in axes[1]:
        ax.set_xlabel('Normalized poloidal flux ψ̂ (0 = axis, 1 = boundary)')
    axes[0, 0].legend(frameon=False, fontsize=9)
    fig.suptitle('Profiles vs normalized poloidal flux (magnitudes; the files use different sign conventions)', fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, 'profiles_tmtx_vs_geqdsk.png'), dpi=140)
    plt.close(fig)

    # TORAX kinetic profiles, one panel per quantity, one line per TM+TX run
    fig, axes = plt.subplots(1, 4, figsize=(16, 3.8))
    for ax, (n, lab, sc) in zip(axes, (('n_e', r'$n_e$ [$10^{20}$ m$^{-3}$]', 1e-20), ('T_e', r'$T_e$ [keV]', 1.0),
                                       ('T_i', r'$T_i$ [keV]', 1.0), ('j_total', r'$j_{total}$ [MA/m²]', 1e-6))):
        for k in cases[1:]:
            p = J[k]['tx_profiles'].get(n, {})
            if 'x' in p:
                ax.plot(p['x'], np.asarray(p['y'], dtype=float)*sc, color=COLOR[k], ls=STYLE[k], lw=2, label=SHORT[k])
        ax.set_xlabel('ρ (normalized toroidal flux radius)'); ax.set_ylabel(lab); style(ax)
    axes[0].legend(frameon=False, fontsize=8)
    fig.suptitle('TORAX final profiles. A gEQDSK has no density or temperature to compare.', fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, 'torax_final_profiles.png'), dpi=140)
    plt.close(fig)

    if len(coils) == 2:
        names = list(coils['tmtx'])
        x = np.arange(len(names))
        fig, ax = plt.subplots(figsize=(13, 4.5))
        for i, k in enumerate(('tmtx', 'tmtx_sec')):
            ax.bar(x + (i - 0.5)*0.38, [coils[k][c] for c in names], width=0.36, color=COLOR[k], label=LABEL[k],
                   edgecolor='white', linewidth=0.6)
        for y in (COIL_LIMIT_MAT, -COIL_LIMIT_MAT):
            ax.axhline(y, color=INK2, ls='--', lw=1)
        ax.text(len(names) - 0.5, COIL_LIMIT_MAT, 'notebook coil limit ±50 MA-turns', ha='right', va='bottom',
                fontsize=8, color=INK2)
        ax.axhline(0, color=INK2, lw=0.8)
        ax.set_xticks(x); ax.set_xticklabels(names)
        ax.set_ylabel('Coil current [MA-turns]'); style(ax)
        ax.legend(frameon=False, fontsize=9, loc='lower right')
        ax.set_title('Coil currents at the end of the coupled run', fontsize=10)
        fig.tight_layout()
        fig.savefig(os.path.join(FIG, 'coil_currents_tmtx_cases.png'), dpi=140)
        plt.close(fig)
    print('\n'.join(L))


if __name__ == '__main__':
    import sys
    sys.path.insert(0, HERE)
    main()
