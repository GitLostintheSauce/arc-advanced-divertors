"""Compare the TokaMaker+TORAX ARC V3A flattop with the public ARC V3A gEQDSK.

Both equilibria are read as gEQDSK files (the TM+TX one is written by run_flattop.py on the same
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
EXTRACT = os.path.join(SCRATCH, 'run', 'extract')
TMTX_EQDSK = os.path.join(EXTRACT, 'tmtx_flattop_final.geqdsk')
TMTX_JSON = os.path.join(EXTRACT, 'tmtx_flattop.json')
REF_EQDSK = os.path.join(SCRATCH, 'inputs', 'geqdsk-ARCv3a')
REF_URL = 'https://zenodo.org/records/19498373/files/geqdsk-ARCv3a'
REF_SHA256 = None   # filled on first download; see fetch_reference()

MU0 = 4e-7*np.pi
R0_DEVICE = 4.62
INK, INK2, GRID = '#0b0b0b', '#52514e', '#e4e3df'
C_REF, C_TMTX = '#2a78d6', '#eb6834'     # gEQDSK, TM+TX
LABEL = {'ref': 'gEQDSK (Zenodo, FreeGS)', 'tmtx': 'TokaMaker+TORAX flattop'}


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


def main():
    os.makedirs(RES, exist_ok=True)
    os.makedirs(FIG, exist_ok=True)
    ref_sha = fetch_reference()
    ref, tm = load(REF_EQDSK), load(TMTX_EQDSK)
    with open(TMTX_JSON) as fh:
        tj = json.load(fh)
    M = {}
    M['ref'], lc_ref = measure(ref)
    M['tmtx'], lc_tm = measure(tm)
    st = tj['tm_stats']
    tx = tj['tx_scalars']
    state = tj['tmtx_state']
    last = lambda k: (state.get(k) or [None])[-1] if isinstance(state.get(k), list) else None

    # Cross-check of the measuring functions against TokaMaker's own numbers for the TM+TX equilibrium
    check = {
        'kappa': (M['tmtx']['kappa'], st['kappa']), 'q95': (M['tmtx']['q95'], abs(st['q_95'])),
        'l_i': (M['tmtx']['l_i'], st['l_i']), 'beta_p': (M['tmtx']['beta_p'], st['beta_pol']/100),
        'beta_N': (M['tmtx']['beta_N'], st.get('beta_n')), 'volume_m3': (M['tmtx']['volume_m3'], st['vol']),
        'Ip_MA (boundary integral vs header)': (M['tmtx']['Ip_from_boundary_MA'], M['tmtx']['Ip_MA']),
    }
    tmtx_context = {
        'Q_fusion': tx.get('Q_fusion'), 'P_fusion_MW': (tx.get('P_fusion') or 0)/1e6 if tx.get('P_fusion') else None,
        'H98': tx.get('H98'), 'tau_E_s': tx.get('tau_E'), 'beta_N_torax': tx.get('beta_N'),
        'li3_torax': tx.get('li3'), 'Ip_torax_MA': (tx.get('Ip') or 0)/1e6 if tx.get('Ip') else None,
        'q95_torax': tx.get('q95'),
        'TM_vs_TX_final_tm_time': {k: last(k) for k in ('Ip_tm', 'Ip_tx', 'l_i_tm', 'l_i_tx', 'beta_N_tm',
                                                         'beta_N_tx', 'q95', 'q95_tm', 'pax_tm', 'vloop_tm', 'vloop_tx')},
    }
    out = {'reference_geqdsk': {'source': REF_URL, 'sha256': ref_sha, 'header': str(ref.get('case', ''))},
           'tmtx_notebook_commit': tj.get('notebook_commit'),
           'measured': M, 'measuring_check_vs_tokamaker_get_stats': check, 'tmtx_context': tmtx_context}
    with open(os.path.join(RES, 'comparison.json'), 'w') as fh:
        json.dump(out, fh, indent=1, default=float)

    # ---------------- table ----------------
    rows = [('Plasma current Ip [MA]', 'Ip_MA', '{:.2f}'), ('Vacuum field at R = 4.62 m [T]', 'B0_at_4.62m_T', '{:.2f}'),
            ('Magnetic axis R [m]', 'R_axis', '{:.3f}'), ('Geometric centre R_geo [m]', 'R_geo', '{:.3f}'),
            ('Minor radius a [m]', 'a', '{:.3f}'), ('Elongation κ', 'kappa', '{:.3f}'),
            ('Triangularity δ (mean, from boundary)', 'delta', '{:.3f}'), ('δ upper / lower', None, None),
            ('Triangularity δ (from X-points: (R_geo − R_x)/a)', 'delta_from_xpoints', '{:.3f}'),
            ('Plasma volume [m³]', 'volume_m3', '{:.1f}'), ('q0', 'q0', '{:.2f}'), ('q95', 'q95', '{:.2f}'),
            ('Axis pressure [MPa]', 'p_axis_MPa', '{:.2f}'), ('Volume-average pressure [MPa]', None, None),
            ('β_N', 'beta_N', '{:.2f}'), ('β_p', 'beta_p', '{:.3f}'), ('l_i', 'l_i', '{:.3f}'),
            ('Upper X-point (R, Z) [m]', None, None), ('Lower X-point (R, Z) [m]', None, None),
            ('Outer-divertor null near (5.0, ±3.7) m', None, None),
            ('Poloidal flux axis→boundary [Wb/rad]', 'psi_axis_to_boundary_Wb_per_rad', '{:.2f}')]
    L = ['| Quantity | gEQDSK (Zenodo) | TokaMaker+TORAX | Difference |', '|---|---|---|---|']
    for name, key, f in rows:
        r, t = M['ref'], M['tmtx']
        if key:
            d = t[key] - r[key]
            rel = f' ({100*d/r[key]:+.0f} %)' if r[key] else ''
            L.append(f'| {name} | {f.format(r[key])} | {f.format(t[key])} | {f.format(d) if f.startswith("{:.") else d}{rel} |')
        elif name.startswith('δ upper'):
            L.append(f"| {name} | {r['deltaU']:.3f} / {r['deltaL']:.3f} | {t['deltaU']:.3f} / {t['deltaL']:.3f} | |")
        elif name.startswith('Volume-average'):
            L.append(f"| {name} | {r['p_vol_avg_Pa']/1e6:.3f} | {t['p_vol_avg_Pa']/1e6:.3f} | "
                     f"{(t['p_vol_avg_Pa']-r['p_vol_avg_Pa'])/1e6:+.3f} ({100*(t['p_vol_avg_Pa']/r['p_vol_avg_Pa']-1):+.0f} %) |")
        elif name.startswith('Outer-divertor null'):
            def sn(m):
                parts = []
                for x in m['secondary_nulls']:
                    parts.append(f"({x['R']:.2f}, {x['Z']:+.2f}) ψ̂={x['psi_n']:.3f}" if x['is_null']
                                 else f"none (min B_p {1e3*x['Bp_T']:.0f} mT)")
                return '; '.join(parts)
            L.append(f"| {name} | {sn(r)} | {sn(t)} | |")
        elif 'X-point' in name:
            i = 0 if name.startswith('Upper') else 1
            xr, xt = r['x_points'][i], t['x_points'][i]
            L.append(f"| {name}, ψ̂ there | ({xr['R']:.3f}, {xr['Z']:+.3f}), {xr['psi_n']:.3f} | "
                     f"({xt['R']:.3f}, {xt['Z']:+.3f}), {xt['psi_n']:.3f} | "
                     f"{100*np.hypot(xt['R']-xr['R'], xt['Z']-xr['Z']):.1f} cm apart |")
    L += ['', 'Check of the measuring code against TokaMaker\'s own `get_stats()` for the TM+TX equilibrium:', '',
          '| Quantity | Measured here | TokaMaker | Relative difference |', '|---|---|---|---|']
    for k, (a, b) in check.items():
        L.append(f'| {k} | {fmt(a)} | {fmt(b)} | {fmt(None if (a is None or not b) else 100*(a/b-1), "{:+.1f} %")} |')
    L += ['', 'TokaMaker+TORAX run context (final state):', '', '| Quantity | Value |', '|---|---|']
    for k, v in tmtx_context.items():
        if k == 'TM_vs_TX_final_tm_time':
            for kk, vv in v.items():
                L.append(f'| {kk} (last TokaMaker time) | {fmt(vv, "{:.4g}")} |')
        else:
            L.append(f'| {k} | {fmt(v, "{:.4g}")} |')
    with open(os.path.join(RES, 'comparison_table.md'), 'w') as fh:
        fh.write('\n'.join(L) + '\n')

    # ---------------- figures ----------------
    def style(ax):
        for s in ('top', 'right'):
            ax.spines[s].set_visible(False)
        ax.grid(color=GRID, lw=0.6)
        ax.set_axisbelow(True)

    # 1) boundary and separatrix overlay
    fig, axes = plt.subplots(1, 2, figsize=(12, 8), gridspec_kw={'width_ratios': [1, 1.1]})
    for ax, (xl, yl, title) in zip(axes, [((2.4, 7.0), (-5.0, 5.0), 'Full view'),
                                          ((3.2, 5.6), (1.4, 4.4), 'Upper divertor (zoom)')]):
        lim = ref['rzlim']
        ax.plot(lim[:, 0], lim[:, 1], color=INK2, lw=1.2, label='limiter (from gEQDSK)')
        for key, e, lc, c, ls in (('ref', ref, lc_ref, C_REF, '-'), ('tmtx', tm, lc_tm, C_TMTX, '--')):
            for i, seg in enumerate(separatrix_pieces(e)):
                ax.plot(seg[:, 0], seg[:, 1], color=c, ls=ls, lw=1.8, label=f'{LABEL[key]} ψ̂ = 1' if i == 0 else None)
            for xp in M[key]['x_points']:
                ax.plot(xp['R'], xp['Z'], 'x', color=c, ms=10, mew=2.5)
            ax.plot(e['raxis'], e['zaxis'], '+', color=c, ms=12, mew=2)
        ax.set_xlim(*xl); ax.set_ylim(*yl); ax.set_aspect('equal')
        ax.set_xlabel('R [m]'); ax.set_ylabel('Z [m]'); ax.set_title(title, fontsize=10)
        style(ax)
    axes[0].legend(fontsize=8, frameon=False, loc='lower left')
    fig.suptitle('ARC V3A: TokaMaker+TORAX flattop vs public gEQDSK (× X-points, + magnetic axis)', fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, 'separatrix_tmtx_vs_geqdsk.png'), dpi=140)
    plt.close(fig)

    # 2) profiles vs normalized poloidal flux, one quantity per panel
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True)
    panels = [('pres', 1e-6, 'Pressure p [MPa]'), ('qpsi', 1.0, 'Safety factor q'),
              ('ffprim', 1.0, "FF′ [T² m² / (Wb/rad)]"), ('pprime', 1e-6, "p′ [MPa / (Wb/rad)]")]
    for ax, (k, s, lab) in zip(axes.flat, panels):
        for key, e, c, ls in (('ref', ref, C_REF, '-'), ('tmtx', tm, C_TMTX, '--')):
            y = np.abs(np.asarray(e[k]))*s
            ax.plot(e['psi_1d'], y, color=c, ls=ls, lw=2, label=LABEL[key])
        ax.set_ylabel(lab); style(ax)
    for ax in axes[1]:
        ax.set_xlabel('Normalized poloidal flux ψ̂ (0 = axis, 1 = boundary)')
    axes[0, 0].legend(frameon=False, fontsize=9)
    fig.suptitle('Profiles vs normalized poloidal flux (magnitudes; the two files use different sign conventions)', fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, 'profiles_tmtx_vs_geqdsk.png'), dpi=140)
    plt.close(fig)

    # 3) TORAX kinetic profiles (no gEQDSK equivalent: a gEQDSK stores only total pressure)
    prof = tj['tx_profiles']
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.8))
    for ax, names, lab, s in ((axes[0], ['n_e'], r'$n_e$ [$10^{20}$ m$^{-3}$]', 1e-20),
                              (axes[1], ['T_e', 'T_i'], 'T [keV]', 1.0),
                              (axes[2], ['j_total', 'j_bootstrap', 'j_ohmic'], 'j [MA/m²]', 1e-6)):
        for n, c in zip(names, [C_TMTX, '#1baf7a', '#6250d6']):
            p = prof.get(n, {})
            if 'x' in p:
                ax.plot(p['x'], np.asarray(p['y'], dtype=float)*s, color=c, lw=2, label=n)
        ax.set_xlabel('ρ (normalized toroidal flux radius)'); ax.set_ylabel(lab); style(ax)
        ax.legend(frameon=False, fontsize=8)
    fig.suptitle('TORAX final profiles (TokaMaker+TORAX flattop). A gEQDSK has no density or temperature to compare.', fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, 'torax_final_profiles.png'), dpi=140)
    plt.close(fig)
    print('\n'.join(L))


if __name__ == '__main__':
    main()
