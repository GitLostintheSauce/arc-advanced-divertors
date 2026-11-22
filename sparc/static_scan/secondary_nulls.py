"""Where are the extra field nulls in the SPARC static scan, and how close are they to the plasma?

Re-solves each scan case exactly as sparc_static_scan.py does, then:
  1. runs the divertor metrics from OpenFUSIONToolkit PR #403 ("TokaMaker: Add divertor geometry
     metrics", fork commit d567ae8) on the solved equilibrium. That gives, for every X-point
     TokaMaker reports, its flux (psi_n), whether it is inside the wall, its distance to the wall
     and where its separatrix crosses the outboard midplane (OMP);
  2. adds plasma-distance measures: distance to the LCFS, to the primary separatrix (all branches
     inside the wall) and to the nearest main X-point, plus |B_p| and the saddle check at the null;
  3. runs an independent null search over the whole mesh (1 cm grid, then Newton refinement), to
     check whether a null exists in cases where TokaMaker reports none.

PR #403 is not in the pinned OFT 26.9 wheel. Its divertor.py is pure Python and only uses
equilibrium methods that 26.9 already has, so the one file is downloaded at a fixed commit into
work/sparc_static/inputs/ (SHA-256 checked) and loaded on its own. No pins change.

Usage (from the repo root):
    .venv/bin/python sparc/static_scan/secondary_nulls.py               # all cases + summary
    .venv/bin/python sparc/static_scan/secondary_nulls.py --case ip_m20  # one case
    .venv/bin/python sparc/static_scan/secondary_nulls.py --summary     # table + plot only

Outputs: results/secondary_nulls.json, results/secondary_nulls.md, results/secondary_nulls.png.
Per-case grids for the plot go to work/sparc_static/nulls/.
"""
import argparse
import hashlib
import importlib.util
import json
import os
import subprocess
import sys
import urllib.request

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import sparc_static_scan as sss  # noqa: E402

NULLS = os.path.join(sss.SCRATCH, 'nulls')
PR403_COMMIT = 'd567ae8d879a73aac328d73cd5b5bdbcab3f76a9'
DIVERTOR_PY = os.path.join(sss.INPUTS, f'divertor_pr403_{PR403_COMMIT[:7]}.py')
DIVERTOR_URL = ('https://raw.githubusercontent.com/GitLostintheSauce/OpenFUSIONToolkit/'
                f'{PR403_COMMIT}/src/python/OpenFUSIONToolkit/TokaMaker/divertor.py')
DIVERTOR_SHA = '21157132a8c5af9c917fe6bf5f7fbb4e6f12d2e12c5e2bc7b792d9bc2034f9e0'
GRID_DX = 0.01      # null-search grid spacing [m]
# |B_p| [T] below which a refined point counts as a null. Not machine zero: B comes from the
# gradient of the order-2 FE psi, which jumps slightly across element edges, so a null lying on an
# edge (e.g. a material interface) only refines to ~1 mT. Typical |B_p| nearby is ~0.4 T.
NULL_TOL = 5.E-3
NEAR_LCFS = 0.6    # extra saddles closer than this to the LCFS go in the table [m]


def load_divertor():
    """Download PR #403's divertor.py if needed, check it, and import it as a standalone module."""
    os.makedirs(sss.INPUTS, exist_ok=True)
    if not os.path.exists(DIVERTOR_PY):
        print('Downloading divertor.py from PR #403', flush=True)
        urllib.request.urlretrieve(DIVERTOR_URL, DIVERTOR_PY)
    with open(DIVERTOR_PY, 'rb') as fh:
        got = hashlib.sha256(fh.read()).hexdigest()
    if got != DIVERTOR_SHA:
        raise RuntimeError(f'{DIVERTOR_PY}: SHA-256 {got} does not match expected {DIVERTOR_SHA}')
    spec = importlib.util.spec_from_file_location('divertor_pr403', DIVERTOR_PY)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def min_dist(pt, curve):
    """Distance from a point to a polyline [n,2]."""
    if curve is None or len(curve) < 2:
        return np.nan
    a, b = curve[:-1], curve[1:]
    v = b - a
    t = np.clip(((pt - a)*v).sum(1)/np.maximum((v**2).sum(1), 1.E-30), 0.0, 1.0)
    return float(np.hypot(*(a + t[:, None]*v - pt).T).min())


def seg_dist(pt, segs):
    """Distance from a point to segments [n,2,2]."""
    if len(segs) == 0:
        return np.nan
    a, v = segs[:, 0], segs[:, 1] - segs[:, 0]
    t = np.clip(((pt - a)*v).sum(1)/np.maximum((v**2).sum(1), 1.E-30), 0.0, 1.0)
    return float(np.hypot(*(a + t[:, None]*v - pt).T).min())


def hessian_det(psi, pt, h=2.E-3):
    """det of the psi Hessian by central differences (< 0 for a saddle / X-point)."""
    R, Z = pt
    p = psi(np.array([[R+h, Z], [R-h, Z], [R, Z+h], [R, Z-h],
                      [R+h, Z+h], [R+h, Z-h], [R-h, Z+h], [R-h, Z-h], [R, Z]]))
    prr = (p[0] - 2*p[8] + p[1])/h**2
    pzz = (p[2] - 2*p[8] + p[3])/h**2
    prz = (p[4] - p[5] - p[6] + p[7])/(4*h**2)
    return float(prr*pzz - prz**2)


def find_nulls(field, bbox):
    """All points with B_p = 0 on the mesh: grid local minima of |B_p|, refined by least squares."""
    from scipy.optimize import least_squares
    R = np.arange(bbox[0], bbox[1] + GRID_DX/2, GRID_DX)
    Z = np.arange(bbox[2], bbox[3] + GRID_DX/2, GRID_DX)
    RR, ZZ = np.meshgrid(R, Z)
    B = field.B(np.column_stack((RR.ravel(), ZZ.ravel())))
    Bp = np.hypot(B[:, 0], B[:, 2]).reshape(RR.shape)
    psi = field.psi(np.column_stack((RR.ravel(), ZZ.ravel()))).reshape(RR.shape)
    nulls = []
    for i in range(1, Bp.shape[0] - 1):
        for j in range(1, Bp.shape[1] - 1):
            win = Bp[i-1:i+2, j-1:j+2]
            if not np.all(np.isfinite(win)) or Bp[i, j] > win.min():
                continue

            def res(x):
                b = field.B(x)[0]
                return np.array([b[0], b[2]]) if np.all(np.isfinite(b)) else np.array([1.E3, 1.E3])
            sol = least_squares(res, [RR[i, j], ZZ[i, j]], x_scale=0.01, xtol=1.E-12, ftol=1.E-14,
                                bounds=([RR[i, j] - 2*GRID_DX, ZZ[i, j] - 2*GRID_DX],
                                        [RR[i, j] + 2*GRID_DX, ZZ[i, j] + 2*GRID_DX]))
            bp = float(np.hypot(*res(sol.x)))
            if bp < NULL_TOL and not any(np.hypot(*(sol.x - n['RZ'])) < 5.E-3 for n in nulls):
                nulls.append({'RZ': sol.x, 'Bp': bp})
    return nulls, R, Z, psi


def mesh_regions():
    """Return a function giving (region name, allow_xpoints) for a point on the SPARC mesh."""
    from matplotlib.tri import Triangulation
    from OpenFUSIONToolkit.TokaMaker.meshing import load_gs_mesh
    pts, lc, reg, coil_dict, cond_dict = load_gs_mesh(sss.MESH)
    names = {1: ('plasma region', True)}
    for dct in (coil_dict, cond_dict):
        for k, v in dct.items():
            names[v['reg_id']] = (k, bool(v.get('allow_xpoints', False)))
    finder = Triangulation(pts[:, 0], pts[:, 1], lc).get_trifinder()

    def region_of(pt):
        cell = int(finder(pt[0], pt[1]))
        return ('outside mesh', False) if cell < 0 else names.get(int(reg[cell]), (f'region {reg[cell]}', False))
    return region_of


def analyse_case(name):
    """Solve one case and record every null with its distances to the plasma."""
    div = load_divertor()
    myOFT, mygs, coil_dict, rec = sss.setup_and_solve(name)
    out = {'case': name, 'mod': rec['mod'], 'converged': rec['converged'], 'pr403_commit': PR403_COMMIT}
    if not rec['converged']:
        _write(name, out)
        return
    equil = mygs._tMaker_equil
    field = div._TokaMakerField(equil)
    walls = [np.vstack((c, c[:1])) if np.linalg.norm(c[0] - c[-1]) > 0 else c for c in mygs.lim_contours]
    segs = div._wall_segments(walls)
    xpts, _ = mygs.get_xpoints()
    region_of = mesh_regions()

    # PR #403 divertor metrics
    try:
        met = div.compute_divertor_metrics(mygs)
        out['pr403_error'] = None
    except Exception as err:   # report rather than hide a metrics failure
        met = None
        out['pr403_error'] = f'{type(err).__name__}: {err}'
    psi_sep = met['psi_sep'] if met else float(field.psi(xpts[-1])[0])
    R_omp_sep = met['R_OMP'] if met else np.nan
    ipri = met['primary_xpoint'] if met else xpts.shape[0] - 1
    o_pt = np.array(mygs.o_point)

    # Main X-points: the reported null nearest each saddle target
    tgts = np.array(rec['saddle_targets'])
    main_idx = sorted({int(np.argmin(np.hypot(*(xpts - t).T))) for t in tgts})

    # LCFS and the primary separatrix (all contour branches of psi_sep on a grid inside the mesh)
    lcfs = mygs.trace_surf(0.999)
    lim_all = np.vstack(walls)
    bbox_wall = [lim_all[:, 0].min(), lim_all[:, 0].max(), lim_all[:, 1].min(), lim_all[:, 1].max()]
    mesh_r = mygs.r
    bbox_mesh = [max(mesh_r[:, 0].min(), 1.E-3), mesh_r[:, 0].max(), mesh_r[:, 1].min(), mesh_r[:, 1].max()]
    nulls, gR, gZ, gpsi = find_nulls(field, bbox_mesh)
    import contourpy
    cg = contourpy.contour_generator(gR, gZ, np.ma.masked_invalid(gpsi))
    # Primary separatrix pieces inside the wall, as segments [n,2,2]
    sep_segs = []
    for ln in cg.lines(psi_sep):
        inside = np.array([any(div._point_in_contour(p, w) for w in walls) for p in ln])
        keep = inside[:-1] & inside[1:]
        sep_segs.append(np.stack((ln[:-1][keep], ln[1:][keep]), axis=1))
    sep_segs = np.concatenate(sep_segs) if sep_segs else np.zeros((0, 2, 2))
    s_wall, _ = div._first_wall_hit(o_pt, o_pt + np.r_[1.E3, 0.0], segs)
    R_wall_omp = o_pt[0] + s_wall*1.E3
    psi_n_wall_omp = float(field.psi_norm(field.psi(np.r_[R_wall_omp, o_pt[1]])[0]))

    def describe(pt, src_idx=None):
        psi_x = float(field.psi(pt)[0])
        b = field.B(pt)[0]
        R_omp = div._omp_crossing(field, psi_x, o_pt, R_wall_omp)
        d = {
            'R': float(pt[0]), 'Z': float(pt[1]),
            'tokamaker_index': src_idx,
            'role': ('primary' if src_idx == ipri else 'main' if src_idx in main_idx else 'extra'),
            'psi': psi_x,
            'psi_n': float(field.psi_norm(psi_x)),
            'Bp_at_null_T': float(np.hypot(b[0], b[2])),
            'saddle': hessian_det(field.psi, pt) < 0.0,
            'inside_wall': any(div._point_in_contour(pt, w) for w in walls),
            'mesh_region': region_of(pt)[0],
            'xpoints_allowed_in_region': region_of(pt)[1],
            'distance_to_wall_m': float(div._distance_to_segments(pt, segs)),
            'distance_to_lcfs_m': min_dist(pt, lcfs),
            'distance_to_separatrix_m': seg_dist(pt, sep_segs),
            'distance_to_nearest_main_xpt_m': float(min(np.hypot(*(xpts[main_idx] - pt).T))),
            'R_OMP_of_its_separatrix_m': float(R_omp),
            'dR_OMP_from_primary_mm': float(1.E3*(R_omp - R_omp_sep)),
        }
        return d

    out['tokamaker_xpoints'] = [describe(xpts[i], i) for i in range(xpts.shape[0])]
    # Independent search: label each null found and flag any TokaMaker did not report
    out['grid_nulls'] = []
    for n in nulls:
        dists = np.hypot(*(xpts - n['RZ']).T)
        idx = int(np.argmin(dists)) if dists.min() < 5.E-3 else None
        d = describe(n['RZ'], idx)
        if np.hypot(*(n['RZ'] - o_pt)) < 5.E-3:
            d['role'] = 'magnetic axis'
        elif idx is None:
            d['role'] = 'not reported by TokaMaker'
        out['grid_nulls'].append(d)
    out['R_OMP_primary_m'] = float(R_omp_sep)
    out['R_wall_OMP_m'] = float(R_wall_omp)
    out['psi_n_wall_OMP'] = psi_n_wall_omp
    out['psi_sep'] = float(psi_sep)
    if met:
        # Strike-point-to-null distances from PR #403 for the primary legs
        out['pr403_targets'] = {k: {'R_target': float(t['R_target']), 'Z_target': float(t['Z_target']),
                                    'secondary_xpoint': int(t['secondary_xpoint']),
                                    'secondary_xpoint_distance_m': float(t['secondary_xpoint_distance'])}
                                for k, t in met['targets'].items() if t['primary']}
    os.makedirs(NULLS, exist_ok=True)
    np.savez(os.path.join(NULLS, f'{name}.npz'), R=gR, Z=gZ, psi=gpsi, psi_sep=psi_sep,
             lcfs=lcfs, wall=walls[0], xpts=xpts)
    _write(name, out)
    print(f'[{name}] {len(xpts)} TokaMaker nulls, {len(nulls)} grid nulls', flush=True)


def _write(name, out):
    os.makedirs(NULLS, exist_ok=True)
    with open(os.path.join(NULLS, f'{name}.json'), 'w') as fh:
        json.dump(out, fh, indent=1, default=float)


def summarise():
    """Collect per-case results into a JSON, a Markdown table and a plot."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    recs = {}
    for name in sss.CASES:
        path = os.path.join(NULLS, f'{name}.json')
        if os.path.exists(path):
            with open(path) as fh:
                recs[name] = json.load(fh)
    with open(os.path.join(sss.OUT, 'secondary_nulls.json'), 'w') as fh:
        json.dump(recs, fh, indent=1)

    def nearby(n):
        """Extra saddle within NEAR_LCFS of the LCFS (not the axis or a main X-point)."""
        return (n['role'] not in ('magnetic axis', 'main', 'primary') and n['saddle']
                and n['distance_to_lcfs_m'] < NEAR_LCFS)

    lines = ['| case | null (R, Z) [m] | reported by TokaMaker | ψ_n | inside wall | to wall [cm] | to LCFS [cm] '
             '| to separatrix [cm] | to main X-pt [cm] | mesh region | \\|B_p\\| left [mT] | other nulls |',
             '|---|---|---|---|---|---|---|---|---|---|---|---|']
    for name, rec in recs.items():
        if not rec['converged']:
            lines.append(f'| {name} | not converged | | | | | | | | | | |')
            continue
        rows = [n for n in rec['grid_nulls'] if nearby(n)]
        n_other = sum(1 for n in rec['grid_nulls'] if n['role'] not in ('magnetic axis', 'main', 'primary')
                      and not nearby(n))
        if not rows:
            lines.append(f'| {name} | none | | | | | | | | | | {n_other} |')
        for n in sorted(rows, key=lambda n: -n['Z']):
            lines.append(
                f"| {name} | ({n['R']:.3f}, {n['Z']:+.3f}) | {'yes' if n['tokamaker_index'] is not None else 'no'} | "
                f"{n['psi_n']:.3f} | {'yes' if n['inside_wall'] else 'no'} | {100*n['distance_to_wall_m']:.1f} | "
                f"{100*n['distance_to_lcfs_m']:.1f} | {100*n['distance_to_separatrix_m']:.1f} | "
                f"{100*n['distance_to_nearest_main_xpt_m']:.1f} | {n['mesh_region']} | "
                f"{1.E3*n['Bp_at_null_T']:.1f} | {n_other} |")
    md = '\n'.join(lines)
    with open(os.path.join(sss.OUT, 'secondary_nulls.md'), 'w') as fh:
        fh.write(f'Extra X-points (saddles of psi) within {100*NEAR_LCFS:.0f} cm of the LCFS, per case. '
                 'The magnetic axis and the two main X-points are left out. "other nulls" counts the '
                 'remaining nulls on the mesh: extrema inside coils and saddles further out. '
                 'Generated by secondary_nulls.py.\n\n' + md + '\n')

    # Plot: upper divertor of the baseline and of each case where TokaMaker reported an extra null
    from matplotlib.tri import Triangulation
    from OpenFUSIONToolkit.TokaMaker.meshing import load_gs_mesh
    pts, lc, reg, coil_dict, cond_dict = load_gs_mesh(sss.MESH)
    tri = Triangulation(pts[:, 0], pts[:, 1], lc)
    shade = {}
    for k in ('VESSEL_INNER', 'VESSEL_OUTER', 'VS_COVER_UPPER', 'VS_COVER_LOWER'):
        shade[cond_dict[k]['reg_id']] = 1.0     # conductors: X-points allowed
    shade[cond_dict['VAC_VESSEL']['reg_id']] = 0.3
    for k, v in coil_dict.items():
        shade[v['reg_id']] = 2.0
    cell_val = np.array([shade.get(int(r), np.nan) for r in reg])
    show = ['baseline'] + [n for n, r in recs.items() if n != 'baseline' and r['converged'] and any(
        g['role'] == 'extra' for g in r['tokamaker_xpoints'])]
    fig, axes = plt.subplots(1, len(show), figsize=(3.4*len(show), 3.9), sharey=True)
    for ax, name in zip(np.atleast_1d(axes), show):
        z = np.load(os.path.join(NULLS, f'{name}.npz'))
        rec = recs[name]
        ax.tripcolor(tri, np.ma.masked_invalid(cell_val), cmap='Greys', vmin=0.0, vmax=3.0, shading='flat')
        ax.plot(z['wall'][:, 0], z['wall'][:, 1], 'k-', lw=1.2)
        ax.contour(z['R'], z['Z'], np.ma.masked_invalid(z['psi']), levels=[float(z['psi_sep'])],
                   colors='tab:red', linewidths=1.0)
        ax.plot(z['lcfs'][:, 0], z['lcfs'][:, 1], 'tab:red', lw=2)
        for n in rec['grid_nulls']:
            if n['Z'] < 0.3:
                continue
            if n['role'] in ('main', 'primary'):
                ax.plot(n['R'], n['Z'], 'kx', ms=8, mew=2)
            elif nearby(n):
                rep = n['tokamaker_index'] is not None
                ax.plot(n['R'], n['Z'], 'o', color='tab:blue', ms=8, mew=2,
                        mfc='tab:blue' if rep else 'none')
                ax.contour(z['R'], z['Z'], np.ma.masked_invalid(z['psi']), levels=[n['psi']],
                           colors='tab:blue', linewidths=0.8, linestyles='--')
                near_dv1 = n['Z'] > 1.4
                ax.annotate(f"{100*n['distance_to_lcfs_m']:.0f} cm to LCFS\n{100*n['distance_to_wall_m']:.0f} cm behind wall",
                            (n['R'], n['Z']), xytext=(4, 30) if near_dv1 else (-20, -38),
                            textcoords='offset points', fontsize=7, ha='left' if near_dv1 else 'center',
                            bbox={'boxstyle': 'round,pad=0.2', 'fc': 'white', 'ec': 'none', 'alpha': 0.8},
                            arrowprops={'arrowstyle': '-', 'lw': 0.6})
        ax.set_xlim(1.2, 2.35)
        ax.set_ylim(0.85, 1.8)
        ax.set_aspect('equal')
        ax.set_title(f"{name} ({len(rec['tokamaker_xpoints'])} nulls reported)", fontsize=9)
        ax.set_xlabel('R (m)')
    np.atleast_1d(axes)[0].set_ylabel('Z (m)')
    fig.suptitle('Upper divertor. Red: LCFS and primary separatrix. Black line: first wall. Black x: main '
                 'X-point. Blue o: extra X-point (filled = reported by TokaMaker), dashed = its separatrix.\n'
                 'Grey: vacuum-vessel gap (light), vessel conductors (mid), coils (dark). '
                 'TokaMaker only reports X-points in the plasma region and the conductors.', fontsize=7.5)
    fig.tight_layout()
    fig.savefig(os.path.join(sss.OUT, 'secondary_nulls.png'), dpi=140)
    plt.close(fig)
    print(md)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--case', choices=list(sss.CASES))
    parser.add_argument('--summary', action='store_true')
    args = parser.parse_args()
    if args.case:
        sss.fetch_inputs()
        analyse_case(args.case)
        return
    if not args.summary:
        sss.fetch_inputs()
        load_divertor()
        os.makedirs(os.path.join(sss.SCRATCH, 'logs'), exist_ok=True)
        for name in sss.CASES:   # one TokaMaker per process
            log = os.path.join(sss.SCRATCH, 'logs', f'nulls_{name}.log')
            with open(log, 'w') as fh:
                proc = subprocess.run([sys.executable, __file__, '--case', name],
                                      stdout=fh, stderr=subprocess.STDOUT)
            print(f'{name}: exit {proc.returncode} (log: {log})', flush=True)
    summarise()


if __name__ == '__main__':
    main()
