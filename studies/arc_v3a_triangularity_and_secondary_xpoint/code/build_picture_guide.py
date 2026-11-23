"""Build studies/arc_v3a_triangularity_and_secondary_xpoint/picture_guide.html, the interactive picture guide.

Takes the real simulation results (per-case JSON, flux arrays, mesh coil positions), turns them
into compact drawing data, and inserts it into picture_guide_template.html. The output is one
self-contained HTML file with no external dependencies, so it works offline.

Usage (from the repo root, after run_equilibria.py and check_coil_cancellation.py have run):
    .venv/bin/python studies/arc_v3a_triangularity_and_secondary_xpoint/code/build_picture_guide.py
    .venv/bin/python studies/arc_v3a_triangularity_and_secondary_xpoint/code/build_picture_guide.py --hosted PATH   # also write a copy for a
                                                                              # host that adds its own
                                                                              # <html>/<head>/<body>
"""
import argparse
import json
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from run_equilibria import ARR, BASE, MESH, RES  # noqa: E402

STUDY = os.path.dirname(HERE)
OUT = os.path.join(STUDY, 'picture_guide.html')
TEMPLATE = os.path.join(HERE, 'picture_guide_template.html')
CASES = ['A_baseline_delta065', 'B_delta060', 'C_secondary_xpoint_out_5cm']


def r3(a):
    return np.round(np.asarray(a, dtype=float), 3).tolist()


def load(n):
    with open(os.path.join(RES, f'{n}.json')) as fh:
        return json.load(fh)


def resample_closed(pts, n=180):
    """Resample a closed curve to n points by arc length, starting at max R, counter-clockwise."""
    pts = np.asarray(pts)
    if np.allclose(pts[0], pts[-1]):
        pts = pts[:-1]
    i0 = int(np.argmax(pts[:, 0]))
    pts = np.roll(pts, -i0, axis=0)
    # make counter-clockwise (positive signed area in R-Z)
    area = 0.5*np.sum(pts[:, 0]*np.roll(pts[:, 1], -1) - np.roll(pts[:, 0], -1)*pts[:, 1])
    if area < 0:
        pts = np.vstack([pts[:1], pts[1:][::-1]])
    closed = np.vstack([pts, pts[:1]])
    s = np.concatenate([[0], np.cumsum(np.hypot(*np.diff(closed, axis=0).T))])
    t = np.linspace(0, s[-1], n, endpoint=False)
    return np.column_stack([np.interp(t, s, closed[:, 0]), np.interp(t, s, closed[:, 1])])


def separatrix_legs(npz, lim, step=0.03):
    """psi_hat = 1 contour pieces inside the limiter, outside the core boundary, decimated."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.path import Path
    fig, ax = plt.subplots()
    cs = ax.tricontour(npz['r'][:, 0], npz['r'][:, 1], npz['lc'], npz['psi'], levels=[1.0])
    plt.close(fig)
    inside = Path(lim)
    pieces = []
    for path in cs.get_paths():
        v, codes = path.vertices, path.codes
        starts = [0] if codes is None else list(np.where(codes == 1)[0])
        starts.append(len(v))
        for a, b in zip(starts[:-1], starts[1:]):
            seg = v[a:b]
            keep = inside.contains_points(seg, radius=1e-3)
            # split into runs that are inside the vessel
            run = []
            for p, k in zip(seg, keep):
                if k:
                    run.append(p)
                elif len(run) > 1:
                    pieces.append(np.array(run)); run = []
                else:
                    run = []
            if len(run) > 1:
                pieces.append(np.array(run))
    out = []
    for p in pieces:
        d = np.concatenate([[0], np.cumsum(np.hypot(*np.diff(p, axis=0).T))])
        if d[-1] < 0.05:
            continue
        t = np.arange(0, d[-1], step)
        t = np.append(t, d[-1])
        out.append(r3(np.column_stack([np.interp(t, d, p[:, 0]), np.interp(t, d, p[:, 1])])))
    return out


def coil_boxes():
    """Bounding box of each coil region in the mesh."""
    from OpenFUSIONToolkit.TokaMaker.meshing import load_gs_mesh
    pts, lc, reg, coil_dict, _ = load_gs_mesh(MESH)
    boxes = {}
    for name, c in coil_dict.items():
        cells = lc[reg == c['reg_id']]
        p = pts[np.unique(cells)]
        boxes[name] = r3([p[:, 0].min(), p[:, 1].min(), p[:, 0].max(), p[:, 1].max()])
    return boxes


def requested_outline(delta):
    from OpenFUSIONToolkit.TokaMaker.util import create_isoflux_xpts
    c = BASE
    pts = create_isoflux_xpts(50, c['R0'], 0, c['a'], kappa_upper=c['kappa'], delta_upper=delta,
                              kappa_lower=c['kappa'], delta_lower=delta)
    return resample_closed(pts, 180)


def main(hosted=None):
    recs = {n: load(n) for n in CASES}
    arrs = {n: np.load(os.path.join(ARR, f'{n}.npz'), allow_pickle=True) for n in CASES}
    lim = np.asarray(arrs['A_baseline_delta065']['lim'])
    null = json.load(open(os.path.join(RES, 'coil_cancellation_check.json')))
    shp = json.load(open(os.path.join(RES, 'shape_measured_from_boundary.json')))

    data = {
        'limiter': r3(lim[::2] if len(lim) > 400 else lim),
        'coils': coil_boxes(),
        'requested': {f'{d:.2f}': r3(requested_outline(d)) for d in (0.55, 0.60, 0.65, 0.70)},
        'cases': {},
        'null_space': {
            'whole_coil_set': null['A_coils']['spread_on_targets_Wb_per_rad'],
            'B_minus_A_total': null['B_minus_A_coils']['spread_on_targets_Wb_per_rad'],
            'B_minus_A_by_pair': {p: v['spread_on_targets_Wb_per_rad']
                                  for p, v in null['B_minus_A_by_pair'].items()},
            'C_minus_A_total': null['C_minus_A_coils']['spread_on_targets_Wb_per_rad'],
        },
        'strike_targets': r3(recs['A_baseline_delta065']['strike_targets']),
    }
    for n in CASES:
        r = recs[n]
        data['cases'][n] = {
            'lcfs': r3(resample_closed(arrs[n]['lcfs'], 180)),
            'legs': separatrix_legs(arrs[n], lim),
            'x_targets': r3(r['x_point_targets']),
            'nulls': r3([[q['R'], q['Z']] for q in r['all_nulls']]),
            'strikes': r3(r['strike_points']),
            'coil_MAt': {k: round(v/1e6, 2) for k, v in r['coil_Aturns'].items()},
            'coil_kA': {k: round(v/1e3, 1) for k, v in r['coil_I_per_turn'].items()},
            'delta_req': r['delta_requested'],
            'delta_oft': round(r['delta'], 3),
            'delta_boundary': round(shp[n]['delta'], 3),
            'kappa': round(r['kappa'], 3),
            'q95': round(r['q95'], 2),
            'its': r['nl_its'],
            'converged': r['converged'],
            'strike_miss_cm': [],
        }
        st = np.array(r['strike_points'])
        for tgt in np.array(r['strike_targets']):
            d = np.hypot(*(st - tgt).T)
            data['cases'][n]['strike_miss_cm'].append(round(float(100*d.min()), 1))
    with open(TEMPLATE) as fh:
        html = fh.read()
    html = html.replace('/*__GUIDE_DATA__*/null', json.dumps(data, separators=(',', ':')))
    with open(OUT, 'w') as fh:
        fh.write(html)
    print(f'wrote {OUT} ({os.path.getsize(OUT)/1024:.0f} kB)')
    if hosted:
        body = re.sub(r'<!doctype html>\s*|</?html[^>]*>\s*|</?head>\s*|</?body>\s*|'
                      r'<meta charset="utf-8">\s*|<meta name="viewport"[^>]*>\s*', '', html)
        with open(hosted, 'w') as fh:
            fh.write(body)
        print(f'wrote {hosted}')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--hosted', help='also write a copy without <html>/<head>/<body> wrappers')
    main(ap.parse_args().hosted)
