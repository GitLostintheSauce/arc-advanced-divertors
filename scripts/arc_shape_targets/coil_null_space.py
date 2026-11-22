"""Diagnostic: is the large A->B coil-current swing a near-"invisible" coil combination?

Takes the solved coil currents from A_baseline.json, B_delta060.json and C_secR505.json,
runs vacuum solves (no plasma) and measures how much the flux from each coil set varies
along the baseline isoflux targets (shape + strike + X-point points). A uniform flux
offset along the boundary does not change the shape, so the shape-relevant number is the
spread (max - min) of flux along the targets.

Writes docs/arc_shape_targets/results/coil_null_space.json.

Usage (from the repo root, after shape_study.py has run):
    .venv/bin/python scripts/arc_shape_targets/coil_null_space.py
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from shape_study import BASE, MESH, RES, build_shape_targets  # noqa: E402


def main():
    from matplotlib.tri import Triangulation, LinearTriInterpolator
    from OpenFUSIONToolkit import OFT_env
    from OpenFUSIONToolkit.TokaMaker import TokaMaker
    from OpenFUSIONToolkit.TokaMaker.meshing import load_gs_mesh

    A = json.load(open(os.path.join(RES, 'A_baseline.json')))
    B = json.load(open(os.path.join(RES, 'B_delta060.json')))
    C = json.load(open(os.path.join(RES, 'C_secR505.json')))
    targets = build_shape_targets(dict(BASE))['iso']
    xpts = np.array(BASE['primary_x_points'] + BASE['secondary_x_points'])

    myOFT = OFT_env(nthreads=1)
    mygs = TokaMaker(myOFT)
    pts, lc, reg, coil_dict, cond_dict = load_gs_mesh(MESH)
    mygs.setup_mesh(pts, lc, reg)
    mygs.setup_regions(cond_dict=cond_dict, coil_dict=coil_dict)
    mygs.setup(order=2, F0=BASE['R0']*BASE['B0'])
    tri = Triangulation(mygs.r[:, 0], mygs.r[:, 1], mygs.lc)

    def flux_from(currents):
        mygs.set_coil_currents(currents)
        psi = np.zeros(mygs.np)
        mygs.vac_solve(psi=psi)
        f = LinearTriInterpolator(tri, psi)
        on_targets = np.array([float(f(*p)) for p in targets])
        # |grad psi| / R = |B_pol| at the X-point targets
        gx, gz = f.gradient(xpts[:, 0], xpts[:, 1])
        bpol = np.hypot(np.asarray(gx), np.asarray(gz))/xpts[:, 0]
        return {'spread_on_targets_Wb_per_rad': float(on_targets.max() - on_targets.min()),
                'std_on_targets_Wb_per_rad': float(on_targets.std()),
                'Bpol_at_xpoint_targets_T': bpol.tolist()}

    IA, IB, IC = A['coil_I_per_turn'], B['coil_I_per_turn'], C['coil_I_per_turn']
    dI = {k: IB[k] - IA[k] for k in IA}
    dIC = {k: IC[k] - IA[k] for k in IA}
    out = {'A_coils': flux_from(IA), 'B_coils': flux_from(IB), 'B_minus_A_coils': flux_from(dI),
           'C_minus_A_coils': flux_from(dIC)}
    # Each coil pair's share of the difference, alone
    pairs = sorted({k[:-1] for k in dI})
    out['B_minus_A_by_pair'] = {p: flux_from({p+'U': dI[p+'U'], p+'L': dI[p+'L']}) for p in pairs}
    out['C_minus_A_by_pair'] = {p: flux_from({p+'U': dIC[p+'U'], p+'L': dIC[p+'L']}) for p in pairs}
    out['dI_per_turn_A'] = dI
    out['dIC_per_turn_A'] = dIC
    out['note'] = ('psi is TokaMaker poloidal flux per radian [Wb/rad]. Spread = max - min of the '
                   'coil-only flux along the baseline isoflux targets.')
    with open(os.path.join(RES, 'coil_null_space.json'), 'w') as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps({k: v['spread_on_targets_Wb_per_rad'] if 'spread_on_targets_Wb_per_rad' in v else
                      {p: round(q['spread_on_targets_Wb_per_rad'], 3) for p, q in v.items()}
                      for k, v in out.items() if isinstance(v, dict) and not k.startswith('dI')}, indent=1))


if __name__ == '__main__':
    main()
