"""Triangularity with the X-points moved along with it (cases D and E) vs pinned X-points (A and B).

Reads results/<case>.json and the flux arrays written by run_equilibria.py for
A_baseline_delta065, B_delta060, D_delta065_xpoints_R3863 and E_delta060_xpoints_R3921. Writes:
    figures/D_vs_E_triangularity_with_xpoints_moved.png   separatrix overlay, full view + upper divertor
    results/delta_with_xpoints_moved.md                   A/B vs D/E side by side
    results/delta_with_xpoints_moved.json                 the same numbers

Usage (from the repo root, after running cases D and E):
    .venv/bin/python studies/arc_v3a_triangularity_and_secondary_xpoint/code/run_equilibria.py --case D_delta065_xpoints_R3863
    .venv/bin/python studies/arc_v3a_triangularity_and_secondary_xpoint/code/run_equilibria.py --case E_delta060_xpoints_R3921
    .venv/bin/python studies/arc_v3a_triangularity_and_secondary_xpoint/code/compare_delta_with_xpoints_moved.py
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import make_figures_and_tables as M  # noqa: E402

A, B = 'A_baseline_delta065', 'B_delta060'
D, E = 'D_delta065_xpoints_R3863', 'E_delta060_xpoints_R3921'
M.LABEL.update({D: 'D: δ=0.65, X-points at R=3.863', E: 'E: δ=0.60, X-points at R=3.921'})
M.COLOR.update({D: '#2a78d6', E: '#eb6834'})
M.STYLE.update({D: '-', E: '--'})


def row(n):
    r, s = M.load(n), M.measured_shape(n)
    cur = r['coil_I_per_turn']
    big = max(cur, key=lambda k: abs(cur[k]))
    xr = r['x_point_targets'][0][0]
    return {
        'delta_requested': r['delta_requested'],
        'primary_x_R': xr,
        'converged': r['converged'],
        'delta_oft': r['delta'], 'delta_boundary': s['delta'],
        'delta_from_xpoint': (s['R_geo'] - xr)/s['a'],
        'kappa_boundary': s['kappa'], 'q95': r['q95'], 'l_i': r['l_i'], 'volume': r['vol'],
        'x_point_miss_mm_max': max(m['miss_mm'] for m in r['x_point_matches']),
        'strike_miss_cm': [m['miss_cm'] for m in M.strike_matches(r)],
        'shape_max_abs_psi_dev': r['isoflux_residuals']['shape']['max_abs_dev'],
        'max_coil': big, 'max_coil_kA_per_turn': cur[big]/1e3,
        'coil_kA_per_turn': {k: v/1e3 for k, v in cur.items()},
    }


def main():
    rows = {n: row(n) for n in (A, B, D, E)}
    with open(os.path.join(M.RES, 'delta_with_xpoints_moved.json'), 'w') as fh:
        json.dump(rows, fh, indent=1)

    M.equilibrium_figure([D, E], 'D_vs_E_triangularity_with_xpoints_moved.png',
                         'Case D (δ 0.65, X-points R = 3.863 m) vs case E (δ 0.60, X-points R = 3.921 m)')

    f3 = lambda v: f'{v:.3f}'
    lines = [
        '| Quantity | A: δ 0.65, X pinned 3.90 | B: δ 0.60, X pinned 3.90 | B − A | '
        'D: δ 0.65, X at 3.863 | E: δ 0.60, X at 3.921 | E − D |',
        '|---|---|---|---|---|---|---|',
    ]

    def add(label, key, fmt=f3, diff=True):
        v = {n: rows[n][key] for n in rows}
        dBA = fmt(v[B] - v[A]) if diff else ''
        dED = fmt(v[E] - v[D]) if diff else ''
        lines.append(f'| {label} | {fmt(v[A])} | {fmt(v[B])} | {dBA} | {fmt(v[D])} | {fmt(v[E])} | {dED} |')

    add('primary X-point target R [m]', 'primary_x_R')
    add('δ achieved (OFT)', 'delta_oft')
    add('δ achieved (from boundary)', 'delta_boundary')
    add('δ from X-point position', 'delta_from_xpoint')
    add('κ (from boundary)', 'kappa_boundary')
    add('q95', 'q95')
    add('l_i', 'l_i')
    add('X-point miss, worst of four [mm]', 'x_point_miss_mm_max', fmt=lambda v: f'{v:.2f}', diff=False)
    add('worst shape-point \\|ψ̂−1\\|', 'shape_max_abs_psi_dev')
    for i, lab in enumerate(['inner strike miss, upper [cm]', 'inner strike miss, lower [cm]',
                             'outer strike miss, upper [cm]', 'outer strike miss, lower [cm]']):
        for n in rows:
            rows[n][f'_s{i}'] = rows[n]['strike_miss_cm'][i]
        add(lab, f'_s{i}', fmt=lambda v: f'{v:.1f}')
    add('largest \\|I\\| per turn [kA]', 'max_coil_kA_per_turn', fmt=lambda v: f'{abs(v):.0f}', diff=False)
    mc = {n: rows[n]['max_coil'] for n in rows}
    lines.append(f'| coil carrying it | {mc[A]} | {mc[B]} |  | {mc[D]} | {mc[E]} |  |')
    lines.append('')
    lines.append('Coil currents per turn [kA], upper coils (lower are the same within 0.5 kA):')
    lines.append('')
    coils = [k for k in rows[A]['coil_kA_per_turn'] if k.endswith('U')]
    lines.append('| Case | ' + ' | '.join(coils) + ' |')
    lines.append('|---|' + '---|'*len(coils))
    for n in (A, B, D, E):
        lines.append(f'| {n} | ' + ' | '.join(f"{rows[n]['coil_kA_per_turn'][k]:+.0f}" for k in coils) + ' |')
    with open(os.path.join(M.RES, 'delta_with_xpoints_moved.md'), 'w') as fh:
        fh.write('\n'.join(lines) + '\n')
    print('\n'.join(lines))


if __name__ == '__main__':
    main()
