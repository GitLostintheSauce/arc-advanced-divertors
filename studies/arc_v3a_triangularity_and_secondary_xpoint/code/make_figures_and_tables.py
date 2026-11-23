"""Figures and tables for the ARC V3A shape-target study.

Reads the per-case JSON (studies/arc_v3a_triangularity_and_secondary_xpoint/results/) and flux arrays (work/arc_v3a_triangularity_and_secondary_xpoint/arrays/)
written by run_equilibria.py. Figures go to studies/arc_v3a_triangularity_and_secondary_xpoint/figures/, tables to studies/arc_v3a_triangularity_and_secondary_xpoint/results/:
    A_vs_B_triangularity_separatrix.png     separatrix overlays, full view + upper divertor zoom
    A_vs_C_secondary_xpoint_separatrix.png
    coil_currents_by_case.png          coil currents per coil (MA-turns and kA per turn), cases A/B/C
    summary_by_case.md               requested vs achieved shape, X-points, flux, strike points, convergence
    coil_currents_by_case.md                  coil currents and differences with units
    shape_measured_from_boundary.json            delta/kappa measured from the traced boundary (see note in report)

Usage (from the repo root):
    .venv/bin/python studies/arc_v3a_triangularity_and_secondary_xpoint/code/make_figures_and_tables.py
"""
import json
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))           # .../code
STUDY = os.path.dirname(HERE)                                  # this study's folder
REPO = os.path.dirname(os.path.dirname(STUDY))
RES = os.path.join(STUDY, 'results')                           # tracked: JSON, tables
FIG = os.path.join(STUDY, 'figures')                           # tracked: PNGs
SCRATCH = os.path.join(REPO, 'work', 'arc_v3a_triangularity_and_secondary_xpoint')   # gitignored
ARR = os.path.join(SCRATCH, 'arrays')

MAIN = ['A_baseline_delta065', 'B_delta060', 'C_secondary_xpoint_out_5cm']
LABEL = {'A_baseline_delta065': 'A: baseline (δ=0.65)', 'B_delta060': 'B: δ=0.60',
         'C_secondary_xpoint_out_5cm': 'C: secondary X R+5 cm'}
COLOR = {'A_baseline_delta065': '#2a78d6', 'B_delta060': '#eb6834', 'C_secondary_xpoint_out_5cm': '#1baf7a'}
STYLE = {'A_baseline_delta065': '-', 'B_delta060': '--', 'C_secondary_xpoint_out_5cm': ':'}
INK, INK2, GRID = '#0b0b0b', '#52514e', '#e4e3df'


def load(n):
    with open(os.path.join(RES, f'{n}.json')) as fh:
        return json.load(fh)


def arrays(n):
    return np.load(os.path.join(ARR, f'{n}.npz'), allow_pickle=True)


def measured_shape(n):
    """kappa and delta from the traced psi_hat=0.9999 boundary, each side measured alike.

    OFT's get_stats() treats the 'active' X-point specially, so in a near-balanced double null
    its upper/lower values swap depending on which null is marginally active.
    """
    lc = arrays(n)['lcfs']
    rmin, rmax = lc[:, 0].min(), lc[:, 0].max()
    R_geo, a = 0.5*(rmin + rmax), 0.5*(rmax - rmin)
    up, lo = lc[np.argmax(lc[:, 1])], lc[np.argmin(lc[:, 1])]
    return {'R_geo': R_geo, 'a': a,
            'kappa': (up[1] - lo[1])/(2*a), 'kappaU': up[1]/a, 'kappaL': -lo[1]/a,
            'deltaU': (R_geo - up[0])/a, 'deltaL': (R_geo - lo[0])/a,
            'delta': 0.5*((R_geo - up[0]) + (R_geo - lo[0]))/a,
            'top': up.tolist(), 'bottom': lo.tolist()}


def strike_matches(rec):
    """Nearest achieved strike point to each strike target."""
    st = np.array(rec['strike_points'])
    out = []
    for tgt in np.array(rec['strike_targets']):
        i = int(np.argmin(np.hypot(*(st - tgt).T)))
        out.append({'target': tgt.tolist(), 'achieved': st[i].tolist(),
                    'miss_cm': float(100*np.hypot(*(st[i] - tgt)))})
    return out


def style_axes(ax):
    ax.set_facecolor('white')
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    ax.tick_params(colors=INK2)


def separatrix_panel(ax, names, xlim=None, ylim=None, legend=True):
    base = arrays(names[0])
    lim = base['lim']
    ax.plot(lim[:, 0], lim[:, 1], color=INK2, lw=1.2)
    for n in names:
        a = arrays(n)
        rec = load(n)
        ax.tricontour(a['r'][:, 0], a['r'][:, 1], a['lc'], a['psi'], levels=[1.0],
                      colors=COLOR[n], linestyles=STYLE[n], linewidths=1.8)
        ax.plot([], [], color=COLOR[n], ls=STYLE[n], lw=1.8, label=f'{LABEL[n]} separatrix')
        xt = np.array(rec['x_point_targets'])
        ax.plot(xt[:, 0], xt[:, 1], 'x', color=COLOR[n], ms=10, mew=2)
        nn = np.array([[q['R'], q['Z']] for q in rec['all_nulls']])
        ax.plot(nn[:, 0], nn[:, 1], 's', mfc='none', mec=COLOR[n], ms=9, mew=1.6)
        st = np.array(rec['strike_points'])
        ax.plot(st[:, 0], st[:, 1], 'o', mfc='none', mec=COLOR[n], ms=8, mew=1.6)
    stt = np.array(load(names[0])['strike_targets'])
    ax.plot(stt[:, 0], stt[:, 1], 'D', mfc='none', mec=INK, ms=7, mew=1.2)
    # legend proxies for marker meanings
    ax.plot([], [], 'x', color=INK2, ms=9, mew=2, ls='none', label='X-point target')
    ax.plot([], [], 's', mfc='none', mec=INK2, ms=8, ls='none', label='null achieved')
    ax.plot([], [], 'o', mfc='none', mec=INK2, ms=8, ls='none', label='strike point achieved')
    ax.plot([], [], 'D', mfc='none', mec=INK, ms=7, ls='none', label='strike target')
    ax.set_aspect('equal')
    if xlim:
        ax.set_xlim(*xlim)
    if ylim:
        ax.set_ylim(*ylim)
    ax.set_xlabel('R [m]')
    ax.set_ylabel('Z [m]')
    style_axes(ax)
    if legend:
        ax.legend(fontsize=7.5, loc='center left', frameon=False)


def equilibrium_figure(pair, fname, title):
    fig, axes = plt.subplots(1, 2, figsize=(13, 7.5), gridspec_kw={'width_ratios': [1, 1.25]})
    separatrix_panel(axes[0], pair, xlim=(2.0, 7.2), ylim=(-4.8, 4.8), legend=True)
    axes[0].set_title('Full view (ψ̂ = 1 contour)', fontsize=10)
    separatrix_panel(axes[1], pair, xlim=(3.2, 5.5), ylim=(1.8, 4.4), legend=False)
    axes[1].set_title('Upper divertor (zoom)', fontsize=10)
    fig.suptitle(title, fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, fname), dpi=140)
    plt.close(fig)


def coil_figure():
    recs = {n: load(n) for n in MAIN}
    coils = list(recs['A_baseline_delta065']['coil_I_per_turn'])
    x = np.arange(len(coils))
    w = 0.26
    fig, axes = plt.subplots(2, 1, figsize=(13, 8.2), sharex=True)
    for ax, key, scale, ylab in ((axes[0], 'coil_Aturns', 1e-6, 'Total current [MA-turns]'),
                                 (axes[1], 'coil_I_per_turn', 1e-3, 'Current per turn [kA]')):
        for i, n in enumerate(MAIN):
            vals = np.array([recs[n][key][c] for c in coils])*scale
            ax.bar(x + (i - 1)*w, vals, width=w - 0.03, color=COLOR[n], label=LABEL[n],
                   edgecolor='white', linewidth=0.6)
        ax.axhline(0, color=INK2, lw=0.8)
        ax.set_ylabel(ylab)
        ax.grid(axis='y', color=GRID, lw=0.6)
        ax.set_axisbelow(True)
        style_axes(ax)
    axes[0].legend(frameon=False, fontsize=9, ncol=3, loc='upper right')
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(coils)
    axes[0].set_title('Coil currents by coil (notebook bound: ±1000 kA per turn on every coil, not drawn)',
                      fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, 'coil_currents_by_case.png'), dpi=140)
    plt.close(fig)


def main():
    os.makedirs(FIG, exist_ok=True)
    recs = {n: load(n) for n in MAIN}
    shp = {n: measured_shape(n) for n in MAIN}
    with open(os.path.join(RES, 'shape_measured_from_boundary.json'), 'w') as fh:
        json.dump(shp, fh, indent=1)

    equilibrium_figure(['A_baseline_delta065', 'B_delta060'], 'A_vs_B_triangularity_separatrix.png',
                       'Case A (δ requested 0.65) vs case B (δ requested 0.60): everything else fixed')
    equilibrium_figure(['A_baseline_delta065', 'C_secondary_xpoint_out_5cm'], 'A_vs_C_secondary_xpoint_separatrix.png',
                       'Case A vs case C (secondary X-point targets R = 5.00 → 5.05 m): everything else fixed')
    coil_figure()

    # ---------------- summary table ----------------
    A = recs['A_baseline_delta065']
    L = []
    L.append('| Quantity | ' + ' | '.join(LABEL[n] for n in MAIN) + ' |')
    L.append('|---|' + '---|'*len(MAIN))

    def row(name, fn):
        L.append(f'| {name} | ' + ' | '.join(fn(n) for n in MAIN) + ' |')

    row('Converged (iterations)', lambda n: f"yes ({recs[n]['nl_its']})" if recs[n]['converged'] else 'NO')
    row('Topology', lambda n: 'diverted' if recs[n]['diverted'] else 'limited')
    row('δ requested', lambda n: f"{recs[n]['delta_requested']:.2f}")
    row('δ achieved, OFT get_stats (U / L)',
        lambda n: f"{recs[n]['delta']:.3f} ({recs[n]['deltaU']:.3f} / {recs[n]['deltaL']:.3f})")
    row('δ achieved, measured from boundary (U / L)',
        lambda n: f"{shp[n]['delta']:.3f} ({shp[n]['deltaU']:.3f} / {shp[n]['deltaL']:.3f})")
    row('δ achieved, from X-point: (R_geo − R_x)/a',
        lambda n: f"{(shp[n]['R_geo'] - recs[n]['x_point_matches'][0]['solved'][0])/shp[n]['a']:.3f}")
    row('κ requested', lambda n: f"{recs[n]['kappa_requested']:.2f}")
    row('κ achieved, OFT get_stats', lambda n: f"{recs[n]['kappa']:.3f}")
    row('κ achieved, measured from boundary', lambda n: f"{shp[n]['kappa']:.3f}")
    row('a, R_geo [m] (boundary)', lambda n: f"{shp[n]['a']:.3f}, {shp[n]['R_geo']:.3f}")
    row('q95, l_i', lambda n: f"{recs[n]['q95']:.3f}, {recs[n]['l_i']:.4f}")
    names = ['upper primary', 'lower primary', 'upper secondary', 'lower secondary']
    for i, nm in enumerate(names):
        row(f'{nm} X-point: target → achieved [m] (miss)',
            lambda n, i=i: (f"{tuple(recs[n]['x_point_matches'][i]['target'])} → "
                            f"({recs[n]['x_point_matches'][i]['solved'][0]:.3f}, "
                            f"{recs[n]['x_point_matches'][i]['solved'][1]:+.3f}) "
                            f"({recs[n]['x_point_matches'][i]['miss_mm']:.1f} mm)"))
    for i, nm in enumerate(names):
        row(f'ψ̂ at {nm} null', lambda n, i=i: f"{recs[n]['x_point_matches'][i]['psi_hat']:.4f}")
    row('ψ̂ gap, secondary − primary (upper / lower)',
        lambda n: (f"{recs[n]['x_point_matches'][2]['psi_hat'] - recs[n]['x_point_matches'][0]['psi_hat']:+.4f} / "
                   f"{recs[n]['x_point_matches'][3]['psi_hat'] - recs[n]['x_point_matches'][1]['psi_hat']:+.4f}"))
    snames = ['upper inner', 'lower inner', 'upper outer', 'lower outer']
    for i, nm in enumerate(snames):
        row(f'{nm} strike: target → nearest achieved [m] (miss)',
            lambda n, i=i: (lambda m: f"{tuple(m['target'])} → ({m['achieved'][0]:.3f}, {m['achieved'][1]:+.3f}) "
                                      f"({m['miss_cm']:.1f} cm)")(strike_matches(recs[n])[i]))
    row('Number of separatrix strike points on wall', lambda n: str(len(recs[n]['strike_points'])))
    for g in ['shape', 'strike', 'primary_xpt', 'secondary_xpt']:
        row(f'Isoflux residual max |ψ̂ − 1|, {g} points',
            lambda n, g=g: f"{recs[n]['isoflux_residuals'][g]['max_abs_dev']:.4f}")
    row('Largest coil |I| / bound', lambda n: (lambda k: f"{recs[n]['coil_bound_fraction'][k]:.2f} ({k})")(
        max(recs[n]['coil_bound_fraction'], key=recs[n]['coil_bound_fraction'].get)))
    row('Σ|coil current| [MA-turns]', lambda n: f"{sum(abs(v) for v in recs[n]['coil_Aturns'].values())/1e6:.1f}")
    with open(os.path.join(RES, 'summary_by_case.md'), 'w') as fh:
        fh.write('\n'.join(L) + '\n')

    # ---------------- coil table ----------------
    coils = list(A['coil_I_per_turn'])
    C = ['| Coil | turns | A [kA/turn] | B [kA/turn] | C [kA/turn] | B−A [kA/turn] | C−A [kA/turn] | '
         'A [MA-turns] | B−A [MA-turns] | C−A [MA-turns] |',
         '|---|---|---|---|---|---|---|---|---|---|']
    for k in coils:
        I = {n: recs[n]['coil_I_per_turn'][k]/1e3 for n in MAIN}
        T = {n: recs[n]['coil_Aturns'][k]/1e6 for n in MAIN}
        C.append(f"| {k} | {A['coil_nturns'][k]} | {I['A_baseline_delta065']:+.1f} | {I['B_delta060']:+.1f} | "
                 f"{I['C_secondary_xpoint_out_5cm']:+.1f} | {I['B_delta060']-I['A_baseline_delta065']:+.1f} | "
                 f"{I['C_secondary_xpoint_out_5cm']-I['A_baseline_delta065']:+.1f} | {T['A_baseline_delta065']:+.2f} | "
                 f"{T['B_delta060']-T['A_baseline_delta065']:+.2f} | {T['C_secondary_xpoint_out_5cm']-T['A_baseline_delta065']:+.2f} |")
    with open(os.path.join(RES, 'coil_currents_by_case.md'), 'w') as fh:
        fh.write('Coil currents. "kA/turn" is the current in one turn (what TokaMaker solves for and what '
                 'the notebook bounds at ±1000 kA/turn); "MA-turns" is per-turn current × turns = total '
                 'current through the coil cross-section.\n\n')
        fh.write('\n'.join(C) + '\n')

    # Reproducibility check against the 4-thread run
    if os.path.exists(os.path.join(RES, 'A_baseline_delta065_4threads.json')):
        a4 = load('A_baseline_delta065_4threads')
        same = all(abs(a4['coil_I_per_turn'][k] - A['coil_I_per_turn'][k]) < 1e-6*max(1, abs(A['coil_I_per_turn'][k]))
                   for k in coils) and abs(a4['delta'] - A['delta']) < 1e-9
        print('A 1-thread vs 4-thread identical:', same)
    print('\n'.join(L))


if __name__ == '__main__':
    main()
