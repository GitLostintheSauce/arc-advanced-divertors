"""Tables and figures from the saved case records (no solving).

Usage (from the repo root, after run_diiid_shapehelper.py and run_iter_hmode.py):
    .venv/bin/python studies/diiid_iter_flux_function_realism/code/summarise.py

Writes results/summary.md and figures/*.png.
"""
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common  # noqa: E402
import run_diiid_shapehelper as diiid  # noqa: E402

# Categorical slots 1-4 of the dataviz reference palette (validated); line styles are a second cue
CAT = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100']
STYLES = ['-', '--', '-.', ':']
# Blue ramp steps 250-700, light -> dark (light end still clears 2:1 on the surface)
BLUE = ['#86b6ef', '#6da7ec', '#5598e7', '#3987e5', '#2a78d6', '#256abf', '#1c5cab', '#184f95', '#104281', '#0d366b']


def seq(n):
    """n evenly spaced steps of the blue ramp, light -> dark."""
    return [BLUE[int(round(i))] for i in np.linspace(0, len(BLUE) - 1, n)]
INK, MUTED, GRID = '#0b0b0b', '#52514e', '#e4e3df'
plt.rcParams.update({'font.size': 9, 'axes.edgecolor': MUTED, 'axes.labelcolor': INK,
                     'xtick.color': MUTED, 'ytick.color': MUTED, 'axes.grid': True,
                     'grid.color': GRID, 'grid.linewidth': 0.6, 'lines.linewidth': 2,
                     'axes.spines.top': False, 'axes.spines.right': False,
                     'figure.facecolor': '#fcfcfb', 'axes.facecolor': '#fcfcfb', 'legend.frameon': False})


def load(name):
    out = common.load_case(name)
    return out if out else (None, None)


def fmt(v, f='.2f'):
    return '' if v is None else (f'{v:{f}}' if isinstance(v, float) else str(v))


def row(rec, label=None):
    if rec is None:
        return None
    if not rec.get('converged', True):
        its = f" ({rec['maxits']} iterations)" if 'maxits' in rec else ''
        return f"| {label or rec['case']} | {rec['Ip_target']/1e6:.2f} | {rec['pax_target']/1e3:.1f} | not converged{its} | | | | | | | |"
    zc = ', '.join(f'{x:.2f}' for x in rec['ffp_zero_crossings']) or 'none'
    return (f"| {label or rec['case']} | {rec['Ip']/1e6:.2f} | {rec['P_ax']/1e3:.1f} | {fmt(rec['q_0'])} | "
            f"{fmt(rec['q_95'])} | {fmt(rec.get('l_i'))} | {fmt(rec.get('beta_N'))} | "
            f"{rec['ffp_axis']:+.2f} | {rec['ffp_min']:+.2f} | {zc} | "
            f"{fmt(rec['pp_edge_peak_over_core_peak'])} |")


HEADER = ('| case | Ip [MA] | P_ax [kPa] | q_0 | q_95 | l_i | beta_N | FF\'(axis) | min FF\' | '
          'FF\' zero crossings (psi_hat) | edge/core P\' peak |\n|---|---|---|---|---|---|---|---|---|---|---|')


def profile_panel(ax_ffp, ax_pp, items, colors, styles=None, title=None):
    for i, (label, arr) in enumerate(items):
        if arr is None:
            continue
        ls = styles[i] if styles else '-'
        ax_ffp.plot(arr['psi'], arr['FFp'], color=colors[i], ls=ls, label=label)
        ax_pp.plot(arr['psi'], arr['Pp']/1e3, color=colors[i], ls=ls, label=label)
    ax_ffp.axhline(0.0, color=MUTED, lw=0.8)
    ax_ffp.set_ylabel(r"FF' [T$^2$ m$^2$ / (Wb/rad)]")
    ax_pp.set_ylabel(r"P' [kPa / (Wb/rad)]")
    if title:
        ax_ffp.set_title(title, loc='left', fontsize=9, color=INK)


def diiid_figures(lines):
    os.makedirs(common.FIGURES, exist_ok=True)
    # 1. Published case and EFIT references
    names = [('shapehelper_as_published', 'shape-helper as published (power-law FF\', P\')'),
             ('efit_profiles_same_shape', 'same shape, EFIT FF\'/P\' shapes'),
             ('efit_192185_baseline', 'TokaMaker with shot 192185 inputs (baseline notebook)'),
             ('efit_192185_gfile', 'EFIT01 gEQDSK, shot 192185')]
    fig, axs = plt.subplots(1, 3, figsize=(12, 3.8))
    items = [(lab, load(n)[1]) for n, lab in names]
    profile_panel(axs[0], axs[1], items, CAT, STYLES)
    for i, (lab, arr) in enumerate(items):
        if arr is not None:
            axs[2].plot(arr['psi_q'], arr['q'], color=CAT[i], ls=STYLES[i], label=lab)
    axs[2].set_ylabel('q')
    for ax in axs:
        ax.set_xlabel(r'$\hat\psi$ (0 = axis, 1 = LCFS)')
        ax.set_xlim(0, 1)
    axs[2].legend(loc='upper left', fontsize=7.5)
    fig.suptitle('DIII-D: solved FF\', P\' and q, published shape-helper case vs EFIT references',
                 x=0.01, ha='left', color=INK)
    fig.tight_layout()
    fig.savefig(os.path.join(common.FIGURES, 'diiid_published_vs_efit.png'), dpi=150)
    plt.close(fig)

    lines += ['## DIII-D: published case and EFIT references', '', HEADER]
    lines += [r for r in (row(load(n)[0], lab) for n, lab in names) if r]
    lines += ['', 'EFIT01 gEQDSK row: Ip, P_ax, q and profiles from the file; l_i and beta_N not in the file.', '']

    # 2. Scans: profiles
    groups = [('Plasma current (pressure 7 kPa)', [k for k in diiid.SCAN if k.startswith('current_') and 'pressure' not in k], 'Ip'),
              ('Axis pressure (current 0.5 MA)', [k for k in diiid.SCAN if k.startswith('pressure_')], 'pax'),
              ('Axis pressure (current 1.2 MA)', [k for k in diiid.SCAN if k.startswith('current_1.20MA_pressure')], 'pax')]
    fig, axs = plt.subplots(len(groups), 2, figsize=(10, 3.2*len(groups)))
    for g, (title, keys, _) in enumerate(groups):
        keys = (['shapehelper_as_published'] if g < 2 else []) + keys
        items = []
        for k in keys:
            rec, arr = load(k)
            if rec is None or arr is None:      # not converged: no profiles (listed in summary.md)
                continue
            lab = f"{rec['Ip_target']/1e6:.2f} MA, {rec['pax_target']/1e3:.0f} kPa"
            items.append((lab if rec.get('converged') else lab + ' (not converged)', arr))
        profile_panel(axs[g, 0], axs[g, 1], items, seq(len(items)),
                      title=title)
        axs[g, 1].legend(fontsize=7, loc='upper right')
    for ax in axs[-1]:
        ax.set_xlabel(r'$\hat\psi$')
    fig.suptitle('DIII-D shape-helper: solved FF\' and P\' as current or pressure changes (profile shapes fixed)',
                 x=0.01, ha='left', color=INK)
    fig.tight_layout()
    fig.savefig(os.path.join(common.FIGURES, 'diiid_scan_profiles.png'), dpi=150)
    plt.close(fig)

    # 3. Pressure scans: scalars against the solved beta_p, so both currents share one axis
    fig, axs = plt.subplots(1, 4, figsize=(13, 3.3))
    for g, (title, keys, var) in enumerate(groups):
        if var == 'Ip':
            continue
        keys = (['shapehelper_as_published'] if g < 2 else []) + keys
        recs = [load(k)[0] for k in keys]
        failed = [r for r in recs if r and not r.get('converged')]
        recs = [r for r in recs if r and r.get('converged')]
        x = [r['beta_pol']/100 for r in recs]
        lab = '0.5 MA' if g == 1 else '1.2 MA'
        ffp_ref = recs[0]['ffp_axis']       # the 7 kPa case at this current
        ys = (('q_0', [r['q_0'] for r in recs]), ('l_i', [r['l_i'] for r in recs]),
              (r'$\beta_N$', [r['beta_N'] for r in recs]),
              ("FF'(axis) / FF'(axis) at 7 kPa", [r['ffp_axis']/ffp_ref for r in recs]))
        for ax, (ylabel, y) in zip(axs, ys):
            ax.plot(x, y, 'o' + STYLES[g-1], color=CAT[g-1], ms=5, label=lab)
            ax.set_ylabel(ylabel)
        if failed:
            lo = min(r['pax_target'] for r in failed)/1e3
            axs[0].text(x[-1], recs[-1]['q_0'], f'  no solution\n  from {lo:.0f} kPa', fontsize=7,
                        color=MUTED, va='center')
    axs[3].axhline(0.0, color=MUTED, lw=0.8)
    for ax in axs:
        ax.set_xlabel(r'$\beta_p$ (solved)')
    axs[1].legend(title='plasma current', fontsize=8)
    fig.suptitle('DIII-D shape-helper: axis pressure scans at fixed profile shapes, plotted against poloidal beta',
                 x=0.01, ha='left', color=INK)
    fig.tight_layout()
    fig.savefig(os.path.join(common.FIGURES, 'diiid_pressure_scan_scalars.png'), dpi=150)
    plt.close(fig)

    lines += ['## DIII-D: one input changed from the published case', '', HEADER]
    lines += [r for r in (row(load(k)[0]) for k in diiid.SCAN) if r]
    reruns = sorted(f[:-5] for f in os.listdir(os.path.join(common.RESULTS, 'cases')) if f.endswith('_maxits300.json'))
    if reruns:
        lines += ['', 'Non-converged cases rerun with 300 nonlinear iterations (the notebook uses the default, 40):', '', HEADER]
        lines += [row(load(k)[0]) for k in reruns]
    lines += ['', 'Cases named current_1.20MA_pressure_* change two inputs (current and pressure) from the published case.', '']


def iter_figures(lines):
    names = [('iter_hmode_wheel_no_isolate', 'OFT 26.9 wheel, isolate_edge_jBS=False'),
             ('iter_hmode_wheel_isolate', 'OFT 26.9 wheel, isolate_edge_jBS=True'),
             ('iter_hmode_oft_8fc4afe_no_isolate', 'bootstrap.py @8fc4afe (PR #316), False'),
             ('iter_hmode_oft_8fc4afe_isolate', 'bootstrap.py @8fc4afe (PR #316), True')]
    items = [(lab, load(n)[1]) for n, lab in names]
    if all(a is None for _, a in items):
        return
    fig, axs = plt.subplots(1, 3, figsize=(12, 3.8))
    profile_panel(axs[0], axs[1], items, CAT, STYLES)
    ref = next(a for _, a in items if a is not None)
    # P' implied by the input kinetic pressure, on the solved equilibrium's flux range
    rec0 = next(load(n)[0] for n, _ in names if load(n)[0])
    dpsi = abs(rec0['psi_axis'] - rec0['psi_lcfs'])
    pp_in = -np.gradient(ref['pressure_input'], ref['psi_input'], edge_order=2)/dpsi
    axs[1].plot(ref['psi_input'], pp_in/1e3, color=MUTED, lw=1, label='-dP/dpsi of input pressure')
    h, lab = axs[1].get_legend_handles_labels()
    fig.legend(h, lab, loc='lower center', ncol=3, fontsize=7.5)
    axs[2].plot(ref['psi_input'], ref['ne']/1e20, color=CAT[0], label=r'$n_e$ [$10^{20}$ m$^{-3}$]')
    axs[2].plot(ref['psi_input'], ref['Te']/1e4, color=CAT[1], ls='--', label=r'$T_e$ [10 keV]')
    axs[2].axvline(0.965 - 0.35/2, color=MUTED, lw=0.8, ls=':')
    axs[2].text(0.965 - 0.35/2, 0.05, ' ne pedestal start\n (0.965 - 0.35/2)', fontsize=7, color=MUTED,
                transform=axs[2].get_xaxis_transform())
    axs[2].legend(fontsize=7)
    for ax in axs:
        ax.set_xlabel(r'$\hat\psi$')
        ax.set_xlim(0, 1)
    fig.suptitle('ITER H-mode example: solved FF\' and P\', and the input kinetic profiles', x=0.01, ha='left', color=INK)
    fig.tight_layout(rect=(0, 0.1, 1, 1))
    fig.savefig(os.path.join(common.FIGURES, 'iter_hmode_profiles.png'), dpi=150)
    plt.close(fig)
    lines += ['## ITER H-mode example', '', HEADER]
    lines += [r for r in (row(load(n)[0], lab) for n, lab in names) if r]
    lines += ['', 'P\' argmax and FF\' zero crossings are in results/cases/*.json.', '']


def main():
    lines = ['# Summary tables', '',
             'Generated by code/summarise.py from results/cases/*.json. FF\' and P\' are TokaMaker\'s solved',
             'profiles (get_profiles), psi_hat = 0 on axis. "edge/core P\' peak" is the largest P\' at psi_hat > 0.85',
             'divided by the largest P\' inside it: above 1 means an edge pedestal dominates P\'.', '']
    diiid_figures(lines)
    iter_figures(lines)
    with open(os.path.join(common.RESULTS, 'summary.md'), 'w') as fh:
        fh.write('\n'.join(lines) + '\n')
    print('\n'.join(lines))


if __name__ == '__main__':
    main()
