"""Analyse the SPARC static scan results and make the report figure.

Reads results/*.json written by sparc_static_scan.py (no TokaMaker needed) and writes
    results/analysis.json   fitted scalings and shape sensitivities quoted in REPORT.md
    results/scaling.png     q95 and beta_p vs Ip, and q95 vs elongation

Usage (from the repo root):
    .venv/bin/python sparc/static_scan/analyze_scan.py
"""
import json
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, 'results')

INK = '#0b0b0b'        # primary text
INK_2 = '#52514e'      # secondary text, reference lines
SERIES = '#2a78d6'     # single data series
SURFACE = '#fcfcfb'

IP_CASES = ['ip_m20', 'ip_m10', 'baseline', 'ip_p10', 'ip_p20']
SHAPE_CASES = ['kappa_m5', 'xpt_Zin5', 'baseline', 'xpt_Zout5', 'kappa_p5']
SHAPE_LABELS = {'kappa_m5': 'shape −5 %', 'xpt_Zin5': 'X-pts in 5 cm', 'baseline': 'baseline',
                'xpt_Zout5': 'X-pts out 5 cm', 'kappa_p5': 'shape +5 %'}


def load(name):
    with open(os.path.join(RES, f'{name}.json')) as fh:
        return json.load(fh)


def main():
    base = load('baseline')

    # Current scan: fit log-log slopes
    ip = np.array([load(n)['Ip'] for n in IP_CASES])/1e6
    q95 = np.array([load(n)['q95'] for n in IP_CASES])
    bp = np.array([load(n)['beta_p'] for n in IP_CASES])
    li = np.array([load(n)['l_i'] for n in IP_CASES])
    slope_q = np.polyfit(np.log(ip), np.log(q95), 1)[0]
    slope_bp = np.polyfit(np.log(ip), np.log(bp), 1)[0]

    # Elongation cases: dq95/dkappa
    kap = np.array([load(n)['kappa'] for n in SHAPE_CASES])
    qk = np.array([load(n)['q95'] for n in SHAPE_CASES])
    dq_dk = np.polyfit(kap, qk, 1)[0]

    # X-point R shifts: d(delta)/dR_x per cm (average of the in and out cases)
    d_in = load('xpt_Rin5')['delta'] - base['delta']
    d_out = load('xpt_Rout5')['delta'] - base['delta']
    ddelta_dR = (d_out - d_in)/10.0  # per cm, from -5 cm to +5 cm

    # Pressure: beta_p response
    bp_ratio = {n: load(n)['beta_p']/base['beta_p'] for n in ['pax_m20', 'pax_p20']}

    # Measurement noise of delta in a balanced double null
    delta_ul_gap = base['deltaU'] - base['deltaL']

    # Total coil effort (sum of |MA-turns|), a crude measure dominated by the CS
    effort = {n: sum(abs(v) for v in load(n)['coil_Aturns'].values())/1e6
              for n in ['baseline'] + IP_CASES[:2] + IP_CASES[3:] + ['kappa_m5', 'kappa_p5', 'xpt_Rin5', 'xpt_Rout5']}

    out = {
        'q95_vs_Ip_loglog_slope': slope_q,
        'beta_p_vs_Ip_loglog_slope': slope_bp,
        'q95_times_Ip_MA': dict(zip(IP_CASES, (q95*ip).tolist())),
        'l_i_range_Ip_scan': [float(li.min()), float(li.max())],
        'dq95_dkappa': dq_dk,
        'ddelta_dRx_per_cm': ddelta_dR,
        'beta_p_ratio_pressure_cases': bp_ratio,
        'baseline_deltaU_minus_deltaL': delta_ul_gap,
        'coil_effort_MAturns': effort,
        'iterations': {n: load(n)['nl_its'] for n in
                       ['baseline', 'ip_m20', 'ip_m10', 'ip_p10', 'ip_p20', 'pax_m20', 'pax_p20',
                        'xpt_Rin5', 'xpt_Rout5', 'xpt_Zin5', 'xpt_Zout5', 'xpt_upper_up5',
                        'kappa_m5', 'kappa_p5']},
    }
    with open(os.path.join(RES, 'analysis.json'), 'w') as fh:
        json.dump(out, fh, indent=1)

    # Figure: three panels, one measure each (no shared or dual axes)
    plt.rcParams.update({'font.size': 10, 'axes.edgecolor': INK_2, 'axes.labelcolor': INK,
                         'xtick.color': INK_2, 'ytick.color': INK_2, 'text.color': INK})
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.0), facecolor=SURFACE)
    ip_ref = np.linspace(ip.min()*0.97, ip.max()*1.03, 50)

    ax = axes[0]
    ax.plot(ip_ref, q95[2]*ip[2]/ip_ref, '--', color=INK_2, lw=1.2, label=r'$\propto 1/I_p$')
    ax.plot(ip, q95, 'o-', color=SERIES, lw=2, ms=8, label='TokaMaker')
    ax.set_xlabel(r'$I_p$ [MA]')
    ax.set_ylabel(r'$q_{95}$')
    ax.set_title(f'(a) $q_{{95}}$ vs current (fitted slope {slope_q:.2f})', fontsize=10)
    ax.legend(frameon=False)

    ax = axes[1]
    ax.plot(ip_ref, bp[2]*(ip[2]/ip_ref)**2, '--', color=INK_2, lw=1.2, label=r'$\propto 1/I_p^2$')
    ax.plot(ip, bp, 'o-', color=SERIES, lw=2, ms=8, label='TokaMaker')
    ax.set_xlabel(r'$I_p$ [MA]')
    ax.set_ylabel(r'$\beta_p$')
    ax.set_title(f'(b) $\\beta_p$ vs current, fixed $p_{{axis}}$ (slope {slope_bp:.2f})', fontsize=10)
    ax.legend(frameon=False)

    ax = axes[2]
    order = np.argsort(kap)
    k_line = np.linspace(kap.min() - 0.01, kap.max() + 0.01, 10)
    fit = np.polyfit(kap, qk, 1)
    ax.plot(k_line, np.polyval(fit, k_line), '--', color=INK_2, lw=1.2,
            label=f'linear fit, {dq_dk:.1f} per unit κ')
    ax.plot(kap[order], qk[order], 'o', color=SERIES, ms=8, label='TokaMaker')
    for n, k, q in zip(SHAPE_CASES, kap, qk):
        ax.annotate(SHAPE_LABELS[n], (k, q), textcoords='offset points', xytext=(6, -12),
                    fontsize=8, color=INK_2)
    ax.set_xlabel(r'elongation $\kappa$')
    ax.set_ylabel(r'$q_{95}$')
    ax.set_title(r'(c) $q_{95}$ vs elongation at fixed $I_p$', fontsize=10)
    ax.legend(frameon=False, loc='upper left')

    for ax in axes:
        ax.set_facecolor(SURFACE)
        ax.grid(color='#e4e3df', lw=0.6)
        for side in ('top', 'right'):
            ax.spines[side].set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(RES, 'scaling.png'), dpi=140, facecolor=SURFACE)
    print(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
