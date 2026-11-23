"""DIII-D shape-helper example: reproduce it, compare with EFIT, and scan one input at a time.

The setup copies the physics cells of OFT's DIIID_shapehelper_ex.ipynb at commit 8fc4afe
(mesh, F0, VSC, Ip/pax targets, shape from xpoints_from_moments/create_isoflux_xpts, coil
regularization, power-law FF' and P', init_psi). Plot cells are left out.

Cases
  shapehelper_as_published   the notebook unchanged
  efit_profiles_same_shape   the notebook's shape and coils, but FF' and P' shapes from the EFIT
                             reconstruction of shot 192185 (same Ip and pax as published)
  efit_192185_baseline       DIIID_baseline_ex.ipynb's inverse solve: EFIT boundary, profiles,
                             Ip, pax and coil targets (TokaMaker with the real shot's inputs)
  scan cases                 one input changed from shapehelper_as_published (see SCAN)

Usage (from the repo root):
    .venv/bin/python studies/diiid_iter_flux_function_realism/code/run_diiid_shapehelper.py           # all cases
    .venv/bin/python studies/diiid_iter_flux_function_realism/code/run_diiid_shapehelper.py --case NAME
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common  # noqa: E402

RCENTR, BCENTR = 1.696, 1.92
PUBLISHED = {'Ip': 5.E5, 'pax': 7.E3, 'ffp': (1.5, 2.0), 'pp': (4.0, 1.0), 'profiles': 'power'}


def _scan():
    cases = {}
    for ip in (0.75, 1.0, 1.25, 1.5):
        cases[f'current_{ip:.2f}MA'] = {'Ip': ip*1.E6}
    for pax in (15, 30, 50, 80, 120):
        cases[f'pressure_{pax}kPa'] = {'pax': pax*1.E3}
    # Pressure scan again at a typical DIII-D H-mode current (two inputs differ from published)
    for pax in (7, 30, 50, 80, 100, 120, 150, 200):
        cases[f'current_1.20MA_pressure_{pax}kPa'] = {'Ip': 1.2E6, 'pax': pax*1.E3}
    for a in (1.0, 2.5, 4.0):
        cases[f'ffp_alpha_{a:g}'] = {'ffp': (a, PUBLISHED['ffp'][1])}
    for g in (1.0, 3.0):
        cases[f'ffp_gamma_{g:g}'] = {'ffp': (PUBLISHED['ffp'][0], g)}
    for a in (1.0, 2.0, 8.0):
        cases[f'pp_alpha_{a:g}'] = {'pp': (a, PUBLISHED['pp'][1])}
    for g in (0.5, 2.0):
        cases[f'pp_gamma_{g:g}'] = {'pp': (PUBLISHED['pp'][0], g)}
    return cases


SCAN = _scan()
CASES = {'shapehelper_as_published': {}, 'efit_profiles_same_shape': {'profiles': 'efit'},
         'efit_192185_baseline': {'baseline': True}, **SCAN}

# Coil targets of shot 192185 [A-turns], from DIIID_baseline_ex.ipynb at 8fc4afe
EFIT_COIL_TARGETS = {
    'ECOILA': -16030.961914062493, 'ECOILB': -15782.150390625,
    'F1A': 1999.5169719827586, 'F2A': 1031.1705953663793, 'F3A': -530.3450296336207,
    'F4A': -691.5127963362069, 'F5A': 12.465672986260776, 'F6A': -2142.877414772727,
    'F7A': 620.4805397727273, 'F8A': -975.8820716594828, 'F9A': 4302.279545454546,
    'F1B': 2213.2513469827586, 'F2B': 1316.7264278017242, 'F3B': -385.2729660560345,
    'F4B': -2295.533943965517, 'F5B': 6.891535265692344, 'F6B': -2792.180113636364,
    'F7B': 754.6947443181818, 'F8B': -886.030239762931, 'F9B': 4588.732102272727,
}


def efit_profile_shapes(eqdsk, npts=50):
    """EFIT FF' and P' resampled on 50 points, as DIIID_baseline_ex.ipynb does (signs as in the file)."""
    x_eq = np.linspace(0.0, 1.0, np.size(eqdsk['ffprim']))
    x = np.linspace(0.0, 1.0, npts)
    return ({'type': 'linterp', 'x': x, 'y': np.interp(x, x_eq, eqdsk['ffprim'])},
            {'type': 'linterp', 'x': x, 'y': np.interp(x, x_eq, eqdsk['pprime'])})


def shapehelper_shape(mygs):
    """Isoflux and saddle targets exactly as in the notebook's shape cells."""
    from OpenFUSIONToolkit.TokaMaker.util import create_isoflux_xpts, xpoints_from_moments
    lim_upper = np.flip(mygs.lim_contour[np.logical_and(mygs.lim_contour[:, 0] > 0.2,
                                                         mygs.lim_contour[:, 1] > 0.06), :], axis=0)

    def max_minor_rad(r0, kappa, delta, lim_pad=0.1):
        rmax, rmin = max(lim_upper[:, 0]), min(lim_upper[:, 0])
        a_from_r = min(r0 - rmin, rmax - r0)
        r_zmax = r0 + a_from_r*np.cos(np.pi/2.0 + np.arcsin(delta))
        zmax_lim = np.interp(r_zmax, lim_upper[:, 0], lim_upper[:, 1]) - lim_pad
        return min(zmax_lim/kappa, a_from_r)

    r0, z0 = RCENTR, 0.05
    ku, du, kl, dl = 1.5, 0.7, 1.8, 0.3
    a = max_minor_rad(r0, 0.5*(ku + kl), 0.5*(du + dl))
    xpoints = xpoints_from_moments(r0, z0, a, ku, du, kl, dl)
    isoflux = create_isoflux_xpts(80, r0, z0, a, ku, du, kl, dl, r_inner_mid=None,
                                  zeta_outer_upper=0.0, zeta_outer_lower=-0.1,
                                  zeta_inner_upper=-0.1, zeta_inner_lower=-0.15)
    mygs.set_isoflux_constraints(isoflux)
    mygs.set_saddle_constraints(xpoints)


def run_case(name, maxits=None, save_as=None):
    from OpenFUSIONToolkit import OFT_env
    from OpenFUSIONToolkit.TokaMaker import TokaMaker
    from OpenFUSIONToolkit.TokaMaker.meshing import load_gs_mesh
    from OpenFUSIONToolkit.TokaMaker.util import create_power_flux_fun, read_eqdsk

    common.fetch_inputs()
    mod = CASES[name]
    inp = {**PUBLISHED, **mod}
    eqdsk = read_eqdsk(common.DIIID_GEQDSK)
    baseline = inp.get('baseline', False)

    myOFT = OFT_env(nthreads=common.NTHREADS)
    mygs = TokaMaker(myOFT)
    if baseline:
        mygs.settings.maxits = 200      # as in DIIID_baseline_ex.ipynb
    if maxits:
        mygs.settings.maxits = maxits   # diagnostic override (the shape-helper notebook uses the default, 40)
    mesh_pts, mesh_lc, mesh_reg, coil_dict, cond_dict = load_gs_mesh(common.DIIID_MESH)
    mygs.setup_mesh(mesh_pts, mesh_lc, mesh_reg)
    mygs.setup_regions(cond_dict=cond_dict, coil_dict=coil_dict)
    F0 = eqdsk['rcentr']*eqdsk['bcentr'] if baseline else RCENTR*BCENTR
    mygs.setup(order=2, F0=F0)
    mygs.set_coil_vsc({'F9A': 1.0, 'F9B': -1.0})

    if baseline:
        Ip, pax = eqdsk['ip'], eqdsk['pres'][0]
        mygs.set_targets(Ip=Ip, pax=pax)
        mygs.set_isoflux_constraints(eqdsk['rzout'].copy())
        targets, w_e = EFIT_COIL_TARGETS, 61.0
    else:
        Ip, pax = inp['Ip'], inp['pax']
        mygs.set_targets(Ip=Ip, pax=pax)
        shapehelper_shape(mygs)
        targets, w_e = {k: 0.0 for k in mygs.coil_sets}, 1.E-2
    reg = []
    for cname, tgt in targets.items():
        if cname.startswith('ECOIL'):
            reg.append(mygs.coil_reg_term({cname: 1.0}, target=tgt, weight=w_e))
        elif cname.startswith('F5'):
            reg.append(mygs.coil_reg_term({cname: 1.0}, target=tgt, weight=1.E2))
        elif cname.startswith('F'):
            reg.append(mygs.coil_reg_term({cname: 1.0}, target=tgt, weight=1.E0))
    reg.append(mygs.coil_reg_term({'#VSC': 1.0}, target=0.0, weight=1.E-2))
    mygs.set_coil_reg(reg_terms=reg)

    if baseline or inp['profiles'] == 'efit':
        ffp_prof, pp_prof = efit_profile_shapes(eqdsk)
    else:
        ffp_prof = create_power_flux_fun(40, *inp['ffp'])
        pp_prof = create_power_flux_fun(40, *inp['pp'])
    mygs.set_profiles(ffp_prof=ffp_prof, pp_prof=pp_prof)

    out = save_as or name
    rec = {'case': out, 'scan_case': name, 'machine': 'DIII-D', 'inputs': {k: v for k, v in inp.items()},
           'Ip_target': float(Ip), 'pax_target': float(pax), 'maxits': int(mygs.settings.maxits)}
    mygs.init_psi(1.7, 0.0, 0.5, 1.3, -0.4)
    try:
        _, its = mygs.solve(return_its=True)
    except ValueError as err:
        rec.update({'converged': False, 'error': str(err)})
        common.save_case(out, rec, None)
        print(f'{out}: NOT CONVERGED: {err}')
        return
    eq, arrays = common.equilibrium_record(mygs)
    rec.update({'converged': True, 'nl_its': int(its), 'diverted': bool(mygs.diverted), **eq})
    common.save_case(out, rec, arrays)
    print(f"{out}: q0 {rec['q_0']:.2f} q95 {rec['q_95']:.2f} li {rec['l_i']:.2f} betaN {rec['beta_N']:.2f} "
          f"FF' zero crossings {rec['ffp_zero_crossings']}")


def efit_record():
    """FF', P' and q straight from the gEQDSK, sign-flipped to TokaMaker's convention (as the
    baseline notebook does when plotting). No solve."""
    from OpenFUSIONToolkit.TokaMaker.util import read_eqdsk
    common.fetch_inputs()
    e = read_eqdsk(common.DIIID_GEQDSK)
    psi = np.linspace(0.0, 1.0, e['nr'])
    ffp, pp = -np.asarray(e['ffprim']), -np.asarray(e['pprime'])
    rec = {'case': 'efit_192185_gfile', 'machine': 'DIII-D', 'source': 'g192185.02440 (EFIT01)',
           'Ip': float(e['ip']), 'P_ax': float(e['pres'][0]), 'q_0': float(e['qpsi'][0]),
           'q_95': float(np.interp(0.95, psi, e['qpsi'])), 'q_min': float(np.min(e['qpsi'])),
           'psi_axis_gfile': float(e['psimag']), 'psi_lcfs_gfile': float(e['psibry']),
           **common.profile_metrics(psi, ffp, pp)}
    arrays = {'psi': psi, 'F': np.asarray(e['fpol']), 'FFp': ffp, 'P': np.asarray(e['pres']),
              'Pp': pp, 'psi_q': psi, 'q': np.asarray(e['qpsi'])}
    common.save_case('efit_192185_gfile', rec, arrays)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--case', choices=list(CASES))
    ap.add_argument('--maxits', type=int, help='override the nonlinear iteration limit (diagnostic)')
    ap.add_argument('--save-as', help='record name for a --maxits rerun (default: the case name)')
    args = ap.parse_args()
    if args.case:
        run_case(args.case, args.maxits, args.save_as)
        return
    efit_record()
    failed = common.run_in_subprocesses(os.path.abspath(__file__), list(CASES))
    if failed:
        sys.exit(f'{len(failed)} case(s) crashed: {failed}')


if __name__ == '__main__':
    main()
