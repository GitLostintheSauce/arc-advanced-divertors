"""ITER H-mode example: reproduce it and record the solved FF' and P'.

The setup copies the physics cells of OFT's ITER_Hmode_ex.ipynb at commit 8fc4afe (mesh, F0,
VSC, coil bounds and regularization, Ip/pax targets, isoflux and X-point targets, Hmode_profiles
kinetic profiles, power-law inductive j_phi, init_psi + solve, then solve_with_bootstrap twice:
isolate_edge_jBS=False, then True starting from that result). Plot cells are left out and
diagnostic_plots is False (it only draws figures).

The notebook calls solve_with_bootstrap from OpenFUSIONToolkit.TokaMaker.bootstrap. Two versions:
  wheel        the pinned OpenFUSIONToolkit 26.9 wheel (no PR #316)
  oft_8fc4afe  bootstrap.py from OFT commit 8fc4afe, which includes PR #316 ("Fix bootstrap
               profile-derivative accuracy and profile/equilibrium grid mismatch"). Loaded as a
               standalone module next to the installed package; no pins change.

Each version writes two records: iter_hmode_<version>_no_isolate and iter_hmode_<version>_isolate.

Usage (from the repo root):
    .venv/bin/python studies/diiid_iter_flux_function_realism/code/run_iter_hmode.py            # both versions
    .venv/bin/python studies/diiid_iter_flux_function_realism/code/run_iter_hmode.py --case wheel
"""
import argparse
import importlib.util
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common  # noqa: E402

CASES = ('wheel', 'oft_8fc4afe')
IP_TARGET, P0_TARGET = 13.0E6, 6.2E5
N_SAMPLE = 257


def load_bootstrap(version):
    if version == 'wheel':
        from OpenFUSIONToolkit.TokaMaker import bootstrap
        return bootstrap
    import OpenFUSIONToolkit.TokaMaker  # noqa: F401  (parent package for the relative imports)
    spec = importlib.util.spec_from_file_location('OpenFUSIONToolkit.TokaMaker.bootstrap_8fc4afe',
                                                  common.BOOTSTRAP_8FC4AFE)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def kinetic_profiles(bs):
    """ne, Te, ni, Ti exactly as in the notebook."""
    xphalf, widthp_Te, widthp_ne = 0.965, 0.1, 0.35
    ne = bs.Hmode_profiles(edge=0.35, ped=0.6, core=1.1, rgrid=N_SAMPLE, expin=1.6, expout=1.6,
                           widthp=widthp_ne, xphalf=xphalf)*1e20
    Te = bs.Hmode_profiles(edge=1500., ped=5000., core=21000., rgrid=N_SAMPLE, expin=1.3, expout=1.7,
                           widthp=widthp_Te, xphalf=xphalf)
    return ne, Te, ne.copy(), Te.copy()


def run_case(version):
    import matplotlib
    matplotlib.use('Agg')
    from OpenFUSIONToolkit import OFT_env
    from OpenFUSIONToolkit.TokaMaker import TokaMaker
    from OpenFUSIONToolkit.TokaMaker.meshing import load_gs_mesh
    from OpenFUSIONToolkit.TokaMaker.util import create_power_flux_fun

    common.fetch_inputs()
    bs = load_bootstrap(version)
    myOFT = OFT_env(nthreads=common.NTHREADS)
    mygs = TokaMaker(myOFT)
    mesh_pts, mesh_lc, mesh_reg, coil_dict, cond_dict = load_gs_mesh(common.ITER_MESH)
    mygs.setup_mesh(mesh_pts, mesh_lc, mesh_reg)
    mygs.setup_regions(cond_dict=cond_dict, coil_dict=coil_dict)
    mygs.settings.maxits = 100
    mygs.setup(order=2, F0=5.3*6.2)
    mygs.set_coil_vsc({'VS': 1.0})
    mygs.set_coil_bounds({key: [-50.E6, 50.E6] for key in mygs.coil_sets})
    mygs.set_targets(Ip=IP_TARGET, pax=P0_TARGET)
    isoflux_pts = np.array([[8.20, 0.41], [8.06, 1.46], [7.51, 2.62], [6.14, 3.78], [4.51, 3.02],
                            [4.26, 1.33], [4.28, 0.08], [4.49, -1.34], [7.28, -1.89], [8.00, -0.68]])
    x_point = np.array([[5.125, -3.4]])
    mygs.set_isoflux_constraints(np.vstack((isoflux_pts, x_point)))
    mygs.set_saddle_constraints(x_point)
    reg = []
    for name in mygs.coil_sets:
        w = 2.E-2 if name.startswith('CS1') else 1.E-2
        if name.startswith(('CS', 'PF', 'VS')):
            reg.append(mygs.coil_reg_term({name: 1.0}, target=0.0, weight=w))
    reg.append(mygs.coil_reg_term({'#VSC': 1.0}, target=0.0, weight=1.E2))
    mygs.set_coil_reg(reg_terms=reg)

    ne, Te, ni, Ti = kinetic_profiles(bs)
    psi_sample = np.linspace(0.0, 1.0, N_SAMPLE)
    pressure_in = 1.602e-19*(ne*Te + ni*Ti)
    inductive_jphi = create_power_flux_fun(N_SAMPLE, 2.25, 2.5)['y']

    mygs.init_psi(6.3, 0.5, 2.0, 1.4, 0.0)
    mygs.solve()
    Zeff = np.linspace(1.5, 1.5, N_SAMPLE)
    mygs.settings.pm = False
    mygs.update_settings()

    for isolate in (False, True):
        name = f"iter_hmode_{version}_{'isolate' if isolate else 'no_isolate'}"
        res = bs.solve_with_bootstrap(mygs, ne, Te, ni, Ti, Zeff, IP_TARGET, inductive_jphi,
                                      scale_jBS=1.0, isolate_edge_jBS=isolate, diagnostic_plots=False)
        eq, arrays = common.equilibrium_record(mygs)
        rec = {'case': name, 'machine': 'ITER', 'bootstrap_module': version, 'isolate_edge_jBS': isolate,
               'Ip_target': IP_TARGET, 'pax_target': P0_TARGET, 'converged': True,
               'diverted': bool(mygs.diverted), **eq}
        arrays.update({'psi_input': psi_sample, 'ne': ne, 'Te': Te, 'pressure_input': pressure_in,
                       'jphi_total': np.asarray(res['total_j_phi']), 'j_BS': np.asarray(res['j_BS']),
                       'j_inductive': np.asarray(res['j_inductive'])})
        if res.get('isolated_j_BS') is not None:
            arrays['isolated_j_BS'] = np.asarray(res['isolated_j_BS'])
        common.save_case(name, rec, arrays)
        print(f"{name}: q0 {rec['q_0']:.2f} q95 {rec['q_95']:.2f} li {rec['l_i']:.2f} betaN {rec['beta_N']:.2f} "
              f"FF' zero crossings {[round(x, 3) for x in rec['ffp_zero_crossings']]} "
              f"P' argmax {rec['pp_argmax_psi']:.3f}", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--case', choices=CASES)
    args = ap.parse_args()
    if args.case:
        run_case(args.case)
        return
    failed = common.run_in_subprocesses(os.path.abspath(__file__), CASES)
    if failed:
        sys.exit(f'{len(failed)} case(s) crashed: {failed}')


if __name__ == '__main__':
    main()
