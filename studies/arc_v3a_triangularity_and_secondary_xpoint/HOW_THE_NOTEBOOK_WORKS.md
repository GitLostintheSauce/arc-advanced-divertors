# How the shape-testing notebook builds the equilibrium

Notebook: `2026-06-08_ARC_diverted_shape_testing.ipynb` by F. Sheehan (Columbia), shared for this study
and published here with permission. The notebook itself is not in this repo;
[`studies/arc_v3a_triangularity_and_secondary_xpoint/code/run_equilibria.py`](../../studies/arc_v3a_triangularity_and_secondary_xpoint/code/run_equilibria.py) reproduces
its equilibrium logic exactly. Cell numbers below are execution counts as saved in the notebook.

## The big picture

TokaMaker solves the **Grad–Shafranov equation**: force balance between plasma pressure and magnetic
forces in an axisymmetric torus. "Free boundary" means the plasma edge is not given. It comes out
of the solution, set by the plasma current *and* the currents in the external coils.

You don't give TokaMaker coil currents. You give it **wishes about the field** (this point should be
on the boundary, this point should be an X-point), and it finds the coil currents that best satisfy
those wishes while also solving force balance. This kind of solve is often called "inverse", because
it works backwards from the desired shape to the coils.

## Where each piece lives

| What | Where | What it does |
|---|---|---|
| Machine, mesh, coils | cell 2 | Loads `2026-07-14_ARC_V3A_public_mesh.h5`, which is byte-identical to the public OFT ARC V3A mesh (PR #298). Sets F0 = R0·B0 = 4.62 × 11.4 T·m (the toroidal field). Makes PF5U/PF5L a vertical-stability pair (`set_coil_vsc`). |
| Coil-current bounds | cell 4 | `set_coil_bounds`: every coil limited to ±1.0×10⁶ **A per turn**. These are *hard* limits the solver must obey. |
| Coil regularization | cell 4 (defined), cell 9 (applied) | `build_min_norm_coil_reg`: a gentle pull of every coil toward 0 A (weight 0.01–0.02), a strong penalty on any difference between upper and lower coils (weight 100, forcing up-down symmetry), and a strong penalty on the vertical-stability "virtual coil" (weight 100). |
| Plasma profiles | cell 5 | `create_power_flux_fun(40, 1.5, 2.0)` for FF′ and `(40, 4.0, 1.0)` for p′. These set the **shape** of the current and pressure profiles across the plasma. Their amplitudes are scaled to hit the targets below. |
| Global targets | cells 6, 9 | Ip = 12 MA and on-axis pressure 0.5 MPa (`set_targets`). These are matched essentially exactly. |
| Shape targets | cell 7 | `a, kappa, delta` → `create_isoflux_xpts(...)` makes 50 points on an ideal D-shaped boundary. |
| X-point targets | cell 7 | Typed in by hand: primary (3.9, ±2.2), secondary (5.0, ±3.8). Saddle weights 100 and 1000. `x_points_new` (from `xpoints_from_moments`) is **only plotted**, never used as a constraint. |
| Strike-point targets | cell 8 | Four points on the wall: (3.51, ±2.36) inner, (4.61, ±3.87) outer. Added as isoflux points, weight 1. |
| Trimming | cell 8 | Drops shape points with \|Z\| above 85 % of the maximum, i.e. the 15 points near the X-point corners. 35 of 50 remain. |
| X-points on the boundary | cell 8 | The four X-points are *also* added as isoflux points (weights 10 primary, 3 secondary), so they are pulled onto the separatrix. |
| Solve | cell 9 | `init_psi` (initial guess), then `solve()`. Afterwards it prints stats, X-point misses with ψ̂, strike points and coil currents. |
| Strike-point picker, export | cells 10–11 | Interactive helper and export for a TORAX config. Not needed for the equilibrium. |

## The two kinds of shape constraint

**Isoflux constraint** (`set_isoflux_constraints`): "these points should all lie on the *same* flux
surface". TokaMaker uses them as boundary points, so they should end up on the separatrix (ψ̂ = 1).
It says nothing about the field direction at those points.

**Saddle constraint** (`set_saddle_constraints`): "the poloidal field should be *zero* here", which
makes the point an X-point (a null). It says nothing about *which* flux surface the null sits on.
A null can sit inside the plasma, on the separatrix, or out in the scrape-off layer. That is why the
notebook adds the X-points as isoflux points too. The comment in cell 8 of the notebook explains this well.

Both kinds are **soft**: they are combined into one weighted least-squares problem. A higher weight
means "this matters more". If two wishes conflict, the solver compromises, leaning toward the
higher-weight one. Only the coil bounds and the Ip/pressure targets are **hard**.

## How the solve works, roughly

1. **Initial guess.** Start from a guess for the flux ψ: a uniform current inside the shape from
   `init_psi`.
2. **Coil fit.** Holding the plasma current fixed, find coil currents that best satisfy all the
   isoflux and saddle wishes, plus the regularization, inside the bounds. This is a small linear
   least-squares problem.
3. **Plasma update.** With those coils, recompute the plasma current from the profiles (p′, FF′),
   scaled to hit Ip and the axis pressure.
4. **Repeat** steps 2–3 until ψ stops changing (tolerance 10⁻⁶). Here that takes 22–25 iterations.

`get_coil_currents()` returns the final current **per turn** in each coil, in A. Multiply by the coil's
turn count from the mesh to get total ampere-turns.

## Things to watch in this notebook

- **The bound comment is wrong.** Cell 4 says "these coils have 100 turns", but the mesh gives
  169–1150 turns per coil. So ±10⁶ A per turn means ±169 to ±1150 MA-turns per coil. These bounds are
  exploratory, not ARC engineering limits.
- **δ does not move the X-points.** `delta` changes the 50 shape points, not the hand-typed X-point
  targets. After trimming, the shape points that would carry most of the δ change are removed.
- **Version difference.** The saved output was made with an OFT v26.6 development build
  (`tmtx-fsa-geometry`) on 4 threads. Re-running with OFT 26.9 (PyPI) reproduces it exactly.
- **Unused imports.** `omfit_classes` is imported but never used, and `ipywidgets` is only used for the
  picker. Neither is installed in this repo's `.venv`; the script skips both.
- **Isoflux reference point.** TokaMaker measures every isoflux point against the *first* one. Here,
  after trimming, the first point is (4.63, −1.73), on the lower shoulder. F. Sheehan's public
  TokaMaker+TORAX project notes that an off-midplane reference point gives an up-down asymmetric
  solve, and rolls the boundary array so the reference is at the outboard midplane. That may explain
  the small up/down differences seen here (ψ̂ 0.9999 vs 1.0011 at the two primary nulls). Not yet tested.
