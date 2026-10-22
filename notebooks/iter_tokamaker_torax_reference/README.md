# ITER TokaMaker+TORAX pulse example: reference run under TORAX 1.4.0

`ITER_TokaMaker_TORAX_pulse_ex_reference_executed.ipynb` is the OpenFUSIONToolkit example notebook
`src/examples/TokaMaker/AdvancedWorkflows/Pulse_Design/TokaMaker_TORAX/ITER/ITER_TokaMaker_TORAX_pulse_ex.ipynb`
(OFT commit 8fc4afe, LGPL-3.0), executed with its code cells unchanged and all stored outputs cleared beforehand.
It ran on 2026-10-21 on Apple M1 (macOS, CPU JAX, float64) with `OFT_env(nthreads=4)`.

## Environment
This notebook does **not** run with this repository's pinned environment. With the pinned `torax==1.3.1`,
`OpenFUSIONToolkit.TokaMaker.pulse_design` (OFT 26.9) fails to import, because it uses `ConfinementMode`, which
first appeared in TORAX 1.4.0.

The run used a separate experimental environment: `requirements-lock.txt` with four changes, plus pytest for tests.

| package | lock | used here |
|---|---|---|
| torax | 1.3.1 | 1.4.0 |
| jax | 0.9.2 | 0.10.0 |
| jaxlib | 0.9.2 | 0.10.0 |
| fusion_surrogates | 0.4.5 | 0.4.6 |

The full `pip freeze` is in `environment_freeze_torax140.txt`. The project pins have not been changed.

## Settings (as published)
- TORAX `tx_dt = 1` s, 25 TokaMaker solve times.
- `fly(run_name='tmp', output_mode='debug', max_loop=2, t_ave_toggle='flattop', t_ave_window=25, relax=True, relax_duration=5)`.
- ADAPTIVE_SOURCE pedestal with the Martin-scaling L-H formation model and a 10 s transition ramp.

## Result
- Every cell ran.
- Both coupling passes completed without TORAX errors. TokaMaker solved 25/25 equilibria in each pass.
- OFT L-H/H-L times: pass 1 136.4 / 499.4 s, pass 2 168.7 / 489.7 s.
- Flux consumption: 170.55 Wb (pass 1), 165.85 Wb (pass 2).

## Limitations
- **Not converged.** With `relax=True`, the inter-loop relax shifts the absolute boundary flux by about 3-4 Wb/rad
  per coupling pass, which drives CS currents toward their limits over more passes. `relax=False` removes this in
  diagnostic tests.
- **L-H timing is timestep-sensitive.** It is about 30 s later at `tx_dt = 0.5`.
- **Not bit-reproducible with 4 OpenMP threads.** Pass-2 flux varied by about 0.1 % between repeat runs.
- **No movie.** The movie cell printed an ffmpeg warning because ffmpeg was not installed.
- OFT's P_LH plots label the Delabie curve "(used)", but TORAX used the configured Martin scaling.
- The notebook outputs contain local absolute paths from the machine it ran on.
