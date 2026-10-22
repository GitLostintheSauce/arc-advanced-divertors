# ITER TokaMaker+TORAX pulse example: reference run under TORAX 1.4.0

This is a reference execution of the OpenFUSIONToolkit example notebook
`src/examples/TokaMaker/AdvancedWorkflows/Pulse_Design/TokaMaker_TORAX/ITER/ITER_TokaMaker_TORAX_pulse_ex.ipynb`
(OFT commit 8fc4afe, LGPL-3.0), executed with its code cells unchanged and all stored outputs cleared beforehand.
It ran on 2026-10-21 on Apple M1 (macOS, CPU JAX, float64) with `OFT_env(nthreads=4)`.

The executed notebook (7.1 MB) is attached to the
[`iter-reference-2026-10-21` release](https://github.com/GitLostintheSauce/arc-advanced-divertors/releases/tag/iter-reference-2026-10-21):
[`ITER_TokaMaker_TORAX_pulse_ex_reference_executed.ipynb`](https://github.com/GitLostintheSauce/arc-advanced-divertors/releases/download/iter-reference-2026-10-21/ITER_TokaMaker_TORAX_pulse_ex_reference_executed.ipynb)
(SHA-256 `df90286b51f2f0af6dbba23d2d890812a80f91a8ee1b1a18bc5ac0c1e616dbb5`).
Before upload, local absolute paths in the outputs were replaced with `<repo>` and `<ipykernel>`, and the kernelspec
was set to `python3`. Nothing else was changed.

CI runs the same example on every push (the "ITER TokaMaker+TORAX pulse example" step) and checks its L-H times
and flux consumption against the results below, so this reference also works as a regression check.

## Environment
The run used the environment this repository pins in `requirements-lock.txt` (`torax==1.4.0`, `jax`/`jaxlib==0.10.0`,
`fusion_surrogates==0.4.6`), plus pytest, pluggy and iniconfig, which the notebook does not use.

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
- OFT's P_LH plots label the Delabie curve "(used)", but TORAX used the configured Martin scaling
  (fix proposed in [OFT PR #401](https://github.com/OpenFUSIONToolkit/OpenFUSIONToolkit/pull/401)).
