# Roadmap

**Question:** how hard is it to *form and hold* each advanced divertor configuration on ARC V3A, given its real coils, shielding and pulse?

The four configurations studied:

| Configuration | Defining feature |
|---|---|
| X-point target (XPT) | Second X-point near the outer target |
| Super-X (SXD) | Outer leg pushed to large major radius |
| X-divertor (XD) | Field lines flare at the target plate |
| Snowflake (SFD) | Second-order null with four legs |

The study runs for about 34 weeks, full time.

## Where code lives

| What | Where |
|---|---|
| ARC scans, configs, figures, paper | This repo |
| Reusable tools (divertor metrics, coil limits, ARC pulse example) | [OFT fork](https://github.com/GitLostintheSauce/OpenFUSIONToolkit), one branch per tool, then a small PR to [OpenFUSIONToolkit](https://github.com/OpenFUSIONToolkit/OpenFUSIONToolkit) |
| TORAX | Used as is, pinned to 1.3.1 |

Rule of thumb: if another tokamak could use it, it goes to OFT. If it's about ARC results, it stays here.

## Stages

### 0 · Set up (week 0)
- Create this repo and a single pinned Python environment (OFT 26.9, TORAX 1.3.1).
- Check the install by reproducing the OFT ARC V3A baseline example.

### 1 · Learn (weeks 1–6)
- **Read** the ARC, STEP and DEMO divertor literature.
- **Run the OFT examples:** `fixed_boundary_ex1`, `ITER_baseline_ex`, `DIIID_shapehelper_ex`, `ARCV3A_control_ex`.
- **Reproduce ARC V3A:** the double-null reference equilibrium, built from the public inputs in Zenodo [19498373](https://zenodo.org/records/19498373).
- **Benchmark** against a published result, e.g. a SPARC XPT at reduced Ip, or an ITER XD within PF coil limits.

### 2 · Build (weeks 7–18)
- **Divertor metrics tool** *(OFT PR)*: strike points, flux expansion, connection length L∥, target angle and ΔR_sep. OFT doesn't have this yet, and every scan needs it.
- **Coil limits tool** *(OFT PR)*: PF/CS forces, current density, stored energy and the shield keep-out zone. Reports which limit binds.
- **ARC flattop scans:** map Ip × shield standoff × coil limits for each configuration, and find the binding limit.
- **l_i / β_p sensitivity:** how far strike points move as l_i and β_p change, and the vertical growth rate for each configuration.

### 3 · Pulse (weeks 19–26)
- **Coil-voltage fix:** follow [OpenFUSIONToolkit#349](https://github.com/OpenFUSIONToolkit/OpenFUSIONToolkit/pull/349). With coupled voltage-driven coils, the current error grows as dt shrinks.
- **ARC pulse with TokaMaker + TORAX** *(OFT example)*: include coil voltage and slew-rate limits, and find when each configuration can form and when it is lost.
- **Score the benefit:** feed the stage 2 divertor geometry into TORAX's extended Lengyel divertor model.
- **Compare with one or two other power plant designs** (e.g. MANTA, which already has an OFT example).
- **Talk and abstract:** a group talk, then APS-DPP.

### 4 · Ship (weeks 27–34)
- **Final figures:** coil cost vs divertor benefit across the pulse, plus the levers that change it (trim coils, lower Ip, thinner shield).
- **Paper:** Nuclear Fusion or PPCF.
- **OFT example:** release the ARC advanced-divertor workflow next to the other OFT `AdvancedWorkflows`, and cite it in the paper.

## TORAX notes

- The OFT TokaMaker+TORAX workflow pins `torax==1.3.1`, so nothing blocks starting.
- Don't upgrade TORAX casually: the workflow imports TORAX internals, and TORAX main has already moved its divertor code.
- Watch [google-deepmind/torax#2383](https://github.com/google-deepmind/torax/pull/2383), a direct TokaMaker → TORAX link that wouldn't need gEQDSK files.
- STEP is spherical, so the default QLKNN transport model is outside its training range there. STEP would need a different transport model.
