Columns: **gEQDSK** = gEQDSK (Zenodo, FreeGS); **TM+TX** = TM+TX flattop (as published); **TM+TX + secondary X-pts** = TM+TX flattop + secondary X-point targets. Percentages are differences from the gEQDSK.

| Quantity | gEQDSK | TM+TX | TM+TX + secondary X-pts |
|---|---|---|---|
| Plasma current Ip [MA] | 12.00 | 12.00 (-0 %) | 12.00 (-0 %) |
| Vacuum field at R = 4.62 m [T] | 11.40 | 11.40 (-0 %) | 11.40 (+0 %) |
| Magnetic axis R [m] | 4.700 | 4.669 (-1 %) | 4.658 (-1 %) |
| Geometric centre R_geo [m] | 4.622 | 4.618 (-0 %) | 4.632 (+0 %) |
| Minor radius a [m] | 1.180 | 1.143 (-3 %) | 1.084 (-8 %) |
| Elongation κ | 1.842 | 1.775 (-4 %) | 1.878 (+2 %) |
| Triangularity δ (mean, from boundary) | 0.589 | 0.580 (-2 %) | 0.550 (-7 %) |
| Triangularity δ (from X-points: (R_geo − R_x)/a) | 0.611 | 0.669 (+9 %) | 0.691 (+13 %) |
| Plasma volume [m³] | 194.3 | 183.4 (-6 %) | 167.3 (-14 %) |
| q0 | 0.89 | 0.55 (-38 %) | 0.52 (-42 %) |
| q95 | 3.20 | 3.29 (+3 %) | 3.26 (+2 %) |
| Axis pressure [MPa] | 3.20 | 1.63 (-49 %) | 1.64 (-49 %) |
| Volume-average pressure [MPa] | 0.720 | 0.630 (-13 %) | 0.655 (-9 %) |
| β_N | 1.56 | 1.32 (-15 %) | 1.31 (-16 %) |
| β_p | 0.882 | 0.701 (-21 %) | 0.700 (-21 %) |
| l_i | 1.343 | 0.954 (-29 %) | 0.960 (-28 %) |
| Poloidal flux axis→boundary [Wb/rad] | 8.80 | 8.02 (-9 %) | 7.81 (-11 %) |
| δ upper / lower (from boundary) | 0.589 / 0.589 | 0.651 / 0.509 | 0.650 / 0.451 |
| Upper main X-point (R, Z) [m], ψ̂ there | (3.901, +2.236), 1.000 | (3.853, +2.129), 1.000 | (3.878, +2.212), 1.000 |
| Lower main X-point (R, Z) [m], ψ̂ there | (3.901, -2.236), 1.000 | (3.853, -2.122), 1.016 | (3.888, -2.163), 1.032 |
| Upper outer-divertor null near (5.0, ±3.7) m | (4.93, +3.64), ψ̂ = 0.999 | none (min B_p 350 mT) | (4.98, +3.60), ψ̂ = 1.082 |
| Lower outer-divertor null near (5.0, ±3.7) m | (4.93, -3.64), ψ̂ = 0.999 | none (min B_p 334 mT) | (4.96, -3.57), ψ̂ = 1.104 |

Check of the measuring code against TokaMaker's own `get_stats()` for each TM+TX run:

| Quantity | TM+TX: here / TokaMaker | TM+TX + secondary X-pts: here / TokaMaker |
|---|---|---|
| kappa | 1.775 / 1.788 | 1.878 / 2.020 |
| q95 | 3.290 / 3.290 | 3.260 / 3.261 |
| l_i | 0.954 / 0.965 | 0.960 / 0.977 |
| beta_p | 0.701 / 0.707 | 0.700 / 0.705 |
| beta_N | 1.322 / 1.323 | 1.307 / 1.300 |
| volume_m3 | 183.408 / 183.416 | 167.270 / 168.189 |
| Ip_MA (boundary integral vs header) | 11.906 / 12.000 | 11.855 / 12.000 |

TM+TX run results (end of the coupled run):

| Quantity | TM+TX | TM+TX + secondary X-pts |
|---|---|---|
| Q_fusion | 10.18 | 9.968 |
| P_fusion_MW | 522.1 | 512 |
| H98 | 0.7471 | 0.7175 |
| beta_N_torax | 1.324 | 1.303 |
| q95_torax | 3.362 | 3.349 |
| li3_torax | 0.7391 | 0.7078 |
| l_i_tm (last time) | 0.7461 | 0.72 |
| l_i_tx (last time) | 0.7391 | 0.7078 |
| beta_N_tm (last time) | 1.323 | 1.3 |
| beta_N_tx (last time) | 1.324 | 1.303 |
| Ip_tm_MA (last time) | 12 | 12 |
| Ip_tx_MA (last time) | 12 | 12 |

Coil currents [MA-turns]. The TM+TX notebook bounds every coil at ±50 MA-turns.

| Coil | TM+TX | TM+TX + secondary X-pts | Change |
|---|---|---|---|
| CS1U | -18.25 | -15.27 | +2.98 |
| CS1L | -15.47 | -8.40 | +7.07 |
| CS2U | +24.96 | -31.76 | -56.72 |
| CS2L | +25.53 | -32.01 | -57.54 |
| CS3U | -49.72 **at limit** | +50.00 **at limit** | +99.72 |
| CS3L | -50.00 **at limit** | +48.53 | +98.53 |
| PF1U | +47.44 | +50.00 **at limit** | +2.56 |
| PF1L | +46.78 | +48.44 | +1.66 |
| PF2U | +42.69 | -4.48 | -47.17 |
| PF2L | +41.10 | -4.51 | -45.61 |
| PF3U | -6.79 | +8.20 | +14.99 |
| PF3L | -6.37 | +7.80 | +14.18 |
| PF4U | +1.11 | -10.79 | -11.90 |
| PF4L | +1.05 | -9.89 | -10.94 |
| PF5U | -1.97 | -0.72 | +1.25 |
| PF5L | -2.17 | -1.19 | +0.98 |
| Σ\|I\| | 381.4 | 332.0 | |
