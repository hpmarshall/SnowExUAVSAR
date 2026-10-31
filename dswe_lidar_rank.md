# dSWE vs lidar: correlation matrix and pair ranking

Pearson r at 90 m between each pair's dSWE retrieval (mean over HH/HV/VH/VV;
per-pol values in dswe_lidar_rank.csv) and BOTH QSI lidar snow depths.
Rank criterion: r vs the season-matched lidar, among pairs passing the gates
(SNOTEL interval dSWE >= +10 mm AND n90 >= 500).

| pair | line | dates | SNOTEL (mm) | r vs 2020 lidar | r vs 2021 lidar | n90 | cross-pol | gates | rank |
|---|---|---|---|---|---|---|---|---|---|
| 0 | 05208 | 2021-03-16->2021-03-22 | +27.9 | +0.074 | +0.174 | 3158 | +0.852 | pass  | 11 |
| 1 | 23205 | 2021-03-16->2021-03-22 | +27.9 | +0.388 | +0.612 | 3595 | +0.915 | pass  | 3 |
| 2 | 23205 | 2021-03-10->2021-03-22 | +27.9 | +0.385 | +0.638 | 3666 | +0.952 | pass  | 2 |
| 3 | 05208 | 2021-03-10->2021-03-16 | +0.0 | +0.025 | -0.164 | 4415 | +0.957 | no snowfall  |  |
| 4 | 23205 | 2021-03-10->2021-03-16 | +0.0 | +0.342 | +0.557 | 3813 | +0.966 | no snowfall  |  |
| 5 | 23205 | 2021-03-03->2021-03-16 | +12.7 | +0.293 | +0.525 | 3715 | +0.938 | pass  | 4 |
| 6 | 05208 | 2021-03-03->2021-03-10 | +12.7 | +0.320 | +0.476 | 4413 | +0.926 | pass  | 5 |
| 7 | 23205 | 2021-03-03->2021-03-10 | +12.7 | -0.026 | +0.204 | 3668 | +0.928 | pass  | 9 |
| 8 | 05208 | 2021-02-10->2021-03-03 | +172.8 | +0.665 | +0.668 | 3220 | +0.901 | pass  | 1 |
| 9 | 23205 | 2021-02-10->2021-03-03 | +172.8 | +0.044 | +0.278 | 3717 | +0.857 | pass custom unwrap | 6 |
| 10 | 05208 | 2021-02-03->2021-02-10 | +15.2 | +0.252 | -0.062 | 4394 | +0.961 | pass  | 16 |
| 11 | 23205 | 2021-02-03->2021-02-10 | +15.2 | +0.174 | +0.052 | 3160 | +0.916 | pass  | 14 |
| 12 | 23205 | 2021-01-27->2021-02-10 | +58.4 | +0.090 | -0.011 | 3410 | +0.412 | pass custom unwrap | 15 |
| 13 | 23205 | 2021-01-27->2021-02-03 | +43.2 | +0.193 | +0.254 | 652 | +0.881 | pass  | 8 |
| 14 | 23205 | 2021-01-20->2021-02-03 | +61.0 | -0.153 | -0.097 | 3181 | +0.735 | pass custom unwrap | 17 |
| 15 | 23205 | 2021-01-20->2021-01-27 | +17.8 | +0.103 | +0.052 | 3708 | +0.766 | pass  | 13 |
| 16 | 23205 | 2021-01-15->2021-01-27 | +17.8 | -0.320 | -0.499 | 3198 | +0.925 | pass  | 19 |
| 17 | 23205 | 2021-01-15->2021-01-20 | +0.0 | -0.334 | -0.528 | 3683 | +0.856 | no snowfall  |  |
| 18 | 23205 | 2020-02-21->2020-03-11 | +30.5 | +0.269 | +0.590 | 2796 | +0.956 | pass  | 7 |
| 19 | 23205 | 2020-02-13->2020-02-21 | +33.0 | +0.140 | +0.371 | 2665 | +0.781 | pass  | 12 |
| 20 | 23205 | 2020-01-31->2020-02-13 | +63.5 | -0.188 | -0.404 | 1660 | +0.589 | pass custom unwrap | 18 |
| 21 | 23205 | 2019-12-20->2020-01-31 | +330.2 | +0.190 | +0.394 | 852 | +0.397 | pass custom unwrap | 10 |

## Top 5 pairs for snow-model precipitation-pattern input

1. **pair 8** (05208 2021-02-10->2021-03-03): r_lidar +0.668, SNOTEL +173 mm, cross-pol +0.90, n90 3220
2. **pair 2** (23205 2021-03-10->2021-03-22): r_lidar +0.638, SNOTEL +28 mm, cross-pol +0.95, n90 3666
3. **pair 1** (23205 2021-03-16->2021-03-22): r_lidar +0.612, SNOTEL +28 mm, cross-pol +0.92, n90 3595
4. **pair 5** (23205 2021-03-03->2021-03-16): r_lidar +0.525, SNOTEL +13 mm, cross-pol +0.94, n90 3715
5. **pair 6** (05208 2021-03-03->2021-03-10): r_lidar +0.476, SNOTEL +13 mm, cross-pol +0.93, n90 4413

## Zero-change example (outside the ranked list)

**pair 4** (23205 2021-03-10->2021-03-16): r vs matched lidar +0.557 with SNOTEL +0.0 mm and n90 3813. A minimal-dSWE
interval bracketing the 2021 lidar whose retrieval still tracks the snowpack
pattern: a useful model test case (near-zero precip input; pattern from
redistribution/settling), excluded from the ranking only by the snowfall gate.

Custom-unwrap caveat (pairs 9/12/14/20/21): SNAPHU connected components carry
independent phase constants, which depresses scene-wide correlations for the
multi-component pairs (esp. 20, 21).

