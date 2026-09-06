# dSWE vs lidar: correlation matrix and pair ranking

Pearson r at 90 m between each pair's dSWE retrieval (mean over HH/HV/VH/VV;
per-pol values in dswe_lidar_rank.csv) and BOTH QSI lidar snow depths.
Rank criterion: r vs the season-matched lidar, among pairs passing the gates
(SNOTEL interval dSWE >= +10 mm AND n90 >= 500).

| pair | line | dates | SNOTEL (mm) | r vs 2020 lidar | r vs 2021 lidar | n90 | cross-pol | gates | rank |
|---|---|---|---|---|---|---|---|---|---|
| 0 | 05208 | 2021-03-16->2021-03-22 | +27.9 | +0.061 | +0.185 | 3157 | +0.820 | pass  | 11 |
| 1 | 23205 | 2021-03-16->2021-03-22 | +27.9 | +0.229 | +0.483 | 3594 | +0.882 | pass  | 3 |
| 2 | 23205 | 2021-03-10->2021-03-22 | +27.9 | +0.362 | +0.577 | 3685 | +0.947 | pass  | 2 |
| 3 | 05208 | 2021-03-10->2021-03-16 | +0.0 | -0.210 | -0.292 | 4415 | +0.948 | no snowfall  |  |
| 4 | 23205 | 2021-03-10->2021-03-16 | +0.0 | +0.370 | +0.486 | 3819 | +0.951 | no snowfall  |  |
| 5 | 23205 | 2021-03-03->2021-03-16 | +12.7 | +0.238 | +0.466 | 3729 | +0.899 | pass  | 4 |
| 6 | 05208 | 2021-03-03->2021-03-10 | +12.7 | +0.335 | +0.454 | 4413 | +0.939 | pass  | 5 |
| 7 | 23205 | 2021-03-03->2021-03-10 | +12.7 | -0.028 | +0.170 | 3679 | +0.911 | pass  | 12 |
| 8 | 05208 | 2021-02-10->2021-03-03 | +172.8 | +0.650 | +0.637 | 3220 | +0.894 | pass  | 1 |
| 9 | 23205 | 2021-02-10->2021-03-03 | +172.8 | +0.257 | +0.318 | 3723 | +0.871 | pass custom unwrap | 6 |
| 10 | 05208 | 2021-02-03->2021-02-10 | +15.2 | +0.123 | -0.142 | 4394 | +0.957 | pass  | 18 |
| 11 | 23205 | 2021-02-03->2021-02-10 | +15.2 | +0.187 | +0.109 | 3157 | +0.936 | pass  | 14 |
| 12 | 23205 | 2021-01-27->2021-02-10 | +58.4 | +0.289 | +0.140 | 3408 | +0.408 | pass custom unwrap | 13 |
| 13 | 23205 | 2021-01-27->2021-02-03 | +43.2 | +0.251 | +0.255 | 653 | +0.892 | pass  | 8 |
| 14 | 23205 | 2021-01-20->2021-02-03 | +61.0 | -0.027 | -0.046 | 3191 | +0.674 | pass custom unwrap | 16 |
| 15 | 23205 | 2021-01-20->2021-01-27 | +17.8 | -0.050 | +0.002 | 3727 | +0.718 | pass  | 15 |
| 16 | 23205 | 2021-01-15->2021-01-27 | +17.8 | -0.402 | -0.467 | 3209 | +0.924 | pass  | 19 |
| 17 | 23205 | 2021-01-15->2021-01-20 | +0.0 | -0.309 | -0.434 | 3695 | +0.823 | no snowfall  |  |
| 18 | 23205 | 2020-02-21->2020-03-11 | +30.5 | +0.313 | +0.550 | 2795 | +0.941 | pass  | 7 |
| 19 | 23205 | 2020-02-13->2020-02-21 | +33.0 | +0.189 | +0.389 | 2658 | +0.706 | pass  | 10 |
| 20 | 23205 | 2020-01-31->2020-02-13 | +63.5 | -0.080 | -0.146 | 1657 | +0.533 | pass custom unwrap | 17 |
| 21 | 23205 | 2019-12-20->2020-01-31 | +330.2 | +0.222 | +0.291 | 858 | +0.349 | pass custom unwrap | 9 |

## Top 5 pairs for snow-model precipitation-pattern input

1. **pair 8** (05208 2021-02-10->2021-03-03): r_lidar +0.637, SNOTEL +173 mm, cross-pol +0.89, n90 3220
2. **pair 2** (23205 2021-03-10->2021-03-22): r_lidar +0.577, SNOTEL +28 mm, cross-pol +0.95, n90 3685
3. **pair 1** (23205 2021-03-16->2021-03-22): r_lidar +0.483, SNOTEL +28 mm, cross-pol +0.88, n90 3594
4. **pair 5** (23205 2021-03-03->2021-03-16): r_lidar +0.466, SNOTEL +13 mm, cross-pol +0.90, n90 3729
5. **pair 6** (05208 2021-03-03->2021-03-10): r_lidar +0.454, SNOTEL +13 mm, cross-pol +0.94, n90 4413

## Zero-change example (outside the ranked list)

**pair 4** (23205 2021-03-10->2021-03-16): r vs matched lidar +0.486 with SNOTEL +0.0 mm and n90 3819. A minimal-dSWE
interval bracketing the 2021 lidar whose retrieval still tracks the snowpack
pattern: a useful model test case (near-zero precip input; pattern from
redistribution/settling), excluded from the ranking only by the snowfall gate.

Custom-unwrap caveat (pairs 9/12/14/20/21): SNAPHU connected components carry
independent phase constants, which depresses scene-wide correlations for the
multi-component pairs (esp. 20, 21).

