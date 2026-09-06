"""Correlate every pair's dSWE retrieval with both QSI lidar snow-depth maps, and
rank the pairs most useful for defining a snow model's sub-km precipitation pattern.

For each pair (22) x pol (4): Pearson r at 90 m between the dSWE map and BOTH lidar
snow depths (2020-02-09 and 2021-03-15). Pearson is shift-invariant, so the two
unanchored pairs (20, 21) are comparable; the custom-unwrapped pairs' independent
per-component phase constants are a real caveat and are footnoted.

Ranking (fully shown in the table so it is auditable):
  r_lidar  = mean over 4 pols of r vs the season-matched lidar
  gates    = SNOTEL interval dSWE >= +10 mm (an interval with no snowfall carries no
             precipitation-pattern information) and n90 >= 500 valid cells
  rank     = r_lidar among gated pairs; cross-pol consistency, coverage, SNOTEL dSWE
             and baseline reported as tiebreaker context.
Plus one "zero-change example" outside the ranked list (user request): pair 4
(2021-03-10 -> 03-16, SNOTEL +0.0 mm, brackets the 2021 lidar) -- a minimal-dSWE
interval whose retrieval still tracks the snowpack pattern is a useful model test
case even though it fails the snowfall gate.

Outputs: dswe_lidar_rank.md (summary + top 5) and dswe_lidar_rank.csv (per-pol).

Run:  PROJ_DATA=/opt/anaconda3/envs/myenv/share/proj \
      /opt/anaconda3/envs/myenv/bin/python rank_dswe_pairs.py
"""

import itertools               # cross-pol combinations
import numpy as np
import pandas as pd
import xarray as xr

from compute_phase_sd_correlations import block_mean, pearson     # shared 90 m helpers

ZARR_PATH = "ZARR/mores_creek_summit.zarr"
POLS = ["HH", "HV", "VH", "VV"]
LIDARS = {"2020": "2020-02-09", "2021": "2021-03-15"}
GATE_SNOTEL_M = 0.010          # >= +10 mm SNOTEL dSWE: interval must contain snowfall
GATE_N90 = 500                 # >= 500 valid 90 m cells: statistical support
CUSTOM_UNWRAPPED = [9, 12, 14, 20, 21]                            # SNAPHU pipeline provenance
ZERO_CHANGE_EXAMPLE = 4        # user-requested extra: minimal-dSWE bracketing pair
# letter labels for the selected pairs, used consistently across the SNOTEL figure,
# the report tables, the subset NetCDF's rank_note, and the subset GUI
PAIR_LETTERS = {8: "A", 2: "B", 1: "C", 5: "D", 6: "E", 4: "F"}

if __name__ == "__main__":
    ds = xr.open_zarr(ZARR_PATH)
    sno = pd.read_csv("snotel_mcs_swe.csv", index_col="date", parse_dates=True)
    sd90 = {k: block_mean(ds["sd"].sel(sd_time=d).squeeze().values) for k, d in LIDARS.items()}

    rows = []
    for pi in range(ds.sizes["pair"]):
        line = str(ds["line_id"][pi].values)
        d1, d2 = str(ds["pair_date1"][pi].values)[:10], str(ds["pair_date2"][pi].values)[:10]
        season = d2[:4]                                           # matched lidar by winter
        dswe_sno = (sno.loc[d2, "swe_mm"] - sno.loc[d1, "swe_mm"]) / 1000.0
        maps = [block_mean(ds["dswe"][pi, p].values) for p in range(4)]
        # cross-pol consistency (same metric as dswe_validation section 2)
        xp = [pearson(maps[a], maps[b])[0] for a, b in itertools.combinations(range(4), 2)]
        xp = [v for v in xp if np.isfinite(v)]
        for pol_i, pol in enumerate(POLS):
            r20, n20 = pearson(sd90["2020"], maps[pol_i])
            r21, n21 = pearson(sd90["2021"], maps[pol_i])
            rows.append(dict(pair=pi, line=line, date1=d1, date2=d2, pol=pol, season=season,
                             snotel_dswe_mm=round(1000 * dswe_sno, 1),
                             r_lidar2020=r20, n2020=n20, r_lidar2021=r21, n2021=n21,
                             crosspol_mean=np.mean(xp) if xp else np.nan,
                             custom_unwrap=pi in CUSTOM_UNWRAPPED))
    df = pd.DataFrame(rows).round({"r_lidar2020": 3, "r_lidar2021": 3, "crosspol_mean": 3})
    df.to_csv("dswe_lidar_rank.csv", index=False)

    # verification: must reproduce dswe_validation.md's pair-4-HH-vs-2021 value (+0.507)
    chk = df[(df.pair == 4) & (df.pol == "HH")].r_lidar2021.iloc[0]
    assert abs(chk - 0.507) < 0.002, f"cross-check vs dswe_validation failed: {chk}"

    # per-pair summary: mean over pols vs each lidar; matched-season stats drive the rank
    summ = []
    for pi, sub in df.groupby("pair"):
        r0 = sub.iloc[0]
        matched = "r_lidar2021" if r0.season == "2021" else "r_lidar2020"
        n_matched = "n2021" if r0.season == "2021" else "n2020"
        r_lidar = sub[matched].mean()
        n90 = int(sub[n_matched].max())
        gated = (r0.snotel_dswe_mm >= 1000 * GATE_SNOTEL_M) and (n90 >= GATE_N90)
        summ.append(dict(pair=pi, line=r0.line, dates=f"{r0.date1}->{r0.date2}",
                         snotel=r0.snotel_dswe_mm,
                         r2020=sub.r_lidar2020.mean(), r2021=sub.r_lidar2021.mean(),
                         r_lidar=r_lidar, n90=n90, crosspol=r0.crosspol_mean,
                         gated=gated, custom=r0.custom_unwrap))
    S = pd.DataFrame(summ).round({"r2020": 3, "r2021": 3, "r_lidar": 3, "crosspol": 3})
    ranked = S[S.gated].sort_values("r_lidar", ascending=False)
    top5 = ranked.head(5)

    md = ["# dSWE vs lidar: correlation matrix and pair ranking", "",
          "Pearson r at 90 m between each pair's dSWE retrieval (mean over HH/HV/VH/VV;",
          "per-pol values in dswe_lidar_rank.csv) and BOTH QSI lidar snow depths.",
          "Rank criterion: r vs the season-matched lidar, among pairs passing the gates",
          f"(SNOTEL interval dSWE >= +{1000*GATE_SNOTEL_M:.0f} mm AND n90 >= {GATE_N90}).", "",
          "| pair | line | dates | SNOTEL (mm) | r vs 2020 lidar | r vs 2021 lidar | n90 | cross-pol | gates | rank |",
          "|---|---|---|---|---|---|---|---|---|---|"]
    rank_of = {int(r.pair): i + 1 for i, (_, r) in enumerate(ranked.iterrows())}
    for _, r in S.iterrows():
        note = "custom unwrap" if r.custom else ""
        gate_txt = "pass" if r.gated else ("no snowfall" if r.snotel < 1000 * GATE_SNOTEL_M else "low n")
        rk = str(rank_of.get(int(r.pair), ""))
        md.append(f"| {int(r.pair)} | {r.line} | {r.dates} | {r.snotel:+.1f} | {r.r2020:+.3f} "
                  f"| {r.r2021:+.3f} | {r.n90} | {r.crosspol:+.3f} | {gate_txt} {note} | {rk} |")

    md += ["", "## Top 5 pairs for snow-model precipitation-pattern input", ""]
    for i, (_, r) in enumerate(top5.iterrows(), 1):
        md.append(f"{i}. **pair {int(r.pair)}** ({r.line} {r.dates}): r_lidar {r.r_lidar:+.3f}, "
                  f"SNOTEL {r.snotel:+.0f} mm, cross-pol {r.crosspol:+.2f}, n90 {r.n90}"
                  + (" *(custom unwrap)*" if r.custom else ""))
    z = S[S.pair == ZERO_CHANGE_EXAMPLE].iloc[0]
    md += ["", "## Zero-change example (outside the ranked list)", "",
           f"**pair {ZERO_CHANGE_EXAMPLE}** ({z.line} {z.dates}): r vs matched lidar "
           f"{z.r_lidar:+.3f} with SNOTEL {z.snotel:+.1f} mm and n90 {z.n90}. A minimal-dSWE",
           "interval bracketing the 2021 lidar whose retrieval still tracks the snowpack",
           "pattern: a useful model test case (near-zero precip input; pattern from",
           "redistribution/settling), excluded from the ranking only by the snowfall gate.", "",
           "Custom-unwrap caveat (pairs 9/12/14/20/21): SNAPHU connected components carry",
           "independent phase constants, which depresses scene-wide correlations for the",
           "multi-component pairs (esp. 20, 21).", ""]
    open("dswe_lidar_rank.md", "w").write("\n".join(md) + "\n")
    print("\n".join(md[md.index("## Top 5 pairs for snow-model precipitation-pattern input") - 1:]))
    print("wrote dswe_lidar_rank.md and dswe_lidar_rank.csv")
