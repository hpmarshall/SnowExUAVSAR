"""Generate the data-table fragments the two LaTeX deliverables \\input.

Reads live sources (the Zarr store, dswe_lidar_rank.csv, dswe_retrieval_log.csv) and
writes tex fragments so the PDFs stay regenerable when numbers change:
  tex_fragments/zarr_dims.tex        dimension/coordinate table (user guide)
  tex_fragments/zarr_vars.tex        data-variable table (user guide)
  tex_fragments/zarr_counts.tex      per-variable 2-D slice counts (user guide)
  tex_fragments/rank_top.tex         top-5 + zero-change table (summary report)
  tex_fragments/rank_appendix.tex    all-22-pairs appendix table (summary report)

The prose lives directly in mcs_zarr_userguide.tex / mcs_summary_report.tex.

Run:  PROJ_DATA=/opt/anaconda3/envs/myenv/share/proj \
      /opt/anaconda3/envs/myenv/bin/python write_deliverable_docs.py
"""

import os
import numpy as np
import pandas as pd
import xarray as xr
import zarr

ZARR_PATH = "ZARR/mores_creek_summit.zarr"
FRAG = "tex_fragments"
TOP_PAIRS = [8, 2, 1, 5, 6]
ZERO_PAIR = 4
from rank_dswe_pairs import PAIR_LETTERS  # shared letter labels (A-F)

COORD_MEANING = {  # one-line meaning per coordinate variable
    "x": "UTM easting of cell centers (m, EPSG:32611)",
    "y": "UTM northing of cell centers (m), stored north to south",
    "pol": "polarization (HH, HV, VH, VV)",
    "line": "UAVSAR flight line id (05208 NE-looking, 23205 SW-looking)",
    "pair": "InSAR pair index (order of uavsar\\_pair\\_manifest\\_mcs.json)",
    "pair_date1": "pair reference (earlier) acquisition date",
    "pair_date2": "pair repeat (later) acquisition date",
    "line_id": "flight line of each pair",
    "heading": "aircraft heading of each pair (deg)",
    "temporal_baseline_days": "days between the pair's acquisitions",
    "flight": "unique flight index (line, date)",
    "flight_date": "acquisition date of each flight",
    "flight_line_id": "flight line of each flight",
    "sd_time": "snow-depth acquisition date",
    "sd_source": "snow-depth product (SNEX20\\_QSI\\_SD or SNEX\\_MCS\\_Lidar)",
    "lidar_time": "MCS\\_Lidar acquisition date (DTM/DSM/CHM)",
}
VAR_NOTES = {  # units + provenance per data variable
    "wrapped_phase": ("complex (unitless)", "ASF INTERFEROMETRY\\_GRD, complex-domain regrid", "build\\_mcs\\_datacube.py"),
    "unwrapped_phase": ("radians", "ASF .unw.grd; pairs 9/12/14/20/21 SNAPHU-unwrapped in-house", "build\\_mcs\\_datacube.py, unwrap\\_missing\\_pairs.py"),
    "coherence": ("unitless [0,1]", "ASF .cor.grd", "build\\_mcs\\_datacube.py"),
    "sd": ("m", "QSI + MCS\\_Lidar snow depth, 0.5\\,m block-averaged to 3\\,m", "build\\_mcs\\_datacube.py"),
    "dtm": ("m", "MCS\\_Lidar bare-earth DTM", "build\\_mcs\\_datacube.py"),
    "dsm": ("m", "MCS\\_Lidar surface DSM", "build\\_mcs\\_datacube.py"),
    "chm": ("m", "MCS\\_Lidar canopy height", "build\\_mcs\\_datacube.py"),
    "lia": ("degrees", "JPL .llh/.lkv look vectors + 2023-02-09 DTM", "compute\\_local\\_incidence\\_angle.py"),
    "atm_delay": ("radians", "ERA5 delay below aircraft altitude, per flight", "compute\\_atm\\_delay.py"),
    "atm_delay_diff": ("radians", "per-pair delay difference (date2$-$date1)", "compute\\_atm\\_delay\\_diff.py"),
    "dswe": ("m w.e.", "Guneriussen-2001 inversion, $\\rho$=250, SNOTEL-anchored, pair 8 deramped", "compute\\_dswe.py"),
}


def esc(s):
    """Escape LaTeX specials in plain strings."""
    return str(s).replace("_", "\\_").replace("%", "\\%").replace("&", "\\&")


def dt(dtype):
    """LaTeX-safe, human-friendly dtype name ('<U10' -> 'str(10)')."""
    t = str(dtype)
    if t.startswith("<U"):
        return f"str({t[2:]})"
    return {"datetime64[ns]": "datetime"}.get(t, t)


if __name__ == "__main__":
    os.makedirs(FRAG, exist_ok=True)
    ds = xr.open_zarr(ZARR_PATH)
    zg = zarr.open(ZARR_PATH, mode="r")

    # ---- dimensions + coordinates table ----
    rows = [r"\begin{tabular}{llll}", r"\hline",
            r"name & size / dims & dtype & meaning \\", r"\hline"]
    for d, n in ds.sizes.items():
        rows.append(f"{esc(d)} (dim) & {n} & -- & {COORD_MEANING.get(d, '')} \\\\")
    for c in sorted(ds.coords):
        if c in ds.sizes:                                        # dimension coords already listed
            meaning = COORD_MEANING.get(c, "")
            rows.append(f"{esc(c)} & coord on {esc(c)} & {dt(ds[c].dtype)} & {meaning} \\\\")
        else:
            rows.append(f"{esc(c)} & coord on {esc(ds[c].dims[0])} & {dt(ds[c].dtype)} "
                        f"& {COORD_MEANING.get(c, '')} \\\\")
    rows += [r"\hline", r"\end{tabular}"]
    open(f"{FRAG}/zarr_dims.tex", "w").write("\n".join(rows) + "\n")

    # ---- data variables table ----
    rows = [r"\begin{tabular}{p{2.6cm}p{3.1cm}p{1.5cm}p{2.0cm}p{5.6cm}}", r"\hline",
            r"variable & dims & dtype & chunks & units / source \\", r"\hline"]
    for v in sorted(ds.data_vars):
        units, src, script = VAR_NOTES[v]
        rows.append(f"{esc(v)} & {esc('(' + ', '.join(ds[v].dims) + ')')} & {dt(ds[v].dtype)} "
                    f"& {esc(str(zg[v].chunks))} & {units}; {src} ({script}) \\\\")
    rows += [r"\hline", r"\end{tabular}"]
    open(f"{FRAG}/zarr_vars.tex", "w").write("\n".join(rows) + "\n")

    # ---- dataset counts: number of 2-D (y,x) slices per variable ----
    rows = [r"\begin{tabular}{lrl}", r"\hline",
            r"variable & 2-D maps & composition \\", r"\hline"]
    for v in sorted(ds.data_vars):
        extra = [d for d in ds[v].dims if d not in ("y", "x")]
        n = int(np.prod([ds.sizes[d] for d in extra])) if extra else 1
        comp = r" $\times$ ".join(f"{ds.sizes[d]} {esc(d)}" for d in extra) if extra else "single map"
        rows.append(f"{esc(v)} & {n} & {comp} \\\\")
    rows += [r"\hline", r"\end{tabular}"]
    open(f"{FRAG}/zarr_counts.tex", "w").write("\n".join(rows) + "\n")

    # ---- ranking tables from dswe_lidar_rank.csv ----
    df = pd.read_csv("dswe_lidar_rank.csv")
    summ = []
    for pi, sub in df.groupby("pair"):
        r0 = sub.iloc[0]
        matched = "r_lidar2021" if r0.season == 2021 else "r_lidar2020"
        summ.append(dict(pair=pi, line=r0.line, d1=r0.date1, d2=r0.date2,
                         snotel=r0.snotel_dswe_mm, r20=sub.r_lidar2020.mean(),
                         r21=sub.r_lidar2021.mean(), rm=sub[matched].mean(),
                         n90=int(sub[["n2020", "n2021"]].max().max()),
                         xp=r0.crosspol_mean, custom=bool(r0.custom_unwrap)))
    S = pd.DataFrame(summ).set_index("pair")

    rows = [r"\begin{tabular}{cclllrrrr}", r"\hline",
            r"label & rank & pair & line & dates & SNOTEL $\Delta$SWE (mm) & $r$ vs lidar & cross-pol & $n_{90}$ \\",
            r"\hline"]
    for i, pi in enumerate(TOP_PAIRS, 1):
        r = S.loc[pi]
        rows.append(f"\\textbf{{{PAIR_LETTERS[pi]}}} & {i} & {pi} & {r.line:05d} & {r.d1} $\\rightarrow$ {r.d2} & {r.snotel:+.0f} "
                    f"& {r.rm:+.3f} & {r.xp:+.2f} & {r.n90} \\\\")
    r = S.loc[ZERO_PAIR]
    rows.append(r"\hline")
    rows.append(f"\\textbf{{{PAIR_LETTERS[ZERO_PAIR]}}} & -- & {ZERO_PAIR} & {r.line:05d} & {r.d1} $\\rightarrow$ {r.d2} & {r.snotel:+.0f} "
                f"& {r.rm:+.3f} & {r.xp:+.2f} & {r.n90} \\\\")
    rows += [r"\hline", r"\end{tabular}"]
    open(f"{FRAG}/rank_top.tex", "w").write("\n".join(rows) + "\n")

    # appendix: every pair, both lidar correlations, gates, provenance
    rows = [r"\begin{tabular}{cllrrrrrl}", r"\hline",
            r"pair & line & dates & SNOTEL & $r$ 2020 & $r$ 2021 & $n_{90}$ & cross-pol & notes \\",
            r"\hline"]
    ranked = S[(S.snotel >= 10) & (S.n90 >= 500)].sort_values("rm", ascending=False)
    rank_of = {pi: i + 1 for i, pi in enumerate(ranked.index)}
    for pi, r in S.iterrows():
        notes = []
        if pi in PAIR_LETTERS: notes.append(f"\\textbf{{Pair {PAIR_LETTERS[pi]}}}")
        if pi in rank_of: notes.append(f"rank {rank_of[pi]}")
        if pi == ZERO_PAIR: notes.append("zero-change example")
        if r.snotel < 10: notes.append("no snowfall")
        if r.n90 < 500: notes.append("low $n$")
        if r.custom: notes.append("SNAPHU unwrap")
        rows.append(f"{pi} & {r.line:05d} & {r.d1} $\\rightarrow$ {r.d2} & {r.snotel:+.0f} "
                    f"& {r.r20:+.3f} & {r.r21:+.3f} & {r.n90} & {r.xp:+.2f} & {', '.join(notes)} \\\\")
    rows += [r"\hline", r"\end{tabular}"]
    open(f"{FRAG}/rank_appendix.tex", "w").write("\n".join(rows) + "\n")
    print(f"wrote 5 fragments to {FRAG}/")
