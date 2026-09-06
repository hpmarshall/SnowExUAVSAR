"""Spatial correlation of lidar snow depth vs InSAR phase, per MCS UAVSAR pair.

Screens the 22 cataloged pairs for snow-modeling usability: for each pair x
polarization, correlates QSI lidar snow depth (2020-02-09 for winter-2020 pairs,
2021-03-15 for 2021 pairs) against
  - unwrapped phase:  plain signed Pearson r, and
  - wrapped phase:    Mardia's circular-linear correlation R in [0,1]
                      (SD vs cos(phi) and sin(phi) jointly -- wrap-proof, unlike
                      Pearson on the bare angle).
Both fields are block-averaged from 3 m to 90 m first (30x30 pixels; the wrapped
phase is averaged in the complex domain) so speckle doesn't crush the statistic;
native 3 m values are also written to the CSV for reference.

Outputs: phase_sd_correlations.csv (every pair x pol, both scales)
         phase_sd_correlations.md  (90 m summary table, one row per pair)

Run:  PROJ_DATA=/opt/anaconda3/envs/myenv/share/proj \
      /opt/anaconda3/envs/myenv/bin/python compute_phase_sd_correlations.py [--check]
"""

import sys                     # --check flag
import numpy as np             # correlation math
import pandas as pd            # CSV/markdown output
import xarray as xr            # datacube access

ZARR_PATH = "ZARR/mores_creek_summit.zarr"
BLOCK = 30                     # 30 x 3 m = 90 m averaging blocks
POLS = ["HH", "HV", "VH", "VV"]
LIDAR_FOR_YEAR = {2020: "2020-02-09", 2021: "2021-03-15"}        # QSI depth per winter


def block_mean(a, k=BLOCK):
    """NaN-aware k x k block average (trims edge remainder); complex input allowed."""
    ny, nx = (a.shape[0] // k) * k, (a.shape[1] // k) * k        # trim to whole blocks
    blocks = a[:ny, :nx].reshape(ny // k, k, nx // k, k)
    if np.iscomplexobj(a):                                        # complex mean: Re/Im separately
        return (np.nanmean(blocks.real, axis=(1, 3)) + 1j * np.nanmean(blocks.imag, axis=(1, 3)))
    return np.nanmean(blocks, axis=(1, 3))


def pearson(x, y):
    """Pearson r over jointly-finite samples; (nan, 0) when degenerate."""
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 10 or np.nanstd(x[m]) == 0 or np.nanstd(y[m]) == 0:
        return np.nan, int(m.sum())
    return float(np.corrcoef(x[m], y[m])[0, 1]), int(m.sum())


def circ_lin_R(x, phi):
    """Mardia's circular-linear correlation R in [0,1] between linear x and angle phi."""
    m = np.isfinite(x) & np.isfinite(phi)
    if m.sum() < 10:
        return np.nan, int(m.sum())
    x, c, s = x[m], np.cos(phi[m]), np.sin(phi[m])
    rxc, _ = pearson(x, c)                                       # x vs cos component
    rxs, _ = pearson(x, s)                                       # x vs sin component
    rcs, _ = pearson(c, s)                                       # cos vs sin
    if np.isnan(rxc) or np.isnan(rxs) or np.isnan(rcs) or abs(rcs) >= 1:
        return np.nan, int(m.sum())
    R2 = (rxc**2 + rxs**2 - 2 * rxc * rxs * rcs) / (1 - rcs**2)
    return float(np.sqrt(max(R2, 0.0))), int(m.sum())


if __name__ == "__main__":
    ds = xr.open_zarr(ZARR_PATH)
    sd = {y: ds["sd"].sel(sd_time=d).squeeze().values for y, d in LIDAR_FOR_YEAR.items()}
    sd90 = {y: block_mean(v) for y, v in sd.items()}             # pre-average once per lidar date

    rows = []
    for pi in range(ds.sizes["pair"]):
        line = str(ds["line_id"][pi].values)
        d1, d2 = str(ds["pair_date1"][pi].values)[:10], str(ds["pair_date2"][pi].values)[:10]
        year = int(d2[:4])                                        # winter season keyed by the later date
        base = int(ds["temporal_baseline_days"][pi])
        med_coh = float(np.nanmedian(ds["coherence"][pi, 0, ::4, ::4].values))  # HH context
        print(f"[{pi:2d}] {line} {d1} -> {d2} ({base}d, lidar {LIDAR_FOR_YEAR[year]})")

        for pol_i, pol in enumerate(POLS):
            unw = ds["unwrapped_phase"][pi, pol_i].values
            wrp = ds["wrapped_phase"][pi, pol_i].values
            # 90 m block averages: plain for unwrapped, complex-domain for wrapped
            unw90 = block_mean(unw)
            phi90 = np.angle(block_mean(wrp))
            phi90[~np.isfinite(block_mean(np.abs(wrp)))] = np.nan  # keep NaN where no wrapped data
            r_unw90, n90 = pearson(sd90[year], unw90)
            R_wrap90, nw90 = circ_lin_R(sd90[year], phi90)
            # native 3 m values (CSV reference only)
            r_unw3, _ = pearson(sd[year], unw)
            phi3 = np.where(np.isfinite(wrp.real), np.angle(wrp), np.nan)
            R_wrap3, _ = circ_lin_R(sd[year], phi3)
            rows.append(dict(pair=pi, line=line, date1=d1, date2=d2, baseline_days=base,
                             lidar=LIDAR_FOR_YEAR[year], pol=pol, med_coh_hh=round(med_coh, 3),
                             n90_unw=n90, n90_wrap=nw90,
                             r_unw_90m=r_unw90, R_wrap_90m=R_wrap90,
                             r_unw_3m=r_unw3, R_wrap_3m=R_wrap3))

    df = pd.DataFrame(rows).round({"r_unw_90m": 3, "R_wrap_90m": 3, "r_unw_3m": 3, "R_wrap_3m": 3})
    df.to_csv("phase_sd_correlations.csv", index=False)
    print("wrote phase_sd_correlations.csv")

    # sanity checks promised in the plan
    assert df["r_unw_90m"].dropna().between(-1, 1).all()
    assert df["R_wrap_90m"].dropna().between(0, 1).all()

    if "--check" in sys.argv:                                     # independent spot re-computation, pair 3 HH
        pi, pol_i, year = 3, 0, 2021
        a = block_mean(ds["unwrapped_phase"][pi, pol_i].values)
        b = sd90[year]
        m = np.isfinite(a) & np.isfinite(b)
        r = np.corrcoef(a[m], b[m])[0, 1]
        stored = df[(df.pair == pi) & (df.pol == "HH")].r_unw_90m.iloc[0]
        print(f"--check pair 3 HH: independent r={r:.3f}, table r={stored:.3f}")
        assert abs(r - stored) < 5e-3

    # markdown summary: one row per pair, 90 m stats across the 4 pols
    lines_out = ["# Snow depth vs phase correlation, per UAVSAR pair (90 m)", "",
                 "Lidar snow depth (QSI 2020-02-09 for winter-2020 pairs, 2021-03-15 for 2021 pairs)",
                 "vs unwrapped phase (Pearson r, signed) and wrapped phase (Mardia circular-linear R, in [0,1]).",
                 "Fields block-averaged to 90 m before correlating; 3 m values in phase_sd_correlations.csv.", "",
                 "| # | line | date1 | date2 | dt (d) | med coh HH | n (90m) | "
                 + " | ".join(f"r_unw {p}" for p in POLS) + " | "
                 + " | ".join(f"R_wrap {p}" for p in POLS) + " |",
                 "|---|------|-------|-------|--------|-----------|---------|"
                 + "|".join(["--------"] * 8) + "|"]
    fmt = lambda v: "" if pd.isna(v) else f"{v:+.3f}"            # blank cell where no data
    fmtu = lambda v: "" if pd.isna(v) else f"{v:.3f}"            # unsigned for R
    for pi in range(ds.sizes["pair"]):
        sub = df[df.pair == pi].set_index("pol")
        r0 = sub.iloc[0]
        lines_out.append(
            f"| {pi} | {r0.line} | {r0.date1} | {r0.date2} | {r0.baseline_days} | {r0.med_coh_hh:.3f} "
            f"| {r0.n90_unw} | " + " | ".join(fmt(sub.loc[p, 'r_unw_90m']) for p in POLS)
            + " | " + " | ".join(fmtu(sub.loc[p, 'R_wrap_90m']) for p in POLS) + " |")
    lines_out += ["", "Pairs 9, 12, 14, 20, 21: ASF ships no .unw.grd; their unwrapped phase comes from",
                  "our SNAPHU pipeline (unwrap_missing_pairs.py, 5x5 looks, smooth cost).",
                  "n (90m): jointly-finite 90 m cells in the unwrapped comparison (HH)."]
    open("phase_sd_correlations.md", "w").write("\n".join(lines_out) + "\n")
    print("wrote phase_sd_correlations.md")
