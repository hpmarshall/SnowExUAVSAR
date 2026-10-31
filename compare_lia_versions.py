"""Validate and compare the two LIA deliverables (lia_pytools.py vs lia_corrected.py).

Checks, per line (05208, 23205):
  1. the corrected GeoTIFF matches the zarr's lia (to float noise) -- the zarr was
     regenerated with the corrected formula on 2026-10-30 (before that it held the
     uavsar_pytools estimate, which this check validated against the pytools tif);
  2. over near-flat terrain (<2 deg slope) the two versions nearly agree -- both tend to
     arccos of the vertical look component at zero slope (~1 deg residual remains because
     calc_inc_angle's sign error acts even at small slopes);
  3. the mean ground->sensor look azimuth matches the line geometry (UAVSAR is
     left-looking: heading 052 -> ~142 deg azimuth, heading 232 -> ~322 deg);
  4. corrected - pytools difference statistics, binned by terrain slope;
  5. the decisive physics check: slopes facing the radar (downslope aspect toward the
     sensor) must show LIA ~ flat incidence - slope, away-facing ~ flat + slope.
Writes LIA_PRODUCTS/lia_comparison.png (maps + difference-vs-slope summary).

Result (2026-10-29 run): the corrected version passes check 5 to within ~1 deg; the
pytools version misses by 20+ deg both ways -- calc_inc_angle pairs/signs the horizontal
look components inconsistently with the DEM gradient in addition to skipping the
surface-normal normalization, so lia_corrected.py is the scientifically right product.

Run:  PROJ_DATA=/opt/anaconda3/envs/myenv/share/proj \
        /opt/anaconda3/envs/myenv/bin/python compare_lia_versions.py
"""

import os                                                         # paths
import numpy as np                                                # array math
import rasterio as rio                                            # GeoTIFF reader
import xarray as xr                                               # zarr reader
import matplotlib                                                 # headless backend
matplotlib.use("Agg")                                             # no display needed
import matplotlib.pyplot as plt                                   # figure
from lia_pytools import DEM_FP, OUT_DIR, LINES, PIXEL_M, COMPS, read_tif  # shared constants/reader

EXPECTED_AZ = {"05208": 142.0, "23205": 322.0}                    # ground->sensor azimuth = heading - 90 + 180

dem = read_tif(DEM_FP)                                            # terrain
d_row, d_col = np.gradient(dem, PIXEL_M)                          # slope components (m/m)
zx, zy = d_col, -d_row                                            # dz/dEast, dz/dNorth (rows run south)
slope = np.degrees(np.arctan(np.hypot(zx, zy)))                   # terrain slope in degrees
aspect = np.degrees(np.arctan2(-zx, -zy)) % 360                   # downslope azimuth (direction the slope faces)
zarr_lia = xr.open_zarr("ZARR/mores_creek_summit.zarr")["lia"]    # the zarr's lia (corrected since 2026-10-30)

fig, axes = plt.subplots(2, 4, figsize=(20, 9))                   # one row per line
for row, line in enumerate(LINES):                                # 05208 then 23205
    pyt = read_tif(os.path.join(OUT_DIR, f"lia_{line}_pytools.tif"))      # uavsar_pytools LIA
    cor = read_tif(os.path.join(OUT_DIR, f"lia_{line}_corrected.tif"))    # textbook LIA
    ref = zarr_lia.sel(line=line).values                          # zarr reference
    both = np.isfinite(cor) & np.isfinite(ref)                    # common valid pixels

    print(f"=== line {line} ===")
    print(f"  [1] corrected tif vs zarr lia: max |diff| = {np.abs(cor[both] - ref[both]).max():.6f} deg "
          f"on {both.sum()} px (valid px: tif {np.isfinite(cor).sum()}, zarr {np.isfinite(ref).sum()})")

    ok = np.isfinite(pyt) & np.isfinite(cor)                      # pixels valid in both versions
    diff = np.where(ok, cor - pyt, np.nan)                        # corrected minus pytools
    flat = ok & (slope < 2)                                       # near-flat terrain
    print(f"  [2] flat (<2 deg) agreement: mean |diff| = {np.nanmean(np.abs(diff[flat])):.3f} deg, "
          f"p99 = {np.nanpercentile(np.abs(diff[flat]), 99):.3f} deg ({flat.sum()} px)")

    e = read_tif(os.path.join(OUT_DIR, f"lkv_{line}_east.tif"))   # geocoded look-vector components,
    n = read_tif(os.path.join(OUT_DIR, f"lkv_{line}_north.tif"))  # for the geometry sanity check
    az = (np.degrees(np.arctan2(np.nanmean(-e), np.nanmean(-n)))) % 360  # mean ground->sensor azimuth
    print(f"  [3] mean look azimuth {az:.1f} deg (expected ~{EXPECTED_AZ[line]:.0f} for this heading)")

    for lo, hi in ((0, 10), (10, 20), (20, 30), (30, 90)):        # difference by slope class
        sel = ok & (slope >= lo) & (slope < hi)                   # pixels in this slope bin
        print(f"  [4] slope {lo:2d}-{hi:2d} deg: mean diff {np.nanmean(diff[sel]):+.2f}, "
              f"mean |diff| {np.nanmean(np.abs(diff[sel])):.2f}, p99 |diff| "
              f"{np.nanpercentile(np.abs(diff[sel]), 99):.2f} deg  ({sel.sum()} px)")

    flat_inc = np.nanmean(cor[ok & (slope < 2)])                  # reference incidence on flat ground
    band = ok & (slope > 20) & (slope < 25)                       # one steep, well-populated slope class
    for name, center, sign in (("facing", EXPECTED_AZ[line], -1),               # faces radar: aspect toward sensor
                               ("away", (EXPECTED_AZ[line] + 180) % 360, +1)):  # faces away: aspect opposite
        sel = band & (np.abs((aspect - center + 180) % 360 - 180) < 30)         # within 30 deg of that aspect
        print(f"  [5] {name:6s} 20-25 deg slopes: expected ~{flat_inc + sign * 22.5:5.1f}, "
              f"corrected {np.nanmean(cor[sel]):5.1f}, pytools {np.nanmean(pyt[sel]):5.1f} deg  ({sel.sum()} px)")

    for col, (arr, title, cmap, vr) in enumerate(                 # three maps + one summary panel
            [(pyt, f"{line} pytools LIA", "viridis", (10, 90)),
             (cor, f"{line} corrected LIA", "viridis", (10, 90)),
             (diff, f"{line} corrected - pytools", "RdBu_r", (-10, 10))]):
        im = axes[row, col].imshow(arr, cmap=cmap, vmin=vr[0], vmax=vr[1])  # map panel
        axes[row, col].set_title(title); axes[row, col].axis("off")         # label, no ticks
        fig.colorbar(im, ax=axes[row, col], shrink=0.7, label="deg")        # scale
    sel = ok & (slope < 45)                                       # summary: |diff| vs slope, 1-deg bins
    bins = np.arange(0, 46)                                       # slope bin edges
    idx = np.digitize(slope[sel], bins)                           # bin index per pixel
    med = [np.nanmedian(np.abs(diff[sel])[idx == i]) for i in range(1, len(bins))]  # median |diff| per bin
    axes[row, 3].plot(bins[:-1] + 0.5, med)                       # the curve
    axes[row, 3].set(xlabel="terrain slope (deg)", ylabel="median |corrected - pytools| (deg)",
                     title=f"{line}: version difference vs slope")  # labels

fig.tight_layout()                                                # pack panels
fig.savefig(os.path.join(OUT_DIR, "lia_comparison.png"), dpi=130)  # quicklook deliverable
print("wrote", os.path.join(OUT_DIR, "lia_comparison.png"))       # confirm
