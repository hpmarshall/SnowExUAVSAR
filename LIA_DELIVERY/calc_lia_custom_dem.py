"""Compute local incidence angle (LIA) for a Mores Creek Summit UAVSAR line with YOUR
own DEM, using the delivered geocoded look vectors. The DEM may be in any CRS and any
resolution: it is reprojected/resampled onto the 3 m UTM-11N look-vector grid first,
then the textbook LIA = arccos(l_hat . n_hat) is computed (see lia_corrected.py).

Usage:  python calc_lia_custom_dem.py <your_dem.tif> <line: 05208|23205> [out.tif]
Needs:  LIA_PRODUCTS/lkv_<line>_{east,north,up}.tif  (in the delivery package)
        lia_pytools.py + lia_corrected.py next to this script
"""

import sys                                                        # command-line arguments
import numpy as np                                                # array math
import rasterio as rio                                            # GeoTIFF I/O
from rasterio.warp import reproject, Resampling                   # DEM -> look-vector grid
from lia_pytools import OUT_DIR, PIXEL_M, read_tif                # delivery-package constants/reader
from lia_corrected import corrected_lia                           # the textbook LIA formula

dem_fp, line = sys.argv[1], sys.argv[2]                           # required arguments
out_fp = sys.argv[3] if len(sys.argv) > 3 else f"{OUT_DIR}/lia_{line}_customdem.tif"  # optional output path

with rio.open(f"{OUT_DIR}/lkv_{line}_east.tif") as ref:           # the look-vector grid is the target grid
    profile = {**ref.profile, "compress": "deflate"}              # reuse its georeferencing for the output
with rio.open(dem_fp) as src:                                     # the user's DEM, any CRS/resolution
    dem = np.full((profile["height"], profile["width"]), np.nan, dtype="float32")  # target array
    reproject(rio.band(src, 1), dem,                              # warp DEM onto the look-vector grid
              dst_transform=profile["transform"], dst_crs=profile["crs"],
              dst_nodata=np.nan, resampling=Resampling.bilinear)  # bilinear: terrain is continuous

lkv = {c: read_tif(f"{OUT_DIR}/lkv_{line}_{c}.tif") for c in ("east", "north", "up")}  # ENU look vector
lia = corrected_lia(dem, lkv["east"], lkv["north"], lkv["up"], PIXEL_M)  # LIA in degrees

with rio.open(out_fp, "w", **profile) as dst:                     # write float32 GeoTIFF, NaN nodata
    dst.write(lia.astype("float32"), 1)
v = lia[np.isfinite(lia)]                                         # stats on valid pixels
print(f"wrote {out_fp}: {v.size} valid px, range {v.min():.1f}-{v.max():.1f}, mean {v.mean():.1f} deg")
