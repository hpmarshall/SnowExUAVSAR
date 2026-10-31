"""Local incidence angle (LIA) for the Mores Creek Summit UAVSAR lines 05208 and 23205,
computed with the textbook formula: the angle between the ground-to-sensor unit look
vector and the UNIT terrain surface normal,

    LIA = arccos( l_hat . n_hat ),   n_hat = (-dz/dE, -dz/dN, 1) / sqrt(1 + |grad z|^2)

This differs from uavsar_pytools' calc_inc_angle (lia_pytools.py), which skips the
normal's normalization, inflating cos(LIA) over sloped terrain. Flat terrain is identical.

Ships as a pair with lia_pytools.py: the segment geocoding helpers are imported from it,
so keep both files together. Same inputs; the geocoded look-vector cache is shared.
Outputs: LIA_PRODUCTS/lia_<line>_corrected.tif  (LIA in degrees, float32, NaN nodata)

Run:  PROJ_DATA=/opt/anaconda3/envs/myenv/share/proj \
        /opt/anaconda3/envs/myenv/bin/python lia_corrected.py
"""

import os                                                         # paths
import numpy as np                                                # array math
import rasterio as rio                                            # GeoTIFF I/O
from rasterio.warp import transform_bounds                        # grid bounds -> lat/lon
from lia_pytools import DEM_FP, OUT_DIR, LINES, PIXEL_M, geocode_or_load  # shared grid + geocoding


def corrected_lia(dem, east, north, up, pixel_m):
    """Textbook LIA (degrees) from a DEM and the east/north/up look-vector components."""
    d_row, d_col = np.gradient(dem, pixel_m)                      # terrain slope along grid axes (m/m)
    zx = d_col                                                    # dz/dEast: columns advance east (+3 m)
    zy = -d_row                                                   # dz/dNorth: rows advance SOUTH, so negate
    mag = np.sqrt(east**2 + north**2 + up**2)                     # look-vector length
    mag[mag == 0] = np.nan                                        # avoid division by zero
    le, ln, lu = -east / mag, -north / mag, -up / mag             # unit look vector, ground -> sensor
    #                                                               (.lkv points sensor -> ground, so negate)
    n_norm = np.sqrt(1.0 + zx**2 + zy**2)                         # |surface normal (-zx, -zy, 1)|
    cos_lia = (-zx * le - zy * ln + lu) / n_norm                  # l_hat . n_hat, both unit vectors
    return np.degrees(np.arccos(np.clip(cos_lia, -1.0, 1.0)))    # degrees; >90 = slope faces away (kept)


if __name__ == "__main__":
    with rio.open(DEM_FP) as src:                                 # the DEM defines the output grid
        dem = src.read(1)                                         # terrain heights (m)
        grid = {"shape": dem.shape, "crs": src.crs, "bounds": tuple(src.bounds),  # grid description
                "bounds_ll": transform_bounds(src.crs, "EPSG:4326", *src.bounds),  # for overlap tests
                "profile": {**src.profile, "compress": "deflate"}}  # output GeoTIFF profile

    for line in LINES:                                            # 05208 then 23205
        print(f"=== line {line} ===")
        lkv = geocode_or_load(line, grid)                         # cached by lia_pytools.py if run first
        lia = corrected_lia(dem, lkv["east"], lkv["north"], lkv["up"], PIXEL_M)  # the textbook angle
        out_fp = os.path.join(OUT_DIR, f"lia_{line}_corrected.tif")  # deliverable grid
        with rio.open(out_fp, "w", **grid["profile"]) as dst:     # float32, EPSG:32611, NaN nodata
            dst.write(lia.astype("float32"), 1)
        v = lia[np.isfinite(lia)]                                 # quick stats on valid pixels
        print(f"  wrote {out_fp}: {v.size} valid px, range {v.min():.1f}-{v.max():.1f}, mean {v.mean():.1f} deg")
