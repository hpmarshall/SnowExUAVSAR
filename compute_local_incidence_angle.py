"""Compute per-line local incidence angle (LIA) for Mores Creek Summit UAVSAR lines,
using JPL Stack .llh/.lkv look-vector geometry and our regridded lidar DTM.

uavsar_pytools.georeference.geolocate_uavsar() hardcodes segment "_1_" when reading
row/col shapes from the .ann file, so it silently mis-shapes segments 2/3 whenever
their pixel counts differ from segment 1 (true here for both lines' segment 3). This
script reimplements just the shape-reading + GDAL geolocation-array warp steps per
segment (reusing the library's own geocodeUsingGdalWarp, which is segment-agnostic),
then mosaics the 3 non-overlapping along-track segments and warps straight onto our
canonical 3 m grid in one step, instead of geolocate_uavsar's native-grid-then-reproject.
"""

import os                                                       # paths
from glob import glob                                           # find a representative .ann per line
import numpy as np                                               # array math
import rasterio as rio                                           # read back warped GeoTIFFs
from rasterio.warp import transform_bounds                       # canonical grid bounds -> lat/lon, for overlap checks
import xarray as xr                                               # datacube I/O
from uavsar_pytools.convert.tiff_conversion import read_annotation  # parse .ann key/value pairs
from uavsar_pytools.georeference import geocodeUsingGdalWarp     # GDAL geolocation-array warp (segment-agnostic)
from uavsar_pytools.incidence_angle import calc_inc_angle        # LIA from DEM + look vector

from build_mcs_datacube import GRID_CRS, GRID_TRANSFORM, GRID_WIDTH, GRID_HEIGHT, GRID_Y, GRID_X  # canonical grid

LKV_DIR = "UAVSAR_LKV"                                           # downloaded JPL Stack files (git-ignored)
SCRATCH_DIR = os.path.join(LKV_DIR, "tmp_geocode")               # scratch GeoTIFFs for GDAL warp, cleaned up after
ZARR_PATH = "ZARR/mores_creek_summit.zarr"
LINE_SUFFIX = {"05208": "BU", "23205": "BC"}                     # JPL Stack filename suffix per line
PIXEL_SIZE_M = 3.0                                               # canonical grid spacing

# canonical grid bounds, derived the same way the grid itself was defined
GRID_BOUNDS = (GRID_TRANSFORM.c, GRID_TRANSFORM.f + GRID_TRANSFORM.e * GRID_HEIGHT,
               GRID_TRANSFORM.c + GRID_TRANSFORM.a * GRID_WIDTH, GRID_TRANSFORM.f)
GRID_BOUNDS_LATLON = transform_bounds(GRID_CRS, "EPSG:4326", *GRID_BOUNDS)  # for cheap per-segment overlap checks


def write_scratch_tif(arr, out_fp):
    """Plain (non-georeferenced) single-band float32 GeoTIFF, as GDAL geolocation warp expects."""
    profile = {"driver": "GTiff", "dtype": "float32", "count": 1,
               "width": arr.shape[1], "height": arr.shape[0], "nodata": np.nan}
    with rio.open(out_fp, "w", **profile) as dst:
        dst.write(arr.astype("float32"), 1)


def segment_shape(desc, kind, seg):
    """Row/col count for one .llh or .lkv segment, read directly from the .ann (not hardcoded to segment 1)."""
    return desc[f"{kind}_{seg}_2x8.set_rows"]["value"], desc[f"{kind}_{seg}_2x8.set_cols"]["value"]


def geocode_segment(llh_fp, lkv_fp, nrows, ncols, out_prefix):
    """Warp one segment's 3 look-vector components onto the canonical grid; returns 3 arrays (raw0, raw1, raw2).
    Each UAVSAR line is split into 3 along-track segments spanning tens of km; only the one segment that
    actually passes over Mores Creek Summit overlaps our small target footprint -- skip the other two outright
    (each full segment is 400-600 MB and a GDAL geolocation warp costs ~80s, so this isn't just cosmetic)."""
    llh = np.fromfile(llh_fp, dtype="<f4")                       # flat lat/lon/height, interleaved by 3
    lat, lon = llh[::3].reshape(nrows, ncols), llh[1::3].reshape(nrows, ncols)
    lon_min, lat_min, lon_max, lat_max = GRID_BOUNDS_LATLON
    if lon.max() < lon_min or lon.min() > lon_max or lat.max() < lat_min or lat.min() > lat_max:
        empty = np.full((GRID_HEIGHT, GRID_WIDTH), np.nan, dtype="float32")
        return [empty, empty, empty]
    lkv = np.fromfile(lkv_fp, dtype="<f4")                       # flat look-vector components, interleaved by 3
    os.makedirs(SCRATCH_DIR, exist_ok=True)
    lat_fp, lon_fp = f"{out_prefix}_lat.tif", f"{out_prefix}_lon.tif"
    write_scratch_tif(lat, lat_fp)                               # geolocation-array inputs for GDAL warp
    write_scratch_tif(lon, lon_fp)
    raws = []
    for comp in range(3):                                        # raw0/raw1/raw2 -> caller maps to y/x/z
        comp_fp = f"{out_prefix}_c{comp}.tif"
        warped_fp = f"{out_prefix}_c{comp}_geocoded.tif"
        write_scratch_tif(lkv[comp::3].reshape(nrows, ncols), comp_fp)
        geocodeUsingGdalWarp(infile=comp_fp, latfile=lat_fp, lonfile=lon_fp, outfile=warped_fp,
                              outsrs=GRID_CRS, bounds=GRID_BOUNDS, spacing=[PIXEL_SIZE_M, PIXEL_SIZE_M],
                              method="bilinear")                  # continuous vector component, no wrapping concern
        with rio.open(warped_fp) as src:
            arr = src.read(1)
            arr = np.where(arr == src.nodata, np.nan, arr)
        raws.append(arr)
    return raws


def geocode_line(line):
    """Mosaic all 3 segments of one line's look vector onto the canonical grid -> (lkv_x, lkv_y, lkv_z)."""
    suffix = LINE_SUFFIX[line]
    ann_fp = sorted(glob(os.path.join(LKV_DIR, f"lowman_{line}_*_L090HH_01_{suffix}.ann")))[0]
    desc = read_annotation(ann_fp)                                # segment shapes are identical across flights (Step 0 finding)
    merged = [np.full((GRID_HEIGHT, GRID_WIDTH), np.nan, dtype="float32") for _ in range(3)]
    for seg in (1, 2, 3):
        nrows, ncols = segment_shape(desc, "llh", seg)
        llh_fp = os.path.join(LKV_DIR, f"lowman_{line}_01_{suffix}_s{seg}_2x8.llh")
        lkv_fp = os.path.join(LKV_DIR, f"lowman_{line}_01_{suffix}_s{seg}_2x8.lkv")
        seg_raws = geocode_segment(llh_fp, lkv_fp, nrows, ncols, os.path.join(SCRATCH_DIR, f"{line}_s{seg}"))
        for i in range(3):                                        # segments don't overlap in azimuth -> first-valid fill
            fill = np.isnan(merged[i]) & ~np.isnan(seg_raws[i])
            merged[i][fill] = seg_raws[i][fill]
    # raw0/raw1/raw2 -> lkv_y/lkv_x/lkv_z: matches uavsar_pytools' own (confusingly swapped but self-consistent)
    # convention, where calc_inc_angle's "x" pairs with the DEM's row-axis gradient and "y" with its column-axis one
    return {"lkv_y": merged[0], "lkv_x": merged[1], "lkv_z": merged[2]}


if __name__ == "__main__":
    ds = xr.open_zarr(ZARR_PATH)
    dem = ds["dtm"].sel(lidar_time="2023-02-09").squeeze().values  # shared terrain surface (Track A plan, Step 3)

    lia_by_line = {}
    for line in LINE_SUFFIX:
        print(f"=== geocoding look vector: line {line} ===")
        lkv = geocode_line(line)
        print(f"  computing LIA for line {line}")
        lia = calc_inc_angle(dem=dem, lkv_x=lkv["lkv_x"], lkv_y=lkv["lkv_y"], lkv_z=lkv["lkv_z"],
                              pixel_size=PIXEL_SIZE_M)
        valid = lia[np.isfinite(lia)]
        print(f"  line {line}: {valid.size} valid px, range {valid.min():.1f}-{valid.max():.1f} deg, "
              f"mean {valid.mean():.1f} deg")
        lia_by_line[line] = lia

    lia_da = xr.DataArray(np.stack([lia_by_line[l] for l in LINE_SUFFIX], axis=0),
                           dims=("line", "y", "x"),
                           coords={"line": list(LINE_SUFFIX), "y": GRID_Y, "x": GRID_X},
                           name="lia", attrs={"units": "degrees", "long_name": "Local incidence angle"})
    xr.Dataset({"lia": lia_da}).to_zarr(ZARR_PATH, mode="a")     # add alongside existing variables
    print("Stored 'lia' in", ZARR_PATH)
