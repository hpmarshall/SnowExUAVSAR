"""Local incidence angle (LIA) for the Mores Creek Summit UAVSAR lines 05208 and 23205:
the uavsar_pytools estimate, i.e. calc_inc_angle applied to the JPL InSAR-Stack look
vectors geocoded onto the 3 m lidar-DTM grid, exactly as stored in the MCS datacube.
Note this library routine does not normalize the terrain surface-normal vector, so on
steep slopes it deviates from the textbook LIA -- see lia_corrected.py for that version.

Inputs : LIA_PRODUCTS/mcs_dtm_3m.tif                       terrain; defines the output grid
         UAVSAR_LKV/lowman_<line>_01_<sfx>_s<1-3>_2x8.llh  per-pixel lat/lon/height (JPL Stack)
         UAVSAR_LKV/lowman_<line>_01_<sfx>_s<1-3>_2x8.lkv  look vector east/north/up (JPL Stack)
         UAVSAR_LKV/lowman_<line>_*_L090HH_01_<sfx>.ann    annotation with segment shapes
Outputs: LIA_PRODUCTS/lia_<line>_pytools.tif               LIA in degrees, float32, NaN nodata
         LIA_PRODUCTS/lkv_<line>_{east,north,up}.tif       geocoded look vector (cached, reused)

Run:  PROJ_DATA=/opt/anaconda3/envs/myenv/share/proj \
        /opt/anaconda3/envs/myenv/bin/python lia_pytools.py
"""

import os                                                         # paths
import shutil                                                     # scratch-dir cleanup
from glob import glob                                             # find one .ann per line
import numpy as np                                                # array math
import rasterio as rio                                            # GeoTIFF I/O
from rasterio.warp import transform_bounds                        # grid bounds -> lat/lon for overlap test
from uavsar_pytools.convert.tiff_conversion import read_annotation   # parse .ann key/value pairs
from uavsar_pytools.georeference import geocodeUsingGdalWarp      # GDAL geolocation-array warp
from uavsar_pytools.incidence_angle import calc_inc_angle         # the uavsar_pytools LIA routine

DEM_FP = "LIA_PRODUCTS/mcs_dtm_3m.tif"                            # terrain input (export_lia_inputs.py)
LKV_DIR = "UAVSAR_LKV"                                            # JPL Stack files
OUT_DIR = "LIA_PRODUCTS"                                          # all outputs land here
TMP_DIR = os.path.join(OUT_DIR, "tmp_geocode")                    # scratch for GDAL warp, removed at end
LINES = {"05208": "BU", "23205": "BC"}                            # line number -> JPL Stack filename suffix
PIXEL_M = 3.0                                                     # grid spacing in metres
COMPS = ("east", "north", "up")                                   # .lkv stores these 3, interleaved per pixel


def write_plain_tif(arr, fp):
    """Non-georeferenced float32 GeoTIFF: the input format GDAL's geolocation warp expects."""
    profile = {"driver": "GTiff", "dtype": "float32", "count": 1,  # single band
               "width": arr.shape[1], "height": arr.shape[0], "nodata": np.nan}  # no CRS on purpose
    with rio.open(fp, "w", **profile) as dst:                     # write it
        dst.write(arr.astype("float32"), 1)


def read_tif(fp):
    """Read band 1 of a GeoTIFF as float with nodata -> NaN."""
    with rio.open(fp) as src:                                     # open
        arr = src.read(1).astype("float32")                       # band 1
        if src.nodata is not None and not np.isnan(src.nodata):   # map numeric nodata
            arr[arr == src.nodata] = np.nan                       # to NaN
        return arr


def geocode_segment(line, sfx, seg, nrows, ncols, grid):
    """Warp one segment's 3 look-vector components onto the DEM grid -> [east, north, up] arrays.
    Each line has 3 along-track segments; only one covers MCS, so segments whose lat/lon
    never enter the DEM bounds are skipped outright (each costs ~minutes of GDAL warping)."""
    llh = np.fromfile(os.path.join(LKV_DIR, f"lowman_{line}_01_{sfx}_s{seg}_2x8.llh"), dtype="<f4")  # lat/lon/height
    lat, lon = llh[::3].reshape(nrows, ncols), llh[1::3].reshape(nrows, ncols)  # de-interleave
    lon_min, lat_min, lon_max, lat_max = grid["bounds_ll"]        # DEM footprint in lat/lon
    if lon.max() < lon_min or lon.min() > lon_max or lat.max() < lat_min or lat.min() > lat_max:  # no overlap?
        return None                                               # skip this segment
    lkv = np.fromfile(os.path.join(LKV_DIR, f"lowman_{line}_01_{sfx}_s{seg}_2x8.lkv"), dtype="<f4")  # look vectors
    pre = os.path.join(TMP_DIR, f"{line}_s{seg}")                 # scratch filename prefix
    write_plain_tif(lat, f"{pre}_lat.tif")                        # geolocation arrays
    write_plain_tif(lon, f"{pre}_lon.tif")                        # for the GDAL warp
    out = []                                                      # one warped array per component
    for i, comp in enumerate(COMPS):                              # east, north, up
        write_plain_tif(lkv[i::3].reshape(nrows, ncols), f"{pre}_{comp}.tif")  # de-interleave component i
        geocodeUsingGdalWarp(infile=f"{pre}_{comp}.tif", latfile=f"{pre}_lat.tif", lonfile=f"{pre}_lon.tif",
                             outfile=f"{pre}_{comp}_geo.tif", outsrs=grid["crs"], bounds=grid["bounds"],
                             spacing=[PIXEL_M, PIXEL_M], method="bilinear")  # continuous field -> bilinear
        out.append(read_tif(f"{pre}_{comp}_geo.tif"))             # read warped result back
    return out                                                    # [east, north, up] on the DEM grid


def geocode_or_load(line, grid):
    """Return {east, north, up} look-vector arrays on the DEM grid, geocoding only on first run
    (results are cached as LIA_PRODUCTS/lkv_<line>_*.tif and reused by lia_corrected.py too)."""
    cache = {c: os.path.join(OUT_DIR, f"lkv_{line}_{c}.tif") for c in COMPS}  # cache file per component
    if all(os.path.exists(fp) for fp in cache.values()):          # already geocoded?
        return {c: read_tif(fp) for c, fp in cache.items()}       # just load
    sfx = LINES[line]                                             # JPL Stack filename suffix
    ann = sorted(glob(os.path.join(LKV_DIR, f"lowman_{line}_*_L090HH_01_{sfx}.ann")))[0]  # any flight's .ann works
    desc = read_annotation(ann)                                   # segment shapes are flight-invariant
    os.makedirs(TMP_DIR, exist_ok=True)                           # scratch space
    merged = {c: np.full(grid["shape"], np.nan, dtype="float32") for c in COMPS}  # start empty
    for seg in (1, 2, 3):                                         # per-segment shapes read from the .ann --
        nrows = desc[f"llh_{seg}_2x8.set_rows"]["value"]          # NOT hardcoded to segment 1 (that is the
        ncols = desc[f"llh_{seg}_2x8.set_cols"]["value"]          # uavsar_pytools bug the datacube build worked around)
        arrs = geocode_segment(line, sfx, seg, nrows, ncols, grid)  # warp or skip
        if arrs is None:                                          # segment doesn't cover MCS
            print(f"  line {line} segment {seg}: no overlap, skipped")
            continue                                              # next segment
        for c, arr in zip(COMPS, arrs):                           # segments don't overlap in azimuth,
            fill = np.isnan(merged[c]) & ~np.isnan(arr)           # so a first-valid mosaic
            merged[c][fill] = arr[fill]                           # fills each pixel once
    for c in COMPS:                                               # cache for reuse
        with rio.open(cache[c], "w", **grid["profile"]) as dst:   # georeferenced this time
            dst.write(merged[c], 1)
    shutil.rmtree(TMP_DIR)                                        # drop warp scratch files
    return merged                                                 # {east, north, up}


if __name__ == "__main__":
    with rio.open(DEM_FP) as src:                                 # the DEM defines the output grid
        dem = src.read(1)                                         # terrain heights (m)
        grid = {"shape": dem.shape, "crs": src.crs, "bounds": tuple(src.bounds),  # grid description
                "bounds_ll": transform_bounds(src.crs, "EPSG:4326", *src.bounds),  # for overlap tests
                "profile": {**src.profile, "compress": "deflate"}}  # output GeoTIFF profile

    for line in LINES:                                            # 05208 then 23205
        print(f"=== line {line} ===")
        lkv = geocode_or_load(line, grid)                         # look vector on the DEM grid
        # uavsar_pytools convention: its "lkv_x" pairs with the DEM row-axis gradient and
        # "lkv_y" with the column-axis one, so north goes in as x and east as y -- this swapped
        # but self-consistent pairing is exactly what built the datacube's lia variable.
        lia = calc_inc_angle(dem=dem, lkv_x=lkv["north"], lkv_y=lkv["east"], lkv_z=lkv["up"],
                             pixel_size=PIXEL_M)                  # LIA in degrees
        out_fp = os.path.join(OUT_DIR, f"lia_{line}_pytools.tif") # deliverable grid
        with rio.open(out_fp, "w", **grid["profile"]) as dst:     # float32, EPSG:32611, NaN nodata
            dst.write(lia.astype("float32"), 1)
        v = lia[np.isfinite(lia)]                                 # quick stats on valid pixels
        print(f"  wrote {out_fp}: {v.size} valid px, range {v.min():.1f}-{v.max():.1f}, mean {v.mean():.1f} deg")
