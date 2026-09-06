"""Build the Mores Creek Summit lidar + UAVSAR InSAR pilot datacube.

Regrids all local MCS lidar products (0.5m native, block-averaged to 3m) and
every UAVSAR interferometric pair overlapping the site (all 4 polarizations;
wrapped phase regridded in the complex domain to respect phase wrapping) onto
one shared 3m grid, and writes the result to ZARR/mores_creek_summit.zarr.
"""

import glob                    # find MCS_Lidar files per date/product
import json                    # pair manifest
import netrc                   # Earthdata credentials
import os                      # paths
import time                    # retry backoff
import numpy as np             # array math
import rasterio                # raster IO + reprojection
import xarray as xr            # datacube assembly + zarr output
from rasterio.transform import Affine        # site grid definition
from rasterio.warp import reproject, Resampling  # regridding
import asf_search as asf                      # authenticated UAVSAR downloads
from uavsar_pytools.convert.tiff_conversion import grd_tiff_convert  # .grd -> GeoTIFF

# ---- canonical 3 m site grid (from the MCS lidar's native 0.5 m grid, EPSG:32611) ----
GRID_CRS = "EPSG:32611"
GRID_TRANSFORM = Affine(3.0, 0, 601558.0, 0, -3.0, 4870872.5)
GRID_WIDTH, GRID_HEIGHT = 2624, 2801                             # 0.5m 15747x16810 // 6, trimmed
GRID_Y = GRID_TRANSFORM.f + GRID_TRANSFORM.e * (np.arange(GRID_HEIGHT) + 0.5)  # cell-center coords
GRID_X = GRID_TRANSFORM.c + GRID_TRANSFORM.a * (np.arange(GRID_WIDTH) + 0.5)

LIDAR_ROOT = "LIDAR/CMR"
UAVSAR_WORK = "UAVSAR_INSAR"                                     # scratch dir, cleaned per pair
ZARR_OUT = "ZARR/mores_creek_summit.zarr"


def regrid_band(path, resampling, band=1):
    """Reproject one band from a local raster onto the shared site grid, NaN outside coverage."""
    with rasterio.open(path) as src:
        arr = src.read(band)
        if src.nodata is not None and not np.iscomplexobj(arr):  # mask real nodata sentinels (e.g. DTM/DSM -3.4e38)
            arr = np.where(arr == src.nodata, np.nan, arr)
        if np.iscomplexobj(arr):                                # wrapped phase: resample Re/Im independently
            re = np.full((GRID_HEIGHT, GRID_WIDTH), np.nan, dtype="float32")
            im = np.full((GRID_HEIGHT, GRID_WIDTH), np.nan, dtype="float32")
            reproject(np.nan_to_num(arr.real, nan=np.nan), re, src_transform=src.transform, src_crs=src.crs,
                      dst_transform=GRID_TRANSFORM, dst_crs=GRID_CRS, resampling=resampling,
                      src_nodata=np.nan, dst_nodata=np.nan)
            reproject(np.nan_to_num(arr.imag, nan=np.nan), im, src_transform=src.transform, src_crs=src.crs,
                      dst_transform=GRID_TRANSFORM, dst_crs=GRID_CRS, resampling=resampling,
                      src_nodata=np.nan, dst_nodata=np.nan)
            return re + 1j * im
        dst = np.full((GRID_HEIGHT, GRID_WIDTH), np.nan, dtype="float32")
        reproject(arr.astype("float32"), dst, src_transform=src.transform, src_crs=src.crs,
                  dst_transform=GRID_TRANSFORM, dst_crs=GRID_CRS, resampling=resampling,
                  src_nodata=np.nan, dst_nodata=np.nan)
        return dst


def build_lidar_dataset():
    """Regrid every MCS lidar acquisition (0.5m -> 3m, block-average) into an xarray Dataset."""
    qsi = {"2020-02-09": "SNEX20_QSI_SD/SNEX20_QSI_SD_0.5M_USIDMC_20200209_20200209.tif",
           "2021-03-15": "SNEX20_QSI_SD/SNEX20_QSI_SD_0.5M_USIDMC_20210315_20210315.tif"}
    mcs_dates = ["20220217", "20220317", "20220407", "20230209", "20230316", "20230405",
                 "20240115", "20240213", "20240315", "20240418", "20250113", "20250129", "20250404"]

    sd_stack, sd_dates, sd_source = [], [], []                   # 'sd' spans both QSI and MCS_Lidar
    for date, rel in sorted(qsi.items()):                        # QSI dates first
        print(f"  lidar sd (QSI) {date}")
        sd_stack.append(regrid_band(os.path.join(LIDAR_ROOT, rel), Resampling.average))
        sd_dates.append(date); sd_source.append("SNEX20_QSI_SD")

    dtm_stack, dsm_stack, chm_stack = [], [], []                 # only available from MCS_Lidar
    mcs_iso = [f"{d[:4]}-{d[4:6]}-{d[6:]}" for d in mcs_dates]
    for d, iso in zip(mcs_dates, mcs_iso):
        print(f"  lidar (MCS_Lidar) {iso}")
        base = os.path.join(LIDAR_ROOT, "SNEX_MCS_Lidar", f"SNEX_MCS_Lidar_{d}")
        sd_stack.append(regrid_band(f"{base}_SD_V01.0.tif", Resampling.average))
        sd_dates.append(iso); sd_source.append("SNEX_MCS_Lidar")
        dtm_stack.append(regrid_band(f"{base}_DTM_V01.0.tif", Resampling.average))
        dsm_stack.append(regrid_band(f"{base}_DSM_V01.0.tif", Resampling.average))
        chm_stack.append(regrid_band(f"{base}_CHM_V01.0.tif", Resampling.average))

    ds = xr.Dataset(
        {"sd": (("sd_time", "y", "x"), np.stack(sd_stack)),
         "dtm": (("lidar_time", "y", "x"), np.stack(dtm_stack)),
         "dsm": (("lidar_time", "y", "x"), np.stack(dsm_stack)),
         "chm": (("lidar_time", "y", "x"), np.stack(chm_stack))},
        coords={"y": GRID_Y, "x": GRID_X,
                "sd_time": np.array(sd_dates, dtype="datetime64[ns]"),
                "sd_source": ("sd_time", sd_source),
                "lidar_time": np.array(mcs_iso, dtype="datetime64[ns]")})
    return ds


def download_zip(url, dest, retries=4):
    """Resumable authenticated download of one InSAR pair bundle (Range-resume pattern).

    Verifies the final size against the server's own Content-Length before declaring success --
    a mid-transfer connection break (IncompleteRead) can otherwise leave a file that looks
    plausible but has a corrupt/missing zip central directory."""
    s = asf.ASFSession()
    u, _, p = netrc.netrc().authenticators("urs.earthdata.nasa.gov")
    s.auth_with_creds(u, p)
    expected = int(s.head(url, allow_redirects=True).headers.get("content-length", 0)) or None

    if os.path.exists(dest) and expected and os.path.getsize(dest) == expected:
        return                                                    # already complete

    for attempt in range(retries):
        try:
            offset = os.path.getsize(dest) if os.path.exists(dest) else 0
            if expected and offset > expected:                     # a prior corrupt run overshot -- start clean
                offset = 0
            headers = {"Range": f"bytes={offset}-"} if offset else {}
            r = s.get(url, headers=headers, stream=True, timeout=120)
            mode = "ab" if offset else "wb"
            with open(dest, mode) as f:
                for chunk in r.iter_content(chunk_size=1 << 20):
                    f.write(chunk)
            got = os.path.getsize(dest)
            if expected and got != expected:                       # silently-truncated response, no exception raised
                raise IOError(f"incomplete download: got {got}, expected {expected}")
            return
        except Exception as e:
            print(f"    retry {attempt + 1}/{retries} downloading {os.path.basename(dest)}: {e}")
            time.sleep(2 ** attempt)
    raise IOError(f"failed to download {url}")


PAIR_CACHE = os.path.join(UAVSAR_WORK, "cache")                 # cropped/regridded results, one file per pair


def process_pair(rec):
    """Download one InSAR pair, convert+crop+regrid all 4 pols' phase/coherence, clean up scratch files.

    Caches the (small, already-cropped) result to disk so a killed/restarted run never redoes a
    pair it already finished -- only the download+conversion of an in-progress pair is repeated."""
    scene = rec["sceneName"]
    cache_path = os.path.join(PAIR_CACHE, f"{scene}.npz")
    if os.path.exists(cache_path):
        print(f"pair {scene}: cached, skipping")
        d = np.load(cache_path)
        return {"wrapped_phase": d["wrapped_phase"], "unwrapped_phase": d["unwrapped_phase"],
                "coherence": d["coherence"]}

    zip_path = os.path.join(UAVSAR_WORK, f"{scene}.zip")
    print(f"pair {scene} ({rec['date1']} -> {rec['date2']}, line {rec['line_id']})")
    download_zip(rec["url"], zip_path)

    import zipfile
    zf = zipfile.ZipFile(zip_path)
    extract_dir = os.path.join(UAVSAR_WORK, scene)
    os.makedirs(extract_dir, exist_ok=True)
    tiff_dir = os.path.join(extract_dir, "tiffs")
    os.makedirs(tiff_dir, exist_ok=True)

    wrapped, unwrapped, coherence = {}, {}, {}
    for pol in ["HH", "HV", "VH", "VV"]:
        names = [n for n in zf.namelist() if f"L090{pol}_01" in n and n.endswith((".int.grd", ".unw.grd", ".cor.grd", ".ann"))]
        for n in names:
            zf.extract(n, extract_dir)
        for ext, store in [("int", wrapped), ("unw", unwrapped), ("cor", coherence)]:
            cands = [os.path.join(extract_dir, n) for n in names if n.endswith(f".{ext}.grd")]
            if not cands:                                          # e.g. unwrapping failed/wasn't produced for this pair
                print(f"    no .{ext}.grd for {pol} -- filling NaN")
                fill_dtype = "complex64" if ext == "int" else "float32"
                store[pol] = np.full((GRID_HEIGHT, GRID_WIDTH), np.nan, dtype=fill_dtype)
                continue
            grd = cands[0]
            grd_tiff_convert(grd, tiff_dir, overwrite=True)
            tif = os.path.join(tiff_dir, os.path.basename(grd) + ".tiff")
            resampling = Resampling.average if ext == "int" else Resampling.bilinear
            store[pol] = regrid_band(tif, resampling)
            os.remove(tif)                                       # reclaim disk immediately
        for n in names:                                          # drop the raw .grd/.ann now that it's converted
            os.remove(os.path.join(extract_dir, n))

    os.rmdir(tiff_dir); os.rmdir(extract_dir)
    os.remove(zip_path)                                          # done with the ~5GB bundle

    pols = ["HH", "HV", "VH", "VV"]
    result = {"wrapped_phase": np.stack([wrapped[p] for p in pols]),
              "unwrapped_phase": np.stack([unwrapped[p] for p in pols]),
              "coherence": np.stack([coherence[p] for p in pols])}
    os.makedirs(PAIR_CACHE, exist_ok=True)
    np.savez_compressed(cache_path, **result)                    # persist before this pair can be lost to a kill
    return result


def build_uavsar_dataset(manifest):
    """Process every pair in the manifest into a (pair, pol, y, x) xarray Dataset.

    A pair that fails for an unanticipated reason is skipped (logged, not fatal) so one bad
    pair can't take down a run that's otherwise made hours of progress on the others."""
    wrapped_all, unwrapped_all, coherence_all, kept = [], [], [], []
    for rec in manifest:
        try:
            result = process_pair(rec)
        except Exception as e:
            print(f"  SKIPPING {rec['sceneName']}: {e}")
            continue
        wrapped_all.append(result["wrapped_phase"])
        unwrapped_all.append(result["unwrapped_phase"])
        coherence_all.append(result["coherence"])
        kept.append(rec)
    manifest = kept

    date1 = np.array([r["date1"] for r in manifest], dtype="datetime64[ns]")
    date2 = np.array([r["date2"] for r in manifest], dtype="datetime64[ns]")
    baseline = ((date2 - date1) / np.timedelta64(1, "D")).astype("int32")
    ds = xr.Dataset(
        {"wrapped_phase": (("pair", "pol", "y", "x"), np.stack(wrapped_all)),
         "unwrapped_phase": (("pair", "pol", "y", "x"), np.stack(unwrapped_all)),
         "coherence": (("pair", "pol", "y", "x"), np.stack(coherence_all))},
        coords={"y": GRID_Y, "x": GRID_X, "pol": ["HH", "HV", "VH", "VV"],
                "pair_date1": ("pair", date1), "pair_date2": ("pair", date2),
                "line_id": ("pair", [r["line_id"] for r in manifest]),
                "heading": ("pair", [r["heading"] for r in manifest]),
                "temporal_baseline_days": ("pair", baseline)})
    return ds


if __name__ == "__main__":
    os.makedirs(UAVSAR_WORK, exist_ok=True)
    os.makedirs("ZARR", exist_ok=True)

    print("=== lidar ===")
    lidar_cache = os.path.join("UAVSAR_INSAR", "lidar_ds.nc")    # cache: avoid redoing 15 reprojections on restart
    if os.path.exists(lidar_cache):
        print("lidar dataset cached, loading")
        lidar_ds = xr.load_dataset(lidar_cache)
    else:
        lidar_ds = build_lidar_dataset()
        lidar_ds.to_netcdf(lidar_cache)

    print("=== UAVSAR ===")
    manifest = json.load(open("uavsar_pair_manifest_mcs.json"))
    uavsar_ds = build_uavsar_dataset(manifest)

    ds = xr.merge([lidar_ds, uavsar_ds])
    ds.attrs["crs"] = GRID_CRS
    ds.attrs["site"] = "Mores Creek Summit"
    chunks = {"y": 512, "x": 512, "sd_time": 1, "lidar_time": 1, "pair": 1, "pol": 1}
    ds.chunk(chunks).to_zarr(ZARR_OUT, mode="w")
    print(f"wrote {ZARR_OUT}")
