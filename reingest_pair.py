"""Re-ingest one InSAR pair into the existing MCS Zarr store, without a full rebuild.

Built for the 23205 2021-02-03 -> 02-10 pair, whose _L090_02 bundle was silently
skipped by the old hardcoded L090{pol}_01 filter in build_mcs_datacube.py (fixed in
commit 9dd99e3), leaving all-NaN slices in the store AND a poisoned all-NaN .npz
cache. A full rebuild would work but mode="w" wipes the separately-appended lia /
atm_delay / atm_delay_diff variables -- so this script re-processes just the one
pair and region-writes its slices into both the Zarr and the GUI's NetCDF mirror.

Run:  PROJ_DATA=/opt/anaconda3/envs/myenv/share/proj \
      /opt/anaconda3/envs/myenv/bin/python reingest_pair.py [sceneName]
"""

import json                    # pair manifest
import os                      # cache deletion
import sys                     # optional sceneName argument
import numpy as np             # amplitude/angle for the NetCDF mirror
import xarray as xr            # zarr region write
import netCDF4                 # in-place patch of ZARR/mcs_gui.nc

from build_mcs_datacube import process_pair, PAIR_CACHE, ZARR_OUT  # reuse the pipeline verbatim

NC_PATH = "ZARR/mcs_gui.nc"
DEFAULT_SCENE = "UA_lowman_23205_21009-004_21012-000_0007d_s01_L090_02"

if __name__ == "__main__":
    scene = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_SCENE

    manifest = json.load(open("uavsar_pair_manifest_mcs.json"))  # same order as the store's 'pair' dim
    idx = next(i for i, r in enumerate(manifest) if r["sceneName"] == scene)
    rec = manifest[idx]
    print(f"pair index {idx}: {scene} ({rec['date1']} -> {rec['date2']}, line {rec['line_id']})")

    cache_fp = os.path.join(PAIR_CACHE, f"{scene}.npz")          # poisoned all-NaN cache from the buggy run
    if os.path.exists(cache_fp):
        os.remove(cache_fp)
        print(f"removed stale cache {cache_fp}")

    result = process_pair(rec)                                   # download ~5GB bundle, convert, regrid, re-cache
    for name, arr in result.items():                             # refuse to write another all-NaN result
        frac = float(np.isfinite(arr if not np.iscomplexobj(arr) else arr.real).mean())
        print(f"  {name}: {frac:.1%} finite")
        assert frac > 0, f"{name} came back all-NaN -- fix did not take, not writing"

    # ---- region-write the three variables into the existing Zarr at this pair index ----
    ds = xr.open_zarr(ZARR_OUT)
    patch = xr.Dataset(
        {v: (("pair", "pol", "y", "x"), result[v][None]) for v in
         ["wrapped_phase", "unwrapped_phase", "coherence"]},
        coords={"y": ds["y"], "x": ds["x"], "pol": ds["pol"]})
    # region writes reject non-region coordinate variables -- drop them, keep pure dims
    patch = patch.drop_vars([c for c in patch.coords if c not in ("pair",)])
    patch.chunk({"pair": 1, "pol": 1, "y": 512, "x": 512}).to_zarr(
        ZARR_OUT, mode="r+", region={"pair": slice(idx, idx + 1)})
    print(f"zarr: wrote pair {idx} into {ZARR_OUT}")

    # ---- patch the GUI NetCDF mirror in place (avoids the ~25 min full re-export) ----
    # netCDF file stores (pair, pol, y, x) in C order, same as the arrays here
    with netCDF4.Dataset(NC_PATH, "a") as nc:
        nc["coherence"][idx] = result["coherence"]
        nc["unwrapped_phase"][idx] = result["unwrapped_phase"]
        nc["wrapped_amplitude"][idx] = np.abs(result["wrapped_phase"]).astype("float32")
        nc["wrapped_angle"][idx] = np.angle(result["wrapped_phase"]).astype("float32")
    print(f"netcdf: patched pair {idx} in {NC_PATH}")

    # ---- verify both stores: this pair now valid, a neighbor pair untouched ----
    ds = xr.open_zarr(ZARR_OUT)                                  # reopen post-write
    nc = xr.open_dataset(NC_PATH)
    for v in ["coherence", "unwrapped_phase"]:
        zfrac = float(np.isfinite(ds[v].isel(pair=idx)).mean())
        print(f"  verify zarr {v}: {zfrac:.1%} finite")
        assert zfrac > 0
        match = np.allclose(ds[v].isel(pair=idx).values, nc[v].isel(pair=idx).values, equal_nan=True)
        print(f"  verify nc {v} == zarr: {match}")
        assert match
    other = int(np.isfinite(ds["coherence"].isel(pair=idx - 1, pol=0, y=slice(0, 512), x=slice(0, 512))).sum())
    print(f"  neighbor pair {idx-1} coherence sample still has {other} finite px (untouched)")
    assert other > 0
    assert "lia" in ds and "atm_delay" in ds and "atm_delay_diff" in ds  # appended vars survived
    print("re-ingest complete and verified")
