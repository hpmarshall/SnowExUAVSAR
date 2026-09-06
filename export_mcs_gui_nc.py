"""Export the MCS Zarr datacube to a NetCDF4 companion file for the MATLAB GUI.

MATLAB R2025a on this machine has no zarrread, and its ncread cannot touch complex
(<c8) arrays at all -- so this script mirrors ZARR/mores_creek_summit.zarr into
ZARR/mcs_gui.nc (HDF5-backed NetCDF4, chunked so ncread(...,start,count) pulls one
2D slice cheaply), splitting the complex wrapped phase into amplitude and angle.
The Zarr remains the single source of truth; the .nc is a regenerable view.

Run:  PROJ_DATA=/opt/anaconda3/envs/myenv/share/proj \
      /opt/anaconda3/envs/myenv/bin/python export_mcs_gui_nc.py [--decimate N]
"""

import argparse                # optional --decimate flag for a lighter file
import numpy as np             # amplitude/angle math
import xarray as xr            # zarr in, netcdf out

ZARR_PATH = "ZARR/mores_creek_summit.zarr"
NC_PATH = "ZARR/mcs_gui.nc"
DEM_DATE = "2023-02-09"        # same lidar DTM used by the LIA and atm-delay tracks


def date_strings(da):
    """datetime64 coordinate -> ISO 'YYYY-MM-DD' strings (MATLAB reads char arrays cleanly)."""
    return np.array([str(v)[:10] for v in da.values])


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--decimate", type=int, default=1, metavar="N",
                    help="keep every Nth pixel in y/x (default 1 = full 3 m resolution)")
    args = ap.parse_args()

    src = xr.open_zarr(ZARR_PATH)
    if args.decimate > 1:                                        # optional coarsening for a lighter file
        src = src.isel(y=slice(None, None, args.decimate), x=slice(None, None, args.decimate))

    out = xr.Dataset()
    out["dem"] = src["dtm"].sel(lidar_time=DEM_DATE).squeeze(drop=True)  # terrain panel base layer
    out["atm_delay"] = src["atm_delay"].astype("float32")        # per-flight delay (f8 -> f4, plenty for display)
    out["lia"] = src["lia"].astype("float32")                    # per-line local incidence angle (deg)
    out["coherence"] = src["coherence"]
    out["unwrapped_phase"] = src["unwrapped_phase"]
    wrapped = src["wrapped_phase"]                               # complex64 -- split for MATLAB
    out["wrapped_amplitude"] = np.abs(wrapped).astype("float32")
    out["wrapped_angle"] = xr.apply_ufunc(np.angle, wrapped, dask="parallelized",
                                          output_dtypes=[np.float32])
    out["atm_delay_diff"] = src["atm_delay_diff"].astype("float32")  # per-pair delay difference
    out["sd"] = src["sd"]                                        # snow depth (QSI + MCS_Lidar dates)
    out["dtm"] = src["dtm"]
    out["dsm"] = src["dsm"]
    out["chm"] = src["chm"]

    # numeric coordinates come along automatically; rewrite the datetime/string ones as
    # plain arrays MATLAB's ncread returns directly
    out = out.assign_coords(
        sd_time=date_strings(src["sd_time"]), lidar_time=date_strings(src["lidar_time"]),
        pair_date1=("pair", date_strings(src["pair_date1"])),
        pair_date2=("pair", date_strings(src["pair_date2"])),
        flight_date=("flight", src["flight_date"].values.astype(str)),
        flight_line_id=("flight", src["flight_line_id"].values.astype(str)),
        line_id=("pair", src["line_id"].values.astype(str)),
        sd_source=("sd_time", src["sd_source"].values.astype(str)),
        heading=("pair", src["heading"].values.astype("int32")),
        temporal_baseline_days=("pair", src["temporal_baseline_days"].values.astype("int32")))

    out.attrs.update(crs=src.attrs.get("crs", "EPSG:32611"), site=src.attrs.get("site", ""),
                     source_zarr=ZARR_PATH, dem_source=f"MCS lidar DTM {DEM_DATE}",
                     lia_units="degrees", phase_units="radians", atm_delay_units="radians")

    # per-variable compression + slice-friendly chunking; string coords get no encoding
    cy, cx = min(256, out.sizes["y"]), min(256, out.sizes["x"])
    enc = {}
    for name, da in out.data_vars.items():
        chunks = tuple(1 for _ in da.dims[:-2]) + (cy, cx)       # one 2D tile per chunk
        enc[name] = {"zlib": True, "complevel": 4, "chunksizes": chunks}
    out.to_netcdf(NC_PATH, engine="netcdf4", encoding=enc)
    print(f"wrote {NC_PATH}")
    for name, da in out.data_vars.items():
        print(f"  {name:18s} {da.dims} {da.dtype}")
