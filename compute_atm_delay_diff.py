"""Per-pair atmospheric phase delay difference, directly comparable to 'unwrapped_phase':
atm_delay_diff = atm_delay(date2) - atm_delay(date1), for each of the 22 cataloged InSAR pairs.
"""

import json                    # pair manifest (line_id/date1/date2, same order as the store's 'pair' dim)
import numpy as np             # array stacking
import pandas as pd            # date string formatting to match atm_delay's flight_date coord
import xarray as xr            # datacube I/O

ZARR_PATH = "ZARR/mores_creek_summit.zarr"
MANIFEST_PATH = "uavsar_pair_manifest_mcs.json"

if __name__ == "__main__":
    ds = xr.open_zarr(ZARR_PATH)
    atm_delay = ds["atm_delay"].load()                            # 20 flights, small enough to hold in memory
    flight_line_id = atm_delay["flight_line_id"].values
    flight_date = atm_delay["flight_date"].values                 # already 'YYYY-MM-DD' strings (Track B Step 1)
    flight_index = {(l, d): i for i, (l, d) in enumerate(zip(flight_line_id, flight_date))}

    manifest = json.load(open(MANIFEST_PATH))                     # same order as the store's existing 'pair' dim
    diffs = []
    for rec in manifest:
        i1 = flight_index[(rec["line_id"], rec["date1"])]
        i2 = flight_index[(rec["line_id"], rec["date2"])]
        diffs.append((atm_delay.isel(flight=i2) - atm_delay.isel(flight=i1)).values)
        print(f"  {rec['line_id']} {rec['date1']} -> {rec['date2']}: "
              f"median {np.nanmedian(diffs[-1]):.3f} rad, std {np.nanstd(diffs[-1]):.3f} rad")

    atm_delay_diff = xr.DataArray(np.stack(diffs, axis=0), dims=("pair", "y", "x"),
                                   coords={"y": ds["y"], "x": ds["x"]}, name="atm_delay_diff",
                                   attrs={"units": "radians",
                                          "long_name": "Atmospheric delay difference (date2 - date1), "
                                                        "directly comparable to unwrapped_phase"})
    xr.Dataset({"atm_delay_diff": atm_delay_diff}).to_zarr(ZARR_PATH, mode="a")
    print(f"Stored 'atm_delay_diff' ({len(diffs)} pairs) in {ZARR_PATH}")
