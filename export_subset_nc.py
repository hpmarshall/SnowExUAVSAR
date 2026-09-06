"""Export the snow-modeling-team subset NetCDF: the 6 selected MCS pairs + 2 lidar surveys.

Contents (Mores Creek Summit only -- the source Zarr holds no other site):
  - pairs, in rank order (rank_dswe_pairs.py): 8, 2, 1, 5, 6, then 4 as the
    zero-change example; a `rank_note` coordinate labels each;
  - per pair: dswe (4 pols), coherence (4 pols), unwrapped_phase (4 pols),
    atm_delay_diff -- enough to use and quality-screen the retrieval;
  - site layers: sd for the two QSI surveys (2020-02-09, 2021-03-15), dem
    (2023-02-09 lidar DTM), lia (both headings);
  - full provenance in global attrs (CRS/grid, density, anchor, sign convention,
    deramp, custom unwrapping, companion documents).

Output: ZARR/mcs_top_pairs_subset.nc (NetCDF4/zlib, 256x256 chunks; git-ignored,
regenerable). Ships with mcs_zarr_userguide.pdf and mcs_summary_report.pdf.

Run:  PROJ_DATA=/opt/anaconda3/envs/myenv/share/proj \
      /opt/anaconda3/envs/myenv/bin/python export_subset_nc.py
"""

import numpy as np
import xarray as xr

ZARR_PATH = "ZARR/mores_creek_summit.zarr"
OUT = "ZARR/mcs_top_pairs_subset.nc"
PAIRS = [8, 2, 1, 5, 6, 4]                                       # rank order, then the bonus example
RANK_NOTE = ["rank 1", "rank 2", "rank 3", "rank 4", "rank 5", "zero-change example"]
SD_DATES = ["2020-02-09", "2021-03-15"]                          # the two MCS QSI surveys
DEM_DATE = "2023-02-09"

if __name__ == "__main__":
    src = xr.open_zarr(ZARR_PATH)

    out = xr.Dataset()
    sel = src.isel(pair=PAIRS)                                   # rank order preserved by isel
    out["dswe"] = sel["dswe"]
    out["coherence"] = sel["coherence"]
    out["unwrapped_phase"] = sel["unwrapped_phase"]
    out["atm_delay_diff"] = sel["atm_delay_diff"].astype("float32")
    out["sd"] = src["sd"].sel(sd_time=SD_DATES)                  # only the two MCS QSI surveys
    out["dem"] = src["dtm"].sel(lidar_time=DEM_DATE).squeeze(drop=True)
    out["lia"] = src["lia"].astype("float32")

    # string/date coordinates in MATLAB-friendly form (same convention as the GUI mirror)
    out = out.assign_coords(
        rank_note=("pair", np.array(RANK_NOTE)),
        pair_date1=("pair", np.array([str(v)[:10] for v in sel["pair_date1"].values])),
        pair_date2=("pair", np.array([str(v)[:10] for v in sel["pair_date2"].values])),
        line_id=("pair", sel["line_id"].values.astype(str)),
        heading=("pair", sel["heading"].values.astype("int32")),
        temporal_baseline_days=("pair", sel["temporal_baseline_days"].values.astype("int32")),
        sd_time=np.array(SD_DATES),
        sd_source=("sd_time", src["sd_source"].sel(sd_time=SD_DATES).values.astype(str)))

    out.attrs.update(
        title="Mores Creek Summit UAVSAR dSWE subset for snow-model precipitation-pattern testing",
        site="Mores Creek Summit, Idaho", crs="EPSG:32611",
        grid="3 m; transform (3,0,601558 / 0,-3,4870872.5); 2801 rows x 2624 cols",
        source_zarr=ZARR_PATH,
        pair_selection="ranked by dSWE-vs-lidar pattern correlation (rank_dswe_pairs.py); "
                       "pair index order in this file: " + ", ".join(map(str, PAIRS)) +
                       " (indices in the master Zarr); see rank_note coordinate",
        dswe_units="m w.e.; SWE change date1 -> date2",
        dswe_method="Guneriussen-2001 inversion of atmosphere-corrected unwrapped phase; "
                    "density 250 kg/m3; local incidence angle from JPL look vectors (lia); "
                    "masked where coherence < 0.35 or cos(LIA) <= 0.05; each map anchored so "
                    "the 90 m window at Mores Creek Summit SNOTEL (637:ID:SNTL) matches the "
                    "station interval dSWE",
        sign_convention="positive dswe = SWE gain over the interval",
        deramp="master-Zarr pair 8 (rank 1 here): planar ramp a+bx+cy removed (fit jointly "
               "with an elevation term that was kept)",
        custom_unwrapping="none of the pairs in this file required it (all carry ASF .unw.grd)",
        atmospheric_correction="ERA5 via phase_o_matic, delay integrated below aircraft "
                               "altitude, differenced per pair (atm_delay_diff, radians); "
                               "already subtracted in dswe, provided for transparency",
        lidar="sd: QSI snow depth (m), 0.5 m block-averaged to 3 m; dem: 2023-02-09 MCS "
              "lidar DTM (m), the surface used for incidence angle and atmospheric delay",
        companion_documents="mcs_zarr_userguide.pdf, mcs_summary_report.pdf",
        contact="HP Marshall, Boise State University")

    cy, cx = 256, 256
    enc = {v: {"zlib": True, "complevel": 4,
               "chunksizes": tuple(1 for _ in out[v].dims[:-2]) + (cy, cx)}
           for v in out.data_vars}
    out.to_netcdf(OUT, engine="netcdf4", encoding=enc)
    print(f"wrote {OUT}")

    # ---- verification: reopen and bit-compare representative slices against the Zarr ----
    chk = xr.open_dataset(OUT)
    assert list(chk["rank_note"].values) == RANK_NOTE
    assert list(chk["sd_time"].values) == SD_DATES and chk.sizes["pair"] == len(PAIRS)
    sl = dict(y=slice(1000, 1200), x=slice(1000, 1200))
    for v, zv, kw in [("dswe", "dswe", dict(pair=8, pol=0)), ("coherence", "coherence", dict(pair=2, pol=3)),
                      ("unwrapped_phase", "unwrapped_phase", dict(pair=6, pol=1)),
                      ("atm_delay_diff", "atm_delay_diff", dict(pair=4))]:
        i = PAIRS.index(kw.pop("pair"))
        a = chk[v].isel(pair=i, **kw, **sl).values
        b = src[zv].isel(pair=PAIRS[i], **kw, **sl).values.astype("float32")
        assert np.allclose(a, b, equal_nan=True), f"{v} mismatch"
    for j, d in enumerate(SD_DATES):
        assert np.allclose(chk["sd"].isel(sd_time=j, **sl).values,
                           src["sd"].sel(sd_time=d).isel(**sl).values, equal_nan=True)
    import os
    print(f"verified; size {os.path.getsize(OUT)/1e9:.2f} GB")
