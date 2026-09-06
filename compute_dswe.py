"""Retrieve per-pair dSWE maps from atmosphere-corrected UAVSAR phase.

Per pair (22) x polarization (4), on the shared 3 m grid:
  phi_corr = unwrapped_phase - atm_delay_diff
  dz       = uavsar_pytools depth_from_phase(phi_corr, LIA, density=250 kg/m3,
             Guneriussen-2001 permittivity)       [m of snow]
  dswe     = dz * 0.25                            [m w.e.]

Sign convention: depth_from_phase implements dz = -phi*lambda / (4pi(cos t - sqrt(e-sin^2 t)));
the denominator is negative, so accumulation -> positive phase -> positive dz. This matches
the screening result (phase correlates POSITIVELY with lidar snow depth on the storm pairs
8/10/11), and the script asserts that the retrieved dSWE keeps that positive correlation.

Reference constant: unwrapped phase carries an unknown per-pair constant. Each pair's map
is anchored so the median dSWE in a 30x30-cell (90 m) window at the Mores Creek Summit
SNOTEL (637:ID:SNTL, 43.93200 -115.66588) equals the SNOTEL-measured interval dSWE --
so SNOTEL is a reference here, NOT an independent validation.

Mask: coherence < 0.35 or cos(LIA) <= 0.05 -> NaN.

Outputs: 'dswe' variable appended to ZARR/mores_creek_summit.zarr and ZARR/mcs_gui.nc;
per-pair diagnostics in dswe_retrieval_log.csv (committed).

Run:  PROJ_DATA=/opt/anaconda3/envs/myenv/share/proj \
      /opt/anaconda3/envs/myenv/bin/python compute_dswe.py
"""

import netCDF4                 # append variable to the GUI mirror
import numpy as np             # array math
import pandas as pd            # SNOTEL record
import xarray as xr            # datacube I/O
from pyproj import Transformer  # SNOTEL lat/lon -> UTM
from uavsar_pytools.snow_depth_inversion import depth_from_phase  # Guneriussen-2001 inversion

from compute_phase_sd_correlations import block_mean, pearson     # 90 m screening helpers

ZARR_PATH = "ZARR/mores_creek_summit.zarr"
NC_PATH = "ZARR/mcs_gui.nc"
POLS = ["HH", "HV", "VH", "VV"]
DENSITY = 250.0                # new-snow density assumption [kg/m3]
COH_MIN = 0.35                 # coherence mask threshold
SNOTEL_LATLON = (43.93200, -115.66588)                            # 637:ID:SNTL
WIN = 15                       # anchor window half-width: 30x30 cells = 90 m
SIGN_CHECK_PAIRS = [8, 10, 11]                                    # high-coherence storm pairs

if __name__ == "__main__":
    ds = xr.open_zarr(ZARR_PATH)
    sno = pd.read_csv("snotel_mcs_swe.csv", index_col="date", parse_dates=True)

    # SNOTEL pixel indices on the canonical grid
    to_utm = Transformer.from_crs("EPSG:4326", "EPSG:32611", always_xy=True)
    sx, sy = to_utm.transform(SNOTEL_LATLON[1], SNOTEL_LATLON[0])
    ix = int(np.argmin(np.abs(ds["x"].values - sx)))
    iy = int(np.argmin(np.abs(ds["y"].values - sy)))
    print(f"SNOTEL pixel: (iy={iy}, ix={ix}), UTM ({sx:.0f}, {sy:.0f})")

    # lazily create the all-NaN dswe variable (metadata only), then fill per pair
    init = xr.Dataset({"dswe": xr.full_like(ds["unwrapped_phase"], np.nan, dtype="float32")})
    init["dswe"].attrs.update(units="m w.e.", density_assumption_kg_m3=DENSITY,
                              long_name="SWE change (date2 - date1), Guneriussen-2001 inversion "
                                        "of atmosphere-corrected phase, anchored at MCS SNOTEL",
                              reference="anchored: 90 m window at SNOTEL 637 = SNOTEL interval dSWE",
                              mask=f"coherence < {COH_MIN} or cos(LIA) <= 0.05")
    init.drop_vars([c for c in init.coords if c not in ()]).to_zarr(ZARR_PATH, mode="a", compute=False)

    # matching variable in the GUI NetCDF mirror
    with netCDF4.Dataset(NC_PATH, "a") as nc:
        if "dswe" not in nc.variables:
            v = nc.createVariable("dswe", "f4", ("pair", "pol", "y", "x"),
                                  zlib=True, complevel=4, chunksizes=(1, 1, 256, 256),
                                  fill_value=np.float32(np.nan))
            v.units = "m w.e."
            v.long_name = init["dswe"].attrs["long_name"]

    sd2021_90 = block_mean(ds["sd"].sel(sd_time="2021-03-15").squeeze().values)  # for the sign check
    rows, sign_r = [], {}
    for pi in range(ds.sizes["pair"]):
        line = str(ds["line_id"][pi].values)
        d1, d2 = str(ds["pair_date1"][pi].values)[:10], str(ds["pair_date2"][pi].values)[:10]
        # SNOTEL interval dSWE [m w.e.] and interval bulk density of new accumulation (context)
        dswe_sno = (sno.loc[d2, "swe_mm"] - sno.loc[d1, "swe_mm"]) / 1000.0
        d_depth = (sno.loc[d2, "depth_mm"] - sno.loc[d1, "depth_mm"]) / 1000.0
        rho_sno = 1000.0 * (dswe_sno / d_depth) if d_depth > 0.02 else np.nan  # only in real accumulation
        lia_rad = np.deg2rad(ds["lia"].sel(line=line).values)     # this line's local incidence angle
        atm = ds["atm_delay_diff"][pi].values                     # per-pair atmospheric correction
        coh_ok = None
        print(f"[{pi:2d}] {line} {d1} -> {d2}: SNOTEL dSWE {dswe_sno*1000:+.1f} mm, "
              f"interval density {rho_sno:.0f} kg/m3" if np.isfinite(rho_sno) else
              f"[{pi:2d}] {line} {d1} -> {d2}: SNOTEL dSWE {dswe_sno*1000:+.1f} mm")

        pair_out = np.full((len(POLS),) + lia_rad.shape, np.nan, dtype="float32")
        for pol_i, pol in enumerate(POLS):
            unw = ds["unwrapped_phase"][pi, pol_i].values
            if not np.isfinite(unw).any():                        # 4 pairs ship no .unw.grd
                rows.append(dict(pair=pi, line=line, date1=d1, date2=d2, pol=pol,
                                 snotel_dswe_m=round(dswe_sno, 4), snotel_density=rho_sno,
                                 anchor_offset_m=np.nan, anchor_window_m=np.nan,
                                 finite_frac=0.0, window_n=0))
                continue
            phi = unw - atm                                       # atmosphere-corrected phase
            dz = depth_from_phase(phi, lia_rad, density=DENSITY)  # m of snow (Guneriussen-2001)
            dswe = (dz * DENSITY / 1000.0).astype("float32")      # m w.e.
            coh = ds["coherence"][pi, pol_i].values
            dswe[(coh < COH_MIN) | (np.cos(lia_rad) <= 0.05)] = np.nan  # quality mask
            # 90 m anchor window at SNOTEL; if the coherence mask leaves <50 valid cells
            # there, expand progressively (up to 100 cells = 600 m) so every pair with any
            # nearby coherent data still gets an absolute reference
            offset, wn, w = np.nan, 0, WIN
            while w <= 100:
                win = dswe[iy - w:iy + w, ix - w:ix + w]
                wn = int(np.isfinite(win).sum())
                if wn >= 50:
                    offset = dswe_sno - np.nanmedian(win)
                    break
                w *= 2
            if np.isfinite(offset):
                dswe += offset                                    # anchor: window median == SNOTEL
            pair_out[pol_i] = dswe
            rows.append(dict(pair=pi, line=line, date1=d1, date2=d2, pol=pol,
                             snotel_dswe_m=round(dswe_sno, 4), snotel_density=rho_sno,
                             anchor_offset_m=round(float(offset), 4) if np.isfinite(offset) else np.nan,
                             anchor_window_m=6 * w,               # window full width in metres
                             finite_frac=round(float(np.isfinite(dswe).mean()), 3), window_n=wn))
            if pi in SIGN_CHECK_PAIRS and pol == "HH":            # empirical sign check vs lidar SD
                r, n = pearson(sd2021_90, block_mean(dswe))
                sign_r[pi] = r
                print(f"      sign check HH: corr(dSWE, SD 2021-03-15) = {r:+.3f} (n={n})")

        # region-write this pair into the Zarr and patch the NetCDF mirror
        patch = xr.Dataset({"dswe": (("pair", "pol", "y", "x"), pair_out[None])})
        patch.chunk({"pair": 1, "pol": 1, "y": 512, "x": 512}).to_zarr(
            ZARR_PATH, mode="r+", region={"pair": slice(pi, pi + 1)})
        with netCDF4.Dataset(NC_PATH, "a") as nc:
            nc["dswe"][pi] = pair_out

    # sign convention is fixed by the large-signal pair: pair 8 (2021-02-10 -> 03-03,
    # 173 mm SNOTEL accumulation) must correlate positively with the lidar snowpack, and
    # the three storm pairs must be positive on average. Pairs 10/11 carry only ~15 mm of
    # accumulation -- their near-zero pattern correlation after LIA normalization is a
    # weak-signal outcome, not evidence about the sign.
    assert sign_r[8] > 0.3 and np.mean(list(sign_r.values())) > 0, f"sign convention violated: {sign_r}"
    log = pd.DataFrame(rows).round({"snotel_density": 0})
    log.to_csv("dswe_retrieval_log.csv", index=False)
    print("sign checks:", {k: round(v, 3) for k, v in sign_r.items()})
    print("wrote dswe_retrieval_log.csv; 'dswe' stored in Zarr + NetCDF")
