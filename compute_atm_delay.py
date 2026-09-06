"""Compute per-flight atmospheric phase delay for Mores Creek Summit UAVSAR flights,
restricted to the atmosphere below UAVSAR's flight altitude (~12.5 km), using ERA5
reanalysis via phase_o_matic.

phase_o_matic is built for satellites: its refractivity profile is referenced to the
top of the loaded ERA5 pressure-level stack (~48 km), and its own get_delay() samples
that profile at the ground -- giving delay for a signal traveling from the model top
down to the surface. UAVSAR only flies to ~12.5 km, so the delay it actually
accumulates is delay(ground) - delay(aircraft_altitude): refractivity is a cumulative
quantity relative to one fixed reference, so that reference cancels out in the
difference regardless of what it was. get_delay_below_aircraft() below implements
exactly that, reusing every other phase_o_matic step unmodified.
"""

import json                    # flight timestamp cache
import os                      # paths
import numpy as np             # array math
import pandas as pd            # timestamps
import xarray as xr            # datacube I/O
import rioxarray               # noqa: F401 -- registers the .rio accessor
import asf_search as asf       # precise flight acquisition times
import scipy.integrate         # phase_o_matic calls the pre-1.14 scipy.integrate.cumtrapz, removed upstream

if not hasattr(scipy.integrate, "cumtrapz"):
    def _cumtrapz_compat(y, x=None, initial=None, **kw):
        # phase_o_matic always passes initial=np.nan as a length-preserving placeholder,
        # then immediately discards that first element itself -- new scipy only allows
        # initial=None/0, and 0 is an equally fine placeholder since it's discarded anyway.
        if initial is not None:
            initial = 0
        return scipy.integrate.cumulative_trapezoid(y, x=x, initial=initial, **kw)
    scipy.integrate.cumtrapz = _cumtrapz_compat
from phase_o_matic.download import download_era
from phase_o_matic.preprocess import (get_vapor_partial_pressure, convert_pressure_to_pascals,
                                       interpolate_to_heights, geopotential_to_geopotential_heights)
from phase_o_matic.phase_delay import calculate_refractive_indexes
from glob import glob                                           # find a representative .ann per line
from uavsar_pytools.convert.tiff_conversion import read_annotation  # parse .ann key/value pairs
from build_mcs_datacube import GRID_CRS, GRID_TRANSFORM, GRID_WIDTH, GRID_HEIGHT, GRID_Y, GRID_X  # canonical grid

ZARR_PATH = "ZARR/mores_creek_summit.zarr"
ATM_WORK_DIR = "ATM_DELAY"                                      # ERA5 scratch (git-ignored)
FLIGHT_TIMES_CACHE = "atm_delay_flight_times.json"
UAVSAR_WAVELENGTH_M = 0.2379                                    # L-band, SnowEx UAVSAR tech reference Table 1
MCS_CENTER_POINT = "POINT(-115.685 43.946)"                     # for resolving each flight's COMPLEX scene


def unique_flights(manifest_path="uavsar_pair_manifest_mcs.json"):
    """Every unique (line_id, date) contributing to the 22 cataloged InSAR pairs."""
    manifest = json.load(open(manifest_path))
    flights = set()
    for r in manifest:
        flights.add((r["line_id"], r["date1"]))
        flights.add((r["line_id"], r["date2"]))
    return sorted(flights)


def load_flight_times():
    """Precise acquisition UTC timestamp per unique flight, cached (Step 1)."""
    if os.path.exists(FLIGHT_TIMES_CACHE):
        return json.load(open(FLIGHT_TIMES_CACHE))
    import datetime
    results = []
    for line_id, date in unique_flights():
        d = datetime.date.fromisoformat(date)
        window = (str(d), str(d + datetime.timedelta(days=1)))
        r = asf.geo_search(platform=asf.PLATFORM.UAVSAR, intersectsWith=MCS_CENTER_POINT,
                           start=window[0], end=window[1], processingLevel="COMPLEX")
        match = [p for p in r if line_id in p.properties["sceneName"]]
        if not match:
            print(f"  MISSING flight metadata: {line_id} {date}")
            continue
        results.append({"line_id": line_id, "date": date, "start_time_utc": match[0].properties["startTime"]})
    json.dump(results, open(FLIGHT_TIMES_CACHE, "w"), indent=1)
    return results


def load_reference_dem():
    """The same 2023-02-09 MCS DTM used for LIA, reprojected to lat/lon (Step 2)."""
    ds = xr.open_zarr(ZARR_PATH)
    dtm = ds["dtm"].sel(lidar_time="2023-02-09").squeeze().rio.write_crs("EPSG:32611")
    return dtm.rio.reproject("EPSG:4326").rename({"x": "longitude", "y": "latitude"})


def line_aircraft_altitude(line_id):
    """Global Average Altitude (m) for one line -- bit-identical across every flight of that
    line (confirmed exhaustively in Track A Step 0), so any one flight's .ann suffices."""
    ann_fp = sorted(glob(f"UAVSAR_LKV/lowman_{line_id}_*_L090HH_01_*.ann"))[0]
    return read_annotation(ann_fp)["global average altitude"]["value"]


def line_lia_radians(line_id, dem):
    """Track A's per-line LIA, reprojected from the canonical UTM grid onto this dem's lat/lon grid."""
    lia_deg = xr.open_zarr(ZARR_PATH)["lia"].sel(line=line_id).rio.write_crs(GRID_CRS)
    lia_deg = lia_deg.rio.reproject_match(dem.rio.write_crs("EPSG:4326"))
    return np.deg2rad(lia_deg).rename({"x": "longitude", "y": "latitude"})


def regrid_to_canonical(da_latlon):
    """Reproject a lat/lon result back onto the canonical 3 m UTM grid used by the rest of the datacube."""
    da_latlon = da_latlon.rio.write_crs("EPSG:4326")
    return da_latlon.rio.reproject(GRID_CRS, transform=GRID_TRANSFORM, shape=(GRID_HEIGHT, GRID_WIDTH))


def get_delay_below_aircraft(refractivity_ds, dem, inc, aircraft_alt_m, wavelength=UAVSAR_WAVELENGTH_M):
    """Delay accumulated only between the ground and aircraft_alt_m (see module docstring).

    Mirrors phase_o_matic.phase_delay.get_delay()'s DEM-sampling step, but additionally
    samples the same (full-height-range) refractivity profile at the constant aircraft
    altitude and subtracts it, instead of referencing the model top."""
    assert "N" in refractivity_ds.data_vars
    assert isinstance(inc, xr.DataArray) or inc < 2 * np.pi, "inc must be radians"

    # sample at the aircraft's cruise altitude (same value everywhere -- constant height)
    n_aircraft = refractivity_ds["N"].interp(height=aircraft_alt_m)
    n_aircraft = n_aircraft.interp(latitude=dem.latitude, longitude=dem.longitude)

    # sample at the ground surface, same fine height interpolation get_delay() itself uses
    ds = refractivity_ds.interp(height=np.round(np.arange(float(dem.min()) - 2, float(dem.max()) + 2)),
                                 method="linear")
    n_ground = ds["N"].interp(latitude=dem.latitude, longitude=dem.longitude, height=dem)

    cos_inc = np.cos(inc)
    delay = (n_ground - n_aircraft) * 4.0 * np.pi / (wavelength * cos_inc)
    # 1/cos(inc) diverges at grazing incidence and is unphysical past 90 deg (radar shadow /
    # back-facing slopes, LIA's own arccos legitimately returns >90 deg there) -- mask both out
    # the same way any real InSAR product would exclude shadow/layover terrain
    delay = delay.where(cos_inc > 0.05)
    delay.attrs.update(units="radians" if wavelength != 4 * np.pi else "meters",
                        long_name="Atmospheric delay below aircraft altitude")
    return delay.transpose(..., "latitude", "longitude")        # match dem's dim order


def _adapt_era5_naming(era):
    """CDS's current backend names fields differently than phase_o_matic (written against an
    older CDS API) expects: 'pressure_level' instead of 'level', 'valid_time' instead of
    'time', and reports pressure units as 'hPa' rather than the exact string 'millibars'
    the library checks for (numerically identical -- 1 hPa = 1 millibar)."""
    era = era.rename({"pressure_level": "level", "valid_time": "time"})
    if era["level"].attrs.get("units") == "hPa":
        era["level"].attrs["units"] = "millibars"
    return era


def refractivity_for_flight(start_time_utc, dem, aircraft_alt_m, work_dir=ATM_WORK_DIR):
    """ERA5 download + preprocessing for one flight's acquisition hour (reused from phase_o_matic)."""
    os.makedirs(work_dir, exist_ok=True)
    era_dir = os.path.join(work_dir, "ERA5")
    os.makedirs(era_dir, exist_ok=True)
    from shapely.geometry import box
    # ERA5's native grid is ~0.25 deg; our MCS crop (~0.1 x 0.08 deg) is smaller than one grid
    # cell, so CDS's area-crop finds no points and errors ("non-empty area crop/mask"). Pad the
    # requested area well past one grid spacing -- costs a handful of extra (free) grid points.
    pad = 0.5
    lon_min, lat_min = float(dem.longitude.min()), float(dem.latitude.min())
    lon_max, lat_max = float(dem.longitude.max()), float(dem.latitude.max())
    subset = box(lon_min - pad, lat_min - pad, lon_max + pad, lat_max + pad)
    era_fp = download_era(pd.Timestamp(start_time_utc), out_dir=era_dir, subset=subset,
                          humid_param="specific_humidity")
    era = xr.open_dataset(era_fp)
    era = _adapt_era5_naming(era)                                # bridge CDS's current field names to phase_o_matic's expected ones
    era = convert_pressure_to_pascals(era)
    era = get_vapor_partial_pressure(era)
    era = geopotential_to_geopotential_heights(era)
    # must span both the ground (for n_ground) and the aircraft's cruise altitude (for
    # n_aircraft) -- capping max_alt at the DEM's own (much lower) elevation left no data
    # near 12.5 km, so get_delay_below_aircraft()'s interp(height=aircraft_alt_m) returned NaN
    ds = interpolate_to_heights(era, min_alt=float(dem.min()), max_alt=aircraft_alt_m + 500)
    return calculate_refractive_indexes(ds)


RESULTS_DIR = os.path.join(ATM_WORK_DIR, "results")             # per-flight cache -- resumable if interrupted


def compute_flight_delay(flight, dem, lia_rad_by_line):
    """One flight's atm_delay on the canonical grid, cached to disk so a long CDS queue can be resumed."""
    out_fp = os.path.join(RESULTS_DIR, f"{flight['line_id']}_{flight['date']}.nc")
    if os.path.exists(out_fp):
        return xr.open_dataset(out_fp)["atm_delay"]
    aircraft_alt = line_aircraft_altitude(flight["line_id"])
    refr = refractivity_for_flight(flight["start_time_utc"], dem, aircraft_alt)
    delay = get_delay_below_aircraft(refr, dem, lia_rad_by_line[flight["line_id"]], aircraft_alt_m=aircraft_alt)
    delay = regrid_to_canonical(delay.squeeze("time", drop=True))
    # dem's own scalar 'lidar_time' coord (and lia's scalar 'line' coord) ride along through every
    # op above -- drop them so this doesn't collide with the store's existing 'lidar_time' dimension
    delay = delay.drop_vars([c for c in delay.coords if c not in ("x", "y")])
    delay.name = "atm_delay"
    os.makedirs(RESULTS_DIR, exist_ok=True)
    delay.to_netcdf(out_fp)                                      # cache before returning, in case a later flight fails
    return delay


if __name__ == "__main__":
    print("=== Step 1: flight timestamps ===")
    flight_times = load_flight_times()
    print(f"{len(flight_times)} unique flights resolved")

    print("=== Step 2: reference DEM ===")
    dem = load_reference_dem()
    print(f"DEM ready: {dem.shape}, elevation {float(dem.min()):.0f}-{float(dem.max()):.0f} m")

    if not os.path.exists(os.path.expanduser("~/.cdsapirc")):
        print("\n~/.cdsapirc not found yet -- stopping here (Steps 1-2 done, "
              "Step 3 needs CDS credentials, Step 4 needs Track A's lia).")
    else:
        print("=== Step 4: per-flight delay using Track A's real LIA ===")
        lia_rad_by_line = {line: line_lia_radians(line, dem) for line in ("05208", "23205")}
        delays, line_ids, dates = [], [], []
        for i, flight in enumerate(flight_times):
            print(f"  [{i+1}/{len(flight_times)}] {flight['line_id']} {flight['date']}")
            try:
                delays.append(compute_flight_delay(flight, dem, lia_rad_by_line))
                line_ids.append(flight["line_id"])
                dates.append(flight["date"])
            except Exception as e:
                print(f"    FAILED: {e}")                        # skip and continue -- rerun later to retry just this flight

        print("=== Step 5: store atm_delay in the Zarr store ===")
        atm_delay = xr.concat(delays, dim="flight").assign_coords(
            flight_line_id=("flight", line_ids), flight_date=("flight", dates))
        atm_delay.name = "atm_delay"
        atm_delay.attrs.update(units="radians", long_name="Atmospheric delay below aircraft altitude")
        xr.Dataset({"atm_delay": atm_delay}).to_zarr(ZARR_PATH, mode="a")
        print(f"Stored 'atm_delay' ({len(delays)} flights) in {ZARR_PATH}")
