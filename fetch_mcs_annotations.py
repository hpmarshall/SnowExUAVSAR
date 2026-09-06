"""Fetch map annotations for the MCS GUI: Idaho Highway 21 and the two named peaks.

Queries OpenStreetMap's Overpass API for the ID 21 polyline (Ponderosa Pine Scenic
Route) and the Pilot Peak / Freeman Peak summit nodes, reprojects everything to the
datacube's UTM grid (EPSG:32611), clips to the grid bounds, and writes
mcs_annotations.json -- small enough to commit, so the GUI never needs the network.

Run:  PROJ_DATA=/opt/anaconda3/envs/myenv/share/proj \
      /opt/anaconda3/envs/myenv/bin/python fetch_mcs_annotations.py
"""

import json                    # output file
import numpy as np             # peak sanity-check against the DEM
import requests                # Overpass API
import xarray as xr            # DEM for the sanity check
from pyproj import Transformer  # lat/lon -> UTM

from build_mcs_datacube import GRID_CRS, GRID_TRANSFORM, GRID_WIDTH, GRID_HEIGHT  # canonical grid

OVERPASS = "https://overpass-api.de/api/interpreter"
OUT_PATH = "mcs_annotations.json"
ZARR_PATH = "ZARR/mores_creek_summit.zarr"

# canonical grid bounds in UTM (same derivation as compute_local_incidence_angle.py)
X0, X1 = GRID_TRANSFORM.c, GRID_TRANSFORM.c + GRID_TRANSFORM.a * GRID_WIDTH
Y0, Y1 = GRID_TRANSFORM.f + GRID_TRANSFORM.e * GRID_HEIGHT, GRID_TRANSFORM.f

QUERY = """
[out:json][timeout:50];
(
  way["ref"="ID 21"](43.89,-115.76,44.00,-115.61);
  node["natural"="peak"]["name"~"Pilot Peak|Freeman Peak"](43.85,-115.85,44.05,-115.55);
);
out geom;
"""

if __name__ == "__main__":
    # Overpass 406's the default python-requests user agent -- identify the tool instead
    r = requests.post(OVERPASS, data={"data": QUERY}, timeout=90,
                      headers={"User-Agent": "SnowExUAVSAR-MCS-annotations/1.0"})
    r.raise_for_status()
    elements = r.json()["elements"]
    to_utm = Transformer.from_crs("EPSG:4326", GRID_CRS, always_xy=True)

    # ---- Highway 21: one segment per OSM way, clipped to the grid bounds ----
    segments = []
    for e in elements:
        if e["type"] != "way":
            continue
        lon = [p["lon"] for p in e["geometry"]]
        lat = [p["lat"] for p in e["geometry"]]
        x, y = to_utm.transform(lon, lat)
        seg = [[round(xi, 1), round(yi, 1)] for xi, yi in zip(x, y)
               if X0 <= xi <= X1 and Y0 <= yi <= Y1]              # keep only in-grid vertices
        if len(seg) >= 2:
            segments.append(seg)
    n_vertices = sum(len(s) for s in segments)
    print(f"Highway 21: {len(segments)} in-grid segments, {n_vertices} vertices")

    # ---- Peaks: named summit nodes, sanity-checked against the DEM ----
    dem = xr.open_zarr(ZARR_PATH)["dtm"].sel(lidar_time="2023-02-09").squeeze()
    peaks = []
    for e in elements:
        if e["type"] != "node":
            continue
        name, ele = e["tags"]["name"], float(e["tags"].get("ele", "nan"))
        x, y = to_utm.transform(e["lon"], e["lat"])
        # DEM elevation in a 15-pixel (45 m) window around the node -- should be near the OSM elevation
        win = dem.sel(x=slice(x - 45, x + 45), y=slice(y + 45, y - 45)).values
        dem_max = float(np.nanmax(win)) if np.isfinite(win).any() else float("nan")
        print(f"{name}: OSM ele {ele:.0f} m, DEM max within 45 m = {dem_max:.1f} m")
        assert abs(dem_max - ele) < 30, f"{name}: DEM disagrees with OSM elevation by >30 m"
        peaks.append({"name": name, "x": round(x, 1), "y": round(y, 1), "ele_m": ele})

    json.dump({"crs": GRID_CRS,
               "highway21": {"name": "Idaho Highway 21", "segments": segments},
               "peaks": peaks},
              open(OUT_PATH, "w"), indent=1)
    print(f"wrote {OUT_PATH}")
