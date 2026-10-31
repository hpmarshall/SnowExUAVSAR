"""Export the 3 m lidar DTM (2023-02-09 survey) from the MCS zarr to a GeoTIFF.
This is the terrain input the standalone LIA scripts (lia_pytools.py, lia_corrected.py)
read, so they need no zarr or repo imports and can be shared with the LKV files alone.

Run:  PROJ_DATA=/opt/anaconda3/envs/myenv/share/proj \
        /opt/anaconda3/envs/myenv/bin/python export_lia_inputs.py
"""

import os                                                         # output-dir creation
import numpy as np                                                # float32 cast / NaN nodata
import rasterio as rio                                            # GeoTIFF writer
import xarray as xr                                               # zarr reader
from build_mcs_datacube import GRID_CRS, GRID_TRANSFORM, GRID_WIDTH, GRID_HEIGHT  # canonical 3 m grid

OUT_DIR = "LIA_PRODUCTS"                                          # deliverable folder (git-ignored)
OUT_FP = os.path.join(OUT_DIR, "mcs_dtm_3m.tif")                  # the DEM the LIA scripts read

os.makedirs(OUT_DIR, exist_ok=True)                               # create the folder on first run
dtm = xr.open_zarr("ZARR/mores_creek_summit.zarr")["dtm"].sel(lidar_time="2023-02-09").squeeze()  # datacube terrain surface
profile = {"driver": "GTiff", "dtype": "float32", "count": 1,     # single-band float32 GeoTIFF
           "width": GRID_WIDTH, "height": GRID_HEIGHT,            # canonical grid shape
           "crs": GRID_CRS, "transform": GRID_TRANSFORM,          # EPSG:32611, 3 m pixels
           "nodata": np.nan, "compress": "deflate"}               # NaN = no lidar coverage
with rio.open(OUT_FP, "w", **profile) as dst:                     # write it
    dst.write(dtm.values.astype("float32"), 1)
print("wrote", OUT_FP)                                            # confirm
