# SnowExUAVSAR
this repo summarizes the UAVSAR data collected for SnowEx, and starts by building a database for Mores Creek Summit area. See README.md for a complete file-by-file index.

## Tools to use
Use asf_search for accessing UAVSAR data
Use uavsar_pytools for reading the binary UAVSAR data

## Environment
All scripts here are run with the conda environment **`myenv`** -- it is the only
environment on this machine that has zarr, asf_search, uavsar_pytools and
phase_o_matic together. The repo's default `python` (anaconda3 base) does NOT have
zarr and cannot open the datacube.

    /opt/anaconda3/envs/myenv/bin/python build_mcs_datacube.py

(or `conda activate myenv` first). Recreate it elsewhere with
`conda env create -f environment.yml`; `environment-myenv-full.yml` is the exact
pinned freeze of the environment that produced the pilot datacube.

Credentials the scripts expect: Earthdata login in `~/.netrc` (UAVSAR downloads via
asf_search) and a CDS API key in `~/.cdsapirc` (ERA5 for phase_o_matic).

Known gotcha: pyproj/rasterio in myenv can pick up the anaconda base env's older
`proj.db`, breaking every EPSG lookup. Run geospatial scripts with
`PROJ_DATA=/opt/anaconda3/envs/myenv/share/proj` set.

## MCS pilot GUI (MATLAB)
`mcs_gui.m` is a four-panel viewer (Terrain/Atmosphere, LIA, UAVSAR, lidar) over the
pilot datacube. MATLAB here has no zarrread, so it reads `ZARR/mcs_gui.nc`, a
regenerable NetCDF mirror of the Zarr:

    PROJ_DATA=/opt/anaconda3/envs/myenv/share/proj \
      /opt/anaconda3/envs/myenv/bin/python export_mcs_gui_nc.py   # rebuild ZARR/mcs_gui.nc
    /Applications/MATLAB_R2025a.app/bin/matlab                    # then >> mcs_gui

`mcs_annotations.json` (Highway 21 + peak labels, committed) is rebuilt with
`fetch_mcs_annotations.py`. `mcs_gui_test.m` drives every control headlessly and
writes screenshots to `gui_test_png/`:

    matlab -sd <repo> -batch mcs_gui_test

## Area of interest
Mores Creek Summit, MCS_domain.kml

## First task
Start by querying the UAVSAR data on ASF using asf_search, and find all UAVSAR flights during the periods January 1, 2020-April 1, 2020 and January 1, 2021-April 1, 2021, in the Western U.S.

Make a table of all flights with location, state, date, and bearing.

## Plan Mode
Always use plan mode and create a plan for tasks, and ask for permission before executing.

## Code style
Keep code as concise as possible, and as simple to understand.  Comment every line. 

## References
Primary reference for SnowEx UAVSAR data: https://nsidc.org/sites/default/files/documents/technical-reference/snex_uavsar-v001-techref.pdf

Primary references for SnowEx LiDAR:
https://nsidc.org/data/snex20_qsi_sd_3m/versions/1
https://nsidc.org/data/snex23_lidar/versions/1
https://nsidc.org/sites/default/files/documents/user-guide/snex21_ps_dsm-v001-userguide.pdf
https://nsidc.org/sites/default/files/documents/user-guide/snex20_gm_swe_sd.pdf

## Extended search 
SnowEx 2017 (Grand Mesa, Telluride)  and flights in 2014-2016 in mountainous areas

## Task 2
Create kml of all lidar acquisitions associated with SnowEx.  These should come from the references for SnowEx LiDAR above, and the Airborne Snow Observatory.  We are interested in any airborne lidar flown during the months January-April that overlaps with the UAVSAR flights from 2014-2021.
