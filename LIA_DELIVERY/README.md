# UAVSAR local incidence angle (LIA) — Mores Creek Summit, lines 05208 & 23205

Local incidence angle grids for the SnowEx UAVSAR lines lowman_05208 (heading 052°) and
lowman_23205 (heading 232°) over Mores Creek Summit, Idaho, computed two ways, plus the
scripts that produced them. All grids: 3 m pixels, UTM zone 11N (EPSG:32611), float32,
degrees, NaN = no data, on the same grid as the MCS UAVSAR/lidar datacube.

Contact: HP Marshall (hpmarshall@boisestate.edu), Boise State University.

## Contents

| file | what it is |
|---|---|
| `lia_pytools.py` | Script #1 — the **uavsar_pytools estimate**: `uavsar_pytools.incidence_angle.calc_inc_angle` applied to the JPL InSAR-Stack look vectors geocoded onto the lidar-DTM grid. Note: this library routine does not normalize the terrain surface normal and mis-pairs the horizontal look components with the DEM gradient, so it is biased on slopes (kept for provenance — it is what MCS products generated before 2026-10-30 used; the datacube has since been regenerated with the corrected version). |
| `lia_corrected.py` | Script #2 — the **corrected (textbook) LIA**: LIA = arccos(l̂ · n̂) with the unit ground→sensor look vector and the unit surface normal. Imports the geocoding helpers from `lia_pytools.py`, so keep the two files together. |
| `calc_lia_custom_dem.py` | Script #3 — compute LIA from the delivered look vectors with **your own DEM** (any CRS/resolution): `python calc_lia_custom_dem.py your_dem.tif 05208` (see *Look vectors* below) |
| `LIA_PRODUCTS/lia_<line>_pytools.tif` | LIA GeoTIFFs, uavsar_pytools estimate (2 files) |
| `LIA_PRODUCTS/lia_<line>_corrected.tif` | LIA GeoTIFFs, corrected version (2 files) |
| `LIA_PRODUCTS/mcs_dtm_3m.tif` | The terrain input: 2023-02-09 lidar DTM, block-averaged to 3 m (defines the grid; needed to rerun the scripts) |
| `LIA_PRODUCTS/lkv_<line>_{east,north,up}.tif` | **The radar look vectors** (6 files, ~120 MB) — a core deliverable, see *Look vectors* below. Also what lets the scripts rerun in seconds without the raw ~GB JPL Stack `.llh`/`.lkv` files. |
| `LIA_PRODUCTS/lia_comparison.png` | Maps of both versions + their difference, and difference vs terrain slope |
| `LIA_PRODUCTS/dswe_lia_sensitivity.png` | Phase vs ΔSWE curves by incidence angle, and the retrieval sensitivity K(θ) with the validation classes marked |

## Look vectors

`lkv_<line>_{east,north,up}.tif` hold the UAVSAR look vector per pixel so LIA can be
recomputed against **any DEM** (different lidar date, snow-on surface, coarser national
DEM, …):

- Components: **East / North / Up** (one GeoTIFF each), on the same 3 m EPSG:32611 grid
  as everything else, bilinearly geocoded from the JPL InSAR-Stack `.llh`/`.lkv` files.
- Direction: **sensor → ground** (the Up component is negative). Normalize before use —
  the magnitude is the slant range in metres (mean ≈ 16.5 km for 05208, ≈ 24.3 km for
  23205), itself occasionally useful.
- Validity exceeds the LIA maps (no DEM dependence): 7.35 M px (05208, full grid) and
  5.98 M px (23205).
- Geometry check (encoded in `compare_lia_versions.py`): the mean ground→sensor azimuth
  is 142° (line 05208, heading 052°) and 322° (line 23205, heading 232°) — UAVSAR is
  left-looking.

To compute LIA with your own DEM (any CRS, any resolution — it is resampled onto the
look-vector grid first):

    python calc_lia_custom_dem.py  path/to/your_dem.tif  05208   # or 23205
    # -> LIA_PRODUCTS/lia_05208_customdem.tif

The angle math is the corrected formula in `lia_corrected.py::corrected_lia()`; the
helper script only adds the DEM reprojection.

## Results

Validation classes are 20–25° slopes whose downslope aspect points toward ("facing") or
away from the radar (UAVSAR is left-looking; ground→sensor azimuth is 142° for 05208 and
322° for 23205). "Expected" is flat-ground incidence ∓ 22.5°. The corrected version
matches the expected geometry to ~1°; the uavsar_pytools estimate misses by 15–25°.
The last column is the resulting error in retrieved ΔSWE (Guneriussen-2001 inversion,
ρ = 250 kg/m³, λ = 23.84 cm) if the uavsar_pytools LIA is used instead of the corrected
one, per 100 mm of true ΔSWE:

| line | class | pixels | expected LIA | corrected LIA | uavsar_pytools LIA | ΔSWE error per 100 mm |
|---|---|---|---|---|---|---|
| 05208 | facing 20–25° slopes | 158,800 | ~28.8° | 28.7° | 44.3° | −15 mm (−15%) |
| 05208 | away 20–25° slopes | 91,425 | ~73.8° | 72.0° | 54.4° | +41 mm (+41%) |
| 05208 | all valid pixels | 3,799,117 | — | 49.8° | 48.6° | +6 mm (+6%) |
| 23205 | facing 20–25° slopes | 68,978 | ~42.7° | 43.3° | 59.1° | −21 mm (−21%) |
| 23205 | away 20–25° slopes | 145,956 | ~87.7° | 87.2° | 70.1° | +53 mm (+53%) |
| 23205 | all valid pixels | 3,347,603 | — | 70.7° | 64.1° | +24 mm (+24%) |

Recommendation: use the `*_corrected.tif` grids; the `*_pytools.tif` grids document what
the uavsar_pytools-based processing used.

## Rerunning

Requires Python with numpy, rasterio, GDAL and uavsar_pytools (e.g.
`conda install -c conda-forge numpy rasterio gdal; pip install uavsar_pytools`).
From this folder:

    python lia_pytools.py      # rewrites LIA_PRODUCTS/lia_<line>_pytools.tif
    python lia_corrected.py    # rewrites LIA_PRODUCTS/lia_<line>_corrected.tif

Both load the geocoded look vectors from `LIA_PRODUCTS/lkv_*.tif` if present; otherwise
they re-geocode from the JPL InSAR-Stack files, which must then sit in `UAVSAR_LKV/`
(`lowman_<line>_01_<BU|BC>_s<1-3>_2x8.{llh,lkv}` plus one `L090HH_01` `.ann` per line,
from the JPL UAVSAR InSAR Stack portal). If EPSG lookups fail, point `PROJ_DATA` at your
environment's `share/proj` directory.
