# DRAFT GitHub issue for SnowEx/uavsar_pytools — review, then paste title + body
# (this file is a local draft, not part of the repo deliverables; delete after posting)

**Title:** `calc_inc_angle`: surface normal never normalized and DEM gradient mis-paired with look-vector components — large incidence-angle errors on slopes

---

## Summary

`uavsar_pytools.incidence_angle.calc_inc_angle` returns incorrect local incidence
angles wherever the terrain is not flat. Two independent defects in the geometry
compound: errors reach 15–25° on 20–25° slopes with real UAVSAR data, and in a clean
synthetic case the function either ignores the terrain entirely or inverts the
facing/away response, depending on argument order. Anything derived from the LIA
downstream (e.g. SWE retrieval from repeat-pass phase) inherits the error with the
wrong sign pattern across aspect. A 30-second, data-free reproduction is below, along
with a suggested fix that keeps the existing function signature.

Found while validating the SnowEx Mores Creek Summit (lowman 05208/23205) InSAR
products at Boise State; happy to open a PR with the fix plus a regression test.

## The two defects

Both are in the final dot product
([`incidence_angle.py` L75](https://github.com/SnowEx/uavsar_pytools/blob/main/uavsar_pytools/incidence_angle.py)):

```python
inc_cos = unit_lkv['x']*dx + unit_lkv['y']*dy + unit_lkv['z']
```

**1. The surface normal is never normalized.** This is the dot product of a *unit*
look vector with the *unnormalized* terrain normal `(dx, dy, 1)`, whose magnitude is
`sqrt(1 + dx² + dy²) > 1` on any slope. The result is therefore not a cosine — it is
too large by exactly that factor, and the error grows with slope (~1° of LIA error
per degree of slope near mid incidence). The helper `arccos_theta` then handles the
out-of-range values with an ad-hoc fold (`arccos(2 - v)` for `v > 1`), which masks
the symptom instead of fixing it and introduces a discontinuous, non-physical mapping.
The correct expression divides by the normal's magnitude, after which no folding is
needed.

**2. The horizontal pairing/sign between look vector and DEM gradient is wrong.**
`np.gradient(dem, pixel_size)` returns derivatives along *(rows, cols)*. On a
standard north-up raster, rows run **south** and columns run **east**, so
`dx, dy = np.gradient(...)` is actually `(-dz/dNorth, +dz/dEast)`. The correct
cosine, with ground→sensor unit vector `(l_e, l_n, l_u)` and unit normal
`(-zx, -zy, 1)/sqrt(1 + zx² + zy²)` (where `zx = dz/dEast`, `zy = dz/dNorth`), is

```
cos(LIA) = (-zx*l_e - zy*l_n + l_u) / sqrt(1 + zx² + zy²)
```

No assignment of `(lkv_x, lkv_y)` makes L75 match this:

- pass `(east, north)` → each horizontal term pairs a look component with the
  *orthogonal* gradient component (`l_e` with the north–south derivative, `l_n` with
  the east–west one);
- pass `(north, east)` → the north term comes out right, but the east term has the
  **wrong sign**.

The docstring specifies neither what `lkv_x`/`lkv_y` mean nor the look-vector sign
convention (the code's negation at L72 implies the JPL sensor→ground convention),
so users can't even work around it consistently.

## Minimal reproduction (no data files)

A uniform plane tilted 20°, radar at 45° flat-ground incidence due east of the
target, so the exact answer is known analytically: **25°** when the plane dips toward
the sensor, **65°** when it dips away.

```python
import numpy as np
from uavsar_pytools.incidence_angle import calc_inc_angle

PIX, N = 5.0, 200                                   # 5 m pixels, 200x200 grid
THETA, SLOPE = np.radians(45), np.radians(20)       # flat incidence, plane slope

cols = np.arange(N) * PIX                           # easting of each column
le, ln, lu = np.sin(THETA), 0.0, np.cos(THETA)      # ground->sensor unit vector (E,N,U)
rng = 20000.0                                       # arbitrary slant range [m]
lkv_e = np.full((N, N), -le * rng)                  # sensor->ground, as shipped in .lkv
lkv_n = np.full((N, N), -ln * rng)
lkv_u = np.full((N, N), -lu * rng)

def textbook_lia(dem):
    """LIA = arccos(l_hat . n_hat), unit surface normal, north-up grid."""
    d_row, d_col = np.gradient(dem, PIX)            # rows run SOUTH, cols run EAST
    zx, zy = d_col, -d_row                          # dz/dEast, dz/dNorth
    n_norm = np.sqrt(1 + zx**2 + zy**2)             # |(-zx, -zy, 1)|
    cos_lia = (-zx*le - zy*ln + lu) / n_norm
    return np.degrees(np.arccos(np.clip(cos_lia, -1, 1)))

for label, sign, exact in [("facing (dips toward sensor)", -1, 25.0),
                           ("away   (dips away)         ", +1, 65.0)]:
    dem = sign * np.tan(SLOPE) * np.tile(cols, (N, 1))          # tilted plane
    inc_xy = calc_inc_angle(dem, lkv_e, lkv_n, lkv_u, PIX)      # lkv_x=east, lkv_y=north
    inc_yx = calc_inc_angle(dem, lkv_n, lkv_e, lkv_u, PIX)      # lkv_x=north, lkv_y=east
    fix = textbook_lia(dem)
    c = slice(5, -5), slice(5, -5)                              # ignore gradient edge rows
    print(f"{label}  exact {exact:5.1f} | calc_inc_angle(E,N) {np.nanmean(inc_xy[c]):6.2f}"
          f" | calc_inc_angle(N,E) {np.nanmean(inc_yx[c]):6.2f} | textbook {np.nanmean(fix[c]):6.2f}")
```

Output (uavsar_pytools 0.7.1):

```
facing (dips toward sensor)  exact  25.0 | calc_inc_angle(E,N)  45.00 | calc_inc_angle(N,E)  63.27 | textbook  25.00
away   (dips away)           exact  65.0 | calc_inc_angle(E,N)  45.00 | calc_inc_angle(N,E)  15.32 | textbook  65.00
```

With `(east, north)` ordering the terrain has **no effect at all** for this geometry
(45.00° on a 20° slope); with `(north, east)` the facing/away response is
**inverted** (a slope facing the radar reads 63° instead of 25°). The textbook
formula is exact.

## Magnitude on real data

On the SnowEx Mores Creek Summit lines (lowman 05208 & 23205, JPL InSAR-Stack look
vectors geocoded onto a 3 m lidar DTM), comparing `calc_inc_angle` against the
textbook formula, with an independent geometric expectation (flat-ground incidence
∓ slope) as referee:

| line | 20–25° slope class | expected LIA | textbook LIA | calc_inc_angle |
|---|---|---|---|---|
| 05208 | facing radar | ~28.8° | 28.7° | 44.3° |
| 05208 | facing away | ~73.8° | 72.0° | 54.4° |
| 23205 | facing radar | ~42.7° | 43.3° | 59.1° |
| 23205 | facing away | ~87.7° | 87.2° | 70.1° |

Propagated through a standard dry-snow SWE inversion (Guneriussen & Høgda 2001),
the LIA error biases retrieved ΔSWE by −15% to −21% on slopes facing the radar and
+41% to +53% on slopes facing away — i.e. it imprints a spurious aspect signal on
any slope-dependent retrieval.

## Suggested fix

Keep the signature, define the conventions explicitly in the docstring
(`lkv_x` = East, `lkv_y` = North, `lkv_z` = Up, sensor→ground as shipped in JPL
`.lkv` files, DEM on a north-up grid), and replace the geometry:

```python
d_row, d_col = np.gradient(dem_arr, pixel_size)   # rows run south on a north-up grid
zx, zy = d_col, -d_row                            # dz/dEast, dz/dNorth
mag = np.sqrt(lkv_e**2 + lkv_n**2 + lkv_u**2)
mag[mag == 0] = np.nan
le, ln, lu = -lkv_e/mag, -lkv_n/mag, -lkv_u/mag   # ground->sensor unit vector
n_norm = np.sqrt(1 + zx**2 + zy**2)               # |(-zx, -zy, 1)|
inc = np.degrees(np.arccos(np.clip((-zx*le - zy*ln + lu)/n_norm, -1, 1)))
```

`arccos_theta` becomes unnecessary (the clipped argument is a true cosine). The
tilted-plane case above makes a good regression test since its answer is analytic.
Note this will change results for every user's sloped pixels, so it probably
warrants a version bump and a changelog warning.

## Environment

uavsar_pytools 0.7.1, numpy 2.2.5, Python 3.10.14, macOS.
