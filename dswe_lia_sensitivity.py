"""How much does the calc_inc_angle LIA error change retrieved dSWE?

Uses the datacube's own inversion physics (compute_dswe.py / uavsar_pytools):
    phi = -(4 pi / lambda) * (cos t - sqrt(eps - sin^2 t)) * dz,   dz = dSWE/(rho/1000)
with lambda = 0.238403545 m (UAVSAR L-band), rho = 250 kg/m3 and Guneriussen-2001
permittivity eps = 1 + 1.6 rho' + 1.8 rho'^3 (rho' in g/cm3). The retrieved dSWE per
radian of phase is K(t) = -(lambda rho / 4 pi) / (cos t - sqrt(eps - sin^2 t)) [mm/rad],
so retrieving with the wrong LIA multiplies dSWE by K(t_wrong)/K(t_right).

Outputs: LIA_PRODUCTS/dswe_lia_sensitivity.png   (phase-vs-dSWE curves + K(t) panel)
         a markdown table on stdout: the facing/away LIA validation classes from
         compare_lia_versions.py, extended with the dSWE error they imply.

Run:  PROJ_DATA=/opt/anaconda3/envs/myenv/share/proj \
        /opt/anaconda3/envs/myenv/bin/python dswe_lia_sensitivity.py
"""

import numpy as np                                                # math
import matplotlib                                                 # headless backend
matplotlib.use("Agg")                                             # no display
import matplotlib.pyplot as plt                                   # figure
from matplotlib import cm                                         # sequential ramp for the theta family
from uavsar_pytools.snow_depth_inversion import depth_from_phase  # datacube inversion, for a round-trip check
from lia_pytools import DEM_FP, OUT_DIR, LINES, PIXEL_M, read_tif  # shared grid + reader

LAMBDA = 0.238403545                                              # UAVSAR L-band wavelength [m]
RHO = 250.0                                                       # datacube new-snow density [kg/m3]
EPS = 1 + 1.6 * (RHO / 1000) + 1.8 * (RHO / 1000) ** 3            # Guneriussen-2001 permittivity (1.428)
EXPECTED_AZ = {"05208": 142.0, "23205": 322.0}                    # ground->sensor azimuth per line


def k_mm_per_rad(theta_deg):
    """Retrieved dSWE per radian of phase, K(theta) [mm w.e. / rad]."""
    t = np.deg2rad(theta_deg)                                     # to radians
    return -LAMBDA * RHO / (4 * np.pi * (np.cos(t) - np.sqrt(EPS - np.sin(t) ** 2)))  # denominator < 0


for t in (20.0, 45.0, 70.0):                                      # round-trip check against the library's own
    lib = float(depth_from_phase(np.array([1.0]), np.deg2rad(t), density=RHO)) * RHO  # library: mm dSWE for 1 rad
    assert abs(lib - k_mm_per_rad(t)) < 1e-9, (t, lib)            # must match exactly

# ---- terrain classes: identical definitions to compare_lia_versions.py ----
dem = read_tif(DEM_FP)                                            # terrain
d_row, d_col = np.gradient(dem, PIXEL_M)                          # gradient along grid axes
zx, zy = d_col, -d_row                                            # dz/dEast, dz/dNorth (rows run south)
slope = np.degrees(np.arctan(np.hypot(zx, zy)))                   # slope [deg]
aspect = np.degrees(np.arctan2(-zx, -zy)) % 360                   # downslope azimuth [deg]

print("| line | pixels | corrected LIA | pytools LIA | dSWE ratio pyt/corr | dSWE error per 100 mm |")
print("|---|---|---|---|---|---|")
for line in LINES:                                                # 05208 then 23205
    cor = read_tif(f"{OUT_DIR}/lia_{line}_corrected.tif")         # right LIA
    pyt = read_tif(f"{OUT_DIR}/lia_{line}_pytools.tif")           # uavsar_pytools (wrong) LIA
    ok = np.isfinite(cor) & np.isfinite(pyt)                      # common valid pixels
    ratio = k_mm_per_rad(pyt) / k_mm_per_rad(cor)                 # per-pixel dSWE multiplier if pytools LIA is used
    band = ok & (slope > 20) & (slope < 25)                       # the validation slope class
    classes = [("facing 20-25° slopes", band, EXPECTED_AZ[line]),           # downslope aspect toward sensor
               ("away 20-25° slopes", band, (EXPECTED_AZ[line] + 180) % 360),  # opposite
               ("all valid pixels", ok, None)]                    # scene-wide net effect
    for name, base, center in classes:                            # one table row per class
        sel = base if center is None else base & (np.abs((aspect - center + 180) % 360 - 180) < 30)  # aspect gate
        r = np.nanmean(ratio[sel])                                # mean per-pixel dSWE multiplier
        print(f"| {line} {name} | {sel.sum()} | {np.nanmean(cor[sel]):.1f}° | {np.nanmean(pyt[sel]):.1f}° "
              f"| {r:.2f} | {100 * (r - 1):+.0f} mm ({100 * (r - 1):+.0f}%) |")

# ---- figure: forward curves + the sensitivity function ----
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.5, 5))           # two panels
dswe = np.linspace(0, 200, 201)                                   # dSWE range [mm w.e.]
thetas = np.arange(20, 81, 10)                                    # the incidence-angle family
ramp = cm.Blues(np.linspace(0.35, 0.95, len(thetas)))             # one-hue sequential: light=20deg, dark=80deg
dodge = {20: -5, 30: 3}                                           # nudge the two near-identical end labels apart
for th, c in zip(thetas, ramp):                                   # one curve per theta
    phi = dswe / k_mm_per_rad(th)                                 # forward model: phase for this dSWE
    ax1.plot(dswe, phi, color=c, lw=2)                            # the curve
    ax1.annotate(f"{th}°", (dswe[-1], phi[-1]), xytext=(4, dodge.get(th, 0)),  # direct label at the right end
                 textcoords="offset points", color=c, va="center", fontsize=10)
ax1.set(xlabel="ΔSWE (mm w.e.)", ylabel="interferometric phase (rad)", xlim=(0, 215),
        title=f"Phase vs ΔSWE by local incidence angle\n(ρ = {RHO:.0f} kg m⁻³, "
              f"Guneriussen 2001, λ = {100 * LAMBDA:.1f} cm)")
ax1.grid(alpha=0.25)                                              # recessive grid

tt = np.linspace(15, 88, 200)                                     # smooth theta axis
ax2.plot(tt, k_mm_per_rad(tt), color="#2b5f8a", lw=2)             # the sensitivity function
for line, marker in (("05208", "o"), ("23205", "s")):             # mark the validation classes on the curve
    cor = read_tif(f"{OUT_DIR}/lia_{line}_corrected.tif")         # reread (cheap) for class means
    pyt = read_tif(f"{OUT_DIR}/lia_{line}_pytools.tif")
    band = np.isfinite(cor) & np.isfinite(pyt) & (slope > 20) & (slope < 25)  # same slope class
    for name, center in (("facing", EXPECTED_AZ[line]), ("away", (EXPECTED_AZ[line] + 180) % 360)):
        sel = band & (np.abs((aspect - center + 180) % 360 - 180) < 30)       # aspect gate
        tc, tp = np.nanmean(cor[sel]), np.nanmean(pyt[sel])       # class-mean LIAs
        ax2.plot(tc, k_mm_per_rad(tc), marker, color="#2b5f8a", ms=8)          # where retrieval SHOULD sit
        ax2.plot(tp, k_mm_per_rad(tp), marker, color="#c0392b", ms=8, mfc="none")  # where pytools puts it
        ax2.annotate("", (tp, k_mm_per_rad(tp)), (tc, k_mm_per_rad(tc)),       # arrow: corrected -> pytools
                     arrowprops=dict(arrowstyle="->", color="#999999", lw=1))
        ax2.annotate(f"{line} {name}", ((tc + tp) / 2, k_mm_per_rad((tc + tp) / 2)),  # name the class pair
                     xytext=(8, 8), textcoords="offset points", fontsize=8, color="#555555")
ax2.plot([], [], "o", color="#2b5f8a", label="corrected LIA (filled)")         # legend proxies
ax2.plot([], [], "o", color="#c0392b", mfc="none", label="pytools LIA (open)")
ax2.legend(frameon=False, fontsize=9)                             # identify the two marker sets
ax2.set(xlabel="local incidence angle (°)", ylabel="retrieved ΔSWE per phase (mm w.e. / rad)",
        title="Retrieval sensitivity K(θ), with the 20-25°-slope\nfacing/away classes of both lines")
ax2.grid(alpha=0.25)                                              # recessive grid

fig.tight_layout()                                                # pack
fig.savefig(f"{OUT_DIR}/dswe_lia_sensitivity.png", dpi=150)       # deliverable figure
print(f"\nwrote {OUT_DIR}/dswe_lia_sensitivity.png")              # confirm
