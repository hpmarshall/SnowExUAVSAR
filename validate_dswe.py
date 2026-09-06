"""Validate the dSWE retrieval four ways -> dswe_validation.md (+ PNG figures).

1. Pair-chain closure: triangles (d1->d2)+(d2->d3)-(d1->d3) on atmosphere-corrected
   PHASE (per-pair reference constants cancel in the triangle; the SNOTEL anchor would
   not), reported in radians and dSWE-equivalent.
2. Cross-polarization consistency: pairwise correlation + RMSD between pols' dSWE at 90 m.
3. SNOTEL anchor diagnostics: the per-pair anchor offset (from dswe_retrieval_log.csv) --
   how far the scene sat from the station before anchoring -- plus interval density.
4. Lidar: 90 m corr of QSI snow depth vs dSWE -- 2021-03-15 vs bracketing pairs 3/4
   (clean test), 2020-02-09 vs adjacent-window pairs 19/18 (proxy only; the bracketing
   2020 pair (20) has no unwrapped phase).

Run:  PROJ_DATA=/opt/anaconda3/envs/myenv/share/proj \
      /opt/anaconda3/envs/myenv/bin/python validate_dswe.py
"""

import json                    # manifest for triangle search
import itertools               # pol combinations
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr

from compute_phase_sd_correlations import block_mean, pearson     # shared 90 m helpers
from compute_dswe import DENSITY                                  # same density assumption
from uavsar_pytools.snow_depth_inversion import depth_from_phase  # rad -> dSWE-equivalent scale

ZARR_PATH = "ZARR/mores_creek_summit.zarr"
POLS = ["HH", "HV", "VH", "VV"]
OUT_MD = "dswe_validation.md"

if __name__ == "__main__":
    ds = xr.open_zarr(ZARR_PATH)
    manifest = json.load(open("uavsar_pair_manifest_mcs.json"))
    key = {(r["line_id"], r["date1"], r["date2"]): i for i, r in enumerate(manifest)}
    md = ["# dSWE retrieval validation", "",
          f"Retrieval: Guneriussen-2001, density {DENSITY:.0f} kg/m3, atmosphere-corrected,",
          "coherence >= 0.35 mask, anchored at MCS SNOTEL (see compute_dswe.py).", ""]

    # ---- 1. closure triangles on corrected phase ----
    triangles = [(i, j, k) for (l1, a, b), i in key.items() for (l2, b2, c), j in key.items()
                 if l1 == l2 and b2 == b for (l3, a3, c3), k in key.items()
                 if l3 == l1 and a3 == a and c3 == c]
    md += ["## 1. Pair-chain closure (corrected phase; constants cancel)", "",
           "| triangle | line | pol | median (rad) | RMS (rad) | RMS dSWE-equiv (mm) | n (90m) |",
           "|---|---|---|---|---|---|---|"]
    print(f"{len(triangles)} closure triangles found")
    mean_lia = {ln: np.deg2rad(float(np.nanmedian(ds["lia"].sel(line=ln).values)))
                for ln in ["05208", "23205"]}                     # dask nanmedian needs loaded values
    for (i, j, k) in triangles:
        line = manifest[i]["line_id"]
        tri_lbl = f"{manifest[i]['date1']}->{manifest[i]['date2']}->{manifest[j]['date2']}"
        for pol_i, pol in enumerate(POLS):
            phi = [ds["unwrapped_phase"][p, pol_i].values - ds["atm_delay_diff"][p].values
                   for p in (i, j, k)]
            closure = block_mean(phi[0]) + block_mean(phi[1]) - block_mean(phi[2])  # 90 m
            v = closure[np.isfinite(closure)]
            if v.size < 10:
                continue
            v = v - np.median(v)                                  # per-pair constants leave one net offset
            rms = float(np.sqrt(np.mean(v**2)))
            # dSWE-equivalent of the RMS closure at the line's median LIA
            dswe_eq = abs(float(depth_from_phase(np.array([rms]), mean_lia[line],
                                                  density=DENSITY)[0])) * DENSITY  # m snow -> mm w.e.
            md.append(f"| {tri_lbl} | {line} | {pol} | {np.median(v):+.3f} | {rms:.3f} "
                      f"| {dswe_eq:.1f} | {v.size} |")
    md.append("")

    # ---- 2. cross-pol consistency at 90 m ----
    md += ["## 2. Cross-polarization consistency (dSWE, 90 m)", "",
           "| pair | line | dates | mean corr | min corr | mean RMSD (mm) |", "|---|---|---|---|---|---|"]
    for pi in range(ds.sizes["pair"]):
        maps = [block_mean(ds["dswe"][pi, p].values) for p in range(4)]
        if not any(np.isfinite(m).any() for m in maps):
            continue
        cors, rmsds = [], []
        for a, b in itertools.combinations(range(4), 2):
            r, n = pearson(maps[a], maps[b])
            m = np.isfinite(maps[a]) & np.isfinite(maps[b])
            if n >= 10:
                cors.append(r)
                rmsds.append(1000 * float(np.sqrt(np.mean((maps[a][m] - maps[b][m])**2))))
        if not cors:
            continue
        d1, d2 = str(ds["pair_date1"][pi].values)[:10], str(ds["pair_date2"][pi].values)[:10]
        md.append(f"| {pi} | {str(ds['line_id'][pi].values)} | {d1}->{d2} "
                  f"| {np.mean(cors):+.3f} | {np.min(cors):+.3f} | {np.mean(rmsds):.1f} |")
    md.append("")

    # ---- 3. SNOTEL anchor diagnostics ----
    log = pd.read_csv("dswe_retrieval_log.csv")
    hh = log[log.pol == "HH"].copy()
    md += ["## 3. SNOTEL anchor diagnostics (HH)", "",
           "SNOTEL is the retrieval's reference, so this is a consistency diagnostic, not",
           "independent validation: anchor_offset is how far the pre-anchor scene sat from",
           "the station value (large |offset| = large unknown phase constant / ramp).", "",
           "| pair | dates | SNOTEL dSWE (mm) | anchor offset (mm) | anchor window (m) | interval density (kg/m3) |",
           "|---|---|---|---|---|---|"]
    for _, r in hh.iterrows():
        if r.window_n == 0 and r.finite_frac == 0:
            continue
        rho = f"{r.snotel_density:.0f}" if np.isfinite(r.snotel_density) and r.snotel_density < 1000 else ""
        off = f"{1000*r.anchor_offset_m:+.1f}" if np.isfinite(r.anchor_offset_m) else "unanchored"
        w = f"{r.anchor_window_m:.0f}" if np.isfinite(r.anchor_offset_m) else "-"
        md.append(f"| {r.pair} | {r.date1}->{r.date2} | {1000*r.snotel_dswe_m:+.1f} | {off} | {w} | {rho} |")
    md.append("")

    # ---- 4. lidar comparisons ----
    md += ["## 4. Lidar snow depth vs dSWE (90 m Pearson r)", "",
           "| pair | dates | lidar | type | " + " | ".join(POLS) + " |",
           "|---|---|---|---|---|---|---|---|"]
    cases = [(3, "2021-03-15", "bracketing"), (4, "2021-03-15", "bracketing"),
             (19, "2020-02-09", "proxy (adjacent window)"), (18, "2020-02-09", "proxy (adjacent window)")]
    fig, axes = plt.subplots(1, len(cases), figsize=(4.2 * len(cases), 4))
    for ax, (pi, ld, typ) in zip(axes, cases):
        sd90 = block_mean(ds["sd"].sel(sd_time=ld).squeeze().values)
        rs = []
        for pol_i in range(4):
            r, n = pearson(sd90, block_mean(ds["dswe"][pi, pol_i].values))
            rs.append(f"{r:+.3f}" if np.isfinite(r) else "")
        d1, d2 = str(ds["pair_date1"][pi].values)[:10], str(ds["pair_date2"][pi].values)[:10]
        md.append(f"| {pi} | {d1}->{d2} | {ld} | {typ} | " + " | ".join(rs) + " |")
        dw = block_mean(ds["dswe"][pi, 0].values)                 # HH scatter figure
        m = np.isfinite(sd90) & np.isfinite(dw)
        ax.plot(sd90[m], 1000 * dw[m], ".", ms=2, alpha=0.3)
        ax.set_xlabel(f"lidar SD {ld} (m)"); ax.set_ylabel("dSWE HH (mm w.e.)")
        ax.set_title(f"pair {pi} {d1}->{d2}\n{typ}, r={rs[0]}", fontsize=9)
        ax.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig("dswe_vs_lidar.png", dpi=130)
    md += ["", "Scatter (HH): dswe_vs_lidar.png. The pair bracketing the 2020 lidar (20,",
           "2020-01-31->02-13) has no unwrapped phase from ASF, so no dSWE map exists for it;",
           "pairs 19/18 compare a fixed snapshot against a different window (proxy only).", ""]

    open(OUT_MD, "w").write("\n".join(md) + "\n")
    print(f"wrote {OUT_MD} and dswe_vs_lidar.png")
