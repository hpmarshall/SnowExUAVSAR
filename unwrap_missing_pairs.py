"""Unwrap the 5 InSAR pairs ASF ships no .unw.grd for, using SNAPHU.

Pairs 9, 12, 14, 20, 21 have wrapped phase + coherence in the datacube but no
unwrapped product -- including pair 20 (23205 2020-01-31 -> 02-13), the only pair
bracketing the 2020-02-09 QSI lidar. Per pair x pol:

  1. complex 5x5 block-mean of the 3 m wrapped phase -> 15 m interferogram
     (multilooking: coherence ~0.3 at 3 m is too noisy to unwrap directly);
  2. snaphu.unwrap(smooth cost, MCF init) with 5x5-averaged coherence;
  3. back to 3 m by snapping to the original wrapped phase:
     unw3 = angle(w3) + 2*pi*round((upsample(unw15) - angle(w3)) / 2*pi)
     -- full 3 m detail preserved, rewrap consistency exact by construction;
  4. QC: rewrap residual asserted ~0; connected components + coverage logged;
     one wrapped-vs-unwrapped PNG per pair (unwrap_qc_pair<NN>.png).

Results are region-written into 'unwrapped_phase' in both the Zarr and the GUI
NetCDF; provenance recorded as the 'custom_unwrapped_pairs' attr on the Zarr array.

Run:  PROJ_DATA=/opt/anaconda3/envs/myenv/share/proj \
      /opt/anaconda3/envs/myenv/bin/python unwrap_missing_pairs.py
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import netCDF4                 # GUI mirror patch
import numpy as np
import snaphu                  # SNAPHU 2.x via the snaphu-py wheel
import xarray as xr
import zarr                    # provenance attr + metadata re-consolidation

from compute_phase_sd_correlations import block_mean              # NaN-aware, complex-capable

ZARR_PATH = "ZARR/mores_creek_summit.zarr"
NC_PATH = "ZARR/mcs_gui.nc"
POLS = ["HH", "HV", "VH", "VV"]
MISSING = [9, 12, 14, 20, 21]  # pairs with no ASF .unw.grd (audit + phase_sd_correlations)
K = 5                          # 5x5 looks: 3 m -> 15 m before unwrapping


def upsample_nearest(a, k, shape):
    """Nearest-neighbor kxk upsample, edge-padded to `shape` (block_mean trimmed edges)."""
    up = np.kron(a, np.ones((k, k), dtype=a.dtype))               # each 15 m cell -> 5x5 block
    pad_y, pad_x = shape[0] - up.shape[0], shape[1] - up.shape[1]
    return np.pad(up, ((0, pad_y), (0, pad_x)), mode="edge")


if __name__ == "__main__":
    ds = xr.open_zarr(ZARR_PATH)
    ny, nx = ds.sizes["y"], ds.sizes["x"]

    for pi in MISSING:
        d1, d2 = str(ds["pair_date1"][pi].values)[:10], str(ds["pair_date2"][pi].values)[:10]
        line = str(ds["line_id"][pi].values)
        assert not np.isfinite(ds["unwrapped_phase"][pi, 0, ::8, ::8]).any(), \
            f"pair {pi} already has unwrapped phase -- refusing to overwrite an ASF product"
        print(f"[pair {pi}] {line} {d1} -> {d2}")

        pair_out = np.full((len(POLS), ny, nx), np.nan, dtype="float32")
        fig, axs = plt.subplots(2, len(POLS), figsize=(4 * len(POLS), 7))
        for pol_i, pol in enumerate(POLS):
            w3 = ds["wrapped_phase"][pi, pol_i].values            # complex, 3 m
            coh3 = ds["coherence"][pi, pol_i].values
            igram = block_mean(w3, K)                             # complex multilook -> 15 m
            corr = np.clip(block_mean(coh3, K), 0, 1)
            mask = np.isfinite(igram) & np.isfinite(corr)         # unwrap only real data
            unw15, cc = snaphu.unwrap(np.where(mask, np.nan_to_num(igram), 0),
                                      np.where(mask, np.nan_to_num(corr), 0),
                                      nlooks=K * K, cost="smooth", init="mcf", mask=mask)
            unw15 = np.where(mask & (cc > 0), unw15, np.nan)      # keep only unwrapped components
            ncomp = len(np.unique(cc[cc > 0]))
            # back to 3 m: nearest-upsample chooses the 2*pi branch, the wrapped phase
            # supplies the detail -- rewrap consistency is exact by construction
            up = upsample_nearest(unw15, K, (ny, nx))
            w_ang = np.angle(w3)
            unw3 = np.where(np.isfinite(w3.real) & np.isfinite(up),
                            w_ang + 2 * np.pi * np.round((up - w_ang) / (2 * np.pi)), np.nan)
            resid = np.angle(np.exp(1j * (unw3 - w_ang)))         # QC: must be ~0 where defined
            max_resid = float(np.nanmax(np.abs(resid))) if np.isfinite(resid).any() else 0.0
            assert max_resid < 1e-4, f"rewrap residual {max_resid} -- branch snapping failed"
            cov15, cov3 = float(np.isfinite(unw15).mean()), float(np.isfinite(unw3).mean())
            print(f"  {pol}: {ncomp} connected components, coverage 15m {cov15:.1%} / 3m {cov3:.1%}, "
                  f"rewrap residual {max_resid:.1e}")
            pair_out[pol_i] = unw3

            axs[0, pol_i].imshow(np.angle(igram), cmap="hsv", vmin=-np.pi, vmax=np.pi)
            axs[0, pol_i].set_title(f"{pol} wrapped (15 m)"); axs[0, pol_i].axis("off")
            im = axs[1, pol_i].imshow(unw15, cmap="RdBu_r")
            axs[1, pol_i].set_title(f"{pol} unwrapped, {ncomp} comp"); axs[1, pol_i].axis("off")
            fig.colorbar(im, ax=axs[1, pol_i], shrink=0.8)
        fig.suptitle(f"pair {pi}: {line} {d1} -> {d2} (SNAPHU smooth/MCF, 5x5 looks)")
        fig.tight_layout(); fig.savefig(f"unwrap_qc_pair{pi:02d}.png", dpi=120)
        plt.close(fig)

        # region-write into the Zarr and patch the GUI NetCDF, same pattern as reingest_pair.py
        patch = xr.Dataset({"unwrapped_phase": (("pair", "pol", "y", "x"), pair_out[None])})
        patch.chunk({"pair": 1, "pol": 1, "y": 512, "x": 512}).to_zarr(
            ZARR_PATH, mode="r+", region={"pair": slice(pi, pi + 1)})
        with netCDF4.Dataset(NC_PATH, "a") as nc:
            nc["unwrapped_phase"][pi] = pair_out
        print(f"  written to Zarr + NetCDF")

    # provenance: mark which pairs carry our unwrapping rather than the ASF product
    zg = zarr.open(ZARR_PATH, mode="r+")
    zg["unwrapped_phase"].attrs["custom_unwrapped_pairs"] = MISSING
    zg["unwrapped_phase"].attrs["custom_unwrap_method"] = \
        "SNAPHU smooth/MCF on 5x5-multilooked wrapped phase, 2pi-snapped back to 3 m"
    zarr.consolidate_metadata(ZARR_PATH)                          # xarray reads consolidated metadata
    print(f"done: pairs {MISSING} unwrapped and recorded in array attrs")
