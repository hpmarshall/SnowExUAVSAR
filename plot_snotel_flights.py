"""Plot Mores Creek Summit SNOTEL SWE with UAVSAR flight and lidar timing.

Two panels (winter 2019-20 and 2020-21): daily SNOTEL SWE (snotel_mcs_swe.csv),
vertical lines at every UAVSAR flight over MCS (line 23205 and 05208 in distinct
colors, dates from the datacube's flight coordinates), and the QSI lidar
acquisition dates. Output: snotel_swe_flights.png (committed).

Run:  PROJ_DATA=/opt/anaconda3/envs/myenv/share/proj \
      /opt/anaconda3/envs/myenv/bin/python plot_snotel_flights.py
"""

import matplotlib
matplotlib.use("Agg")          # headless
import matplotlib.pyplot as plt
import pandas as pd            # SNOTEL CSV
import xarray as xr           # flight dates from the datacube

ZARR_PATH = "ZARR/mores_creek_summit.zarr"
OUT = "snotel_swe_flights.png"
LIDAR_DATES = {"2020-02-09": "QSI lidar", "2021-03-15": "QSI lidar"}
LINE_COLOR = {"23205": "tab:red", "05208": "tab:orange"}         # per UAVSAR line
WINDOWS = [("2019-11-01", "2020-05-01", "Winter 2019–20"),       # panel date ranges
           ("2020-11-01", "2021-05-01", "Winter 2020–21")]

if __name__ == "__main__":
    sno = pd.read_csv("snotel_mcs_swe.csv", index_col="date", parse_dates=True)
    ds = xr.open_zarr(ZARR_PATH)
    flights = sorted(zip(ds["flight_line_id"].values, pd.to_datetime(ds["flight_date"].values)))

    fig, axes = plt.subplots(2, 1, figsize=(11, 7), sharey=True)
    for ax, (t0, t1, title) in zip(axes, WINDOWS):
        win = sno.loc[t0:t1]
        ax.plot(win.index, win["swe_mm"], color="tab:blue", lw=2, label="SNOTEL 637 SWE")
        seen = set()                                             # one legend entry per line id
        # both lines fly the same days in 2021 -- draw 23205 solid first, 05208 dotted on
        # top so coincident flights show as red-with-orange-dots instead of hiding one line
        for line, d in sorted(flights, key=lambda f: f[0], reverse=True):
            if not (pd.Timestamp(t0) <= d <= pd.Timestamp(t1)):
                continue
            lbl = f"UAVSAR {line}" if line not in seen else None
            seen.add(line)
            style = {"lw": 1.4, "ls": "-"} if line == "23205" else {"lw": 2.0, "ls": ":"}
            ax.axvline(d, color=LINE_COLOR[line], alpha=0.9, label=lbl, **style)
        for d, name in LIDAR_DATES.items():
            d = pd.Timestamp(d)
            if pd.Timestamp(t0) <= d <= pd.Timestamp(t1):
                ax.axvline(d, color="tab:green", lw=2.2, ls="--", label=name)
        ax.set_title(title)
        ax.set_ylabel("SWE (mm)")
        ax.grid(alpha=0.3)
        ax.legend(loc="upper left", fontsize=9)
    axes[1].set_xlabel("Date")
    fig.suptitle("Mores Creek Summit SNOTEL SWE with UAVSAR flights and QSI lidar", y=0.995)
    fig.tight_layout()
    fig.savefig(OUT, dpi=150)
    print(f"wrote {OUT}")
