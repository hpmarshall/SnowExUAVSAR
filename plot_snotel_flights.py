"""Plot Mores Creek Summit SNOTEL SWE + air temperature with UAVSAR/lidar timing
and the selected InSAR pairs.

Two panels (winter 2019-20 and 2020-21): daily SNOTEL SWE (left axis) and daily mean
air temperature in Celsius (right axis, with a 0 degC line), vertical lines at every
UAVSAR flight (23205 solid red, 05208 dotted orange) and QSI lidar date (dashed green),
and -- bottom panel, where all six selected pairs fall -- double-headed arrows spanning
each selected pair's interval, labeled "Pair A, r=+0.64" etc. Letters come from
rank_dswe_pairs.PAIR_LETTERS and the r values (polarization-mean, season-matched lidar)
are read live from dswe_lidar_rank.csv, so the figure always matches the report tables.

Output: snotel_swe_flights.png (committed).

Run:  PROJ_DATA=/opt/anaconda3/envs/myenv/share/proj \
      /opt/anaconda3/envs/myenv/bin/python plot_snotel_flights.py
"""

import matplotlib
matplotlib.use("Agg")          # headless
import matplotlib.pyplot as plt
import pandas as pd            # SNOTEL CSV
import xarray as xr           # flight/pair dates from the datacube

from rank_dswe_pairs import PAIR_LETTERS, ZERO_CHANGE_EXAMPLE     # shared letter labels

ZARR_PATH = "ZARR/mores_creek_summit.zarr"
OUT = "snotel_swe_flights.png"
LIDAR_DATES = {"2020-02-09": "QSI lidar", "2021-03-15": "QSI lidar"}
LINE_COLOR = {"23205": "tab:red", "05208": "tab:orange"}         # per UAVSAR line
WINDOWS = [("2019-11-01", "2020-05-01", "Winter 2019–20"),       # panel date ranges
           ("2020-11-01", "2021-05-01", "Winter 2020–21")]
ARROW_ROW_MM = 48              # vertical spacing between stacked pair arrows (mm SWE units)




if __name__ == "__main__":
    sno = pd.read_csv("snotel_mcs_swe.csv", index_col="date", parse_dates=True)
    ds = xr.open_zarr(ZARR_PATH)
    flights = sorted(zip(ds["flight_line_id"].values, pd.to_datetime(ds["flight_date"].values)))

    # selected pairs: letter, interval, and season-matched pol-mean r from the ranking CSV
    rank = pd.read_csv("dswe_lidar_rank.csv")
    sel = []
    for pi, letter in PAIR_LETTERS.items():
        sub = rank[rank.pair == pi]
        matched = "r_lidar2021" if sub.season.iloc[0] == 2021 else "r_lidar2020"
        sel.append(dict(letter=letter, pair=pi,
                        t0=pd.Timestamp(sub.date1.iloc[0]), t1=pd.Timestamp(sub.date2.iloc[0]),
                        r=sub[matched].mean(), zero=pi == ZERO_CHANGE_EXAMPLE))
    sel.sort(key=lambda s: s["letter"])

    fig, axes = plt.subplots(2, 1, figsize=(11, 8.5))
    for ax, (t0, t1, title) in zip(axes, WINDOWS):
        win = sno.loc[t0:t1]
        ax.plot(win.index, win["swe_mm"], color="tab:blue", lw=2, label="SNOTEL 637 SWE", zorder=3)
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

        # right axis: daily mean air temperature (degC), thin so SWE stays dominant
        axr = ax.twinx()
        axr.plot(win.index, win["tavg_c"], color="tab:purple", lw=0.9, alpha=0.6,
                 label="air temperature")
        axr.axhline(0, color="gray", lw=1, ls="--", alpha=0.8)   # the 0 degC line
        axr.set_ylabel("Air temperature (°C)", color="tab:purple")
        axr.tick_params(axis="y", labelcolor="tab:purple")
        h1, l1 = ax.get_legend_handles_labels()
        h2, l2 = axr.get_legend_handles_labels()
        ax.legend(h1 + h2, l1 + l2, loc="upper left", fontsize=9)

        # selected-pair arrows: all six fall in winter 2020-21 (bottom panel).
        # One row per pair (letter order) so labels never collide horizontally.
        in_panel = [s for s in sel if pd.Timestamp(t0) <= s["t0"] and s["t1"] <= pd.Timestamp(t1)]
        band = ARROW_ROW_MM * (len(in_panel) + 0.8) if in_panel else 20
        ax.set_ylim(-band, 780)                                  # reserve a band below SWE=0
        # keep the temperature curve above the arrow band: map its display range onto
        # the panel fraction the SWE data area occupies
        tmin, tmax = -14.0, 46.0                                 # temp drawn in lower ~half of data area
        frac = band / (band + 780.0)
        axr.set_ylim(tmin - (tmax - tmin) * frac / (1 - frac), tmax)
        for ri, s in enumerate(in_panel):
            y = -ARROW_ROW_MM * (ri + 1)
            color = "dimgray" if s["zero"] else "black"
            ax.annotate("", xy=(s["t0"], y), xytext=(s["t1"], y),
                        arrowprops=dict(arrowstyle="<->", color=color, lw=1.4))
            ax.text(s["t1"] + pd.Timedelta(days=2), y,
                    f"Pair {s['letter']}, r={s['r']:+.2f}", ha="left", va="center",
                    fontsize=8.5, color=color, fontweight="bold")
    axes[1].set_xlabel("Date")
    fig.suptitle("Mores Creek Summit SNOTEL SWE and air temperature, with UAVSAR flights, "
                 "QSI lidar, and the selected InSAR pairs", y=0.995)
    fig.tight_layout()
    fig.savefig(OUT, dpi=150)
    print(f"wrote {OUT}")
