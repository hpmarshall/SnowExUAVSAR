"""Fetch Mores Creek Summit SNOTEL (637:ID:SNTL) daily SWE, snow depth and air temperature.

Pulls WTEQ (SWE, inches), SNWD (depth, inches) and TAVG (daily mean air
temperature, degF) from the NRCS AWDB REST API for
water years 2020 and 2021, converts to metric, and writes snotel_mcs_swe.csv
(committed -- it is small and the retrieval + plots depend on it).

Station metadata (from the same API): lat 43.93200, lon -115.66588, 6090 ft,
inside the MCS canonical grid.

Run:  /opt/anaconda3/envs/myenv/bin/python fetch_snotel_mcs.py
"""

import pandas as pd            # timeseries assembly + CSV
import requests                # AWDB REST API

AWDB = "https://wcc.sc.egov.usda.gov/awdbRestApi/services/v1/data"
TRIPLET = "637:ID:SNTL"                                          # Mores Creek Summit
BEGIN, END = "2019-10-01", "2021-06-30"                          # WY2020 + WY2021
OUT = "snotel_mcs_swe.csv"
IN2MM = 25.4                                                     # inches -> mm

if __name__ == "__main__":
    r = requests.get(AWDB, params={"stationTriplets": TRIPLET, "elements": "WTEQ,SNWD,TAVG",
                                    "beginDate": BEGIN, "endDate": END, "duration": "DAILY"},
                     headers={"User-Agent": "SnowExUAVSAR-research/1.0"}, timeout=60)
    r.raise_for_status()
    series = {}
    for el in r.json()[0]["data"]:                               # one entry per element
        code = el["stationElement"]["elementCode"]
        vals = {v["date"]: v["value"] for v in el["values"]}
        series[code] = pd.Series(vals, dtype=float)

    df = pd.DataFrame({"swe_mm": series["WTEQ"] * IN2MM,         # SWE inches -> mm
                       "depth_mm": series["SNWD"] * IN2MM,       # depth inches -> mm
                       "tavg_c": (series["TAVG"] - 32) * 5 / 9})  # degF -> degC
    df.index.name = "date"
    df = df.round(1).sort_index()
    df.to_csv(OUT)
    print(f"wrote {OUT}: {len(df)} days {df.index[0]} -> {df.index[-1]}, "
          f"peak SWE {df.swe_mm.max():.0f} mm")
