"""Find all SnowEx-associated UAVSAR flights (winter 2020 & 2021, Western U.S.) via asf_search.

Builds a table of flights (location, state, date, bearing) and cross-checks the dates
against Table 2 of the NSIDC SnowEx UAVSAR technical reference (snex_uavsar-v001-techref.pdf).
"""

import json                    # read/write the local query cache
import os                      # check whether the cache file exists
import re                      # regex for parsing UAVSAR scene names
import pandas as pd            # table handling and CSV output

# query windows: Jan-Apr 2014/2015/2016 (pre-SnowEx search), the 2017/2020/2021 SnowEx winters,
# and the 2024-2026 seasons (idahos campaign, March 2025)
WINDOWS = [("2014-01-01", "2014-05-01"), ("2015-01-01", "2015-05-01"), ("2016-01-01", "2016-05-01"),
           ("2016-11-01", "2017-05-01"), ("2019-11-01", "2020-05-01"), ("2020-11-01", "2021-05-01"),
           ("2024-11-01", "2026-05-01")]

# cache filename includes the windows so changing them automatically triggers a fresh query
CACHE = "asf_uavsar_cache_" + "_".join(f"{s}_{e}" for s, e in WINDOWS) + ".json"

# Western U.S. bounding box (lon -125..-102, lat 31..49) as WKT for the spatial filter
WESTERN_US = "POLYGON((-125 31,-102 31,-102 49,-125 49,-125 31))"

# JPL campaign code -> (SnowEx site name, state) from Table 2 of the NSIDC tech reference
SNOWEX_SITES = {
    "grmesa": ("Grand Mesa", "CO"),
    "irnton": ("Senator Beck Basin (Ironton)", "CO"),
    "tellur": ("Senator Beck Basin (Telluride)", "CO"),  # 2017 lines; omitted from the tech ref
    "Slumgl": ("Slumgullion Landslide (Lake City)", "CO"),  # not SnowEx, included at user request
    "GrMesa": ("Grand Mesa (pre-SnowEx 2015)", "CO"),       # single 2015-04-28 acquisition, 2 lines
    "tuolum": ("Tuolumne Meadows (Sierra Nevada)", "CA"),   # April 2016 Sierra snow flights, 3 lines
    "SMAP13": ("SMAP cal/val (Sierra Nevada crossing)", "CA"),  # 2014/15; line 13500 crosses Tuolumne
    "idahos": ("Idaho 2025 (Boise Mtns/Sawtooths/McCall)", "ID"),  # Mar 2025; line 27441 covers Mores Creek
    "peeler": ("East River (Peeler Peak)", "CO"),
    "fraser": ("Fraser Experimental Forest", "CO"),
    "rockmt": ("Cameron Pass (Rocky Mtn NP)", "CO"),
    "lowman": ("Boise River Basin (Lowman)", "ID"),
    "silver": ("Reynolds Creek (Silver City)", "ID"),
    "stlake": ("Little Cottonwood Canyon (Salt Lake City)", "UT"),
    "uticam": ("Central Agricultural Research Center (Utica)", "MT"),
    "alamos": ("Jemez River (Los Alamos)", "NM"),
    "dorado": ("American River Basin (Eldorado NF)", "CA"),
    "donner": ("Sagehen Creek (Donner)", "CA"),
    "sierra": ("Lakes Basin (Sierra NF)", "CA"),
}

# acquisition dates per site listed in Table 2 of the tech reference (within the query windows)
TECHREF_DATES = {
    "grmesa": ["2017-02-06", "2017-02-22", "2017-02-25", "2017-03-08", "2017-03-31",  # 2017 dates from tech-ref text (sec. 4)
               "2020-02-01", "2020-02-12", "2020-02-19", "2020-02-26", "2020-03-12",
               "2021-01-27", "2021-02-03", "2021-02-10", "2021-03-03", "2021-03-10", "2021-03-16", "2021-03-22"],
    "irnton": ["2020-02-12", "2020-02-19", "2020-02-26", "2020-03-12",
               "2021-01-15", "2021-01-21", "2021-01-28", "2021-02-04", "2021-02-11",
               "2021-02-23", "2021-03-03", "2021-03-10", "2021-03-16", "2021-03-22"],
    "peeler": ["2019-12-20", "2020-02-12", "2020-02-19", "2020-02-26", "2020-03-12"],
    "fraser": ["2020-02-12", "2020-02-19", "2020-02-26", "2020-03-12",
               "2021-01-15", "2021-01-20", "2021-01-27", "2021-02-03", "2021-02-23",
               "2021-03-03", "2021-03-10", "2021-03-16", "2021-03-22"],
    "rockmt": ["2020-02-12", "2020-02-19", "2020-02-26", "2020-03-12",
               "2021-01-15", "2021-01-20", "2021-01-27", "2021-02-03", "2021-02-23",
               "2021-03-03", "2021-03-10", "2021-03-16", "2021-03-22"],
    "lowman": ["2019-12-20", "2020-01-31", "2020-02-13", "2020-02-21", "2020-03-11",
               "2021-01-15", "2021-01-20", "2021-01-27", "2021-02-03", "2021-02-10",
               "2021-02-23", "2021-03-03", "2021-03-10", "2021-03-16", "2021-03-22"],
    "silver": ["2020-01-31", "2020-02-13", "2020-02-21", "2020-03-11"],
    "stlake": ["2020-01-31", "2020-02-13", "2020-02-21", "2020-03-12",
               "2021-01-15", "2021-01-21", "2021-01-28", "2021-02-03", "2021-02-10",
               "2021-02-23", "2021-03-03", "2021-03-10", "2021-03-16", "2021-03-22"],
    "uticam": ["2021-01-15", "2021-01-20", "2021-02-23"],
    "alamos": ["2020-02-12", "2020-02-19", "2020-02-26"],
    "dorado": ["2020-01-31", "2020-02-12", "2020-02-19", "2020-02-26", "2020-03-11"],
    "donner": ["2019-12-20", "2020-01-31", "2020-02-12", "2020-02-19", "2020-02-26", "2020-03-11"],
    "sierra": ["2020-01-31", "2020-02-12", "2020-02-19", "2020-02-26", "2020-03-12"],
}

# regex for single-acquisition (PolSAR/CX) scene names: UA_<campaign>_<HHHnn>_<YYnnn>_<nnn>_<YYMMDD>_...
# (line counter <nn> may contain letters, e.g. tellur_1701T / tellur_3503F)
CX_NAME = re.compile(r"^UA_(\w+?)_(\d{3})(\w{2})_(\d{5})_(\d{3})_(\d{6})_")

def load_flights():
    """Return parsed UAVSAR flights as a list of dicts, querying ASF only if no local cache exists."""
    if os.path.exists(CACHE):                                   # cache hit: skip the remote query
        return json.load(open(CACHE))                           # load previously parsed results
    import asf_search as asf                                    # ASF client (imported only when querying)
    rows = []                                                   # accumulator for parsed flights
    for start, end in WINDOWS:                                  # query each SnowEx winter window
        results = asf.geo_search(platform=asf.PLATFORM.UAVSAR,  # UAVSAR platform only
                                 intersectsWith=WESTERN_US,     # Western U.S. spatial filter
                                 processingLevel="COMPLEX",     # one CX product per flight/line
                                 start=start, end=end)          # temporal filter
        print(f"{start}..{end}: {len(results)} COMPLEX products")  # progress report
        for p in results:                                       # loop over returned products
            m = CX_NAME.match(p.properties["sceneName"])        # parse the scene name
            if not m:                                           # warn on any unparseable name (never drop silently)
                print("WARNING: could not parse", p.properties["sceneName"])
                continue
            camp, hdg, cnt, flt, acq, ymd = m.groups()          # campaign, heading, counter, flight, acq#, date
            date = f"20{ymd[:2]}-{ymd[2:4]}-{ymd[4:6]}"         # YYMMDD -> ISO date
            rows.append({"campaign": camp, "date": date, "bearing_deg": int(hdg),  # parsed fields
                         "line_id": f"{hdg}{cnt}", "scene": p.properties["sceneName"],
                         "geometry": p.geometry})               # GeoJSON footprint for mapping
    json.dump(rows, open(CACHE, "w"))                           # save cache for future runs
    return rows                                                 # hand parsed flights to the caller


if __name__ == "__main__":                                      # run the table build only as a script
    df = pd.DataFrame(load_flights()).drop_duplicates(["campaign", "date", "line_id"])  # unique flights only

    snowex = df[df.campaign.isin(SNOWEX_SITES)].copy()          # keep SnowEx campaigns for the main table
    snowex["location"] = snowex.campaign.map(lambda c: SNOWEX_SITES[c][0])  # add SnowEx site name
    snowex["state"] = snowex.campaign.map(lambda c: SNOWEX_SITES[c][1])     # add state
    snowex = snowex.sort_values(["state", "location", "date", "bearing_deg"])  # tidy ordering
    table = snowex[["location", "state", "date", "bearing_deg", "campaign", "line_id"]]  # final columns
    table.to_csv("snowex_uavsar_flights.csv", index=False)      # save the flight table
    print(f"\n{len(table)} SnowEx flights saved to snowex_uavsar_flights.csv")  # summary

    other = df[~df.campaign.isin(SNOWEX_SITES)]                 # non-SnowEx flights, for completeness
    print(f"{len(other)} non-SnowEx UAVSAR flights in window: "  # report their campaign codes
          f"{sorted(other.campaign.unique())}")

    print("\n=== Cross-check vs NSIDC tech reference (Table 2) ===")  # header for the comparison
    for camp, ref_dates in TECHREF_DATES.items():               # compare each site
        asf_dates = set(snowex[snowex.campaign == camp].date)   # dates found on ASF
        ref = set(ref_dates)                                    # dates listed in the tech reference
        missing = sorted(ref - asf_dates)                       # in tech ref but not on ASF
        extra = sorted(asf_dates - ref)                         # on ASF but not in tech ref
        status = "OK" if not missing and not extra else f"missing={missing} extra={extra}"  # verdict
        print(f"{camp:8s} ({SNOWEX_SITES[camp][0][:30]:30s}): "  # per-site result line
              f"{len(asf_dates):2d} ASF / {len(ref):2d} ref -> {status}")
