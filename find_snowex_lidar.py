"""Build a KML + table of SnowEx-associated airborne lidar acquisitions (Jan-April)
that overlap the 2014-2021 UAVSAR flight lines, per AGENTS.md Task 2.

Sources: NSIDC/CMR granules (ASO 2013-19, QSI 2020-21, Prairie Station UAV-lidar,
Mores Creek Summit lidar) plus documented ASO SnowEx 2017/2020 flights absent from CMR.
"""

import json                    # cache read/write
import os                      # cache existence check
import requests                # CMR REST queries
from shapely.geometry import Polygon, box  # footprint overlap tests
from find_snowex_flights import load_flights, SNOWEX_SITES  # cached UAVSAR footprints + site names

CACHE = "cmr_lidar_cache.json"                                  # local cache of CMR granule metadata
CMR = "https://cmr.earthdata.nasa.gov/search/granules.json"     # CMR granule search endpoint

# lidar collections to query at CMR (short name -> dataset label)
COLLECTIONS = {"ASO_50M_SWE": "ASO lidar (NSIDC 2013-19)",
               "SNEX20_QSI_SD": "QSI lidar 0.5m/3m (SNEX20_QSI_SD[_3m])",
               "SNEX20_GM_Lidar": "QSI lidar Grand Mesa IOP (SNEX20_GM_Lidar)",
               "SNEX21_PS_DSM": "Prairie Station UAV-lidar (SNEX21_PS_DSM)",
               "SNEX_MCS_Lidar": "Mores Creek Summit lidar (SNEX_MCS_Lidar)"}

# site/basin code -> (name, state) for codes appearing in granule filenames
LIDAR_SITES = {"USIDMC": ("Mores Creek Summit", "ID"), "USIDBS": ("Banner Summit", "ID"),
               "USIDDC": ("Dry Creek", "ID"), "USCOFR": ("Fraser Experimental Forest", "CO"),
               "USCOCP": ("Cameron Pass", "CO"), "USUTLC": ("Little Cottonwood Canyon", "UT"),
               "USCATB": ("Tuolumne Basin", "CA"), "USCAMB": ("Merced Basin", "CA"),
               "USCASJ": ("San Joaquin Basin", "CA"), "USCOUB": ("Uncompahgre Basin", "CO"),
               "USCOCB": ("Crested Butte / East River", "CO"), "USCOGE": ("Gunnison East / East River", "CO"),
               "USCOGT": ("Gunnison Taylor", "CO"), "USCORG": ("Rio Grande", "CO"),
               "USCOCJ": ("Conejos", "CO"), "USCOBR": ("Blue River", "CO"),
               "USCOCM": ("Castle/Maroon", "CO"), "USWAOL": ("Olympic Mtns", "WA"),
               "USCACE": ("Cherry-Eleanor", "CA"), "USCAJW": ("San Joaquin West (JW)", "CA"),
               "USCAKC": ("Kings Canyon", "CA"), "USCAKN": ("Kern", "CA"),
               "USCAKW": ("Kaweah", "CA"), "USCALB": ("Lakes Basin (Mammoth)", "CA"),
               "USCALV": ("Lee Vining", "CA"), "USCARC": ("Rush Creek", "CA"),
               "USCASF": ("San Joaquin South Fork (SF)", "CA"), "USCATE": ("Tuolumne East (TE)", "CA"),
               "CASC": ("Central Agricultural Research Center (CARC)", "MT"),
               "MCS": ("Mores Creek Summit", "ID")}


def cmr_granules(short_name):
    """Return all granule entries for a CMR collection (paged, public API)."""
    out, page = [], 1                                           # results and page counter
    while True:                                                 # page until fewer than page_size returned
        r = requests.get(CMR, params={"short_name": short_name, "page_size": 200, "page_num": page})
        entries = r.json()["feed"]["entry"]                     # granule entries on this page
        out += entries                                          # accumulate
        if len(entries) < 200:                                  # last page reached
            return out
        page += 1                                               # next page


def load_lidar():
    """Return raw CMR granule metadata for all lidar collections, cached locally."""
    if os.path.exists(CACHE):                                   # cache hit: skip remote queries
        return json.load(open(CACHE))
    data = {sn: cmr_granules(sn) for sn in COLLECTIONS}         # query every collection
    json.dump(data, open(CACHE, "w"))                           # save cache for future runs
    return data


def poly_from_entry(e):
    """Build a shapely polygon from a CMR granule's polygons (lat lon pairs) or bounding box."""
    if e.get("polygons"):                                       # preferred: real footprint polygon
        nums = [float(x) for x in e["polygons"][0][0].split()]  # flat list: lat lon lat lon ...
        return Polygon([(nums[i+1], nums[i]) for i in range(0, len(nums), 2)])  # -> (lon, lat) ring
    s, w, n, ee = [float(x) for x in e["boxes"][0].split()]     # fallback: S W N E box
    return box(w, s, ee, n)                                     # box as polygon


lidar_raw = load_lidar()                                        # all granule metadata (cached)

acqs = {}                                                       # (dataset, site, date) -> acquisition record
for sn, entries in lidar_raw.items():                           # walk each collection
    for e in entries:                                           # each granule
        gid = e.get("producer_granule_id") or e["title"]        # granule filename
        date = e["time_start"][:10]                             # acquisition date (ISO)
        if int(date[5:7]) not in (1, 2, 3, 4):                  # keep January-April only
            continue
        site = next((c for c in LIDAR_SITES if c in gid), "?")  # site code from filename
        if sn == "SNEX20_GM_Lidar":                             # GM IOP granules carry no site code
            site = "USCOGM"
        note = "snow-off" if "snowoff" in gid else ""           # flag snow-off DSMs
        if "orthomosaic" in gid:                                # skip non-lidar orthomosaic product
            continue
        poly = poly_from_entry(e)                               # granule footprint
        if "USCOUB" in gid:                                     # CMR bug: USCOUB granule carries a CA polygon
            poly = box(-107.9, 37.85, -107.55, 38.1)            # approx upper Uncompahgre (Senator Beck area)
            note = "approx footprint (bad CMR metadata)"
        acqs.setdefault((sn, site, date), {"poly": poly, "note": note})  # dedupe resolutions

def zip_footprint(zip_glob):
    """Return the WGS84 footprint of a snow-depth GeoTIFF inside a downloaded ASO zip, or None."""
    import glob                                                 # find the local zip
    zips = glob.glob(os.path.join("ASO_data", zip_glob))        # match the product zip
    if not zips:                                                # data not downloaded -> caller falls back
        return None
    import zipfile, rasterio                                    # read tif header inside the zip
    from rasterio.warp import transform_bounds                  # UTM -> WGS84 bounds
    tifs = [n for n in zipfile.ZipFile(zips[0]).namelist()      # prefer a 50m snow-depth tif (small)
            if n.endswith(".tif") and "swe" in n.lower()] or \
           [n for n in zipfile.ZipFile(zips[0]).namelist() if n.endswith(".tif")]
    with rasterio.open(f"/vsizip/{zips[0]}/{tifs[0]}") as src:  # open without extracting
        return box(*transform_bounds(src.crs, "EPSG:4326", *src.bounds))  # footprint in lon/lat


# documented ASO SnowEx flights absent from CMR: (site code, date, local zip pattern, fallback poly, note)
gm_poly = poly_from_entry(lidar_raw["SNEX20_GM_Lidar"][0])      # Grand Mesa footprint (from QSI granule)
aso_poly = {s: poly_from_entry(e) for sn in ["ASO_50M_SWE"] for e in lidar_raw[sn]  # ASO basin footprints
            for s in [next((c for c in LIDAR_SITES if c in (e.get("producer_granule_id") or "")), "?")]}
rcew = box(-116.87, 43.03, -116.63, 43.27)                      # approximate Reynolds Creek watershed box
MANUAL = [("USCOGM", "2017-02-08", None, gm_poly, "ASO SnowEx17"),
          ("USCOGM", "2017-02-16", None, gm_poly, "ASO SnowEx17"),
          ("USCOGM", "2017-02-25", None, gm_poly, "ASO SnowEx17"),
          ("USCOGM", "2020-02-01", "*GrandMesa*2020Feb1-2*.zip", gm_poly, "ASO SnowEx20 (flown Feb 1-2)"),
          ("USCOGM", "2020-02-13", "*GrandMesa*2020Feb13*.zip", gm_poly, "ASO SnowEx20"),
          ("USCOCB", "2020-02-14", "*EastRiver*2020Feb14-20*.zip", aso_poly["USCOCB"], "ASO SnowEx20 (East River, flown Feb 14-20)"),
          ("USCOGT", "2020-02-20", "*TaylorRiver*2020Feb20*.zip", aso_poly["USCOGT"], "ASO SnowEx20"),
          ("USIDRC", "2020-02-18", "*Reynolds*2020Feb18-19*.zip", rcew, "ASO SnowEx20 (Reynolds Creek, flown Feb 18-19)")]
MANUAL = [(s, d, zip_footprint(z) if z else None, f, n) for s, d, z, f, n in MANUAL]  # try local data first
MANUAL = [(s, d, p if p is not None else f,                     # real footprint, else fallback
           n + ("" if p is not None else "; approx footprint")) for s, d, p, f, n in MANUAL]
LIDAR_SITES["USCOGM"] = ("Grand Mesa", "CO")                    # add codes used only by manual entries
LIDAR_SITES["USIDRC"] = ("Reynolds Creek", "ID")
for site, date, poly, note in MANUAL:                           # add manual flights (footprints approximate)
    acqs[("ASO_SnowEx", site, date)] = {"poly": poly, "note": note}

uav = {}                                                        # UAVSAR (campaign, line_id) -> shapely polygon
for f in load_flights():                                        # walk cached UAVSAR flights
    if f["campaign"] in SNOWEX_SITES:                           # catalogued lines only
        uav.setdefault((f["campaign"], f["line_id"]), Polygon(f["geometry"]["coordinates"][0]))

rows, dropped = [], set()                                       # kept acquisitions and dropped sites
for (sn, site, date), a in sorted(acqs.items()):                # test every lidar acquisition
    hits = sorted({c for (c, l), p in uav.items() if p.intersects(a["poly"])})  # overlapping UAVSAR campaigns
    if not hits:                                                # no UAVSAR overlap -> exclude
        dropped.add(site)
        continue
    name, state = LIDAR_SITES.get(site, (site, "?"))            # readable site name
    rows.append({"location": name, "state": state, "date": date,  # record for table/KML
                 "dataset": COLLECTIONS.get(sn, a["note"] or sn), "note": a["note"],
                 "uavsar_overlap": ",".join(hits), "site": site,
                 "flag": "" if "2014" <= date <= "2021-12-31" else "outside 2014-2021",
                 "poly": a["poly"]})

import pandas as pd                                             # tabular output
df = pd.DataFrame([{k: v for k, v in r.items() if k != "poly"} for r in rows])  # drop geometry column
df.to_csv("snowex_lidar_flights.csv", index=False)              # save the lidar flight table
print(f"{len(df)} lidar acquisitions kept; sites with no UAVSAR overlap: {sorted(dropped)}")

kml = ['<?xml version="1.0" encoding="UTF-8"?>',                # KML header
       '<kml xmlns="http://www.opengis.net/kml/2.2"><Document>',
       '<name>SnowEx lidar acquisitions</name>']
for code in sorted({r["site"] for r in rows}):                  # one folder per lidar site
    name, state = LIDAR_SITES.get(code, (code, "?"))
    kml.append(f'<Folder><name>{name}, {state} (lidar)</name>')
    for r in [r for r in rows if r["site"] == code]:            # each acquisition of this site
        coords = " ".join(f"{x},{y},0" for x, y in r["poly"].exterior.coords)  # KML coordinate string
        extra = " ".join(t for t in [r["note"], r["flag"]] if t)  # optional annotations
        kml.append(                                             # outline-heavy placemark, toggleable
            f'<Placemark><name>{name} {r["date"]}</name>'
            f'<description>{r["dataset"]} | overlaps UAVSAR: {r["uavsar_overlap"]}'
            f'{" | " + extra if extra else ""}</description>'
            f'<Style><LineStyle><color>ffff00ff</color><width>3</width></LineStyle>'
            f'<PolyStyle><color>2dff00ff</color></PolyStyle></Style>'
            f'<Polygon><outerBoundaryIs><LinearRing><coordinates>{coords}</coordinates>'
            f'</LinearRing></outerBoundaryIs></Polygon></Placemark>')
    kml.append('</Folder>')                                     # close site folder
kml.append('</Document></kml>')                                 # close document
open("snowex_lidar.kml", "w").write("\n".join(kml))             # write the KML
print(f"{len(rows)} placemarks -> snowex_lidar.kml")
