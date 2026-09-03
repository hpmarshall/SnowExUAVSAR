"""Build a KML of lidar/UAVSAR overlap: the spatial intersection of each lidar
footprint with each UAVSAR flight line, kept only where the two were acquired
within 3 months (92 days) of each other.
"""

from datetime import date                                       # date arithmetic for the 3-month test
import find_snowex_lidar as fsl                                 # runs the lidar catalog (cached); gives rows + uav polys
from find_snowex_flights import load_flights                    # cached UAVSAR flights (for dates per line)

MAX_DAYS = 92                                                   # "within 3 months"

d = lambda s: date(*map(int, s.split("-")))                     # ISO string -> date

uav_dates = {}                                                  # (campaign, line_id) -> sorted flight dates
for f in load_flights():                                        # walk cached UAVSAR flights
    if f["campaign"] in fsl.SNOWEX_SITES:                       # catalogued lines only
        uav_dates.setdefault((f["campaign"], f["line_id"]), set()).add(f["date"])

pairs = {}                                                      # (lidar site, campaign, line) -> overlap record
for r in fsl.rows:                                              # every lidar acquisition (with polygon)
    for (camp, lid), poly in fsl.uav.items():                   # every UAVSAR flight line
        if not poly.intersects(r["poly"]):                      # must overlap spatially
            continue
        close = sorted(u for u in uav_dates[(camp, lid)]        # UAVSAR dates within 3 months of the lidar date
                       if abs((d(u) - d(r["date"])).days) <= MAX_DAYS)
        if not close:                                           # no temporally close UAVSAR flight -> skip
            continue
        key = (r["site"], camp, lid)                            # one placemark per site/line pairing
        rec = pairs.setdefault(key, {"geom": poly.intersection(r["poly"]),  # intersection polygon
                                     "matches": [], "loc": r["location"], "state": r["state"]})
        rec["matches"].append(f'lidar {r["date"]} ~ UAVSAR {", ".join(close)}')  # record the date pairing

def rings(geom):
    """Yield exterior rings for a Polygon or MultiPolygon."""
    for g in getattr(geom, "geoms", [geom]):                    # MultiPolygon -> parts; Polygon -> itself
        yield g.exterior.coords                                 # outer ring coordinates

kml = ['<?xml version="1.0" encoding="UTF-8"?>',                # KML header
       '<kml xmlns="http://www.opengis.net/kml/2.2"><Document>',
       '<name>Lidar-UAVSAR overlap (within 3 months)</name>']
for (site, camp, lid), rec in sorted(pairs.items()):            # one placemark per lidar-site x UAVSAR-line pair
    polys = "".join(                                            # KML polygon(s) for the intersection
        f'<Polygon><outerBoundaryIs><LinearRing><coordinates>'
        f'{" ".join(f"{x},{y},0" for x, y in ring)}'
        f'</coordinates></LinearRing></outerBoundaryIs></Polygon>' for ring in rings(rec["geom"]))
    kml.append(                                                 # orange filled overlap region, toggleable
        f'<Placemark><name>{rec["loc"]} x {camp} {lid[:3]}° ({lid})</name>'
        f'<description>{rec["loc"]}, {rec["state"]} lidar ∩ UAVSAR {camp} line {lid} | '
        f'{"; ".join(rec["matches"])}</description>'
        f'<Style><LineStyle><color>ff0080ff</color><width>2</width></LineStyle>'
        f'<PolyStyle><color>7f0080ff</color></PolyStyle></Style>'
        f'<MultiGeometry>{polys}</MultiGeometry></Placemark>')
kml.append('</Document></kml>')                                 # close document
open("snowex_overlap.kml", "w").write("\n".join(kml))           # write the KML
print(f"{len(pairs)} overlap regions -> snowex_overlap.kml")

import pandas as pd                                             # tabular output of the same regions
rows = []                                                       # one row per overlap placemark
for (site, camp, lid), rec in sorted(pairs.items()):            # walk the KML placemarks
    lidar_dates = sorted({m.split(" ~ ")[0].replace("lidar ", "") for m in rec["matches"]})  # lidar side
    uav_close = sorted({u for m in rec["matches"] for u in m.split("UAVSAR ")[1].split(", ")})  # radar side
    rows.append({"location": rec["loc"], "state": rec["state"],  # table row
                 "uavsar_line": f'{camp} {lid} ({int(lid[:3])} deg)',
                 "lidar_dates": " ".join(lidar_dates), "uavsar_dates_within_3mo": " ".join(uav_close)})
df = pd.DataFrame(rows)                                         # assemble the table
df.to_csv("snowex_overlap.csv", index=False)                    # CSV version
with open("overlap_table.md", "w") as f:                        # markdown version
    f.write("# Lidar-UAVSAR Overlap Regions (within 3 months)\n\n")
    f.write(f"{len(df)} overlap regions in snowex_overlap.kml: each is the spatial intersection of a\n")
    f.write("lidar footprint and a UAVSAR flight line where acquisitions are within 92 days.\n\n")
    f.write(df.to_markdown(index=False))
print(f"{len(df)} rows -> snowex_overlap.csv, overlap_table.md")
