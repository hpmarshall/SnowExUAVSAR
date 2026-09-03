"""Build a Google Earth KML of every SnowEx UAVSAR flight-line footprint.

Each flight line is a separate Placemark (individually toggleable), grouped in one
Folder per SnowEx site, using the actual ASF swath footprint polygons.
"""

from find_snowex_flights import load_flights, SNOWEX_SITES  # cached ASF query + site lookup

# one distinct KML color (aabbggrr hex) per state so overlapping sites are distinguishable
STATE_COLORS = {"CO": "7f0000ff", "ID": "7f00aaff", "UT": "7f00ff00",
                "MT": "7fff0000", "NM": "7fff00ff", "CA": "7fffff00"}

lines = {}                                                      # unique flight lines keyed by (campaign, line_id)
for f in load_flights():                                        # walk cached flights
    if f["campaign"] not in SNOWEX_SITES:                       # keep SnowEx campaigns only
        continue
    key = (f["campaign"], f["line_id"])                         # one entry per flight line
    lines.setdefault(key, {**f, "dates": []})["dates"].append(f["date"])  # collect acquisition dates

kml = ['<?xml version="1.0" encoding="UTF-8"?>',                # KML header
       '<kml xmlns="http://www.opengis.net/kml/2.2"><Document>',
       '<name>SnowEx UAVSAR flight lines</name>']
for camp in sorted({c for c, _ in lines}):                      # one folder per SnowEx site
    site, state = SNOWEX_SITES[camp]                            # site name and state
    kml.append(f'<Folder><name>{site}, {state}</name>')         # open the site folder
    for (c, lid), f in sorted(lines.items()):                   # each flight line of this site
        if c != camp:                                           # skip other sites' lines
            continue
        ring = f["geometry"]["coordinates"][0]                  # outer ring of the GeoJSON footprint
        coords = " ".join(f"{lon},{lat},0" for lon, lat in ring)  # KML lon,lat,alt coordinate string
        dates = ", ".join(sorted(set(f["dates"])))              # acquisition dates for the description
        kml.append(                                             # one toggleable placemark per line
            f'<Placemark><name>{site} {f["bearing_deg"]:03d}° ({lid})</name>'
            f'<description>State: {state} | Campaign: {c} | Line: {lid} | '
            f'{len(set(f["dates"]))} acquisitions: {dates}</description>'
            f'<Style><LineStyle><color>ff{STATE_COLORS[state][2:]}</color><width>2</width></LineStyle>'
            f'<PolyStyle><color>{STATE_COLORS[state]}</color></PolyStyle></Style>'
            f'<Polygon><outerBoundaryIs><LinearRing><coordinates>{coords}</coordinates>'
            f'</LinearRing></outerBoundaryIs></Polygon></Placemark>')
    kml.append('</Folder>')                                     # close the site folder
kml.append('</Document></kml>')                                 # close the document

open("snowex_flightlines.kml", "w").write("\n".join(kml))       # write the KML file
print(f"{len(lines)} flight lines from {len({c for c, _ in lines})} sites -> snowex_flightlines.kml")
