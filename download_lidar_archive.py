"""Download a local copy of every lidar acquisition in snowex_lidar_flights.csv.

Resumable by design: each file streams to a .part temp name and is only renamed to
its final name once the byte count matches the expected size. A killed/interrupted
run (laptop sleep, network change, Ctrl-C) can just be re-started -- complete files
are skipped, partial ones are redone. Safe to run from cron/background repeatedly.
"""

import os                      # paths, sizes
import time                    # retry backoff
import datetime                # pad the temporal search window
import requests                # streaming HTTP downloads
import pandas as pd            # read the lidar catalog
import earthaccess             # NSIDC/CMR authenticated downloads

ROOT = "LIDAR"                                                  # all lidar data lives here (git-ignored)
CMR_ROOT = os.path.join(ROOT, "CMR")                             # NSIDC/CMR-hosted granules, by collection
BUCKET_ROOT = os.path.join(ROOT, "ASO_bucket")                   # public ASO S3 zips, by campaign

earthaccess.login(strategy="netrc")                              # authenticate once for all CMR downloads
EDL_SESSION = earthaccess.get_requests_https_session()           # reused authenticated session

# dataset label (as it appears in snowex_lidar_flights.csv) -> CMR short_name
DATASET_TO_SN = {
    "ASO lidar (NSIDC 2013-19)": "ASO_50M_SWE",
    "Mores Creek Summit lidar (SNEX_MCS_Lidar)": "SNEX_MCS_Lidar",
    "Prairie Station UAV-lidar (SNEX21_PS_DSM)": "SNEX21_PS_DSM",
    "QSI lidar 0.5m/3m (SNEX20_QSI_SD[_3m])": "SNEX20_QSI_SD",    # 3m sibling handled separately below
    "QSI lidar Grand Mesa IOP (SNEX20_GM_Lidar)": "SNEX20_GM_Lidar",
    "ASO snow-off DTM (bare-earth reference, ASO_3M_PCDTM)": "ASO_3M_PCDTM",
}

# ASO public S3 bucket: (site code, date) -> full "AllData_and_Reports" zip URL
BUCKET_BASE = "https://asopublic.s3-us-west-1.amazonaws.com"
BUCKET_URLS = {
    ("USCOGM", "2020-02-01"): f"{BUCKET_BASE}/USCO/GM/2020/0201/ASO_GrandMesa_mosaic_2020Feb1-2_AllData_and_Reports.zip",
    ("USCOGM", "2020-02-13"): f"{BUCKET_BASE}/USCO/GM/2020/0213/ASO_GrandMesa_mosaic_2020Feb13_AllData_and_Reports.zip",
    ("USCOCB", "2020-02-14"): f"{BUCKET_BASE}/USCO/GE/2020/0214/ASO_EastRiver_mosaic_2020Feb14-20_AllData_and_Reports.zip",
    ("USCOGT", "2020-02-20"): f"{BUCKET_BASE}/USCO/GT/2020/0220/ASO_TaylorRiver_mosaic_2020Feb20_AllData_and_Reports.zip",
    ("USIDRC", "2020-02-18"): f"{BUCKET_BASE}/USID/RY/2020/0218/ASO_Reynolds_mosaic_2020Feb18-19_AllData_and_Reports.zip",
    ("USCOUB", "2015-04-30"): f"{BUCKET_BASE}/USCO/UB/2015/0430/ASO_Uncompahgre_mosaic_2015Apr30_AllData_and_Reports.zip",
    ("USCOAN", "2021-04-19"): f"{BUCKET_BASE}/USCO/AN/2021/0419/ASO_Animas_mosaic_2021Apr19_AllData_and_Reports.zip",
    ("USCODL", "2021-04-20"): f"{BUCKET_BASE}/USCO/DL/2021/0420/ASO_Dolores_mosaic_2021Apr20-21_AllData_and_Reports.zip",
}

# rows with no downloadable product at all (already archived locally by the user)
LOCAL_ONLY = {("USCOSB", "2017-02-16"), ("USCOUB", "2015-04-29"), ("USCOUB", "2016-06-04")}


def download(url, dest, expected_mb=None, exact_bytes=None, session=None, retries=4):
    """Stream url to dest atomically. Skip if already complete; retry with backoff on failure.

    CMR's catalog `granule_size` can be stale for a given object (observed: several ASO_3M_SD
    granules report a size a few % larger than the file the server actually sends -- a genuine
    metadata/data mismatch, not a truncated transfer). So it is only used for a loose pre-download
    skip check; the authoritative size for validating an actual transfer is the live response's
    own Content-Length header. Pass exact_bytes (e.g. from a HEAD request) instead of expected_mb
    when the precise per-file size is already known (multi-file granules), for a tight skip check."""
    os.makedirs(os.path.dirname(dest), exist_ok=True)            # ensure the destination folder exists
    if exact_bytes is not None:                                   # precise size known -- tight tolerance
        catalog_bytes, loose_tol = exact_bytes, 1000
    else:                                                          # CMR/earthaccess report granule size in MiB
        catalog_bytes = int(expected_mb * 1024 * 1024) if expected_mb else None
        loose_tol = max(5_000_000, int((catalog_bytes or 0) * 0.10))  # generous -- catalog size is only a hint
    if os.path.exists(dest):                                      # already have a final (non-.part) file
        if catalog_bytes is None or abs(os.path.getsize(dest) - catalog_bytes) < loose_tol:
            return "exists"                                        # roughly matches (or unknown) -> trust it
    part = dest + ".part"                                          # temp name while streaming
    s = session or requests                                        # authenticated session for CMR, else plain
    for attempt in range(retries):                                 # retry transient network failures
        try:
            with s.get(url, stream=True, timeout=120) as r:
                r.raise_for_status()
                content_length = r.headers.get("content-length")   # authoritative size for THIS transfer
                server_bytes = int(content_length) if content_length else None
                with open(part, "wb") as f:                         # stream to disk, never fully in memory
                    for chunk in r.iter_content(chunk_size=1 << 20):
                        f.write(chunk)
                got = os.path.getsize(part)
                if server_bytes is not None:                        # trust the server over stale catalog metadata
                    if got != server_bytes:
                        raise IOError(f"truncated: got {got}, server said {server_bytes}")
                elif catalog_bytes and abs(got - catalog_bytes) > loose_tol:  # no header -- fall back to catalog
                    raise IOError(f"size mismatch: got {got}, expected ~{catalog_bytes}")
                os.replace(part, dest)                               # atomic: only now does dest "exist"
                return "downloaded"
        except Exception as e:                                      # network hiccup, sleep/wake, etc.
            print(f"    retry {attempt + 1}/{retries} for {os.path.basename(dest)}: {e}")
            time.sleep(2 ** attempt)                                 # exponential backoff
    return "failed"


def cmr_granule(short_name, site, date, retries=4):
    """Find the earthaccess granule for a (site, date) in a CMR collection, or None.

    Retries on transient network failure (DNS drop, connection reset, etc.) instead
    of letting the whole multi-hour run die on one flaky request."""
    ymd = date.replace("-", "")                                    # YYYY-MM-DD -> YYYYMMDD
    d = datetime.date.fromisoformat(date)                           # a same-day (date, date) range matches
    window = (str(d), str(d + datetime.timedelta(days=1)))          # nothing in CMR/earthaccess -- pad by 1 day
    for attempt in range(retries):
        try:
            results = earthaccess.search_data(short_name=short_name, temporal=window)
            if short_name == "SNEX20_GM_Lidar":                       # single-site collection; filenames don't
                by_date = [g for g in results if ymd in g["umm"]["GranuleUR"]]  # contain "USCOGM" at all
                full = [g for g in by_date if "subset" not in g["umm"]["GranuleUR"].lower()]  # prefer the
                return (full or by_date or [None])[0]                # full SD product over the smaller SDsubset
            exact = [g for g in results if site in g["umm"]["GranuleUR"] and ymd in g["umm"]["GranuleUR"]]
            if exact:                                                # site code + date both in the filename
                return exact[0]
            by_site = [g for g in results if site in g["umm"]["GranuleUR"]]  # e.g. "..._snowoff_V01.0.tif"
            if len(by_site) == 1:                                    # no date in the name, but CMR's own
                return by_site[0]                                    # temporal window already narrowed to 1
            return None                                             # search succeeded, just no match
        except Exception as e:                                      # network drop during the search itself
            print(f"    retry {attempt + 1}/{retries} for CMR search {short_name} {site} {date}: {e}")
            time.sleep(2 ** attempt)
    print(f"    giving up on CMR search {short_name} {site} {date} after {retries} attempts")
    return "search_failed"                                          # distinguish from a genuine "not found"


counts = {"downloaded": 0, "exists": 0, "failed": 0, "skipped_local": 0, "skipped_no_product": 0, "row_errors": 0}
df = pd.read_csv("snowex_lidar_flights.csv")                        # the full lidar catalog to archive

for _, row in df.iterrows():                                        # walk every cataloged acquisition
    key = (row.site, row.date)                                      # (site code, ISO date)

    if key in LOCAL_ONLY:                                            # already archived by the user, no download
        counts["skipped_local"] += 1
        continue

    try:                                                             # nothing about one row may kill the whole run
        # (site, date) is NOT a unique key -- e.g. Grand Mesa 2020-02-01 has both a QSI CMR row and a
        # separate ASO bucket row for the same day. Resolve the CMR/dataset mapping first so it can
        # never be shadowed by an unrelated bucket entry that happens to share the same (site, date).
        sn = DATASET_TO_SN.get(row.dataset)                          # CMR-hosted row?
        if sn is None:                                               # e.g. GM 2017 proxy rows, handled below
            sn = "ASO_3M_SD" if "SnowEx17" in str(row.dataset) and row.site == "USCOGM" else None

        if sn is None and key in BUCKET_URLS:                        # ASO public S3, no auth needed
            url = BUCKET_URLS[key]
            dest = os.path.join(BUCKET_ROOT, os.path.basename(url))
            print(f"[bucket] {row.location} {row.date}: {os.path.basename(url)}")
            counts[download(url, dest)] += 1
            continue
        if sn is None:
            counts["skipped_no_product"] += 1
            print(f"[skip] {row.location} {row.date}: no known downloadable source ({row.dataset})")
            continue

        g = cmr_granule(sn, row.site, row.date)
        if g == "search_failed":                                     # CMR unreachable after retries -- try again later
            counts["row_errors"] += 1
            continue
        if g is None:
            counts["skipped_no_product"] += 1
            print(f"[missing] {row.location} {row.date}: no {sn} granule found in CMR")
            continue
        links = g.data_links()                                       # a granule can bundle >1 file (e.g. MCS Lidar:
        if len(links) > 1:                                            # CHM/DSM/DTM/SD are 4 separate products)
            print(f"[CMR:{sn}] {row.location} {row.date}: {g['umm']['GranuleUR']} ({len(links)} files)")
            for url in links:                                          # each file gets its own destination + size
                fname = url.rsplit("/", 1)[-1]
                dest = os.path.join(CMR_ROOT, sn, fname)
                head = EDL_SESSION.head(url, allow_redirects=True, timeout=30)
                exact = int(head.headers["content-length"]) if "content-length" in head.headers else None
                print(f"  [{fname}] ~{(exact or 0) / 1e6:.0f} MB")
                counts[download(url, dest, exact_bytes=exact, session=EDL_SESSION)] += 1
        else:
            dest = os.path.join(CMR_ROOT, sn, g["umm"]["GranuleUR"])
            print(f"[CMR:{sn}] {row.location} {row.date}: {g['umm']['GranuleUR']} (~{g.size():.0f} MB)")
            counts[download(links[0], dest, expected_mb=g.size(), session=EDL_SESSION)] += 1

        if sn == "SNEX20_QSI_SD":                                    # also fetch the 3m sibling product
            g3 = cmr_granule("SNEX20_QSI_SD_3m", row.site, row.date)
            if g3 and g3 != "search_failed":
                dest3 = os.path.join(CMR_ROOT, "SNEX20_QSI_SD_3m", g3["umm"]["GranuleUR"])
                counts[download(g3.data_links()[0], dest3, expected_mb=g3.size(), session=EDL_SESSION)] += 1

        if row.dataset == "ASO lidar (NSIDC 2013-19)":               # also fetch the full-res 3m companion
            g3 = cmr_granule("ASO_3M_SD", row.site, row.date)
            if g3 and g3 != "search_failed":
                dest3 = os.path.join(CMR_ROOT, "ASO_3M_SD", g3["umm"]["GranuleUR"])
                print(f"  [3m companion] {g3['umm']['GranuleUR']} (~{g3.size():.0f} MB)")
                counts[download(g3.data_links()[0], dest3, expected_mb=g3.size(), session=EDL_SESSION)] += 1
            elif g3 is None:
                print(f"  [3m companion] none published for {row.site} {row.date}")
    except Exception as e:                                           # unexpected failure -- log and keep going
        counts["row_errors"] += 1
        print(f"[error] {row.location} {row.date}: {e}")

print("\n=== summary ===")
for k, v in counts.items():
    print(f"{k}: {v}")
if counts["row_errors"]:
    print(f"\n{counts['row_errors']} row(s) hit a network/error interruption -- just re-run this script to retry them.")
