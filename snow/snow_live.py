"""Real snow on the ground for the web twin (served by 05_WebTwin/serve.py as GET /api/snow).

- NOAA NOHRSC SNODAS (NSIDC G02158, no login): today's 1 km snow depth and snow water equivalent grid. Only the two
  members we need are read from the daily tar with HTTP range requests (about 1.6 MB, a few MB in winter), and only the
  rows over the park are decompressed. New grids are posted around 13:00 UTC; until then yesterday's is used.
- Sentinel-2 L2A scene classification (ESA Copernicus via Microsoft Planetary Computer, no login): where the satellite saw snow,
  bare ground or nothing (cloud, shadow) at 20 m, resampled by the Planetary Computer onto the twin's own 30 m grid. Newest
  clear pass of the last 10 days (every few days). Written to 05_WebTwin/data/live/s2_snow_<date>.png: 0 unknown, 128 bare, 255 snow.
  SNODAS for that day is included too, so the browser can tell snow that fell after the pass.
- CDEC (California Data Exchange Center) snow pillows and depth sensors in and around the twin: hourly readings.
  Each sensor's September median (the ground is bare then) is subtracted, since several read 1-3 inches on bare ground.
Results are cached in 00_Pipeline/cache (SNODAS per day, Sentinel-2 for 6 hours, stations for an hour)."""
import datetime as dt, io, json, math, os, statistics, threading, time, urllib.request, zlib

PIPE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(PIPE, "cache")
E0, N0 = 288000.0, 4189000.0                       # the twin's origin (UTM 11N): x = E - E0, z = -(N - N0)
SNODAS = "https://noaadata.apps.nsidc.org/NOAA/G02158/masked/{y}/{m:02d}_{mon}/SNODAS_{y}{m:02d}{d:02d}.tar"
GRID = {"cols": 6935, "rows": 3351, "lon0": -124.73333333333333, "lat0": 52.875, "step": 1 / 120}   # masked grid, NW corner
BOX = {"lon": (-119.75, -119.15), "lat": (37.66, 38.00)}          # the twin and a little around it
CDEC = "https://cdec.water.ca.gov/dynamicapp/req/JSONDataServlet?Stations={s}&SensorNums={n}&dur_code={d}&Start={a}&End={b}"
STATIONS = [   # id, name, lat, lon, elevation m (CDEC station pages); all report hourly snow depth (sensor 18)
    ("TUM", "Tuolumne Meadows", 37.876215, -119.349968, 2591),
    ("DAN", "Dana Meadows", 37.898045, -119.258808, 2975),
    ("TES", "Tioga Pass entrance", 37.910870, -119.258507, 3031),
    ("TNY", "Tenaya Lake", 37.837250, -119.450185, 2460),
    ("WHW", "White Wolf", 37.859445, -119.651607, 2408),
    ("FLV", "Vogelsang", 37.794571, -119.347404, 3077),
    ("ERY", "Ellery Lake", 37.934967, -119.232201, 2940),
]
PC_STAC = "https://planetarycomputer.microsoft.com/api/stac/v1/search"
PC_GRID = ("https://planetarycomputer.microsoft.com/api/data/v1/item/bbox/266000,4177000,305000,4201000/1300x800.npy"
           "?collection=sentinel-2-l2a&item={id}&assets=SCL&coord_crs=EPSG:32611&dst_crs=EPSG:32611&resampling=nearest")
LIVE = os.path.join(os.path.dirname(PIPE), "05_WebTwin", "data", "live")
UA = {"User-Agent": "YosemiteDigitalTwin/1.0 (local research viewer)"}
_lock = threading.Lock()
_state = {"snodas": None, "stations": None, "stations_at": 0.0, "s2": None, "s2_at": 0.0}


def _get(url, rng=None, timeout=60):
    h = dict(UA)
    if rng:
        h["Range"] = f"bytes={rng[0]}-{rng[1]}"
    with urllib.request.urlopen(urllib.request.Request(url, headers=h), timeout=timeout) as r:
        return r.read()


def latlon_to_xz(lat, lon):
    """WGS84/NAD83 lat, lon -> the twin's x, z (UTM zone 11N, Krueger series; about 1 mm accurate here)."""
    a, f, k0 = 6378137.0, 1 / 298.257222101, 0.9996
    n = f / (2 - f)
    A = a / (1 + n) * (1 + n * n / 4 + n ** 4 / 64)
    al = [n / 2 - 2 * n * n / 3 + 5 * n ** 3 / 16, 13 * n * n / 48 - 3 * n ** 3 / 5, 61 * n ** 3 / 240]
    phi, lam = math.radians(lat), math.radians(lon + 117)
    e = 2 * math.sqrt(n) / (1 + n)
    t = math.sinh(math.atanh(math.sin(phi)) - e * math.atanh(e * math.sin(phi)))
    xi, eta = math.atan2(t, math.cos(lam)), math.atanh(math.sin(lam) / math.sqrt(1 + t * t))
    E = 500000 + k0 * A * (eta + sum(al[j] * math.cos(2 * (j + 1) * xi) * math.sinh(2 * (j + 1) * eta) for j in range(3)))
    N = k0 * A * (xi + sum(al[j] * math.sin(2 * (j + 1) * xi) * math.cosh(2 * (j + 1) * eta) for j in range(3)))
    return E - E0, -(N - N0)


# ---------------------------------------------------------------------------------------------- SNODAS
def _members(url):
    """Walk the tar's 512-byte headers with range requests: {name: (data offset, size)}."""
    out, off = {}, 0
    while True:
        h = _get(url, (off, off + 511))
        if len(h) < 512 or not h[:100].strip(b"\0"):
            return out
        name = h[:100].split(b"\0")[0].decode()
        size = int(h[124:136].split(b"\0")[0].strip() or b"0", 8)
        out[name] = (off + 512, size)
        off += 512 + (size + 511) // 512 * 512


def _rows(blob, r0, r1):
    """Decompress a gzip'ed big-endian int16 grid only up to row r1 and return rows r0..r1 (inclusive)."""
    need, d, buf = (r1 + 1) * GRID["cols"] * 2, zlib.decompressobj(16 + zlib.MAX_WBITS), bytearray()
    for i in range(0, len(blob), 1 << 16):
        buf += d.decompress(blob[i:i + (1 << 16)])
        if len(buf) >= need:
            break
    return memoryview(bytes(buf))[r0 * GRID["cols"] * 2:need]


def _snodas_day(day):
    url = SNODAS.format(y=day.year, m=day.month, mon=day.strftime("%b"), d=day.day)
    mem = _members(url)
    pick = lambda code: next(k for k in mem if code in k and k.endswith(".dat.gz"))
    s, g = GRID["step"], GRID
    c0, c1 = int((BOX["lon"][0] - g["lon0"]) / s), int(math.ceil((BOX["lon"][1] - g["lon0"]) / s))
    r0, r1 = int((g["lat0"] - BOX["lat"][1]) / s), int(math.ceil((g["lat0"] - BOX["lat"][0]) / s))
    out = {}
    for key, code in (("depth_mm", "ssmv11036tS__"), ("swe_mm", "ssmv11034tS__")):
        o, n = mem[pick(code)]
        raw = _rows(_get(url, (o, o + n - 1)), r0, r1)
        vals = []
        for r in range(r1 - r0 + 1):
            for c in range(c0, c1 + 1):
                i = (r * g["cols"] + c) * 2
                v = int.from_bytes(raw[i:i + 2], "big", signed=True)
                vals.append(v if 0 <= v < 20000 else -1)          # depth and SWE are m x 1000 = mm; -9999 = no data
        out[key] = vals
    w, h = c1 - c0 + 1, r1 - r0 + 1
    # cell centers -> the twin's frame: a plane fit (exact to a few metres over 50 km) turns x, z into grid column, row
    pts = [(c, r, *latlon_to_xz(g["lat0"] - (r0 + r + 0.5) * s, g["lon0"] + (c0 + c + 0.5) * s)) for r in (0, h - 1) for c in (0, w - 1)]
    pts.append((w / 2 - 0.5, h / 2 - 0.5, *latlon_to_xz(g["lat0"] - (r0 + h / 2) * s, g["lon0"] + (c0 + w / 2) * s)))
    fit = lambda k: _plane([(p[2], p[3], p[k]) for p in pts])
    return {"date": day.isoformat(), "valid": f"{day.isoformat()}T06:00Z", "w": w, "h": h, "col": fit(0), "row": fit(1),
            "source": "NOAA NOHRSC SNODAS, 1 km (NSIDC G02158)", **out}


def _plane(p):
    """Least squares v = a + b*x + c*z."""
    n = len(p); sx = sum(q[0] for q in p); sz = sum(q[1] for q in p); sv = sum(q[2] for q in p)
    sxx = sum(q[0] ** 2 for q in p); szz = sum(q[1] ** 2 for q in p); sxz = sum(q[0] * q[1] for q in p)
    sxv = sum(q[0] * q[2] for q in p); szv = sum(q[1] * q[2] for q in p)
    M = [[n, sx, sz], [sx, sxx, sxz], [sz, sxz, szz]]; y = [sv, sxv, szv]
    det = lambda m: (m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
                     + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]))
    D = det(M)
    return [det([[y[i] if j == k else M[i][j] for j in range(3)] for i in range(3)]) / D for k in range(3)]


def snodas_for(day):
    """One day's grid, from the cache or NSIDC."""
    path = os.path.join(CACHE, f"snodas_{day:%Y%m%d}.json")
    if os.path.exists(path):
        return json.load(open(path))
    os.makedirs(CACHE, exist_ok=True)
    g = _snodas_day(day)
    json.dump(g, open(path, "w"), separators=(",", ":"))
    return g


def snodas():
    now = dt.datetime.now(dt.timezone.utc)
    newest = now.date() if now.hour * 60 + now.minute >= 13 * 60 + 15 else now.date() - dt.timedelta(days=1)
    cur = _state["snodas"]
    for back in range(4):                              # newest first; a day can be missing
        day = newest - dt.timedelta(days=back)
        if cur and cur["date"] >= day.isoformat():
            return cur
        try:
            _state["snodas"] = snodas_for(day)
            return _state["snodas"]
        except Exception:   # noqa: BLE001 - not posted yet or a network error: try the day before
            continue
    return cur


# ---------------------------------------------------------------------------------------------- Sentinel-2
def sentinel2():
    """Newest clear Sentinel-2 view of the twin from the last 10 days, as a PNG on the twin's 30 m grid."""
    if _state["s2"] and time.time() - _state["s2_at"] < 6 * 3600:
        return _state["s2"]
    import numpy as np
    from PIL import Image
    today = dt.datetime.now(dt.timezone.utc).date()
    body = {"collections": ["sentinel-2-l2a"], "bbox": [-119.70, 37.70, -119.18, 37.98], "limit": 40,
            "datetime": f"{today - dt.timedelta(days=10)}T00:00:00Z/{today + dt.timedelta(days=1)}T00:00:00Z",
            "sortby": [{"field": "datetime", "direction": "desc"}]}
    req = urllib.request.Request(PC_STAC, data=json.dumps(body).encode(), headers={**UA, "Content-Type": "application/json"})
    items = json.load(urllib.request.urlopen(req, timeout=60))["features"]
    by_day = {}
    for f in items:
        p = f["properties"]
        by_day.setdefault(p["datetime"][:10], {})[p.get("s2:mgrs_tile")] = (f["id"], p.get("eo:cloud_cover", 100))
    day = next((d for d in sorted(by_day, reverse=True) if "11SKB" in by_day[d] and by_day[d]["11SKB"][1] < 50), None)
    if day is None:
        raise RuntimeError("no clear Sentinel-2 pass in the last 10 days")
    png = os.path.join(LIVE, f"s2_snow_{day.replace('-', '')}.png")
    if not os.path.exists(png):
        def scl(item):
            a = np.load(io.BytesIO(_get(PC_GRID.format(id=item), timeout=180)))
            return a[0], a[1] > 0                      # classes, and where the tile has data
        cls, ok = scl(by_day[day]["11SKB"][0])
        if not ok.all() and "11SKC" in by_day[day]:    # the northern strip is in the next tile
            c2, ok2 = scl(by_day[day]["11SKC"][0])
            fill = ~ok & ok2
            cls[fill], ok = c2[fill], ok | ok2
        out = np.zeros(cls.shape, np.uint8)            # 0: not seen (cloud, shadow, no data)
        out[np.isin(cls, (4, 5, 6)) & ok] = 128        # vegetation, bare ground, water
        out[(cls == 11) & ok] = 255                    # snow and ice
        os.makedirs(LIVE, exist_ok=True)
        for old in os.listdir(LIVE):
            if old.startswith("s2_snow_"):
                os.remove(os.path.join(LIVE, old))
        Image.fromarray(out).save(png, optimize=True)
    a = np.asarray(Image.open(png))
    then = None
    try:
        then = snodas_for(dt.date.fromisoformat(day))["depth_mm"]
    except Exception:   # noqa: BLE001
        pass
    _state["s2"] = {"date": day, "item": by_day[day]["11SKB"][0], "cloud_pct": round(by_day[day]["11SKB"][1], 1),
                    "png": f"data/live/{os.path.basename(png)}", "seen_pct": round(float((a > 0).mean() * 100), 1),
                    "snow_pct": round(float((a == 255).mean() * 100), 2), "snodas_depth_mm": then,
                    "source": "Copernicus Sentinel-2 L2A scene classification (ESA), via Microsoft Planetary Computer"}
    _state["s2_at"] = time.time()
    return _state["s2"]


# ---------------------------------------------------------------------------------------------- CDEC stations
def _cdec(sensors, dur, start, end):
    ids = ",".join(s[0] for s in STATIONS)
    return json.loads(_get(CDEC.format(s=ids, n=sensors, d=dur, a=start, b=end), timeout=40))


def stations():
    now = dt.datetime.now()
    if _state["stations"] and now.timestamp() - _state["stations_at"] < 3600:
        return _state["stations"]
    by = now.year if now.month >= 10 else now.year - 1             # the last September that is over
    base_path = os.path.join(CACHE, f"cdec_baseline2_{by}.json")
    if os.path.exists(base_path):
        base = json.load(open(base_path))
    else:   # bare-ground reading of each depth (18) and water-content (3) sensor: its median over September
        rows = _cdec("18,3", "D", f"{by}-09-01", f"{by}-09-30")
        if any(not any(r["stationId"] == s[0] and r["SENSOR_NUM"] == 18 and r["value"] != -9999 for r in rows) for s in STATIONS):
            rows += _cdec("18", "H", f"{by}-09-01", f"{by}-09-30")   # some only keep hourly depth
        vals = {}
        for r in rows:
            if r["value"] not in (-9999, None) and -5 <= r["value"] < (30 if r["SENSOR_NUM"] == 18 else 5):
                vals.setdefault(f"{r['stationId']}/{r['SENSOR_NUM']}", []).append(r["value"])
        base = {k: statistics.median(v) for k, v in vals.items()}
        if base:                                                    # never keep an empty baseline
            os.makedirs(CACHE, exist_ok=True)
            json.dump(base, open(base_path, "w"))
        else:                                                       # fall back to the newest one on disk
            old = sorted(f for f in os.listdir(CACHE) if f.startswith("cdec_baseline2_")) if os.path.isdir(CACHE) else []
            base = json.load(open(os.path.join(CACHE, old[-1]))) if old else {}
    rows = _cdec("18,3", "H", (now - dt.timedelta(days=2)).strftime("%Y-%m-%d"), (now + dt.timedelta(days=1)).strftime("%Y-%m-%d"))
    last = {}
    for r in rows:   # JSON times are Pacific Standard Time; keep the newest valid reading of each sensor
        v = r["value"]
        if v in (-9999, None) or not (-5 < v < 300):
            continue
        last[(r["stationId"], r["SENSOR_NUM"])] = (r["date"], v)
    out = []
    for sid, name, lat, lon, elev in STATIONS:
        x, z = latlon_to_xz(lat, lon)
        d, w = last.get((sid, 18)), last.get((sid, 3))
        depth = None if d is None else max(0.0, d[1] - base.get(f"{sid}/18", 0)) * 2.54   # inches -> cm, minus the bare-ground reading
        swe = None if w is None else max(0.0, w[1] - base.get(f"{sid}/3", 0)) * 25.4      # inches -> mm, likewise
        if depth is not None and depth > 50 and swe is not None and swe / (depth * 10) < 0.05:
            depth = None                                                                  # far too light for snow: a sensor fault
        if depth is not None and depth < 3 and (swe is None or swe < 5):
            depth, swe = 0.0, None if swe is None else 0                                  # within the sensors' noise: bare ground
        out.append({"id": sid, "name": name, "lat": lat, "lon": lon, "elev_m": elev, "x": round(x, 1), "z": round(z, 1),
                    "depth_cm": None if depth is None else round(depth, 1), "swe_mm": None if swe is None else round(swe),
                    "time_pst": (d or w or (None,))[0]})
    _state["stations"], _state["stations_at"] = out, now.timestamp()
    return out


def snapshot():
    with _lock:
        res = {"updated": dt.datetime.now(dt.timezone.utc).isoformat(timespec="minutes")}
        for key, fn in (("snodas", snodas), ("sentinel2", sentinel2), ("stations", stations)):
            try:
                res[key] = fn()
            except Exception as e:   # noqa: BLE001 - the twin falls back to elevation rules without it
                res[key], res[key + "_error"] = None, str(e)[:200]
        return res


if __name__ == "__main__":
    s = snapshot()
    g = s["snodas"]
    if g:
        snow = [v for v in g["depth_mm"] if v > 0]
        print(g["date"], g["w"], "x", g["h"], "cells,", len(snow), "with snow, max", max(snow, default=0), "mm")
    print(json.dumps(s["stations"], indent=1)[:1500], s.get("snodas_error"), s.get("stations_error"))
