"""Today's real snow for the online twin, run by a GitHub Action in the yosemite-twin repository every 6 hours
(.github/workflows/snow.yml). Same data as the laptop's GET /api/snow (00_Pipeline/snow_live.py, copied next to this file):
NOAA SNODAS depth, the newest clear Sentinel-2 view, CDEC snow sensors. Writes snow/latest.json and the Sentinel-2 PNG,
which GitHub Pages serves at https://amirkhazaeii.github.io/yosemite-twin/snow/ for the twin (window.TWIN_CONFIG.snowUrl)."""
import datetime as dt, json, os, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import snow_live  # noqa: E402

BASE_URL = "https://amirkhazaeii.github.io/yosemite-twin/snow/"
snow_live.LIVE = HERE                                  # the Sentinel-2 PNG goes next to latest.json (old ones are removed)
snow_live.CACHE = os.path.join(HERE, "cache")          # kept in the repository: the sensors' September baseline

snap = snow_live.snapshot()
s2 = snap.get("sentinel2")
if s2 and s2.get("png"):
    s2["png"] = BASE_URL + os.path.basename(s2["png"])
json.dump(snap, open(os.path.join(HERE, "latest.json"), "w"), separators=(",", ":"))
# SNODAS day caches older than 3 days are not needed again
if os.path.isdir(snow_live.CACHE):
    for f in os.listdir(snow_live.CACHE):
        p = os.path.join(snow_live.CACHE, f)
        if f.startswith("snodas_") and time.time() - os.path.getmtime(p) > 3 * 86400:
            os.remove(p)
print(dt.datetime.now(dt.timezone.utc).isoformat(timespec="minutes"),
      {k: ("ok" if snap.get(k) else snap.get(k + "_error", "none")) for k in ("snodas", "sentinel2", "stations")})
