"""Atualiza cells.csv com o CH4 mais recente do Sentinel-5P (Google Earth Engine).
Sem getInfo gigante: a grelha é dividida em blocos pequenos e cada um é descarregado por URL.
Credenciais: variáveis EE_SERVICE_ACCOUNT_JSON e EE_PROJECT (GitHub Actions) ou o teu login local."""
import ee, os, io, json, time, urllib.request, concurrent.futures as cf
import pandas as pd
from datetime import date, timedelta

PROJECT = os.environ.get("EE_PROJECT", "banded-lexicon-509920-b4")
if os.environ.get("EE_SERVICE_ACCOUNT_JSON"):
    k = json.loads(os.environ["EE_SERVICE_ACCOUNT_JSON"])
    ee.Initialize(ee.ServiceAccountCredentials(k["client_email"], key_data=os.environ["EE_SERVICE_ACCOUNT_JSON"]), project=PROJECT)
else:
    ee.Initialize(project=PROJECT)

END = date.today(); START = END - timedelta(days=540); RECENT = END - timedelta(days=60)
BAND = "CH4_column_volume_mixing_ratio_dry_air"
col = ee.ImageCollection("COPERNICUS/S5P/OFFL/L3_CH4").select(BAND)
base = col.filterDate(str(START), str(END))
t0 = ee.Date(str(START))
fit = base.map(lambda i: i.addBands(ee.Image(i.date().difference(t0, "day")).float().rename("t"))).select(["t", BAND]).reduce(ee.Reducer.linearFit())
stack = (base.mean().rename("mean").addBands(fit.select("scale").rename("slope"))
         .addBands(col.filterDate(str(RECENT), str(END)).mean().rename("recent")))

def tile(w, s, size=20):
    def c_lon(x):
        x = ee.Number(x)
        return ee.List.sequence(s, s + size - 1).map(lambda la: ee.Feature(
            ee.Geometry.Rectangle([x, ee.Number(la).subtract(.5), x.add(1), ee.Number(la).add(.5)], None, False),
            {"lon": x.add(.5), "lat": la}))
    fc = ee.FeatureCollection(ee.List.sequence(w, w + size - 1).map(c_lon).flatten())
    res = stack.reduceRegions(fc, ee.Reducer.mean(), scale=7000, tileScale=4).filter(ee.Filter.notNull(["mean", "slope"]))
    url = res.getDownloadURL("CSV", ["lat", "lon", "mean", "slope", "recent"])
    for n in range(4):
        try: return pd.read_csv(io.StringIO(urllib.request.urlopen(url, timeout=300).read().decode()))
        except Exception as e: time.sleep(10 * (n + 1)); err = e
    print("falhou", w, s, err); return pd.DataFrame()

jobs = [(w, s) for w in range(-180, 180, 20) for s in range(-60, 80, 20)]
with cf.ThreadPoolExecutor(6) as ex: parts = list(ex.map(lambda a: tile(*a), jobs))
df = pd.concat(parts).drop_duplicates(["lat", "lon"])
if len(df) < 20000: raise SystemExit(f"só {len(df)} células, a manter os dados antigos")
df.round(5).to_csv("cells.csv", index=False)
json.dump({"start": str(START), "end": str(END)}, open("meta.json", "w"))
print("ok", len(df), "células", START, END)
