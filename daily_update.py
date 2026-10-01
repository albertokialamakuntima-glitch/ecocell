"""Atualiza cells.csv com o CH4 mais recente do Sentinel-5P (Google Earth Engine).
O cálculo pesado corre como TAREFAS em segundo plano (exportação para assets do Earth Engine),
sem limite de 5 minutos. Depois só se descarregam os resultados já calculados.
Credenciais: EE_SERVICE_ACCOUNT_JSON e EE_PROJECT (GitHub Actions) ou o teu login local."""
import ee, os, io, json, time, urllib.request
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

ROOT = f"projects/{PROJECT}/assets/ecocell"
try: ee.data.getAsset(ROOT)
except Exception: ee.data.createAsset({"type": "FOLDER"}, ROOT)

def grid(w, e, s, n):
    def c_lon(x):
        x = ee.Number(x)
        return ee.List.sequence(s, n - 1).map(lambda la: ee.Feature(
            ee.Geometry.Rectangle([x, ee.Number(la).subtract(.5), x.add(1), ee.Number(la).add(.5)], None, False),
            {"lon": x.add(.5), "lat": la}))
    return ee.FeatureCollection(ee.List.sequence(w, e - 1).map(c_lon).flatten())

W = 40
jobs = [(w, s, min(s + 70, 80)) for w in range(-180, 180, W) for s in (-60, 10)]
tasks = []
for w, s, n in jobs:
    name = f"t_{w}_{s}".replace("-", "m"); aid = f"{ROOT}/{name}"
    try: ee.data.deleteAsset(aid)
    except Exception: pass
    res = stack.reduceRegions(grid(w, w + W, s, n), ee.Reducer.mean(), scale=7000, tileScale=4).filter(ee.Filter.notNull(["mean", "slope"]))
    t = ee.batch.Export.table.toAsset(collection=res, description=name, assetId=aid)
    t.start(); tasks.append((w, s, n, aid, t)); print("enviada", name, flush=True)

deadline = time.time() + 130 * 60
while time.time() < deadline:
    st = [t.status()["state"] for *_, t in tasks]
    print({x: st.count(x) for x in set(st)}, flush=True)
    if all(x in ("COMPLETED", "FAILED", "CANCELLED") for x in st): break
    time.sleep(60)

def read(aid):
    url = ee.FeatureCollection(aid).getDownloadURL("CSV", ["lat", "lon", "mean", "slope", "recent"])
    for n in range(4):
        try: return pd.read_csv(io.StringIO(urllib.request.urlopen(url, timeout=300).read().decode()))
        except Exception as e: time.sleep(10 * (n + 1)); err = e
    print("falhou a leitura", aid, err); return None

df = pd.read_csv("cells.csv") if os.path.exists("cells.csv") else pd.DataFrame(columns=["lat", "lon", "mean", "slope", "recent"])
ok = 0
for w, s, n, aid, t in tasks:
    s_ = t.status()
    if s_["state"] != "COMPLETED": print("tarefa sem sucesso", aid, s_["state"], s_.get("error_message", "")); continue
    part = read(aid)
    if part is None or part.empty: continue
    m = (df.lon >= w) & (df.lon < w + W) & (df.lat >= s) & (df.lat < n)
    df = pd.concat([df[~m], part[["lat", "lon", "mean", "slope", "recent"]]]); ok += 1
print(f"{ok}/{len(tasks)} blocos atualizados")
if ok == 0 or len(df) < 20000: raise SystemExit("nada de novo para guardar, a manter os dados antigos")
df.drop_duplicates(["lat", "lon"]).round(5).to_csv("cells.csv", index=False)
json.dump({"start": str(START), "end": str(END), "blocos_ok": f"{ok}/{len(tasks)}"}, open("meta.json", "w"))
print("ok", len(df), "células", START, END)
