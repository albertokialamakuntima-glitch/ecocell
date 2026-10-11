"""Validação independente: tendência do metano do Sentinel-5P (satélite diferente do AIRS) entre 2019-01 e 2026-08,
calculada no Earth Engine com regressão harmónica (tendência + ciclo sazonal), por célula de 1°.
Corre como TAREFAS em segundo plano (como o daily_update.py) e depois descarrega os resultados.

Uso:  py s5p_tendencia_2019.py          (demora cerca de 30 a 60 minutos; pode deixar a janela aberta)
Saída: s5p_slope_2019_2026.csv  (lat, lon, slope_ppb_ano)
Precisa de: o teu login do Earth Engine (o mesmo do export_world.py).
"""
import ee, io, time, math, urllib.request
import pandas as pd

PROJECT = "banded-lexicon-509920-b4"
ee.Initialize(project=PROJECT)
BAND = "CH4_column_volume_mixing_ratio_dry_air"
T0 = ee.Date("2019-01-01")
col = ee.ImageCollection("COPERNICUS/S5P/OFFL/L3_CH4").filterDate("2019-01-01", "2026-09-01").select(BAND)

# Médias MENSAIS (92 imagens) em vez de ~2800 imagens diárias: muito menos memória e tempo
N_MESES = 92  # 2019-01 a 2026-08
def mes_img(i):
    d0 = T0.advance(ee.Number(i), "month")
    y = col.filterDate(d0, d0.advance(1, "month")).mean().select(BAND).rename("y")
    t = ee.Image.constant(d0.advance(15, "day").difference(T0, "year")).float().rename("t")
    w = t.multiply(2 * math.pi)
    return ee.Image.cat([ee.Image.constant(1).rename("c"), t, w.sin().rename("s"), w.cos().rename("co"), y]).updateMask(y.mask())
mensais = ee.ImageCollection.fromImages(ee.List.sequence(0, N_MESES - 1).map(mes_img))
reg = mensais.reduce(ee.Reducer.linearRegression(4, 1))
slope = reg.select("coefficients").arrayProject([0]).arrayFlatten([["c", "t", "s", "co"]]).select("t").rename("slope")   # ppb por ano

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

W = 40; tasks = []
for w in range(-180, 180, W):
    for s, n in ((-60, 10), (10, 80)):
        name = f"s5p_t_{w}_{s}".replace("-", "m"); aid = f"{ROOT}/{name}"
        try: ee.data.deleteAsset(aid)
        except Exception: pass
        res = slope.reduceRegions(grid(w, w + W, s, n), ee.Reducer.mean(), scale=25000, tileScale=8).filter(ee.Filter.notNull(["mean"]))
        t = ee.batch.Export.table.toAsset(collection=res, description=name, assetId=aid); t.start(); tasks.append((aid, t)); print("enviada", name, flush=True)

vistos = set()
while True:
    st = [t.status()["state"] for _, t in tasks]; print({x: st.count(x) for x in set(st)}, flush=True)
    for aid, t in tasks:
        e = t.status()
        if e["state"] == "FAILED" and aid not in vistos: vistos.add(aid); print("  ERRO", aid.split("/")[-1], ":", str(e.get("error_message", ""))[:200], flush=True)
    if all(x in ("COMPLETED", "FAILED", "CANCELLED") for x in st): break
    time.sleep(60)

partes = []
for aid, t in tasks:
    if t.status()["state"] != "COMPLETED": print("sem sucesso:", aid, t.status().get("error_message", "")); continue
    url = ee.FeatureCollection(aid).getDownloadURL("CSV", ["lat", "lon", "mean"])
    for k in range(4):
        try: partes.append(pd.read_csv(io.StringIO(urllib.request.urlopen(url, timeout=300).read().decode()))); break
        except Exception as e: time.sleep(10 * (k + 1))
df = pd.concat(partes).rename(columns={"mean": "slope_ppb_ano"}).drop_duplicates(["lat", "lon"])
df.round(4).to_csv("s5p_slope_2019_2026.csv", index=False)
print(f"\nGravado s5p_slope_2019_2026.csv com {len(df)} células. Mediana da tendência: {df.slope_ppb_ano.median():.2f} ppb/ano")
