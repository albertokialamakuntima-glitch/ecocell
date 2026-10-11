"""2.ª variável da NASA: precipitação anual (GPM IMERG, missão da NASA/JAXA), por célula de 1°, de 2003 a 2024.
Serve para o tema do desafio: variáveis que sobem numa região e descem noutra (a chuva faz isso) e para ligar ao metano
(as zonas húmidas e a chuva influenciam o metano).
Corre como TAREFAS no Earth Engine (como o s5p_tendencia_2019.py) e grava gpm_anual_1deg.csv  (lat, lon, y2003 ... y2024, em mm/ano).
Uso:  py gpm_exportar.py          (cerca de 20 a 60 minutos)
"""
import ee, io, time, urllib.request
import pandas as pd

PROJECT = "banded-lexicon-509920-b4"
ee.Initialize(project=PROJECT)
ANOS = list(range(2003, 2025))        # anos completos; o GPM mensal V07 cobre este período
col = ee.ImageCollection("NASA/GPM_L3/IMERG_MONTHLY_V07").select("precipitation")   # mm/hora (média do mês)

def total_ano(y):
    def mm(img):                        # mm/hora x horas do mês = mm no mês
        dias = img.date().advance(1, "month").difference(img.date(), "day")
        return img.multiply(dias.multiply(24))
    return col.filterDate(f"{y}-01-01", f"{y + 1}-01-01").map(mm).sum().rename(f"y{y}")
ann = ee.Image.cat([total_ano(y) for y in ANOS])

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
        name = f"gpm_t_{w}_{s}".replace("-", "m"); aid = f"{ROOT}/{name}"
        try: ee.data.deleteAsset(aid)
        except Exception: pass
        res = ann.reduceRegions(grid(w, w + W, s, n), ee.Reducer.mean(), scale=25000, tileScale=8).filter(ee.Filter.notNull([f"y{ANOS[0]}"]))
        t = ee.batch.Export.table.toAsset(collection=res, description=name, assetId=aid); t.start(); tasks.append((aid, t)); print("enviada", name, flush=True)

vistos = set()
while True:
    st = [t.status()["state"] for _, t in tasks]; print({x: st.count(x) for x in set(st)}, flush=True)
    for aid, t in tasks:
        e = t.status()
        if e["state"] == "FAILED" and aid not in vistos: vistos.add(aid); print("  ERRO", aid.split("/")[-1], ":", str(e.get("error_message", ""))[:200], flush=True)
    if all(x in ("COMPLETED", "FAILED", "CANCELLED") for x in st): break
    time.sleep(60)

cols = ["lat", "lon"] + [f"y{y}" for y in ANOS]; partes = []
for aid, t in tasks:
    if t.status()["state"] != "COMPLETED": continue
    url = ee.FeatureCollection(aid).getDownloadURL("CSV", cols)
    for k in range(4):
        try: partes.append(pd.read_csv(io.StringIO(urllib.request.urlopen(url, timeout=300).read().decode()))); break
        except Exception: time.sleep(10 * (k + 1))
df = pd.concat(partes).drop_duplicates(["lat", "lon"])
df.round(1).to_csv("gpm_anual_1deg.csv", index=False)
print(f"\nGravado gpm_anual_1deg.csv com {len(df)} células. Precipitação média: {df[[c for c in cols[2:]]].mean().mean():.0f} mm/ano")
