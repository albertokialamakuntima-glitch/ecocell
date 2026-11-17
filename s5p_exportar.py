"""Metano do Sentinel-5P (ESA Copernicus, TROPOMI) pelo Google Earth Engine: médias mensais em blocos de 1°, desde 2019,
até ao último mês que o Earth Engine já tem. Corre todos os dias (daily_update.py) e grava os «anos móveis» (anos_moveis.py):
  ch4s5p_anual_f32.bin + ch4s5p_anual.json
Coleção: COPERNICUS/S5P/OFFL/L3_CH4, banda CH4_column_volume_mixing_ratio_dry_air_bias_corrected (ppb, coluna total).
Contém dados Copernicus Sentinel modificados.

Uso: py s5p_exportar.py        Precisa: earthengine-api numpy e o login do Earth Engine (o mesmo do variaveis_exportar.py)
"""
import io, os, json, time, datetime, urllib.request
import numpy as np
import ee
from anos_moveis import grava

PROJECT = os.environ.get("EE_PROJECT") or "banded-lexicon-509920-b4"
if os.environ.get("EE_SERVICE_ACCOUNT_JSON"):
    k = json.loads(os.environ["EE_SERVICE_ACCOUNT_JSON"])
    ee.Initialize(ee.ServiceAccountCredentials(k["client_email"], key_data=os.environ["EE_SERVICE_ACCOUNT_JSON"]), project=PROJECT)
else:
    ee.Initialize(project=PROJECT)

BANDA = "CH4_column_volume_mixing_ratio_dry_air_bias_corrected"
col = ee.ImageCollection("COPERNICUS/S5P/OFFL/L3_CH4").select(BANDA)
hoje = datetime.date.today()


def mensal(a, m, nmin=4):
    ini = ee.Date.fromYMD(a, m, 1); sub = col.filterDate(ini, ini.advance(1, "month"))
    nome = f"m{a}{m:02d}"
    return ee.Image(ee.Algorithms.If(sub.size().gte(nmin), sub.mean().rename(nome), ee.Image.constant(-9999).rename(nome).updateMask(0)))


def descarrega(img, tentativas=4):
    url = img.unmask(-9999).toFloat().getDownloadURL({"format": "NPY", "crs": "EPSG:4326",
                                                     "crs_transform": [1, 0, -180, 0, -1, 90], "dimensions": [360, 180]})
    for k in range(tentativas):
        try: return np.load(io.BytesIO(urllib.request.urlopen(url, timeout=600).read()))
        except Exception as e: print("   falhou, nova tentativa:", str(e)[:120], flush=True); time.sleep(15 * (k + 1))
    raise RuntimeError("não consegui descarregar")


meses, cubo = [], []
for a in range(2019, hoje.year + 1):                       # um pedido por ano (12 bandas mensais)
    ms = [(a, m) for m in range(1, 13) if (a, m) >= (2019, 2) and (a, m) <= (hoje.year, hoje.month)]
    t0 = time.time(); arr = descarrega(ee.Image.cat([mensal(x, y) for x, y in ms]))
    for x, y in ms:
        v = arr[f"m{x}{y:02d}"].astype("float32"); v[(v <= -9998) | (v < 1000) | (v > 3000)] = np.nan
        meses.append(f"{x}-{y:02d}"); cubo.append(v)
    print(f"{a}: {len(ms)} meses em {time.time() - t0:.0f} s", flush=True)
cubo = np.array(cubo)
ok = [i for i in range(len(meses)) if np.isfinite(cubo[i]).any()]
if not ok: raise SystemExit("Sentinel-5P: sem dados")
cubo, meses = cubo[: ok[-1] + 1], meses[: ok[-1] + 1]
grava("ch4s5p", meses, cubo, min_meses=6, fonte="ESA Copernicus Sentinel-5P TROPOMI L3 OFFL (Google Earth Engine)")
