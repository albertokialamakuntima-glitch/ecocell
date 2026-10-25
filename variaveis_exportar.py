"""Seletor de variáveis: descarrega do Google Earth Engine as MÉDIAS ANUAIS de várias variáveis das missões e modelos da NASA,
numa grelha global de 1° (360 x 180), e grava tudo em variaveis_anuais.npz.

Cada variável é pedida de uma só vez (um pedido por variável, sem tarefas): ao pedir a imagem a 1°, o Earth Engine usa as versões
reduzidas (pirâmides) dos dados, por isso é rápido mesmo para produtos de 1 km.

Uso:   py variaveis_exportar.py              (todas)
       py variaveis_exportar.py lst ndvi     (só algumas)
Precisa: earthengine-api numpy  e o teu login do Earth Engine (o mesmo dos outros scripts).
"""
import sys, io, os, json, time, urllib.request
import numpy as np
import ee

PROJECT = os.environ.get("EE_PROJECT") or "banded-lexicon-509920-b4"
if os.environ.get("EE_SERVICE_ACCOUNT_JSON"):                   # GitHub Actions: conta de serviço guardada nos segredos
    k = json.loads(os.environ["EE_SERVICE_ACCOUNT_JSON"])
    ee.Initialize(ee.ServiceAccountCredentials(k["client_email"], key_data=os.environ["EE_SERVICE_ACCOUNT_JSON"]), project=PROJECT)
else:
    ee.Initialize(project=PROJECT)
import datetime
ANOS = list(range(2003, datetime.datetime.now(datetime.timezone.utc).year))   # de 2003 até ao último ano completo (avança sozinho em janeiro)

def mod(col, banda, f, mascara=None):
    def prep(i):
        x = i.select(banda)
        if mascara: x = x.updateMask(mascara(x))
        return f(x).rename("v").copyProperties(i, ["system:time_start"])
    return ee.ImageCollection(col).map(prep)

# cada variável: coleção, conversão para unidades físicas e número mínimo de imagens por ano para o ano contar
VARS = {
 "lst":   dict(c=lambda: mod("MODIS/061/MOD11A2", "LST_Day_1km", lambda x: x.multiply(0.02).subtract(273.15)), n=20),
 "ndvi":  dict(c=lambda: mod("MODIS/061/MOD13A2", "NDVI", lambda x: x.multiply(0.0001)), n=12),
 "et":    dict(c=lambda: mod("MODIS/061/MOD16A2GF", "ET", lambda x: x.multiply(0.1 * 46), lambda x: x.lt(32700)), n=20),   # kg/m²/8 dias -> mm/ano
 "chuva": dict(c=lambda: mod("NASA/GPM_L3/IMERG_MONTHLY_V07", "precipitation", lambda x: x.multiply(8766)), n=10),   # mm/h (média do mês) -> mm/ano
 "tar":   dict(c=lambda: mod("NASA/FLDAS/NOAH01/C/GL/M/V001", "Tair_f_tavg", lambda x: x.subtract(273.15)), n=10),
 "solo":  dict(c=lambda: mod("NASA/FLDAS/NOAH01/C/GL/M/V001", "SoilMoi00_10cm_tavg", lambda x: x), n=10),
 "agua":  dict(c=lambda: mod("NASA/GRACE/MASS_GRIDS_V04/MASCON_CRI", "lwe_thickness", lambda x: x), n=6),
 "sst":   dict(c=lambda: mod("NASA/OCEANDATA/MODIS-Aqua/L3SMI", "sst", lambda x: x), n=100),
 "clor":  dict(c=lambda: mod("NASA/OCEANDATA/MODIS-Aqua/L3SMI", "chlor_a", lambda x: x.log10(), lambda x: x.gt(0)), n=100),   # log10(mg/m³)
}

def anual(col, ano, nmin):
    sub = col.filterDate(f"{ano}-01-01", f"{ano + 1}-01-01")
    m = sub.mean().rename(f"y{ano}")
    return ee.Image(ee.Algorithms.If(sub.size().gte(nmin), m, ee.Image.constant(-9999).rename(f"y{ano}").updateMask(0)))

def descarrega(img, tentativas=4):
    url = img.unmask(-9999).toFloat().getDownloadURL({"format": "NPY", "crs": "EPSG:4326",
                                                     "crs_transform": [1, 0, -180, 0, -1, 90], "dimensions": [360, 180]})
    for k in range(tentativas):
        try:
            return np.load(io.BytesIO(urllib.request.urlopen(url, timeout=600).read()))
        except Exception as e:
            print("   falhou, nova tentativa:", str(e)[:120], flush=True); time.sleep(15 * (k + 1))
    raise RuntimeError("não consegui descarregar")

escolha = sys.argv[1:] or list(VARS)
saida = "variaveis_anuais.npz"
dados = dict(np.load(saida)) if os.path.exists(saida) else {}
for v in escolha:
    t0 = time.time(); print(f"{v}: a pedir {len(ANOS)} anos ao Earth Engine...", flush=True)
    col = VARS[v]["c"]()
    img = ee.Image.cat([anual(col, a, VARS[v]["n"]) for a in ANOS])
    arr = descarrega(img)                                   # array estruturado: um campo por ano
    cubo = np.stack([arr[f"y{a}"] for a in ANOS]).astype("float32")
    cubo[cubo <= -9998] = np.nan
    dados[v] = cubo; dados["anos"] = np.array(ANOS)
    np.savez_compressed(saida, **dados)
    ok = np.isfinite(cubo).any(0).sum()
    print(f"   ok em {time.time() - t0:.0f} s: {ok} células com dados | média {np.nanmean(cubo):.3f}", flush=True)
print("Gravado", saida)
