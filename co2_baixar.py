"""Dióxido de carbono (CO2) da NASA: OCO-2 assimilado no modelo GEOS, Level 3 mensal, v10r (OCO2_GEOS_L3CO2_MONTH).
Descarrega os meses de 2015 a 2021 (~3 MB cada), passa a grelha de 0,5° x 0,625° para blocos de 1° e grava as médias anuais
(só anos com os 12 meses) em:
  co2_anual_f32.bin   float32, (anos, 180, 360), Norte em cima, em ppm (coluna total, XCO2)
  co2_anual.json      anos e descrição
O produto da NASA termina em fevereiro de 2022, por isso basta correr uma vez (o daily_update.py só o corre se faltar o resultado).

Uso (com o token definido): py co2_baixar.py        Precisa: numpy earthaccess netCDF4
"""
import os, re, glob, json, time
import numpy as np, earthaccess, netCDF4

ANOS = list(range(2015, 2022))
CACHE, TMP = "airs_tmp/co2_cache", "airs_tmp/co2_dl"                # dentro de airs_tmp/, que a atualização automática nunca envia para o GitHub
os.makedirs(CACHE, exist_ok=True); os.makedirs(TMP, exist_ok=True)


def para_1grau(ds):
    """XCO2 (mol/mol) numa grelha de 0,5° x 0,625° -> média em blocos de 1° (ppm), Norte em cima"""
    nome = next((v for v in ds.variables if v.upper() == "XCO2"), None) or next((v for v in ds.variables if "XCO2" in v.upper()), None)
    if not nome: raise KeyError("sem XCO2; variáveis: " + ", ".join(ds.variables))
    la = np.asarray(ds.variables[next(v for v in ds.variables if v.lower() in ("lat", "latitude"))][:], dtype="float64")
    lo = np.asarray(ds.variables[next(v for v in ds.variables if v.lower() in ("lon", "longitude"))][:], dtype="float64")
    x = np.ma.filled(ds.variables[nome][:].astype("float64"), np.nan).squeeze()
    if x.shape != (la.size, lo.size): x = x.reshape(la.size, lo.size)
    if np.nanmedian(x) < 1: x = x * 1e6                                   # mol/mol -> ppm
    x[(x < 300) | (x > 500)] = np.nan
    lo = np.where(lo >= 180, lo - 360, lo)
    LO, LA = np.meshgrid(lo, la)
    r = np.clip(np.floor(90 - LA), 0, 179).astype(int); c = np.clip(np.floor(LO + 180), 0, 359).astype(int)
    ok = np.isfinite(x); S = np.zeros((180, 360)); N = np.zeros((180, 360))
    np.add.at(S, (r[ok], c[ok]), x[ok]); np.add.at(N, (r[ok], c[ok]), 1)
    return np.where(N > 0, S / np.maximum(N, 1), np.nan).astype("float32")


if __name__ == "__main__":
    if os.environ.get("EARTHDATA_TOKEN"): earthaccess.login(strategy="environment")
    else: earthaccess.login(strategy="interactive", persist=True)
    res = earthaccess.search_data(short_name="OCO2_GEOS_L3CO2_MONTH", version="10r", temporal=(f"{ANOS[0]}-01-01", f"{ANOS[-1]}-12-31"), count=-1)
    print("ficheiros mensais encontrados:", len(res))
    for i, g in enumerate(res, 1):
        m = re.search(r"month_(\d{4})(\d{2})", g.data_links()[0])
        if not m: continue
        out = f"{CACHE}/{m.group(1)}-{m.group(2)}.npy"
        if os.path.exists(out): continue
        for tentativa in range(5):
            try:
                f = str(earthaccess.download([g], TMP)[0])
                with netCDF4.Dataset(f) as ds: np.save(out, para_1grau(ds))
                os.remove(f); print(f"[{i}/{len(res)}] {m.group(1)}-{m.group(2)} ok", flush=True); break
            except Exception as e:
                print(f"  falhou {m.group(1)}-{m.group(2)} ({tentativa + 1}/5):", str(e)[:120], flush=True); time.sleep(15 * (tentativa + 1))
    cub, anos = [], []
    for y in ANOS:
        fs = sorted(glob.glob(f"{CACHE}/{y}-*.npy"))
        if len(fs) < 12: print(f"{y}: só {len(fs)} meses, fica de fora"); continue
        cub.append(np.nanmean([np.load(f) for f in fs], 0)); anos.append(y)
    if len(anos) < 5: raise SystemExit("CO2: anos completos insuficientes")
    np.asarray(cub, dtype="<f4").tofile("co2_anual_f32.bin")
    json.dump(dict(anos=anos, forma=[len(anos), 180, 360], unidade="ppm", produto="OCO2_GEOS_L3CO2_MONTH v10r", variavel="XCO2"), open("co2_anual.json", "w", encoding="utf-8"))
    print(f"co2_anual_f32.bin: {anos[0]}-{anos[-1]} | média global {np.nanmean(cub[0]):.1f} ppm em {anos[0]} -> {np.nanmean(cub[-1]):.1f} ppm em {anos[-1]}")
