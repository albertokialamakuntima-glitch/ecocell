"""Compara a tendência do AIRS (NASA) com a do Sentinel-5P (ESA) no MESMO período (2019-01 a 2026-08), célula a célula.
Se dois instrumentos independentes concordam em quais regiões sobem mais ou menos, o resultado ganha credibilidade.

Entradas (na mesma pasta): airs_ch4_400hPa.npz, s5p_slope_2019_2026.csv, world.json
Uso: py comparar_s5p_airs.py        Saídas: validacao.json, validacao.png      Precisa: numpy pandas scipy matplotlib shapely
"""
import json
import numpy as np, pandas as pd
from scipy import stats
from shapely.geometry import shape
from shapely import contains_xy
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

z = np.load("airs_ch4_400hPa.npz"); meses = [str(m) for m in z["meses"]]; A = z["A"].astype("float64")
lat_a = (89.5 - np.arange(180)) if bool(z["norte_em_cima"]) else (-89.5 + np.arange(180)); lon_a = -179.5 + np.arange(360)
sel = [i for i, m in enumerate(meses) if "2019-01" <= m <= "2026-08"]
A = A[sel]; ano = np.array([int(meses[i][:4]) + (int(meses[i][5:]) - 1) / 12 for i in sel]); mes = np.array([int(meses[i][5:]) - 1 for i in sel])
for m in range(12): A[mes == m] -= np.nanmean(A[mes == m], axis=0)                      # tira o ciclo sazonal
tc = ano - ano.mean(); ok = np.isfinite(A)
num = np.nansum(np.where(ok, A * tc[:, None, None], 0), 0); den = np.maximum((ok * tc[:, None, None] ** 2).sum(0), 1e-9)
slope_a = np.where(ok.sum(0) > 0.7 * len(sel), num / den, np.nan)                      # ppb por ano, em cada célula do AIRS

s5 = pd.read_csv("s5p_slope_2019_2026.csv")
ia = {round(float(v), 1): k for k, v in enumerate(lat_a)}; ja = {round(float(v), 1): k for k, v in enumerate(lon_a)}
def airs_em(la, lo):                                                                     # a grelha do S5P tem latitudes inteiras: média das 2 linhas vizinhas do AIRS
    j = ja.get(round(lo, 1));  v = [slope_a[ia[round(la + d, 1)], j] for d in (-0.5, 0.5) if round(la + d, 1) in ia] if j is not None else []
    return np.nanmean(v) if v and np.isfinite(v).any() else np.nan
s5["airs"] = [airs_em(a, b) for a, b in zip(s5.lat, s5.lon)]
s5 = s5.dropna(); s5 = s5[(s5.slope_ppb_ano.abs() < 60)]
w = np.cos(np.radians(s5.lat.values)); ga = np.average(s5.airs, weights=w); gs = np.average(s5.slope_ppb_ano, weights=w)
s5["rel_airs"] = s5.airs - ga; s5["rel_s5p"] = s5.slope_ppb_ano - gs
r = stats.pearsonr(s5.rel_airs, s5.rel_s5p)[0]; rho = stats.spearmanr(s5.rel_airs, s5.rel_s5p)[0]
conc = float((np.sign(s5.rel_airs) == np.sign(s5.rel_s5p)).mean())
print(f"{len(s5)} células comparadas. Tendência global: AIRS {ga:+.2f} | Sentinel-5P {gs:+.2f} ppb/ano")
print(f"Diferenças regionais (face à média global): Pearson r = {r:.2f} | Spearman = {rho:.2f} | mesmo sinal em {conc:.0%} das células")

paises = {}
LON, LAT = np.meshgrid(s5.lon.unique(), s5.lat.unique())
for f in json.load(open("world.json", encoding="utf-8"))["features"]:
    nome = f["properties"]["name"]
    if nome == "Antarctica": continue
    m = contains_xy(shape(f["geometry"]), s5.lon.values, s5.lat.values)
    if m.sum() >= 5: paises[nome] = dict(airs=round(float(np.average(s5.rel_airs[m], weights=w[m])), 2), s5p=round(float(np.average(s5.rel_s5p[m], weights=w[m])), 2), n=int(m.sum()))
pa = np.array([[v["airs"], v["s5p"]] for v in paises.values()]); rp = stats.pearsonr(pa[:, 0], pa[:, 1])[0]
print(f"Por país ({len(paises)}): r = {rp:.2f}")
json.dump(dict(periodo=["2019-01", "2026-08"], n_celulas=len(s5), global_airs=round(ga, 2), global_s5p=round(gs, 2), pearson=round(r, 3), spearman=round(rho, 3),
               concordancia_sinal=round(conc, 3), pearson_paises=round(rp, 3), paises=paises), open("validacao.json", "w", encoding="utf-8"), ensure_ascii=False)
fig, ax = plt.subplots(1, 2, figsize=(10, 4.2)); sm = s5.sample(min(len(s5), 4000), random_state=0)
ax[0].scatter(sm.rel_airs, sm.rel_s5p, s=4, alpha=.4); ax[0].axhline(0, c="#999", lw=.6); ax[0].axvline(0, c="#999", lw=.6)
ax[0].set(xlabel="AIRS (NASA): tendência face à global (ppb/ano)", ylabel="Sentinel-5P (ESA)", title=f"Células: r = {r:.2f}")
ax[1].scatter(pa[:, 0], pa[:, 1], s=14, c="#1f7a4d"); ax[1].set(xlabel="AIRS, por país", ylabel="Sentinel-5P, por país", title=f"Países: r = {rp:.2f}")
fig.tight_layout(); fig.savefig("validacao.png", dpi=130); print("Gravado validacao.json e validacao.png")
