"""Camada de células do metano SÓ com dados da NASA (substitui o Sentinel-5P da ESA).
A partir do metano mensal do NASA Aqua AIRS (400 hPa, airs_ch4_400hPa.npz do airs_baixar.py) calcula, por célula de 1°:
  mean   = média dos últimos 3 anos completos (36 meses)
  slope  = tendência em ppb/dia (regressão harmónica: constante + tendência + seno + cosseno anuais, para tirar a estação)
  recent = média + anomalia dos 2 últimos meses face ao esperado (tendência + estação)
Grava cells.csv e meta.json no formato que o build_site.py já usa.
Uso: py airs_celulas.py            Precisa: numpy
"""
import json
import numpy as np

z = np.load("airs_ch4_400hPa.npz"); ms = [str(m) for m in z["meses"]]; A = z["A"].astype("float64")
if not bool(z["norte_em_cima"]): A = A[:, ::-1, :]
fim = len(ms) - 1
while not ms[fim].endswith("-12") and fim > 0 and np.isfinite(A[fim]).mean() < .2: fim -= 1   # último mês com dados
ini = fim - 35; sel = list(range(ini, fim + 1)); n = len(sel)
Y = A[sel].reshape(n, -1)
t = np.arange(n) / 12.0; w = 2 * np.pi * t
X = np.column_stack([np.ones(n), t, np.sin(w), np.cos(w)])
ok = np.isfinite(Y); val = ok.sum(0) >= 24
coef = np.full((4, Y.shape[1]), np.nan)
for j in np.where(val)[0]:                                # mínimos quadrados só com os meses válidos de cada célula
    m = ok[:, j]; coef[:, j] = np.linalg.lstsq(X[m], Y[m, j], rcond=None)[0]
media = np.nanmean(Y, 0)
ajuste = X @ np.nan_to_num(coef)
res = np.where(ok, Y - ajuste, np.nan)
recente = media + np.nanmean(res[-2:], 0)
lat = 89.5 - np.arange(180); lon = -179.5 + np.arange(360); LON, LAT = np.meshgrid(lon, lat)
linhas = ["lat,lon,mean,slope,recent"]
for j in np.where(val & np.isfinite(media))[0]:
    rr = recente[j] if np.isfinite(recente[j]) else media[j]
    linhas.append(f"{LAT.ravel()[j]:.1f},{LON.ravel()[j]:.1f},{media[j]:.3f},{coef[1, j] / 365.25:.6f},{rr:.3f}")
open("cells.csv", "w", encoding="utf-8").write("\n".join(linhas) + "\n")
a0, a1 = ms[ini], ms[fim]
json.dump({"start": f"{a0}-01", "end": f"{a1}-28", "fonte": "NASA Aqua AIRS L3 v7.0, CH4 400 hPa (diurno)", "metodo": "harmonica_36meses"}, open("meta.json", "w"))
print(f"cells.csv: {len(linhas) - 1} células | janela {a0} a {a1} | média {np.nanmean(media):.1f} ppb | tendência média {np.nanmedian(coef[1]):+.2f} ppb/ano")
