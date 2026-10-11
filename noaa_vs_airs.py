"""Compara o crescimento do metano medido pelo AIRS (NASA, satélite) com o da NOAA (rede de superfície).
Antes: descarrega a série MENSAL GLOBAL de metano da NOAA (CSV) em https://gml.noaa.gov/ccgg/trends_ch4/
       ("Globally averaged marine surface monthly mean data") e guarda como  ch4_mm_gl.csv  na mesma pasta.
Entradas: airs_ch4_400hPa.npz, ch4_mm_gl.csv       Uso: py noaa_vs_airs.py       Saída: noaa_vs_airs.json
"""
import json, re
import numpy as np
from scipy import stats

z = np.load("airs_ch4_400hPa.npz"); meses = [str(m) for m in z["meses"]]; A = z["A"].astype("float64")
lat = (89.5 - np.arange(180)) if bool(z["norte_em_cima"]) else (-89.5 + np.arange(180)); mi = np.array([int(m[5:]) - 1 for m in meses])
for m in range(12): A[mi == m] -= np.nanmean(A[mi == m], axis=0)
w = np.cos(np.radians(lat))[:, None] * np.ones((1, 360)); v = np.isfinite(A)
g = (np.where(v, A, 0) * w).sum((1, 2)) / np.maximum((v * w).sum((1, 2)), 1e-9)
airs = dict(zip(meses, g))

noaa = {}
for ln in open("ch4_mm_gl.csv", encoding="utf-8", errors="ignore"):
    if ln.startswith("#") or not ln.strip(): continue
    p = [x.strip() for x in ln.split(",")]
    try: noaa[f"{int(p[0])}-{int(p[1]):02d}"] = float(p[3])          # colunas: ano, mês, data decimal, média (ppb)
    except (ValueError, IndexError): pass
if not noaa: raise SystemExit("Não consegui ler o ficheiro da NOAA. Manda-me as 5 primeiras linhas dele.")

comum = [m for m in meses if m in noaa]
def cresc(d):                                # crescimento = diferença a 12 meses (tira a sazonalidade)
    return {m: d[m] - d[prev] for m in comum for prev in [f"{int(m[:4]) - 1}-{m[5:]}"] if prev in d}
ca, cn = cresc(airs), cresc(noaa); ks = [m for m in comum if m in ca and m in cn]
a = np.array([ca[m] for m in ks]); n = np.array([cn[m] for m in ks])
r = stats.pearsonr(a, n)[0]; sl, ic = np.polyfit(n, a, 1)
lags = {L: float(np.corrcoef(a[max(0, L):len(a) + min(0, L)], n[max(0, -L):len(n) - max(0, L)])[0, 1]) for L in range(-8, 9)}
L = max(lags, key=lags.get)
anual = {}
for m in ks: anual.setdefault(m[:4], []).append((ca[m], cn[m]))
anual = {y: [round(float(np.mean([p[0] for p in x])), 2), round(float(np.mean([p[1] for p in x])), 2)] for y, x in anual.items() if len(x) == 12}
res = dict(n_meses=len(ks), periodo=[ks[0], ks[-1]], correlacao_crescimento=round(float(r), 3), razao_airs_noaa=round(float(sl), 2),
           desfasamento_meses=int(L), correlacao_no_desfasamento=round(lags[L], 3), crescimento_anual_airs_noaa=anual)
json.dump(res, open("noaa_vs_airs.json", "w"), ensure_ascii=False)
print(f"{len(ks)} meses comparados ({ks[0]} a {ks[-1]}).")
print(f"Crescimento a 12 meses, AIRS vs NOAA: r = {r:.2f}  |  o AIRS mostra {sl:.2f}x o crescimento da NOAA  |  melhor desfasamento: {L:+d} meses (r = {lags[L]:.2f})")
print("Crescimento anual (AIRS, NOAA) ppb/ano:", {y: v for y, v in list(anual.items())[-6:]})
