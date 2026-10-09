"""Tendências de metano (NASA AIRS) para o mundo inteiro: por célula de 1° e por país.
Responde às 4 perguntas do desafio: o que muda (direção), onde (célula/país), quanto (ppb/ano, IC 95%)
e se é significativo (Mann-Kendall corrigido para autocorrelação + correção FDR de Benjamini-Hochberg).

Entrada : airs_ch4_400hPa.npz (do airs_baixar.py), world.json (do teu repositório) e continentes.json
Saída   : airs_trends.json
Uso     : py airs_tendencias_mundo.py            (completo, uns minutos)
          py airs_tendencias_mundo.py 3          (teste rápido: 1 célula em cada 3)
Precisa : pip install numpy scipy shapely
"""
import sys, json
import numpy as np
from scipy import stats
from shapely.geometry import shape
from shapely import contains_xy
from shapely.ops import unary_union
from shapely.validation import make_valid

PASSO = int(sys.argv[1]) if len(sys.argv) > 1 else 1
z = np.load("airs_ch4_400hPa.npz")
meses = [str(m) for m in z["meses"]]; A = z["A"].astype("float64")          # (tempo, 180, 360), diurno
nt = len(meses); ano = np.array([int(m[:4]) + (int(m[5:]) - 1) / 12 for m in meses]); mes = np.array([int(m[5:]) - 1 for m in meses])
lat = (89.5 - np.arange(180)) if bool(z["norte_em_cima"]) else (-89.5 + np.arange(180)); lon = -179.5 + np.arange(360)

def sem_sazonal(x):                      # x: (tempo, ...) -> tira a média de cada mês do ano
    out = x.copy()
    for m in range(12):
        s = mes == m
        if s.any(): out[s] = x[s] - np.nanmean(x[s], axis=0)
    return out

def tendencia(t, y):
    """Sen + Mann-Kendall com variância corrigida pela autocorrelação (n efetivo)."""
    ok = ~np.isnan(y); t, y = t[ok], y[ok]; n = len(y)
    if n < 0.6 * nt or n < 36: return None
    sl, ic, lo, hi = stats.theilslopes(y, t, 0.95)
    tau = stats.kendalltau(t, y)[0]; S = tau * n * (n - 1) / 2; var = n * (n - 1) * (2 * n + 5) / 18
    r = y - (ic + sl * t); r1 = float(np.corrcoef(r[:-1], r[1:])[0, 1]); r1 = min(max(r1, 0.0), 0.9)
    var *= (1 + r1) / (1 - r1)                                              # inflaciona a variância se houver memória
    zz = (S - np.sign(S)) / np.sqrt(var); p = 2 * (1 - stats.norm.cdf(abs(zz)))
    return sl, lo, hi, p, n, r1

def fdr(p):                               # Benjamini-Hochberg
    p = np.asarray(p); o = np.argsort(p); q = np.empty_like(p); m = len(p)
    q[o] = np.minimum.accumulate((p[o] * m / (np.arange(m) + 1))[::-1])[::-1]
    return np.minimum(q, 1)

anom = sem_sazonal(A)
LON, LAT = np.meshgrid(lon, lat); W = np.cos(np.radians(LAT))
FEATS = [f for f in json.load(open("world.json", encoding="utf-8"))["features"] if f["properties"]["name"] != "Antarctica"]
terra = unary_union([make_valid(shape(f["geometry"]).simplify(0.2)).buffer(2) for f in FEATS])
mascara = contains_xy(terra, LON, LAT)
print(f"{nt} meses ({meses[0]} a {meses[-1]}). Células de terra e costa: {int(mascara.sum())}. A calcular (passo {PASSO})...")

# tendência global do planeta (média ponderada de TODAS as células)
v_ = ~np.isnan(anom); gser = np.nansum(np.where(v_, anom, 0) * W, axis=(1, 2)) / np.maximum((v_ * W).sum(axis=(1, 2)), 1e-9)
rg = tendencia(ano, gser); glob = dict(slope=rg[0], lo=rg[1], hi=rg[2], p=rg[3], r1=rg[5]); print(f"Tendência global: {rg[0]:+.2f} ppb/ano")
anos = sorted(set(int(a) for a in ano)); iy = [(ano >= a) & (ano < a + 1) for a in anos]

cel = []
for i in range(0, 180, PASSO):
    for j in range(0, 360, PASSO):
        if not mascara[i, j]: continue
        r = tendencia(ano, anom[:, i, j])
        if r: cel.append((i, j) + r)
    if i % 20 == 0: print("  linha", i, "de 180", flush=True)
p_c = np.array([c[5] for c in cel]); q_c = fdr(p_c)
cells = []
for c, q in zip(cel, q_c):
    ser = [None if not np.isfinite(np.nanmean(anom[m, c[0], c[1]])) else int(round(float(np.nanmean(anom[m, c[0], c[1]])))) for m in iy]
    cells.append([round(float(lat[c[0]]), 1), round(float(lon[c[1]]), 1), round(c[2], 3), round(c[3], 3), round(c[4], 3), round(float(c[5]), 4),
                  round(float(q), 4), int(q < 0.05 and (c[3] > 0 or c[4] < 0)), round(c[7], 2), ser])
print(f"{len(cells)} células analisadas; {sum(c[7] for c in cells)} com tendência significativa (q<0,05).")

# ---- países (média ponderada por cos(latitude) das anomalias das células dentro do país)
import os
CONT = json.load(open("continentes.json", encoding="utf-8")) if os.path.exists("continentes.json") else {}
mk = {}
def resumir(m):
    w = W[m]; X = anom[:, m]; v = ~np.isnan(X)
    ser = np.where(v.sum(1) > 0, (np.nansum(np.where(v, X, 0) * w, axis=1) / np.maximum((v * w).sum(1), 1e-9)), np.nan)
    r = tendencia(ano, ser)
    if not r: return None
    base = float(np.nanmean(A[:, m])); med = [float(np.nanmean(ser[m_])) if m_.any() else None for m_ in iy]
    return dict(slope=r[0], lo=r[1], hi=r[2], p=r[3], n=r[4], r1=r[5], celulas=int(m.sum()), media=base,
                pct_decada=r[0] * 10 / base * 100, anos=anos, serie_anual=[None if x is None else round(x, 2) for x in med])
paises = {}
for f in FEATS:
    nome = f["properties"]["name"]
    if nome == "Antarctica": continue
    m = contains_xy(shape(f["geometry"]), LON, LAT)
    if nome in CONT: mk[CONT[nome]] = mk.get(CONT[nome], False) | m
    if m.sum() < 3: continue
    w = W[m]; X = anom[:, m]; v = ~np.isnan(X)
    ser = np.where(v.sum(1) > 0, (np.nansum(np.where(v, X, 0) * w, axis=1) / np.maximum((v * w).sum(1), 1e-9)), np.nan)
    r = tendencia(ano, ser)
    if not r: continue
    base = float(np.nanmean(A[:, m]))
    med = [float(np.nanmean(ser[m_])) if m_.any() else None for m_ in iy]
    paises[nome] = dict(slope=r[0], lo=r[1], hi=r[2], p=r[3], n=r[4], r1=r[5], celulas=int(m.sum()), media=base,
                        pct_decada=r[0] * 10 / base * 100, anos=anos, serie_anual=[None if x is None else round(x, 2) for x in med])
nomes = list(paises); qp = fdr([paises[k]["p"] for k in nomes])
for k, q in zip(nomes, qp):
    d = paises[k]; d["q"] = float(q); d["significativo"] = int(q < 0.05 and (d["lo"] > 0 or d["hi"] < 0))
    for c in ("slope", "lo", "hi", "pct_decada", "media"): d[c] = round(d[c], 3)
    d["p"] = round(d["p"], 5); d["q"] = round(d["q"], 5); d["r1"] = round(d["r1"], 2)

cont_out = {k_: d_ for k_, m_ in mk.items() for d_ in [resumir(m_)] if d_}
for (k_, d_), q_ in zip(cont_out.items(), fdr([d_["p"] for d_ in cont_out.values()]) if cont_out else []):
    d_["q"] = float(q_); d_["significativo"] = int(q_ < 0.05 and (d_["lo"] > 0 or d_["hi"] < 0))
    for c_ in ("slope", "lo", "hi", "pct_decada", "media"): d_[c_] = round(d_[c_], 3)
    d_["p"] = round(d_["p"], 5); d_["q"] = round(d_["q"], 5); d_["r1"] = round(d_["r1"], 2)
serie_global = [round(float(np.nanmean(gser[m_])), 2) for m_ in iy]
print(f"{len(cont_out)} continentes calculados.")
json.dump(dict(periodo=[meses[0], meses[-1]], nivel_hPa=400, instrumento="NASA Aqua AIRS L3 mensal v7.0 (diurno)",
               metodo="Sen + Mann-Kendall (variância corrigida por autocorrelação lag-1), FDR Benjamini-Hochberg q<0.05",
               anos_celulas=anos, globais=glob, serie_global=serie_global, continentes=cont_out, celulas=cells, paises=paises), open("airs_trends.json", "w", encoding="utf-8"), ensure_ascii=False)
sig = [k for k in nomes if paises[k]["significativo"]]
print(f"\n{len(paises)} países; {len(sig)} com tendência significativa. Gravado airs_trends.json")
for k in sorted(nomes, key=lambda k: -abs(paises[k]["slope"]))[:8]:
    d = paises[k]; print(f"  {k:28s} {d['slope']:+.2f} ppb/ano [{d['lo']:+.2f}, {d['hi']:+.2f}]  q={d['q']:.3g}  {'SIGNIFICATIVO' if d['significativo'] else 'não significativo'}")
