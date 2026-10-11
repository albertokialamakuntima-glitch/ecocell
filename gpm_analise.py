"""Tendência da precipitação anual (NASA GPM IMERG), 2003-2024: por célula de 1° e por país. Sen + Mann-Kendall + FDR."""
import json, numpy as np, pandas as pd
from scipy import stats
from shapely.geometry import shape
from shapely import contains_xy
d = pd.read_csv("gpm_anual_1deg.csv"); cols = [c for c in d.columns if c.startswith("y")]; anos = np.array([int(c[1:]) for c in cols])
Y = d[cols].values.astype(float); print(len(d), "células;", len(anos), "anos;", anos[0], "a", anos[-1])
def tend(t, y):
    ok = np.isfinite(y)
    if ok.sum() < 18: return None
    t, y = t[ok], y[ok]; sl, ic, lo, hi = stats.theilslopes(y, t, 0.95); p = stats.kendalltau(t, y)[1]
    return sl, lo, hi, float(p), float(y.mean())
def fdr(p):
    p = np.asarray(p); o = np.argsort(p); q = np.empty_like(p); m = len(p); q[o] = np.minimum.accumulate((p[o] * m / (np.arange(m) + 1))[::-1])[::-1]; return np.minimum(q, 1)
R = [tend(anos, y) for y in Y]; idx = [i for i, r in enumerate(R) if r and r[4] > 50]            # ignora células quase sem chuva (<50 mm/ano)
sl = np.array([R[i][0] for i in idx]); p = np.array([R[i][3] for i in idx]); q = fdr(p); mean = np.array([R[i][4] for i in idx]); lat = d.lat.values[idx]; lon = d.lon.values[idx]
pct = sl * 10 / mean * 100; sig = q < 0.05
print(f"{len(idx)} células com chuva. Significativas (q<0,05): {sig.sum()} ({sig.mean():.0%}) | mais húmidas: {(sig & (sl > 0)).sum()} | mais secas: {(sig & (sl < 0)).sum()}")
print(f"Mediana da tendência: {np.median(pct):+.2f} %/década | células com sinal positivo {np.mean(sl > 0):.0%}, negativo {np.mean(sl < 0):.0%}")
w = np.cos(np.radians(d.lat.values)); LON, LAT = d.lon.values, d.lat.values
paises = {}
for f in json.load(open("world.json", encoding="utf-8"))["features"]:
    nome = f["properties"]["name"]
    if nome == "Antarctica": continue
    m = contains_xy(shape(f["geometry"]), LON, LAT)
    if m.sum() < 4: continue
    ser = np.array([np.nansum(Y[m, k] * w[m]) / np.sum(w[m] * np.isfinite(Y[m, k])) for k in range(len(anos))]); r = tend(anos, ser)
    if r and r[4] > 50: paises[nome] = dict(slope=r[0], lo=r[1], hi=r[2], p=r[3], media=r[4], n=int(m.sum()), serie=[round(float(v), 1) for v in ser])
nomes = list(paises); qp = fdr([paises[k]["p"] for k in nomes])
for k, qq in zip(nomes, qp):
    c = paises[k]; c["q"] = round(float(qq), 4); c["pct_decada"] = round(c["slope"] * 10 / c["media"] * 100, 2); c["sig"] = int(qq < 0.05 and (c["lo"] > 0 or c["hi"] < 0))
    for a in ("slope", "lo", "hi", "media"): c[a] = round(c[a], 2)
    c["p"] = round(c["p"], 4)
mais = sorted([k for k in nomes if paises[k]["sig"] and paises[k]["slope"] > 0], key=lambda k: -paises[k]["pct_decada"]); menos = sorted([k for k in nomes if paises[k]["sig"] and paises[k]["slope"] < 0], key=lambda k: paises[k]["pct_decada"])
print(f"\nPaíses: {len(nomes)} | significativamente mais húmidos: {len(mais)} | mais secos: {len(menos)}")
print("Mais húmidos:", [(k, paises[k]["pct_decada"]) for k in mais[:8]]); print("Mais secos:", [(k, paises[k]["pct_decada"]) for k in menos[:8]])
# ligação exploratória ao metano: chuva tropical (terra, 30S-30N) vs crescimento anual da NOAA, sem tendência
tr = np.isfinite(Y).all(1) & (np.abs(LAT) <= 30); terra = np.zeros(len(d), bool)
for f in json.load(open("world.json", encoding="utf-8"))["features"]:
    if f["properties"]["name"] != "Antarctica": terra |= contains_xy(shape(f["geometry"]), LON, LAT)
tm = tr & terra; chuva_trop = np.array([np.average(Y[tm, k], weights=w[tm]) for k in range(len(anos))])
NO = json.load(open("noaa_vs_airs.json"))["crescimento_anual_airs_noaa"]; ys = [int(a) for a in anos if str(a) in NO]
cn = np.array([NO[str(a)][1] for a in ys]); cp = np.array([chuva_trop[list(anos).index(a)] for a in ys])
dt = lambda x: x - np.polyval(np.polyfit(np.arange(len(x)), x, 1), np.arange(len(x)))
r0 = stats.pearsonr(dt(cp), dt(cn)); r1 = stats.pearsonr(dt(cp[:-1]), dt(cn[1:])); rs = stats.spearmanr(dt(cp), dt(cn))
print(f"\nChuva tropical em terra vs crescimento do metano (NOAA), {len(ys)} anos, sem tendência: Pearson r = {r0[0]:.2f} (p = {r0[1]:.2f}) | com 1 ano de atraso r = {r1[0]:.2f} (p = {r1[1]:.2f}) | Spearman {rs[0]:.2f}")
# séries regionais de exemplo: células significativamente mais húmidas ou mais secas dentro de um país (anomalia % face à média)
exemplos = []
sigmask = np.zeros(len(d), bool); sigsl = np.zeros(len(d)); sigmask[idx] = sig; sigsl[idx] = sl; meanall = np.nanmean(Y, axis=1)
for f in json.load(open("world.json", encoding="utf-8"))["features"]:
    nome = f["properties"]["name"]
    for alvo, sinal, rot in (("China", 1, "mais húmidas"), ("Brazil", -1, "mais secas"), ("Angola", -1, "mais secas"), ("Russia", 1, "mais húmidas"), ("Zambia", -1, "mais secas")):
        if nome != alvo: continue
        m = contains_xy(shape(f["geometry"]), LON, LAT) & sigmask & (sigsl * sinal > 0) & (meanall >= 300) & (np.abs(LAT) <= 60)
        if m.sum() < 3: continue
        sr = np.nanmean(Y[m], axis=0); pc = (sr / sr.mean() - 1) * 100; k = np.polyfit(anos, pc, 1)[0]
        exemplos.append(dict(pais=alvo, tipo=rot, n=int(m.sum()), anos=[int(a) for a in anos], pct=[round(float(x), 1) for x in pc], pct_ano=round(float(k), 2)))
print("Séries regionais:", [(e["pais"], e["tipo"], e["n"], e["pct_ano"]) for e in exemplos])
out = dict(exemplos=exemplos, periodo=[int(anos[0]), int(anos[-1])], fonte="NASA GPM IMERG mensal V07, totais anuais, células de 1°", n_celulas=len(idx), sig_pct=round(float(sig.mean()), 3),
           mais_humidas=int((sig & (sl > 0)).sum()), mais_secas=int((sig & (sl < 0)).sum()), mediana_pct_decada=round(float(np.median(pct)), 2),
           paises=paises, ligacao=dict(n=len(ys), pearson=round(float(r0[0]), 3), p=round(float(r0[1]), 3), pearson_atraso1=round(float(r1[0]), 3), p_atraso1=round(float(r1[1]), 3), spearman=round(float(rs[0]), 3),
                                      anos=ys, chuva_tropical=[round(float(x), 1) for x in cp], crescimento_noaa=[round(float(x), 2) for x in cn]),
           celulas=[[float(a), float(b), round(float(c), 2), int(s), int(m)] for a, b, c, s, m in zip(lat, lon, pct, sig, mean)])
json.dump(out, open("gpm_trends.json", "w", encoding="utf-8"), ensure_ascii=False); print("Gravado gpm_trends.json")
