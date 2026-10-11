"""QUANDO mudou? Tendência global do metano (NASA AIRS) por períodos, aceleração e faixas de latitude.
Entradas (mesma pasta): airs_ch4_400hPa.npz   (opcional: s5p_slope_2019_2026.csv para confirmar as faixas com o Sentinel-5P)
Uso: py airs_periodos.py          Saídas: periodos.json, periodos.png         Precisa: numpy scipy matplotlib pandas
"""
import json, os
import numpy as np, pandas as pd
from scipy import stats
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

z = np.load("airs_ch4_400hPa.npz"); meses = [str(m) for m in z["meses"]]; A = z["A"].astype("float64"); nt = len(meses)
lat = (89.5 - np.arange(180)) if bool(z["norte_em_cima"]) else (-89.5 + np.arange(180))
ano = np.array([int(m[:4]) + (int(m[5:]) - 1) / 12 for m in meses]); mes = np.array([int(m[5:]) - 1 for m in meses])
for m in range(12): A[mes == m] -= np.nanmean(A[mes == m], axis=0)               # tira o ciclo sazonal
W = np.cos(np.radians(lat))[:, None] * np.ones((1, 360))

def serie(mask_lat):
    X = A[:, mask_lat, :]; w = W[mask_lat]; v = np.isfinite(X)
    return (np.where(v, X, 0) * w).sum((1, 2)) / np.maximum((v * w).sum((1, 2)), 1e-9)

def ar1(r): return float(min(max(np.corrcoef(r[:-1], r[1:])[0, 1], 0), 0.9))

def tendencia(t, y):
    """Mínimos quadrados com erro-padrão corrigido por autocorrelação (n efetivo) + Mann-Kendall corrigido."""
    n = len(y); b, a = np.polyfit(t, y, 1); r = y - (a + b * t); r1 = ar1(r); f = np.sqrt((1 + r1) / (1 - r1))
    se = np.sqrt(((r ** 2).sum() / (n - 2)) / ((t - t.mean()) ** 2).sum()) * f; h = 1.96 * se
    tau = stats.kendalltau(t, y)[0]; S = tau * n * (n - 1) / 2; var = n * (n - 1) * (2 * n + 5) / 18 * (1 + r1) / (1 - r1)
    p = 2 * (1 - stats.norm.cdf(abs((S - np.sign(S)) / np.sqrt(var))))
    return dict(slope=round(b, 2), lo=round(b - h, 2), hi=round(b + h, 2), p=round(float(p), 5), n=n, r1=round(r1, 2))

g = serie(np.ones(180, bool)); periodos = [("2002-09", "2008-12"), ("2009-01", "2018-12"), ("2019-01", "2026-08"), ("2002-09", "2026-08")]
def por_periodo(y):
    out = {}
    for a, b in periodos:
        s = np.array([a <= m <= b for m in meses]); out[f"{a} a {b}"] = tendencia(ano[s], y[s])
    return out
res = dict(global_=por_periodo(g))
print("TENDÊNCIA GLOBAL por período (ppb/ano, IC 95% corrigido por autocorrelação):")
for k, v in res["global_"].items(): print(f"  {k}: {v['slope']:+.2f}  [{v['lo']:+.2f}, {v['hi']:+.2f}]  p={v['p']:.3g}")

# aceleração: ajuste quadrático + bootstrap por blocos (24 meses)
t0 = ano - ano.mean(); c = np.polyfit(t0, g, 2); fit = np.polyval(c, t0); r = g - fit; rng = np.random.default_rng(0); acc = []
for _ in range(1000):
    idx = np.concatenate([np.arange(s, s + 24) % nt for s in rng.integers(0, nt, nt // 24 + 1)])[:nt]
    acc.append(np.polyfit(t0, fit + r[idx], 2)[0] * 2)
res["aceleracao"] = dict(valor=round(float(c[0] * 2), 3), lo=round(float(np.percentile(acc, 2.5)), 3), hi=round(float(np.percentile(acc, 97.5)), 3),
                         significativa=bool(np.percentile(acc, 2.5) > 0 or np.percentile(acc, 97.5) < 0))
print(f"\nACELERAÇÃO: {res['aceleracao']['valor']:+.3f} ppb/ano² [{res['aceleracao']['lo']:+.3f}, {res['aceleracao']['hi']:+.3f}] "
      f"-> {'significativa' if res['aceleracao']['significativa'] else 'NÃO significativa'}")

# tendência numa janela móvel de 60 meses
roll = []
for i in range(0, nt - 60 + 1, 6):
    tt = tendencia(ano[i:i + 60], g[i:i + 60]); roll.append([round(float(ano[i:i + 60].mean()), 2), tt["slope"], tt["lo"], tt["hi"]])
res["janela_5_anos"] = roll

# faixas de latitude
faixas = {"Norte (30°N-90°N)": lat >= 30, "Trópicos (30°S-30°N)": (lat < 30) & (lat > -30), "Sul (90°S-30°S)": lat <= -30}
res["faixas"] = {k: por_periodo(serie(m)) for k, m in faixas.items()}
print("\nFAIXAS DE LATITUDE (2019-2026):")
if os.path.exists("s5p_slope_2019_2026.csv"):
    s5 = pd.read_csv("s5p_slope_2019_2026.csv"); s5 = s5[s5.slope_ppb_ano.abs() < 60]
for k, m in faixas.items():
    a = res["faixas"][k]["2019-01 a 2026-08"]; txt = f"  {k}: AIRS {a['slope']:+.2f} [{a['lo']:+.2f}, {a['hi']:+.2f}]"
    if os.path.exists("s5p_slope_2019_2026.csv"):
        mm = (s5.lat >= 30) if k.startswith("Norte") else (s5.lat <= -30) if k.startswith("Sul") else ((s5.lat < 30) & (s5.lat > -30))
        v = float(np.average(s5.slope_ppb_ano[mm], weights=np.cos(np.radians(s5.lat[mm])))); res["faixas"][k]["s5p_2019_2026"] = round(v, 2); txt += f" | Sentinel-5P {v:+.2f}"
    print(txt)
res["serie_anual"] = [[int(a), round(float(np.nanmean(g[(ano >= a) & (ano < a + 1)])), 2)] for a in range(int(ano[0]), int(ano[-1]) + 1)]
json.dump(res, open("periodos.json", "w", encoding="utf-8"), ensure_ascii=False)

fig, ax = plt.subplots(1, 2, figsize=(11, 4)); ax[0].plot(ano, g, c="#1f7a4d", lw=1)
cores = ["#2166ac", "#e6a117", "#b2182b"]
for (a, b), cor in zip(periodos[:3], cores):
    s = np.array([a <= m <= b for m in meses]); k = np.polyfit(ano[s], g[s], 1); ax[0].plot(ano[s], np.polyval(k, ano[s]), c=cor, lw=2.5, label=f"{a[:4]}-{b[:4]}: {k[0]:+.1f} ppb/ano")
ax[0].set(title="Metano global (AIRS, 400 hPa), sem sazonalidade", ylabel="anomalia (ppb)"); ax[0].legend(fontsize=8)
R = np.array(roll); ax[1].plot(R[:, 0], R[:, 1], c="#1f7a4d"); ax[1].fill_between(R[:, 0], R[:, 2], R[:, 3], alpha=.25, color="#1f7a4d")
ax[1].axhline(0, c="#999", lw=.6); ax[1].set(title="Tendência em janelas de 5 anos (IC 95%)", ylabel="ppb/ano", xlabel="ano (centro da janela)")
fig.tight_layout(); fig.savefig("periodos.png", dpi=130); print("\nGravado periodos.json e periodos.png")
