"""«Quando mudou?» do metano global, SÓ com dados da NASA (Aqua AIRS, 400 hPa).
A partir do metano mensal (airs_ch4_400hPa.npz do airs_baixar.py) calcula a média global pesada pela área, tira a estação
(média de cada mês do ano) e grava periodos.json no formato que o build_site.py e a janela «Quando mudou?» usam:
  global_       tendência de Sen (ppb/ano) por período, com IC 95 %, p de Mann-Kendall (corrigido pela autocorrelação) e r1
  aceleracao    termo quadrático (ppb/ano em cada ano) com IC 95 %
  janela_5_anos tendência em janelas de 60 meses, de 6 em 6 meses: [centro, ppb/ano, IC baixo, IC alto]
  faixas        as mesmas tendências por faixa de latitude (Norte, Trópicos, Sul)
  serie_anual   média anual da anomalia global (ppb)
Grava também as mesmas tendências por zona, que a janela mostra conforme a zona escolhida no mapa:
  periodos_zonas.json  continentes e países
  periodos_blocos.bin  blocos de 1°
Uso: py airs_periodos.py            Precisa: numpy shapely, world.json e continentes.json
"""
import json, math
import numpy as np

z = np.load("airs_ch4_400hPa.npz"); ms = [str(m) for m in z["meses"]]; A = z["A"].astype("float64")
if not bool(z["norte_em_cima"]): A = A[:, ::-1, :]
lat = 89.5 - np.arange(180)
w = np.cos(np.radians(lat))[:, None] * np.ones((1, 360))


def media_faixa(la0, la1):
    """média mensal pesada pela área entre as latitudes la0 e la1 (meses com pouca cobertura ficam em falta)"""
    m = (lat >= la0) & (lat < la1); W = w[m]; S = []
    for k in range(len(ms)):
        X = A[k][m]; ok = np.isfinite(X)
        S.append(np.nan if ok.mean() < .3 else float((X[ok] * W[ok]).sum() / W[ok].sum()))
    return np.array(S)


def sem_estacao(s):
    """anomalia: tira a média de cada mês do ano"""
    mes = np.array([int(x[5:7]) for x in ms]); a = s.copy()
    for k in range(1, 13):
        sel = (mes == k) & np.isfinite(s)
        if sel.any(): a[mes == k] -= np.nanmean(s[sel])
    return a


T = np.array([int(x[:4]) + (int(x[5:7]) - .5) / 12 for x in ms])   # tempo em anos (meio do mês)


def r1(y):
    y = y - y.mean(); return float((y[1:] * y[:-1]).sum() / (y * y).sum()) if len(y) > 3 else 0.0


def mq(t, y, grau, r_=False):
    """mínimos quadrados com erro padrão inflacionado pela autocorrelação: devolve (coeficiente mais alto, erro[, r1])"""
    ok = np.isfinite(y); t, y = t[ok] - t[ok].mean(), y[ok]; X = np.vander(t, grau + 1)
    b, *_ = np.linalg.lstsq(X, y, rcond=None); res = y - X @ b; r = r1(res); n = len(y)
    s2 = (res @ res) / (n - grau - 1) * max(1.0, (1 + r) / (1 - r)); cov = s2 * np.linalg.inv(X.T @ X)
    return (float(b[0]), math.sqrt(cov[0, 0]), r) if r_ else (float(b[0]), math.sqrt(cov[0, 0]))


def sen(t, y):
    """tendência de Sen (por ano); IC 95 % e p a partir do erro da regressão, corrigido pela autocorrelação (como nas janelas de 5 anos)"""
    ok = np.isfinite(y); t, y = t[ok], y[ok]; n = len(y)
    i, j = np.triu_indices(n, 1); s = float(np.median((y[j] - y[i]) / (t[j] - t[i])))
    _, e, r = mq(t, y, 1, True); p = math.erfc(abs(s / e) / math.sqrt(2)) if e > 0 else 0.0
    return dict(slope=round(s, 2), lo=round(s - 1.96 * e, 2), hi=round(s + 1.96 * e, 2), p=round(p, 5), n=n, r1=round(r, 2))


fim = len(ms) - 1
G = sem_estacao(media_faixa(-90, 90))
while fim > 0 and not np.isfinite(G[fim]): fim -= 1
ini = int(np.argmax(np.isfinite(G)))
nome = lambda a, b: f"{ms[a]} a {ms[b]}"
i09, i19 = ms.index("2009-01"), ms.index("2019-01")
PER = [(ini, i09 - 1), (i09, i19 - 1), (i19, fim), (ini, fim)]


def tendencias(s):
    return {nome(a, b): sen(T[a:b + 1], s[a:b + 1]) for a, b in PER}


out = {"global_": tendencias(G)}
c, e = mq(T[ini:fim + 1], G[ini:fim + 1], 2)                       # y = c·t² + … → aceleração = 2c
out["aceleracao"] = dict(valor=round(2 * c, 3), lo=round(2 * (c - 1.96 * e), 3), hi=round(2 * (c + 1.96 * e), 3), significativa=bool(abs(c) > 1.96 * e))
jan = []
for a in range(ini, fim - 59 + 1, 6):
    b, eb = mq(T[a:a + 60], G[a:a + 60], 1)
    jan.append([round(float(T[a:a + 60].mean()), 2), round(b, 2), round(b - 1.96 * eb, 2), round(b + 1.96 * eb, 2)])
out["janela_5_anos"] = jan
out["faixas"] = {k: tendencias(sem_estacao(media_faixa(a, b))) for k, (a, b) in
                 {"Norte (30°N-90°N)": (30, 90), "Trópicos (30°S-30°N)": (-30, 30), "Sul (90°S-30°S)": (-90, -30)}.items()}
anos = sorted({int(x[:4]) for x in ms[ini:fim + 1]})
out["serie_anual"] = [[y, round(float(np.nanmean(G[[k for k in range(ini, fim + 1) if ms[k].startswith(str(y))]])), 2)] for y in anos
                      if sum(1 for k in range(ini, fim + 1) if ms[k].startswith(str(y)) and np.isfinite(G[k])) >= 6]
json.dump(out, open("periodos.json", "w", encoding="utf-8"), ensure_ascii=False)
pk = max(jan, key=lambda r: r[1])
print(f"periodos.json: {ms[ini]} a {ms[fim]} | global {out['global_'][nome(ini, fim)]['slope']:+.2f} ppb/ano | aceleração {out['aceleracao']['valor']:+.3f} | pico janela {pk[0]:.1f}: {pk[1]:+.2f}")


# ---------- o mesmo por continente e por país (periodos_zonas.json) ----------
def zona(mask):
    """resumo de uma zona: g = [ppb/ano, IC baixo, IC alto, p] por período; ac = aceleração; j = janelas de 5 anos, de ano a ano"""
    Wm = w[mask]; X = A[:, mask]; ok = np.isfinite(X); peso = (ok * Wm).sum(1)
    s = np.where(ok.mean(1) >= .3, (np.nan_to_num(X) * Wm).sum(1) / np.where(peso > 0, peso, 1), np.nan)
    a = sem_estacao(s)
    if np.isfinite(a[ini:fim + 1]).sum() < 60: return None
    g = []
    for p, q in PER:
        r = sen(T[p:q + 1], a[p:q + 1]) if np.isfinite(a[p:q + 1]).sum() >= 24 else None
        g.append([r["slope"], r["lo"], r["hi"], r["p"]] if r else None)
    c_, e_ = mq(T[ini:fim + 1], a[ini:fim + 1], 2)
    J = []
    for p in range(ini, fim - 59 + 1, 12):
        if np.isfinite(a[p:p + 60]).sum() < 36: continue
        b_, eb = mq(T[p:p + 60], a[p:p + 60], 1)
        J.append([round(float(T[p:p + 60].mean()), 1), round(b_, 2), round(b_ - 1.96 * eb, 2), round(b_ + 1.96 * eb, 2)])
    return dict(g=g, ac=[round(2 * c_, 3), round(2 * (c_ - 1.96 * e_), 3), round(2 * (c_ + 1.96 * e_), 3)], j=J)


import shapely
from shapely.geometry import MultiPolygon, box, shape
from shapely.ops import transform
from shapely.strtree import STRtree


def corta_linha_data(s):
    """países que atravessam ±180° (Rússia, Fiji...) viriam como uma faixa à volta do mundo: divide-os em duas partes (como no build_site.py)"""
    partes = list(s.geoms) if s.geom_type == "MultiPolygon" else [s]; res = []
    for p in partes:
        b = p.bounds
        if b[2] - b[0] <= 180: res.append(p); continue
        q = transform(lambda x, y, z=None: (np.where(np.asarray(x) < 0, np.asarray(x) + 360, x), y), p).buffer(0)
        for bx, sh in ((box(0, -90, 180, 90), 0), (box(180, -90, 360, 90), -360)):
            r = q.intersection(bx)
            if r.is_empty: continue
            r = transform(lambda x, y, z=None: (np.asarray(x) + sh, y), r)
            res += [x for x in (r.geoms if hasattr(r, "geoms") else [r]) if x.geom_type == "Polygon" and x.area > 0]
    return MultiPolygon(res)


feats = json.load(open("world.json", encoding="utf-8"))["features"]
geoms = [corta_linha_data(shape(f["geometry"])) for f in feats]; nomes = [f["properties"]["name"] for f in feats]
lon = -179.5 + np.arange(360); LON, LAT = np.meshgrid(lon, lat)
a_, b_ = STRtree(geoms).query(shapely.points(LON.ravel(), LAT.ravel()), predicate="within")
pais = np.full(LAT.size, -1); pais[a_] = b_; pais = pais.reshape(LAT.shape)
paises = {}
for k, nm in enumerate(nomes):
    m = pais == k
    if not m.any():                                                # país mais pequeno do que um bloco: usa o bloco onde ele fica
        rp = geoms[k].representative_point(); m[min(179, int(89.5 - rp.y + .5)), min(359, int(rp.x + 179.5 + .5))] = True
    r = zona(m)
    if r: paises[nm] = r
try:
    cont_de = json.load(open("continentes.json", encoding="utf-8"))
except OSError:
    cont_de = {}
continentes = {}
for ct in sorted(set(cont_de.values())):
    idx = [k for k, nm in enumerate(nomes) if cont_de.get(nm) == ct]
    r = zona(np.isin(pais, idx))
    if r: continentes[ct] = r
json.dump(dict(per=[nome(a, b) for a, b in PER], continentes=continentes, paises=paises), open("periodos_zonas.json", "w", encoding="utf-8"), ensure_ascii=False, separators=(",", ":"))


# ---------- o mesmo por bloco de 1° (periodos_blocos.bin) ----------
# regressão linear com o mesmo erro corrigido pela autocorrelação (Sen em 64 800 blocos demoraria minutos todos os dias)
# formato: Int16 por bloco (linha·360 + coluna, com o Norte primeiro), 2 valores por período [ppb/ano·100, meia largura do IC 95 %·100]; -32768 = sem dados
mes = np.array([int(x[5:7]) for x in ms])
An = A.copy()
for k in range(1, 13):
    sel = mes == k
    if sel.any():
        with np.errstate(invalid="ignore"):
            import warnings
            with warnings.catch_warnings():
                warnings.simplefilter("ignore"); An[sel] -= np.nanmean(A[sel], 0)
B = np.full((180, 360, len(PER) * 2), -32768, dtype=np.int16)
for n, (p, q) in enumerate(PER):
    Y = An[p:q + 1]; t = T[p:q + 1]; ok = np.isfinite(Y); cnt = ok.sum(0); Z = np.where(ok, Y, 0.0)
    tt = (t[:, None, None] * ok); tm = tt.sum(0) / np.maximum(cnt, 1); ym = Z.sum(0) / np.maximum(cnt, 1)
    dt = (t[:, None, None] - tm) * ok; dy = (Z - ym) * ok
    sxx = (dt * dt).sum(0); bb = (dt * dy).sum(0) / np.where(sxx > 0, sxx, 1)
    res = (dy - bb * dt) * ok
    ac0 = (res * res).sum(0); ac1 = (res[1:] * res[:-1]).sum(0)
    rr = np.clip(ac1 / np.where(ac0 > 0, ac0, 1), 0, .95)
    s2 = ac0 / np.maximum(cnt - 2, 1) * (1 + rr) / (1 - rr)
    ee = np.sqrt(s2 / np.where(sxx > 0, sxx, 1)) * 1.96
    bom = (cnt >= max(24, .6 * len(t))) & (sxx > 0)
    B[..., n * 2][bom] = np.clip(np.round(bb[bom] * 100), -32000, 32000)
    B[..., n * 2 + 1][bom] = np.clip(np.round(ee[bom] * 100), 0, 32000)
B.astype("<i2").tofile("periodos_blocos.bin")
print(f"periodos_zonas.json: {len(continentes)} continentes, {len(paises)} de {len(nomes)} países | periodos_blocos.bin: {(B[..., 6] != -32768).sum()} blocos")
