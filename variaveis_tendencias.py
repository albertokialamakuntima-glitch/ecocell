"""Seletor de variáveis: calcula, para TODAS as variáveis, o que muda, onde, quanto e se é significativo — com o mesmo método.

Método (igual para todas): médias anuais por célula de 1°; tendência de Sen com IC 95 %; Mann-Kendall com a variância
corrigida pela autocorrelação lag-1 dos resíduos; correção FDR de Benjamini-Hochberg (q < 0,05) em cada variável.
Também por país (média ponderada pela área das células do país) e para o conjunto (terra, mar ou mundo).

Entradas (mesma pasta): variaveis_anuais.npz (do variaveis_exportar.py), airs_ch4_400hPa.npz, gpm_anual_1deg.csv, world.json
Saídas:  variaveis.json        (resumo que vai dentro do site: mapas, países, rankings)
         vars/<id>.json        (estatísticas e séries anuais de cada célula; o site só as carrega quando precisas)
Uso:     py variaveis_tendencias.py
Precisa: numpy scipy shapely
"""
import json, os, base64
import numpy as np
from scipy import stats
from shapely.geometry import shape
from shapely.strtree import STRtree
import shapely

# ---------- metadados das variáveis (nomes em pt/en, unidades, missão, onde faz sentido) ----------
META = {
 "ch4":  dict(pt="Metano (troposfera média)", en="Methane (mid-troposphere)", un="ppb", tr="ppb/década", tr_en="ppb/decade", fonte="NASA Aqua AIRS (400 hPa)", tipo="satélite", dom="tudo", d=1),
 "chuva":dict(pt="Chuva anual", en="Annual rainfall", un="mm/ano", tr="%/década", tr_en="%/decade", fonte="NASA/JAXA GPM IMERG", tipo="satélite", dom="tudo", d=1, pct=True),
 "tar":  dict(pt="Temperatura do ar (2 m)", en="Air temperature (2 m)", un="°C", tr="°C/década", tr_en="°C/decade", fonte="Modelo NASA FLDAS (Noah, forçado por MERRA-2/CHIRPS)", tipo="modelo", dom="terra", d=2),
 "lst":  dict(pt="Temperatura da superfície (dia)", en="Land surface temperature (day)", un="°C", tr="°C/década", tr_en="°C/decade", fonte="NASA Terra MODIS (MOD11A2)", tipo="satélite", dom="terra", d=2, ate=2020,
              nota="Só até 2020: desde então o satélite Terra deixou de corrigir a órbita e passa mais cedo, numa hora mais fresca, o que criaria um falso arrefecimento.",
              nota_en="Only up to 2020: since then the Terra satellite stopped correcting its orbit and passes earlier, at a cooler hour, which would create a false cooling."),
 "ndvi": dict(pt="Vegetação (NDVI)", en="Vegetation (NDVI)", un="NDVI", tr="NDVI/década", tr_en="NDVI/decade", fonte="NASA Terra MODIS (MOD13A2)", tipo="satélite", dom="terra", d=3),
 "et":   dict(pt="Evapotranspiração", en="Evapotranspiration", un="mm/ano", tr="mm/ano por década", tr_en="mm/yr per decade", fonte="NASA Terra MODIS (MOD16A2GF)", tipo="satélite", dom="terra", d=1),
 "solo": dict(pt="Humidade do solo (0–10 cm)", en="Soil moisture (0–10 cm)", un="m³/m³", tr="m³/m³ por década", tr_en="m³/m³ per decade", fonte="Modelo NASA FLDAS (Noah)", tipo="modelo", dom="terra", d=4),
 "agua": dict(pt="Água total armazenada", en="Total water storage", un="cm de água", tr="cm/década", tr_en="cm/decade", fonte="NASA/DLR GRACE e GRACE-FO (mascons JPL)", tipo="satélite", dom="terra", d=2),
 "sst":  dict(pt="Temperatura do mar", en="Sea surface temperature", un="°C", tr="°C/década", tr_en="°C/decade", fonte="NASA Aqua MODIS (L3SMI)", tipo="satélite", dom="mar", d=2),
 "clor": dict(pt="Clorofila do oceano (fitoplâncton)", en="Ocean chlorophyll (phytoplankton)", un="mg/m³", tr="%/década", tr_en="%/decade", fonte="NASA Aqua MODIS (L3SMI)", tipo="satélite", dom="mar", d=1, pct_log=True),
}

lat = 89.5 - np.arange(180); lon = -179.5 + np.arange(360)
LON, LAT = np.meshgrid(lon, lat); W = np.cos(np.radians(LAT))

# ---------- dados: tudo para cubos (anos, 180, 360), norte em cima ----------
cubos = {}
z = np.load("variaveis_anuais.npz"); anos = [int(a) for a in z["anos"]]
for k in ["tar", "lst", "ndvi", "et", "solo", "agua", "sst", "clor"]:   # (a chuva é lida mais abaixo)
    if k in z.files:
        fim = META[k].get("ate", anos[-1]); s = [i for i, a in enumerate(anos) if a <= fim]
        cubos[k] = ([anos[i] for i in s], z[k][s].astype("float64"))
if os.path.exists("airs_ch4_400hPa.npz"):                       # metano: meses -> médias anuais (anos com 10+ meses)
    a = np.load("airs_ch4_400hPa.npz"); ms = [str(m) for m in a["meses"]]; A = a["A"].astype("float64")
    if not bool(a["norte_em_cima"]): A = A[:, ::-1, :]
    aa = sorted({int(m[:4]) for m in ms}); cub = []
    for y in aa:
        s = [i for i, m in enumerate(ms) if m.startswith(str(y))]
        cub.append(np.nanmean(A[s], 0) if len(s) >= 10 else np.full((180, 360), np.nan))
    ok = [i for i, y in enumerate(aa) if len([m for m in ms if m.startswith(str(y))]) >= 10]
    cubos["ch4"] = ([aa[i] for i in ok], np.array([cub[i] for i in ok]))
if "chuva" in z.files: cubos["chuva"] = (anos, z["chuva"].astype("float64"))   # chuva do GPM pelo Earth Engine (variaveis_exportar.py)
elif os.path.exists("gpm_anual_1deg.csv"):                      # ou da tabela antiga (lat, lon, y2003...) -> cubo
    import csv
    rd = list(csv.reader(open("gpm_anual_1deg.csv", encoding="utf-8"))); cab = rd[0]; ya = [int(c[1:]) for c in cab[2:]]
    C = np.full((len(ya), 180, 360), np.nan)
    for r in rd[1:]:
        la, lo = float(r[0]), float(r[1]); i = int(np.clip(np.floor(90 - la), 0, 179)); j = int(np.clip(np.floor(lo + 180), 0, 359))
        C[:, i, j] = [float(x) if x not in ("", "nan") else np.nan for x in r[2:]]
    cubos["chuva"] = (ya, C)
for k, (ya, cubo) in list(cubos.items()):   # tira os anos finais sem dados (a fonte ainda não tem o ano completo), para não anunciar anos que não existem
    cobre = np.isfinite(cubo).reshape(len(ya), -1).sum(1); n = len(ya)
    while n > 1 and cobre[n - 1] < 0.5 * cobre[:n - 1].max(): n -= 1
    cubos[k] = (ya[:n], cubo[:n])

# ---------- países: a que país pertence o centro de cada célula ----------
from shapely.geometry import MultiPolygon, box
from shapely.ops import transform
def corta_linha_data(s):
    """Países que atravessam ±180° (Rússia, Fiji...) viriam como uma faixa à volta do mundo: divide-os em duas partes (como no build_site.py)."""
    partes = list(s.geoms) if s.geom_type == "MultiPolygon" else [s]; out = []
    for p in partes:
        b = p.bounds
        if b[2] - b[0] <= 180: out.append(p); continue
        q = transform(lambda x, y, z=None: (np.where(np.asarray(x) < 0, np.asarray(x) + 360, x), y), p).buffer(0)
        for bx, sh in ((box(0, -90, 180, 90), 0), (box(180, -90, 360, 90), -360)):
            r = q.intersection(bx)
            if r.is_empty: continue
            r = transform(lambda x, y, z=None: (np.asarray(x) + sh, y), r)
            out += [x for x in (r.geoms if hasattr(r, "geoms") else [r]) if x.geom_type == "Polygon" and x.area > 0]
    return MultiPolygon(out)
feats = json.load(open("world.json", encoding="utf-8"))["features"]
geoms = [corta_linha_data(shape(f["geometry"])) for f in feats]; nomes = [f["properties"]["name"] for f in feats]
pts = shapely.points(LON.ravel(), LAT.ravel()); arv = STRtree(geoms)
pais = np.full(LAT.size, -1); a_, b_ = arv.query(pts, predicate="within"); pais[a_] = b_
pais = pais.reshape(LAT.shape); terra = pais >= 0

# ---------- estatística vetorizada (todas as células de uma vez) ----------
def tendencias(t, Y):
    """t: (n,) anos; Y: (n, m) séries. Devolve Sen, IC95 baixo/alto, p de Mann-Kendall corrigido, n válidos."""
    n, m = Y.shape; i, j = np.triu_indices(n, 1)
    d = Y[j] - Y[i]; dt = (t[j] - t[i])[:, None]
    sl = d / dt                                                        # declives de todos os pares
    sen = np.nanmedian(sl, 0)
    S = np.nansum(np.sign(d), 0); nv = np.isfinite(Y).sum(0)
    var = nv * (nv - 1) * (2 * nv + 5) / 18.0
    ic = np.nanmedian(Y - sen * t[:, None], 0); res = Y - (ic + sen * t[:, None])   # resíduos da reta de Sen
    r0, r1v = res[:-1], res[1:]; ok = np.isfinite(r0) & np.isfinite(r1v)
    a0 = np.where(ok, r0 - np.nanmean(res, 0), 0); a1 = np.where(ok, r1v - np.nanmean(res, 0), 0)
    r1 = np.clip(np.sum(a0 * a1, 0) / np.maximum(np.sum(np.where(np.isfinite(res), (res - np.nanmean(res, 0)) ** 2, 0), 0), 1e-12), 0, .9)
    var = var * (1 + r1) / (1 - r1)
    zz = (S - np.sign(S)) / np.sqrt(np.maximum(var, 1e-12)); p = 2 * (1 - stats.norm.cdf(np.abs(zz)))
    # IC 95 % de Sen (método por ordem dos declives, com a variância corrigida)
    N = np.isfinite(sl).sum(0); Ca = 1.96 * np.sqrt(np.maximum(var, 0)); srt = np.sort(sl, 0)       # NaN vão para o fim
    lo_i = np.clip(np.floor((N - Ca) / 2).astype(int), 0, None); hi_i = np.clip(np.ceil((N + Ca) / 2).astype(int), 0, np.maximum(N - 1, 0))
    col = np.arange(m); lo = srt[np.minimum(lo_i, sl.shape[0] - 1), col]; hi = srt[np.minimum(hi_i, sl.shape[0] - 1), col]
    return sen, lo, hi, p, nv

def fdr(p):
    p = np.asarray(p, float); q = np.full_like(p, np.nan); ok = np.isfinite(p); pv = p[ok]; o = np.argsort(pv); mm = len(pv)
    if mm: qq = np.minimum.accumulate((pv[o] * mm / (np.arange(mm) + 1))[::-1])[::-1]; tmp = np.empty(mm); tmp[o] = np.minimum(qq, 1); q[ok] = tmp
    return q

def b64(a): return base64.b64encode(np.ascontiguousarray(a).tobytes()).decode()

os.makedirs("vars", exist_ok=True); saida = {"anos_ref": anos, "vars": []}
for k, M in META.items():
    if k not in cubos: print("sem dados:", k); continue
    ya, cubo = cubos[k]; t = np.array(ya, float); n = len(ya)
    dom = terra if M["dom"] == "terra" else (~terra) if M["dom"] == "mar" else np.ones_like(terra)
    validos = (np.isfinite(cubo).sum(0) >= max(8, int(.7 * n))) & dom
    if k in ("chuva",): validos &= (np.nanmean(cubo, 0) >= 100)        # em desertos extremos a percentagem não faz sentido
    idx = np.where(validos.ravel())[0]; Y = cubo.reshape(n, -1)[:, idx]
    sen, lo, hi, p, nv = tendencias(t, Y); q = fdr(p)
    media = np.nanmean(Y, 0)
    # «quanto»: por década; nas variáveis em % é relativo à média (clorofila: está em log10)
    def conv(x, med):
        if M.get("pct"): return x * 10 / np.maximum(med, 1e-9) * 100
        if M.get("pct_log"): return (10 ** (x * 10) - 1) * 100
        return x * 10
    S, L, H = conv(sen, media), conv(lo, media), conv(hi, media)
    sig = (q < 0.05) & ((L > 0) | (H < 0)); sinal = np.where(sig, np.sign(S), 0)
    wv = W.ravel()[idx]; tot = wv.sum()
    sobe = float((wv * (sinal > 0)).sum() / tot * 100); desce = float((wv * (sinal < 0)).sum() / tot * 100)
    esc = float(np.nanpercentile(np.abs(S), 95)) or 1.0
    R = np.full(180 * 360, -128, np.int8); R[idx] = np.clip(np.round(S / esc * 120), -127, 127).astype(np.int8)
    MK = np.zeros(180 * 360, bool); MK[idx] = sig; bits = np.packbits(MK)
    # séries do conjunto (terra, mar ou mundo) e por país: média ponderada das anomalias
    def serie(mask_idx):
        Yv = Y[:, mask_idx]; wv_ = wv[mask_idx]; an = Yv - np.nanmean(Yv, 0); v = np.isfinite(an)
        return np.where(v.sum(1) > 0, (np.where(v, an, 0) * wv_).sum(1) / np.maximum((v * wv_).sum(1), 1e-12), np.nan), float(np.nanmean(Yv))
    def resumo(sr, med):
        s_, l_, h_, p_, _ = tendencias(t, sr[:, None]); s_, l_, h_ = conv(s_, np.array([med])), conv(l_, np.array([med])), conv(h_, np.array([med]))
        return [round(float(s_[0]), 4), round(float(l_[0]), 4), round(float(h_[0]), 4), float(p_[0])]
    g, gm = serie(np.arange(len(idx))); G = resumo(g, gm)
    paises = {}; pi = pais.ravel()[idx]; pp = []
    for c in np.unique(pi[pi >= 0]):
        mi = np.where(pi == c)[0]
        if len(mi) < 3: continue
        sr, md = serie(mi)
        if np.isfinite(sr).sum() < max(8, int(.7 * n)): continue
        r = resumo(sr, md); paises[nomes[c]] = r + [round(md, M["d"] + 2), [None if not np.isfinite(x) else round(float(x), M["d"] + 2) for x in sr]]; pp.append(r[3])
    for (nm, v), qq in zip(list(paises.items()), fdr(pp)):
        v.insert(4, round(float(qq), 5)); v.insert(5, int(qq < 0.05 and (v[1] > 0 or v[2] < 0)) * int(np.sign(v[0])))   # [s, lo, hi, p, q, sinal, média, série]
    saida["vars"].append(dict(id=k, nome=dict(pt=M["pt"], en=M["en"]), un=M["un"], tr=dict(pt=M["tr"], en=M["tr_en"]), fonte=M["fonte"], tipo=M["tipo"], dom=M["dom"], nota=dict(pt=M.get("nota", ""), en=M.get("nota_en", "")),
        anos=[ya[0], ya[-1]], n=len(idx), sobe=round(sobe, 1), desce=round(desce, 1), esc=round(esc, 5), r=b64(R), m=b64(bits),
        glob=dict(s=G[0], lo=G[1], hi=G[2], p=round(G[3], 6), serie=[None if not np.isfinite(x) else round(float(x), M["d"] + 2) for x in g], media=round(gm, M["d"] + 2)), paises=paises))
    # ficheiro por célula (carregado só quando preciso): estatísticas + séries quantizadas em int16
    prec = 10 ** (-(M["d"] + 1)); base = np.where(np.isfinite(media), media, 0)
    Q = np.where(np.isfinite(Y), np.clip(np.round((Y - base) / prec), -32767, 32767), -32768).astype(np.int16)
    json.dump(dict(id=k, anos=ya, prec=prec, cel=idx.tolist(), s=np.round(S, 5).tolist(), lo=np.round(L, 5).tolist(), hi=np.round(H, 5).tolist(),
                   q=np.round(q, 5).tolist(), med=np.round(media, M["d"] + 2).tolist(), serie=b64(Q.T)), open(f"vars/{k}.json", "w", encoding="utf-8"), separators=(",", ":"))
    print(f"{k:6s} {len(idx):6d} células | sobe {sobe:5.1f} % | desce {desce:5.1f} % da área (significativo) | global {G[0]:+.3f} {M['tr']} (p={G[3]:.3g}) | {len(paises)} países")
# variáveis sem dados nesta execução (ex.: o metano do AIRS na atualização automática, que não tem o ficheiro local): mantém as anteriores
if os.path.exists("variaveis.json"):
    feitos = {v["id"] for v in saida["vars"]}
    for v in json.load(open("variaveis.json", encoding="utf-8")).get("vars", []):
        if v["id"] not in feitos and v["id"] in META and os.path.exists(f"vars/{v['id']}.json"): saida["vars"].append(v); print("mantida a versão anterior de", v["id"])
    ordem = list(META); saida["vars"].sort(key=lambda v: ordem.index(v["id"]))
json.dump(saida, open("variaveis.json", "w", encoding="utf-8"), ensure_ascii=False, separators=(",", ":"))
print("Gravado variaveis.json e vars/*.json")
