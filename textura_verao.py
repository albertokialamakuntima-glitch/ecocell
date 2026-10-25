"""Converte a textura Blue Marble (com neve de inverno) numa Terra de "verão", como as imagens de referência:
tira a neve das terras de latitude média e alta (exceto Gronelândia, Antártida e Ártico profundo),
reaviva o verde e o castanho e aclara um pouco o oceano.  Uso: py verao.py entrada.jpg world.json saida.jpg"""
import sys, json
import numpy as np
from PIL import Image, ImageDraw, ImageFilter

ent, wj, sai = sys.argv[1:4]
im = Image.open(ent).convert("RGB"); W, H = im.size
A = np.asarray(im).astype(np.float32)
lat = 90 - (np.arange(H) + .5) / H * 180
LAT = np.repeat(lat[:, None], W, 1)

def px(lon, la): return ((lon + 180) / 360 * W, (90 - la) / 180 * H)
terra, gelo = Image.new("L", (W, H)), Image.new("L", (W, H))
for f in json.load(open(wj, encoding="utf-8"))["features"]:
    g = f["geometry"]; polys = [g["coordinates"]] if g["type"] == "Polygon" else g["coordinates"]
    alvo = gelo if f["properties"]["name"] in ("Greenland", "Antarctica") else terra
    for p in polys:
        for k, anel in enumerate(p):
            lo = np.array([c[0] for c in anel]); la = [c[1] for c in anel]
            if lo.max() - lo.min() > 180: lo = np.where(lo < 0, lo + 360, lo)      # anel que cruza a linha de data
            for sh in ((0, -360) if lo.max() > 180 else (0,)):
                ImageDraw.Draw(alvo).polygon([px(a + sh, b) for a, b in zip(lo, la)], fill=0 if k else 255)
dil = lambda m, r: np.asarray(m.filter(ImageFilter.MaxFilter(r))).astype(np.float32) / 255
gelo = dil(gelo, 3); terra = dil(terra, 5) * (1 - gelo)           # um pouco de margem na costa

mx, mn = A.max(2), A.min(2); sat = mx - mn
neve = np.clip((mn - 110) / 60, 0, 1) * np.clip((45 - sat) / 25, 0, 1)
faixa = np.clip((np.abs(LAT) - 28) / 8, 0, 1) * np.clip((84 - np.abs(LAT)) / 8, 0, 1)   # 36°–76°, com transições
peso = (neve * faixa * terra)[..., None]

# cor de substituição: média local das zonas sem neve (cor real da região) misturada com uma paleta por latitude
livre = ((1 - neve) * terra)[..., None]
def blur(x, r): return np.asarray(Image.fromarray(np.clip(x, 0, 255).astype(np.uint8)).filter(ImageFilter.GaussianBlur(r))).astype(np.float32)
num = sum(blur(A[..., k] * livre[..., 0], 18)[..., None] * (k == np.arange(3)) for k in range(3))
den = blur(livre[..., 0] * 255, 18)[..., None] / 255
local = num / np.maximum(den, 1e-3)
al = np.abs(LAT)[..., None]
pal = np.stack([np.interp(al[..., 0], [40, 52, 60, 66, 74], v) for v in   # floresta mista -> taiga -> tundra, sem degraus
                ([66, 58, 42, 72, 100], [86, 82, 64, 76, 94], [42, 38, 32, 46, 62])], -1).astype(np.float32)
conf = np.clip(den * 3, 0, 1) * .35
subs = conf * local + (1 - conf) * pal
rng = np.random.default_rng(3)
ruido = np.asarray(Image.fromarray((rng.random((H // 4, W // 4)) * 255).astype(np.uint8)).resize((W, H), Image.BICUBIC)).astype(np.float32)[..., None] / 255
subs *= .85 + .3 * ruido                                         # textura para não ficar liso
B = A * (1 - peso) + subs * peso

# afinação geral: terra um pouco mais viva, oceano um pouco mais claro e azul (como as referências)
t3 = terra[..., None]; cinza = B.mean(2, keepdims=True)
B = B * (1 - t3) + (cinza + (B - cinza) * 1.18) * 1.06 * t3
mar = ((1 - terra) * (1 - gelo) * np.clip((70 - np.abs(LAT)) / 10, 0, 1))[..., None]   # sem mexer no gelo marinho
B = B * (1 - mar) + (B * [1.10, 1.18, 1.22] + [4, 10, 18]) * mar
Image.fromarray(np.clip(B, 0, 255).astype(np.uint8)).save(sai, quality=86, optimize=True, progressive=True)
print("gravado", sai, f"| neve retirada em {float((peso[..., 0] > .5).mean()):.1%} da imagem")
