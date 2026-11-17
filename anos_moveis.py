"""Comum ao metano do Sentinel-5P e ao CO2 do OCO-2/OCO-3: de médias mensais por bloco de 1° para «anos móveis».

Cada «ano» são 12 meses seguidos que acabam no último mês com dados (ex.: out/2025 a set/2026 conta como 2026), por isso
a tendência usa sempre os dados mais recentes, sem o erro de comparar anos incompletos (a estação do ano pesa muito).
Grava <id>_anual_f32.bin (float32, anos x 180 x 360, Norte em cima) e <id>_anual.json, que o variaveis_tendencias.py lê.
"""
import json, warnings
import numpy as np


def grava(vid, meses, cubo, min_meses=6, min_anos=6, fonte=""):
    """meses: ['2019-02', ...] por ordem; cubo: (meses, 180, 360) médias mensais (NaN = sem dados)"""
    cobre = np.isfinite(cubo).reshape(len(meses), -1).sum(1)
    fim = len(meses) - 1
    while fim > 0 and cobre[fim] < 0.5 * np.median(cobre[cobre > 0]): fim -= 1      # último mês com cobertura normal
    # tira a estação do ano (média de cada mês, por faixa de latitude): um bloco com poucos meses deixa de ficar enviesado
    mes = np.array([int(x[5:7]) for x in meses[: fim + 1]]); cubo = np.array(cubo[: fim + 1], dtype="float64")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        linha = np.nanmean(cubo, axis=(0, 2))                                  # média de cada latitude
        for k in range(1, 13):
            sel = mes == k
            if sel.any(): cubo[sel] -= np.nan_to_num(np.nanmean(cubo[sel], axis=(0, 2)) - linha)[None, :, None]
    blocos, rot = [], []
    a = fim
    while a - 11 >= 0:
        X = cubo[a - 11: a + 1]; n = np.isfinite(X).sum(0)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore"); m = np.nanmean(X, 0)
        blocos.append(np.where(n >= min_meses, m, np.nan)); rot.append(int(meses[a][:4]))
        a -= 12
    blocos, rot = blocos[::-1], rot[::-1]
    if len(blocos) < min_anos:
        print(f"{vid}: só {len(blocos)} anos móveis (precisa de {min_anos}): ainda não grava"); return False
    np.asarray(blocos, dtype="<f4").tofile(f"{vid}_anual_f32.bin")
    ini = meses[fim - 12 * len(blocos) + 1]
    json.dump(dict(anos=rot, forma=[len(blocos), 180, 360], desde=ini, ate=meses[fim], fonte=fonte,
                   regra=f"anos de 12 meses que acabam em {meses[fim]}; bloco conta se tiver {min_meses}+ meses"),
              open(f"{vid}_anual.json", "w", encoding="utf-8"))
    print(f"{vid}: {len(blocos)} anos móveis, de {ini} a {meses[fim]} | média do último ano {np.nanmean(blocos[-1]):.2f}")
    return True
