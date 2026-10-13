"""Motor de prioridade: pontuação explicável (0 a 100) para ordenar zonas a investigar.

NÃO é uma probabilidade nem prova de emissão: é uma média ponderada de quatro critérios entre 0 e 1.
Os pesos são configuráveis. Equivalente em Python da pontuação do mapa (verificado por testes contra os valores do JavaScript).
Critérios: anomalia (percentil entre todas as células), persistência (percentil da tendência, só se positiva, fiável e |lat| <= 60),
qualidade dos dados (tendência fiável 0,5; |lat| <= 60 0,3; terra 0,2 ou mar 0,1) e proximidade a instalações conhecidas (1 - d/300 km).
"""
import math
import numpy as np

REF = [("Refinaria de Luanda", -8.78, 13.3), ("Refinaria de Cabinda (Malembo)", -5.33, 12.18), ("Ras Tanura", 26.64, 50.16),
       ("Jamnagar", 22.35, 69.83), ("Roterdão", 51.89, 4.32), ("Port Arthur", 29.88, -93.94), ("Ulsan", 35.49, 129.33)]
PESOS_DEFAULT = {"an": 40.0, "pe": 25.0, "ql": 20.0, "px": 15.0}
RAIO_KM = 300.0

def haversine_km(a, b, c, d):
    r = math.pi / 180
    x = math.sin((c - a) * r / 2) ** 2 + math.cos(a * r) * math.cos(c * r) * math.sin((d - b) * r / 2) ** 2
    return 12742 * math.asin(math.sqrt(x))

def dist_instalacao(lat, lon, refs=REF):
    return min(haversine_km(lat, lon, la, lo) for _, la, lo in refs)

def percentil(arr_ordenado, x):
    """Percentil por posição média: (n(<x) + n(=x)/2) / n."""
    a = np.asarray(arr_ordenado)
    if a.size == 0: return 0.0
    esq = np.searchsorted(a, x, "left"); dir_ = np.searchsorted(a, x, "right")
    return float((esq + (dir_ - esq) / 2) / a.size)

def fiavel(cell): return cell["ch"] is not None and abs(cell["lat"]) <= 60

def criterios(cell, A_ord, H_ord, refs=REF):
    dk = dist_instalacao(cell["lat"], cell["lon"], refs)
    return {"an": percentil(A_ord, cell["an"]),
            "pe": percentil(H_ord, cell["ch"]) if (fiavel(cell) and cell["ch"] > 0) else 0.0,
            "ql": (0.5 if cell["ch"] is not None else 0.0) + (0.3 if abs(cell["lat"]) <= 60 else 0.0) + (0.1 if cell["sea"] else 0.2),
            "px": max(0.0, 1.0 - dk / RAIO_KM)}, dk

def pontuacao(cell, A_ord, H_ord, pesos=None, refs=REF):
    p = pesos or PESOS_DEFAULT
    k, dk = criterios(cell, A_ord, H_ord, refs); S = sum(p.values())
    if S <= 0: raise ValueError("A soma dos pesos tem de ser positiva.")
    return 100.0 * sum(p[x] * k[x] for x in ("an", "pe", "ql", "px")) / S, k, dk

def preparar_populacao(cells):
    A = sorted(c["an"] for c in cells); H = sorted(c["ch"] for c in cells if fiavel(c))
    return A, H

def recomendacao(cell, A_ord, H_ord, pesos=None, dados_industriais=False):
    """Texto explicável. Sem dados industriais, recomenda recolha de dados em vez de estimar o potencial."""
    s, k, dk = pontuacao(cell, A_ord, H_ord, pesos)
    falta = ["caudal e composição do gás da fonte", "confirmação da fonte (inventário e fogachos)", "outro período ou instrumento para confirmar a anomalia"]
    return {"pontuacao": round(s, 1), "criterios": k, "instalacao_mais_proxima_km": round(dk), "dados_em_falta": [] if dados_industriais else falta,
            "acao": "estimar potencial com os dados fornecidos (não validados)" if dados_industriais else "recolher dados antes de estimar o potencial"}
