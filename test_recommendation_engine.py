import sys, os, json, pytest
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
from recommendation_engine import *

def pop(n=200):
    import random; random.seed(1)
    return [{"lat": random.uniform(-50, 50), "lon": random.uniform(-100, 100), "an": random.gauss(0, 10), "ch": random.choice([None, random.gauss(1, 1)]), "sea": random.choice([0, 1])} for _ in range(n)]

def test_pontuacao_entre_0_e_100():
    c = pop(); A, H = preparar_populacao(c)
    for x in c: assert 0 <= pontuacao(x, A, H)[0] <= 100

def test_maior_anomalia_maior_pontuacao_com_peso_so_na_anomalia():
    c = pop(); A, H = preparar_populacao(c); w = {"an": 1, "pe": 0, "ql": 0, "px": 0}
    a = max(c, key=lambda x: x["an"]); b = min(c, key=lambda x: x["an"])
    assert pontuacao(a, A, H, w)[0] > pontuacao(b, A, H, w)[0]

def test_proximidade_zero_alem_de_300_km():
    cell = {"lat": 0.0, "lon": -150.0, "an": 0.0, "ch": 1.0, "sea": 0}
    assert criterios(cell, [0.0], [1.0])[0]["px"] == 0.0

def test_proximidade_maxima_na_propria_instalacao():
    cell = {"lat": -8.78, "lon": 13.3, "an": 0.0, "ch": 1.0, "sea": 0}
    assert criterios(cell, [0.0], [1.0])[0]["px"] == pytest.approx(1.0, abs=1e-6)

def test_tendencia_nao_fiavel_nao_conta_para_persistencia():
    cell = {"lat": 10.0, "lon": 10.0, "an": 5.0, "ch": None, "sea": 0}
    assert criterios(cell, [0.0, 5.0], [1.0])[0]["pe"] == 0.0

def test_mudar_pesos_muda_a_ordem():
    c = pop(300); A, H = preparar_populacao(c)
    ordem = lambda w: sorted(range(len(c)), key=lambda i: -pontuacao(c[i], A, H, w)[0])[:10]
    assert ordem({"an": 100, "pe": 0, "ql": 0, "px": 0}) != ordem({"an": 0, "pe": 0, "ql": 0, "px": 100})

def test_pesos_invalidos():
    with pytest.raises(ValueError): pontuacao({"lat": 0, "lon": 0, "an": 0, "ch": 0, "sea": 0}, [0.0], [0.0], {"an": 0, "pe": 0, "ql": 0, "px": 0})

def test_sem_dados_industriais_recomenda_recolher_dados():
    c = pop(); A, H = preparar_populacao(c); r = recomendacao(c[0], A, H)
    assert "recolher dados" in r["acao"] and len(r["dados_em_falta"]) > 0

def test_igual_ao_javascript_do_mapa():
    d = json.load(open(os.path.join(os.path.dirname(__file__), "data", "js_scores.json")))
    A = sorted(a for a, h, la in d["populacao"]); H = sorted(h for a, h, la in d["populacao"] if h is not None and abs(la) <= 60)
    dif = 0.0
    for x in d["amostra"]:
        s, _, dk = pontuacao(x, A, H, d["pesos"]); dif = max(dif, abs(s - x["sc"]), abs(dk - x["dk"]))
    assert dif < 1e-3          # arredondamentos de 6 casas no ficheiro de referência
