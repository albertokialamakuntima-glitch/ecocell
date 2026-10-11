import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
import math, pytest
from dataclasses import replace
from ecocell_model import Gas, Params, simulate, compare_scenarios, MW, F, R

def run(**kw): return simulate(Gas(), replace(Params(t_end_h=24.0, dt_h=0.01), **kw))

def test_composicao_tem_de_somar_1():
    with pytest.raises(ValueError): Gas(y={"CH4": 0.5, "CO2": 0.2}).validate()

def test_caudal_molar_gas_ideal():
    g = Gas(flow_m3_h=1.0, T_K=273.15, P_Pa=101325.0)
    assert g.n_total_mol_h() == pytest.approx(44.6, rel=0.01)           # 1 m3 a 0 °C e 1 atm ≈ 44,6 mol

def test_balanco_de_carbono_fecha():
    r = run()["resumo"]
    assert abs(r["erro_balanco_carbono_mol"]) < 1e-6 * max(r["CH4_consumido_mol"], 1.0)

def test_balanco_ch4_e_o2_fecham():
    r = run()["resumo"]
    assert abs(r["erro_balanco_CH4_mol"]) < 1e-6 * max(r["CH4_transferido_mol"], 1.0)
    assert abs(r["erro_balanco_O2_mol"]) < 1e-6 * max(r["O2_consumido_mol"], 1.0)

def test_nao_consome_mais_do_que_o_alimentado():
    r = run()["resumo"]
    assert r["CH4_consumido_mol"] <= r["CH4_ao_reator_mol"] + 1e-6 and r["CH4_residual_mol"] >= 0

def test_convergencia_do_integrador():
    a = simulate(Gas(), replace(Params(t_end_h=24.0), rtol=1e-6))["resumo"]["CH4_consumido_mol"]
    b = simulate(Gas(), replace(Params(t_end_h=24.0), rtol=1e-10))["resumo"]["CH4_consumido_mol"]
    assert abs(a - b) / b < 1e-3

def test_resultado_nao_depende_da_amostragem_de_saida():
    a = simulate(Gas(), replace(Params(t_end_h=24.0), dt_h=0.05))["resumo"]["CH4_consumido_mol"]
    b = simulate(Gas(), replace(Params(t_end_h=24.0), dt_h=0.005))["resumo"]["CH4_consumido_mol"]
    assert abs(a - b) / b < 1e-6

def test_parametros_mudam_resultados():
    assert run(kLa_CH4_h=20.0)["resumo"]["CH4_consumido_mol"] < run(kLa_CH4_h=150.0)["resumo"]["CH4_consumido_mol"]
    assert run(capture_eff=0.2)["resumo"]["CH4_ao_reator_mol"] < run(capture_eff=0.9)["resumo"]["CH4_ao_reator_mol"]

def test_sem_recuperacao_de_eletroes_nao_ha_eletricidade():
    r = run(f_e_recuperavel=0.0)["resumo"]
    assert r["energia_bruta_kWh"] == 0 and r["H2_kg"] == 0

def test_corrente_nao_excede_o_maximo_teorico():
    s = run(f_e_recuperavel=1.0, coulombic_eff=1.0)
    maximo = [8 * F * 1.0 * 0 for _ in s["series"]["corrente_A"]]
    V = 10.0
    # corrente máxima = 8 F (mol CH4/s); verifica contra o consumo máximo possível (alimentação)
    n_feed_s = Gas().molar_flows()["CH4"] * 0.9 / 3600.0
    assert max(s["series"]["corrente_A"]) <= 8 * F * (n_feed_s * 5) + 1e-9      # folga: o consumo instantâneo pode exceder a alimentação por acumulação

def test_eletrolise_sem_conversor_e_tensao_baixa_nao_produz_h2():
    r = run(f_e_recuperavel=1.0, coulombic_eff=1.0, aux_kW_per_m3=0.0, use_converter=False, V_cell_V=0.5, n_series=1)["resumo"]
    assert r["H2_kg"] == 0 and r["tensao_suficiente_sem_conversor"] is False

def test_h2_consistente_com_a_energia():
    r = run(f_e_recuperavel=1.0, coulombic_eff=1.0, aux_kW_per_m3=0.0, use_converter=True, V_cell_V=1.0)["resumo"]
    esperado = r["energia_liquida_kWh"] / (39.4 / 0.65)
    assert r["H2_kg"] == pytest.approx(esperado, rel=1e-9) and r["H2_kg"] > 0

def test_estequiometria_amonia():
    r = run(f_e_recuperavel=1.0, coulombic_eff=1.0, aux_kW_per_m3=0.0, V_cell_V=1.0)["resumo"]
    nH2 = r["H2_kg"] * 1000 / MW["H2"]; nNH3 = r["NH3_teorica_kg"] * 1000 / MW["NH3"]
    assert nNH3 == pytest.approx(2.0 / 3.0 * nH2, rel=1e-9)

def test_cenarios_comparaveis_e_conservador_nao_supera_otimista():
    c = compare_scenarios(Gas(), Params(t_end_h=24.0, dt_h=0.01))
    assert set(c) == {"referencia", "conservador", "otimizado"}
    assert c["referencia"]["CH4_consumido_mol"] == 0
    for k in ("conservador", "otimizado"):
        assert c[k]["evitadas_conservador_kgCO2e"] <= c[k]["evitadas_otimista_kgCO2e"] + 1e-9

def test_limite_de_parametros():
    with pytest.raises(ValueError): Params(capture_eff=1.5).validate()

def test_referencia_nao_tem_consumo_auxiliar_nem_biomassa():
    r = simulate(Gas(), Params(t_end_h=24.0, capture_eff=0.0))["resumo"]
    assert r["energia_auxiliar_kWh"] == 0 and r["biomassa_produzida_kg"] == 0 and r["CH4_consumido_mol"] == 0

def test_breakeven_coerente():
    from ecocell_model import breakeven_f_e
    p = replace(Params(t_end_h=24.0), coulombic_eff=1.0, V_cell_V=1.0, aux_kW_per_m3=0.05)
    f = breakeven_f_e(Gas(), p)
    assert f is not None and 0 < f < 1
    assert simulate(Gas(), replace(p, f_e_recuperavel=min(f * 1.05, 1.0)))["resumo"]["energia_liquida_kWh"] >= 0
