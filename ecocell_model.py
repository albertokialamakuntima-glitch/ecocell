"""EcoCell Engine: modelo conceptual executável (balanços de massa e energia).

IMPORTANTE: todos os parâmetros marcados como [hipótese] são valores de DEMONSTRAÇÃO editáveis, não medições.
O desempenho do EcoCell (seletividade em metanol, recuperação de eletrões, eficiências) ainda não foi validado
experimentalmente; este modelo serve para explorar cenários e identificar que dados faltam medir.

Etapas: A) gás; B) biorreator (Monod duplo + transferência gás-líquido; integração rígida LSODA, porque o RK4 explícito é instável aqui); C) bioeletroquímica;
D) eletrólise; E) amónia (limite estequiométrico); F) balanço global e emissões.
"""
from dataclasses import dataclass, field, replace, asdict
import math
import numpy as np
from scipy.integrate import solve_ivp

R = 8.314462618          # J/(mol K)
F = 96485.33212          # C/mol
GWP100_CH4 = 29.8        # IPCC AR6, metano fóssil (kg CO2e por kg CH4)
MW = {"CH4": 16.043, "CO2": 44.009, "O2": 31.998, "N2": 28.014, "H2": 2.016, "NH3": 17.031, "CH3OH": 32.042}  # g/mol
M_C = 0.012011           # kg/mol
HHV_H2_KWH_KG = 39.4     # kWh/kg (poder calorífico superior)

PROVENIENCIA = {
    "kH_CH4_298": "literatura (constante de Henry, a validar)", "kH_O2_298": "literatura (a validar)",
    "mu_max_h": "[hipótese] ordem de grandeza de metanotróficos", "Y_XS_kg_mol": "[hipótese] rendimento típico 0,3 a 0,6 g/g",
    "kLa_CH4_h": "[hipótese] depende do reator", "sigma_metanol": "[hipótese] não medido para o EcoCell",
    "f_e_recuperavel": "[hipótese] não medido", "coulombic_eff": "[hipótese] não medido", "V_cell_V": "[hipótese] não medido",
    "eta_eletrolisador": "[hipótese] típico de eletrolisadores comerciais", "capture_eff": "[hipótese] depende da instalação",
}

@dataclass
class Gas:
    flow_m3_h: float = 20.0               # m3/h nas condições T_K e P_Pa
    T_K: float = 303.15
    P_Pa: float = 101325.0
    y: dict = field(default_factory=lambda: {"CH4": 0.60, "CO2": 0.35, "O2": 0.00, "N2": 0.05})   # frações molares
    def validate(self):
        if abs(sum(self.y.values()) - 1.0) > 1e-9: raise ValueError("A composição molar tem de somar 1.")
        if any(v < 0 for v in self.y.values()): raise ValueError("Frações molares negativas.")
        if not (250 <= self.T_K <= 400) or self.flow_m3_h < 0 or self.P_Pa <= 0: raise ValueError("T, P ou caudal fora dos limites.")
    def n_total_mol_h(self): return self.P_Pa * self.flow_m3_h / (R * self.T_K)       # gás ideal
    def molar_flows(self): return {k: v * self.n_total_mol_h() for k, v in self.y.items()}

@dataclass
class Params:
    capture_eff: float = 0.90             # fração do gás disponível que é captada
    routed_fraction: float = 1.0          # fração do captado encaminhada ao reator
    V_liq_m3: float = 10.0                # volume líquido do biorreator
    kLa_CH4_h: float = 150.0              # 1/h  coeficiente volumétrico de transferência de CH4
    kLa_O2_h: float = 200.0               # 1/h
    o2_partial_atm: float = 0.21          # O2 fornecido por ar (borbulhamento)
    kH_CH4_298: float = 1.4               # mol/(m3 atm) a 298,15 K;  d ln kH / d(1/T) = 1900 K
    kH_O2_298: float = 1.3                # mol/(m3 atm);  1700 K
    mu_max_h: float = 0.20                # 1/h
    Ks_CH4: float = 0.02                  # mol/m3
    Ko_O2: float = 0.01                   # mol/m3
    Y_XS_kg_mol: float = 0.008            # kg biomassa por mol de CH4 (crescimento)
    m_s_mol_kg_h: float = 0.5             # manutenção, mol CH4 por kg biomassa por hora
    k_d_h: float = 0.01                   # decaimento 1/h
    c_X: float = 0.5                      # kg C por kg de biomassa
    o2_per_ch4: float = 1.6               # mol O2 consumido por mol de CH4 consumido
    sigma_metanol: float = 0.02           # fração do CH4 consumido recuperável como metanol (seletividade)
    k_ox_metanol_h: float = 0.05          # oxidação posterior do metanol, 1/h
    X0_kg_m3: float = 0.05
    # bioeletroquímica
    f_e_recuperavel: float = 0.01         # fração dos 8 eletrões por CH4 que chega ao elétrodo
    coulombic_eff: float = 0.7
    V_cell_V: float = 0.35                # tensão por célula
    n_series: int = 1
    # energia
    aux_kW_per_m3: float = 0.3            # consumo auxiliar (agitação, aeração) por m3 de líquido
    use_converter: bool = True
    dcdc_eff: float = 0.85
    eta_eletrolisador: float = 0.65       # sobre o PCS do H2
    nh3_yield: float = 1.0                # fator sobre o máximo teórico (1 = limite estequiométrico)
    # emissões
    ref_flare_eff: float = 0.0            # cenário de referência: 0 = gás libertado para a atmosfera
    grid_ef_kgco2_kwh: float = 0.4
    t_end_h: float = 72.0
    dt_h: float = 0.01                    # só define a amostragem de saída
    rtol: float = 1e-9                    # tolerância relativa do integrador
    def validate(self):
        for n in ("capture_eff", "routed_fraction", "sigma_metanol", "f_e_recuperavel", "coulombic_eff", "dcdc_eff", "eta_eletrolisador", "ref_flare_eff"):
            if not (0 <= getattr(self, n) <= 1): raise ValueError(f"{n} tem de estar entre 0 e 1.")
        if self.V_liq_m3 <= 0 or self.dt_h <= 0 or self.t_end_h <= 0: raise ValueError("Volume, dt e tempo têm de ser positivos.")
        if self.Y_XS_kg_mol * self.c_X / M_C > 1: raise ValueError("Rendimento incompatível com o balanço de carbono (>1 C por CH4).")

def henry(k298, B, T): return k298 * math.exp(B * (1.0 / T - 1.0 / 298.15))

def _rates(state, p, g, n_feed_reactor):
    S, O, X, M = state[:4]
    Sst = henry(p.kH_CH4_298, 1900.0, g.T_K) * g.y["CH4"] * (g.P_Pa / 101325.0)
    Ost = henry(p.kH_O2_298, 1700.0, g.T_K) * p.o2_partial_atm
    S = max(S, 0.0); O = max(O, 0.0)
    mu = p.mu_max_h * S / (p.Ks_CH4 + S) * O / (p.Ko_O2 + O)
    r_S = mu * X / p.Y_XS_kg_mol + p.m_s_mol_kg_h * X * S / (p.Ks_CH4 + S)  # mol CH4 /m3/h (manutenção também limitada por substrato)
    r_CX = mu * X * p.c_X / M_C                                            # carbono para biomassa, mol/m3/h
    r_CX = min(r_CX, r_S)
    r_M = max(0.0, min(p.sigma_metanol * r_S, r_S - r_CX))
    decay_C = p.k_d_h * X * p.c_X / M_C
    co2 = r_S - r_CX - r_M + p.k_ox_metanol_h * M + decay_C
    tr = max(0.0, p.kLa_CH4_h * (Sst - S) * p.V_liq_m3)                    # mol/h, limitado pela alimentação
    tr = min(tr, n_feed_reactor)
    tr_O2 = p.kLa_O2_h * (Ost - O)
    d = np.empty(9)
    d[0] = tr / p.V_liq_m3 - r_S
    d[1] = tr_O2 - p.o2_per_ch4 * r_S
    d[2] = (mu - p.k_d_h) * X
    d[3] = r_M - p.k_ox_metanol_h * M
    d[4] = r_S * p.V_liq_m3                  # CH4 consumido acumulado (mol)
    d[5] = p.o2_per_ch4 * r_S * p.V_liq_m3   # O2 consumido acumulado (mol)
    d[6] = co2 * p.V_liq_m3                  # CO2 produzido acumulado (mol)
    d[7] = tr                                # CH4 transferido acumulado (mol)
    d[8] = tr_O2 * p.V_liq_m3                # O2 transferido acumulado (mol)
    return d, r_S, Sst, Ost

def simulate(g: Gas, p: Params):
    """Integra o biorreator por Runge-Kutta de 4.ª ordem e devolve séries e balanço global."""
    g.validate(); p.validate()
    if p.capture_eff * p.routed_fraction == 0:          # EcoCell desligado (cenário de referência): sem reator, sem consumo auxiliar
        p = replace(p, X0_kg_m3=0.0, aux_kW_per_m3=0.0)
    flows = g.molar_flows()
    n_feed = flows["CH4"] * p.capture_eff * p.routed_fraction              # mol/h de CH4 que chega ao reator
    _, _, Sst, Ost = _rates(np.array([0, 0, p.X0_kg_m3, 0, 0, 0, 0, 0, 0.0]), p, g, n_feed)
    y = np.array([0.0, Ost, p.X0_kg_m3, 0.0, 0, 0, 0, 0, 0.0])     # CH4 dissolvido inicial = 0: sem alimentação não há consumo
    npts = max(50, min(2000, int(round(p.t_end_h / p.dt_h))))
    sol = solve_ivp(lambda t, yy: _rates(yy, p, g, n_feed)[0], (0.0, p.t_end_h), y, method="LSODA", rtol=p.rtol, atol=1e-12,
                    t_eval=np.linspace(0.0, p.t_end_h, npts + 1))
    if not sol.success: raise RuntimeError("Integração falhou: " + sol.message)
    T = sol.t; Ys = sol.y.T; rS = np.array([_rates(row, p, g, n_feed)[1] for row in Ys])
    V = p.V_liq_m3
    # --- C) bioeletroquímica
    I_A = 8 * F * p.f_e_recuperavel * p.coulombic_eff * rS * V / 3600.0           # corrente (A) = e- recuperados
    V_stack = p.V_cell_V * p.n_series
    P_W = V_stack * I_A
    E_gross_kWh = float(np.trapezoid(P_W, T) / 1000.0)
    E_aux_kWh = p.aux_kW_per_m3 * V * p.t_end_h
    E_conv = E_gross_kWh * (p.dcdc_eff if p.use_converter else 1.0)
    grid_kWh = max(E_aux_kWh - E_conv, 0.0); E_avail = max(E_conv - E_aux_kWh, 0.0)
    # --- D) hidrogénio
    volt_ok = (V_stack >= 1.23) or p.use_converter
    spec = HHV_H2_KWH_KG / p.eta_eletrolisador
    H2_kg = E_avail / spec if volt_ok else 0.0
    O2_el_kg = H2_kg * (MW["O2"] / (2 * MW["H2"]))
    # --- E) amónia (limite teórico)
    nH2 = H2_kg * 1000 / MW["H2"]; nNH3 = (2.0 / 3.0) * nH2 * p.nh3_yield
    # --- F) balanço global
    H = p.t_end_h
    avail = {k: v * H for k, v in flows.items()}
    ch4_cons, o2_cons, co2_prod, ch4_tr, o2_tr = Ys[-1, 4], Ys[-1, 5], Ys[-1, 6], Ys[-1, 7], Ys[-1, 8]
    X_end, M_end = Ys[-1, 2], Ys[-1, 3]
    bio_kg = (X_end - p.X0_kg_m3) * V; meoh_kg = M_end * V * MW["CH3OH"] / 1000.0
    ch4_resid = avail["CH4"] - ch4_cons
    ref_emit = avail["CH4"] * (1 - p.ref_flare_eff)
    kg = lambda mol, s: mol * MW[s] / 1000.0
    co2_flare_ref = avail["CH4"] * p.ref_flare_eff
    ref_co2e = kg(ref_emit, "CH4") * GWP100_CH4 + kg(co2_flare_ref, "CO2")
    # EcoCell: CH4 residual + CO2 respirado (otimista: carbono de produtos fica retido) ou todo o carbono consumido (conservador)
    eco_opt = kg(ch4_resid, "CH4") * GWP100_CH4 + kg(co2_prod, "CO2") + grid_kWh * p.grid_ef_kgco2_kwh
    carbon_products = kg((X_end - p.X0_kg_m3) * V * p.c_X / M_C + M_end * V, "CO2")
    eco_cons = eco_opt + carbon_products
    resumo = {
        "CH4_disponivel_mol": avail["CH4"], "CH4_captado_mol": avail["CH4"] * p.capture_eff, "CH4_ao_reator_mol": n_feed * H,
        "CH4_transferido_mol": ch4_tr, "CH4_consumido_mol": ch4_cons, "CH4_residual_mol": ch4_resid,
        "CO2_alimentado_mol": avail["CO2"], "CO2_produzido_respiracao_mol": co2_prod, "O2_consumido_mol": o2_cons,
        "biomassa_produzida_kg": bio_kg, "metanol_recuperavel_kg": meoh_kg,
        "energia_bruta_kWh": E_gross_kWh, "energia_auxiliar_kWh": E_aux_kWh, "energia_liquida_kWh": E_conv - E_aux_kWh,
        "eletricidade_da_rede_kWh": grid_kWh, "tensao_pilha_V": V_stack, "tensao_suficiente_sem_conversor": V_stack >= 1.23,
        "H2_kg": H2_kg, "O2_eletrolise_kg": O2_el_kg, "NH3_teorica_kg": nNH3 * MW["NH3"] / 1000.0,
        "emissoes_referencia_kgCO2e": ref_co2e, "emissoes_ecocell_otimista_kgCO2e": eco_opt, "emissoes_ecocell_conservador_kgCO2e": eco_cons,
        "evitadas_otimista_kgCO2e": ref_co2e - eco_opt, "evitadas_conservador_kgCO2e": ref_co2e - eco_cons,
        "conversao_CH4_pct": 100 * ch4_cons / avail["CH4"] if avail["CH4"] else 0.0,
    }
    # verificações de balanço (valores devem ser ~0)
    c_bio = (X_end - p.X0_kg_m3) * V * p.c_X / M_C
    resumo["erro_balanco_carbono_mol"] = ch4_cons - (c_bio + M_end * V + co2_prod)
    resumo["erro_balanco_CH4_mol"] = ch4_tr - (ch4_cons + (Ys[-1, 0] - Ys[0, 0]) * V)
    resumo["erro_balanco_O2_mol"] = o2_tr - (o2_cons + (Ys[-1, 1] - Ys[0, 1]) * V)
    series = {"t_h": T, "S_CH4_mol_m3": Ys[:, 0], "O2_mol_m3": Ys[:, 1], "X_kg_m3": Ys[:, 2], "metanol_mol_m3": Ys[:, 3],
              "corrente_A": I_A, "potencia_W": P_W}
    return {"resumo": resumo, "series": series}

CENARIOS = {
    "referencia": {"capture_eff": 0.0, "f_e_recuperavel": 0.0},
    "conservador": {"capture_eff": 0.85, "kLa_CH4_h": 100.0, "kLa_O2_h": 200.0, "f_e_recuperavel": 0.01, "sigma_metanol": 0.02, "coulombic_eff": 0.6},
    "otimizado": {"capture_eff": 0.95, "kLa_CH4_h": 200.0, "kLa_O2_h": 400.0, "f_e_recuperavel": 0.05, "sigma_metanol": 0.10, "coulombic_eff": 0.8},
}

def compare_scenarios(g: Gas, base: Params, cenarios=CENARIOS):
    """Gémeo digital: executa os cenários (hipóteses claramente identificadas) e devolve o resumo de cada um."""
    return {nome: simulate(g, replace(base, **alt))["resumo"] for nome, alt in cenarios.items()}


def breakeven_f_e(g: Gas, p: Params, tol=1e-4):
    """Fração de eletrões recuperados necessária para a eletricidade bruta (após conversor) igualar o consumo auxiliar.
    Devolve None se for preciso mais de 100% (impossível). É a pergunta de engenharia mais útil: que medição falta?"""
    def net(f): return simulate(g, replace(p, f_e_recuperavel=f))["resumo"]["energia_liquida_kWh"]
    if net(1.0) < 0: return None
    lo, hi = 0.0, 1.0
    while hi - lo > tol:
        mid = (lo + hi) / 2
        if net(mid) >= 0: hi = mid
        else: lo = mid
    return hi
