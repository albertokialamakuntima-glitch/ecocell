"""Executa os três cenários do gémeo digital e grava resultados_cenarios.json e relatorio_cenarios.md.
TODOS os resultados são SIMULADOS com parâmetros de demonstração (hipóteses editáveis em ecocell_model.py)."""
import json
from dataclasses import replace, asdict
from ecocell_model import Gas, Params, CENARIOS, compare_scenarios, breakeven_f_e, PROVENIENCIA

g, base = Gas(), Params()
res = compare_scenarios(g, base)
eq = {n: breakeven_f_e(g, replace(base, **CENARIOS[n])) for n in ("conservador", "otimizado")}
json.dump({"aviso": "Resultados simulados com hipóteses de demonstração; não são medições.", "gas": asdict(g), "parametros_base": asdict(base),
           "cenarios_parametros": CENARIOS, "resultados": res, "f_e_de_equilibrio": eq, "proveniencia": PROVENIENCIA},
          open("resultados_cenarios.json", "w"), ensure_ascii=False, indent=1, default=float)
linhas = [("CH4 disponível (mol)", "CH4_disponivel_mol"), ("CH4 consumido (mol)", "CH4_consumido_mol"), ("Conversão de CH4 (%)", "conversao_CH4_pct"),
          ("Biomassa produzida (kg)", "biomassa_produzida_kg"), ("Metanol recuperável (kg)", "metanol_recuperavel_kg"),
          ("Energia bruta (kWh)", "energia_bruta_kWh"), ("Energia auxiliar (kWh)", "energia_auxiliar_kWh"), ("Energia líquida (kWh)", "energia_liquida_kWh"),
          ("H2 potencial (kg)", "H2_kg"), ("NH3 teórica (kg)", "NH3_teorica_kg"),
          ("Emissões evitadas, otimista (kg CO2e)", "evitadas_otimista_kgCO2e"), ("Emissões evitadas, conservador (kg CO2e)", "evitadas_conservador_kgCO2e")]
md = ["# Cenários do gémeo digital EcoCell (SIMULADO)\n", "> Hipóteses de demonstração, não medições. Horizonte de 72 h, 20 m3/h de gás com 60% de CH4.\n",
      "| Indicador | " + " | ".join(res) + " |", "|---|" + "---|" * len(res)]
for nome, k in linhas: md.append(f"| {nome} | " + " | ".join(f"{res[c][k]:,.1f}" for c in res) + " |")
md += ["", "**Fração de eletrões recuperados necessária para a eletricidade bruta igualar o consumo auxiliar:** " +
       ", ".join(f"{n}: {v*100:.0f}%" if v is not None else f"{n}: impossível" for n, v in eq.items()) + ".",
       "", "Com estas hipóteses o sistema consome mais eletricidade do que produz. O que falta medir: recuperação de eletrões, seletividade em metanol e consumo auxiliar real."]
open("relatorio_cenarios.md", "w", encoding="utf-8").write("\n".join(md)); print("\n".join(md))
