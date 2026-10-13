# EcoCell Engine

Núcleo científico, em Python, do EcoCell Earth Intelligence: balanços de massa e energia de um sistema conceptual de tratamento de gases com biorreator metanotrófico, bioeletroquímica, eletrólise e amónia.

**Aviso essencial.** O desempenho do EcoCell ainda não foi validado experimentalmente. Todos os parâmetros marcados como hipótese são valores de demonstração editáveis. Os resultados são simulações e nunca medições. O modelo serve para explorar cenários e para dizer que dados faltam medir.

## Executar
```
pip install -r requirements-engine.txt
python -m pytest -q          # 26 testes
python demo_cenarios.py      # grava resultados_cenarios.json e relatorio_cenarios.md
```

## Etapas
- **A) Gás:** caudal molar (gás ideal), balanço por componente, captado, encaminhado, consumido e residual.
- **B) Biorreator:** Monod duplo (CH4 e O2), manutenção limitada por substrato, decaimento, transferência gás-líquido (lei de Henry com correção de temperatura), seletividade em metanol e oxidação posterior do metanol. O balanço de carbono fecha (verificado pelos testes). A integração é rígida (LSODA): o RK4 explícito mostrou-se instável para estes parâmetros.
- **C) Bioeletroquímica:** I = 8 F x fração de eletrões recuperados x eficiência coulómbica x consumo de CH4; P = V I; energia por integração no tempo. Distingue energia bruta, auxiliar e líquida.
- **D) Hidrogénio:** consumo específico = PCS do H2 / eficiência; sem conversor, tensões abaixo de 1,23 V não alimentam o eletrolisador.
- **E) Amónia:** só o limite estequiométrico (2 mol de NH3 por 3 mol de H2). Não ocorre no biorreator.
- **F) Balanço global:** emissões evitadas em duas contabilidades: otimista (carbono dos produtos retido) e conservadora (todo o carbono consumido volta como CO2). GWP100 do CH4 = 29,8 (IPCC AR6, fóssil).

## Motor de prioridade (`recommendation_engine.py`)
Pontuação explicável de 0 a 100 para ordenar zonas a investigar, com pesos configuráveis (anomalia 40, persistência 25, qualidade dos dados 20, proximidade a fontes conhecidas 15). **Não é uma probabilidade nem prova de emissão.** Sem dados industriais da fonte, recomenda recolher dados antes de estimar o potencial. Os testes verificam que o resultado coincide com o do mapa (JavaScript) em 208 células.

## Gémeo digital
`compare_scenarios` corre referência (EcoCell desligado), conservador e otimizado. `breakeven_f_e` calcula a recuperação de eletrões necessária para a eletricidade bruta igualar o consumo auxiliar.

## Limitações
- O CO2 alimentado passa pelo reator sem ser assimilado (a assimilação de CO2 não está modelada).
- A estequiometria de O2 é um parâmetro único (mol de O2 por mol de CH4).
- Não distingue concentração atmosférica de caudal de emissão: o caudal de gás tem de vir de dados da instalação.
- Os parâmetros biológicos são ordens de grandeza, a substituir por dados experimentais.

## Motor de prioridade (recommendation_engine.py)
Duas pontuações de propósito. A **prioridade de investigação** usa só satélite (anomalia, tendência fiável, qualidade dos dados, proximidade a instalações conhecidas), com pesos editáveis e uma banda de incerteza heurística de ±5 a ±30 pontos. A **prioridade de aplicação do EcoCell** só existe se houver dados da instalação (caudal e composição do gás); sem eles devolve "dados insuficientes" e a recomendação é recolher dados. A pontuação é um índice de priorização, não uma probabilidade científica.
