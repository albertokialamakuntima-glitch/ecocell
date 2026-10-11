# EcoCell Earth Intelligence: reprodução dos números

Este diretório contém os scripts que geram todos os resultados do mapa. Cada passo produz um ficheiro que o passo seguinte usa.

## Dados usados
- **Metano (NASA Aqua AIRS)**: produto AIRS3STM v7.0 (nível 3 mensal, 1°), obtido no GES DISC, nível de 400 hPa, 2002-09 a 2026-08.
- **Chuva (NASA/JAXA GPM)**: IMERG mensal V07 (DOI 10.5067/GPM/IMERG/3B-MONTH/07), totais anuais 2003–2024, lido pelo Google Earth Engine.
- **Verificação externa (não são variáveis do desafio)**: ESA Sentinel-5P (via Earth Engine) e NOAA GML (média global mensal de metano).

## Ordem de execução
1. `airs_baixar.py` descarrega o AIRS e cria `airs_ch4_400hPa.npz` (precisa de um token do Earthdata em `EARTHDATA_TOKEN`).
2. `airs_tendencias_mundo.py` (precisa de `world.json` e `continentes.json`) cria `airs_trends.json`: tendência de Sen, Mann-Kendall com correção de autocorrelação e FDR, por célula, país e continente.
3. `airs_periodos.py` cria `periodos.json`: tendência por período, aceleração (bootstrap) e faixas de latitude.
4. `gpm_exportar.py` cria `gpm_anual_1deg.csv` (precisa do Earth Engine); `gpm_analise.py` cria `gpm_trends.json`.
5. `s5p_tendencia_2019.py` cria `s5p_slope_2019_2026.csv`; `comparar_s5p_airs.py` cria `validacao.json`.
6. `noaa_vs_airs.py` (precisa de `ch4_mm_gl.csv`, descarregado em https://gml.noaa.gov/ccgg/trends_ch4/) cria `noaa_vs_airs.json`.
7. `build_site.py` junta tudo e gera `index.html`. `daily_update.py` atualiza só a camada de apoio do Sentinel-5P.

## Requisitos
Python 3.12, `numpy`, `scipy`, `pandas`, `shapely>=2`, `matplotlib`, `earthaccess`, `netCDF4`, `xarray`, `pydap`, `earthengine-api`.

## Limites
O AIRS capta cerca de 45% da amplitude do crescimento medido à superfície; as diferenças regionais do AIRS não foram confirmadas pelo Sentinel-5P; no IMERG, oceanos e latitudes altas são menos fiáveis. Ver a janela "Métodos e fontes" do mapa.

## Como citar
Ver `CITATION.cff` na raiz do repositório.
