# Seletor de variáveis (EcoCell)

O painel **Variáveis da Terra** (primeiro botão da barra) permite escolher 10 variáveis de missões e modelos da NASA e ver,
para cada uma, **o que muda, onde, quanto e se é significativo**, com o mesmo método estatístico.

| Variável | Fonte | Período |
|---|---|---|
| Metano (troposfera média) | NASA Aqua AIRS (400 hPa) | 2003–2025 |
| Chuva anual | NASA/JAXA GPM IMERG | 2003–2024 |
| Temperatura do ar (2 m) | Modelo NASA FLDAS | 2003–2024 |
| Temperatura da superfície (dia) | NASA Terra MODIS | 2003–2020 (*) |
| Vegetação (NDVI) | NASA Terra MODIS | 2003–2024 |
| Evapotranspiração | NASA Terra MODIS | 2003–2024 |
| Humidade do solo (0–10 cm) | Modelo NASA FLDAS | 2003–2024 |
| Água total armazenada | NASA/DLR GRACE e GRACE-FO | 2003–2024 |
| Temperatura do mar | NASA Aqua MODIS | 2003–2024 |
| Clorofila do oceano | NASA Aqua MODIS | 2003–2024 |

(*) Desde 2020 o satélite Terra deixou de corrigir a órbita e passa mais cedo, numa hora mais fresca: incluir esses anos criaria um falso arrefecimento.

**Método (igual para todas):** médias anuais por célula de 1°; tendência de Sen com IC 95 %; Mann-Kendall com a variância corrigida
pela autocorrelação; correção FDR de Benjamini-Hochberg (q < 0,05) sobre todas as células; o mesmo por país.

## Funções do painel
1. **Variável:** escolher a variável; mapa sobe/desce (cores fortes = significativo) no mapa 2D e no globo 3D; resumo global; países que mais sobem e descem.
2. **Que variáveis mudam mais?:** ranking pela % da área com mudança significativa.
3. **Ligar variáveis:** correlação entre as tendências de duas variáveis nas mesmas células e nos países.
4. **Clicar numa célula** (com uma variável escolhida): série anual e tendência dessa célula.
5. **Exportar CSV:** por país ou por célula (filtra pelo continente ou país aberto).

## Ficheiros a pôr no GitHub (na raiz do repositório)
- `template.html` e `build_site.py` (atualizados)
- `variaveis.json` (resumo que vai dentro do site)
- a pasta **`vars/`** inteira (séries de cada célula; o site só as carrega quando clicas ou exportas)
- `variaveis_exportar.py`, `variaveis_tendencias.py` (para refazer os dados; opcional no site)

No GitHub: **Add file → Upload files**, arrasta os ficheiros e a pasta `vars`, e **Commit changes**.

## Refazer os dados (no teu computador)
```
py variaveis_exportar.py        (Earth Engine: médias anuais, ~2 minutos)
py variaveis_tendencias.py      (precisa também de airs_ch4_400hPa.npz, gpm_anual_1deg.csv e world.json)
```
