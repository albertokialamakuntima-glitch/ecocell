# EcoCell Earth Intelligence — atualização diária

O site fica em `docs/index.html` e é reconstruído todos os dias com o CH4 mais recente do Sentinel-5P.

## Montar (uma vez)
1. Cria um repositório no GitHub e sobe esta pasta inteira.
2. Google Cloud: cria uma *service account*, ativa a Earth Engine API, gera uma chave JSON e regista a conta em https://code.earthengine.google.com/register (projeto não comercial).
3. No GitHub, em Settings > Secrets > Actions, cria `EE_SERVICE_ACCOUNT_JSON` (o conteúdo do JSON) e `EE_PROJECT` (o teu Project ID).
4. Settings > Pages: fonte "Deploy from a branch", pasta `/docs`.
5. Actions > "Atualização diária" > Run workflow, para testar já.

## Horário
O GitHub corre às 23:00 UTC, que são 00:00 em Luanda. Pode atrasar uns minutos. Os dados do Sentinel-5P (OFFL) chegam ao Earth Engine com alguns dias de atraso, por isso "atualizar todos os dias" traz a última data disponível.

## Alternativa só no teu PC (Windows)
`schtasks /create /tn EcoCell /tr "py C:\caminho\daily_update.py" /sc daily /st 00:00` e depois `py build_site.py`. O PC tem de estar ligado.
