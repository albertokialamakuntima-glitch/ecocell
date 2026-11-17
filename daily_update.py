"""Atualização automática do EcoCell — com dados da NASA (e o metano do Sentinel-5P, da ESA), todos os dias às 00:00 de Luanda (ver daily.yml).

  1. airs_baixar.py          metano mensal do NASA Aqua AIRS (NASA Earthdata, precisa do segredo EARTHDATA_TOKEN);
                             só descarrega os meses novos (os antigos ficam na cache do GitHub Actions)
  2. airs_celulas.py         camada de células do metano (últimos 36 meses) -> cells.csv, meta.json
     airs_periodos.py        «Quando mudou?» do metano (mundo, continentes, países e blocos) -> periodos*.json
     co2_baixar.py           CO2 dos satélites OCO-2 e OCO-3 da NASA (dias novos + o passado aos poucos) -> co2_anual_*
     s5p_exportar.py         metano do Sentinel-5P (ESA Copernicus) pelo Earth Engine -> ch4s5p_anual_*
  3. variaveis_exportar.py   médias anuais das variáveis da NASA no Google Earth Engine
                             (MODIS, GPM IMERG, GRACE/GRACE-FO, FLDAS, MODIS Aqua oceano), de 2003 ao último ano completo
  4. variaveis_tendencias.py tendências, significância, países e séries por célula
Se uma fonte falhar (rede, quota), mantêm-se os dados anteriores dessa parte e o resto continua.
"""
import os, subprocess, sys

def corre(*args):
    print(f"\n=== {' '.join(args)} ===", flush=True)
    return subprocess.run([sys.executable, *args]).returncode == 0

airs_ok = False
if os.environ.get("EARTHDATA_TOKEN"):
    airs_ok = corre("airs_baixar.py") and os.path.exists("airs_ch4_400hPa.npz") and corre("airs_celulas.py")
    if airs_ok and not corre("airs_periodos.py"): print("::warning::«Quando mudou?» do metano não atualizado: mantém-se a versão anterior.")
    if not airs_ok: print("::warning::Metano AIRS não atualizado: mantém-se a versão anterior.")
else:
    print("::warning::Falta o segredo EARTHDATA_TOKEN: o metano AIRS não foi atualizado hoje.")

co2_novo = False
if os.environ.get("EARTHDATA_TOKEN"):
    co2_novo = corre("co2_baixar.py") and os.path.exists("co2_anual.json")
    if not co2_novo: print("::warning::CO2 (OCO-2/OCO-3) ainda não tem anos suficientes ou falhou: mantém-se o anterior.")

s5p_ok = corre("s5p_exportar.py") and os.path.exists("ch4s5p_anual.json")
if not s5p_ok: print("::warning::Metano do Sentinel-5P não atualizado: mantém-se o anterior.")

ee_ok = corre("variaveis_exportar.py")
if not ee_ok: print("::warning::Earth Engine falhou: as variáveis mantêm os dados anteriores.")

if (ee_ok or airs_ok or co2_novo or s5p_ok) and os.path.exists("variaveis_anuais.npz"):
    if not corre("variaveis_tendencias.py"): print("::warning::Falhou o cálculo das tendências: mantêm-se os dados anteriores.")
print("Atualização concluída.")
