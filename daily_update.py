"""Atualização automática do EcoCell — SÓ com dados da NASA, todos os dias às 00:00 de Luanda (ver daily.yml).

  1. airs_baixar.py          metano mensal do NASA Aqua AIRS (NASA Earthdata, precisa do segredo EARTHDATA_TOKEN);
                             só descarrega os meses novos (os antigos ficam na cache do GitHub Actions)
  2. airs_celulas.py         camada de células do metano (últimos 36 meses) -> cells.csv, meta.json
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
    if not airs_ok: print("::warning::Metano AIRS não atualizado: mantém-se a versão anterior.")
else:
    print("::warning::Falta o segredo EARTHDATA_TOKEN: o metano AIRS não foi atualizado hoje.")

ee_ok = corre("variaveis_exportar.py")
if not ee_ok: print("::warning::Earth Engine falhou: as variáveis mantêm os dados anteriores.")

if (ee_ok or airs_ok) and os.path.exists("variaveis_anuais.npz"):
    if not corre("variaveis_tendencias.py"): print("::warning::Falhou o cálculo das tendências: mantêm-se os dados anteriores.")
print("Atualização concluída (só dados da NASA).")
