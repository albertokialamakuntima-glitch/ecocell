# --- arranque: instala sozinho o que faltar (funciona mesmo que o requirements.txt ou o workflow estejam errados) ---
import sys, subprocess, importlib
def _garantir(modulo, pacote):
    try:
        importlib.import_module(modulo)
    except ImportError:
        print(f"[arranque] a instalar {pacote} ...", flush=True)
        subprocess.check_call([sys.executable, "-m", "pip", "install", "--quiet", pacote])
for _m, _p in (("ee", "earthengine-api"), ("numpy", "numpy"), ("pandas", "pandas")):
    _garantir(_m, _p)
# --- fim do arranque ---
"""Atualiza cells.csv com o CH4 mais recente do Sentinel-5P (Google Earth Engine).
O cálculo pesado corre como TAREFAS em segundo plano (exportação para assets do Earth Engine),
sem limite de 5 minutos. Depois só se descarregam os resultados já calculados.
Credenciais: EE_SERVICE_ACCOUNT_JSON e EE_PROJECT (GitHub Actions) ou o teu login local."""
import ee, os, io, json, time, urllib.request
import pandas as pd
from datetime import date, timedelta

PROJECT = os.environ.get("EE_PROJECT", "banded-lexicon-509920-b4")
if os.environ.get("EE_SERVICE_ACCOUNT_JSON"):
    k = json.loads(os.environ["EE_SERVICE_ACCOUNT_JSON"])
    ee.Initialize(ee.ServiceAccountCredentials(k["client_email"], key_data=os.environ["EE_SERVICE_ACCOUNT_JSON"]), project=PROJECT)
else:
    ee.Initialize(project=PROJECT)

# --- controlo de carga: o Earth Engine tem uma quota mensal gratuita; ao exceder entra em modo restrito (fila lenta) ---
FORCAR = os.environ.get("GITHUB_EVENT_NAME") in (None, "", "workflow_dispatch")      # execução manual: ignora os intervalos
DIAS_ENTRE_ATUALIZACOES = int(os.environ.get("EE_DIAS_ATUALIZACAO", 7))   # o Sentinel-5P muda pouco de um dia para o outro
DIAS_ENTRE_TENTATIVAS = int(os.environ.get("EE_DIAS_TENTATIVA", 3))       # se falhar, espera antes de tentar de novo
ESPERA_MAX_MIN = float(os.environ.get("EE_ESPERA_MAX_MIN", 70))           # o workflow tem limite de 90 minutos
SEM_ARRANQUE_MIN = float(os.environ.get("EE_SEM_ARRANQUE_MIN", 15))       # se nenhuma tarefa arrancar neste tempo, desiste
PAUSA_S = float(os.environ.get("EE_PAUSA_S", 60))

def _ler(p):
    try: return json.load(open(p, encoding="utf-8"))
    except Exception: return None
def _aviso(msg): print(f"::warning::{msg}", flush=True)
hoje = date.today()
def _fim(res): json.dump({"ultima_tentativa": str(hoje), "resultado": res}, open("estado_ee.json", "w"))
def _cancelar_pendentes():
    """Cancela as tarefas deste programa (descrição t_...) que ainda estejam na fila ou a correr."""
    n = 0
    try:
        for t in ee.batch.Task.list():
            st = t.status()
            if st.get("state") in ("READY", "RUNNING") and str(st.get("description", "")).startswith("t_"):
                try: t.cancel(); n += 1
                except Exception: pass
    except Exception as e:
        print("não foi possível listar as tarefas antigas:", e)
    return n
print("tarefas antigas canceladas:", _cancelar_pendentes(), flush=True)
_meta, _est = _ler("meta.json"), _ler("estado_ee.json")
if not FORCAR:
    if _meta and (hoje - date.fromisoformat(_meta["end"])).days < DIAS_ENTRE_ATUALIZACOES:
        print(f"Dados com menos de {DIAS_ENTRE_ATUALIZACOES} dias (até {_meta['end']}): nada a fazer."); sys.exit(0)
    if _est and (hoje - date.fromisoformat(_est["ultima_tentativa"])).days < DIAS_ENTRE_TENTATIVAS:
        print(f"Última tentativa em {_est['ultima_tentativa']} ({_est.get('resultado', '?')}): a aguardar."); sys.exit(0)
_fim("em curso")


END = date.today(); START = END - timedelta(days=540); RECENT = END - timedelta(days=60)
BAND = "CH4_column_volume_mixing_ratio_dry_air"
col = ee.ImageCollection("COPERNICUS/S5P/OFFL/L3_CH4").select(BAND)
base = col.filterDate(str(START), str(END))
t0 = ee.Date(str(START))
fit = base.map(lambda i: i.addBands(ee.Image(i.date().difference(t0, "day")).float().rename("t"))).select(["t", BAND]).reduce(ee.Reducer.linearFit())
stack = (base.mean().rename("mean").addBands(fit.select("scale").rename("slope"))
         .addBands(col.filterDate(str(RECENT), str(END)).mean().rename("recent")))

ROOT = f"projects/{PROJECT}/assets/ecocell"
try: ee.data.getAsset(ROOT)
except Exception: ee.data.createAsset({"type": "FOLDER"}, ROOT)

def grid(w, e, s, n):
    def c_lon(x):
        x = ee.Number(x)
        return ee.List.sequence(s, n - 1).map(lambda la: ee.Feature(
            ee.Geometry.Rectangle([x, ee.Number(la).subtract(.5), x.add(1), ee.Number(la).add(.5)], None, False),
            {"lon": x.add(.5), "lat": la}))
    return ee.FeatureCollection(ee.List.sequence(w, e - 1).map(c_lon).flatten())

W = 40
jobs = [(w, s, min(s + 70, 80)) for w in range(-180, 180, W) for s in (-60, 10)]
tasks = []
for w, s, n in jobs:
    name = f"t_{w}_{s}".replace("-", "m"); aid = f"{ROOT}/{name}"
    try: ee.data.deleteAsset(aid)
    except Exception: pass
    res = stack.reduceRegions(grid(w, w + W, s, n), ee.Reducer.mean(), scale=7000, tileScale=4).filter(ee.Filter.notNull(["mean", "slope"]))
    t = ee.batch.Export.table.toAsset(collection=res, description=name, assetId=aid)
    t.start(); tasks.append((w, s, n, aid, t)); print("enviada", name, flush=True)

t_ini = time.time(); deadline = t_ini + ESPERA_MAX_MIN * 60; arrancou = False
while time.time() < deadline:
    st = [t.status()["state"] for *_, t in tasks]
    print({x: st.count(x) for x in set(st)}, flush=True)
    if any(x != "READY" for x in st): arrancou = True
    if all(x in ("COMPLETED", "FAILED", "CANCELLED") for x in st): break
    if not arrancou and time.time() - t_ini > SEM_ARRANQUE_MIN * 60:
        _aviso("O Earth Engine não arrancou nenhuma tarefa (provável modo restrito por quota). Tarefas canceladas; os dados anteriores mantêm-se.")
        _cancelar_pendentes(); _fim("sem arranque (fila bloqueada)"); sys.exit(0)
    time.sleep(PAUSA_S)
_cancelar_pendentes()      # o que não terminou a tempo não fica na fila a gastar quota

def read(aid):
    url = ee.FeatureCollection(aid).getDownloadURL("CSV", ["lat", "lon", "mean", "slope", "recent"])
    for n in range(4):
        try: return pd.read_csv(io.StringIO(urllib.request.urlopen(url, timeout=300).read().decode()))
        except Exception as e: time.sleep(10 * (n + 1)); err = e
    print("falhou a leitura", aid, err); return None

df = pd.read_csv("cells.csv") if os.path.exists("cells.csv") else pd.DataFrame(columns=["lat", "lon", "mean", "slope", "recent"])
ok = 0
for w, s, n, aid, t in tasks:
    s_ = t.status()
    if s_["state"] != "COMPLETED": print("tarefa sem sucesso", aid, s_["state"], s_.get("error_message", "")); continue
    part = read(aid)
    if part is None or part.empty: continue
    m = (df.lon >= w) & (df.lon < w + W) & (df.lat >= s) & (df.lat < n)
    df = pd.concat([df[~m], part[["lat", "lon", "mean", "slope", "recent"]]]); ok += 1
print(f"{ok}/{len(tasks)} blocos atualizados")
if ok == 0 or len(df) < 20000:
    _aviso("Nenhum bloco foi atualizado: a manter os dados anteriores."); _fim("sem blocos atualizados"); sys.exit(0)
df.drop_duplicates(["lat", "lon"]).round(5).to_csv("cells.csv", index=False)
json.dump({"start": str(START), "end": str(END), "blocos_ok": f"{ok}/{len(tasks)}"}, open("meta.json", "w"))
_fim("ok"); print("ok", len(df), "células", START, END)
