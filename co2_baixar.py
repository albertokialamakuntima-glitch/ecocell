"""Dióxido de carbono (CO2) da NASA: satélites OCO-2 e OCO-3, ficheiros diários L2 Lite (v11), desde 2019.
Para cada dia lê SÓ latitude, longitude, XCO2 e a bandeira de qualidade (pelo OPeNDAP: ~2 MB em vez de 60-90 MB), junta as
medições boas (bandeira 0) em médias mensais por bloco de 1° e grava os «anos móveis» (anos_moveis.py):
  co2_anual_f32.bin + co2_anual.json
Corre todos os dias (daily_update.py): começa pelos dias mais recentes e, com o tempo que sobra, vai preenchendo o passado
(nos meses com mais de 13 meses só lê 1 dia em cada 3, que chega para a média do mês). Guarda o progresso em airs_cache/oco/,
a pasta que o GitHub Actions já guarda entre execuções, por isso cada noite continua onde a anterior ficou.

Uso (com o token definido): py co2_baixar.py        Precisa: numpy earthaccess xarray pydap netCDF4
Variável opcional: CO2_MINUTOS (tempo máximo por execução; 40 por omissão).
"""
import os, re, glob, json, time, datetime
import numpy as np, earthaccess
from anos_moveis import grava

PASTA = "airs_cache/oco"; os.makedirs(PASTA, exist_ok=True); os.makedirs("airs_tmp/oco", exist_ok=True)
FEITOS = f"{PASTA}/feitos.txt"
LIMITE = float(os.environ.get("CO2_MINUTOS", "40")) * 60
INICIO = datetime.date(2019, 1, 1)
PRODUTOS = [("OCO2_L2_Lite_FP", "11.3r"), ("OCO2_L2_Lite_FP", "11.2r"), ("OCO3_L2_Lite_FP", "11r")]
t0 = time.time()


def acumula(mes, la, lo, x):
    """soma e contagem do mês, por bloco de 1° (Norte em cima)"""
    f = f"{PASTA}/{mes}.npz"
    S, N = (np.load(f)["S"], np.load(f)["N"]) if os.path.exists(f) else (np.zeros((180, 360)), np.zeros((180, 360)))
    r = np.clip(np.floor(90 - la), 0, 179).astype(int); c = np.clip(np.floor(((lo + 180) % 360)), 0, 359).astype(int)
    np.add.at(S, (r, c), x); np.add.at(N, (r, c), 1)
    np.savez_compressed(f, S=S, N=N)


def le_opendap(g, sessao):
    import xarray as xr
    url = next(u["URL"] for u in g["umm"]["RelatedUrls"] if "OPENDAP" in str(u.get("Subtype", "")).upper() or "opendap" in u["URL"])
    pref = "dap4://" if "opendap.earthdata.nasa.gov" in url else "dap2://"
    ds = xr.open_dataset(pref + url.split("://", 1)[1], engine="pydap", session=sessao)
    v = {k: np.asarray(ds[k].values) for k in ("latitude", "longitude", "xco2", "xco2_quality_flag")}
    ds.close(); return v


def le_ficheiro(g):
    import netCDF4
    f = str(earthaccess.download([g], "airs_tmp/oco")[0])
    with netCDF4.Dataset(f) as ds: v = {k: np.asarray(ds.variables[k][:]) for k in ("latitude", "longitude", "xco2", "xco2_quality_flag")}
    os.remove(f); return v


if __name__ == "__main__":
    if os.environ.get("EARTHDATA_TOKEN"): earthaccess.login(strategy="environment")
    else: earthaccess.login(strategy="interactive", persist=True)
    feitos = set(open(FEITOS).read().split()) if os.path.exists(FEITOS) else set()
    hoje = datetime.date.today(); recente = hoje - datetime.timedelta(days=400)
    gr, vistos = [], set()
    for sn, ver in PRODUTOS:                                   # OCO-2 na versão mais nova que exista para cada dia; OCO-3
        try: res = earthaccess.search_data(short_name=sn, version=ver, temporal=(str(INICIO), str(hoje)), count=-1)
        except Exception as e: print(sn, ver, "falhou a pesquisa:", str(e)[:100]); continue
        for g in res:
            nome = g.data_links()[0].rsplit("/", 1)[-1]; m = re.search(r"(oco\d)_LtCO2_(\d{6})", nome)
            if not m: continue
            chave = f"{m.group(1)}_{m.group(2)}"
            if chave in vistos: continue
            d = datetime.datetime.strptime(m.group(2), "%y%m%d").date()
            if d < recente and d.toordinal() % 3: continue      # passado: 1 dia em cada 3
            vistos.add(chave); gr.append((d, chave, g))
    gr.sort(key=lambda x: x[0], reverse=True)                  # mais recentes primeiro
    falta = [x for x in gr if x[1] not in feitos]
    print(f"OCO-2/OCO-3: {len(gr)} dias a usar, {len(falta)} por ler", flush=True)
    sessao = earthaccess.get_requests_https_session(); lidos = 0
    for d, chave, g in falta:
        if time.time() - t0 > LIMITE: print("tempo esgotado nesta execução: continua amanhã"); break
        try:
            try: v = le_opendap(g, sessao)
            except Exception as e1: print(f"  OPeNDAP falhou em {chave} ({str(e1)[:80]}); ficheiro inteiro...", flush=True); v = le_ficheiro(g)
            ok = (v["xco2_quality_flag"] == 0) & np.isfinite(v["xco2"]) & (v["xco2"] > 350) & (v["xco2"] < 480)
            acumula(d.strftime("%Y-%m"), v["latitude"][ok].astype(float), v["longitude"][ok].astype(float), v["xco2"][ok].astype(float))
            with open(FEITOS, "a") as fh: fh.write(chave + "\n")
            lidos += 1
            if lidos % 25 == 0: print(f"  {lidos} dias lidos ({chave})", flush=True)
        except Exception as e:
            print(f"  falhou {chave}: {str(e)[:100]}", flush=True); time.sleep(5)
    print(f"lidos {lidos} dias nesta execução")
    fs = sorted(glob.glob(f"{PASTA}/????-??.npz"))
    if not fs: raise SystemExit("CO2: ainda sem meses")
    meses = [os.path.basename(f)[:7] for f in fs]
    todos = [f"{a}-{m:02d}" for a in range(int(meses[0][:4]), int(meses[-1][:4]) + 1) for m in range(1, 13)]
    todos = todos[todos.index(meses[0]): todos.index(meses[-1]) + 1]
    cubo = np.full((len(todos), 180, 360), np.nan, dtype="float32")
    for f, m in zip(fs, meses):
        z = np.load(f); cubo[todos.index(m)] = np.where(z["N"] >= 3, z["S"] / np.maximum(z["N"], 1), np.nan)
    grava("co2", todos, cubo, min_meses=3, fonte="NASA OCO-2 e OCO-3 L2 Lite v11 (GES DISC)")
