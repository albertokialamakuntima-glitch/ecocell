"""Descarrega o metano mensal do AIRS (NASA) do mundo inteiro, mês a mês, e guarda só o necessário
(o nível de 400 hPa, diurno e noturno). Cada ficheiro da NASA tem ~200 MB; é descarregado, lido e apagado.

Uso (no cmd, com o token definido):
    py airs_baixar.py 2019-01 2019-12     <- TESTE com 12 meses (~2,5 GB de rede)
    py airs_baixar.py                     <- TUDO: 2002-09 até hoje (agora só ~250 KB por mês: minutos)
É retomável: se parar, corre outra vez e continua onde ficou (guarda tudo em airs_cache/).
No fim cria airs_ch4_400hPa.npz.
"""
import os, re, sys, glob
import numpy as np, earthaccess, netCDF4

NIVEL = 6  # posição do nível de 400 hPa na lista (1000, 925, 850, 700, 600, 500, 400, ...)
os.makedirs("airs_cache", exist_ok=True); os.makedirs("airs_tmp", exist_ok=True)

def extrair(path):
    ds = netCDF4.Dataset(path)
    def nivel(nome):
        v = np.ma.filled(ds.variables[nome][NIVEL].astype("float32"), np.nan)
        v[(v <= 0) | (v > 5000)] = np.nan
        return v
    A, D = nivel("CH4_VMR_A"), nivel("CH4_VMR_D")
    nf = True
    try:
        m = re.search(r"UpperLeftPointMtrs=\(([-\d.]+),([-\d.]+)\)", ds.getncattr("StructMetadata.0"))
        nf = float(m.group(2)) > 0          # norte em cima?
    except Exception: pass
    ds.close()
    return A, D, nf

def escalar(x):
    """O AIRS guarda o metano como fração de mistura (~1.8e-6). Converte para ppb (~1800) e limpa valores inválidos."""
    x = np.asarray(x, dtype="float64")
    v = x[np.isfinite(x) & (x > 0)]
    if v.size and np.median(v) < 1: x = x * 1e9
    return np.where((x > 0) & (x < 5000), x, np.nan).astype("float32")

def juntar():
    fs = sorted(glob.glob("airs_cache/*.npz"))
    meses = [os.path.basename(f)[:7] for f in fs]; z = [np.load(f) for f in fs]
    A = np.stack([escalar(a["A"]) for a in z]); D = np.stack([escalar(a["D"]) for a in z])
    np.savez_compressed("airs_ch4_400hPa.npz", meses=np.array(meses), A=A, D=D, norte_em_cima=bool(z[0]["nf"]))
    print(f"\nGravado airs_ch4_400hPa.npz com {len(meses)} meses ({meses[0]} a {meses[-1]}).")
    print("Média global do 1.º mês, diurno (deve rondar 1750-1850 ppb):", round(float(np.nanmean(A[0])), 1))
    print("Média global do último mês, diurno:", round(float(np.nanmean(A[-1])), 1), "| células com dados no 1.º mês:", f"{np.isfinite(A[0]).mean():.0%}")

if __name__ == "__main__":
    ini = (sys.argv[1] + "-01") if len(sys.argv) > 2 else "2002-09-01"
    fim = (sys.argv[2] + "-28") if len(sys.argv) > 2 else "2030-01-01"
    if os.environ.get("EARTHDATA_TOKEN"): earthaccess.login(strategy="environment")
    else: earthaccess.login(strategy="interactive", persist=True)
    res = earthaccess.search_data(short_name="AIRS3STM", version="7.0", temporal=(ini, fim), count=-1)
    print("ficheiros mensais encontrados:", len(res))
    import time
    def extrair_opendap(g, sessao):
        """Lê SÓ o nível de 400 hPa pelo OPeNDAP (poucos KB em vez de 200 MB)."""
        import xarray as xr
        url = next(r["URL"] for r in g["umm"]["RelatedUrls"] if "OPENDAP" in str(r.get("Subtype", "")).upper())
        ds = xr.open_dataset("dap2://" + url.split("://")[1], engine="pydap", session=sessao)
        res_ = {}
        for suf, grp in (("A", "ascending"), ("D", "descending")):
            v = ds[f"/{grp}/Data_Fields/CH4_VMR_{suf}"]
            lev = np.asarray(ds[f"StdPressureLev_{grp}"].values, dtype="float64")
            k = int(np.argmin(np.abs(lev - 400.0)))
            x = np.asarray(v.isel({v.dims[0]: k}).values, dtype="float32")
            x[(x <= 0) | (x > 5000)] = np.nan
            res_[suf] = x; lat = np.asarray(ds[f"YDim_{grp}"].values)
        ds.close()
        return res_["A"], res_["D"], bool(lat[0] > lat[-1]), float(lev[k])

    def processar(i, g):
        link = g.data_links()[0]; m = re.search(r"AIRS\.(\d{4})\.(\d{2})\.\d{2}\.L3", link)
        if not m: return True
        out = f"airs_cache/{m.group(1)}-{m.group(2)}.npz"
        if os.path.exists(out): return True
        for tentativa in range(6):
            try:
                try:
                    A, D, nf, lv = extrair_opendap(g, earthaccess.get_requests_https_session()); via = f"opendap, nível {lv:g} hPa"
                except Exception as e1:
                    print(f"  OPeNDAP falhou ({str(e1)[:90]}); a tentar o ficheiro inteiro...", flush=True)
                    f = str(earthaccess.download([g], "airs_tmp")[0]); A, D, nf = extrair(f); os.remove(f); via = "ficheiro inteiro"
                np.savez_compressed(out, A=A, D=D, nf=nf)
                print(f"[{i}/{len(res)}] {m.group(1)}-{m.group(2)} ok ({via})", flush=True); return True
            except Exception as e:
                espera = 20 * (tentativa + 1)
                print(f"  falhou {m.group(1)}-{m.group(2)} ({tentativa+1}/6), espero {espera}s:", str(e)[:80], flush=True)
                for x in glob.glob("airs_tmp/*"):
                    try: os.remove(x)
                    except Exception: pass
                time.sleep(espera)
                try:
                    earthaccess.login(strategy="environment") if os.environ.get("EARTHDATA_TOKEN") else earthaccess.login(strategy="interactive", persist=True)
                except Exception: pass
        return False
    falhados = [g for i, g in enumerate(res, 1) if not processar(i, g)]
    if falhados:
        print(f"\n{len(falhados)} mês(es) falharam. Espero 2 minutos e tento outra vez...")
        time.sleep(120)
        falhados = [g for g in falhados if not processar(0, g)]
        if falhados: print("Ainda faltam meses. Corre o mesmo comando outra vez mais tarde: ele continua de onde ficou.")
    juntar()
