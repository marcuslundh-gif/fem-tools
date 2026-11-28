"""Materialdata från materials.json (master, se docs/projektstandard.md).

    from femtools import material
    m = material("6082-T6", t=4)           # tjockleksberoende minimivärden
    m = material("PA12 SLS (EOS PA 2200)") # typvärden; E/Rm ges för Z (konservativt) om riktning saknas
"""
import json
from importlib import resources

_DATA = None


def materials():
    global _DATA
    if _DATA is None:
        _DATA = json.loads(resources.files("femtools").joinpath("materials.json").read_text(encoding="utf-8"))
    return {k: v for k, v in _DATA.items() if not k.startswith("_")}


def material(name, t=None, riktning="z"):
    """Returnerar dict med E, nu, rho, alpha, typ, kalla (+ Rp02/Rm eller Rm för polymerer).
    Metaller (rader med t_min/t_max, även TYPVARDE) kräver tjocklek t [mm]. Polymerer: riktning 'xy' eller 'z' (default z = konservativt)."""
    db = materials()
    if name not in db:
        raise KeyError(f"Material '{name}' saknas i materials.json. Finns: {', '.join(db)}. Fråga, gissa inte.")
    rows = db[name]
    if "t_min" in rows[0]:                       # metall (MINIMI eller TYPVARDE): tjockleksuppslag
        if t is None:
            raise ValueError(f"{name}: tjocklek t [mm] krävs (värdena är tjockleksberoende).")
        for r in rows:
            if r["t_min"] < t <= r["t_max"]:
                return dict(r, namn=name)
        raise ValueError(f"{name}: t = {t} mm ligger utanför tabellen ({rows[0]['t_min']}–{rows[-1]['t_max']} mm).")
    r = dict(rows[0], namn=name)
    k = "xy" if riktning == "xy" else "z"
    r["E"], r["Rm"], r["riktning"] = r[f"E_{k}"], r[f"Rm_{k}"], k
    return r
