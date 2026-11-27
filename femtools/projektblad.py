"""Projektblad (inputs.json) med stoppkontroll enligt docs/projektstandard.md.

    from femtools import projektblad
    pb = projektblad.load("inputs.json")   # stoppar om SF, lastfall eller material saknas

Laster, säkerhetsfaktorer, material och upplag får aldrig antas. Fält som saknas
eller är null ger ProjektbladError med en lista över vad agenten ska fråga om.
"""
import json

KRAV = {
    "case": "namn på beräkningsfallet (YYMMDD_namn)",
    "geometri.step": "STEP-fil eller geometribeskrivning",
    "material": "materialnamn enligt femtools/materials.json",
    "sakerhetsfaktorer.SF_y": "säkerhetsfaktor mot sträckgräns",
    "lastfall": "minst ett lastfall med storlek, riktning, angreppsyta och källa",
    "randvillkor": "minst ett upplag med låsta frihetsgrader och motivering",
    "utvardering.storhet": "vilken storhet som utvärderas",
    "utvardering.plats": "var storheten utvärderas",
}


class ProjektbladError(ValueError):
    pass


def _get(d, dotted):
    for k in dotted.split("."):
        if not isinstance(d, dict) or d.get(k) in (None, "", [], {}):
            return None
        d = d[k]
    return d


KRAV_5 = {   # från standard 5.0
    "betraktarvy.fram": "produktens framsida: riktning från produkten mot betraktaren i monterat läge, t.ex. [0, 1, 0]",
    "betraktarvy.upp": "uppåt i monterat läge, t.ex. [0, 0, 1]",
}


def _major(pb):
    try:
        return int(str(pb.get("standardversion") or "0").split(".")[0])
    except ValueError:
        return 0


def check(pb):
    """Lista över saknade obligatoriska fält (tom lista = OK). Betraktarvy krävs från standardversion 5.0."""
    krav = dict(KRAV, **(KRAV_5 if _major(pb) >= 5 else {}))
    saknas = [f"{k}: {v}" for k, v in krav.items() if _get(pb, k) is None]
    for i, lf in enumerate(pb.get("lastfall") or []):
        for k in ("namn", "storlek", "riktning", "yta", "kalla"):
            if lf.get(k) in (None, ""):
                saknas.append(f"lastfall[{i}].{k}")
    for i, rv in enumerate(pb.get("randvillkor") or []):
        for k in ("rand", "typ", "lasta_dof", "motivering"):
            if rv.get(k) in (None, "", []):
                saknas.append(f"randvillkor[{i}].{k}")
    return saknas


def load(path):
    with open(path, encoding="utf-8") as f:
        pb = json.load(f)
    saknas = check(pb)
    if saknas:
        raise ProjektbladError("Projektbladet saknar värden, fråga Marcus innan modellering:\n  - "
                               + "\n  - ".join(saknas))
    return pb
