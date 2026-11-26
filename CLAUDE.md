# fem-tools

Verifierade FEM-verktyg (gmsh + scikit-fem P2) och regler för alla FEM-beräkningar. Repot är master för reglerna och är publikt: lägg aldrig projektblad, laster, STEP, resultat eller rapporter från kund- eller jobbprojekt här.

## Regler
- Standard: `docs/projektstandard.md` (läs bara de avsnitt som behövs). Register över verifierade verktyg: `docs/verifiering.md`.
- Arbetsgång: skillen fem-berakning. Rapport: skillen fem-rapport. Handräkning med enheter: skillen enhetskontroll.
- Gissa aldrig last, material, upplag eller säkerhetsfaktor: fråga. Projektbladet fylls i före modellering och läses med `femtools.projektblad.load()`, som stoppar om något saknas.
- Importera femtools, skriv inte om verifierad kod. Nytt verktyg används först när det har ett test i `femtools/verify.py` och en rad i `docs/verifiering.md`.
- Enheter mm, N, MPa, t, s. Materialvärden bara via `femtools.material()`.

## Uppstart (varje ny container)
```bash
apt-get update -qq && apt-get install -y -qq libglu1-mesa libxcursor1 libxinerama1 libxft2   # gmsh kräver libGLU
pip install -q -e .            # i repot; annars: pip install git+https://github.com/marcuslundh-gif/fem-tools
python -m femtools.verify      # måste sluta med "femtools verifierad OK"
```

## Beräkningsfall
- Fallmapp utanför repot: `~/fem-cases/<YYMMDD_case>/` med `inputs.json` (kopia av `mallar/projektblad.json`), `src/` (skript och `<case>_handcalc.py`), `results/`, `fig/`, `report.md`.
- Leverans: PDF-rapport (stående A4, `YYMMDD <namn>.pdf`) och resultat till projektets mapp på Drive, `FEM/<YYMMDD_case>/`. Kunddata bara i projektets `Kunddata/`.
- Rapportens början anger standardversion, `femtools.__version__`, commit-hash och projektbladet.
- Avsluta varje fall med förslag till revision av standarden. Ändringen görs som pull request i detta repo; Marcus godkänner genom att merga.

## Prestanda
- Molnmiljön: 4 kärnor, 15 GB (uppmätt 2026-11-26). DOF-taket är inte uppmätt här: mät minnet för de två grövsta näten och extrapolera, stoppa och fråga över 70 % av RAM.
- Körningar över cirka 10 minuter i bakgrunden. Spara resultat till fil efter varje nät.

## Kontext och tokens
- Skriv aldrig ut arrayer eller hela loggar: max 20 rader sammanfattning, resten till fil.
- Granska max 2 figurer per steg, kontrollera resten via filstorlek och nyckelvärden.
- Läs aldrig STEP som text. Använd `help(femtools.<funktion>)` i stället för att läsa källkoden.
