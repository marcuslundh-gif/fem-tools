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
apt-get install -y -qq openfoam                                                              # turbulent CFD (v1912), bara vid behov
pip install -q -e .            # i repot; annars: pip install git+https://github.com/marcuslundh-gif/fem-tools
python -m femtools.verify      # måste sluta med "femtools verifierad OK"
```

## Beräkningsfall
- Fallmapp utanför repot: `~/fem-cases/<YYMMDD_case>/` med `inputs.json` (kopia av `mallar/projektblad.json`), `src/` (skript och `<case>_handcalc.py`), `results/`, `fig/`, `report.md`.
- Leverans: PDF-rapport (stående A4, `YYMMDD <namn>.pdf`) och resultat till projektets mapp på Drive, `FEM/<YYMMDD_case>/`. Kunddata bara i projektets `Kunddata/`.
- Figurer med `femtools.report` (L2-projektion, snitt, enhetlig stil). Fig. 1–7 i fem-rapport är minimum: lägg till de figurer som behövs för resultatet eller beslutet och motivera varje extra figur i bildtexten.
- Interaktiv 3D-visare är standardleverans bredvid PDF: `femtools.webreport` med nät, alla lastfall, termik och CFD (`Report.flow` + `openfoam.read_case`), samma tabeller. Publiceras som artifact; tabellerna kopieras tabbseparerade till PPT/Excel. Spara nodala fält till fil i varje körning.
- Betraktarvy: alla vyer (visare, lastskiss, resultatbilder, CAD-vyer, skisser) framifrån som en människa ser produkten i monterat läge. Framsidan står i projektbladet (`betraktarvy.fram`, `betraktarvy.upp`); fråga om den inte är självklar.
- Rapportens början anger standardversion, `femtools.__version__`, commit-hash och projektbladet.
- Avsluta varje fall med förslag till revision av standarden. Ändringen görs som pull request i detta repo; Marcus godkänner genom att merga.

## Prestanda
- Molnmiljön: 4 kärnor, 15 GB (uppmätt 2026-11-26). Cirka 25 GB per miljon DOF, tak cirka 400 000 DOF för slanka detaljer, lägre för kompakta. Mät alltid minnet för de två grövsta näten och extrapolera, stoppa och fråga över 70 % av RAM.
- Assembleringen tar ungefär dubbla lösningstiden: kör flera lastfall på samma matris.
- Körningar över cirka 10 minuter i bakgrunden. Spara resultat till fil efter varje nät.

## CFD
- Laminärt (Re under omslaget, rör cirka 2300): `femtools.flow` (scikit-fem). Turbulent: `femtools.openfoam` (OpenFOAM, SI-enheter).
- Tunga turbulenta 3D-fall körs på smartfish (WSL2, OpenFOAM v2406, 16 trådar, 16 GB): fallet i `~/fem-cases/` i WSL, aldrig i Drive-synkad mapp eller på /mnt/c.
- RANS-modellens fel (k-ω SST) är cirka ±3–5 % i friktion. Ange det i rapporten och använd säkerhetsmarginal därefter.
- CFD-nätstudie: minst två nät (kvot ≥ 1,5) för varje storhet i ett beslut; spara och återanvänd nätet (HXT ger olika nät, Δp ±7 % vid samma h).
- STEP → OpenFOAM-nät: gmsh-tetraeder av vätskekroppen med fysiska grupper (inlet, outlet, walls, fluid), `Mesh.MshFileVersion 2.2`, `gmshToFoam`, `transformPoints -scale '(0.001 0.001 0.001)'`, sätt sedan `walls` till `type wall` i `constant/polyMesh/boundary`. Inloppshastighet = Q / nätets inloppsarea.
- OpenFOAM v1912 (apt): inga function objects (surfaceFieldValue kraschar), efterbehandla i Python ur polyMesh och tidskataloger (`writeFormat ascii`). Parallellt: `decomposePar` med `method hierarchical` (scotch saknas), `mpirun --allow-run-as-root -np 4 simpleFoam -parallel`, `reconstructPar`.
- Tetraedernät med SIMPLE stannar ofta vid residualer 1e-3–1e-2: redovisa Δp som medel och min–max över andra halvan av iterationerna.

## Kontext och tokens
- Skriv aldrig ut arrayer eller hela loggar: max 20 rader sammanfattning, resten till fil.
- Granska max 2 figurer per steg, kontrollera resten via filstorlek och nyckelvärden.
- Läs aldrig STEP som text. Använd `help(femtools.<funktion>)` i stället för att läsa källkoden.
