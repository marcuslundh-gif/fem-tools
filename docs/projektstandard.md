# Projektstandard för FEM-beräkningar

Version 5.1, 2026-11-27. Ägare: Marcus Lundh.

Master: denna fil i GitHub-repot `fem-tools` (`docs/projektstandard.md`). Kopia i Google Drive-mappen "Projektstandard för FEM-beräkningar" synkas efter varje godkänd ändring.

## Revisionslogg

| Version | Datum | Ändring |
| :--- | :--- | :--- |
| 1.0 | 2026-11-12 | Första version: enheter, två material, randvillkor, acceptanskriterier. |
| 2.0 | 2026-11-17 | Säkerhetsfaktor flyttad till projektblad per projekt. Materialdata tjockleksberoende och med HAZ. Nya avsnitt: meshning, lösare och prestanda, verifiering, utmattning/egenfrekvens/knäckning, verktygsregister, erfarenhetsbank och revisionsprocess. RBE2/RBE3 översatt till scikit-fem. |
| 3.0 | 2026-11-19 | Rapportkrav skärpta: lastskiss med randvillkor först i rapporten, obligatoriska resultatbilder (deformation, spänning, snitt, konvergens) med förklarande bildtext. Detaljerade krav flyttade till skillen fem-rapport, standarden har kärnreglerna. |
| 3.1 | 2026-11-19 | femtools.report 1.1.0 verifierad och införd i verktygsregistret (öppen punkt stängd). Ny erfarenhetspost om apt-get update. |
| 3.2 | 2026-11-19 | Värmeutvidgning α införd för alla material. Nya material med typvärden (ej minimivärden): PA12 SLS (EOS PA 2200), PA12 GB SLS (EOS PA 3200 GF), PA12 MJF (HP 3D HR PA 12), Alumide. Värden kontrollerade mot EOS, HP, Materialise och 3D Systems. |
| 4.0 | 2026-11-26 | Master flyttad till GitHub-repot fem-tools så att alla agenter (Claude Code i molnet, lokalt och på jobbet) kör samma regler. claude.ai-specifika delar borttagna (2 kärnor/7 GB, kör klart inom svaret, synk till claude.ai-projektet). Verktygsregistret flyttat till `docs/verifiering.md` och rättat mot koden. Materialdata maskinläsbar i `femtools/materials.json`. Projektblad med stoppkontroll (`femtools.projektblad`). Regel för ANTAGET förtydligad, resonansmarginal 20 % även i skillen, sekretessregel för publikt repo. Revision via pull request. |
| 4.1 | 2026-11-26 | Nya material: 1.4301 (rör, plattstång kapad och varmvalsad) och 42CrMo4+QT. Nya erfarenhetsposter: minnestak för kompakta solidmodeller, randfrihetsgrader vid summering av laster, solida rör i CAD, singulär topp i svetsrot och linjäriserad hot-spot. Regel: lastskiss och projektblad före beräkning gäller även handräkning och tippning. |
| 4.2 | 2026-11-26 | femtools 1.3.0: `static_ext` (tryck, egenvikt, rotation, termisk last, randvillkor per komponent) `heat`/`heat_transient` (värmeledning med fast temperatur, flöde, konvektion, värmekälla) `buckling` (linjär knäckning, geometrisk styvhet från förspänningen), `flow` (laminär Stokes/Navier-Stokes, Taylor-Hood) `openfoam` (turbulent CFD med OpenFOAM, k-ω SST) och `webreport` (interaktiv HTML-rapport som komplement till PDF). Verifierade mot Lamé, roterande rör, hängande stång, rör med konvektion, Timoshenkos termiska rörspänning, inspänd stång, transient platta, Euler-pelare fall 1 och 3, Kovasznay, Hagen-Poiseuille och turbulent rör mot Colebrook; se `docs/verifiering.md`. Bygger på 4.1 (PR #2). |
| 4.3 | 2026-11-26 | Erfarenheter från ett fall med STEP från Fusion 360, termik, superponerade lastfall och CFD i OpenFOAM: STEP-enhet, gmsh-fel på periodiska ytor, STEP → OpenFOAM-nät, efterbehandling utan surfaceFieldValue, minnestak och konvergens med liten nätkvot, AMG-prestanda. 1.4301 rör med temperaturberoende Rp0,2, λ och c_p i `materials.json`. Bygger på 4.2 (PR #3). |
| 5.0 | 2026-11-27 | Nya regler: betraktarvy framifrån för alla vyer (3D-visare, lastskiss, resultatbilder, CAD-vyer och skisser) med framsidan i projektbladet (`betraktarvy`, krävs av `projektblad.check` från 5.0); interaktiv 3D-visare (`femtools.webreport`) är standardleverans i varje beräkning, med CFD-fält när CFD ingår; CFD kräver nätstudie med minst två nät och sparat nät. femtools 1.4.0: `webreport` med betraktarvy, snittplan och CFD-lager (`flow`, `slice_tets`, `streamlines`), `openfoam.read_case`. |
| 5.1 | 2026-11-27 | femtools 1.5.0: `shell` (tunnväggiga rotationsskal i CalculiX, S8R, icke-linjär kollaps med NLGEOM och plasticitet, Koiter-imperfektion), `webreport.surface` (skalytor), typvärden för metaller i `material()`. Nya material: AA3104-H19 och AA5182-H48 (burkplåt). CalculiX med MKL Pardiso (`tools/build_ccx_pardiso.sh`). Erfarenheter från knäckning av tunt skal: `*BUCKLE` med S8R ger falska moder, falsens hårnål kräver grovt nät, SPOOLES långsam. Tillägg 2026-11-28: strukturerat nät med S8R och krav på fallande gren för skalkollaps, `orient()` för normalernas riktning, reaktion i alla upplag för olinjära körningar, ingen Richardson när kollapsformen byter plats. |

## Användning

Standarden gäller alla FEM- och hållfasthetsberäkningar som görs med femtools, oavsett miljö: Claude Code i molnet, lokalt och agenten på jobbet. Arbetsgången står i skillen fem-berakning, rapportkraven i skillen fem-rapport och den praktiska uppstarten i `CLAUDE.md`. Projektspecifika värden (säkerhetsfaktorer, lastfall, krav) står i projektbladet för varje beräkningsfall, inte här. Om projektbladet och standarden säger olika gäller projektbladet för just det projektet. Standarden revideras enligt avsnittet Revisionsprocess.

Sekretess: repot är publikt. Här ligger bara regler, verktyg, generiska verifieringsfall och generiska lärdomar. Projektblad, laster, STEP-filer, resultat och rapporter från kund- eller jobbprojekt läggs aldrig i repot, utan i projektets mapp på Drive (Kunddata/ när det är kunddata) eller i Claude.

## Enhetssystem

All beräkning och redovisning sker i mm, N, MPa, ton, s. Densitet anges i ton/mm³, vilket ger egenfrekvenser direkt i Hz. Kontrollera STEP-filens enhet mot bounding box innan meshning.

## Projektblad (per projekt)

Varje beräkningsfall har ett projektblad, `inputs.json` i fallmappen (utanför repot, se Rapport och spårbarhet), som fylls i innan modellering. Mall: `mallar/projektblad.json`. Skripten läser projektbladet med `femtools.projektblad.load()`, som stoppar om säkerhetsfaktorer, lastfall eller material saknas. Saknas ett värde ska agenten fråga.

Projektbladet gäller hela beräkningsfallet, även handräkningar, tippning och egenfrekvenser utanför FEM. Laster som bara används i en handräkning (till exempel personlast mot en bänk) får inte antas utan att stå i projektbladet.

Vad som får antas:

- Aldrig: laster, lastfall, säkerhetsfaktorer, material, upplag. Fråga.
- Får antas konservativt och märkas [ANTAGET]: geometridetaljer utanför utvärderingsområdet (radier, fasningar), friktionstal, imperfektioner, ytfinhet för utmattning. Antagandet listas i rapporten med påverkan på resultatet.

| Parameter | Exempel | Kommentar |
| :--- | :--- | :--- |
| SF mot sträckgräns, SF_y | 1,5 | Bestäms per projekt. |
| SF mot brott, SF_u | 2,0 | Bestäms per projekt, används om relevant. |
| SF mot utmattning, SF_f | 1,35 | Eller partialkoefficient enligt tillämplig standard. |
| SF mot knäckning, SF_b | 3,0 | Gäller linjär buckling-faktor. |
| Deformationskrav | L/300 | Anges med referenslängd L och vilken punkt som utvärderas. |
| Lastfall och lastfaktorer | drift, prov, olycka | Varje lastfall med storlek, riktning, angreppsyta och källa. |
| Livslängd och lastcykler | 10⁶ cykler | Krävs för utmattningsbedömning. |
| Drivfrekvenser | 10 rpm, 50 Hz | Krävs för egenfrekvensbedömning. |
| Temperatur och miljö | 20 °C, inomhus | Påverkar materialdata och korrosion. |
| Tillämplig standard | EN 1993, EN 1999, IIW, kundkrav | Styr partialkoefficienter och utmattningsklasser. |

## Materialdata

### Minimivärden (metaller)

Maskinläsbar master för alla materialvärden: `femtools/materials.json` (läses med `femtools.material(namn, t)`). Tabellerna nedan är en läsbar sammanställning. Ändras ett värde ändras båda i samma pull request; `python -m femtools.verify` kontrollerar att JSON-filen är komplett.

Värden är minimivärden enligt angiven standard. Använd materialcertifikat när det finns. Nytt material läggs till med källa, kontrollerat mot tillverkarens datablad och minst en annan tillverkare.

| Material | Standard och produkt | Tjocklek t [mm] | E [MPa] | ν [-] | Rp0,2 / ReH [MPa] | Rm [MPa] | ρ [ton/mm³] | α [1/K] |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| Stål S355 | EN 10025-2 | t ≤ 16 | 210 000 | 0,30 | 355 | 470–630 | 7,85e-9 | 1,2e-5 |
| Stål S355 | EN 10025-2 | 16 < t ≤ 40 | 210 000 | 0,30 | 345 | 470–630 | 7,85e-9 | 1,2e-5 |
| Stål S355 | EN 10025-2 | 40 < t ≤ 63 | 210 000 | 0,30 | 335 | 470–630 | 7,85e-9 | 1,2e-5 |
| Alu 6082-T6 | EN 1999-1-1 tab 3.2b, strängpressat | t ≤ 5 | 70 000 | 0,33 | 250 | 290 | 2,70e-9 | 2,3e-5 |
| Alu 6082-T6 | EN 1999-1-1 tab 3.2b, strängpressat | 5 < t ≤ 15 | 70 000 | 0,33 | 260 | 310 | 2,70e-9 | 2,3e-5 |
| Alu 6082-T6, svetspåverkad zon (HAZ) | EN 1999-1-1 tab 3.2b, MIG, t ≤ 15 | t ≤ 15 | 70 000 | 0,33 | 125 | 185 | 2,70e-9 | 2,3e-5 |
| Rostfritt 1.4301, rör | Tibnor Rostfria konstruktionsrör | t ≤ 8 | 200 000\* | 0,30\* | 210 | 520–720 | 7,9e-9\* | 1,6e-5\* |
| Rostfritt 1.4301, plattstång kapad från plåt | Tibnor 1.4301/1.4307 | 3 ≤ t ≤ 12 | 200 000\* | 0,30\* | 210 | 520–720 | 7,9e-9\* | 1,6e-5\* |
| Rostfritt 1.4301, plattstång varmvalsad | Tibnor 1.4301/1.4307 | 8 ≤ t ≤ 25 | 200 000\* | 0,30\* | 190 | 500–700 | 7,9e-9\* | 1,6e-5\* |
| 42CrMo4+QT, rund (t = diameter) | Swiss Steel, EN ISO 683-2 | 16 < d ≤ 40 | 210 000\* | 0,30\* | 750 | 1 000–1 200 | 7,85e-9\* | 1,2e-5\* |
| 42CrMo4+QT, rund (t = diameter) | Swiss Steel, EN ISO 683-2 | 40 < d ≤ 80 | 210 000\* | 0,30\* | 650 | 900–1 100 | 7,85e-9\* | 1,2e-5\* |

\* ANTAGET: typvärden enligt EN 10088-1 respektive allmän stållitteratur, inte kontrollerade mot källa. Kontrollera innan de används i en termisk analys eller egenfrekvensberäkning med liten marginal. Härdat och anlöpt stål: dimensionerande diameter är ämnets, inte den färdiga detaljens (en axel med fläns Ø140 svarvas ur ett grövre ämne).

α för stål enligt EN 1993-1-1 avsnitt 3.2.6, för aluminium enligt EN 1999-1-1 avsnitt 3.2.5.

Version 1.0 angav 6082-T6 med 260/310 MPa utan tjockleksgräns. För väggar under 5 mm gäller 250/290 MPa. Svetsade aluminiumdetaljer utvärderas alltid med HAZ-värden inom HAZ-zonen.

### Typvärden (3D-printade polymerer)

TYPVÄRDE: leverantörens medelvärden, inte minimivärden. Rapporten ska ange att typvärden använts. Ingen sträckgräns anges av tillverkarna: utvärdera mot Rm med säkerhetsfaktor enligt projektbladet. Använd Z-värdet när lastriktningen mot byggriktningen är okänd. Värden markerade \* är ANTAGET (saknas i tillverkarens datablad). α gäller rumstemperatur och ökar över cirka 50 °C.

| Material | Källa | E XY / Z [MPa] | ν [-] | Rm XY / Z [MPa] | ρ [ton/mm³] | α [1/K] |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| PA12 SLS (EOS PA 2200) | EOS Material Data Sheet PA 2200 Balance; α: EOS Product information PA 2200 (ISO 11359) | 1 650 / 1 650 | 0,35\* | 48 / 42 | 9,30e-10 | 1,09e-4 |
| PA12 GB SLS (EOS PA 3200 GF, 40 vikt-% glaskulor) | EOS Material Data Sheet PA 3200 GF | 3 200 / 2 500 | 0,32\* | 51 / 47 | 1,22e-9 | 9,0e-5\* |
| PA12 MJF (HP 3D HR PA 12) | HP datablad HP 3D HR PA 12 | 1 700 / 1 800 | 0,35\* | 48 / 48 | 1,01e-9 | 1,1e-4\* |
| Alumide (EOS PA12-Al, 50 vikt-% Al) | EOS Alumide; Materialise PA-AF | 3 800 / 3 200 | 0,34\* | 48 / 40 | 1,36e-9 | 8,5e-5\* |

Kontroll mot andra tillverkare (2026-11-19): Materialise PA12 SLS 1 650 MPa / 48 MPa, 3D Systems DuraForm ProX PA 1 770 / 47 MPa, Materialise PA-GF 3 200 / 51 MPa, Materialise PA-AF 3 800 / 48 MPa. α för Alumide uppskattad med blandningsregel (cirka 27 vol-% Al) och måste ligga mellan aluminium (2,3e-5) och PA12 (1,1e-4); 3D Systems DuraForm ProX AF+ anger 1,45e-4 för 0–145 °C. Termiska analyser med Alumide kräver mätt värde eller känslighetsstudie.


### Typvärden (metaller)

TYPVÄRDE även för metall när bara uppmätta värden finns (t.ex. tunn burkplåt utan EN-minimivärde i tjockleken). `material()` returnerar `typ` = TYPVARDE; rapporten anger det och kollapslasten redovisas med spridning från känslighetsstudie.

| Material | Källa | t [mm] | E [MPa] | ν [-] | Rp0,2 [MPa] | Rm [MPa] | ρ [ton/mm³] |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| AA3104-H19 (burkkropp) | Wędrychowicz m.fl. 2023 (dragprov, typvärde); Steelforce/Signi | 0–0,5 | 69 000 | 0,33 | 293 | 328 | 2,72e-9 |
| AA5182-H48 (burklock) | Steelforce (minimivärden), E/ν/ρ ANTAGET | 0–0,5 | 70 000 | 0,33 | 329 | 370 | 2,65e-9 |

## Randvillkor

- Undvik helt stumma inspänningar (alla 6 DOF låsta) om det inte motsvarar en faktisk fastsvetsning i mycket styv struktur.
- Lås exakt de stelkroppsrörelser som krävs (6 i 3D) och kontrollera via reaktionskrafterna.
- Skruvhål: fördela lasten realistiskt. I kommersiell FEM används RBE2 eller RBE3. I scikit-fem gäller motsvarande: lasten läggs som tryck eller skjuvning på bricka- eller anliggningsytan (RBE3-motsvarighet), och upplag modelleras som låsning av ringytan under brickan eller skruvskallen, inte hela hålets mantelyta (RBE2-motsvarighet).
- Förspända skruvförband utvärderas enligt VDI 2230. FEM används för klämlängdens styvhet och lastfördelning.
- Osäker inspänning räknas i båda gränsfallen (ledad och fast). Rapportera vilket som är konservativt för den utvärderade storheten.
- Utnyttja symmetriplan när geometri, last och randvillkor är symmetriska. Symmetri får inte användas för egenfrekvenser utan att antisymmetriska moder hanteras.
- Kontakt modelleras bara när den avgör resultatet. Annars Hertz analytiskt eller förenklad lastöverföring.
- Lastskissen (se Rapport och spårbarhet) ritas och granskas mot projektbladet innan lösning. Stämmer den inte löses inget.

## Meshning

Standardinställningar i gmsh (finns i `femtools/mesh.py`, funktionerna `setup` och `refine_and_mesh`):

| Inställning | Värde | Syfte |
| :--- | :--- | :--- |
| 3D-algoritm | HXT (Mesh.Algorithm3D = 10) | Parallell, snabb för stora nät. |
| 2D-algoritm | Frontal-Delaunay (Mesh.Algorithm = 6) | Jämn ytkvalitet. |
| Trådar | General.NumThreads = 0 | Alla kärnor. |
| Storlek från krökning | 24 element per 2π | Hål och radier löses upp automatiskt. |
| Lokal förfining | Distance + Threshold på kritiska ytor | Fint där spänningen utvärderas, grovt annars. |
| Geometrireparation | OCCFixDegenerated, OCCFixSmallEdges, OCCFixSmallFaces | Färre meshfel från STEP. |
| Elementtyp | P1-nät från gmsh, P2-basis i scikit-fem | Kvadratiskt fält. Raka kanter, därför lokal förfining vid krökta kritiska ytor. |

Riktvärden för elementstorlek:

- Minst 8 kvadratiska element per 90° båge i en radie där spänningen utvärderas.
- Minst 2 kvadratiska element genom tjockleken i böjbelastade tunnväggiga delar, helst 3. Är tjockleken under ungefär 1/20 av spännvidden: överväg skalmodell eller balkmodell.
- Lokal storlek h_loc vid kritisk yta, övergång till global storlek h_glob mellan 2·h_loc och 10·h_loc från ytan.
- Konvergensstudien skalar h_glob och h_loc med samma faktor, minst 1,5 mellan nät.
- Begränsar minnet finaste nätet: behåll h_glob och skala h_loc med 1,5 nedåt från det finaste nät som ryms. Ett nät med kvot under 1,2 ger brusig Richardson-extrapolering. Redovisa alla nivåer i tabell och ange avvikelsen från kravet.

Kvalitetskrav innan lösning:

- Elementkvalitet gamma (gmsh) minst 0,1 i hela nätet och minst 0,3 i utvärderingsområdet. Sämre element redovisas med position.
- Inga lösa noder eller dubbla volymer. Kontrollera att nätets volym stämmer med CAD-volymen inom 0,5 %.

STEP från CAD:

- Kontrollera enheten med nätets bounding box mot ritningen innan nätning. Fusion 360 skriver alltid STEP i cm, även när dokumentet är i mm. Konvertera med OCCT (`STEPCAFControl_Writer`, `write.step.unit MM`) eller sätt `Geometry.OCCTargetUnit = "MM"` vid import.
- gmsh läser inte kroppsnamnen i STEP. Identifiera solider via ordning och volym och redovisa volymen mot CAD.
- Ger gmsh "Impossible to mesh periodic surface" (sömmar på cylindrar som korsas av radier): dela soliden i CAD i kvartar eller halvor genom symmetriplanen, exportera och slå ihop med `occ.fragment`, som ger ett konformt nät över snittytorna.

Förenkling av geometri (defeaturing):

- Ta bort små radier, fasningar, hål, gängor och text som ligger utanför lastvägen och utvärderingsområdet. Dokumentera varje borttagning i rapporten.
- Behåll alla detaljer i lastvägen och i utvärderingsområdet. Inre skarpa hörn i utvärderingsområdet modelleras med verklig radie från ritningen.

Effektivitet:

- Använd symmetri när det är tillåtet.
- Spara nät per nivå (`results/<case>/mesh_<nivå>.msh`) så att de kan återanvändas mellan lastfall.
- Kör flera lastfall på samma nät med samma faktorisering av styvhetsmatrisen.

## Lösare och prestanda

- Direktlösare: pypardiso (Intel MKL, flertrådad). scipy spsolve (SuperLU, en tråd) endast under cirka 100 000 DOF.
- Mycket stora modeller: konjugerade gradienter med algebraisk multigrid (pyamg) som förkonditionerare.
- Logga DOF, meshtid, lösningstid, minnestopp och antal trådar för varje körning.
- Mät minnet för de två grövsta näten och extrapolera innan det finaste körs. Stoppa och fråga om uppskattningen överstiger 70 % av tillgängligt RAM.
- Långa körningar (över cirka 10 minuter) körs i bakgrunden med resultat sparat till fil efter varje nät, så att en avbruten studie kan återupptas.
- Styvhetsmatrisen integreras med intorder 2 (exakt för P2 på raka tetraedrar). Spänningar utvärderas med intorder 4. Egenfrekvenser löses med shift-invert där pypardiso faktoriserar.
- Uppmätt 2026-11-26 i Claude Code-molnet (4 kärnor, 15 GB), P2-statik på slank balkgeometri:

| DOF | Assemblering [s] | Lösning [s] | Minnestopp [GB] |
| :--- | :--- | :--- | :--- |
| 49 000 | 3,0 | 0,9 | 0,85 |
| 112 000 | 7,2 | 3,0 | 2,2 |
| 213 000 | 19,7 | 8,3 | 4,6 |
| 362 000 | 38,5 | 18,2 | 8,9 |

- Minnet växer ungefär 25 GB per miljon DOF. Med 70 %-regeln ger det ett tak kring 400 000 DOF i molnet. Kompakta detaljer ger mer fill-in än en slank balk och därmed lägre tak: mät alltid de två grövsta näten. Större modeller körs lokalt.
- Assembleringen tar mer tid än lösningen. Kör därför flera lastfall på samma assemblerade matris i stället för att anropa `static` en gång per lastfall.
- Den tidigare claude.ai-containern (2 kärnor, 7 GB) klarade cirka 300 000 DOF.
- CalculiX (skal, icke-linjärt): apt-versionen använder SPOOLES. Bygg `ccx_pardiso` med `tools/build_ccx_pardiso.sh` (MKL Pardiso, 2,3 gånger snabbare på samma modell). `femtools.shell.ccx_exe()` väljer den när den finns.
- Kompakt svetsad detalj (rör, flänsar, konsoler, 5 superponerade lastfall med spänningar sparade): 206 000 DOF 7,9 GB, 265 000 DOF 8,0 GB, 382 000 DOF 11,4 GB, 90–210 s per nät.
- pyamg (smoothed aggregation med stelkroppsmoder, CG, tol 1e-10) löste 265 000 DOF på 374 s per lastfall mot 20 s med pypardiso: 15–30 gånger långsammare. Bara för ett fåtal lastfall när pypardiso inte ryms i minnet; verifiera mot pypardiso på ett grövre nät.

## Verifiering

Obligatoriskt för varje beräkningsfall:

- Jämvikt: summan av reaktionskrafter plus pålagda laster ska vara under 0,1 % av pålagd last, per komponent.
- Förskjutningsstyrda och olinjära körningar: skriv ut reaktionen i alla upplag i varje körning (CalculiX `*NODE PRINT,NSET=...,TOTALS=ONLY` per upplag), så att jämvikten kan kontrolleras utan omkörning.
- Nätkonvergens: minst 3 nät. Godkänt när utvärderingsstorheten ändras mindre än 3 % mellan de två finaste. Redovisa tabell och Richardson-extrapolerat värde.
- Handräkning eller analytisk lösning för storleksordning, i ett eget skript. FEM och handräkning ska stämma inom cirka 10 % i en jämförbar punkt, annars förklaras skillnaden.
- Rimlighet: deformationsform (skalad plott), tecken och var maxvärden hamnar.
- Singularitet: en topp som växer för varje förfining är en singularitet och rapporteras inte som spänning. Utvärdera med verklig radie, linjariserad spänning eller hot-spot-metod och ange vilken.

## CFD

- Laminärt med `femtools.flow`, turbulent med `femtools.openfoam` (se `CLAUDE.md` för nätning från STEP och efterbehandling).
- Nätstudie krävs för varje storhet som används i ett beslut (Δp, flödesfördelning, värmeövergång): minst två nät med cellstorlek i kvot minst 1,5, helst tre. Redovisa storheten per nät och ange skillnaden mellan de två finaste som osäkerhet. Utan nätstudie redovisas resultatet som osäkert med en uppskattad marginal och motivering.
- Spara CFD-nätet (`fluid.msh` och `constant/polyMesh`) per nivå och återanvänd det. gmsh HXT ger inte samma nät mellan körningar: samma h gav Δp +6–7 % i en upprepad körning.
- Med SIMPLE på tetraedrar som stannar vid residualer 1e-3–1e-2 redovisas storheten som medel och min–max över andra halvan av iterationerna. RANS-modellens fel (k-ω SST cirka ±3–5 % i friktion) läggs till nätosäkerheten.

## Visning och betraktarvy

- Alla vyer visar produkten som en människa ser den när hon står framför produkten i monterat läge: framifrån med uppåt uppåt. En väggmonterad produkt ses från rummet, inte från väggen. En bänk- eller golvstående produkt ses från användarens sida.
- Framsidan anges i projektbladet som `betraktarvy.fram` (riktning från produkten mot betraktaren) och `betraktarvy.upp`. Är den inte självklar ur ritning eller montage frågar agenten. `projektblad.check` stoppar om den saknas (standardversion 5.0 och senare).
- Regeln gäller 3D-visaren (startvy "Fram iso"), lastskissen, resultatbilder och snitt i PDF-rapporten, CAD-vyer (skärmbilder och renderingar från Fusion 360 eller SolidWorks) och skisser. Snitt och detaljvyer från andra håll får läggas till och märks med betraktningsriktning ("sedd från väggen", "underifrån").
- CAD-ritning, även utan beräkning (Fusion 360 via API, SolidWorks, skisser): bestäm framsidan före modellering och fråga om den inte är självklar. Sätt hemvyn (Home/ViewCube i Fusion, vyn Front i SolidWorks) så att produktens framsida är Front och hemvyn är "Fram iso", och ta skärmbilder därifrån. Arbetsgången står i skillen cad-ritning.
- Interaktiv 3D-visare är standardleverans i varje beräkning, bredvid PDF-rapporten: `femtools.webreport` med nät, alla lastfall (von Mises, huvudspänning, förskjutning med skalfaktor), termik när den ingår och CFD-fält (`Report.flow` med `openfoam.read_case`: |U| i snitt, väggtryck, strömlinjer) när CFD ingår. Samma tabeller som rapporten. Publiceras som artifact, utan datornamn eller Drive-sökvägar.
- Spara nodala fält till fil i varje körning (t.ex. pickle med p, t, nodvärden per lastfall och temperatur) så att visaren kan byggas utan omkörning.

## Acceptanskriterier

- Säkerhetsfaktorer enligt projektbladet. Saknas de beräknas inget slutresultat, agenten frågar.
- Max deformation: L/300 om projektbladet inte anger annat.
- Spänningsutvärdering: von Mises för segt stål och aluminium. Lokala singulariteter vid vassa inre hörn ignoreras om de inte är utmattningskritiska. Är de utmattningskritiska utvärderas de med hot-spot- eller kerbspänningsmetod.
- Flagga alltid: von Mises över sträckgränsen (linjär analys ogiltig lokalt), knäckning som kan vara dimensionerande, stora deformationer (över cirka 1/10 av tjockleken).

## Utmattning, egenfrekvens och knäckning

- Utmattning, grundmaterial: Kf = 1 + q(Kt − 1). Goodman eller Gerber med korrigerad utmattningsgräns (yta, storlek, tillförlitlighet). Säkerhetsfaktor mot utmattning redovisas separat från statisk.
- Utmattning, svetsar: nominell spänning eller hot-spot-spänningsvidd mot FAT-klass enligt EN 1993-1-9, EN 1999-1-3 eller IIW. Variabel amplitud med Palmgren-Miner. Ange FAT-klass och källa.
- Egenfrekvens: modalanalys med verkliga randvillkor och tillsatta massor. Marginal mot drivfrekvenser och deras övertoner anges i projektbladet. Riktvärde om inget anges: minst 20 % (gäller även skillen fem-berakning).
- Knäckning: linjär buckling-analys ger lastfaktor. Tunnväggiga eller imperfektionskänsliga konstruktioner kräver icke-linjär analys eller standardens reduktionskurvor.
- Tunna skal (R/t över cirka 100): icke-linjär kollaps med `femtools.shell` (NLGEOM, plasticitet, förskjutningsstyrd last), kollapslast = högsta kraft före första fallande gren. Imperfektion som Koiter-form med halvvåglängden l = π√(Rt)/(12(1−ν²))^¼ och amplitud ξ·t i ett svep (minst ξ 0; 0,1; 0,25; 0,5), med en svag omkretsstörning så att symmetrin bryts. Redovisa jämförelse mot NASA SP-8007 (knockdown γ). Linjär `*BUCKLE` i CalculiX används inte för S8R (falska moder). Nät: strukturerat (meridian × omkrets, `revolve_shell(structured=True)`) med S8R och nätstudie på hela modellen; inte ostrukturerat nät med lokala förfiningszoner. Kollapslasten godtas bara när kraftkurvan har en fallande gren efter toppen (`peak_last` falskt) och vågmönstret inte följer nätet. Byter kollapsformen plats mellan näten (till exempel hals → vägg) används inte Richardson: redovisa ändringen mellan de två finaste och formen per nät. Gränslast där lösaren stannar vid toppen (snäppknäckning) godtas när sista inkrementet är under 0,1 % av förskjutningen och dF/dδ ≈ 0. Projektbladet anger imperfektionsformen; för rotationsskal Koiter axisymmetrisk, inte affin mot linjär knäckform.

## Rapport och spårbarhet

Detaljerade krav och mall finns i skillen fem-rapport. Kärnreglerna:

- Fallmapp: `<projektmapp på Drive>/FEM/<YYMMDD_case>/` med `inputs.json` (projektblad), `src/` (skript, handräkning), `results/` (nät, resultatfiler), `fig/` och `report.md`. Arbetskopia i molnet läggs utanför repot (t.ex. `~/fem-cases/<YYMMDD_case>/`); `cases/` och `results/` är spärrade i `.gitignore`.
- Rapport som PDF (stående A4), filnamn `YYMMDD <namn>.pdf`, när rapporten ska delas. Skisser och diagram som SVG, fältplottar som PNG.
- Lastskiss först (fig. 1): ritad från analyserad geometri med alla laster (storlek, riktning, angreppsyta, lastfall), alla upplag med låsta DOF, symmetriplan, koordinatsystem, mått och reaktionskrafter. Antaganden märks [ANTAGET]. BC-tabell under skissen.
- Figurerna nedan och i skillen fem-rapport (fig. 1–7) är ett minimum, inte ett tak. Agenten lägger själv till de figurer som behövs för att förstå resultatet eller fatta beslutet, till exempel utnyttjandegrad, huvudspänningsriktningar (lastväg), kontakttryck, modformer, känslighetsstudie eller jämförelse före och efter en ändring. Varje extra figur har samma bildtext som övriga plus en mening om varför den behövs.
- Figurerna byggs med byggstenarna i `femtools.report` (L2-projektion, snittuttag, enhetlig stil). Ny figurtyp som blir återkommande läggs in där med test.
- Alla bilder ritas i betraktarvyn (se Visning och betraktarvy). Resultatbilder: deformation (hel modell, förstoringsfaktor angiven), von Mises (skala 0 till Rp0,2, max markerat), detaljvy med nät kring kritisk punkt, minst två snitt (snittplan markerade i översiktsbilden, samma färgskala), konvergensplott. Modform eller bucklingform när relevant.
- Varje resultatbild har bildtext med siffror: vad som visas, var och hur stort (utnyttjandegrad mot gräns), varför, tillförlitlighet (konvergerat, singularitet) och konsekvens.
- Skilj VERIFIERAT från ANTAGET.
- Ange version av denna standard, femtools-version och commit-hash (`femtools.__version__`, `git rev-parse --short HEAD`) och projektbladets innehåll i rapportens början.

## Verktygsregister

Flyttat till `docs/verifiering.md` (v4.0) så att registret ligger bredvid koden och testerna i `femtools/verify.py`. Regeln gäller oförändrat: ett verktyg används i beräkningar först när det är verifierat mot handräkning och infört i registret.

## Erfarenhetsbank

Lärdomar från genomförda beräkningar. Varje post: datum, projekt, lärdom, åtgärd i standarden.

| Datum | Projekt | Lärdom | Åtgärd |
| :--- | :--- | :--- | :--- |
| 2026-11-17 | Tryckkraft skärverktyg | Spänningstopp i skarp slitsrot konvergerar inte. | Krav på verklig radie i utvärderingsområdet. |
| 2026-11-17 | Allmänt | scipy spsolve är enkeltrådad och blir flaskhals för 3D P2. | pypardiso som standardlösare. |
| 2026-11-17 | Allmänt | gmsh i claude.ai-containern kräver systembiblioteket libglu1-mesa (apt-get install libglu1-mesa libxcursor1 libxinerama1 libxft2). | Står i beräkningsmiljön för webben. |
| 2026-11-19 | Allmänt | Rapporter utan lastskiss och snittbilder är svåra att granska och återanvända. | v3.0: obligatorisk lastskiss, resultatbilder och bildtexter. Detaljer i skillen fem-rapport. |
| 2026-11-19 | Allmänt | apt-get install libglu1-mesa misslyckas ("Unable to locate package") i en ny container tills apt-get update körts. | Kör apt-get update före apt-get install i beräkningsmiljön. |
| 2026-11-19 | Allmänt | Nodvärden som medelvärde av elementvärden gav 26,9 MPa i snittet mot 30,0 MPa (balkteori), cirka 10 % för lågt vid ytan. | femtools.report använder L2-projektion till P1. Jämför alltid snittvärden mot handräkning. |
| 2026-11-19 | Allmänt | SolidWorks-materialfiler (.sldmat) från tjänsteleverantörer avvek från tillverkarens datablad: PA 2200 E 2 000 mot 1 650 MPa, PA 3200 GF E 3 600 mot 3 200 MPa, Alumide α 1,2e-5 (fysikaliskt omöjligt, under aluminiums värde). | Materialdata hämtas från tillverkarens datablad och kontrolleras mot minst en annan tillverkare. .sldmat-filerna rättade. |
| 2026-11-26 | Allmänt | gmsh kräver libGLU även i Claude Code-molnet (`OSError: libGLU.so.1`). | Uppstarten i `CLAUDE.md` kör `apt-get update` och `apt-get install` före `pip install`. |
| 2026-11-26 | Allmänt | Registret angav femtools.report 1.1.0 som verifierat, men koden fanns inte i repot (bara 1.0.0). | Registret ligger nu i repot och verify.py måste täcka varje registrerat verktyg. |
| 2026-11-26 | Allmänt | intorder 4 för P2-styvheten gjorde assembleringen 3,4 gånger långsammare än intorder 2 utan noggrannhetsvinst (10,3 mot 3,0 s vid 49 000 DOF). Modal med scipy shift-invert (en tråd) tog 26,7 s mot 6,5 s med pypardiso. | femtools 1.2.0: intorder 2 för K, pypardiso i modal. verify 14 s → 8 s, samma resultat. |
| 2026-11-26 | Allmänt | gmsh HXT med flera trådar ger inte identiskt nät mellan körningar (plåt med hål: 33 651–33 852 DOF, σ_max 32,9–33,4 MPa). | Spara nätet (`.msh`) per nivå och återanvänd det. Krävs exakt reproducerbarhet: `setup(threads=1)`. |
| 2026-11-26 | Allmänt | P2-solidmodell av svetsad rörram: 0,26 M DOF 3,5 GB, 0,50 M DOF 6,7 GB, 0,93 M DOF 11,8 GB (över 70 %-gränsen), 1,4 M DOF överskred 15 GB och containern startades om (alla bakgrundsjobb förlorades). Kompakt axel: 0,14 M DOF 2,5 GB, 0,25 M DOF 4,0 GB. | Minnet växer ungefär 13 GB per miljon DOF för solidmodeller, vilket ger ett tak kring 0,8 M DOF i molnet. Mät de två grövsta näten, extrapolera, och kör aldrig två stora jobb samtidigt. |
| 2026-11-26 | Allmänt | `basis.get_dofs()` utan argument ger bara randens frihetsgrader. Summering av pålagd last med volymkraft blev 10 % fel i jämviktskontrollen. | Summera per komponent över `nodal_dofs[i]` och `edge_dofs[i]` (P2). |
| 2026-11-26 | Allmänt | Rör ritade som solida stänger i CAD (förenklad modell) gav fel massa och styvhet. | Jämför STEP-volym per balk med katalogens kg/m före meshning. Ersätt med verkligt tvärsnitt. |
| 2026-11-26 | Allmänt | Spänningen 4–12 mm från ett inre hörn vid svetsrot växte med förfiningen (111 → 137 MPa) fast den låg utanför hörnet. | Utvärdera svetsdetaljer med linjäriserad spänning per snitt (minsta kvadrat-anpassning a + b·y + c·z) och IIW hot-spot 0,4t/1,0t, inte max i en zon. |
| 2026-11-26 | Allmänt | En fast rapportfunktion med ett givet antal figurer styr agenten mot just de figurerna. | femtools.report byggd som byggstenar; standardfigurerna är minimum och agenten väljer extra figurer med motivering. |
| 2026-11-26 | Allmänt | Fusion 360 exporterar STEP i cm oavsett dokumentets enhet. | Enhetskontroll av STEP före nätning (Meshning, STEP från CAD). |
| 2026-11-26 | Allmänt | gmsh 4.x kunde inte näta en Fusion-solid ("Impossible to mesh periodic surface"); ombyggnad i gmsh OCC med samma radier misslyckades på fillets. Kvartsdelning i CAD + `occ.fragment` fungerade. | Rutin i Meshning, STEP från CAD. |
| 2026-11-26 | Allmänt | Konvergens av max von Mises i ett hörn med R4 sjönk med finare nät (107 → 81 → 70 MPa): grova nät med låg elementkvalitet (γ 0,01–0,03) överskattar. | Avtagande serie: finaste värdet är konservativt. Kontrollera γ i utvärderingsområdet. |
| 2026-11-26 | Allmänt | Toppvärde i integrationspunkt och L2-projicerat nodvärde skiljer upp till 25 % i hörn. | Utvärdera max av båda och redovisa vilket som styr. |
| 2026-11-26 | Allmänt | OpenFOAM v1912 från apt: functionObject `surfaceFieldValue` kraschar ("IOstream sha1"), `#includeFunc` med argument tolkas inte, scotch saknas (decomposePar method hierarchical fungerar). | Efterbehandla Δp och flöde i Python ur polyMesh och tidskataloger; inga function objects i controlDict. |
| 2026-11-26 | Allmänt | STEP → OpenFOAM: gmsh-tetraeder med fysiska grupper inlet/outlet/walls, msh 2.2, `gmshToFoam`, `transformPoints -scale 0.001`; väggpatchen blir typ `patch` och måste ändras till `wall` för väggfunktioner. Inloppshastigheten ska räknas mot nätets facetterade area (1,6 % mindre än cirkeln). | Står i CFD-reglerna i `CLAUDE.md`. |
| 2026-11-26 | Allmänt | Facetterade ytor ger fel resultantkraft om traktion räknas mot teoretisk area (stutsände 179 mot 199 mm², last 180 i stället för 200 N). | Räkna traktion som F / nätets area och kontrollera ΣR per delfall. |
| 2026-11-27 | Allmänt | Nodala fält sparades inte i fallet; den interaktiva visaren krävde omkörning av FEM och CFD. CFD med samma h gav Δp +6–7 % (icke-identiskt nät). Visaren öppnade från väggsidan, inte framifrån. | Fält sparas i varje körning, CFD-nätstudie och sparat nät, betraktarvy i projektbladet (v5.0). |
| 2026-11-27 | Tunnväggigt rotationsskal (dryckesburk) | CalculiX `*BUCKLE` med S8R gav brusformade moder under den klassiska knäcklasten, även med konsistenta randlaster, strukturerat nät och S8. Icke-linjär kollaps på perfekt cylinder gav 2 616 N mot klassiska 2 652 N och Koiter ξ 0,1 +3,3 %. | Kollapslast med icke-linjär analys; `*BUCKLE` för S8R används inte (Knäckning). |
| 2026-11-28 | Tunnväggigt rotationsskal | S8R på ostrukturerat nät med lokal förfining i en hals gav falsk kollaps: lasten sjönk med förfiningen (942 → 684 N), vågantalet följde nätet (37 → 70) och kurvan saknade fallande gren. S8 på samma nät gav "nonpositive jacobian" eller låsning (cylinder +14 % vid h = l/3). Strukturerat nät med S8R gav riktig kollaps som konvergerar (1 599 / 1 389 / 1 300 N vid h 2,25 / 1,5 / 1,0). | Knäckning: strukturerat nät, S8R, krav på fallande gren. `write_inp(quad_type=)`, `orient()` kontrollerar normalernas riktning. |
| 2026-11-27 | Tunnväggigt rotationsskal | Lokal nätförfining som omfattade en hårnålsvikt (två skalytor med motriktade normaler, falsen) stoppade CalculiX direkt: "gen3dnor: size of estimated shell normal ... smaller than 1.e-10". | Förfiningszoner slutar före hårnålar och tätt vikta detaljer; där behåller nätet grundstorleken. |
| 2026-11-27 | Tunnväggigt rotationsskal | CalculiX frd innehöll bara expanderade 3D-noder för skal. | `*NODE FILE, OUTPUT=2D` ger mittytans noder; `femtools.shell.collapse_step` gör det. |
| 2026-11-27 | Allmänt | CalculiX från apt (SPOOLES) tog 2,3 gånger längre tid än samma version byggd med MKL Pardiso. | `tools/build_ccx_pardiso.sh` och `ccx_exe()` (Lösare och prestanda). |
| 2026-11-27 | Allmänt | `material()` hanterade inte typvärden för metaller med tjockleksintervall. | femtools 1.5.0: tjockleksuppslag även för TYPVARDE-metaller. |

## Revisionsprocess

- Efter varje avslutat beräkningsfall föreslår agenten ändringar i standarden: nya lärdomar, nya verifierade verktyg, nya material. Förslaget redovisas som en lista i slutet av rapporten.
- Ändringen görs som pull request i fem-tools. Marcus granskar och godkänner genom att merga. Versionsnumret höjs och revisionsloggen fylls i i samma pull request.
- Mindre ändring (ny post i erfarenhetsbank, verktygsregister eller materialdata) höjer andra siffran, t.ex. 4.0 till 4.1. Ändring av regler eller kriterier höjer första siffran.
- Lärdomar skrivs generiskt (typ av detalj och fenomen), aldrig med kundnamn, artikelnummer eller ritningsnummer.
- Efter merge kopieras filen till Drive-mappen "Projektstandard för FEM-beräkningar". Föregående version flyttas till `_arkiv/` med version och datum i filnamnet.
