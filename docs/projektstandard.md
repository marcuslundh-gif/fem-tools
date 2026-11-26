# Projektstandard för FEM-beräkningar

Version 3.2, 2026-11-19. Ägare: Marcus Lundh.

## Revisionslogg

| Version | Datum | Ändring |
| :--- | :--- | :--- |
| 1.0 | 2026-11-12 | Första version: enheter, två material, randvillkor, acceptanskriterier. |
| 2.0 | 2026-11-17 | Säkerhetsfaktor flyttad till projektblad per projekt. Materialdata tjockleksberoende och med HAZ. Nya avsnitt: meshning, lösare och prestanda, verifiering, utmattning/egenfrekvens/knäckning, verktygsregister, erfarenhetsbank och revisionsprocess. RBE2/RBE3 översatt till scikit-fem. |
| 3.0 | 2026-11-19 | Rapportkrav skärpta: lastskiss med randvillkor först i rapporten, obligatoriska resultatbilder (deformation, spänning, snitt, konvergens) med förklarande bildtext. Detaljerade krav flyttade till skillen fem-rapport, standarden har kärnreglerna. |
| 3.1 | 2026-11-19 | femtools.report 1.1.0 verifierad och införd i verktygsregistret (öppen punkt stängd). Ny erfarenhetspost om apt-get update. |
| 3.2 | 2026-11-19 | Värmeutvidgning α införd för alla material. Nya material med typvärden (ej minimivärden): PA12 SLS (EOS PA 2200), PA12 GB SLS (EOS PA 3200 GF), PA12 MJF (HP 3D HR PA 12), Alumide. Värden kontrollerade mot EOS, HP, Materialise och 3D Systems. |

## Användning

Standarden gäller alla FEM- och hållfasthetsberäkningar i projektet "FEM Agent", lokalt i Claude Code och på claude.ai. Projektspecifika värden (säkerhetsfaktorer, lastfall, krav) står i projektbladet för varje beräkningsfall, inte här. Om projektbladet och standarden säger olika gäller projektbladet för just det projektet. Standarden revideras enligt avsnittet Revisionsprocess.

## Enhetssystem

All beräkning och redovisning sker i mm, N, MPa, ton, s. Densitet anges i ton/mm³, vilket ger egenfrekvenser direkt i Hz. Kontrollera STEP-filens enhet mot bounding box innan meshning.

## Projektblad (per projekt)

Varje beräkningsfall har ett projektblad, `inputs/<case>.json`, som fylls i innan modellering. Saknas ett värde ska agenten fråga. Inga standardvärden antas för säkerhetsfaktorer eller laster.

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

Värden är minimivärden enligt angiven standard. Använd materialcertifikat när det finns. Nytt material läggs till med källa.

| Material | Standard och produkt | Tjocklek t [mm] | E [MPa] | ν [-] | Rp0,2 / ReH [MPa] | Rm [MPa] | ρ [ton/mm³] | α [1/K] |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| Stål S355 | EN 10025-2 | t ≤ 16 | 210 000 | 0,30 | 355 | 470–630 | 7,85e-9 | 1,2e-5 |
| Stål S355 | EN 10025-2 | 16 < t ≤ 40 | 210 000 | 0,30 | 345 | 470–630 | 7,85e-9 | 1,2e-5 |
| Stål S355 | EN 10025-2 | 40 < t ≤ 63 | 210 000 | 0,30 | 335 | 470–630 | 7,85e-9 | 1,2e-5 |
| Alu 6082-T6 | EN 1999-1-1 tab 3.2b, strängpressat | t ≤ 5 | 70 000 | 0,33 | 250 | 290 | 2,70e-9 | 2,3e-5 |
| Alu 6082-T6 | EN 1999-1-1 tab 3.2b, strängpressat | 5 < t ≤ 15 | 70 000 | 0,33 | 260 | 310 | 2,70e-9 | 2,3e-5 |
| Alu 6082-T6, svetspåverkad zon (HAZ) | EN 1999-1-1 tab 3.2b, MIG, t ≤ 15 | t ≤ 15 | 70 000 | 0,33 | 125 | 185 | 2,70e-9 | 2,3e-5 |

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

Standardinställningar i gmsh (finns i `mesh_template.py`):

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

Kvalitetskrav innan lösning:

- Elementkvalitet gamma (gmsh) minst 0,1 i hela nätet och minst 0,3 i utvärderingsområdet. Sämre element redovisas med position.
- Inga lösa noder eller dubbla volymer. Kontrollera att nätets volym stämmer med CAD-volymen inom 0,5 %.

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
- Praktiskt tak i claude.ai-containern (2 kärnor, 7 GB uppmätt 2026-11-17): cirka 300 000 DOF. Större modeller körs lokalt.

## Verifiering

Obligatoriskt för varje beräkningsfall:

- Jämvikt: summan av reaktionskrafter plus pålagda laster ska vara under 0,1 % av pålagd last, per komponent.
- Nätkonvergens: minst 3 nät. Godkänt när utvärderingsstorheten ändras mindre än 3 % mellan de två finaste. Redovisa tabell och Richardson-extrapolerat värde.
- Handräkning eller analytisk lösning för storleksordning, i ett eget skript. FEM och handräkning ska stämma inom cirka 10 % i en jämförbar punkt, annars förklaras skillnaden.
- Rimlighet: deformationsform (skalad plott), tecken och var maxvärden hamnar.
- Singularitet: en topp som växer för varje förfining är en singularitet och rapporteras inte som spänning. Utvärdera med verklig radie, linjariserad spänning eller hot-spot-metod och ange vilken.

## Acceptanskriterier

- Säkerhetsfaktorer enligt projektbladet. Saknas de beräknas inget slutresultat, agenten frågar.
- Max deformation: L/300 om projektbladet inte anger annat.
- Spänningsutvärdering: von Mises för segt stål och aluminium. Lokala singulariteter vid vassa inre hörn ignoreras om de inte är utmattningskritiska. Är de utmattningskritiska utvärderas de med hot-spot- eller kerbspänningsmetod.
- Flagga alltid: von Mises över sträckgränsen (linjär analys ogiltig lokalt), knäckning som kan vara dimensionerande, stora deformationer (över cirka 1/10 av tjockleken).

## Utmattning, egenfrekvens och knäckning

- Utmattning, grundmaterial: Kf = 1 + q(Kt − 1). Goodman eller Gerber med korrigerad utmattningsgräns (yta, storlek, tillförlitlighet). Säkerhetsfaktor mot utmattning redovisas separat från statisk.
- Utmattning, svetsar: nominell spänning eller hot-spot-spänningsvidd mot FAT-klass enligt EN 1993-1-9, EN 1999-1-3 eller IIW. Variabel amplitud med Palmgren-Miner. Ange FAT-klass och källa.
- Egenfrekvens: modalanalys med verkliga randvillkor och tillsatta massor. Marginal mot drivfrekvenser och deras övertoner anges i projektbladet. Riktvärde om inget anges: minst 20 %.
- Knäckning: linjär buckling-analys ger lastfaktor. Tunnväggiga eller imperfektionskänsliga konstruktioner kräver icke-linjär analys eller standardens reduktionskurvor.

## Rapport och spårbarhet

Detaljerade krav och mall finns i skillen fem-rapport. Kärnreglerna:

- Rapport i `results/<YYMMDD_case>/report.md`, figurer i `fig/`. PDF (stående A4) när rapporten ska delas. Skisser och diagram som SVG, fältplottar som PNG.
- Lastskiss först (fig. 1): ritad från analyserad geometri med alla laster (storlek, riktning, angreppsyta, lastfall), alla upplag med låsta DOF, symmetriplan, koordinatsystem, mått och reaktionskrafter. Antaganden märks [ANTAGET]. BC-tabell under skissen.
- Resultatbilder: deformation (hel modell, förstoringsfaktor angiven), von Mises (skala 0 till Rp0,2, max markerat), detaljvy med nät kring kritisk punkt, minst två snitt (snittplan markerade i översiktsbilden, samma färgskala), konvergensplott. Modform eller bucklingform när relevant.
- Varje resultatbild har bildtext med siffror: vad som visas, var och hur stort (utnyttjandegrad mot gräns), varför, tillförlitlighet (konvergerat, singularitet) och konsekvens.
- Skilj VERIFIERAT från ANTAGET.
- Ange version av denna standard, skriptversion (commit-hash när repo finns) och projektbladets innehåll i rapportens början.

## Verktygsregister

Verktyg läggs till här först när de är verifierade mot handräkning. Ange version, verifieringsfall och resultat.

| Verktyg | Version | Verifieringsfall | Resultat |
| :--- | :--- | :--- | :--- |
| fem_template.py (scikit-fem 12.0.2, P2, pypardiso) | 1.0, 2026-11-17 | Konsol 100×10×10 mm, 100 N | Nedböjning 0,1904 mm mot 0,1905 mm (balkteori). σ vid L/2 30,3 mot 30,0 MPa. Jämvikt exakt. |
| modalanalys i fem-berakning-skillen | 1.0, 2026-11-17 | Samma konsol, stål | f1 834 Hz mot 836 Hz (Euler-Bernoulli). |
| mesh_template.py (gmsh 4.15.2, HXT, lokal förfining) | 1.0, 2026-11-17 | Plåt 100×40×5 mm med hål Ø10, 10 MPa drag | σ_max 33,1–33,3 MPa på tre nät (ändring 0,6 %) mot 32,3 MPa (Heywood). Skillnaden 3 % tillskrivs 3D-effekt genom tjockleken. Meshtid under 0,2 s. |
| femtools.report (matplotlib 3.11) | 1.1.0, 2026-11-19 | Samma konsol, snitt vid x = L/2, nät 48×6×6 | Nedböjning 0,1910 mot 0,1905 mm. σ i snitt 30,0 mot 30,0 MPa (L2-projektion till P1). Fem figurer skapas, ingår i python -m femtools.verify. Ej testat: plåt med hål, STEP-geometri, nät över cirka 100 000 element (snittfunktionen loopar i Python). |

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

## Revisionsprocess

- Efter varje avslutat beräkningsfall föreslår agenten ändringar i standarden: nya lärdomar, nya verifierade verktyg, nya material. Förslaget redovisas som en lista i slutet av rapporten.
- Marcus godkänner. Först därefter uppdateras dokumentet, versionsnumret höjs och revisionsloggen fylls i.
- Mindre ändring (ny post i erfarenhetsbank eller verktygsregister) höjer andra siffran, t.ex. 2.0 till 2.1. Ändring av regler eller kriterier höjer första siffran.
- Samma version ska finnas i Google Drive-mappen "Projektstandard för FEM-beräkningar" och som projektdokument i claude.ai-projektet "FEM Agent".
