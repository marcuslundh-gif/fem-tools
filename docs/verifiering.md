# Verktygsregister

Ett verktyg används i beräkningar först när det är verifierat mot handräkning och infört här. Varje rad ska ha ett test i `femtools/verify.py`. Ändringar görs som pull request och höjer andra siffran i standarden (`docs/projektstandard.md`).

| Verktyg | Version | Verifieringsfall | Resultat | Test |
| :--- | :--- | :--- | :--- | :--- |
| `femtools.static` (scikit-fem 12.0.2, P2, pypardiso) | 1.0.0, 2026-11-17 | Konsol 100×10×10 mm, 100 N | Nedböjning 0,1904 mot 0,1905 mm (balkteori). σ vid L/2 30,3 mot 30,0 MPa. Jämvikt 1e-10. | `cantilever()` |
| `femtools.modal` | 1.0.0, 2026-11-17 | Samma konsol, stål | f1 833–834 Hz mot 836 Hz (Euler-Bernoulli). | `cantilever()` |
| `femtools.mesh` (gmsh 4.15.2, HXT, lokal förfining) | 1.0.0, 2026-11-17 | Plåt 100×40×5 mm med hål Ø10, 10 MPa drag | σ_max 33,1–33,4 MPa mot 32,3 MPa (Heywood). Skillnaden 3 % tillskrivs 3D-effekt genom tjockleken. Meshtid under 0,2 s. | `plate_hole()` |
| `femtools.richardson` | 1.0.0 | Formel enligt fem-berakning | Ej eget test. | saknas |
| `femtools.material`, `materials.json` | 1.2.0, 2026-11-26 | Tjockleksuppslag S355/6082-T6, Z-riktning PA 2200 | Värden lika med standardens tabell. | `data()` |
| `femtools.projektblad` | 1.2.0, 2026-11-26 | Tom mall ska stoppas, ifyllt blad ska passera | Mallen ger 15 saknade fält, ifyllt blad 0. | `data()` |

## Ej i repot

| Verktyg | Status | Åtgärd |
| :--- | :--- | :--- |
| `femtools.report` 1.1.0 (matplotlib 3.11) | Registrerad i standard v3.1 (konsol, snitt vid L/2: 0,1910 mot 0,1905 mm, σ 30,0 mot 30,0 MPa, L2-projektion till P1), men koden finns inte i repot. | Lägg in koden med test i `verify.py`, eller skriv ny och verifiera om. Används inte förrän dess. |

## Körningar

| Datum | Miljö | femtools | Resultat |
| :--- | :--- | :--- | :--- |
| 2026-11-17 | claude.ai-container, 2 kärnor, 7 GB | 1.0.0 | OK |
| 2026-11-26 | Claude Code-molnet, 4 kärnor, 15 GB, Python 3.13 | 1.2.0 | OK på 14 s (konsol 0,1904/0,1905 mm, f1 833/836 Hz, plåt 33,4/32,3 MPa) |
