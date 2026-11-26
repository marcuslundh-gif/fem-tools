# Verktygsregister

Ett verktyg används i beräkningar först när det är verifierat mot handräkning och infört här. Varje rad ska ha ett test i `femtools/verify.py`. Ändringar görs som pull request och höjer andra siffran i standarden (`docs/projektstandard.md`).

| Verktyg | Version | Verifieringsfall | Resultat | Test |
| :--- | :--- | :--- | :--- | :--- |
| `femtools.static` (scikit-fem 12.0.2, P2, pypardiso, K med intorder 2) | 1.2.0, 2026-11-26 | Konsol 100×10×10 mm, 100 N | Nedböjning 0,1904 mot 0,1905 mm (balkteori), oförändrad mot 1.0.0. σ vid L/2 30,3 mot 30,0 MPa (1.0.0). Jämvikt 5e-10. | `cantilever()` |
| `femtools.modal` | 1.2.0, 2026-11-26 (shift-invert med pypardiso) | Samma konsol, stål | f1 833 Hz mot 836 Hz (Euler-Bernoulli), samma som 1.0.0. | `cantilever()` |
| `femtools.mesh` (gmsh 4.15.2, HXT, lokal förfining) | 1.0.0, 2026-11-17 | Plåt 100×40×5 mm med hål Ø10, 10 MPa drag | σ_max 32,9–33,4 MPa (HXT-nätet varierar mellan körningar) mot 32,3 MPa (Heywood). Skillnaden 3 % tillskrivs 3D-effekt genom tjockleken. Meshtid under 0,2 s. | `plate_hole()` |
| `femtools.richardson` | 1.0.0 | Formel enligt fem-berakning | Ej eget test. | saknas |
| `femtools.material`, `materials.json` | 1.2.0, 2026-11-26 | Tjockleksuppslag S355/6082-T6, Z-riktning PA 2200 | Värden lika med standardens tabell. | `data()` |
| `femtools.report.nodal` + `sample` (L2-projektion till P1, punktlokalisering) | 1.2.0, 2026-11-26 | Konsol 48×6×6 och 96×6×6 tetraedrar, snitt genom tjockleken vid L/2 | σ_xx 30,03–30,06 mot 30,00 MPa (balkteori), max avvikelse 0,2 % över snittet. Plåt med hål: σ_xx i hålkant 32,6 MPa mot 32,3 MPa (Heywood). 213 000 DOF: nodal 1,5 s, plansnitt 0,5 s. Ersätter report 1.1.0 från claude.ai, vars kod aldrig kom in i repot. | `cantilever()`, `plate_hole()` |
| `femtools.report.plot_*` (matplotlib 3.11) | 1.2.0, 2026-11-26 | Deformation, fält, plansnitt, linje skapas utan fel | Visuellt granskade på konsolen. `plot_bc` och `plot_convergence` granskade men ej i test. | `cantilever()` |
| `femtools.static_ext` – tryck (−p·n), randvillkor per komponent (symmetri) | 1.3.0, 2026-11-26 | Kvarts tjockväggigt rör a 20, b 40, H 8 mm, plan töjning, p 10 MPa inne, ca 20 000 DOF | σ_t,inner 16,51–16,58 mot 16,67 MPa (Lamé), högst −1,0 % (HXT-nätet varierar mellan körningar; P2 på raka kanter mot krökt yta). Jämvikt F + R < 1e-6. | `loads_test()` |
| `femtools.static_ext` – rotation | 1.3.0, 2026-11-26 | Samma rör, stål, 6000 rpm kring z | σ_t,inner 4,40–4,41 mot 4,43 MPa (roterande skiva, ν' = ν/(1−ν)), −0,5 %. Summa F_x 462,6 mot 462,8 N. | `loads_test()` |
| `femtools.static_ext` – egenvikt | 1.3.0, 2026-11-26 | Hängande stång 1000×10×10 mm, stål | σ vid infästning 0,0770 mot 0,0770 MPa (ρgL), förlängning 0,1834 mot 0,1834 µm (ρgL²/2E). | `loads_test()` |
| `femtools.static_ext` – termisk last (α, ΔT konstant eller P2-fält) | 1.3.0, 2026-11-26 | Inspänd stång ΔT 50 K; rör med logaritmisk temperatur från `heat` | Stång −126,0 mot −126,0 MPa (−EαΔT). Rör σ_t inner −5,97 mot −5,94 MPa, ytter 3,77 mot 3,77 MPa (Timoshenko & Goodier, plan töjning), under 0,6 %. | `thermal_test()` |
| `femtools.heat` – stationär värmeledning, fast T, konvektion | 1.3.0, 2026-11-26 | Rör ovan, T_inner 120 °C, konvektion 50 W/(m²·K) mot 20 °C | Värmeflöde 2445,2 mot 2445,5 mW (seriekoppling ledning + konvektion), T_ytter 117,30 mot 117,30 °C. Energibalans 2e-13. Flöde (`flux`) och värmekälla (`Q`) ingår i balansen men saknar eget analytiskt fall. | `thermal_test()` |
| `femtools.heat_transient` – implicit Euler / theta-metod | 1.3.0, 2026-11-26 | Platta 20 mm, båda sidor steg till 100 °C, Fourier-tal 0,1, 800 steg | T i mitten 52,52 mot 52,55 °C (Fourierserie). 200 steg: 52,44 °C. | `thermal_test()` |
| `femtools.projektblad` | 1.2.0, 2026-11-26 | Tom mall ska stoppas, ifyllt blad ska passera | Mallen ger 15 saknade fält, ifyllt blad 0. | `data()` |

## Körningar

| Datum | Miljö | femtools | Resultat |
| :--- | :--- | :--- | :--- |
| 2026-11-17 | claude.ai-container, 2 kärnor, 7 GB | 1.0.0 | OK |
| 2026-11-26 | Claude Code-molnet, 4 kärnor, 15 GB, Python 3.13 | 1.0.0-kod | OK på 13–14 s |
| 2026-11-26 | Samma | 1.2.0 | OK på 8–10 s (konsol 0,1904/0,1905 mm, f1 833/836 Hz, snitt 30,06/30,00 MPa, plåt 33,1/32,3 MPa) |
| 2026-11-26 | Claude Code-molnet, 4 kärnor, 15 GB | 1.3.0 | OK på 15 s (nya fall: laster och termik, se ovan) |
