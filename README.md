# femtools

Verifierade FEM-verktyg: gmsh-meshning (HXT, lokal förfining, kvalitetskontroll), scikit-fem P2-statik (traktion, tryck, egenvikt, rotation, termisk last), modalanalys, linjär knäckning och värmeledning (stationär och transient) med pypardiso, materialdata med källor och projektblad med stoppkontroll. Enheter mm, N, MPa, t, s.

Repot är master för reglerna: [`docs/projektstandard.md`](docs/projektstandard.md) och registret över verifierade verktyg [`docs/verifiering.md`](docs/verifiering.md). Agentinstruktioner: [`CLAUDE.md`](CLAUDE.md).

Installera (gmsh kräver först `apt-get update && apt-get install -y libglu1-mesa libxcursor1 libxinerama1 libxft2` på Linux):

    pip install git+https://github.com/marcuslundh-gif/fem-tools

Verifiera:

    python -m femtools.verify

Användning:

    from femtools import mesh_step, static, modal, richardson, material, projektblad, report
    pb = projektblad.load("inputs.json")
    mat = material(pb["material"], t=pb["tjocklek_mm"])
