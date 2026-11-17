# femtools

Verifierade FEM-mallar för FEM-agenten: gmsh-meshning (HXT, lokal förfining, kvalitetskontroll)
och scikit-fem P2-statik/modal med pypardiso. Enheter mm, N, MPa, t, s.

Installera (claude.ai-container: först `apt-get install -y libglu1-mesa libxcursor1 libxinerama1 libxft2`):

    pip install --break-system-packages git+https://github.com/<konto>/fem-tools

Verifiera:

    python -m femtools.verify

Användning:

    from femtools import mesh_step, static, modal, richardson
