"""Gränssnitt mot OpenFOAM (ESI, v1912 i molnet via apt, v2406 på smartfish i WSL) för turbulent CFD.

OpenFOAM körs i SI-enheter (m, s, kinematisk viskositet m²/s, kinematiskt tryck p/ρ m²/s²), till skillnad
från resten av femtools. Omräkning görs i den här modulen.

    from femtools import openfoam
    openfoam.available()                     # versionssträng eller None
    log = openfoam.run(case, "blockMesh", "simpleFoam")
    r = openfoam.turbulent_pipe("/tmp/rör", Re=1e5)     # verifieringsfall: utvecklad rörströmning, k-ω SST

Uppstart i molnet: `apt-get install -y openfoam` (paketet installerar lösarna i /usr/bin).
På smartfish: kör från Windows med `wsl -e bash -c 'source /usr/lib/openfoam/openfoam2406/etc/bashrc && ...'`
och lägg fallen i Linux-hemmet (~/fem-cases), inte på /mnt/c.
Kör tunga fall parallellt: decomposePar, `mpirun -np N <lösare> -parallel`, reconstructPar.
"""
import os
import re
import shutil
import subprocess
import numpy as np

_BASHRC = ["/usr/share/openfoam/etc/bashrc", *sorted(
    (f"/usr/lib/openfoam/{d}/etc/bashrc" for d in (os.listdir("/usr/lib/openfoam") if os.path.isdir("/usr/lib/openfoam") else [])),
    reverse=True)]

HEADER = """FoamFile
{{
    version     2.0;
    format      ascii;
    class       {cls};
    object      {obj};
}}
"""


def _env_prefix():
    for b in _BASHRC:
        if os.path.isfile(b):
            return f"source {b} >/dev/null 2>&1; "
    return ""


def available():
    """OpenFOAM-version (WM_PROJECT_VERSION) om lösarna finns, annars None."""
    try:
        out = subprocess.run(["bash", "-c", _env_prefix() + 'command -v simpleFoam >/dev/null && echo "$WM_PROJECT_VERSION"'],
                             capture_output=True, text=True, timeout=60).stdout.strip()
    except (OSError, subprocess.TimeoutExpired):
        return None
    return out or None


def write(case, files):
    """Skriv fall: files = {"system/controlDict": (klass, innehåll), ...}. Rubriken läggs till automatiskt."""
    for path, (cls, body) in files.items():
        p = os.path.join(case, path)
        os.makedirs(os.path.dirname(p), exist_ok=True)
        with open(p, "w") as fh:
            fh.write(HEADER.format(cls=cls, obj=os.path.basename(path)) + body)


def run(case, *cmds, timeout=3600):
    """Kör OpenFOAM-kommandon i fallmappen. Loggar till log.<kommando>. Stoppar vid fel."""
    logs = {}
    for c in cmds:
        name = c.split()[0]
        log = os.path.join(case, f"log.{name}")
        r = subprocess.run(["bash", "-c", f"{_env_prefix()}cd '{case}' && {c} > 'log.{name}' 2>&1"],
                           timeout=timeout)
        txt = open(log).read()
        if r.returncode != 0 or "FOAM FATAL" in txt:
            raise RuntimeError(f"{c} misslyckades, se {log}:\n" + "\n".join(txt.splitlines()[-20:]))
        logs[name] = txt
    return logs


def colebrook(Re, eps_D=0.0):
    """Darcy-friktionsfaktor (Colebrook-White; eps_D = 0 ger Prandtls lag för släta rör)."""
    f = 0.02
    for _ in range(100):
        f = (-2.0 * np.log10(eps_D / 3.7 + 2.51 / (Re * np.sqrt(f))))**-2
    return f


def turbulent_pipe(case, Re=1e5, n_r=80, wall_ratio=0.005, angle=5.0, iters=4000):
    """Fullt utvecklad turbulent rörströmning: kilformad (axisymmetrisk) sektor, periodisk i x, drivs av
    meanVelocityForce (Ubar 1 m/s), k-ω SST med upplöst vägg (y+ < 1). D = 1 m, nu = 1/Re.
    Returnerar dict med f (Darcy), f_colebrook, y+ första cell, iterationer."""
    D, U, L = 1.0, 1.0, 0.1
    R, nu = D / 2, U * D / Re
    a = np.radians(angle / 2)
    c, s = np.cos(a), np.sin(a)
    if os.path.isdir(case):
        shutil.rmtree(case)
    k0 = 1.5 * (0.05 * U)**2
    om0 = np.sqrt(k0) / (0.09**0.25 * 0.07 * D)
    files = {
        "system/blockMeshDict": ("dictionary", f"""
convertToMeters 1;
vertices ( (0 0 0) ({L} 0 0) ({L} {R * c} {-R * s}) (0 {R * c} {-R * s}) ({L} {R * c} {R * s}) (0 {R * c} {R * s}) );
blocks ( hex (0 1 2 3 0 1 4 5) (2 {n_r} 1) simpleGrading (1 {wall_ratio} 1) );
edges ( arc 3 5 (0 {R} 0) arc 2 4 ({L} {R} 0) );
boundary
(
    front {{ type wedge; faces ((0 3 2 1)); }}
    back  {{ type wedge; faces ((0 1 4 5)); }}
    wall  {{ type wall;  faces ((3 5 4 2)); }}
    inlet  {{ type cyclic; neighbourPatch outlet; faces ((0 5 3 0)); }}
    outlet {{ type cyclic; neighbourPatch inlet;  faces ((1 2 4 1)); }}
    axis  {{ type empty; faces ((0 1 1 0)); }}
);
"""),
        "system/controlDict": ("dictionary", f"""
application simpleFoam; startFrom startTime; startTime 0; stopAt endTime; endTime {iters}; deltaT 1;
writeControl timeStep; writeInterval {iters}; purgeWrite 1; writeFormat ascii; writePrecision 10;
runTimeModifiable false;
"""),
        "system/fvSchemes": ("dictionary", """
ddtSchemes { default steadyState; }
gradSchemes { default Gauss linear; }
divSchemes { default none; div(phi,U) bounded Gauss linearUpwind grad(U); div(phi,k) bounded Gauss upwind;
             div(phi,omega) bounded Gauss upwind; div((nuEff*dev2(T(grad(U))))) Gauss linear; }
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
wallDist { method meshWave; }
"""),
        "system/fvSolution": ("dictionary", """
solvers
{
    p { solver GAMG; smoother GaussSeidel; tolerance 1e-9; relTol 0.01; }
    "(U|k|omega)" { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-10; relTol 0.1; }
}
SIMPLE { nNonOrthogonalCorrectors 0; consistent yes; pRefCell 0; pRefValue 0; residualControl { p 1e-7; U 1e-8; "(k|omega)" 1e-8; } }
relaxationFactors { equations { U 0.9; ".*" 0.7; } }
"""),
        "constant/fvOptions": ("dictionary", f"""
momentumSource
{{
    type meanVelocityForce; active yes;
    meanVelocityForceCoeffs {{ selectionMode all; fields (U); Ubar ({U} 0 0); }}
}}
"""),
        "constant/transportProperties": ("dictionary", f"transportModel Newtonian;\nnu {nu};\n"),
        "constant/turbulenceProperties": ("dictionary",
                                          "simulationType RAS;\nRAS { RASModel kOmegaSST; turbulence on; printCoeffs off; }\n"),
    }
    cyc = "inlet { type cyclic; } outlet { type cyclic; } front { type wedge; } back { type wedge; } axis { type empty; }"
    fields = {
        "0/U": ("volVectorField", f"dimensions [0 1 -1 0 0 0 0];\ninternalField uniform ({U} 0 0);\n"
                                  f"boundaryField {{ wall {{ type noSlip; }} {cyc} }}\n"),
        "0/p": ("volScalarField", "dimensions [0 2 -2 0 0 0 0];\ninternalField uniform 0;\n"
                                  f"boundaryField {{ wall {{ type zeroGradient; }} {cyc} }}\n"),
        "0/k": ("volScalarField", f"dimensions [0 2 -2 0 0 0 0];\ninternalField uniform {k0};\n"
                                  f"boundaryField {{ wall {{ type fixedValue; value uniform 1e-12; }} {cyc} }}\n"),
        "0/omega": ("volScalarField", f"dimensions [0 0 -1 0 0 0 0];\ninternalField uniform {om0};\n"
                                      f"boundaryField {{ wall {{ type omegaWallFunction; value uniform {om0}; }} {cyc} }}\n"),
        "0/nut": ("volScalarField", "dimensions [0 2 -1 0 0 0 0];\ninternalField uniform 0;\n"
                                    f"boundaryField {{ wall {{ type nutLowReWallFunction; value uniform 0; }} {cyc} }}\n"),
    }
    files.update(fields)
    write(case, files)
    logs = run(case, "blockMesh", "checkMesh", "simpleFoam")
    log = logs["simpleFoam"]
    g = re.findall(r"pressure gradient = ([-\deE.+]+)", log)
    if not g:
        raise RuntimeError("Hittar inget tryckgradientvärde i log.simpleFoam (meanVelocityForce).")
    dpdx = float(g[-1])                                          # kinematiskt, m/s² (= Δp/(ρL))
    f = dpdx * D / (0.5 * U**2)
    it = re.findall(r"^Time = (\d+)", log, re.M)
    # y+ i första cellens mitt från gradering (functionObject yPlus ger IO-fel i Debian-paketet v1912)
    r_ = wall_ratio**(1 / (n_r - 1))
    d1 = R * (1 - r_) / (1 - r_**n_r)                          # cellstorlek vid axeln
    y1 = 0.5 * d1 * r_**(n_r - 1)                              # halva väggcellen
    yplus = np.sqrt(dpdx * D / 4) * y1 / nu                    # u_tau = sqrt(tau_w/rho) = sqrt(dpdx·D/4)
    return dict(f=f, f_colebrook=colebrook(Re), Re=Re, dpdx=dpdx, iterations=int(it[-1]) if it else None,
                yplus=yplus, converged="SIMPLE solution converged" in log)


# ---------------------------------------------------------------------------- efterbehandling (ascii, tetraedernät)
def _foam_list(path):
    """Första listan i en OpenFOAM-fil (polyMesh/points, faces, owner, neighbour) som text efter '('."""
    import re
    t = open(path).read()
    m = re.search(r"^(\d+)\s*\n\(", t, re.M)
    return int(m.group(1)), t[m.end():]


def _internal_field(path, vec):
    import re
    t = open(path).read()
    m = re.search(r"internalField\s+nonuniform\s+List<\w+>\s*(\d+)\s*\(", t)
    if m is None:
        u = re.search(r"internalField\s+uniform\s+(\([^)]*\)|[-\deE.+]+)", t).group(1)
        return None, np.array(u.strip("()").split(), float)
    k, body = int(m.group(1)), t[m.end():]
    if vec:
        a = np.array(re.findall(r"\(([-\deE.+]+) ([-\deE.+]+) ([-\deE.+]+)\)", body)[:k], float)
    else:
        a = np.array(body[:body.index(")")].split()[:k], float)
    return k, a


def read_case(case, time=None, rho=1.0, scale=1000.0, wall="walls"):
    """Läs ett OpenFOAM-fall med rena tetraedrar (gmshToFoam) och writeFormat ascii för visning och efterbehandling.

    Returnerar dict med P (punkter × scale, standard mm), tet (ncell, 4), Uc och pc (cellvärden, p × rho = Pa för
    simpleFoam), U och p (nodvärden = medel av anslutande celler; U = 0 i väggnoder, noSlip), patches
    {namn: (nf, 3) trianglar} och time. Används av webreport.Report.flow."""
    import re
    pm = os.path.join(case, "constant", "polyMesh")
    n, b = _foam_list(os.path.join(pm, "points"))
    P = np.array(re.findall(r"\(([-\deE.+]+) ([-\deE.+]+) ([-\deE.+]+)\)", b)[:n], float) * scale
    n, b = _foam_list(os.path.join(pm, "faces"))
    F = [list(map(int, f.split())) for f in re.findall(r"\d+\(([\d ]+)\)", b)[:n]]
    if any(len(f) != 3 for f in F):
        raise ValueError("read_case kräver rena tetraedernät (alla ytor trianglar)")
    F = np.array(F)
    n, b = _foam_list(os.path.join(pm, "owner")); own = np.array(b[:b.index(")")].split()[:n], int)
    n, b = _foam_list(os.path.join(pm, "neighbour")); nei = np.array(b[:b.index(")")].split()[:n], int)
    bnd = open(os.path.join(pm, "boundary")).read()
    patches = {}
    for name, blk in re.findall(r"(\w+)\s*\{([^}]*)\}", bnd[bnd.index("(") + 1:]):
        nf = re.search(r"nFaces\s+(\d+)", blk); sf = re.search(r"startFace\s+(\d+)", blk)
        if nf and sf:
            s0 = int(sf.group(1)); patches[name] = F[s0:s0 + int(nf.group(1))]
    if time is None:
        ts = [d for d in os.listdir(case) if re.fullmatch(r"[\d.eE+-]+", d) and float(d) > 0
              and os.path.exists(os.path.join(case, d, "U"))]
        time = max(ts, key=float)
    _, Uc = _internal_field(os.path.join(case, str(time), "U"), True)
    _, pc = _internal_field(os.path.join(case, str(time), "p"), False)
    ncell = len(Uc)
    cells = np.concatenate([own, nei]); fp = np.concatenate([F, F[:len(nei)]])
    pairs = np.unique(np.column_stack([np.repeat(cells, 3), fp.ravel()]), axis=0)
    if len(pairs) != 4 * ncell:
        raise ValueError("read_case: cellerna är inte tetraedrar")
    tet = pairs[:, 1].reshape(ncell, 4)
    cnt = np.maximum(np.bincount(tet.ravel(), minlength=len(P)), 1)
    U = np.column_stack([np.bincount(tet.ravel(), weights=np.repeat(Uc[:, k], 4), minlength=len(P)) / cnt for k in range(3)])
    pc = pc * rho
    p = np.bincount(tet.ravel(), weights=np.repeat(pc, 4), minlength=len(P)) / cnt
    if wall in patches:                               # noSlip i väggnoder, utom kanten mot andra patchar
        other = [v for k, v in patches.items() if k != wall]
        keep = np.unique(np.concatenate([o.ravel() for o in other])) if other else np.array([], int)
        U[np.setdiff1d(np.unique(patches[wall]), keep)] = 0.0
    return dict(P=P, tet=tet, Uc=Uc, pc=pc, U=U, p=p, patches=patches, time=str(time))
