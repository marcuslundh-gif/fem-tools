"""Tunnväggiga rotationsskal med CalculiX (ccx, S8R/S6): nät, linjär knäckning och olinjär kollaps.

Enheter mm, N, MPa. Kräver `apt-get install calculix-ccx` (ccx 2.21 verifierad).
Arbetsgång:
    nodes, el, info = revolve_shell(chains, h)          # gmsh, 8-nods fyrhörningar (+ 6-nods trianglar vid axeln)
    write_inp(...) med buckle_step() eller collapse_step()
    run(job, cwd)  ->  read_buckle(), read_frd(), read_totals()
Verifiering: femtools.verify.shell_test() (cylinder mot klassisk knäckspänning 0,605·E·t/R).
"""
import os, re, subprocess, time
import numpy as np


# ---------------------------------------------------------------- nät
def revolve_shell(chains, h, h_func=None, quad=True, structured=False, zones=None):
    """Nät av mittytan från rotationsprofiler.

    chains: lista av (N, 2)-arrayer med (r, z) [mm]; kedjorna ska mötas i ändpunkterna.
    h: elementstorlek [mm]. h_func(r, z) -> h kan ge lokal storlek (långsamt, Python-callback).
    zones: [(z_min, z_max, r_max, h_zon[, kedjor])] lokal storlek i cylindriska zoner (gmsh Box-fält, snabbt);
        kedjor (lista av kedjeindex) begränsar zonen till de kedjornas ytor.
    Returnerar nodes (n, 3), elems {"S8R": (m, 8), "S6": (k, 6)} (0-baserade), info med
    elementens kedjeindex ("chain": {"S8R": ..., "S6": ...}) och tidsåtgång.
    """
    import gmsh
    t0 = time.perf_counter()
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("rev")
    occ = gmsh.model.occ
    pt_cache = {}

    def pt(r, z):
        key = (round(r, 6), round(z, 6))
        if key not in pt_cache:
            pt_cache[key] = occ.addPoint(r, 0, z)
        return pt_cache[key]

    inp, inp_chain = [], []
    for ci, P in enumerate(chains):
        P = np.asarray(P, float)
        ids = [occ.addPoint(r, 0, z) for r, z in P]       # fragment slår ihop kedjornas gemensamma cirklar
        c = occ.addSpline(ids)
        for k in range(2):                                    # två halvvarv: undviker periodisk söm
            rv = occ.revolve([(1, c)], 0, 0, 0, 0, 0, 1, np.pi)
            for d in rv:
                if d[0] == 2:
                    inp.append(d); inp_chain.append(ci)
            occ.rotate([(1, c)], 0, 0, 0, 0, 0, 1, np.pi)
    _, omap = occ.fragment(inp, [])                          # konformt nät över kedjegränser
    occ.synchronize()
    chain_of = {}
    for ci, outs in zip(inp_chain, omap):
        for d in outs:
            if d[0] == 2:
                chain_of[d[1]] = ci
    surfs = [(2, s) for s in chain_of]
    gmsh.option.setNumber("Mesh.MeshSizeMax", h)
    gmsh.option.setNumber("Mesh.MeshSizeMin", h / 20)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    if zones:
        fl = []
        for zz in zones:
            z0, z1, rmax, hz = zz[:4]
            f = gmsh.model.mesh.field.add("Box")
            for k, v in (("VIn", hz), ("VOut", h), ("XMin", -rmax), ("XMax", rmax), ("YMin", -rmax), ("YMax", rmax),
                         ("ZMin", z0), ("ZMax", z1)):
                gmsh.model.mesh.field.setNumber(f, k, v)
            if len(zz) > 4:
                ss = [s_ for s_, c_ in chain_of.items() if c_ in zz[4]]
                cc = sorted({abs(b[1]) for s_ in ss for b in gmsh.model.getBoundary([(2, s_)], oriented=False)})
                fr = gmsh.model.mesh.field.add("Restrict")
                gmsh.model.mesh.field.setNumber(fr, "InField", f)
                gmsh.model.mesh.field.setNumbers(fr, "SurfacesList", ss)
                gmsh.model.mesh.field.setNumbers(fr, "CurvesList", cc)
                f = fr
            fl.append(f)
        fm = gmsh.model.mesh.field.add("Min"); gmsh.model.mesh.field.setNumbers(fm, "FieldsList", fl)
        gmsh.model.mesh.field.setAsBackgroundMesh(fm)
    if h_func is not None:
        gmsh.model.mesh.setSizeCallback(lambda dim, tag, x, y, z, lc: min(h, h_func(np.hypot(x, y), z)))
    if structured:
        # strukturerat nät: samma antal delningar på alla halvcirklar, profilkurvor efter längd
        rmax = max(np.abs(np.asarray(P)[:, 0]).max() for P in chains)
        n_circ = int(np.ceil(np.pi * rmax / h))
        for _, c in gmsh.model.getEntities(1):
            bb = gmsh.model.getBoundingBox(1, c)
            if abs(bb[4] - bb[1]) < 1e-6:                     # i planet y = 0: profilkurva
                L = occ.getMass(1, c)
                gmsh.model.mesh.setTransfiniteCurve(c, int(np.ceil(L / h)) + 1)
            else:
                gmsh.model.mesh.setTransfiniteCurve(c, n_circ + 1)
        for _, s in surfs:
            gmsh.model.mesh.setTransfiniteSurface(s)
            gmsh.model.mesh.setRecombine(2, s)
    elif quad:
        gmsh.option.setNumber("Mesh.Algorithm", 6)            # Frontal-Delaunay + Blossom (alg 8 misslyckas på roterade ytor)
        gmsh.option.setNumber("Mesh.RecombineAll", 1)
        gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 1)
    gmsh.option.setNumber("Mesh.ElementOrder", 2)
    gmsh.option.setNumber("Mesh.SecondOrderIncomplete", 1)
    gmsh.model.mesh.generate(2)
    tags, xyz, _ = gmsh.model.mesh.getNodes()
    idx = {t: i for i, t in enumerate(tags)}
    nodes = xyz.reshape(-1, 3)
    elems = {"S8R": [], "S6": []}
    chain = {"S8R": [], "S6": []}
    for _, s in surfs:
        et, _, en = gmsh.model.mesh.getElements(2, s)
        for t, conn in zip(et, en):
            name = {16: "S8R", 9: "S6"}.get(int(t))
            if name is None:
                raise RuntimeError(f"Oväntad elementtyp {t}")
            nn = 8 if name == "S8R" else 6
            c = np.vectorize(idx.get)(conn).reshape(-1, nn)
            elems[name].append(c)
            chain[name].append(np.full(len(c), chain_of[s]))
    gmsh.finalize()
    for k in list(elems):
        elems[k] = np.vstack(elems[k]) if elems[k] else np.zeros((0, 8 if k == "S8R" else 6), int)
        chain[k] = np.concatenate(chain[k]) if chain[k] else np.zeros(0, int)
    used = np.unique(np.concatenate([e.ravel() for e in elems.values()]))
    remap = -np.ones(len(nodes), int); remap[used] = np.arange(len(used))
    nodes = nodes[used]
    elems = {k: remap[v] for k, v in elems.items()}
    elems, nflip = orient(elems)
    return nodes, elems, dict(chain=chain, time=time.perf_counter() - t0, flipped=nflip)


_REV = {"S8R": [0, 3, 2, 1, 7, 6, 5, 4], "S6": [0, 2, 1, 5, 4, 3]}


def orient(elems):
    """Gör elementens normaler konsistenta över hela ytan (även över kedjegränser och halvvarvsfogar).

    Grannar som delar en kant ska genomlöpa den i motsatt riktning; annars vänds grannen. Motriktade normaler vid
    tangentkontinuerliga fogar ger CalculiX-felet "gen3dnor: size of estimated shell normal ... smaller than 1.e-10"
    och felaktiga medelnormaler. Returnerar (elems, antal vända element). Kräver orienterbar yta.
    """
    from collections import deque
    keys = [k for k in ("S8R", "S6") if len(elems.get(k, []))]
    conn = [(k, i, elems[k][i, :4 if k == "S8R" else 3]) for k in keys for i in range(len(elems[k]))]
    edges = {}
    for j, (_, _, c) in enumerate(conn):
        for a in range(len(c)):
            u, v = int(c[a]), int(c[(a + 1) % len(c)])
            edges.setdefault((min(u, v), max(u, v)), []).append((j, u < v))
    flip = np.zeros(len(conn), bool); seen = np.zeros(len(conn), bool)
    nb = [[] for _ in conn]
    for lst in edges.values():
        if len(lst) == 2:
            (a, da), (b, db) = lst
            nb[a].append((b, da == db)); nb[b].append((a, da == db))
    for s0 in range(len(conn)):
        if seen[s0]:
            continue
        seen[s0] = True; q = deque([s0])
        while q:
            a = q.popleft()
            for b, same in nb[a]:
                want = flip[a] ^ same          # samma riktning på kanten => motsatt orientering
                if not seen[b]:
                    seen[b] = True; flip[b] = want; q.append(b)
                elif flip[b] != want:
                    raise RuntimeError("Ytan är inte orienterbar (eller nätet har en felaktig kant)")
    out = {k: v.copy() for k, v in elems.items()}
    for j, (k, i, _) in enumerate(conn):
        if flip[j]:
            out[k][i] = elems[k][i, _REV[k]]
    return out, int(flip.sum())


def ring_nodes(nodes, r, z, tol=1e-3):
    """Noder på cirkeln (r, z) [0-baserade index]."""
    rr = np.hypot(nodes[:, 0], nodes[:, 1])
    return np.where((np.abs(rr - r) < tol) & (np.abs(nodes[:, 2] - z) < tol))[0]


def elem_centroids(nodes, elems):
    return {k: nodes[v[:, :4 if k == "S8R" else 3]].mean(axis=1) for k, v in elems.items()}


# ---------------------------------------------------------------- inp
def _fmt_ids(ids, per=16):
    ids = [str(int(i) + 1) for i in ids]
    return "\n".join(",".join(ids[i:i + per]) for i in range(0, len(ids), per))


def write_inp(path, nodes, elems, thickness, materials, mat_of, nsets, steps, round_t=4, quad_type="S8R"):
    """Skriver CalculiX-indata.

    thickness: {"S8R": (m,), "S6": (k,)} tjocklek per element [mm]; grupperas per avrundat värde.
    materials: {namn: {"E":, "nu":, "plastic": [(sigma, eps_p), ...] eller None}}
    mat_of: {"S8R": array av materialnamn per element, ...}
    nsets: {namn: 0-baserade nodindex}
    steps: text för *STEP-blocken (se buckle_step, collapse_step).
    quad_type: "S8R" (standard) eller "S8". Använd S8R på strukturerat nät (revolve_shell(structured=True)):
        S8R på ostrukturerat, lokalt förfinat nät gav falska kollapsmoder (vågantal följde nätet, ingen fallande gren,
        lasten sjönk med förfining). S8 (full integration) låser: cylinderkollaps +14 % vid h = l/3, +4 % vid h = l/6.
        Nyckeln i elems heter "S8R" oavsett typ.
    """
    L = ["*HEADING", "femtools.shell", "*NODE"]
    L += [f"{i + 1},{x:.6f},{y:.6f},{z:.6f}" for i, (x, y, z) in enumerate(nodes)]
    eid, groups = 1, {}
    for typ, conn in elems.items():
        if len(conn) == 0:
            continue
        L.append(f"*ELEMENT,TYPE={quad_type if typ == 'S8R' else typ},ELSET=E{typ}")
        for k, c in enumerate(conn):
            L.append(f"{eid}," + ",".join(str(i + 1) for i in c))
            key = (round(float(thickness[typ][k]), round_t), mat_of[typ][k])
            groups.setdefault(key, []).append(eid)
            eid += 1
    for name, ids in nsets.items():
        L += [f"*NSET,NSET={name}", _fmt_ids(ids)]
    for j, ((t, mat), ids) in enumerate(sorted(groups.items())):
        L += [f"*ELSET,ELSET=T{j}", _fmt_ids(np.array(ids) - 1),
              f"*SHELL SECTION,ELSET=T{j},MATERIAL={mat}", f"{t}"]
    for name, m in materials.items():
        L += [f"*MATERIAL,NAME={name}", "*ELASTIC", f"{m['E']},{m['nu']}"]
        if m.get("plastic"):
            L.append("*PLASTIC")
            L += [f"{s},{e}" for s, e in m["plastic"]]
    L.append(steps)
    open(path, "w").write("\n".join(L) + "\n")
    return dict(n_nodes=len(nodes), n_elems=eid - 1, n_sections=len(groups))


def buckle_step(fix_bcs, cloads, nmodes=5):
    """Linjär knäckning. fix_bcs: [(nset, dof_från, dof_till)], cloads: [(nodindex 0-bas, dof, kraft)]."""
    L = ["*BOUNDARY"] + [f"{n},{a},{b}" for n, a, b in fix_bcs]
    L += ["*STEP", "*BUCKLE", f"{nmodes},1e-7,{max(4*nmodes,40)},2000", "*CLOAD"]
    L += [f"{i + 1},{d},{f:.8g}" for i, d, f in cloads]
    L += ["*NODE FILE,OUTPUT=2D", "U", "*END STEP"]
    return "\n".join(L)


def collapse_step(fix_bcs, disp_set, delta, inc0=0.02, inc_min=1e-5, inc_max=0.05, inc_n=400,
                  out_set=None, frequency=1, stress=True):
    """Olinjär (NLGEOM) förskjutningsstyrd belastning: disp_set får uz = -delta linjärt i tiden 0..1."""
    L = ["*BOUNDARY"] + [f"{n},{a},{b}" for n, a, b in fix_bcs]
    L += [f"*STEP,NLGEOM,INC={inc_n}", "*STATIC", f"{inc0},1.0,{inc_min},{inc_max}",
          "*BOUNDARY", f"{disp_set},3,3,{-delta}",
          f"*NODE PRINT,NSET={out_set or disp_set},TOTALS=ONLY", "RF",
          f"*NODE FILE,OUTPUT=2D,FREQUENCY={frequency}", "U"]
    if stress:                                    # mittytans värden (2D-utdata): membranspänning och PEEQ
        L += [f"*EL FILE,OUTPUT=2D,FREQUENCY={frequency}", "S,PEEQ"]
    L += ["*END STEP"]
    return "\n".join(L)


# ---------------------------------------------------------------- körning och läsning
def ccx_exe():
    """ccx_pardiso (egen build med MKL Pardiso, 2–3 gånger snabbare) om den finns, annars ccx (apt, SPOOLES).
    Styrs av miljövariabeln FEMTOOLS_CCX."""
    import shutil
    return os.environ.get("FEMTOOLS_CCX") or ("ccx_pardiso" if shutil.which("ccx_pardiso") else "ccx")


def run(job, cwd, ncpu=4, timeout=None, log=None, exe=None):
    exe = exe or ccx_exe()
    env = dict(os.environ, OMP_NUM_THREADS=str(ncpu), MKL_NUM_THREADS=str(ncpu), CCX_NPROC_EQUATION_SOLVER=str(ncpu))
    t0 = time.perf_counter()
    with open(log or os.path.join(cwd, job + ".log"), "w") as f:
        p = subprocess.run([exe, "-i", job], cwd=cwd, stdout=f, stderr=subprocess.STDOUT, env=env, timeout=timeout)
    txt = open(log or os.path.join(cwd, job + ".log")).read()
    ok = p.returncode == 0 and "*ERROR" not in txt
    return dict(ok=ok, exe=exe, time=time.perf_counter() - t0, returncode=p.returncode,
                errors=[l for l in txt.splitlines() if "ERROR" in l][:5])


def read_buckle(dat):
    """Knäckfaktorer ur .dat."""
    txt = open(dat).read()
    m = re.search(r"B U C K L I N G   F A C T O R   O U T P U T(.*?)(\n\s*\n\s*\n|\Z)", txt, re.S)
    vals = []
    for line in m.group(1).splitlines():
        s = line.split()
        if len(s) == 2 and s[0].isdigit():
            vals.append(float(s[1]))
    return np.array(vals)


def read_totals(dat):
    """Summerade reaktionskrafter per inkrement ur .dat (NODE PRINT TOTALS=ONLY): (tid, (n,3))."""
    t, F = [], []
    cur = None
    for line in open(dat):
        m = re.search(r"(?:total )?forces? \(fx,fy,fz\) for set \S+ and time\s+([-\d.E+]+)", line)
        if m:
            cur = float(m.group(1)); continue
        if cur is not None and line.strip():
            s = line.split()
            try:
                F.append([float(v) for v in s[-3:]]); t.append(cur)
            except ValueError:
                pass
            cur = None
    return np.array(t), np.array(F)


def read_frd(frd, nnodes, fields=False):
    """Nodresultat ur .frd för mittytans noder (1..nnodes).

    fields=False: lista av (steg-värde, U (nnodes, 3)) (bakåtkompatibelt).
    fields=True: lista av dict(value=, DISP=(n,3), STRESS=(n,6) xx,yy,zz,xy,yz,zx, PE=(n,1) ...) per inkrement.
    """
    out, cur, name, arr, ncomp = [], None, None, None, 0
    with open(frd) as f:
        for line in f:
            if line.startswith("  100CL"):
                val = float(line[12:24])
                if cur is None or cur["value"] != val:
                    if cur is not None and cur.get("_done"):
                        out.append(cur)
                    cur = dict(value=val)
                name = None
            elif line.startswith(" -4"):
                s = line.split(); name, ncomp = s[1], int(s[2])
                if name == "DISP":
                    ncomp = 3
                arr = np.zeros((nnodes, ncomp)); cur[name] = arr; cur["_done"] = True
            elif line.startswith(" -5"):
                continue
            elif name is not None and line.startswith(" -1"):
                nid = int(line[3:13])
                if nid <= nnodes:
                    ln = line.rstrip("\n")
                    vals = [float(ln[13 + 12 * k: 25 + 12 * k]) for k in range(min(ncomp, 6, (len(ln) - 13) // 12))]
                    arr[nid - 1, :len(vals)] = vals
            elif name is not None and line.startswith(" -3"):
                name = None
    if cur is not None and cur.get("_done"):
        out.append(cur)
    for c in out:
        c.pop("_done", None)
    if fields:
        return out
    return [(c["value"], c["DISP"]) for c in out if "DISP" in c]


def edge_loads(nodes, elems, ring, total, dof=3):
    """Konsistenta nodkrafter för jämn linjelast längs en ring av kvadratiska elementkanter.

    ring: nodindex på kanten. total: total kraft [N] (tecken = riktning i dof). Hörn 1/6, mitt 4/6 per kant.
    """
    S = set(int(i) for i in ring)
    F = {}
    pairs = {"S8R": [(0, 1, 4), (1, 2, 5), (2, 3, 6), (3, 0, 7)], "S6": [(0, 1, 3), (1, 2, 4), (2, 0, 5)]}
    Ltot = 0.0
    for typ, conn in elems.items():
        for c in conn:
            for a, b, m in pairs[typ]:
                if c[a] in S and c[b] in S and c[m] in S:
                    L = np.linalg.norm(nodes[c[a]] - nodes[c[m]]) + np.linalg.norm(nodes[c[m]] - nodes[c[b]])
                    Ltot += L
                    for n, w in ((c[a], 1 / 6), (c[b], 1 / 6), (c[m], 4 / 6)):
                        F[n] = F.get(n, 0.0) + w * L
    return [(n, dof, total * f / Ltot) for n, f in F.items()], Ltot


def available():
    """CalculiX-version eller None."""
    import shutil
    exe = ccx_exe()
    if not shutil.which(exe):
        return None
    p = subprocess.run([exe, "-v"], capture_output=True, text=True)
    m = re.search(r"Version\s+([\d.]+)", p.stdout + p.stderr)
    return (m.group(1) if m else "okänd") + (" (Pardiso)" if "pardiso" in exe else " (SPOOLES)")


def koiter_axisym(xi, nu):
    """Koiters (1963) knäckfaktor för axisymmetrisk imperfektion w0 = xi·t (klassisk våglängd), lång cylinder."""
    c = 1.5 * np.sqrt(3 * (1 - nu**2)) * abs(xi)
    b = 2 + c
    return (b - np.sqrt(b**2 - 4)) / 2


def cylinder_collapse(workdir, R=29.0, t=0.1, L=20.0, E=69000.0, nu=0.33, h=1.0, xi=0.1, seed=0.01, n_seed=12, ncpu=4):
    """Verifieringsfall: axiellt tryckt cylinder, ledad med radiellt låsta ränder, olinjär förskjutningsstyrning.

    Imperfektion w = t·[xi·cos(π(z−L/2)/l) + seed·cos(n_seed·θ)]·sin(πz/L), l = klassisk axisymmetrisk halvvåglängd.
    Returnerar dict med F_max, P_cl, kvot, Koiters värde, axiell styvhet mot EA/L och tid.
    """
    t0 = time.perf_counter()
    os.makedirs(workdir, exist_ok=True)
    P = np.c_[np.full(41, R), np.linspace(0, L, 41)]
    nodes, el, _ = revolve_shell([P], h, structured=True)
    bot, top = ring_nodes(nodes, R, 0), ring_nodes(nodes, R, L)
    lh = np.pi * np.sqrt(R * t) / (12 * (1 - nu**2)) ** 0.25
    z, th = nodes[:, 2], np.arctan2(nodes[:, 1], nodes[:, 0])
    w = t * (xi * np.cos(np.pi * (z - L / 2) / lh) + seed * np.cos(n_seed * th)) * np.sin(np.pi * z / L)
    X = nodes.copy(); X[:, 0] += w * np.cos(th); X[:, 1] += w * np.sin(th)
    P_cl = 2 * np.pi * E * t**2 / np.sqrt(3 * (1 - nu**2)); k = E * 2 * np.pi * R * t / L
    d = 1.25 * P_cl / k
    st = collapse_step([("BOT", 1, 3), ("TOP", 1, 2)], "TOP", d, inc0=0.1, inc_max=0.1, inc_min=1e-5, inc_n=300)
    tk = {key: np.full(len(v), t) for key, v in el.items()}
    mk = {key: np.array(["AL"] * len(v)) for key, v in el.items()}
    write_inp(os.path.join(workdir, "cyl.inp"), X, el, tk, {"AL": dict(E=E, nu=nu)}, mk, {"BOT": bot, "TOP": top}, st)
    r = run("cyl", workdir, ncpu)
    tt, F = read_totals(os.path.join(workdir, "cyl.dat"))
    Fz = -F[:, 2]
    return dict(ok=r["ok"], F_max=float(Fz.max()), P_cl=P_cl, ratio=float(Fz.max() / P_cl), koiter=float(koiter_axisym(xi, nu)),
                k=float(Fz[0] / (tt[0] * d)), k_hand=k, peak_last=bool(np.argmax(Fz) == len(Fz) - 1),
                n_nodes=len(nodes), time=time.perf_counter() - t0)
