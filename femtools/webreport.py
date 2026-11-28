"""Interaktiv FEM-rapport som fristående HTML (roterbar 3D-modell, lastfallsväljare, kopierbara tabeller).

Komplement till PDF-rapporten (femtools.report + fem-rapport-skillen): samma nodala fält, men läsaren kan
rotera, byta lastfall och fält, skala deformationen och kopiera tabellerna (tabbseparerat, klistras in
direkt i PowerPoint, Excel och Word). Sidan publiceras som artifact eller öppnas lokalt i webbläsaren.

    from femtools import webreport, report
    n = report.nodal(r, E, nu)
    w = webreport.Report("Tryckkärl P1", meta={"Fall": "261126_x", "Standard": "4.2"})
    w.mesh(m)
    w.field("LF1 Provtryck", "von Mises", n["vm"], "MPa", u=n["u"], limit=Rp02, limit_label="Rp0,2")
    w.field("LF1 Provtryck", "Förskjutning", n["umag"], "mm", u=n["u"])
    w.field("LF2 Drift", "Temperatur", webreport.vertex_values(rh["basis"], rh["T"]), "°C")
    w.table("Resultat", ["Storhet", "FEM", "Hand", "Avvikelse"], rows, caption="...")
    w.text("Slutsats", "Utnyttjandegrad 0,62 ...")
    w.write("fig/rapport.html")

Betraktarvy (standard 5.0): ange produktens framsida med front, riktningen från produkten mot betraktaren,
och up. Sidan öppnar i "Fram iso" så att modellen syns som en människa ser den framför produkten, inte från
väggen eller underlaget. Exempel, produkt skruvad mot en vägg i y = −120: front=(0, 1, 0), up=(0, 0, 1).

CFD i samma visare: f = openfoam.read_case(fall, rho=998); w.flow("D2 spolning", f, slices={"y": 0},
seeds=...) lägger vätskedomänen som eget fall med |U| i snitt, väggtryck och strömlinjer. Strukturen visas
då som genomskinlig kontur. Snittplan (x, y, z) finns för alla fall.

Bara ytan skickas till sidan (randtrianglar och deras noder), så sidan blir liten: cirka 60 byte per
ytnod och fält. Över cirka 200 000 ytnoder: grovare nät för visning, eller skriv data som separat fil.
"""
import json
import datetime
import numpy as np

_CDN = ("https://cdn.jsdelivr.net/npm/three@0.128.0/build/three.min.js",
        "https://cdn.jsdelivr.net/npm/three@0.128.0/examples/js/controls/OrbitControls.js")


def vertex_values(basis, x):
    """Värden i hörnnoderna för en skalär P1/P2-lösning (t.ex. temperatur från thermal.heat)."""
    return np.asarray(x)[basis.nodal_dofs[0]]


def slice_tets(P, tet, val, axis, value):
    """Plant snitt genom ett tetraedernät (marching tetrahedra). P (n, 3), tet (nc, 4), val (n,) nodvärden.
    axis 0/1/2 eller "x"/"y"/"z". Returnerar (punkter (m, 3), trianglar (k, 3), värden (m,)), linjärt interpolerat."""
    ax = "xyz".index(axis) if isinstance(axis, str) else int(axis)
    d = P[:, ax] - value
    pos = d[tet] > 0
    nn = pos.sum(1)
    cut = np.where((nn > 0) & (nn < 4))[0]
    pts, vals, tris = [], [], []
    E = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
    o = [i for i in range(3) if i != ax]
    for c in cut:
        v = tet[c]
        q = []
        for a, b in E:
            ia, ib = v[a], v[b]
            if (d[ia] > 0) != (d[ib] > 0):
                w = d[ia] / (d[ia] - d[ib])
                q.append((P[ia] + w * (P[ib] - P[ia]), val[ia] + w * (val[ib] - val[ia])))
        if len(q) == 4:                       # fyrhörning: ordna runt mittpunkten
            X = np.array([x for x, _ in q]); cen = X.mean(0)
            q = [q[i] for i in np.argsort(np.arctan2(X[:, o[1]] - cen[o[1]], X[:, o[0]] - cen[o[0]]))]
        i0 = len(pts)
        tris += [[i0, i0 + 1, i0 + 2]] + ([[i0, i0 + 2, i0 + 3]] if len(q) == 4 else [])
        for x, u in q:
            pts.append(x); vals.append(u)
    if not pts:
        return np.zeros((0, 3)), np.zeros((0, 3), int), np.zeros(0)
    pts = np.array(pts)
    uq, first, inv = np.unique(np.round(pts, 6), axis=0, return_index=True, return_inverse=True)
    return pts[first], inv.ravel()[np.array(tris)], np.asarray(vals)[first]


def streamlines(P, tet, U, seeds, ds=None, max_steps=4000, nudges=300):
    """Strömlinjer i ett nodalt hastighetsfält på tetraedrar (RK2, fast steglängd ds i längdenheten för P).
    Vid väggen (punkten hamnar utanför nätet) knuffas linjen mot senaste cellens mittpunkt, högst `nudges` gånger
    totalt och 8 i följd (då har linjen nått utloppet eller fastnat i en stagnationspunkt).
    Returnerar lista med (punkter (k, 3), fart (k,)). Visualisering: ingen massbalans per linje."""
    from scipy.spatial import cKDTree
    T = P[tet]
    cen = T.mean(1)
    tree = cKDTree(cen)
    Minv = np.linalg.inv(np.transpose(T[:, 1:] - T[:, :1], (0, 2, 1)))
    Um = np.linalg.norm(U, axis=1)
    umin = 1e-4 * Um.max()
    if ds is None:
        ds = 0.5 * np.cbrt(np.abs(np.linalg.det(T[:, 1:] - T[:, :1])).mean() / 6)

    def vel(x, guess):
        for k in (16, 128):
            ids = tree.query(x, k=min(k, len(cen)))[1]
            if guess is not None:
                ids = np.concatenate([[guess], ids])
            for c in ids:
                lam = Minv[c] @ (x - T[c, 0])
                bc = np.array([1 - lam.sum(), *lam])
                if bc.min() > -1e-9:
                    return bc @ U[tet[c]], c
        return None, None

    out = []
    for x in np.asarray(seeds, float):
        u0, g = vel(x, None)
        if u0 is None:
            continue
        X, S, nk, nrow, ok_steps = [x.copy()], [np.linalg.norm(u0)], 0, 0, 0
        for _ in range(max_steps):
            u1, g1 = vel(x, g)
            ok = u1 is not None and np.linalg.norm(u1) > umin
            if ok:
                u2, g2 = vel(x + 0.5 * ds * u1 / np.linalg.norm(u1), g1)
                ok = u2 is not None and np.linalg.norm(u2) > umin
            if ok:
                xn = x + ds * u2 / np.linalg.norm(u2)
                ok = vel(xn, g2)[0] is not None
            if not ok:
                if nk >= nudges or nrow >= 8 or g is None or u1 is None:     # utlopp eller fast i stagnation
                    break
                nk += 1; nrow += 1; ok_steps = 0
                x = x + ds * (cen[g] - x) / max(np.linalg.norm(cen[g] - x), 1e-12)
                continue
            x, g, ok_steps = xn, g2, ok_steps + 1
            if ok_steps >= 3:
                nrow = 0
            X.append(x.copy()); S.append(np.linalg.norm(u2))
        if len(X) > 3:
            out.append((np.array(X), np.array(S)))
    return out


class Report:
    def __init__(self, title, meta=None, subtitle=None, front=(0, -1, 0), up=(0, 0, 1)):
        """front: riktning från produkten mot betraktaren (produktens framsida), up: uppåt. Startvyn är
        "Fram iso" från front. Ange alltid front för produkter med en given monteringsriktning."""
        self.title, self.subtitle = title, subtitle
        f, u = np.asarray(front, float), np.asarray(up, float)
        if abs(np.dot(f, u)) > 1e-6 * np.linalg.norm(f) * np.linalg.norm(u):
            raise ValueError("front och up måste vara vinkelräta")
        self.front, self.up = (f / np.linalg.norm(f)).tolist(), (u / np.linalg.norm(u)).tolist()
        self.flows = []
        self.meta = dict(meta or {})
        self.meta.setdefault("Skapad", datetime.date.today().isoformat())
        try:
            from . import __version__
            self.meta.setdefault("femtools", __version__)
        except ImportError:  # pragma: no cover
            pass
        self.cases, self.blocks = {}, []
        self._surf = None

    # ------------------------------------------------------------------ modell
    def mesh(self, m):
        """Ytan av ett tetraedernät (skfem MeshTet). Anropas en gång, före field()."""
        bf = m.boundary_facets()
        tri = m.facets[:, bf].copy()                                 # (3, ntri), globala nodnummer
        # orientera utåt (normalen bort från elementets tyngdpunkt) så att skuggningen blir rätt
        P0, P1, P2 = (m.p[:, tri[i]] for i in range(3))
        nrm = np.cross((P1 - P0).T, (P2 - P0).T).T
        cen = m.p[:, m.t[:, m.f2t[0, bf]]].mean(axis=1)
        flip = np.einsum("ij,ij->j", nrm, P0 - cen) < 0
        tri[1, flip], tri[2, flip] = tri[2, flip], tri[1, flip].copy()
        used = np.unique(tri)
        remap = -np.ones(m.p.shape[1], dtype=np.int64)
        remap[used] = np.arange(len(used))
        self._surf = dict(nodes=used, p=m.p[:, used].T, tri=remap[tri].T)
        return self

    def surface(self, p, tri, quads=None):
        """Skalyta (t.ex. mittytan från femtools.shell): p (N, 3) noder, tri (k, 3) trianglar (0-baserade).
        Ritas dubbelsidigt. quads (m, 4) ges för att nätlinjerna ska visa fyrhörningarna i stället för trianglar.
        Anropas en gång, före field(); nodala fält ges per nod i p."""
        p = np.asarray(p, float); tri = np.asarray(tri, dtype=np.int64)
        self._surf = dict(nodes=np.arange(len(p)), p=p, tri=tri, double=True,
                          quads=None if quads is None else np.asarray(quads, dtype=np.int64))
        return self

    def field(self, case, name, values, unit, u=None, limit=None, limit_label=None, digits=4):
        """Nodalt fält (värde per hörnnod i hela nätet) i lastfallet `case`. u: (3, N) förskjutning för
        deformerad form (samma fall delar u). limit: gränsvärde som markeras i färgskalan."""
        if self._surf is None:
            raise RuntimeError("Anropa mesh(m) före field().")
        idx = self._surf["nodes"]
        v = np.asarray(values, float)[idx]
        c = self.cases.setdefault(case, dict(name=case, u=None, fields=[]))
        if u is not None and c["u"] is None:
            c["u"] = np.round(np.asarray(u, float)[:, idx].T.ravel(), 7).tolist()
        imax = int(np.argmax(v))
        c["fields"].append(dict(name=name, unit=unit, values=_sig(v, digits), limit=limit,
                                limit_label=limit_label, min=float(v.min()), max=float(v.max()),
                                at=self._surf["p"][imax].round(3).tolist()))
        return self

    def transient(self, case, name, t, frames, unit, curves=None, t_unit="min", t_div=60.0, digits=3,
                  limit=None, limit_label=None, curve_unit=None, curve_title=None):
        """Transient nodalt fält med tidsreglage: frames är en lista med nodvärden (en per sparad tid t [s],
        samma längd som nätets noder). Färgskalan är fast över alla tider. curves {namn: (t [s], y)} ritas
        som kurvdiagram under 3D-vyn med en markör vid vald tid (t.ex. svalningskurva). t visas som t/t_div i t_unit."""
        if self._surf is None:
            raise RuntimeError("Anropa mesh(m) före transient().")
        idx = self._surf["nodes"]
        F = np.asarray([np.asarray(v, float)[idx] for v in frames])
        c = self.cases.setdefault(case, dict(name=case, u=None, fields=[]))
        imax = np.unravel_index(int(np.argmax(F)), F.shape)
        c["fields"].append(dict(name=name, unit=unit, values=_sig(F[0], digits), limit=limit, limit_label=limit_label,
                                min=float(F.min()), max=float(F.max()), at=self._surf["p"][imax[1]].round(3).tolist(),
                                frames=[_sig(f, digits) for f in F], t=(np.asarray(t, float) / t_div).round(4).tolist(),
                                t_unit=t_unit, curve_unit=curve_unit or unit, curve_title=curve_title,
                                curves=[dict(name=k, t=(np.asarray(v[0], float) / t_div).round(4).tolist(),
                                             y=_sig(v[1], 4)) for k, v in (curves or {}).items()]))
        return self

    def flow(self, case, f, slices=None, seeds=None, ds=None, wall="walls", unitU="m/s", unitp="Pa", digits=4):
        """CFD-fält som eget fall: f från openfoam.read_case (P, tet, U, p, patches) eller en dict med samma nycklar.
        slices {"y": 0.0, ...}: plana snitt med |U|. seeds (k, 3): startpunkter för strömlinjer (t.ex. på inloppet).
        Väggpatchen ritas med nodalt tryck p."""
        P, tet, U = np.asarray(f["P"], float), np.asarray(f["tet"]), np.asarray(f["U"], float)
        Um = np.linalg.norm(U, axis=1)
        fl = dict(name=case, slices=[], wall=None, stream=[], unitU=unitU, unitp=unitp,
                  Umax=float(Um.max()), time=f.get("time"))
        for ax, v in (slices or {}).items():
            pts, tri, val = slice_tets(P, tet, Um, ax, v)
            fl["slices"].append(dict(name=f"{ax} = {v:g}", p=np.round(pts, 3).ravel().tolist(),
                                     tri=tri.ravel().astype(int).tolist(), values=_sig(val, digits)))
        wt = f.get("patches", {}).get(wall)
        if wt is not None and f.get("p") is not None:
            used = np.unique(wt); rm = -np.ones(len(P), int); rm[used] = np.arange(len(used))
            pv = np.asarray(f["p"], float)[used]
            fl["wall"] = dict(p=np.round(P[used], 3).ravel().tolist(), tri=rm[wt].ravel().astype(int).tolist(),
                              values=_sig(pv, digits), min=float(pv.min()), max=float(pv.max()))
        if seeds is not None:
            for X, S in streamlines(P, tet, U, seeds, ds=ds):
                X, S = X[::2], S[::2]
                fl["stream"].append(dict(p=np.round(X, 2).ravel().tolist(), values=_sig(S, 3)))
        self.flows.append(fl)
        return self

    # ------------------------------------------------------------------ text och tabeller
    def table(self, title, columns, rows, caption=None):
        """Tabell; rader som listor (tal formateras med svensk decimalkomma vid visning och kopiering)."""
        self.blocks.append(dict(kind="table", title=title, columns=list(columns),
                                rows=[[_cell(x) for x in r] for r in rows], caption=caption))
        return self

    def text(self, heading, body):
        """Textavsnitt. Stycken separeras med tom rad."""
        self.blocks.append(dict(kind="text", title=heading, body=body))
        return self

    # ------------------------------------------------------------------ utdata
    def data(self):
        s = self._surf
        return dict(title=self.title, subtitle=self.subtitle, meta=self.meta, front=self.front, up=self.up, flows=self.flows,
                    mesh=None if s is None else dict(p=np.round(s["p"], 4).ravel().tolist(),
                                                     tri=s["tri"].ravel().astype(int).tolist(),
                                                     double=bool(s.get("double", False)),
                                                     quads=None if s.get("quads") is None else s["quads"].ravel().astype(int).tolist()),
                    cases=list(self.cases.values()), blocks=self.blocks)

    def html(self):
        d = json.dumps(self.data(), ensure_ascii=False, separators=(",", ":")).replace("</", "<\\/")
        return (_TEMPLATE.replace("__TITLE__", _esc(self.title))
                .replace("__THREE__", _CDN[0]).replace("__ORBIT__", _CDN[1]).replace("__DATA__", d))

    def write(self, path):
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(self.html())
        return path


def _sig(v, digits):
    a = np.asarray(v, float)
    scale = np.abs(a).max() or 1.0
    q = 10.0 ** (np.floor(np.log10(scale)) - digits)
    return (np.round(a / q) * q).round(12).tolist()


def _cell(x):
    if isinstance(x, (np.floating, np.integer)):
        x = x.item()
    return x


def _esc(s):
    return str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


_TEMPLATE = r"""<title>__TITLE__</title>
<meta name="claude:viewport-fit-sides" content="cover">
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500&family=IBM+Plex+Sans+Condensed:wght@500;600&family=IBM+Plex+Sans:wght@400;500;600&display=swap">
<style>
/* Layout: rapportblad — huvud med fallets metadata, visare (3D + reglage) i full bredd, sedan tabeller och text i läsbredd */
:root{
  --bg:#f5f6f7; --paper:#ffffff; --ink:#16202a; --muted:#5b6774; --rule:#d7dce1;
  --accent:#1f5f8b; --accent-ink:#ffffff; --canvas:#eef1f4; --mesh:#20303f; --warn:#b4491b;
  --display:"IBM Plex Sans Condensed","Arial Narrow",system-ui,sans-serif;
  --body:"IBM Plex Sans",system-ui,-apple-system,"Segoe UI",sans-serif;
  --mono:"IBM Plex Mono",ui-monospace,"SFMono-Regular",Menlo,Consolas,monospace;
}
@media (prefers-color-scheme: dark){:root:not([data-theme="light"]){
  --bg:#0f151b; --paper:#161e26; --ink:#e3e8ed; --muted:#93a1ae; --rule:#2a3642;
  --accent:#6fb1df; --accent-ink:#0b1218; --canvas:#1b252f; --mesh:#c9d4de; --warn:#f08a5d; color-scheme:dark}}
:root[data-theme="dark"]{
  --bg:#0f151b; --paper:#161e26; --ink:#e3e8ed; --muted:#93a1ae; --rule:#2a3642;
  --accent:#6fb1df; --accent-ink:#0b1218; --canvas:#1b252f; --mesh:#c9d4de; --warn:#f08a5d; color-scheme:dark}
*{box-sizing:border-box}
body{background:var(--bg);color:var(--ink);font:15px/1.55 var(--body);
  padding-left:max(16px,env(safe-area-inset-left,0px));padding-right:max(16px,env(safe-area-inset-right,0px));padding-block:24px 48px}
.wrap{max-width:1180px;margin:0 auto;display:grid;gap:28px}
header{display:grid;gap:10px;border-bottom:2px solid var(--ink);padding-bottom:14px}
.kicker{font:500 12px/1 var(--mono);letter-spacing:.08em;text-transform:uppercase;color:var(--muted)}
h1{font:600 clamp(26px,4vw,38px)/1.1 var(--display);margin:0;text-wrap:balance}
.sub{color:var(--muted);margin:0;max-width:65ch}
dl.meta{display:flex;flex-wrap:wrap;gap:6px 22px;margin:4px 0 0;font:13px/1.4 var(--mono)}
dl.meta div{display:flex;gap:6px} dl.meta dt{color:var(--muted)} dl.meta dd{margin:0}
.viewer{display:grid;grid-template-columns:minmax(0,1fr) 260px;gap:0;background:var(--paper);border:1px solid var(--rule)}
@media (max-width:820px){.viewer{grid-template-columns:minmax(0,1fr)}}
.stage{position:relative;min-width:0;aspect-ratio:16/10;max-width:100%;background:var(--canvas)}
.stage canvas{display:block;width:100%;height:100%;touch-action:none}
.legend{position:absolute;right:12px;top:34px;bottom:12px;display:flex;gap:8px;align-items:stretch;pointer-events:none;
  font:12px/1 var(--mono);color:var(--ink)}
.legend .bar{width:14px;border:1px solid var(--rule)}
.legend .ticks{position:relative;width:64px}
.legend .ticks span{position:absolute;left:0;transform:translateY(50%);white-space:nowrap;font-variant-numeric:tabular-nums;
  text-shadow:0 0 3px var(--canvas),0 0 3px var(--canvas)}
.legend .ticks span.lim{color:var(--warn);font-weight:500}
.legend .unit{position:absolute;top:-14px;left:-22px;transform:translateY(-100%);color:var(--muted);white-space:nowrap}
.readout{position:absolute;left:12px;bottom:12px;font:12px/1.4 var(--mono);color:var(--muted);pointer-events:none;max-width:60%}
.readout b{color:var(--ink);font-weight:500}
.controls{border-left:1px solid var(--rule);padding:16px;display:grid;gap:14px;align-content:start;min-width:0}
@media (max-width:820px){.controls{border-left:0;border-top:1px solid var(--rule);grid-template-columns:repeat(auto-fit,minmax(200px,1fr))}}
.controls label{display:grid;gap:5px;font:500 12px/1.2 var(--body);letter-spacing:.04em;text-transform:uppercase;color:var(--muted)}
select,input[type=range]{width:100%;font:14px var(--body);color:var(--ink)}
select{padding:7px 8px;background:var(--bg);border:1px solid var(--rule);border-radius:3px}
.row{display:flex;gap:6px;flex-wrap:wrap}
button{font:500 13px/1 var(--body);padding:8px 10px;border:1px solid var(--rule);background:var(--bg);color:var(--ink);border-radius:3px;cursor:pointer}
button:hover{border-color:var(--accent)} button:focus-visible,select:focus-visible,input:focus-visible{outline:2px solid var(--accent);outline-offset:2px}
button.primary{background:var(--accent);color:var(--accent-ink);border-color:var(--accent)}
.check{display:flex!important;grid-auto-flow:column;align-items:center;gap:8px;text-transform:none!important;letter-spacing:0!important;font-size:14px!important;color:var(--ink)!important}
.scaleval{font:13px var(--mono);color:var(--ink);text-transform:none;letter-spacing:0}
section{display:grid;gap:10px;max-width:900px}
h2{font:600 20px/1.2 var(--display);margin:0;letter-spacing:.01em}
section p{margin:0;max-width:68ch}
.tablebox{overflow-x:auto;border:1px solid var(--rule);background:var(--paper)}
table{border-collapse:collapse;width:100%;font-size:14px}
th,td{padding:7px 12px;text-align:left;border-bottom:1px solid var(--rule);white-space:nowrap}
th{font:600 13px/1.2 var(--body);color:var(--muted);background:var(--bg)}
td.num{text-align:right;font-family:var(--mono);font-variant-numeric:tabular-nums}
tr:last-child td{border-bottom:0}
.caption{display:flex;justify-content:space-between;gap:12px;align-items:baseline;flex-wrap:wrap;color:var(--muted);font-size:13px}
.toast{font:12px var(--mono);color:var(--accent)}
.nodata{padding:24px;color:var(--muted)}
@media print{.controls,.caption button{display:none}.viewer{grid-template-columns:1fr}body{background:#fff}}
@media (prefers-reduced-motion:reduce){*{scroll-behavior:auto}}
[hidden]{display:none!important}
.tl{grid-column:1/-1;border-top:1px solid var(--rule);padding:12px 16px;display:grid;gap:8px}
.tl .row{align-items:center} .tl input[type=range]{flex:1;min-width:160px}
.tl svg{width:100%;height:auto;display:block;font:12px var(--mono)}
.tl .ax{stroke:var(--rule)} .tl .lab{fill:var(--muted)} .tl .cur{stroke:var(--warn);stroke-width:1.5}
</style>
<script src="__THREE__"></script>
<script src="__ORBIT__"></script>
<div class="wrap">
  <header>
    <div class="kicker">FEM-rapport · interaktiv</div>
    <h1 id="title"></h1>
    <p class="sub" id="subtitle" hidden></p>
    <dl class="meta" id="meta"></dl>
  </header>
  <div class="viewer" id="viewer">
    <div class="stage" id="stage">
      <canvas id="cv" aria-label="3D-modell med resultatfält"></canvas>
      <div class="legend" id="legend"><div class="bar" id="bar"></div><div class="ticks" id="ticks"><span class="unit" id="unit"></span></div></div>
      <div class="readout" id="readout"></div>
    </div>
    <div class="controls">
      <label for="case">Lastfall<select id="case"></select></label>
      <label for="fld">Resultat<select id="fld"></select></label>
      <label for="scale">Deformationsskala <span class="scaleval" id="scalev"></span><input type="range" id="scale" min="0" max="100" value="0"></label>
      <label class="check" for="edges"><input type="checkbox" id="edges"> Visa nät</label>
      <label class="check" for="maxm"><input type="checkbox" id="maxm" checked> Markera max</label>
      <label class="check" for="strm" id="strmrow" hidden><input type="checkbox" id="strm" checked> Strömlinjer</label>
      <label for="clipax">Snittplan<select id="clipax"><option value="">Av</option><option value="0">x</option><option value="1">y</option><option value="2">z</option></select></label>
      <label for="clip" id="cliprow" hidden>Läge <span class="scaleval" id="clipv"></span><input type="range" id="clip" min="0" max="1000" value="500"></label>
      <label class="check" for="clipf" id="clipfrow" hidden><input type="checkbox" id="clipf"> Vänd snittsida</label>
      <div class="row" role="group" aria-label="Vy">
        <button type="button" data-view="iso" class="primary">Fram iso</button><button type="button" data-view="front">Fram</button>
        <button type="button" data-view="top">Topp</button><button type="button" data-view="side">Sida</button><button type="button" data-view="back">Bak</button>
      </div>
    </div>
    <div class="tl" id="tl" hidden>
      <div class="row"><button type="button" id="play">▶ Spela</button>
        <input type="range" id="tsl" min="0" max="0" value="0" aria-label="Tid"><span class="scaleval" id="tval"></span></div>
      <svg id="chart" viewBox="0 0 900 250" role="img" aria-label="Kurvor över tiden"></svg>
    </div>
  </div>
  <div id="blocks" class="wrap" style="gap:28px;margin:0;justify-items:start"></div>
</div>
<script>
const D = __DATA__;
const $ = id => document.getElementById(id);
const fmt = (x, d=3) => {
  if (typeof x !== "number" || !isFinite(x)) return String(x ?? "");
  const a = Math.abs(x);
  let s;
  if (a !== 0 && (a < 1e-3 || a >= 1e7)) s = x.toExponential(d - 1);
  else if (a >= 10 ** d) s = String(Math.round(x));
  else s = String(Number(x.toPrecision(d)));
  return s.replace(".", ",");
};
const css = n => getComputedStyle(document.documentElement).getPropertyValue(n).trim();

// Turbo-färgskala (Google AI, polynomapproximation), 0..1 -> rgb
function turbo(t){ t=Math.min(1,Math.max(0,t));
  const r=0.13572138+t*(4.61539260+t*(-42.66032258+t*(132.13108234+t*(-152.94239396+t*59.28637943))));
  const g=0.09140261+t*(2.19418839+t*(4.84296658+t*(-14.18503333+t*(4.27729857+t*2.82956604))));
  const b=0.10667330+t*(12.64194608+t*(-60.58204836+t*(110.36276771+t*(-89.90310912+t*27.34824973))));
  return [Math.min(1,Math.max(0,r)),Math.min(1,Math.max(0,g)),Math.min(1,Math.max(0,b))]; }

// ---------- rubrik, metadata, tabeller, text
$("title").textContent = D.title;
if (D.subtitle){ $("subtitle").textContent = D.subtitle; $("subtitle").hidden = false; }
for (const [k,v] of Object.entries(D.meta||{})){ const d=document.createElement("div");
  d.innerHTML = "<dt></dt><dd></dd>"; d.children[0].textContent=k; d.children[1].textContent=v; $("meta").appendChild(d); }
D.blocks.forEach((b,i) => {
  const s=document.createElement("section"); const h=document.createElement("h2"); h.textContent=b.title; s.appendChild(h);
  if (b.kind==="text"){ String(b.body).split(/\n\s*\n/).forEach(p=>{const e=document.createElement("p"); e.textContent=p.trim(); s.appendChild(e);}); }
  else {
    const box=document.createElement("div"); box.className="tablebox"; const t=document.createElement("table");
    const hr=t.createTHead().insertRow(); b.columns.forEach(c=>{const th=document.createElement("th"); th.textContent=c; hr.appendChild(th);});
    const tb=t.createTBody(); b.rows.forEach(r=>{const tr=tb.insertRow(); r.forEach(c=>{const td=tr.insertCell();
      if (typeof c==="number"){ td.className="num"; td.textContent=fmt(c,4);} else td.textContent=c ?? ""; });});
    box.appendChild(t); s.appendChild(box);
    const cap=document.createElement("div"); cap.className="caption";
    const ct=document.createElement("span"); ct.textContent=b.caption||""; const btn=document.createElement("button");
    btn.type="button"; btn.id="copy"+i; btn.textContent="Kopiera tabell"; const toast=document.createElement("span"); toast.className="toast";
    btn.addEventListener("click", async () => {
      const tsv=[b.columns, ...b.rows.map(r=>r.map(c=>typeof c==="number"?fmt(c,4):(c??"")))].map(r=>r.join("\t")).join("\n");
      try { await navigator.clipboard.writeText(tsv); toast.textContent="Kopierad, klistra in i PowerPoint eller Excel"; }
      catch(e){ const r=document.createRange(); r.selectNodeContents(t); const sel=getSelection(); sel.removeAllRanges(); sel.addRange(r);
        toast.textContent="Markerad, tryck Ctrl+C"; }
      setTimeout(()=>toast.textContent="",3500);
    });
    const right=document.createElement("span"); right.className="row"; right.style.alignItems="center"; right.append(toast,btn);
    cap.append(ct,right); s.appendChild(cap);
  }
  $("blocks").appendChild(s);
});

// ---------- 3D
if (!D.mesh || !window.THREE){ $("stage").innerHTML='<div class="nodata">3D-vyn kunde inte laddas (inget nät eller three.js blockerat). Tabellerna nedan gäller ändå.</div>'; }
else (function(){
  const P = new Float32Array(D.mesh.p), N = P.length/3, T = new Uint32Array(D.mesh.tri);
  const cv=$("cv"), renderer=new THREE.WebGLRenderer({canvas:cv,antialias:true,preserveDrawingBuffer:true});
  renderer.setPixelRatio(Math.min(2,window.devicePixelRatio||1)); renderer.localClippingEnabled=true;
  const clipP=new THREE.Plane(new THREE.Vector3(0,0,-1),1e9), CL=[clipP];
  const scene=new THREE.Scene(), cam=new THREE.PerspectiveCamera(35,1.6,1e-3,1e7);
  scene.add(new THREE.HemisphereLight(0xffffff,0x6b7480,0.85)); const dl=new THREE.DirectionalLight(0xffffff,0.45); cam.add(dl); scene.add(cam);
  const geo=new THREE.BufferGeometry(); const pos=new Float32Array(P); const col=new Float32Array(N*3);
  geo.setAttribute("position",new THREE.BufferAttribute(pos,3)); geo.setAttribute("color",new THREE.BufferAttribute(col,3));
  geo.setIndex(new THREE.BufferAttribute(T,1)); geo.computeBoundingSphere();
  const mat=new THREE.MeshLambertMaterial({vertexColors:true,side:D.mesh.double?THREE.DoubleSide:THREE.FrontSide,clippingPlanes:CL,polygonOffset:true,polygonOffsetFactor:1,polygonOffsetUnits:1});
  const mesh=new THREE.Mesh(geo,mat); scene.add(mesh);
  const wireGeo=()=>{ if(!D.mesh.quads) return new THREE.WireframeGeometry(geo);
    const Q=D.mesh.quads, E=[]; for(let i=0;i<Q.length;i+=4) for(let k=0;k<4;k++) E.push(Q[i+k],Q[i+(k+1)%4]);
    const g=new THREE.BufferGeometry(); g.setAttribute("position",geo.attributes.position); g.setIndex(E); return g; };
  const wire=new THREE.LineSegments(wireGeo(),new THREE.LineBasicMaterial({transparent:true,opacity:0.35,clippingPlanes:CL}));
  wire.visible=false; scene.add(wire);
  const bs=geo.boundingSphere, R=bs.radius||1;
  const marker=new THREE.Mesh(new THREE.SphereGeometry(R*0.025,16,12),new THREE.MeshBasicMaterial({color:0x000000})); scene.add(marker);
  const ctl=new THREE.OrbitControls(cam,cv); ctl.target.copy(bs.center); ctl.enableDamping=true;
  // betraktarvy: F = mot betraktaren (produktens framsida), Uv = upp, S = betraktarens högra sida
  const F=new THREE.Vector3(...(D.front||[0,-1,0])), Uv=new THREE.Vector3(...(D.up||[0,0,1])), S=new THREE.Vector3().crossVectors(Uv,F);
  function view(k){ const c=bs.center, d=R*3.2;
    const v={iso:F.clone().add(S.clone().multiplyScalar(-0.45)).add(Uv.clone().multiplyScalar(0.45)), front:F.clone(),
             top:Uv.clone().add(F.clone().multiplyScalar(0.02)), side:S.clone(), back:F.clone().negate()}[k].normalize();
    cam.position.set(c.x+d*v.x,c.y+d*v.y,c.z+d*v.z); cam.up.copy(k==="top"?F.clone().negate():Uv); ctl.update(); }
  document.querySelectorAll("[data-view]").forEach(b=>b.addEventListener("click",()=>{view(b.dataset.view);
    document.querySelectorAll("[data-view]").forEach(x=>x.classList.toggle("primary",x===b));}));
  view("iso");

  const sc=$("case"), sf=$("fld"), ss=$("scale");
  D.cases.forEach((c,i)=>sc.add(new Option(c.name,i)));
  // CFD-fall: snitt (|U|), vägg (p), strömlinjer; strukturen visas som genomskinlig kontur
  const ghost=new THREE.MeshLambertMaterial({color:0x9aa7b4,transparent:true,opacity:0.16,depthWrite:false,side:THREE.DoubleSide,clippingPlanes:CL});
  function tmesh(o,basic){ const g=new THREE.BufferGeometry(); const p=new Float32Array(o.p); g.setAttribute("position",new THREE.BufferAttribute(p,3));
    g.setAttribute("color",new THREE.BufferAttribute(new Float32Array(p.length),3)); g.setIndex(new THREE.BufferAttribute(new Uint32Array(o.tri),1));
    g.computeVertexNormals(); g.userData.v=o.values; return new THREE.Mesh(g, basic? new THREE.MeshBasicMaterial({vertexColors:true,side:THREE.DoubleSide,clippingPlanes:CL})
      : new THREE.MeshLambertMaterial({vertexColors:true,side:THREE.DoubleSide,clippingPlanes:CL})); }
  const FL=(D.flows||[]).map((f,i)=>{ sc.add(new Option("CFD: "+f.name,"f"+i)); const grp=new THREE.Group(); grp.visible=false; scene.add(grp);
    const sl=f.slices.map(o=>{const m=tmesh(o,true); grp.add(m); return m;}); const wall=f.wall? tmesh(f.wall,false):null; if(wall) grp.add(wall);
    const pos=[], v=[]; f.stream.forEach(L=>{ for(let k=0;k<L.values.length-1;k++){ pos.push(...L.p.slice(3*k,3*k+6)); v.push(L.values[k],L.values[k+1]); }});
    const g=new THREE.BufferGeometry(); g.setAttribute("position",new THREE.Float32BufferAttribute(pos,3)); g.setAttribute("color",new THREE.Float32BufferAttribute(new Float32Array(pos.length),3));
    g.userData.v=v; const st=new THREE.LineSegments(g,new THREE.LineBasicMaterial({vertexColors:true,clippingPlanes:CL})); grp.add(st);
    return {f,grp,sl,wall,st}; });
  const isFlow=()=>String(sc.value).startsWith("f");
  function fillFields(){ const prev=sf.selectedOptions[0]?.text; sf.innerHTML="";
    if (isFlow()){ const o=FL[+sc.value.slice(1)].f; if(o.slices.length) sf.add(new Option("|U| i snitt ["+o.unitU+"]","U")); if(o.wall) sf.add(new Option("p på vägg ["+o.unitp+"]","p")); }
    else D.cases[+sc.value].fields.forEach((f,i)=>sf.add(new Option(f.name+" ["+f.unit+"]",i)));
    const j=[...sf.options].findIndex(o=>o.text===prev); if(j>=0) sf.selectedIndex=j; }
  function autoScale(c){ if(!c.u) return 0; let um=0; for(let i=0;i<c.u.length;i+=3) um=Math.max(um,Math.hypot(c.u[i],c.u[i+1],c.u[i+2]));
    return um>0? 0.08*2*R/um : 0; }
  function themeColors(){ if(!$("tl").hidden) update(); renderer.setClearColor(new THREE.Color(css("--canvas")||"#eef1f4"));
    wire.material.color.set(css("--mesh")||"#20303f"); marker.material.color.set(css("--ink")||"#16202a"); }
  function paint(g,vals,lo,hi){ const c=g.attributes.color.array; for(let i=0;i<vals.length;i++) c.set(turbo((vals[i]-lo)/(hi-lo)),3*i); g.attributes.color.needsUpdate=true; }
  function legend(lo,hi,unit,limit,label){
    const stops=[]; for(let k=0;k<=10;k++){const [r,g,b]=turbo(k/10); stops.push(`rgb(${r*255|0},${g*255|0},${b*255|0}) ${k*10}%`);}
    $("bar").style.background=`linear-gradient(to top, ${stops.join(",")})`;
    const tk=$("ticks"); tk.querySelectorAll("span:not(.unit)").forEach(e=>e.remove()); $("unit").textContent=unit;
    for(let k=0;k<=5;k++){ const v=lo+(hi-lo)*k/5, e=document.createElement("span"); e.style.bottom=(k*20)+"%"; e.textContent=fmt(v,3); tk.appendChild(e); }
    if (limit!=null && limit>=lo && limit<=hi){ const e=document.createElement("span"); e.className="lim";
      e.style.bottom=(100*(limit-lo)/(hi-lo))+"%"; e.textContent="◂ "+(label||"gräns")+" "+fmt(limit,3); tk.appendChild(e); } }
  function clip(){ const ax=$("clipax").value; ["cliprow","clipfrow"].forEach(id=>$(id).hidden=ax==="");
    if(ax===""){ clipP.normal.set(0,0,-1); clipP.constant=1e9; return; }
    const a=+ax, bb=new THREE.Box3().setFromBufferAttribute(geo.attributes.position); const lo=bb.min.getComponent(a), hi=bb.max.getComponent(a);
    const v=lo+(hi-lo)*(+$("clip").value)/1000, n=[0,0,0]; n[a]=$("clipf").checked?1:-1;
    // standardsida: behåll halvan bort från betraktaren så att snittet syns framifrån
    if (Math.abs([F.x,F.y,F.z][a])>0.5 && [F.x,F.y,F.z][a]<0) n[a]=-n[a];
    clipP.normal.set(...n); clipP.constant=-n[a]*v+1e-3*(hi-lo); $("clipv").textContent="xyz"[a]+" = "+fmt(v,4); }
  const tsl=$("tsl"); let playT=null;
  const PAL=["--accent","--warn","--ink","--muted"];
  function chart(f,ti){ const svg=$("chart"); if(!f.curves||!f.curves.length){ svg.innerHTML=""; svg.style.display="none"; return; }
    svg.style.display=""; const W=900,H=250,l=58,r=150,tp=30,b=34; let t0=Infinity,t1=-Infinity,y0=Infinity,y1=-Infinity;
    f.curves.forEach(c=>{c.t.forEach(v=>{t0=Math.min(t0,v);t1=Math.max(t1,v);}); c.y.forEach(v=>{y0=Math.min(y0,v);y1=Math.max(y1,v);});});
    const pad=(y1-y0)*0.05||1; y0-=pad; y1+=pad; const X=v=>l+(W-l-r)*(v-t0)/((t1-t0)||1), Y=v=>H-b-(H-tp-b)*(v-y0)/((y1-y0)||1);
    let g=""; for(let k=0;k<=5;k++){ const yv=y0+(y1-y0)*k/5, tv=t0+(t1-t0)*k/5;
      g+=`<line class="ax" x1="${l}" x2="${W-r}" y1="${Y(yv)}" y2="${Y(yv)}"/><text class="lab" x="${l-6}" y="${Y(yv)+4}" text-anchor="end">${fmt(yv,3)}</text>`;
      g+=`<text class="lab" x="${X(tv)}" y="${H-b+16}" text-anchor="middle">${fmt(tv,3)}</text>`; }
    g+=`<text class="lab" x="${(l+W-r)/2}" y="${H-4}" text-anchor="middle">tid [${f.t_unit}]</text><text class="lab" x="4" y="${tp-14}">${f.curve_title?f.curve_title+" ":""}[${f.curve_unit}]</text>`;
    const tc=f.t[ti];
    f.curves.forEach((c,i)=>{ const col=css(PAL[i%PAL.length]); let d=""; c.t.forEach((v,j)=>{ d+=(j?"L":"M")+X(v).toFixed(1)+","+Y(c.y[j]).toFixed(1); });
      let j=0; while(j<c.t.length-1 && c.t[j]<tc) j++; const yc=c.y[j];
      g+=`<path d="${d}" fill="none" stroke="${col}" stroke-width="2"/><circle cx="${X(c.t[j])}" cy="${Y(yc)}" r="4" fill="${col}"/>`;
      g+=`<text x="${W-r+8}" y="${tp+14+18*i}" fill="${col}">${c.name}: ${fmt(yc,3)}</text>`; });
    g+=`<line class="cur" x1="${X(tc)}" x2="${X(tc)}" y1="${tp}" y2="${H-b}"/>`; svg.innerHTML=g; }
  $("play").addEventListener("click",()=>{ if(playT){ clearInterval(playT); playT=null; $("play").textContent="▶ Spela"; return; }
    $("play").textContent="❚❚ Paus"; if(+tsl.value>=+tsl.max) tsl.value=0;
    playT=setInterval(()=>{ if(+tsl.value>=+tsl.max){ clearInterval(playT); playT=null; $("play").textContent="▶ Spela"; return; }
      tsl.value=+tsl.value+1; update(); },120); });
  tsl.addEventListener("input",update);
  function update(){
    clip();
    const flow=isFlow(); $("strmrow").hidden=!flow; FL.forEach(o=>o.grp.visible=false);
    if (flow){
      const o=FL[+sc.value.slice(1)], f=o.f; o.grp.visible=true; mesh.material=ghost; marker.visible=false; ss.disabled=true; $("scalev").textContent="–";
      for(let i=0;i<3*N;i++) pos[i]=P[i]; geo.attributes.position.needsUpdate=true; geo.computeVertexNormals();
      const U=sf.value==="U"; o.sl.forEach(m=>m.visible=U); if(o.wall) o.wall.visible=!U; o.st.visible=$("strm").checked;
      let lo=0, hi=f.Umax||1;
      if (U){ o.sl.forEach(m=>paint(m.geometry,m.geometry.userData.v,lo,hi)); legend(lo,hi,f.unitU); }
      else { lo=f.wall.min; hi=f.wall.max>lo?f.wall.max:lo+1; paint(o.wall.geometry,o.wall.geometry.userData.v,lo,hi); legend(lo,hi,f.unitp); }
      paint(o.st.geometry,o.st.geometry.userData.v,0,f.Umax||1);
      $("readout").innerHTML=`${U?"|U| max <b>"+fmt(f.Umax,4)+" "+f.unitU+"</b>":"p <b>"+fmt(f.wall.min,4)+" … "+fmt(f.wall.max,4)+" "+f.unitp+"</b>"}${f.time?", iteration "+f.time:""}`;
      return;
    }
    mesh.material=mat;
    const c=D.cases[+sc.value], f=c.fields[+sf.value]; const lo=f.min, hi=f.max>f.min?f.max:f.min+1;
    const tr=!!f.frames; $("tl").hidden=!tr; let vals=f.values, ti=0;
    if (tr){ if(+tsl.max!==f.frames.length-1){ tsl.max=f.frames.length-1; tsl.value=0; } ti=Math.min(+tsl.value,f.frames.length-1);
      vals=f.frames[ti]; $("tval").textContent="t = "+fmt(f.t[ti],3)+" "+f.t_unit; chart(f,ti); }
    paint(geo,vals,lo,hi);
    const auto=autoScale(c); ss.disabled=!c.u; const s=c.u? auto*(+ss.value/25) : 0;
    $("scalev").textContent = c.u ? (s>0? "× "+fmt(s,3)+(+ss.value===25?" (auto)":"") : "av") : "ingen förskjutning";
    for(let i=0;i<3*N;i++) pos[i]=P[i]+(c.u? s*c.u[i]:0);
    geo.attributes.position.needsUpdate=true; geo.computeVertexNormals();
    if(!D.mesh.quads){ wire.geometry.dispose(); wire.geometry=wireGeo(); }
    let im=0; for(let i=1;i<N;i++) if(vals[i]>vals[im]) im=i; marker.position.set(pos[3*im],pos[3*im+1],pos[3*im+2]);
    marker.visible=$("maxm").checked;
    legend(lo,hi,f.unit,f.limit,f.limit_label);
    const util = f.limit? `, utnyttjande <b>${fmt(f.max/f.limit,3)}</b>`:"";
    if (tr){ let vmx=-Infinity,vmn=Infinity; for(const v of vals){ if(v>vmx) vmx=v; if(v<vmn) vmn=v; }
      $("readout").innerHTML=`t = <b>${fmt(f.t[ti],3)} ${f.t_unit}</b>: max <b>${fmt(vmx,4)} ${f.unit}</b> vid (${[pos[3*im],pos[3*im+1],pos[3*im+2]].map(x=>fmt(x,4)).join("; ")}) mm, min <b>${fmt(vmn,4)}</b>`; }
    else $("readout").innerHTML=`max <b>${fmt(f.max,4)} ${f.unit}</b> vid (${f.at.map(x=>fmt(x,4)).join("; ")}) mm, min <b>${fmt(f.min,4)}</b>${util}`;
  }
  ["clip","clipf","strm"].forEach(id=>$(id).addEventListener("input",update));
  $("clipax").addEventListener("input",()=>{ const ax=$("clipax").value; if(ax!==""){ const bb=new THREE.Box3().setFromBufferAttribute(geo.attributes.position);
    const lo=bb.min.getComponent(+ax), hi=bb.max.getComponent(+ax); $("clip").value = (lo<0&&hi>0)? Math.round(1000*(0-lo)/(hi-lo)) : 500; } update(); });
  sc.addEventListener("change",()=>{fillFields();update();}); sf.addEventListener("change",update);
  ss.addEventListener("input",update); $("edges").addEventListener("change",e=>{wire.visible=e.target.checked;});
  $("maxm").addEventListener("change",update);
  fillFields(); ss.value = 25; update(); themeColors();
  matchMedia("(prefers-color-scheme: dark)").addEventListener("change",themeColors);
  new MutationObserver(themeColors).observe(document.documentElement,{attributes:true,attributeFilter:["data-theme"]});
  function size(){ const w=cv.clientWidth,h=cv.clientHeight; if(cv.width!==w||cv.height!==h){renderer.setSize(w,h,false); cam.aspect=w/h; cam.updateProjectionMatrix();} }
  (function loop(){ size(); ctl.update(); renderer.render(scene,cam); requestAnimationFrame(loop); })();
})();
</script>
"""
