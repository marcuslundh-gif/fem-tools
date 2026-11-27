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


class Report:
    def __init__(self, title, meta=None, subtitle=None):
        self.title, self.subtitle = title, subtitle
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
        return dict(title=self.title, subtitle=self.subtitle, meta=self.meta,
                    mesh=None if s is None else dict(p=np.round(s["p"], 4).ravel().tolist(),
                                                     tri=s["tri"].ravel().astype(int).tolist()),
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
      <div class="row" role="group" aria-label="Vy">
        <button type="button" data-view="iso" class="primary">Iso</button>
        <button type="button" data-view="x">X</button><button type="button" data-view="y">Y</button><button type="button" data-view="z">Z</button>
      </div>
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
  renderer.setPixelRatio(Math.min(2,window.devicePixelRatio||1));
  const scene=new THREE.Scene(), cam=new THREE.PerspectiveCamera(35,1.6,1e-3,1e7);
  scene.add(new THREE.HemisphereLight(0xffffff,0x6b7480,0.85)); const dl=new THREE.DirectionalLight(0xffffff,0.45); cam.add(dl); scene.add(cam);
  const geo=new THREE.BufferGeometry(); const pos=new Float32Array(P); const col=new Float32Array(N*3);
  geo.setAttribute("position",new THREE.BufferAttribute(pos,3)); geo.setAttribute("color",new THREE.BufferAttribute(col,3));
  geo.setIndex(new THREE.BufferAttribute(T,1)); geo.computeBoundingSphere();
  const mat=new THREE.MeshLambertMaterial({vertexColors:true,side:THREE.FrontSide,polygonOffset:true,polygonOffsetFactor:1,polygonOffsetUnits:1});
  const mesh=new THREE.Mesh(geo,mat); scene.add(mesh);
  const wire=new THREE.LineSegments(new THREE.WireframeGeometry(geo),new THREE.LineBasicMaterial({transparent:true,opacity:0.35}));
  wire.visible=false; scene.add(wire);
  const bs=geo.boundingSphere, R=bs.radius||1;
  const marker=new THREE.Mesh(new THREE.SphereGeometry(R*0.025,16,12),new THREE.MeshBasicMaterial({color:0x000000})); scene.add(marker);
  const ctl=new THREE.OrbitControls(cam,cv); ctl.target.copy(bs.center); ctl.enableDamping=true;
  function view(k){ const c=bs.center, d=R*3.2; const v={iso:[1,-1.2,0.9],x:[1,0,0],y:[0,-1,0],z:[0,0,1]}[k];
    const n=Math.hypot(...v); cam.position.set(c.x+d*v[0]/n,c.y+d*v[1]/n,c.z+d*v[2]/n); cam.up.set(0,0,1); if(k==="z") cam.up.set(0,1,0); ctl.update(); }
  document.querySelectorAll("[data-view]").forEach(b=>b.addEventListener("click",()=>{view(b.dataset.view);
    document.querySelectorAll("[data-view]").forEach(x=>x.classList.toggle("primary",x===b));}));
  view("iso");

  const sc=$("case"), sf=$("fld"), ss=$("scale");
  D.cases.forEach((c,i)=>sc.add(new Option(c.name,i)));
  function fillFields(){ const c=D.cases[+sc.value]; const prev=sf.selectedOptions[0]?.text; sf.innerHTML="";
    c.fields.forEach((f,i)=>sf.add(new Option(f.name+" ["+f.unit+"]",i)));
    const j=c.fields.findIndex(f=>f.name+" ["+f.unit+"]"===prev); if(j>=0) sf.value=j; }
  function autoScale(c){ if(!c.u) return 0; let um=0; for(let i=0;i<c.u.length;i+=3) um=Math.max(um,Math.hypot(c.u[i],c.u[i+1],c.u[i+2]));
    return um>0? 0.08*2*R/um : 0; }
  function themeColors(){ renderer.setClearColor(new THREE.Color(css("--canvas")||"#eef1f4"));
    wire.material.color.set(css("--mesh")||"#20303f"); marker.material.color.set(css("--ink")||"#16202a"); }
  function update(){
    const c=D.cases[+sc.value], f=c.fields[+sf.value]; const lo=f.min, hi=f.max>f.min?f.max:f.min+1;
    for(let i=0;i<N;i++){ const rgb=turbo((f.values[i]-lo)/(hi-lo)); col.set(rgb,3*i); }
    geo.attributes.color.needsUpdate=true;
    const auto=autoScale(c); ss.disabled=!c.u; const s=c.u? auto*(+ss.value/25) : 0;
    $("scalev").textContent = c.u ? (s>0? "× "+fmt(s,3)+(+ss.value===25?" (auto)":"") : "av") : "ingen förskjutning";
    for(let i=0;i<3*N;i++) pos[i]=P[i]+(c.u? s*c.u[i]:0);
    geo.attributes.position.needsUpdate=true; geo.computeVertexNormals();
    wire.geometry.dispose(); wire.geometry=new THREE.WireframeGeometry(geo);
    let im=0; for(let i=1;i<N;i++) if(f.values[i]>f.values[im]) im=i; marker.position.set(pos[3*im],pos[3*im+1],pos[3*im+2]);
    marker.visible=$("maxm").checked;
    // färgskala
    const stops=[]; for(let k=0;k<=10;k++){const [r,g,b]=turbo(k/10); stops.push(`rgb(${r*255|0},${g*255|0},${b*255|0}) ${k*10}%`);}
    $("bar").style.background=`linear-gradient(to top, ${stops.join(",")})`;
    const tk=$("ticks"); tk.querySelectorAll("span:not(.unit)").forEach(e=>e.remove()); $("unit").textContent=f.unit;
    for(let k=0;k<=5;k++){ const v=lo+(hi-lo)*k/5, e=document.createElement("span"); e.style.bottom=(k*20)+"%"; e.textContent=fmt(v,3); tk.appendChild(e); }
    if (f.limit!=null && f.limit>=lo && f.limit<=hi){ const e=document.createElement("span"); e.className="lim";
      e.style.bottom=(100*(f.limit-lo)/(hi-lo))+"%"; e.textContent="◂ "+(f.limit_label||"gräns")+" "+fmt(f.limit,3); tk.appendChild(e); }
    const util = f.limit? `, utnyttjande <b>${fmt(f.max/f.limit,3)}</b>`:"";
    $("readout").innerHTML=`max <b>${fmt(f.max,4)} ${f.unit}</b> vid (${f.at.map(x=>fmt(x,4)).join("; ")}) mm, min <b>${fmt(f.min,4)}</b>${util}`;
  }
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
