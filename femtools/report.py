"""Byggstenar för FEM-rapportens figurer (verktygslåda, inte en fast rapport).

Verifierat (se docs/verifiering.md): nodal() med L2-projektion och sample() för snittuttag.
Plottfunktionerna ger enhetlig stil (ljus bakgrund, enheter, axelnamn, färgskala 0..Rp0,2,
förstoringsfaktor). Fig. 1–7 i skillen fem-rapport är ett minimum; agenten lägger till de
figurer som behövs för resultatet eller beslutet och kombinerar byggstenarna fritt.

    from femtools import static, report
    r = static(m, E, nu, traction)
    n = report.nodal(r, E, nu)                       # spänningar och u i hörnnoder (L2-projektion)
    report.plot_field(m, n["vm"], "fig/3_vm.png", "von Mises, LF1", vmax=Rp02)
    x, pts = report.line([0, 5, 0], [0, 5, 10], 50)  # snittlinje genom tjockleken
    report.plot_line(x, {"σ_xx": report.sample(m, n["sxx"], pts)}, "fig/5_snitt.svg", "Snitt A–A")
"""
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from mpl_toolkits.mplot3d.art3d import Poly3DCollection  # noqa: E402

COMP = {"sxx": (0, 0), "syy": (1, 1), "szz": (2, 2), "sxy": (0, 1), "syz": (1, 2), "sxz": (0, 2)}


# ---------------------------------------------------------------- data
def nodal(res, E, nu):
    """Spänningskomponenter, von Mises, huvudspänningar och förskjutning i nätets hörnnoder.
    Spänningen L2-projiceras till P1 (medelvärde av elementvärden ger ca 10 % för lågt vid ytan)."""
    from skfem import Basis, ElementTetP1
    sb = res["basis"]
    m = sb.mesh
    pb = Basis(m, ElementTetP1(), quadrature=sb.quadrature)
    s = res["sigma"]
    out = {k: pb.project(s[i, j]) for k, (i, j) in COMP.items()}
    S = np.empty((m.p.shape[1], 3, 3))
    for k, (i, j) in COMP.items():
        S[:, i, j] = S[:, j, i] = out[k]
    out["principal"] = np.linalg.eigvalsh(S)[:, ::-1].T          # s1 >= s2 >= s3
    out["vm"] = np.sqrt(0.5 * ((out["sxx"] - out["syy"])**2 + (out["syy"] - out["szz"])**2
                               + (out["szz"] - out["sxx"])**2)
                        + 3 * (out["sxy"]**2 + out["syz"]**2 + out["sxz"]**2))
    out["u"] = res["u"][sb.nodal_dofs]                             # (3, N) i hörnnoder
    out["umag"] = np.linalg.norm(out["u"], axis=0)
    return out


def locate(m, pts, k=40):
    """Element och barycentriska koordinater för punkter (3, n). Punkter utanför nätet: element -1."""
    from scipy.spatial import cKDTree
    pts = np.atleast_2d(np.asarray(pts, float))
    if pts.shape[0] != 3:
        pts = pts.T
    P = m.p[:, m.t]                                                # (3, 4, nel)
    T = np.transpose(P[:, 1:] - P[:, :1], (2, 0, 1))               # (nel, 3, 3)
    Tinv = np.linalg.inv(T)
    tree = cKDTree(P.mean(axis=1).T)
    _, cand = tree.query(pts.T, k=min(k, m.t.shape[1]))
    cand = np.atleast_2d(cand)
    n = pts.shape[1]
    elem = -np.ones(n, int)
    lam = np.zeros((4, n))
    todo = np.arange(n)
    for c in range(cand.shape[1]):
        e = cand[todo, c]
        l123 = np.einsum("nij,jn->in", Tinv[e], pts[:, todo] - P[:, 0, e])
        l = np.vstack([1 - l123.sum(0), l123])
        ok = l.min(0) >= -1e-9
        elem[todo[ok]], lam[:, todo[ok]] = e[ok], l[:, ok]
        todo = todo[~ok]
        if todo.size == 0:
            break
    return elem, lam


def sample(m, values, pts):
    """Interpolera nodvärden (N,) eller (c, N) i punkter (3, n). NaN utanför nätet."""
    elem, lam = locate(m, pts)
    v = np.asarray(values)
    inside = elem >= 0
    vt = v[..., m.t[:, np.where(inside, elem, 0)]]                 # (..., 4, n)
    out = np.einsum("...kn,kn->...n", vt, lam)
    out[..., ~inside] = np.nan
    return out


def line(p0, p1, n=50):
    """Punkter längs en snittlinje: (avstånd [mm], punkter (3, n))."""
    p0, p1 = np.asarray(p0, float), np.asarray(p1, float)
    s = np.linspace(0, 1, n)
    return s * np.linalg.norm(p1 - p0), p0[:, None] + (p1 - p0)[:, None] * s


def surface(m):
    """Ytrianglar (3, nf) på nätets rand."""
    return m.facets[:, m.boundary_facets()]


# ---------------------------------------------------------------- stil
def style():
    plt.rcParams.update({"figure.facecolor": "white", "axes.facecolor": "white", "font.size": 9,
                         "axes.grid": True, "grid.alpha": 0.3, "savefig.bbox": "tight"})


def save(fig, path, dpi=200):
    """PNG (fältplottar, minst 1600 px bred) eller SVG (diagram) efter filändelse."""
    fig.savefig(path, dpi=dpi)
    plt.close(fig)
    return path


def _ax3d(m, view, figsize=(8, 5.5)):
    style()
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(projection="3d")
    ax.view_init(*view)
    ax.set_box_aspect(np.ptp(m.p, axis=1) + 1e-9)                  # verkliga proportioner
    ax.set_xlabel("X [mm]"); ax.set_ylabel("Y [mm]"); ax.set_zlabel("Z [mm]")
    ax.grid(False)
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.set_pane_color((1, 1, 1, 0))
        axis.set_major_locator(plt.MaxNLocator(4))
    return fig, ax


def _cmap(over="0.55"):
    c = plt.get_cmap("turbo").copy()
    c.set_over(over)
    return c


# ---------------------------------------------------------------- figurer
def plot_field(m, values, path, title, unit="MPa", vmin=0.0, vmax=None, u=None, scale=0.0,
               view=(25, -60), edges=False, mark_max=True, mask=None, figsize=(8, 5.5)):
    """Nodfält på ytan. vmax = gräns (t.ex. Rp0,2): värden över visas grå. u, scale: deformerad form.
    mask: bool (N,) för noder som ska uteslutas ur max (t.ex. kända singulariteter)."""
    p = m.p + (scale * u if u is not None else 0)
    tri = surface(m)
    fig, ax = _ax3d(m, view, figsize)
    vals = np.asarray(values)
    vmax = float(np.nanmax(vals)) if vmax is None else vmax
    norm = plt.Normalize(vmin, vmax)
    pc = Poly3DCollection(p[:, tri].transpose(2, 1, 0), cmap=_cmap(), norm=norm,
                          edgecolor="k" if edges else "none", linewidth=0.1)
    pc.set_array(vals[tri].mean(axis=0))
    ax.add_collection3d(pc)
    ax.auto_scale_xyz(*[(a.min(), a.max()) for a in p])
    cb = fig.colorbar(pc, ax=ax, shrink=0.7, pad=0.12, extend="max" if vals.max() > vmax else "neither")
    cb.set_label(unit)
    if mark_max:
        v = np.where(mask, -np.inf, vals) if mask is not None else vals
        i = int(np.nanargmax(v))
        ax.scatter(*p[:, i], color="k", s=25, depthshade=False)
        ax.text(*p[:, i], f"  max {vals[i]:.3g} {unit}\n  ({p[0, i]:.1f}, {p[1, i]:.1f}, {p[2, i]:.1f})",
                fontsize=8)
    ax.set_title(title + (f"  (förstoring {scale:.3g}x)" if u is not None and scale else ""))
    return save(fig, path)


def plot_deformation(m, n, path, title, scale=None, view=(25, -60), limit=None):
    """Deformerad form färgad med |u| [mm] över odeformerad kontur. scale=None: max 10 % av största mått.
    limit: deformationskrav [mm] i titeln."""
    umax = float(n["umag"].max())
    if scale is None:
        scale = 0.1 * np.ptp(m.p, axis=1).max() / max(umax, 1e-30)
    t = title + (f", krav {limit:.3g} mm, utnyttjat {umax / limit:.0%}" if limit else "")
    tri = surface(m)
    p = m.p + scale * n["u"]
    fig, ax = _ax3d(m, view)
    ax.add_collection3d(Poly3DCollection(m.p[:, tri].transpose(2, 1, 0), facecolor="0.85",
                                         edgecolor="0.6", linewidth=0.05, alpha=0.15))
    norm = plt.Normalize(0, umax)
    pc = Poly3DCollection(p[:, tri].transpose(2, 1, 0), cmap=_cmap(), norm=norm, edgecolor="none")
    pc.set_array(n["umag"][tri].mean(axis=0))
    ax.add_collection3d(pc)
    allp = np.hstack([m.p, p])
    ax.auto_scale_xyz(*[(a.min(), a.max()) for a in allp])
    fig.colorbar(pc, ax=ax, shrink=0.7, pad=0.12).set_label("|u| [mm]")
    i = int(n["umag"].argmax())
    ax.scatter(*p[:, i], color="k", s=25, depthshade=False)
    ax.text(*p[:, i], f"  max {umax:.3g} mm", fontsize=8)
    ax.set_title(f"{t}  (förstoring {scale:.3g}x)")
    return save(fig, path)


def plot_bc(m, path, title, fixed=(), loads=(), view=(25, -60), sym=()):
    """Skiss över randvillkor och laster från analyserad geometri (fig. 1).
    fixed: [(facet-index, "Ux Uy Uz")], loads: [(facet-index, vektor [N], "LF1 1000 N")],
    sym: [(facet-index, "symmetri y=0")]. Facet-index t.ex. m.boundaries["fix"]."""
    fig, ax = _ax3d(m, view)
    tri = surface(m)
    ax.add_collection3d(Poly3DCollection(m.p[:, tri].transpose(2, 1, 0), facecolor="0.9",
                                         edgecolor="0.7", linewidth=0.05, alpha=0.35))
    L = np.ptp(m.p, axis=1).max()
    for group, col in ((fixed, "tab:blue"), (sym, "tab:green")):
        for f, lab in group:
            ax.add_collection3d(Poly3DCollection(m.p[:, m.facets[:, f]].transpose(2, 1, 0),
                                                 facecolor=col, alpha=0.6, edgecolor="none"))
            c = m.p[:, m.facets[:, f]].mean(axis=(1, 2))
            ax.text(*c, "  " + lab, color=col, fontsize=8)
    for f, vec, lab in loads:
        fc = m.p[:, m.facets[:, f]]
        ax.add_collection3d(Poly3DCollection(fc.transpose(2, 1, 0), facecolor="tab:red", alpha=0.5,
                                             edgecolor="none"))
        c = fc.mean(axis=(1, 2))
        d = np.asarray(vec, float)
        d = 0.25 * L * d / np.linalg.norm(d)
        ax.quiver(*(c - d), *d, color="tab:red", linewidth=2, arrow_length_ratio=0.25)
        ax.text(*(c - d), lab + "  ", color="tab:red", fontsize=8, ha="right")
    ax.auto_scale_xyz(*[(a.min(), a.max()) for a in m.p])
    ax.set_title(title)
    return save(fig, path)


def plot_plane(m, values, path, title, origin, normal, unit="MPa", vmin=0.0, vmax=None, res=200):
    """Plansnitt genom nätet: fältet samplas på ett rutnät i planet (NaN utanför materialet).
    Samma vmin/vmax som översiktsbilden ger jämförbar färgskala."""
    o, nrm = np.asarray(origin, float), np.asarray(normal, float) / np.linalg.norm(normal)
    i1, i2 = sorted(np.argsort(abs(nrm))[:2])                     # två globala axlar närmast planet
    e1 = np.eye(3)[i1] - nrm[i1] * nrm; e1 /= np.linalg.norm(e1)
    e2 = np.eye(3)[i2] - nrm[i2] * nrm - (np.eye(3)[i2] @ e1) * e1; e2 /= np.linalg.norm(e2)
    q = m.p - o[:, None]
    r1, r2 = e1 @ q, e2 @ q
    g1, g2 = np.meshgrid(np.linspace(r1.min(), r1.max(), res), np.linspace(r2.min(), r2.max(), res))
    pts = o[:, None] + e1[:, None] * g1.ravel() + e2[:, None] * g2.ravel()
    v = sample(m, values, pts).reshape(g1.shape)
    style()
    ext = np.ptp(g1) / max(np.ptp(g2), 1e-9)
    fig, ax = plt.subplots(figsize=(8, float(np.clip(8 / ext, 2.5, 8))))
    vmax = float(np.nanmax(v)) if vmax is None else vmax
    im = ax.pcolormesh(g1, g2, np.ma.masked_invalid(v), cmap=_cmap(), vmin=vmin, vmax=vmax, shading="auto")
    fig.colorbar(im, ax=ax, fraction=0.03, extend="max" if np.nanmax(v) > vmax else "neither").set_label(unit)
    def _lab(e):
        t = [("" if abs(abs(c) - 1) < 1e-9 else f"{abs(c):.2g}") + "XYZ"[i] for i, c in enumerate(e) if abs(c) > 1e-9]
        return "+".join(t) + " [mm]"
    ax.set_xlabel(_lab(e1)); ax.set_ylabel(_lab(e2))
    ax.set_aspect("equal")
    i = np.nanargmax(v)
    ax.plot(g1.ravel()[i], g2.ravel()[i], "k+", ms=10)
    ax.set_title(f"{title}  (max {np.nanmax(v):.3g} {unit})")
    return save(fig, path)


def plot_line(dist, curves, path, title, unit="MPa", ref=None, xlabel="Avstånd längs snitt [mm]"):
    """Värden längs snittlinje. curves: {"namn": värden}. ref: {"namn": värden eller konstant} (handräkning)."""
    style()
    fig, ax = plt.subplots(figsize=(7, 4))
    for k, v in curves.items():
        ax.plot(dist, v, lw=2, label=k)
    for k, v in (ref or {}).items():
        ax.plot(dist, np.broadcast_to(v, np.shape(dist)), "k--", lw=1, label=k)
    ax.set_xlabel(xlabel); ax.set_ylabel(unit); ax.set_title(title); ax.legend()
    return save(fig, path)


def plot_convergence(dofs, values, path, title, unit="MPa", f_ext=None, tol=0.03):
    """Utvärderingsstorhet mot DOF, Richardson-värde som linje och ±tol-band."""
    style()
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.semilogx(dofs, values, "o-", lw=2, label="FEM")
    if f_ext is not None:
        ax.axhline(f_ext, color="k", ls="--", lw=1, label=f"Richardson {f_ext:.4g}")
        ax.axhspan(f_ext * (1 - tol), f_ext * (1 + tol), color="0.85", label=f"±{tol:.0%}")
    ax.set_xlabel("DOF"); ax.set_ylabel(unit); ax.set_title(title); ax.legend()
    return save(fig, path)
