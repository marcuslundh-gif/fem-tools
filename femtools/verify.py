"""Verifiering mot handräkning. Kör: python -m femtools.verify"""
import time
import numpy as np
from skfem import MeshTet
import gmsh
from .core import static, modal
from .loads import static_ext
from .thermal import heat, heat_transient
from .mesh import setup, refine_and_mesh, to_skfem
from .materials import material, materials
from . import projektblad, report


def cantilever():
    L, B, H, E, nu, F = 100.0, 10.0, 10.0, 210e3, 0.3, -100.0
    m = MeshTet.init_tensor(np.linspace(0, L, 41), np.linspace(0, B, 5), np.linspace(0, H, 5))
    m = m.with_boundaries({"fix": lambda x: np.isclose(x[0], 0), "load": lambda x: np.isclose(x[0], L)})
    r = static(m, E, nu, [0, 0, F / (B * H)])
    uz = r["u"][r["basis"].get_dofs("load").all("u^3")].mean()
    I = B * H**3 / 12
    d_hand = F * L**3 / (3 * E * I)
    n = report.nodal(r, E, nu)                                     # L2-projektion, snitt vid L/2
    x, pts = report.line([L / 2, B / 2, 0], [L / 2, B / 2, H], 11)
    s_sec = report.sample(m, n["sxx"], pts)
    s_hand = -F * L / 2 * (pts[2] - H / 2) / I
    e_sec = np.abs(s_sec - s_hand).max() / np.abs(s_hand).max()
    import tempfile, os
    with tempfile.TemporaryDirectory() as td:                        # figurer skapas utan fel
        report.plot_deformation(m, n, os.path.join(td, "d.png"), "test")
        report.plot_field(m, n["vm"], os.path.join(td, "f.png"), "test", vmax=355)
        report.plot_plane(m, n["sxx"], os.path.join(td, "p.png"), "test", [0, B / 2, 0], [0, 1, 0], res=60)
        report.plot_line(x, {"FEM": s_sec}, os.path.join(td, "l.svg"), "test", ref={"hand": s_hand})
        nfig = len(os.listdir(td))
    f1 = modal(m, E, nu, 7.85e-9, k=2)[0]
    f_hand = 1.8751**2 / (2 * np.pi) * np.sqrt(E * I / (7.85e-9 * B * H * L**4))
    eq = max(abs(a + b) for a, b in zip(r["R"], r["F"])) / abs(F)
    print(f"Konsol: delta {uz:.4f} mm (hand {d_hand:.4f}), f1 {f1:.0f} Hz (hand {f_hand:.0f}), "
          f"jämvikt {eq:.1e}, DOF {r['DOF']}, lösare {r['solver']}")
    print(f"Snitt L/2 (report, L2-projektion): sigma_xx {s_sec[-1]:.2f} MPa (hand {s_hand[-1]:.2f}), "
          f"max avvikelse {e_sec:.1%}, {nfig} figurer")
    assert abs(uz / d_hand - 1) < 0.01 and abs(f1 / f_hand - 1) < 0.01 and eq < 1e-3
    assert e_sec < 0.03 and nfig == 4


def plate_hole():
    L, W, T, d, s_nom = 100.0, 40.0, 5.0, 10.0, 10.0
    setup()
    box = gmsh.model.occ.addBox(0, 0, 0, L, W, T)
    cyl = gmsh.model.occ.addCylinder(L / 2, W / 2, 0, 0, 0, T, d / 2)
    gmsh.model.occ.cut([(3, box)], [(3, cyl)])
    gmsh.model.occ.synchronize()
    hole = [s for _, s in gmsh.model.getEntities(2)
            if abs(gmsh.model.occ.getMass(2, s) - np.pi * d * T) < 1e-3 * np.pi * d * T]
    t_mesh, q = refine_and_mesh(4.0, 0.8, hole)
    m = to_skfem()
    gmsh.finalize()
    m = m.with_boundaries({"fix": lambda x: np.isclose(x[0], 0), "load": lambda x: np.isclose(x[0], L)})
    from skfem import ElementVector, ElementTetP2, Basis
    ib = Basis(m, ElementVector(ElementTetP2()))
    pts = m.p.T
    nd = ib.nodal_dofs
    p0 = np.argmin(np.linalg.norm(pts - [0, 0, 0], axis=1))
    p1 = np.argmin(np.linalg.norm(pts - [0, W, 0], axis=1))
    p2 = np.argmin(np.linalg.norm(pts - [0, 0, T], axis=1))
    D = np.unique(np.concatenate([ib.get_dofs("fix").all("u^1"), [nd[1, p0], nd[2, p0], nd[2, p1], nd[1, p2]]]))
    r = static(m, 210e3, 0.3, [s_nom, 0, 0], D=D)
    sxx = r["sigma"][0, 0].max()
    hand = (2 + (1 - d / W)**3) * s_nom * W / (W - d)
    n = report.nodal(r, 210e3, 0.3)
    s_node = report.sample(m, n["sxx"], np.array([[L / 2], [W / 2 + d / 2], [T / 2]]))[0]
    print(f"Plåt med hål: sigma_xx,max {sxx:.1f} MPa, i hålkant (L2-nod) {s_node:.1f} MPa (Heywood {hand:.1f}), "
          f"gamma_min {q['gamma_min']:.2f}, mesh {t_mesh:.2f} s, DOF {r['DOF']}")
    assert abs(sxx / hand - 1) < 0.06 and abs(s_node / hand - 1) < 0.06


def _tube(a, b, H, h):
    """Kvarts tjockväggigt rör (axel z) med ränder sx, sy, ends, inner, outer."""
    setup()
    o = gmsh.model.occ.addCylinder(0, 0, 0, 0, 0, H, b)
    i = gmsh.model.occ.addCylinder(0, 0, 0, 0, 0, H, a)
    r, _ = gmsh.model.occ.cut([(3, o)], [(3, i)])
    gmsh.model.occ.intersect(r, [(3, gmsh.model.occ.addBox(0, 0, 0, b, b, H))])
    gmsh.model.occ.synchronize()
    refine_and_mesh(h, n_curv=64)
    m = to_skfem()
    gmsh.finalize()
    rr = lambda x: np.hypot(x[0], x[1])
    side = lambda x: ~(np.isclose(x[0], 0) | np.isclose(x[1], 0) | np.isclose(x[2], 0) | np.isclose(x[2], H))
    return m.with_boundaries({"sx": lambda x: np.isclose(x[0], 0), "sy": lambda x: np.isclose(x[1], 0),
                              "ends": lambda x: np.isclose(x[2], 0) | np.isclose(x[2], H),
                              "inner": lambda x: (rr(x) < a + h / 5) & side(x),
                              "outer": lambda x: (rr(x) > b - h / 5) & side(x)})


def _polar(m, n, R, z):
    """(sigma_r, sigma_theta) vid radie R, 45°, höjd z, från L2-projicerade fält."""
    P = np.array([[R / np.sqrt(2)], [R / np.sqrt(2)], [z]])
    sxx, syy, sxy = (report.sample(m, n[k], P)[0] for k in ("sxx", "syy", "sxy"))
    return (sxx + syy) / 2 + sxy, (sxx + syy) / 2 - sxy


def loads_test():
    """static_ext: tryck (Lamé), rotation (roterande rör), egenvikt (hängande stång)."""
    a, b, H, E, nu = 20.0, 40.0, 8.0, 210e3, 0.3
    m = _tube(a, b, H, 2.5)
    fix = {"sx": [0], "sy": [1], "ends": [2]}                                 # plan töjning
    r = static_ext(m, E, nu, fix=fix, pressure={"inner": 10.0})
    st_p = _polar(m, report.nodal(r, E, nu), a + 0.01, H / 2)[1]
    st_p_hand = 10.0 * (b**2 + a**2) / (b**2 - a**2)
    eq_p = abs(r["F"][0] + r["R"][0]) / (10.0 * a * H)
    rho, om = 7.85e-9, 2 * np.pi * 100                                       # 6000 rpm
    rw = static_ext(m, E, nu, fix=fix, rho=rho, omega=om)
    st_w = _polar(m, report.nodal(rw, E, nu), a + 0.01, H / 2)[1]
    v = nu / (1 - nu)                                                        # plan töjning: nu' = nu/(1-nu)
    st_w_hand = rho * om**2 * ((3 + v) / 4 * b**2 + (1 - v) / 4 * a**2)
    F_w_hand = rho * om**2 * (b**3 - a**3) / 3 * H
    L = 1000.0                                                               # hängande stång 10×10
    mb = MeshTet.init_tensor(np.linspace(0, L, 41), np.linspace(0, 10, 3), np.linspace(0, 10, 3))
    mb = mb.with_boundaries({"top": lambda x: np.isclose(x[0], L), "y0": lambda x: np.isclose(x[1], 0),
                             "z0": lambda x: np.isclose(x[2], 0)})
    r = static_ext(mb, E, nu, fix={"top": [0], "y0": [1], "z0": [2]}, rho=rho, gravity=[-9810.0, 0, 0])
    s_top = report.sample(mb, report.nodal(r, E, nu)["sxx"], np.array([[L - 1e-6], [5], [5]]))[0]
    tip = r["u"][r["basis"].split_indices()[0]].min()
    s_g, tip_hand = rho * 9810.0 * L, -rho * 9810.0 * L**2 / (2 * E)
    print(f"Laster: tryck sigma_t,i {st_p:.2f} MPa (Lamé {st_p_hand:.2f}), rotation sigma_t,i {st_w:.3f} MPa "
          f"(hand {st_w_hand:.3f}), F_x {rw['F'][0]:.1f}/{F_w_hand:.1f} N; egenvikt sigma {s_top:.4f} MPa "
          f"(hand {s_g:.4f}), spets {tip * 1e3:.4f} µm (hand {tip_hand * 1e3:.4f})")
    assert abs(st_p / st_p_hand - 1) < 0.02 and eq_p < 1e-6
    assert abs(st_w / st_w_hand - 1) < 0.02 and abs(rw["F"][0] / F_w_hand - 1) < 0.01
    assert abs(s_top / s_g - 1) < 0.01 and abs(tip / tip_hand - 1) < 0.01 and abs(r["F"][0] + r["R"][0]) < 1e-6


def thermal_test():
    """heat (rör med konvektion), static_ext med dT (termisk spänning), heat_transient (platta)."""
    a, b, H, E, nu, al = 20.0, 40.0, 8.0, 210e3, 0.3, 1.2e-5
    k, h, Ti, Tinf = 50.0, 0.05, 120.0, 20.0                                 # 50 W/(m·K), 50 W/(m²·K)
    m = _tube(a, b, H, 2.5)
    rh = heat(m, k, T_fix={"inner": Ti}, conv={"outer": (h, Tinf)})
    R_l, R_c = 4 * np.log(b / a) / (2 * np.pi * k * H), 4 / (h * 2 * np.pi * b * H)
    Q_hand = (Ti - Tinf) / (R_l + R_c)
    To = Tinf + Q_hand * R_c
    Tb = rh["basis"].interpolator(rh["T"])(np.array([[(b - 0.01) / np.sqrt(2)], [(b - 0.01) / np.sqrt(2)], [H / 2]]))[0]
    r = static_ext(m, E, nu, fix={"sx": [0], "sy": [1], "ends": [2]}, alpha=al, dT=rh["T"] - Tinf)
    n = report.nodal(r, E, nu)
    C = E * al * (Ti - To) / (2 * (1 - nu) * np.log(b / a))                  # Timoshenko, plan töjning
    st_a, st_b = _polar(m, n, a + 0.01, H / 2)[1], _polar(m, n, b - 0.01, H / 2)[1]
    st_a_h = C * (1 - 2 * b**2 / (b**2 - a**2) * np.log(b / a))
    st_b_h = C * (1 - 2 * a**2 / (b**2 - a**2) * np.log(b / a))
    L = 100.0                                                                # inspänd stång, fri tvärkontraktion
    mb = MeshTet.init_tensor(np.linspace(0, L, 11), np.linspace(0, 10, 3), np.linspace(0, 10, 3))
    mb = mb.with_boundaries({"ends": lambda x: np.isclose(x[0], 0) | np.isclose(x[0], L),
                             "y0": lambda x: np.isclose(x[1], 0), "z0": lambda x: np.isclose(x[2], 0)})
    rc = static_ext(mb, E, nu, fix={"ends": [0], "y0": [1], "z0": [2]}, alpha=al, dT=50.0)
    s_c = rc["sigma"][0, 0].mean()
    Ls, rho, c = 20.0, 7.85e-9, 460e6                                        # platta, båda sidor till 100 °C
    ms = MeshTet.init_tensor(np.linspace(0, Ls, 21), np.linspace(0, 2, 2), np.linspace(0, 2, 2))
    ms = ms.with_boundaries({"ends": lambda x: np.isclose(x[0], 0) | np.isclose(x[0], Ls)})
    t_end = 0.1 * Ls**2 / (k / (rho * c))                                    # Fourier-tal 0,1
    rt = heat_transient(ms, k, rho, c, 0.0, t_end, 800, T_fix={"ends": 100.0}, save=3)
    Tc = rt["basis"].interpolator(rt["T"][-1])(np.array([[Ls / 2], [1], [1]]))[0]
    nn = np.arange(1, 400, 2)
    Tc_h = 100 * (1 - np.sum(4 / (nn * np.pi) * np.sin(nn * np.pi / 2) * np.exp(-nn**2 * np.pi**2 * 0.1)))
    print(f"Termik: Q {rh['flows']['T_fix']:.1f} mW (hand {Q_hand:.1f}), "
          f"T_ytter {Tb:.2f} (hand {To:.2f}), balans {rh['balance']:.0e}; sigma_t inner {st_a:.2f}/{st_a_h:.2f}, "
          f"ytter {st_b:.2f}/{st_b_h:.2f} MPa; inspänd {s_c:.1f}/{-E * al * 50:.1f} MPa; "
          f"transient T_mitt {Tc:.2f}/{Tc_h:.2f} °C ({rt['DOF']} DOF, {len(rt['T'])} sparade)")
    assert abs(rh["flows"]["T_fix"] / Q_hand - 1) < 0.005 and abs(Tb - To) < 0.05 and rh["balance"] < 1e-8
    assert abs(st_a / st_a_h - 1) < 0.02 and abs(st_b / st_b_h - 1) < 0.02
    assert abs(s_c / (-E * al * 50) - 1) < 1e-6
    assert abs(Tc - Tc_h) < 0.2


def data():
    """materials.json komplett och tjockleksuppslag rätt; projektbladsmallen ska stoppas."""
    import json, pathlib
    for namn, rows in materials().items():
        for r in rows:
            assert {"nu", "rho", "alpha", "typ", "kalla"} <= set(r), namn
    assert material("S355", t=20)["Rp02"] == 345 and material("6082-T6", t=4)["Rp02"] == 250
    assert material("PA12 SLS (EOS PA 2200)")["Rm"] == 42                     # Z-riktning default
    mall = pathlib.Path(__file__).resolve().parent.parent / "mallar" / "projektblad.json"
    if mall.exists():                                                          # finns i repo, ej i pip-installation
        assert len(projektblad.check(json.loads(mall.read_text(encoding="utf-8")))) >= 5
    print(f"Data: {len(materials())} material OK, projektbladskontroll OK")


if __name__ == "__main__":
    t0 = time.perf_counter()
    cantilever()
    plate_hole()
    loads_test()
    thermal_test()
    data()
    print(f"femtools verifierad OK ({time.perf_counter() - t0:.1f} s)")
