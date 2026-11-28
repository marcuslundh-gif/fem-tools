"""Verifiering mot handräkning. Kör: python -m femtools.verify"""
import time
import numpy as np
from skfem import MeshTet
import gmsh
from .core import static, modal
from .loads import static_ext
from .thermal import heat, heat_transient, heat_transient_coupled, SIGMA
from .buckling import buckling
from .cfd import flow, boundary_flow
from . import openfoam
from .core import richardson
from .mesh import setup, refine_and_mesh, to_skfem
from .materials import material, materials
from . import projektblad, report, webreport, shell


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


def thermal_coupled_test():
    """heat_transient_coupled: klumpad nod + konvektion (2 noder, exakt med expm) och strålning +
    temperaturberoende konvektion (ODE med solve_ivp). Kuben har mycket hög k, så den är isoterm."""
    from scipy.linalg import expm
    from scipy.integrate import solve_ivp
    t0 = time.perf_counter()
    a = 10.0
    mc = MeshTet.init_tensor(*[np.linspace(0, a, 4)] * 3)
    mc = mc.with_boundaries({"x0": lambda x: np.isclose(x[0], 0), "x1": lambda x: np.isclose(x[0], a),
                             "alla": lambda x: np.isclose(x[0], 0) | np.isclose(x[0], a) | np.isclose(x[1], 0)
                             | np.isclose(x[1], a) | np.isclose(x[2], 0) | np.isclose(x[2], a)})
    rho, c, A = 2.4e-9, 830e6, a * a
    Cs, CL, hc, h0, G = rho * c * a**3, 4000.0, 0.2, 0.01, 2.0
    r = heat_transient_coupled(mc, 1e5, rho, c, 20.0, 2000.0, 2000, save=2,
                               conv={"x0": [(hc, "L")], "x1": [(h0, 0.0)]},
                               lumped=dict(C=CL, T0=80.0, loss=lambda TL, t: (G * TL, G)))
    Amat = np.array([[-(hc * A + h0 * A) / Cs, hc * A / Cs], [hc * A / CL, -(hc * A + G) / CL]])
    xs, xl = expm(Amat * 2000.0) @ np.array([20.0, 80.0])
    Ts = r["T"][-1].mean()
    # strålning (eps 0,9) + konvektion h = 1e-3·(Ts − 20)^0,25 från 300 °C mot 20 °C
    rho2, c2, eps = 1e-9, 1e8, 0.9
    C2, A2 = rho2 * c2 * a**3, 6 * a * a
    hfun = lambda Ts: 1e-3 * max(Ts - 20.0, 0.0) ** 0.25
    r2 = heat_transient_coupled(mc, 1e5, rho2, c2, 300.0, 30.0, 300, save=2,
                                conv={"alla": lambda s: [(hfun(s["Ts"]), 20.0)]}, rad={"alla": (eps, 20.0)})
    ode = lambda t, T: [-A2 / C2 * (hfun(T[0]) * (T[0] - 20) + eps * SIGMA * ((T[0] + 273.15)**4 - 293.15**4))]
    Tref2 = solve_ivp(ode, (0, 30.0), [300.0], rtol=1e-10, atol=1e-10).y[0, -1]
    T2 = r2["T"][-1].mean()
    print(f"Termik kopplad: kropp {Ts:.3f}/{xs:.3f} °C, nod {r['TL'][-1]:.3f}/{xl:.3f} °C, balans {r['balance']:.1e}; "
          f"strålning+konv {T2:.3f}/{Tref2:.3f} °C, balans {r2['balance']:.1e}, {r2['refactor']} faktoriseringar, "
          f"{time.perf_counter() - t0:.1f} s")
    assert abs(Ts - xs) < 0.05 and abs(r["TL"][-1] - xl) < 0.05 and r["balance"] < 1e-8
    assert abs(T2 - Tref2) < 0.3 and r2["balance"] < 1e-3


def buckling_test():
    """buckling: Euler-pelare 500×10×10 mm, fall 1 (inspänd-fri) och fall 3 (inspänd-ledad)."""
    E, nu, L, a = 210e3, 0.3, 500.0, 10.0
    I = a**4 / 12
    m = MeshTet.init_tensor(np.linspace(0, L, 51), np.linspace(0, a, 3), np.linspace(0, a, 3))
    m = m.with_boundaries({"bas": lambda x: np.isclose(x[0], 0), "topp": lambda x: np.isclose(x[0], L)})
    t0 = time.perf_counter()
    r1 = buckling(m, E, nu, k=3, fix="bas", traction={"topp": [-1.0 / a**2, 0, 0]})        # 1 N tryck
    r3 = buckling(m, E, nu, k=1, fix={"bas": "all", "topp": [1, 2]}, traction={"topp": [-1.0 / a**2, 0, 0]})
    P1, P3 = np.pi**2 * E * I / (4 * L**2), 20.19 * E * I / L**2                          # 4,4934² = 20,19
    f = r1["factors"]
    print(f"Knäckning: fall 1 {f[0]:.0f} N (Euler {P1:.0f}), mod 3/mod 1 {f[2] / f[0]:.2f} (teori 9,00), "
          f"fall 3 {r3['factors'][0]:.0f} N (Euler {P3:.0f}), {r1['DOF']} DOF, {time.perf_counter() - t0:.1f} s")
    assert abs(f[0] / P1 - 1) < 0.005 and abs(f[0] / f[1] - 1) < 0.001               # två lika moder (kvadrat)
    assert abs(f[2] / f[0] / 9 - 1) < 0.02 and abs(r3["factors"][0] / P3 - 1) < 0.01


def cfd_test():
    """flow: Kovasznay (exakt Navier-Stokes, Re 40) och Hagen-Poiseuille (Stokes, Δp med Richardson)."""
    t0 = time.perf_counter()
    Re = 40.0
    lam = Re / 2 - np.sqrt(Re**2 / 4 + 4 * np.pi**2)
    ue = lambda x: np.array([1 - np.exp(lam * x[0]) * np.cos(2 * np.pi * x[1]),
                             lam / (2 * np.pi) * np.exp(lam * x[0]) * np.sin(2 * np.pi * x[1]), 0 * x[0]])
    pe = lambda x: (1 - np.exp(2 * lam * x[0])) / 2
    m = MeshTet.init_tensor(np.linspace(-0.5, 1, 17), np.linspace(-0.5, 1.5, 22), np.linspace(0, 0.1, 2))
    m = m.with_boundaries({"alla": lambda x: np.ones(x.shape[1], bool)})
    r = flow(m, mu=1 / Re, rho=1.0, u_fix={"alla": ue}, p_ref=([1.0, 0.5, 0.0], float(pe(np.array([1.0])))))
    ub, xq = r["ubasis"], r["ubasis"].mapping.F(r["ubasis"].X)
    uq, pq, dx = ub.interpolate(r["u"]).value, r["pbasis"].interpolate(r["p"]).value, ub.dx
    eu = np.sqrt((((uq - ue(xq))**2).sum(0) * dx).sum() / ((ue(xq)**2).sum(0) * dx).sum())
    ep = np.sqrt((((pq - pe(xq))**2) * dx).sum() / ((pe(xq)**2) * dx).sum())
    D, L, mu, U = 10.0, 30.0, 1.0e-9, 50.0                                   # vatten, mm-t-s-MPa
    R = D / 2
    prof = lambda x: np.array([0 * x[0], 0 * x[0], 2 * U * np.clip(1 - (x[0]**2 + x[1]**2) / R**2, 0, None)])
    dps = []
    for nc in (16, 24, 36):                                                  # element per varv, r = 1,5
        setup()
        gmsh.model.occ.addCylinder(0, 0, 0, 0, 0, L, R)
        gmsh.model.occ.synchronize()
        refine_and_mesh(2 * np.pi * R / nc, n_curv=nc)
        mp = to_skfem()
        gmsh.finalize()
        mp = mp.with_boundaries({"in": lambda x: np.isclose(x[2], 0), "ut": lambda x: np.isclose(x[2], L),
                                 "vagg": lambda x: ~(np.isclose(x[2], 0) | np.isclose(x[2], L))})
        rp = flow(mp, mu, 0.0, u_fix={"in": prof, "vagg": [0, 0, 0]})
        (Qi, pi), (Qo, po) = boundary_flow(mp, rp, "in"), boundary_flow(mp, rp, "ut")
        dps.append(pi - po)
    p_obs, dp_ext = richardson(*dps, 1.5)
    dp_hand = 32 * mu * L * U / D**2
    print(f"CFD laminärt: Kovasznay Re 40 L2-fel u {eu:.2%}, p {ep:.2%} ({r['iterations']} Picard, {r['DOF']} DOF); "
          f"Poiseuille dp {dps[-1] * 1e6:.4f} Pa, Richardson {dp_ext * 1e6:.4f} (p {p_obs:.1f}) mot {dp_hand * 1e6:.4f} Pa, "
          f"Q ut {Qo:.0f}/{np.pi * R**2 * U:.0f} mm3/s, {time.perf_counter() - t0:.1f} s")
    assert eu < 0.002 and ep < 0.01 and abs(dp_ext / dp_hand - 1) < 0.005 and abs(dps[-1] / dp_hand - 1) < 0.02
    assert abs(Qo / (np.pi * R**2 * U) - 1) < 0.002


def openfoam_test():
    """OpenFOAM k-ω SST: utvecklad turbulent rörströmning mot Colebrook (släta rör). Hoppas över utan OpenFOAM."""
    ver = openfoam.available()
    if not ver:
        print("OpenFOAM: saknas i miljön, turbulent CFD är INTE verifierad här (apt-get install -y openfoam)")
        return
    import tempfile, os
    t0 = time.perf_counter()
    res = []
    with tempfile.TemporaryDirectory() as td:
        for Re, n, w in ((1e5, 160, 0.005), (1e6, 200, 0.0005)):
            r = openfoam.turbulent_pipe(os.path.join(td, f"ror{int(Re)}"), Re=Re, n_r=n, wall_ratio=w, iters=6000)
            res.append(r)
    txt = ", ".join(f"Re {r['Re']:.0e}: f {r['f']:.5f} mot {r['f_colebrook']:.5f} ({r['f'] / r['f_colebrook'] - 1:+.1%}, "
                    f"y+ {r['yplus']:.2f})" for r in res)
    print(f"OpenFOAM {ver} k-omega SST rör: {txt}, {time.perf_counter() - t0:.0f} s")
    for r in res:
        assert abs(r["f"] / r["f_colebrook"] - 1) < 0.03 and r["yplus"] < 1


def webreport_test():
    """webreport: ytan utåtorienterad, fält och tabeller hamnar i en fristående HTML-sida."""
    import tempfile, os, json
    m = MeshTet.init_tensor(np.linspace(0, 2, 5), np.linspace(0, 1, 3), np.linspace(0, 1, 4))
    w = webreport.Report("Test", meta={"Fall": "verify"})
    w.mesh(m)
    p, t = w._surf["p"], w._surf["tri"]
    nrm = np.cross(p[t[:, 1]] - p[t[:, 0]], p[t[:, 2]] - p[t[:, 0]])
    out = (np.einsum("ij,ij->i", nrm, p[t].mean(1) - [1, 0.5, 0.5]) > 0).mean()
    vals = m.p[0] * 10
    w.field("LF1", "x", vals, "MPa", u=np.vstack([m.p[0] * 0.01, 0 * m.p[0], 0 * m.p[0]]), limit=15)
    tt = np.array([0.0, 60.0, 120.0])
    w.transient("Transient", "T", tt, [vals * (1 - 0.25 * k) for k in range(3)], "°C",
                curves={"medel": (tt, [vals.mean() * (1 - 0.25 * k) for k in range(3)])})
    w.table("T", ["a", "b"], [["x", 1.5]])
    w.text("Text", "Stycke 1\n\nStycke 2")
    with tempfile.TemporaryDirectory() as td:
        html = open(w.write(os.path.join(td, "r.html")), encoding="utf-8").read()
    d = w.data()
    f = d["cases"][0]["fields"][0]
    print(f"Webbrapport: {len(t)} yttrianglar, {out:.0%} utåt, max {f['max']:.1f} vid {f['at']}, {len(html) // 1024} kB")
    assert out == 1.0 and abs(f["max"] - 20) < 1e-9 and f["at"][0] == 2.0
    assert "<title>Test</title>" in html and "three@0.128.0" in html and json.dumps(d["title"]) in html
    ft = d["cases"][1]["fields"][0]
    assert len(ft["frames"]) == 3 and abs(ft["max"] - 20) < 1e-9 and abs(min(ft["frames"][2]) - 0) < 1e-9
    assert ft["t"] == [0.0, 1.0, 2.0] and abs(ft["curves"][0]["y"][2] - vals.mean() * 0.5) < 1e-3 and 'id="tsl"' in html


def webreport_flow_test():
    """webreport.flow + openfoam.read_case: syntetiskt OpenFOAM-fall (box 2×1×1 tetraedrar, U = (1, 0, 0), p = x).
    Snittet y = 0,5 ska ha arean 2,0 och p/|U| exakt interpolerat; strömlinjen ska gå rakt från x 0 till 2."""
    import tempfile, os
    m = MeshTet.init_tensor(np.linspace(0, 2, 9), np.linspace(0, 1, 5), np.linspace(0, 1, 5))
    P, tet = m.p.T, m.t.T
    fc = np.sort(np.vstack([tet[:, [0, 1, 2]], tet[:, [0, 1, 3]], tet[:, [0, 2, 3]], tet[:, [1, 2, 3]]]), axis=1)
    cell = np.tile(np.arange(len(tet)), 4)
    uq, inv, cnt = np.unique(fc, axis=0, return_inverse=True, return_counts=True)
    inv = inv.ravel()
    own, nei = {}, {}
    for f, c in zip(inv, cell):
        (nei if f in own else own)[f] = c if f not in own else c
    order = [f for f in range(len(uq)) if cnt[f] == 2] + [f for f in range(len(uq)) if cnt[f] == 1]
    nint = int((cnt == 2).sum())
    lst = lambda xs: f"{len(xs)}\n(\n" + "\n".join(xs) + "\n)\n"
    with tempfile.TemporaryDirectory() as td:
        pm = os.path.join(td, "constant", "polyMesh"); os.makedirs(pm); os.makedirs(os.path.join(td, "1"))
        w = lambda p, t: open(os.path.join(td, p), "w").write("FoamFile{}\n" + t)
        w("constant/polyMesh/points", lst([f"({x / 1000} {y / 1000} {z / 1000})" for x, y, z in P]))
        w("constant/polyMesh/faces", lst(["3(%d %d %d)" % tuple(uq[f]) for f in order]))
        w("constant/polyMesh/owner", lst([str(own[f]) for f in order]))
        w("constant/polyMesh/neighbour", lst([str(nei[f]) for f in order[:nint]]))
        w("constant/polyMesh/boundary", f"1\n(\nwalls\n{{\n type wall;\n nFaces {len(uq) - nint};\n startFace {nint};\n}}\n)\n")
        cen = P[tet].mean(1)
        w("1/U", f"internalField nonuniform List<vector> {len(tet)}\n(\n" + "\n".join("(1 0 0)" for _ in tet) + "\n)\n;")
        w("1/p", f"internalField nonuniform List<scalar> {len(tet)}\n(\n" + "\n".join(f"{x:.9g}" for x in cen[:, 0]) + "\n)\n;")
        f = openfoam.read_case(td, rho=2.0, wall="ingen")
    vol = np.abs(np.linalg.det(f["P"][f["tet"]][:, 1:] - f["P"][f["tet"]][:, :1])).sum() / 6
    assert np.allclose(f["P"], P) and abs(vol - 2.0) < 1e-9 and np.allclose(f["pc"], 2 * cen[:, 0])
    pts, tri, val = webreport.slice_tets(P, tet, P[:, 0], "y", 0.5)
    area = 0.5 * np.linalg.norm(np.cross(pts[tri[:, 1]] - pts[tri[:, 0]], pts[tri[:, 2]] - pts[tri[:, 0]]), axis=1).sum()
    L = webreport.streamlines(P, tet, f["U"], [[0.01, 0.37, 0.61]], ds=0.05)[0][0]
    w_ = webreport.Report("Flöde", front=(0, 1, 0)).mesh(m)
    w_.flow("test", f, slices={"y": 0.5, "z": 0.25}, seeds=[[0.01, 0.37, 0.61]], wall="walls")
    d = w_.data()["flows"][0]
    print(f"Flöde i webbrapport: snittarea {area:.6f} (2), x-fel {np.abs(val - pts[:, 0]).max():.1e}, "
          f"strömlinje x {L[0, 0]:.2f}→{L[-1, 0]:.3f}, avdrift {np.abs(L[L[:, 0] < 1.9, 1:] - [0.37, 0.61]).max():.1e}")
    assert abs(area - 2.0) < 1e-9 and np.abs(val - pts[:, 0]).max() < 1e-12
    inner = L[:, 0] < 1.9                     # vid utloppet knuffas linjen (randbeteende, se streamlines)
    assert L[-1, 0] > 1.9 and np.abs(L[inner, 1:] - [0.37, 0.61]).max() < 1e-9
    assert len(d["slices"]) == 2 and d["wall"] is not None and len(d["stream"]) == 1 and w_.front == [0.0, 1.0, 0.0]


def data():
    """materials.json komplett och tjockleksuppslag rätt; projektbladsmallen ska stoppas."""
    import json, pathlib
    for namn, rows in materials().items():
        for r in rows:
            assert {"nu", "rho", "alpha", "typ", "kalla"} <= set(r), namn
    assert material("S355", t=20)["Rp02"] == 345 and material("6082-T6", t=4)["Rp02"] == 250
    assert material("1.4301 plattstång varmvalsad", t=10)["Rp02"] == 190 and material("42CrMo4+QT", t=60)["Rp02"] == 650
    assert material("PA12 SLS (EOS PA 2200)")["Rm"] == 42                     # Z-riktning default
    assert material("AA3104-H19 burkkropp", t=0.095)["Rp02"] == 293 and material("AA5182-H48 burklock", t=0.208)["E"] == 70000
    mall = pathlib.Path(__file__).resolve().parent.parent / "mallar" / "projektblad.json"
    if mall.exists():                                                          # finns i repo, ej i pip-installation
        saknas = projektblad.check(json.loads(mall.read_text(encoding="utf-8")))
        assert len(saknas) >= 5 and any(k.startswith("betraktarvy.fram") for k in saknas)
    print(f"Data: {len(materials())} material OK, projektbladskontroll OK")


def shell_test():
    """shell (CalculiX S8R, olinjärt): cylinder R 29, t 0,1, L 20 mm med Koiters axisymmetriska imperfektion 0,1t.
    Hoppas över utan CalculiX."""
    ver = shell.available()
    if not ver:
        print("Skal: CalculiX saknas i miljön, femtools.shell är INTE verifierad här (apt-get install -y calculix-ccx)")
        return
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        r = shell.cylinder_collapse(td, xi=0.1)
    print(f"Skal ccx {ver}: axiell styvhet {r['k']:.0f} N/mm (EA/L {r['k_hand']:.0f}, {r['k'] / r['k_hand'] - 1:+.1%}), "
          f"kollaps {r['F_max']:.0f} N = {r['ratio']:.3f}·P_cl (Koiter {r['koiter']:.3f}, {r['ratio'] / r['koiter'] - 1:+.1%}), "
          f"{r['n_nodes']} noder, {r['time']:.0f} s")
    assert not r["peak_last"]                     # topp före sista inkrement (ccx avbryter efter kollaps)
    assert abs(r["k"] / r["k_hand"] - 1) < 0.02 and abs(r["ratio"] / r["koiter"] - 1) < 0.06


if __name__ == "__main__":
    t0 = time.perf_counter()
    cantilever()
    plate_hole()
    loads_test()
    thermal_test()
    thermal_coupled_test()
    buckling_test()
    cfd_test()
    openfoam_test()
    shell_test()
    webreport_test()
    webreport_flow_test()
    data()
    print(f"femtools verifierad OK ({time.perf_counter() - t0:.1f} s)")
