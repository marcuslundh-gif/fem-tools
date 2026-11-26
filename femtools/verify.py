"""Verifiering mot handräkning. Kör: python -m femtools.verify"""
import time
import numpy as np
from skfem import MeshTet
import gmsh
from .core import static, modal
from .mesh import setup, refine_and_mesh, to_skfem
from .materials import material, materials
from . import projektblad


def cantilever():
    L, B, H, E, nu, F = 100.0, 10.0, 10.0, 210e3, 0.3, -100.0
    m = MeshTet.init_tensor(np.linspace(0, L, 41), np.linspace(0, B, 5), np.linspace(0, H, 5))
    m = m.with_boundaries({"fix": lambda x: np.isclose(x[0], 0), "load": lambda x: np.isclose(x[0], L)})
    r = static(m, E, nu, [0, 0, F / (B * H)])
    uz = r["u"][r["basis"].get_dofs("load").all("u^3")].mean()
    I = B * H**3 / 12
    d_hand = F * L**3 / (3 * E * I)
    f1 = modal(m, E, nu, 7.85e-9, k=2)[0]
    f_hand = 1.8751**2 / (2 * np.pi) * np.sqrt(E * I / (7.85e-9 * B * H * L**4))
    eq = max(abs(a + b) for a, b in zip(r["R"], r["F"])) / abs(F)
    print(f"Konsol: delta {uz:.4f} mm (hand {d_hand:.4f}), f1 {f1:.0f} Hz (hand {f_hand:.0f}), "
          f"jämvikt {eq:.1e}, DOF {r['DOF']}, lösare {r['solver']}")
    assert abs(uz / d_hand - 1) < 0.01 and abs(f1 / f_hand - 1) < 0.01 and eq < 1e-3


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
    print(f"Plåt med hål: sigma_xx,max {sxx:.1f} MPa (Heywood {hand:.1f}), gamma_min {q['gamma_min']:.2f}, "
          f"mesh {t_mesh:.2f} s, DOF {r['DOF']}")
    assert abs(sxx / hand - 1) < 0.06


def data():
    """materials.json komplett och tjockleksuppslag rätt; projektbladsmallen ska stoppas."""
    import json, pathlib
    for namn, rows in materials().items():
        for r in rows:
            assert {"nu", "rho", "alpha", "typ", "kalla"} <= set(r), namn
    assert material("S355", t=20)["Rp02"] == 345 and material("6082-T6", t=4)["Rp02"] == 250
    assert material("1.4301 plattstång varmvalsad", t=10)["Rp02"] == 190 and material("42CrMo4+QT", t=60)["Rp02"] == 650
    assert material("PA12 SLS (EOS PA 2200)")["Rm"] == 42                     # Z-riktning default
    mall = pathlib.Path(__file__).resolve().parent.parent / "mallar" / "projektblad.json"
    if mall.exists():                                                          # finns i repo, ej i pip-installation
        assert len(projektblad.check(json.loads(mall.read_text(encoding="utf-8")))) >= 5
    print(f"Data: {len(materials())} material OK, projektbladskontroll OK")


if __name__ == "__main__":
    t0 = time.perf_counter()
    cantilever()
    plate_hole()
    data()
    print(f"femtools verifierad OK ({time.perf_counter() - t0:.1f} s)")
