"""Linjär statik med fler lasttyper (mm, N, MPa, t, s), P2-tetraedrar.

    from femtools import static_ext
    r = static_ext(m, E, nu,
                   fix={"sym_x": [0], "sym_y": [1], "bas": "all"},   # randnamn -> komponenter (0=x, 1=y, 2=z)
                   traction={"lock": [0, 0, -2.0]},                   # MPa, global riktning
                   pressure={"inner": 10.0},                          # MPa, verkar mot ytan (−p·n)
                   rho=7.85e-9, gravity=[0, 0, -9810.0],              # t/mm³, mm/s²
                   omega=314.0, axis=([0, 0, 0], [0, 0, 1]),          # rad/s, rotationsaxel (punkt, riktning)
                   alpha=1.2e-5, dT=dT)                               # 1/K; K, konstant eller P2-fält från thermal.heat

Returnerar samma nycklar som `static` (basis, u, vm, sigma, xq, R, F, DOF, solver). sigma är total spänning
C:(ε − αΔT·I), dvs. termisk töjning är borttagen. F är summan av alla pålagda laster (yta + volym).
`static` lämnas orörd; den här funktionen är verifierad separat i verify.loads_test().
"""
import numpy as np
from skfem import ElementVector, ElementTetP2, Basis, FacetBasis, asm, condense, solve, LinearForm
from skfem.helpers import div
from skfem.models.elasticity import linear_elasticity, lame_parameters
from .core import SOLVER, SOLVER_NAME


def _dofs(ib, fix, D):
    """Fixerade DOF från randnamn. fix: str, lista av str (alla komponenter) eller dict namn -> 'all'/[komp]."""
    idx = []
    if fix is not None:
        if isinstance(fix, str):
            fix = {fix: "all"}
        elif not isinstance(fix, dict):
            fix = {f: "all" for f in fix}
        for name, comps in fix.items():
            d = ib.get_dofs(name)
            comps = [0, 1, 2] if comps == "all" else comps
            idx += [d.all(f"u^{c + 1}") for c in comps]
    if D is not None:
        idx.append(np.asarray(D, dtype=np.int64))
    if not idx:
        raise ValueError("Inga randvillkor: ange fix eller D.")
    return np.unique(np.concatenate(idx))


def _dT_field(basis, dT):
    """ΔT i integrationspunkterna för `basis` (konstant eller P2-DOF-vektor på samma nät)."""
    if dT is None:
        return np.zeros((basis.nelems, basis.X.shape[1]))
    if np.isscalar(dT):
        return np.full((basis.nelems, basis.X.shape[1]), float(dT))
    return basis.with_element(ElementTetP2()).interpolate(np.asarray(dT, float)).value


def stress(basis, u, E, nu, alpha=0.0, dT=None):
    """Spänning C:(ε − αΔT·I) och von Mises i integrationspunkter. Returnerar (vm, sigma)."""
    lam, mu = lame_parameters(E, nu)
    G = basis.interpolate(u).grad
    eps = 0.5 * (G + G.transpose(1, 0, 2, 3))
    th = (3 * lam + 2 * mu) * alpha * _dT_field(basis, dT)
    s = (2 * mu * eps + (lam * (eps[0, 0] + eps[1, 1] + eps[2, 2]) - th) * np.eye(3)[:, :, None, None])
    vm = np.sqrt(0.5 * ((s[0, 0] - s[1, 1])**2 + (s[1, 1] - s[2, 2])**2 + (s[2, 2] - s[0, 0])**2)
                 + 3 * (s[0, 1]**2 + s[1, 2]**2 + s[0, 2]**2))
    return vm, s


def static_ext(m, E, nu, fix="fix", D=None, traction=None, pressure=None, rho=None, gravity=None,
               omega=None, axis=((0, 0, 0), (0, 0, 1)), alpha=None, dT=None):
    e = ElementVector(ElementTetP2())
    ib = Basis(m, e, intorder=2)
    lam, mu = lame_parameters(E, nu)
    K = asm(linear_elasticity(lam, mu), ib)
    f = np.zeros(K.shape[0])

    for name, t in (traction or {}).items():
        t = np.asarray(t, float)
        fb = FacetBasis(m, e, facets=m.boundaries[name], intorder=4)
        f += asm(LinearForm(lambda v, w: t[0] * v[0] + t[1] * v[1] + t[2] * v[2]), fb)

    for name, p in (pressure or {}).items():
        fb = FacetBasis(m, e, facets=m.boundaries[name], intorder=4)
        f += asm(LinearForm(lambda v, w: -p * (w.n[0] * v[0] + w.n[1] * v[1] + w.n[2] * v[2])), fb)

    vb = Basis(m, e, intorder=4)
    if gravity is not None or omega is not None:
        if rho is None:
            raise ValueError("rho [t/mm³] krävs för gravitation/rotation.")
        g = np.zeros(3) if gravity is None else np.asarray(gravity, float)
        a0, ax = np.asarray(axis[0], float), np.asarray(axis[1], float)
        ax = ax / np.linalg.norm(ax)
        w2 = 0.0 if omega is None else float(omega)**2

        @LinearForm
        def body(v, w):
            x = [w.x[i] - a0[i] for i in range(3)]
            xa = x[0] * ax[0] + x[1] * ax[1] + x[2] * ax[2]
            b = [rho * (g[i] + w2 * (x[i] - xa * ax[i])) for i in range(3)]     # N/mm³
            return b[0] * v[0] + b[1] * v[1] + b[2] * v[2]

        f += asm(body, vb)

    a = 0.0 if alpha is None else float(alpha)
    if dT is not None and a != 0.0:
        th = (3 * lam + 2 * mu) * a

        @LinearForm
        def thermal(v, w):
            return th * w.dT * div(v)

        f += asm(thermal, vb, dT=_dT_field(vb, dT))

    Dd = _dofs(ib, fix, D)
    u = solve(*condense(K, f, D=Dd), solver=SOLVER)
    r = K @ u - f
    R, F = [], []
    for ci in ib.split_indices():                      # alla DOF per komponent (även inre, för volymlaster)
        R.append(r[ci[np.isin(ci, Dd)]].sum())
        F.append(f[ci].sum())
    vm, s = stress(vb, u, E, nu, a, dT)
    return dict(basis=vb, u=u, vm=vm, sigma=s, xq=vb.mapping.F(vb.X), R=R, F=F, DOF=K.shape[0],
                solver=SOLVER_NAME)
