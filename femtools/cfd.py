"""Laminär inkompressibel strömning (Stokes / stationär Navier-Stokes), Taylor-Hood P2/P1-tetraedrar.

Enheter mm, t, s, MPa (som resten av femtools):
    μ [Pa·s] × 1e-9 -> MPa·s      (vatten 20 °C: 1.0e-3 Pa·s -> 1.0e-9)
    ρ [kg/m³] × 1e-12 -> t/mm³    (vatten: 998 -> 9.98e-10)
    hastighet mm/s, tryck MPa (1 Pa = 1e-6 MPa), volymflöde mm³/s

    from femtools import flow, boundary_flow
    r = flow(m, mu=1.0e-9, rho=9.98e-10,
             u_fix={"inlopp": lambda x: [0, 0, U(x)], "vagg": [0, 0, 0]})   # utlopp: fri rand (do-nothing, p ≈ 0)
    Q, p_mean = boundary_flow(m, r, "inlopp")

Randvillkor: hastighet på namngivna ränder (konstant vektor eller funktion av x, form (3, n)); ränder utan
villkor blir fria utlopp, där μ∂u/∂n − p·n = 0. Saknas fri rand fixeras trycket i en punkt (p_ref).
rho = 0 ger Stokes (krypströmning). Navier-Stokes löses med Picard-iteration (Oseen), som konvergerar för
måttliga Re (cirka < 1000). Laminärt bara: kontrollera Re mot omslaget (rör cirka 2300). Turbulent: OpenFOAM.
"""
import numpy as np
import scipy.sparse as sp
from skfem import ElementVector, ElementTetP2, ElementTetP1, Basis, FacetBasis, asm, BilinearForm, LinearForm
from skfem.helpers import ddot, grad, div, dot
from .core import SOLVER, SOLVER_NAME


def _u_dofs(ub, u_fix):
    D, vals = [], []
    for name, val in u_fix.items():
        d = ub.get_dofs(name)
        for c in range(3):
            dc = d.all(f"u^{c + 1}")
            if callable(val):
                pts = ub.doflocs[:, dc]
                vc = np.asarray(val(pts), float)[c] * np.ones(len(dc))
            else:
                vc = np.full(len(dc), float(val[c]))
            D.append(dc)
            vals.append(vc)
    D = np.concatenate(D)
    D, i = np.unique(D, return_index=True)
    return D, np.concatenate(vals)[i]


def flow(m, mu, rho=0.0, u_fix=None, p_ref=None, body=None, tol=1e-8, maxit=50, u0=None):
    ue = ElementVector(ElementTetP2())
    ub = Basis(m, ue, intorder=4)
    pb = ub.with_element(ElementTetP1())
    A = asm(BilinearForm(lambda u, v, w: mu * ddot(grad(u), grad(v))), ub)
    B = asm(BilinearForm(lambda u, q, w: -q * div(u)), ub, pb)          # (np × nu)
    nu_, np_ = A.shape[0], B.shape[0]
    f = np.zeros(nu_)
    if body is not None:                                                 # volymkraft [N/mm³]
        b = np.asarray(body, float)
        f += asm(LinearForm(lambda v, w: b[0] * v[0] + b[1] * v[1] + b[2] * v[2]), ub)
    D, uD = _u_dofs(ub, u_fix or {})
    free_bnd = set(m.boundaries) - set(u_fix or {})
    allfacets = np.unique(np.concatenate([m.boundaries[k] for k in (u_fix or {})])) if u_fix else []
    fix_p = p_ref is not None or len(allfacets) == len(m.boundary_facets())
    x = np.zeros(nu_ + np_)
    if u0 is not None:
        x[:nu_] = u0
    x[D] = uD
    Dall = D
    if fix_p:                                                            # trycket bestäms bara upp till en konstant
        pt, pv = (np.zeros(3), 0.0) if p_ref is None else (np.asarray(p_ref[0], float), float(p_ref[1]))
        j = np.argmin(np.linalg.norm(m.p.T - pt, axis=1))
        Dall = np.concatenate([D, [nu_ + j]])
        x[nu_ + j] = pv
    I = np.setdiff1d(np.arange(nu_ + np_), Dall)
    it, change = 0, np.inf
    for it in range(1, (maxit if rho else 1) + 1):
        Ak = A
        if rho:
            wv = ub.interpolate(x[:nu_])
            Ak = A + asm(BilinearForm(lambda u, v, w: rho * np.einsum("jeq,ijeq,ieq->eq", w.wv, grad(u), v)),
                         ub, wv=wv)                                      # ρ (w·∇)u · v, w = förra iterationen
        K = sp.bmat([[Ak, B.T], [B, None]], format="csr")
        rhs = np.concatenate([f, np.zeros(np_)]) - K @ x
        xn = x.copy()
        xn[I] = x[I] + SOLVER(K[I][:, I], rhs[I])
        change = np.linalg.norm(xn[:nu_] - x[:nu_]) / max(np.linalg.norm(xn[:nu_]), 1e-300)
        x = xn
        if change < tol:
            break
    if rho and change >= tol:
        raise RuntimeError(f"Navier-Stokes konvergerade inte på {maxit} Picard-iterationer (ändring {change:.1e}).")
    u, p = x[:nu_], x[nu_:]
    return dict(ubasis=ub, pbasis=pb, u=u, p=p, iterations=it, change=change, DOF=nu_ + np_,
                solver=SOLVER_NAME, free_boundaries=sorted(free_bnd), mu=mu, rho=rho)


def boundary_flow(m, r, name):
    """Volymflöde ut genom randen [mm³/s] (negativt = in) och areaviktat medeltryck [MPa]."""
    fb = FacetBasis(m, ElementVector(ElementTetP2()), facets=m.boundaries[name], intorder=4)
    fp = fb.with_element(ElementTetP1())
    from skfem import Functional
    Q = asm(Functional(lambda w: dot(w.u, w.n)), fb, u=fb.interpolate(r["u"]))
    Ar = asm(Functional(lambda w: 1.0 + 0 * w.x[0]), fb)
    pm = asm(Functional(lambda w: w.p), fp, p=fp.interpolate(r["p"])) / Ar
    return Q, pm


def reynolds(rho, U, L, mu):
    """Re = ρUL/μ (enheter enligt ovan, ger dimensionslöst tal)."""
    return rho * U * L / mu
