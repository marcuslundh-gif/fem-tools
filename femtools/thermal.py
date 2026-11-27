"""Värmeledning, stationär och transient, P2-tetraedrar.

Enheter i mm-N-s-systemet (K för temperatur). Effekt i mW = N·mm/s. Omräkning från SI:
    k   [W/(m·K)]   -> samma siffra       (mW/(mm·K))
    h   [W/(m²·K)]  -> × 1e-3             (mW/(mm²·K))
    q   [W/m²]      -> × 1e-3             (mW/mm²)
    Q   [W/m³]      -> × 1e-6             (mW/mm³)
    c   [J/(kg·K)]  -> × 1e6              (N·mm/(t·K)),  rho i t/mm³ som i resten av femtools
Kontroll: stål k 50, rho 7.85e-9, c 460e6 ger temperaturledningstal 13,8 mm²/s (1,38e-5 m²/s).

    from femtools import heat, heat_transient
    r = heat(m, k=50, T_fix={"inner": 120}, conv={"outer": (10e-3, 20)})
    r["T"]          # P2-DOF-vektor, kan ges direkt som dT till static_ext (dT = T - T_ref)
    r["balance"]    # relativ energibalans, ska vara < 1e-6

Randvillkor per randnamn: T_fix (K), flux (mW/mm², positivt = in i kroppen), conv (h, T_inf).
Värden i T_fix får vara funktioner av tiden t i heat_transient.
"""
import numpy as np
from skfem import ElementTetP2, Basis, FacetBasis, asm, condense, solve, BilinearForm, LinearForm
from skfem.helpers import dot, grad
from .core import SOLVER, SOLVER_NAME


def _system(m, k, flux, conv, Q):
    ib = Basis(m, ElementTetP2(), intorder=4)
    K = asm(BilinearForm(lambda u, v, w: k * dot(grad(u), grad(v))), ib)
    f = np.zeros(K.shape[0])
    if Q:
        f += asm(LinearForm(lambda v, w: Q * v), ib)
    for name, q in (flux or {}).items():
        fb = FacetBasis(m, ElementTetP2(), facets=m.boundaries[name], intorder=4)
        f += asm(LinearForm(lambda v, w: q * v), fb)
    for name, (h, Tinf) in (conv or {}).items():
        fb = FacetBasis(m, ElementTetP2(), facets=m.boundaries[name], intorder=4)
        K = K + asm(BilinearForm(lambda u, v, w: h * u * v), fb)
        f += asm(LinearForm(lambda v, w: h * Tinf * v), fb)
    return ib, K.tocsr(), f


def _fixed(ib, T_fix, t=0.0):
    D, vals = [], []
    for name, T in (T_fix or {}).items():
        d = ib.get_dofs(name).all()
        D.append(d)
        vals.append(np.full(len(d), T(t) if callable(T) else float(T)))
    if not D:
        return np.array([], dtype=np.int64), np.array([])
    D, i = np.unique(np.concatenate(D), return_index=True)
    return D, np.concatenate(vals)[i]


def _flows(m, ib, T, K, f, D, flux, conv, Q):
    """Värmeflöden [mW] in i kroppen per randvillkor och volym; summan ska bli noll."""
    out = {}
    if len(D):
        out["T_fix"] = (K @ T - f)[D].sum()               # reaktion i fixerade noder = värme in
    ones = lambda fb: asm(LinearForm(lambda v, w: v), fb)
    for name, q in (flux or {}).items():
        fb = FacetBasis(m, ElementTetP2(), facets=m.boundaries[name], intorder=4)
        out[f"flux:{name}"] = q * ones(fb).sum()
    for name, (h, Tinf) in (conv or {}).items():
        fb = FacetBasis(m, ElementTetP2(), facets=m.boundaries[name], intorder=4)
        out[f"conv:{name}"] = asm(LinearForm(lambda v, w: h * (Tinf - w.T) * v), fb, T=fb.interpolate(T)).sum()
    if Q:
        out["Q"] = Q * ones(ib).sum()
    return out


def heat(m, k, T_fix=None, flux=None, conv=None, Q=0.0):
    """Stationär värmeledning. k [mW/(mm·K)], Q [mW/mm³]. Returnerar dict med basis, T, flows, balance, DOF."""
    ib, K, f = _system(m, k, flux, conv, Q)
    D, TD = _fixed(ib, T_fix)
    if not len(D) and not conv:
        raise ValueError("Temperaturen är obestämd: ange T_fix eller conv.")
    T = np.zeros(K.shape[0])
    T[D] = TD
    T = solve(*condense(K, f, x=T, D=D), solver=SOLVER) if len(D) else solve(K, f, solver=SOLVER)
    flows = _flows(m, ib, T, K, f, D, flux, conv, Q)
    scale = max(abs(v) for v in flows.values()) or 1.0
    return dict(basis=ib, T=T, flows=flows, balance=abs(sum(flows.values())) / scale,
                Tmin=T.min(), Tmax=T.max(), DOF=K.shape[0], solver=SOLVER_NAME)


def heat_transient(m, k, rho, c, T0, t_end, nsteps, T_fix=None, flux=None, conv=None, Q=0.0,
                   theta=1.0, save=None):
    """Transient värmeledning, theta-metod (1 = implicit Euler, robust; 0.5 = Crank-Nicolson).
    rho [t/mm³], c [N·mm/(t·K)] (= J/(kg·K) × 1e6), T0 konstant eller P2-vektor.
    save: antal jämnt fördelade tidpunkter att spara (default alla steg).
    Returnerar dict med basis, t (sparade tider), T (lista med P2-vektorer), DOF."""
    import pypardiso
    ib, K, f = _system(m, k, flux, conv, Q)
    M = asm(BilinearForm(lambda u, v, w: rho * c * u * v), ib).tocsr()
    dt = t_end / nsteps
    A = (M / dt + theta * K).tocsr()
    B = (M / dt - (1 - theta) * K).tocsr()
    T = np.full(K.shape[0], float(T0)) if np.isscalar(T0) else np.asarray(T0, float).copy()
    D, _ = _fixed(ib, T_fix, 0.0)
    I = np.setdiff1d(np.arange(K.shape[0]), D)
    AII = A[I][:, I].tocsr()
    AID = A[I][:, D].tocsr()
    ps = pypardiso.PyPardisoSolver()
    ps.factorize(AII)
    keep = set(range(nsteps + 1)) if save is None else set(np.linspace(0, nsteps, save).round().astype(int))
    ts, Ts = [0.0], [T.copy()]
    for n in range(1, nsteps + 1):
        t = n * dt
        rhs = B @ T + f
        _, TD = _fixed(ib, T_fix, t)
        T[D] = TD
        T[I] = ps.solve(AII, rhs[I] - AID @ TD)
        if n in keep:
            ts.append(t)
            Ts.append(T.copy())
    return dict(basis=ib, t=np.array(ts), T=Ts, dt=dt, DOF=K.shape[0], solver="pypardiso")
