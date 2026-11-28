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


SIGMA = 5.670374419e-14        # Stefan-Boltzmann i mW/(mm²·K⁴) (5,670e-8 W/(m²·K⁴))


def heat_transient_coupled(m, k, rho, c, T0, t_end, nsteps, conv=None, rad=None, lumped=None,
                           T_abs=273.15, save=None, refactor_tol=0.03):
    """Transient värmeledning (implicit Euler) med tids- och temperaturberoende randvillkor,
    lokalt linjariserad strålning och en klumpad nod (t.ex. välblandad vätska) kopplad via randerna.

    conv: {rand: fn(s) -> [(h, Tref), ...]} eller en fast lista [(h, Tref), ...]. h i mW/(mm²·K).
          Tref är en temperatur eller "L" (den klumpade noden). s = dict(t, Ts, TL, A): tid,
          areamedeltemperatur på randen och nodens temperatur från föregående steg, randens area [mm²].
    rad:  {rand: (eps, T_sur)} grå strålning mot omgivning med temperaturen T_sur. Linjariseras per
          integrationspunkt med föregående steg, h_r = eps·σ·(T² + T_sur²)(T + T_sur) i absolut temperatur
          (T + T_abs). Temperaturer i °C med T_abs = 273,15 eller i K med T_abs = 0.
    lumped: dict(C, T0, loss) med C [mW·s/K] = m·c_p, och loss(TL, t) -> (P, dP/dTL): förlust från noden
          till omgivningen [mW], linjariseras kring föregående steg (semi-implicit).
    Systemmatrisen faktoriseras om när någon randkoefficient har ändrats mer än refactor_tol sedan
    senaste faktoriseringen; däremellan tas skillnaden med explicit korrektion.
    save: antal jämnt fördelade sparade fält (default alla steg).

    Returnerar dict: basis, t och T (sparade P2-vektorer), tt, TL, Q (dict per rand och "lumped" med
    flöde ut ur systemet per steg [mW]; "L:<rand>" = från noden in i kroppen), E (lagrad energi relativt
    0 grader [mW·s]), balance (|ΔE + ∫ΣQ_ut dt| / max(|ΔE|, ∫|ΣQ_ut| dt)), DOF, refactor (antal)."""
    import pypardiso
    import scipy.sparse as sp
    ib = Basis(m, ElementTetP2(), intorder=4)
    K0 = asm(BilinearForm(lambda u, v, w: k * dot(grad(u), grad(v))), ib).tocsr()
    M = asm(BilinearForm(lambda u, v, w: rho * c * u * v), ib).tocsr()
    N = K0.shape[0]
    names = sorted(set(conv or {}) | set(rad or {}))
    fb, Mb, bv, A = {}, {}, {}, {}
    for n in names:
        fb[n] = FacetBasis(m, ElementTetP2(), facets=m.boundaries[n], intorder=4)
        Mb[n] = asm(BilinearForm(lambda u, v, w: u * v), fb[n]).tocsr()
        bv[n] = asm(LinearForm(lambda v, w: v + 0 * w.x[0]), fb[n])
        A[n] = bv[n].sum()
    hasL = lumped is not None
    nL = N + 1 if hasL else N
    dt = t_end / nsteps
    T = np.full(N, float(T0)) if np.isscalar(T0) else np.asarray(T0, float).copy()
    TL = float(lumped["T0"]) if hasL else 0.0
    CL = float(lumped["C"]) if hasL else 0.0
    keep = set(range(nsteps + 1)) if save is None else set(np.linspace(0, nsteps, save).round().astype(int))
    Md = M / dt
    one = np.ones(N)
    E = lambda T, TL: one @ (M @ T) + CL * TL
    E0 = E(T, TL)
    ts, Ts_, tt, TLs = [0.0], [T.copy()], [0.0], [TL]
    Q = {}
    fact = None
    ps = pypardiso.PyPardisoSolver()
    nref = 0
    rad_form_K = BilinearForm(lambda u, v, w: w.hr * u * v)
    rad_form_f = LinearForm(lambda v, w: w.hr * w.Tsur * v)
    sumQ = 0.0
    sumabs = 0.0
    for step in range(1, nsteps + 1):
        t = step * dt
        K = K0.copy()
        f = np.zeros(nL)
        cvec = np.zeros(N)
        KLL = 0.0
        fL = 0.0
        coeff = []
        terms = {}                                   # rand -> (lista av (h, Tref)), strålningsmatris
        for n in names:
            s = dict(t=t, Ts=bv[n] @ T / A[n], TL=TL, A=A[n])
            lst = []
            if conv and n in conv:
                lst = conv[n](s) if callable(conv[n]) else conv[n]
            for h, Tref in lst:
                K = K + h * Mb[n]
                if isinstance(Tref, str):
                    cvec -= h * bv[n]
                    KLL += h * A[n]
                else:
                    f[:N] += h * Tref * bv[n]
                coeff.append(h)
            Kr = fr = None
            if rad and n in rad:
                eps, Tsur = rad[n]
                Tp = fb[n].interpolate(T) + T_abs
                Sa = Tsur + T_abs
                hr = eps * SIGMA * (Tp ** 2 + Sa ** 2) * (Tp + Sa)
                Kr = asm(rad_form_K, fb[n], hr=hr).tocsr()
                fr = asm(rad_form_f, fb[n], hr=hr, Tsur=Tsur + 0 * hr)
                K = K + Kr
                f[:N] += fr
                coeff.append(float(np.mean(hr)))
            terms[n] = (lst, Kr, fr)
        if hasL:
            P, dP = lumped["loss"](TL, t)
            KLL += dP
            fL = -P + dP * TL
            coeff.append(dP)
            f[N] = CL / dt * TL + fL
            Amat = sp.bmat([[Md + K, sp.csr_matrix(cvec[:, None])],
                            [sp.csr_matrix(cvec[None, :]), sp.csr_matrix([[CL / dt + KLL]])]], format="csr")
            x0 = np.append(T, TL)
        else:
            Amat = (Md + K).tocsr()
            x0 = T
        f[:N] += Md @ T
        coeff = np.array(coeff)
        if fact is None or len(coeff) != len(fact[1]) or np.max(np.abs(coeff - fact[1]) / np.maximum(np.abs(fact[1]), 1e-30)) > refactor_tol:
            if fact is not None:
                ps.free_memory(everything=True)              # annars läcker MKL-minne vid varje omfaktorisering
            ps.factorize(Amat)
            fact = (ps, coeff, Amat)
            nref += 1
            x = ps.solve(Amat, f)
        else:
            ps, _, Af = fact
            x = ps.solve(Af, f - (Amat - Af) @ x0)
        T = x[:N].copy()
        TLn = x[N] if hasL else 0.0
        # flöden ut ur systemet vid det nya tillståndet
        qsum = 0.0
        qs = {}
        for n in names:
            lst, Kr, fr = terms[n]
            bt = bv[n] @ T
            for i, (h, Tref) in enumerate(lst):
                key = n if not isinstance(Tref, str) else "L:" + n
                q = h * (bt - Tref * A[n]) if not isinstance(Tref, str) else h * (TLn * A[n] - bt)
                qs[key] = qs.get(key, 0.0) + q
                if not isinstance(Tref, str):
                    qsum += q
            if Kr is not None:
                q = (Kr @ T).sum() - fr.sum()
                qs["rad:" + n] = q
                qsum += q
        if hasL:
            q = P + dP * (TLn - TL)
            qs["lumped"] = q
            qsum += q
            TL = TLn
        for key, q in qs.items():
            Q.setdefault(key, []).append(q)
        sumQ += qsum * dt
        sumabs += abs(qsum) * dt
        tt.append(t)
        TLs.append(TL)
        if step in keep:
            ts.append(t)
            Ts_.append(T.copy())
    ps.free_memory(everything=True)
    dE = E(T, TL) - E0
    bal = abs(dE + sumQ) / max(abs(dE), sumabs, 1e-30)
    return dict(basis=ib, t=np.array(ts), T=Ts_, tt=np.array(tt), TL=np.array(TLs),
                Q={k_: np.array(v) for k_, v in Q.items()}, E0=E0, dE=dE, balance=bal, dt=dt,
                DOF=nL, refactor=nref, solver="pypardiso")
