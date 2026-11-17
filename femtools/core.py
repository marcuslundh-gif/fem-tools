import numpy as np
from skfem import ElementVector, ElementTetP2, Basis, FacetBasis, asm, condense, solve, LinearForm, BilinearForm
from skfem.helpers import dot
from skfem.models.elasticity import linear_elasticity, lame_parameters

try:
    import pypardiso
    SOLVER = lambda A, b: pypardiso.spsolve(A.tocsr(), b)
    SOLVER_NAME = "pypardiso"
except ImportError:  # pragma: no cover
    from scipy.sparse.linalg import spsolve
    SOLVER = lambda A, b: spsolve(A.tocsc(), b)
    SOLVER_NAME = "scipy-spsolve (en tråd)"


def von_mises(ib, u, E, nu):
    """von Mises i integrationspunkter, form (nel, nqp)."""
    G = ib.interpolate(u).grad
    eps = 0.5 * (G + G.transpose(1, 0, 2, 3))
    lam, mu = lame_parameters(E, nu)
    s = 2 * mu * eps + lam * (eps[0, 0] + eps[1, 1] + eps[2, 2]) * np.eye(3)[:, :, None, None]
    return np.sqrt(0.5 * ((s[0, 0] - s[1, 1])**2 + (s[1, 1] - s[2, 2])**2 + (s[2, 2] - s[0, 0])**2)
                   + 3 * (s[0, 1]**2 + s[1, 2]**2 + s[0, 2]**2)), s


def static(m, E, nu, traction, fix="fix", load="load", D=None):
    """Linjär statik, P2. traction [MPa] på rand `load`, låsning på rand `fix` (eller egna DOF i D).
    Returnerar dict med basis, u, vm, sigma, xq, R (reaktioner), F (pålagd last), DOF."""
    e = ElementVector(ElementTetP2())
    ib = Basis(m, e, intorder=4)
    K = asm(linear_elasticity(*lame_parameters(E, nu)), ib)
    t = np.asarray(traction, float)
    fb = FacetBasis(m, e, facets=m.boundaries[load], intorder=4)

    @LinearForm
    def trac(v, w):
        return t[0] * v[0] + t[1] * v[1] + t[2] * v[2]

    f = asm(trac, fb)
    Dd = ib.get_dofs(fix) if D is None else D
    u = solve(*condense(K, f, D=Dd), solver=SOLVER)
    r = K @ u - f
    alld = ib.get_dofs()
    R = [r[alld.all(f"u^{i+1}")][np.isin(alld.all(f"u^{i+1}"), np.asarray(Dd.all() if hasattr(Dd, "all") else Dd))].sum()
         for i in range(3)]
    F = [f[alld.all(f"u^{i+1}")].sum() for i in range(3)]
    vm, s = von_mises(ib, u, E, nu)
    return dict(basis=ib, u=u, vm=vm, sigma=s, xq=ib.mapping.F(ib.X), R=R, F=F, DOF=K.shape[0], solver=SOLVER_NAME)


def modal(m, E, nu, rho, k=6, fix="fix"):
    """Egenfrekvenser [Hz] med rho i t/mm³."""
    from scipy.sparse.linalg import eigsh
    ib = Basis(m, ElementVector(ElementTetP2()))
    K = asm(linear_elasticity(*lame_parameters(E, nu)), ib)

    @BilinearForm
    def mass(u, v, w):
        return rho * dot(u, v)

    M = asm(mass, ib)
    I = ib.complement_dofs(ib.get_dofs(fix))
    w2, _ = eigsh(K[I][:, I], k=k, M=M[I][:, I], sigma=0)
    return np.sort(np.sqrt(np.abs(w2)) / (2 * np.pi))


def richardson(f1, f2, f3, r):
    """Observerad ordning p och extrapolerat värde (grov→fin, förfiningsfaktor r)."""
    p = np.log(abs((f1 - f2) / (f2 - f3))) / np.log(r)
    return p, f3 + (f3 - f2) / (r**p - 1)
