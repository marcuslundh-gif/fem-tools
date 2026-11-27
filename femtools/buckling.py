"""Linjär knäckning (egenvärdesknäckning), P2-tetraedrar, mm-N-MPa.

    from femtools import buckling
    r = buckling(m, E, nu, k=4, fix="bas", traction={"topp": [0, 0, -1.0]})
    r["factors"]   # lastfaktorer λ (stigande, positiva); kritisk last = λ × pålagd last
    r["modes"]     # knäckformer (DOF-vektorer, normerade till max |u| = 1)
    r["static"]    # förspänningsfallet från static_ext (spänning, reaktioner)

(K + λ·K_σ)·φ = 0, där K_σ = ∫ σ0_ij ∂_i u_k ∂_j v_k dV byggs av förspänningen σ0 från static_ext med samma
laster och randvillkor (alla nyckelord till static_ext fungerar, även rho/gravity/omega och alpha/dT).
λ skalar ALLA laster i förspänningsfallet. Negativa λ (knäckning vid omvänd last) tas bort.
Linjär knäckning överskattar bärförmågan för imperfektionskänsliga konstruktioner (tunna skal, cylindrar
under axialtryck): använd reduktionsfaktorer eller icke-linjär analys där.
"""
import numpy as np
from skfem import ElementVector, ElementTetP2, Basis, asm, BilinearForm
from skfem.models.elasticity import linear_elasticity, lame_parameters
from scipy.sparse.linalg import eigsh, LinearOperator
import scipy.sparse as sp
from .loads import static_ext, _dofs


def buckling(m, E, nu, k=4, fix="fix", D=None, **loads):
    r0 = static_ext(m, E, nu, fix=fix, D=D, **loads)
    e = ElementVector(ElementTetP2())
    ib = Basis(m, e, intorder=2)
    K = asm(linear_elasticity(*lame_parameters(E, nu)), ib)
    vb = r0["basis"]                                          # intorder 4, samma som spänningen i r0
    s = r0["sigma"]

    @BilinearForm
    def geometric(u, v, w):                                   # skalär: σ0_ij ∂_i u ∂_j v
        return np.einsum("ijeq,ieq,jeq->eq", w.s, u.grad, v.grad, optimize=True)

    # K_σ är blockdiagonal (samma skalära block för varje komponent); vektor-DOF = 3·skalär + komponent
    KGs = asm(geometric, vb.with_element(ElementTetP2()), s=s)
    KG = sp.kron(KGs, sp.identity(3), format="csr")
    Dd = _dofs(ib, fix, D)
    I = np.setdiff1d(np.arange(K.shape[0]), Dd)
    A, G = K[I][:, I].tocsr(), KG[I][:, I].tocsr()
    import pypardiso
    ps = pypardiso.PyPardisoSolver()
    ps.factorize(A)
    assert np.array_equal(ib.split_indices()[1], 3 * np.arange(K.shape[0] // 3) + 1)   # DOF-ordning för kron
    ps.set_phase(33)                                          # bara lösning, faktoriseringen återanvänds
    op = LinearOperator(A.shape, matvec=lambda x: ps._call_pardiso(A, np.ascontiguousarray(x, dtype=float).ravel()),
                        dtype=float)
    # −K_σ φ = μ K φ, μ = 1/λ; största positiva μ ger minsta positiva λ
    from threadpoolctl import threadpool_limits
    with threadpool_limits(limits=1, user_api="blas"):        # BLAS-trådar som spinner i ARPACK tävlar annars
        mu, V = eigsh(-G, k=min(2 * k + 2, A.shape[0] - 2), M=A, Minv=op, which="LA")   # med Pardiso: 2,5x långsammare
    keep = mu > 0
    lam, V = 1.0 / mu[keep], V[:, keep]
    order = np.argsort(lam)[:k]
    modes = []
    for j in order:
        x = np.zeros(K.shape[0])
        x[I] = V[:, j]
        modes.append(x / np.abs(x).max())
    return dict(factors=lam[order], modes=modes, basis=ib, static=r0, DOF=K.shape[0])
