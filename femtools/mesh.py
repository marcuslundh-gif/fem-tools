"""Effektiv meshning med gmsh -> scikit-fem (mm, N, MPa).

Principer:
- HXT 3D-algoritm (parallell), Frontal-Delaunay 2D, flera trådar.
- Global storlek h_glob + lokal förfining h_loc runt utvalda ytor (Distance+Threshold).
- Storlek från krökning (element per 2π) så att radier och hål löses upp.
- Kvalitetskontroll (gamma/SICN) innan lösning.
- P1-nät från gmsh, P2-basis i skfem (raka kanter; förfina vid krökta kritiska ytor).

Användning:
  from femtools import mesh_step          # STEP-fil -> skfem-nät
  python -m femtools.verify               # verifieringsfall: plåt med hål, Kt mot Heywood
"""
import time
import numpy as np
import gmsh


def setup(threads=0, occ_fix=True):
    """occ_fix=False stänger av OCC-reparationen. Krävs för svepta kroppar (t.ex. ett öra svept längs en bana
    med fillets) där OCCFixSmallFaces/-Edges ändrar soliden: volymen sjönk 0,9 % och occ.fragment gav fel delar."""
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.option.setNumber("General.NumThreads", threads)      # 0 = alla kärnor
    gmsh.option.setNumber("Mesh.Algorithm", 6)                # Frontal-Delaunay 2D
    gmsh.option.setNumber("Mesh.Algorithm3D", 10)             # HXT (parallell)
    gmsh.option.setNumber("Mesh.Optimize", 1)
    fix = 1 if occ_fix else 0
    gmsh.option.setNumber("Geometry.OCCFixDegenerated", fix)
    gmsh.option.setNumber("Geometry.OCCFixSmallEdges", fix)
    gmsh.option.setNumber("Geometry.OCCFixSmallFaces", fix)


def refine_and_mesh(h_glob, h_loc=None, refine_surfs=(), dist_min=None, dist_max=None,
                    n_curv=24):
    gmsh.option.setNumber("Mesh.MeshSizeMax", h_glob)
    gmsh.option.setNumber("Mesh.MeshSizeMin", (h_loc or h_glob) / 2)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", n_curv)  # element per 2π
    if h_loc and refine_surfs:
        f_d = gmsh.model.mesh.field.add("Distance")
        gmsh.model.mesh.field.setNumbers(f_d, "SurfacesList", list(refine_surfs))
        f_t = gmsh.model.mesh.field.add("Threshold")
        gmsh.model.mesh.field.setNumber(f_t, "InField", f_d)
        gmsh.model.mesh.field.setNumber(f_t, "SizeMin", h_loc)
        gmsh.model.mesh.field.setNumber(f_t, "SizeMax", h_glob)
        gmsh.model.mesh.field.setNumber(f_t, "DistMin", dist_min or 2 * h_loc)
        gmsh.model.mesh.field.setNumber(f_t, "DistMax", dist_max or 10 * h_loc)
        gmsh.model.mesh.field.setAsBackgroundMesh(f_t)
        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    t0 = time.perf_counter()
    gmsh.model.mesh.generate(3)
    t_mesh = time.perf_counter() - t0
    # Kvalitet: gamma (inradie/omradie-baserad), 1 = perfekt
    _, tets, _ = gmsh.model.mesh.getElements(3)
    q = gmsh.model.mesh.getElementQualities(tets[0], "gamma")
    return t_mesh, dict(n_tet=len(tets[0]), gamma_min=float(q.min()),
                        gamma_p1=float(np.percentile(q, 1)),
                        frac_bad=float((q < 0.1).mean()))


def to_skfem():
    from skfem import MeshTet
    tags, xyz, _ = gmsh.model.mesh.getNodes()
    xyz = xyz.reshape(-1, 3)
    idx = np.zeros(int(tags.max()) + 1, dtype=np.int64)
    idx[tags.astype(np.int64)] = np.arange(len(tags))
    _, ntags = gmsh.model.mesh.getElementsByType(4)            # 4 = tet4
    t = idx[ntags.astype(np.int64)].reshape(-1, 4).T
    used = np.unique(t)                                        # ta bort lösa noder
    remap = -np.ones(len(xyz), dtype=np.int64)
    remap[used] = np.arange(len(used))
    return MeshTet(xyz[used].T, remap[t])


def mesh_step(path, h_glob, h_loc=None, refine_surfs=(), threads=0, occ_fix=True):
    """Läs STEP, kontrollera enheter via bounding box, meshar. Returnerar skfem-nät."""
    setup(threads, occ_fix)
    gmsh.model.occ.importShapes(path)
    gmsh.model.occ.synchronize()
    bb = gmsh.model.getBoundingBox(-1, -1)
    print("bounding box [mm?]:", np.round(bb, 3))            # kontrollera mm vs m
    t, q = refine_and_mesh(h_glob, h_loc, refine_surfs)
    m = to_skfem()
    gmsh.finalize()
    return m, t, q


