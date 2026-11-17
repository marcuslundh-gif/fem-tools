"""femtools – verifierade FEM-mallar (mm, N, MPa, t, s). Verifiering: python -m femtools.verify"""
from .core import static, modal, richardson, von_mises
from .mesh import setup, refine_and_mesh, to_skfem, mesh_step
__version__ = "1.0.0"
