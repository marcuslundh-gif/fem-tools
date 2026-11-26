"""femtools – verifierade FEM-mallar (mm, N, MPa, t, s). Verifiering: python -m femtools.verify
Regler: docs/projektstandard.md. Register över verifierade verktyg: docs/verifiering.md."""
from .core import static, modal, richardson, von_mises
from .mesh import setup, refine_and_mesh, to_skfem, mesh_step
from .materials import material, materials
from . import projektblad
__version__ = "1.2.0"
