"""femtools – verifierade FEM-mallar (mm, N, MPa, t, s). Verifiering: python -m femtools.verify
Regler: docs/projektstandard.md. Register över verifierade verktyg: docs/verifiering.md."""
from .core import static, modal, richardson, von_mises
from .loads import static_ext
from .thermal import heat, heat_transient
from .buckling import buckling
from .mesh import setup, refine_and_mesh, to_skfem, mesh_step
from .materials import material, materials
from . import projektblad, report
__version__ = "1.3.0"
