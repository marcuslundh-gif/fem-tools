"""femtools – verifierade FEM-mallar (mm, N, MPa, t, s). Verifiering: python -m femtools.verify
Regler: docs/projektstandard.md. Register över verifierade verktyg: docs/verifiering.md."""
from .core import static, modal, richardson, von_mises
from .loads import static_ext
from .thermal import heat, heat_transient, heat_transient_coupled
from .buckling import buckling
from .cfd import flow, boundary_flow, reynolds
from . import openfoam
from .mesh import setup, refine_and_mesh, to_skfem, mesh_step
from .materials import material, materials
from . import projektblad, report, webreport, shell
__version__ = "1.6.0"
