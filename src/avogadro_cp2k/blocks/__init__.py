#  This source file is part of the Avogadro project.
#  This source code is released under the 3-Clause BSD License, (see "LICENSE").

"""CP2K output matchers in the style of cp2k-output-tools.

Every matcher has the signature ``(content: str) -> Optional[BlockMatch]`` and
returns the LAST occurrence of its block, since the final step of a multi-step
run (GEO_OPT, CELL_OPT, MD) is the one that matters.
"""

from ..vendored.cp2k_output_tools.blocks.mulliken import (
    match_mulliken_population_analysis,
)
from .cell import match_cell
from .energy import match_total_energy
from .geometry import match_atomic_coordinates, match_atom_count
from .hirshfeld import match_hirshfeld_charges
from .run_info import match_run_info
from .tddfpt import match_tddfpt

builtin_matchers = [
    match_run_info,
    match_cell,
    match_atom_count,
    match_atomic_coordinates,
    match_total_energy,
    match_mulliken_population_analysis,
    match_hirshfeld_charges,
    match_tddfpt,
]

__all__ = ["builtin_matchers"]
