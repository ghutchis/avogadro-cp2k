# avogadro-cp2k: minimal replacement for upstream's blocks/__init__.py, which
# imports every block module and a pint unit registry. Only the Mulliken block
# is vendored; avogadro_cp2k.blocks supplies the rest. See VENDORED.md.
from .mulliken import match_mulliken_population_analysis

builtin_matchers = [match_mulliken_population_analysis]
