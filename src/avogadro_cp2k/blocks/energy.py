#  This source file is part of the Avogadro project.
#  This source code is released under the 3-Clause BSD License, (see "LICENSE").

"""The final ``ENERGY| Total FORCE_EVAL`` line."""

from typing import Optional

import regex as re

from ..vendored.cp2k_output_tools.blocks.common import FLOAT, BlockMatch
from ._last import last_match

# " ENERGY| Total FORCE_EVAL ( QS ) energy [hartree]   -233.7225..."
# Older releases print "energy (a.u.):" instead of "[hartree]".
TOTAL_ENERGY_RE = re.compile(
    rf"""
^[ \t]* ENERGY\| [ \t]+ Total\ FORCE_EVAL [ \t]* \( [ \t]* (?P<method>\w+) [ \t]* \)
[ \t]+ energy [ \t]* (?:\[hartree\] | \(a\.u\.\) [ \t]* :?) [ \t]+
(?P<energy>{FLOAT}) [ \t]* $
""",
    re.VERSION1 | re.MULTILINE | re.VERBOSE,
)


def match_total_energy(content: str) -> Optional[BlockMatch]:
    match = last_match(TOTAL_ENERGY_RE, content)
    if match is None:
        return None
    return BlockMatch(
        {"total energy": {"method": match["method"], "hartree": float(match["energy"])}},
        match.span(),
    )
