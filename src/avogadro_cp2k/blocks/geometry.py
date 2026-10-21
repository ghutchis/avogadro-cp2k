#  This source file is part of the Avogadro project.
#  This source code is released under the 3-Clause BSD License, (see "LICENSE").

"""Quickstep atomic coordinates and the total atom count."""

from typing import Optional

import regex as re

from ..utils import BOHR_TO_ANGSTROM
from ..vendored.cp2k_output_tools.blocks.common import FLOAT, BlockMatch
from ._last import last_match

_FLAGS = re.VERSION1 | re.MULTILINE | re.VERBOSE

# MODULE QUICKSTEP: ATOMIC COORDINATES IN ANGSTROM
#
#   Atom Kind Element         X             Y             Z       Z(eff)     Mass
#      1    1 S    16      1.086966      1.304941      1.378422   6.0000  32.0650
#
# Nine columns: index, kind, element symbol, ATOMIC NUMBER, x, y, z, Z(eff), mass.
ATOMIC_COORDINATES_RE = re.compile(
    rf"""
^[ \t]* MODULE\ QUICKSTEP: [ \t]+ ATOMIC\ COORDINATES\ IN [ \t]+ (?P<unit>\w+) [ \t]* \n
(?:[ \t]*\n)*
^[ \t]* Atom [ \t]+ Kind [ \t]+ Element [^\n]* \n
(?:
  ^[ \t]*
  (?P<atom>\d+) [ \t]+
  (?P<kind>\d+) [ \t]+
  (?P<element>[A-Za-z][^\s]*) [ \t]+
  (?P<number>\d+) [ \t]+
  (?P<x>{FLOAT}) [ \t]+ (?P<y>{FLOAT}) [ \t]+ (?P<z>{FLOAT}) [ \t]+
  (?P<zeff>{FLOAT}) [ \t]+ (?P<mass>{FLOAT}) [ \t]* \n
)+
""",
    _FLAGS,
)

# " Total number of            - Atomic kinds:      4"
# "                            - Atoms:            46"
ATOM_COUNT_RE = re.compile(r"^[ \t]* - [ \t]+ Atoms: [ \t]+ (?P<count>\d+) [ \t]* $", _FLAGS)

_UNIT_TO_ANGSTROM = {"angstrom": 1.0, "bohr": BOHR_TO_ANGSTROM}


def match_atomic_coordinates(content: str) -> Optional[BlockMatch]:
    # avogadro-cp2k: last block, not the first.
    match = last_match(ATOMIC_COORDINATES_RE, content)
    if match is None:
        return None

    scale = _UNIT_TO_ANGSTROM.get(match["unit"].lower())
    if scale is None:
        return None  # unknown unit: refuse rather than guess

    captures = match.capturesdict()
    atoms = []
    for i in range(len(captures["atom"])):
        atoms.append(
            {
                "element": captures["element"][i],
                "kind": int(captures["kind"][i]),
                "atomic number": int(captures["number"][i]),
                # always angstrom in the returned data
                "position": [
                    float(captures["x"][i]) * scale,
                    float(captures["y"][i]) * scale,
                    float(captures["z"][i]) * scale,
                ],
                "z_eff": float(captures["zeff"][i]),
                "mass": float(captures["mass"][i]),
            }
        )
    return BlockMatch({"atomic coordinates": {"unit": "angstrom", "atoms": atoms}}, match.span())


def match_atom_count(content: str) -> Optional[BlockMatch]:
    """The ``- Atoms:`` line of TOTAL NUMBERS AND MAXIMUM NUMBERS, used to
    check that the coordinate table was read completely."""
    match = last_match(ATOM_COUNT_RE, content)
    if match is None:
        return None
    return BlockMatch({"atom count": int(match["count"])}, match.span())
