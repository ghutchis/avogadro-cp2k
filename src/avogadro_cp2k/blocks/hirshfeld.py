#  This source file is part of the Avogadro project.
#  This source code is released under the 3-Clause BSD License, (see "LICENSE").

"""Hirshfeld charges (``Hirshfeld Charges`` block of the population print)."""

from typing import Optional

import regex as re

from ..vendored.cp2k_output_tools.blocks.common import FLOAT, BlockMatch
from ._last import last_match

# Restricted:
#
#                           Hirshfeld Charges
#
#   #Atom  Element  Kind  Ref Charge     Population                    Net charge
#       1       S      1       6.000          6.110                        -0.110
#   ...
#   Total Charge                                                           -0.000
#
# Unrestricted rows carry four values after the reference charge (population
# alpha, population beta, spin moment, net charge). The trailing numbers are
# captured as ONE group per row and counted afterwards; the net charge is the
# last column in both cases.
HIRSHFELD_CHARGES_RE = re.compile(
    rf"""
^[ \t]* Hirshfeld\ Charges [ \t]* \n
(?:[ \t]*\n)*
^[ \t]* \# [^\n]* \n
(?:
  ^[ \t]*
  (?P<atom>\d+) [ \t]+
  (?P<element>\w+) [ \t]+
  (?P<kind>\d+) [ \t]+
  (?P<reference>{FLOAT}) [ \t]+
  (?P<values>{FLOAT}(?:[ \t]+{FLOAT})*) [ \t]* \n
)+
(?:[ \t]*\n)*
^[ \t]* Total\ Charge [ \t]+ (?P<total>{FLOAT}(?:[ \t]+{FLOAT})*) [ \t]* \n
""",
    re.VERSION1 | re.MULTILINE | re.VERBOSE,
)


def match_hirshfeld_charges(content: str) -> Optional[BlockMatch]:
    # avogadro-cp2k: last block (the final step of a multi-step run).
    match = last_match(HIRSHFELD_CHARGES_RE, content)
    if match is None:
        return None

    captures = match.capturesdict()
    per_atom = []
    for i in range(len(captures["atom"])):
        values = [float(v) for v in captures["values"][i].split()]
        entry = {
            "element": captures["element"][i],
            "kind": int(captures["kind"][i]),
            "reference charge": float(captures["reference"][i]),
            "charge": values[-1],  # net charge: always the last column
        }
        if len(values) == 2:  # restricted: population, net charge
            entry["population"] = values[0]
        elif len(values) == 4:  # unrestricted
            entry["population alpha"] = values[0]
            entry["population beta"] = values[1]
            entry["spin"] = values[2]
        else:
            return None  # a layout we do not know; do not guess columns
        per_atom.append(entry)

    total = [float(v) for v in captures["total"][0].split()]
    return BlockMatch(
        {"hirshfeld charges": {"per atom": per_atom, "total charge": total[-1]}},
        match.span(),
    )
