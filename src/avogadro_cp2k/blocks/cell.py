#  This source file is part of the Avogadro project.
#  This source code is released under the 3-Clause BSD License, (see "LICENSE").

"""The ``CELL|`` block (never ``CELL_REF|`` / ``CELL_TOP|``)."""

from typing import Optional

import regex as re

from ..utils import BOHR_TO_ANGSTROM
from ..vendored.cp2k_output_tools.blocks.common import FLOAT, BlockMatch
from ._last import last_match

_FLAGS = re.VERSION1 | re.MULTILINE | re.VERBOSE

# A run of consecutive " CELL| ..." lines. "CELL_REF|" does not match
# because the character after "CELL" must be a literal pipe.
CELL_RUN_RE = re.compile(r"(?:^[ \t]* CELL\| [^\n]* \n?)+", _FLAGS)

# " CELL| Vector a [angstrom]:  9.215  0.000  0.000  |a| =  9.215000"
VECTOR_RE = re.compile(
    rf"""
^[ \t]* CELL\| [ \t]+ Vector [ \t]+ (?P<axis>[abc]) [ \t]+ \[(?P<unit>[^\]\n]*)\] :?
[ \t]+ (?P<x>{FLOAT}) [ \t]+ (?P<y>{FLOAT}) [ \t]+ (?P<z>{FLOAT})
[ \t]+ \|[abc]\| [ \t]* = [ \t]* (?P<norm>{FLOAT}) [ \t]* $
""",
    _FLAGS,
)
# " CELL| Angle (b,c), alpha [degree]:    113.320000"
ANGLE_RE = re.compile(
    rf"""
^[ \t]* CELL\| [ \t]+ Angle [ \t]+ \( [abc],[abc] \), [ \t]+ (?P<name>alpha|beta|gamma)
[ \t]+ \[degree\] :? [ \t]+ (?P<value>{FLOAT}) [ \t]* $
""",
    _FLAGS,
)
PERIODICITY_RE = re.compile(
    r"^[ \t]* CELL\| [ \t]+ Periodicity [ \t]+ (?P<periodicity>\w+) [ \t]* $", _FLAGS
)

_UNIT_TO_ANGSTROM = {"angstrom": 1.0, "bohr": BOHR_TO_ANGSTROM}


def match_cell(content: str) -> Optional[BlockMatch]:
    # avogadro-cp2k: the last CELL| block that actually holds three vectors.
    run = None
    for candidate in CELL_RUN_RE.finditer(content):
        if len(VECTOR_RE.findall(candidate[0])) >= 3:
            run = candidate
    if run is None:
        return None

    text = run[0]
    vectors = {}
    norms = {}
    for m in VECTOR_RE.finditer(text):
        unit = m["unit"].strip().lower()
        scale = _UNIT_TO_ANGSTROM.get(unit)
        if scale is None:
            return None  # an unknown length unit: better no cell than a wrong one
        vectors[m["axis"]] = [float(m["x"]) * scale, float(m["y"]) * scale, float(m["z"]) * scale]
        norms[m["axis"]] = float(m["norm"]) * scale

    if set(vectors) != {"a", "b", "c"}:
        return None

    angles = {m["name"]: float(m["value"]) for m in ANGLE_RE.finditer(text)}
    periodicity = last_match(PERIODICITY_RE, text)

    data = {
        "vectors": [vectors["a"], vectors["b"], vectors["c"]],  # angstrom
        "lengths": [norms["a"], norms["b"], norms["c"]],  # angstrom
        "periodicity": periodicity["periodicity"] if periodicity else None,
    }
    if set(angles) == {"alpha", "beta", "gamma"}:
        data["angles"] = [angles["alpha"], angles["beta"], angles["gamma"]]  # degrees

    return BlockMatch({"cell": data}, run.span())
