#  This source file is part of the Avogadro project.
#  This source code is released under the 3-Clause BSD License, (see "LICENSE").

"""Program banner, GLOBAL run settings and the DFT| spin/charge header."""

from typing import Optional

import regex as re

from ..vendored.cp2k_output_tools.blocks.common import BlockMatch
from ._last import last_match

_FLAGS = re.VERSION1 | re.MULTILINE | re.VERBOSE

VERSION_RE = re.compile(
    r"^[ \t]* CP2K\| [ \t]+ version\ string: [ \t]+ (?P<version>.*?) [ \t]* $",
    _FLAGS,
)
RUN_TYPE_RE = re.compile(
    r"^[ \t]* GLOBAL\| [ \t]+ Run\ type [ \t]+ (?P<run_type>\S+) [ \t]* $", _FLAGS
)
PROJECT_RE = re.compile(
    r"^[ \t]* GLOBAL\| [ \t]+ Project\ name [ \t]+ (?P<project>\S.*?) [ \t]* $",
    _FLAGS,
)
# " DFT| Spin restricted Kohn-Sham (RKS) calculation      RKS"
SPIN_TREATMENT_RE = re.compile(
    r"^[ \t]* DFT\| [ \t]+ Spin\ [^\n]*? \((?P<kind>[RU]KS)\)", _FLAGS
)
MULTIPLICITY_RE = re.compile(
    r"^[ \t]* DFT\| [ \t]+ Multiplicity [ \t]+ (?P<value>[+\-]?\d+) [ \t]* $", _FLAGS
)
CHARGE_RE = re.compile(
    r"^[ \t]* DFT\| [ \t]+ Charge [ \t]+ (?P<value>[+\-]?\d+) [ \t]* $", _FLAGS
)


def match_run_info(content: str) -> Optional[BlockMatch]:
    data = {}
    spans = []

    # avogadro-cp2k: version/run type/project are one-per-file, so the last
    # occurrence is only a safeguard against concatenated outputs.
    for key, pattern, group in (
        ("version", VERSION_RE, "version"),
        ("run type", RUN_TYPE_RE, "run_type"),
        ("project", PROJECT_RE, "project"),
    ):
        match = last_match(pattern, content)
        if match:
            data[key] = match[group]
            spans.append(match.span())

    match = last_match(SPIN_TREATMENT_RE, content)
    if match:
        data["spin treatment"] = match["kind"]
        spans.append(match.span())

    for key, pattern in (("multiplicity", MULTIPLICITY_RE), ("charge", CHARGE_RE)):
        match = last_match(pattern, content)
        if match:
            data[key] = int(match["value"])
            spans.append(match.span())

    if not data:
        return None
    return BlockMatch({"run info": data}, spans)
