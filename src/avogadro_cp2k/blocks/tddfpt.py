#  This source file is part of the Avogadro project.
#  This source code is released under the 3-Clause BSD License, (see "LICENSE").

"""TDDFPT excited-state tables and their consistency check.

The printed table looks like::

   R-TDDFPT states of multiplicity 1
   Transition dipoles calculated using velocity formulation

           State    Excitation        Transition dipole (a.u.)        Oscillator
           number   energy (eV)       x           y           z     strength (a.u.)
           ------------------------------------------------------------------------
   TDDFPT|      1       1.51072   1.1066E-01  3.6585E-02 -7.3801E-03   1.63771E-01

followed (usually) by an "Excitation analysis" table of occupied -> virtual
amplitudes per state. Unrestricted runs print ``U-TDDFPT states``.
"""

import math
from typing import Optional

import regex as re

from ..utils import HARTREE_TO_EV
from ..vendored.cp2k_output_tools.blocks.common import FLOAT, BlockMatch

_FLAGS = re.VERSION1 | re.MULTILINE | re.VERBOSE

TDDFPT_BLOCK_RE = re.compile(
    rf"""
^[ \t]* (?P<spin>[RU])-TDDFPT\ states (?:\ of\ multiplicity [ \t]+ (?P<multiplicity>\d+))? [ \t]* \n
(?:[ \t]*\n)*
(?:^[ \t]* Transition\ dipoles\ calculated\ using [ \t]+ (?P<formulation>[^\n]*?) [ \t]+ formulation [ \t]* \n)?
(?:[ \t]*\n)*
^[ \t]* State [^\n]* \n
^[ \t]* number [^\n]* \n
^[ \t]* -+ [ \t]* \n
(?:
  ^[ \t]* TDDFPT\| [ \t]+
  (?P<state>\d+) [ \t]+
  (?P<energy>{FLOAT}) [ \t]+
  (?P<dx>{FLOAT}) [ \t]+ (?P<dy>{FLOAT}) [ \t]+ (?P<dz>{FLOAT}) [ \t]+
  (?P<oscillator>{FLOAT}) [ \t]* \n
)+
""",
    _FLAGS,
)

EXCITATION_ANALYSIS_MARKER = "Excitation analysis"
_STATE_LINE_RE = re.compile(rf"^[ \t]*(\d+)[ \t]+({FLOAT})[ \t]+eV[ \t]*$")
_TRANSITION_LINE_RE = re.compile(rf"^[ \t]*(\d+)[ \t]+(\d+)[ \t]+({FLOAT})[ \t]*$")


def normalize_formulation(text):
    """'velocity', 'length', 'berry phase', or None if absent/unrecognised."""
    if not text:
        return None
    lowered = text.strip().lower()
    if "velocity" in lowered:
        return "velocity"
    if "berry" in lowered:
        return "berry phase"
    if "length" in lowered:
        return "length"
    return None


def _parse_excitation_analysis(text):
    """Per-state (occupied, virtual, amplitude) lists, or {} if not found.

    Best effort only: any surprise just ends the parse, and nothing here may
    stop the rest of the file from being read.
    """
    start = text.find(EXCITATION_ANALYSIS_MARKER)
    if start < 0:
        return {}
    states = {}
    current = None
    for line in text[start:].splitlines()[1:]:
        stripped = line.strip()
        if stripped.startswith("---"):
            if states:  # the closing rule of the table
                break
            continue
        m = _STATE_LINE_RE.match(line)
        if m:
            current = {"energy": float(m.group(2)), "transitions": []}
            states[int(m.group(1))] = current
            continue
        m = _TRANSITION_LINE_RE.match(line)
        if m and current is not None:
            current["transitions"].append(
                {"occupied": int(m.group(1)), "virtual": int(m.group(2)), "amplitude": float(m.group(3))}
            )
    return states


def _block_from_match(match, analysis_text):
    captures = match.capturesdict()
    analysis = _parse_excitation_analysis(analysis_text)
    states = []
    for i in range(len(captures["state"])):
        number = int(captures["state"][i])
        energy = float(captures["energy"][i])
        state = {
            "number": number,
            "energy": energy,  # eV
            "transition dipole": [
                float(captures["dx"][i]),
                float(captures["dy"][i]),
                float(captures["dz"][i]),
            ],  # a.u.
            "oscillator strength": float(captures["oscillator"][i]),
        }
        entry = analysis.get(number)
        if entry is not None and abs(entry["energy"] - energy) < 5e-4:
            state["transitions"] = entry["transitions"]
        states.append(state)
    return {
        "spin": match["spin"],  # "R" or "U"
        "multiplicity": int(match["multiplicity"]) if match["multiplicity"] else None,
        "formulation": normalize_formulation(match["formulation"]),
        "states": states,
    }


def match_tddfpt(content: str) -> Optional[BlockMatch]:
    """All blocks of the LAST complete TDDFPT calculation in the file.

    A calculation prints one block per spin/multiplicity (for example a singlet
    and a triplet block). A block whose (spin, multiplicity) was already seen
    in the current group belongs to a later step, and starts a new group.
    """
    matches = list(TDDFPT_BLOCK_RE.finditer(content))
    if not matches:
        return None

    groups = []  # each: list of (match, block)
    for index, match in enumerate(matches):
        end = matches[index + 1].start() if index + 1 < len(matches) else len(content)
        block = _block_from_match(match, content[match.end() : end])
        key = (block["spin"], block["multiplicity"])
        if not groups or key in {(b["spin"], b["multiplicity"]) for _, b in groups[-1]}:
            groups.append([])
        groups[-1].append((match, block))

    last = groups[-1]
    return BlockMatch(
        {"tddfpt": {"blocks": [block for _, block in last]}},
        [m.span() for m, _ in last],
    )


def check_oscillator_strengths(block, warn, rel_tol=0.05, min_strength=1e-3):
    """Recompute f from the printed energies and transition dipoles.

    CP2K's printed "transition dipole" is the dipole <r> for the length (and
    Berry phase) formulation but the momentum <p> for the velocity formulation:

        length / Berry phase:  f = (2/3) * dE * |mu|^2
        velocity:              f = 2 * |p|^2 / (3 * dE)

    with dE in hartree. States with f > ``min_strength`` that disagree by more
    than ``rel_tol`` are reported through ``warn``; the printed value is never
    altered. Returns the number of disagreeing states, or None if the
    formulation is unknown and nothing was checked.
    """
    formulation = block["formulation"]
    if formulation not in ("length", "berry phase", "velocity"):
        return None

    bad = 0
    for state in block["states"]:
        printed = state["oscillator strength"]
        if not printed > min_strength:
            continue
        delta_e = state["energy"] / HARTREE_TO_EV
        square = sum(c * c for c in state["transition dipole"])
        if delta_e <= 0.0:
            continue
        if formulation == "velocity":
            computed = 2.0 * square / (3.0 * delta_e)
        else:
            computed = (2.0 / 3.0) * delta_e * square
        if not math.isfinite(computed) or abs(computed - printed) / printed > rel_tol:
            bad += 1
            warn(
                f"TDDFPT state {state['number']} ({formulation} form): printed "
                f"oscillator strength {printed:.5g} but recomputed from the "
                f"printed energy and transition dipole it is {computed:.5g}"
            )
    return bad
