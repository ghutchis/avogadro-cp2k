#  This source file is part of the Avogadro project.
#  This source code is released under the 3-Clause BSD License, (see "LICENSE").

"""Turn CP2K main-output text into a Chemical JSON (CJSON) dictionary."""

import math
import sys

import numpy as np

from .blocks import builtin_matchers
from .blocks.tddfpt import check_oscillator_strengths
from .trajectory import load_trajectory, sibling_path
from .utils import atomic_number_from_label, sanitize
from .vendored.cp2k_output_tools.parser import parse_iter

# Run types for which CP2K writes a per-step trajectory beside the main output
# instead of reprinting the coordinates there.
TRAJECTORY_RUN_TYPES = {
    "GEO_OPT",
    "GEOMETRY_OPTIMIZATION",
    "CELL_OPT",
    "MD",
    "MOLECULAR_DYNAMICS",
}
CELL_RUN_TYPES = {"CELL_OPT"}

# A recomputed frame energy that differs from the output file's final energy by
# more than this (hartree) suggests the trajectory belongs to another run.
ENERGY_MATCH_TOLERANCE = 1.0e-5

# Printed cell vectors carry three decimals but the lengths and angles carry
# six; the two are merged only if they agree to within rounding.
CELL_CONSISTENCY_TOLERANCE = 2.0e-3


def _default_warn(message):
    print(f"avogadro-cp2k: warning: {message}", file=sys.stderr)


def parse_blocks(content):
    """Run every matcher and merge the resulting dictionaries."""
    parsed = {}
    for block in parse_iter(content, matchers=builtin_matchers):
        parsed.update(block)
    return parsed


def standard_orientation_vectors(lengths, angles):
    """Cell vectors (a along x, b in the xy plane) from lengths and angles."""
    a, b, c = lengths
    alpha, beta, gamma = (math.radians(x) for x in angles)
    cos_a, cos_b, cos_g = math.cos(alpha), math.cos(beta), math.cos(gamma)
    sin_g = math.sin(gamma)
    cx = c * cos_b
    cy = c * (cos_a - cos_b * cos_g) / sin_g
    cz_squared = c * c - cx * cx - cy * cy
    if sin_g == 0.0 or cz_squared <= 0.0:
        return None
    return np.array([[a, 0.0, 0.0], [b * cos_g, b * sin_g, 0.0], [cx, cy, math.sqrt(cz_squared)]])


def cell_vectors(cell):
    """Best cell vectors (angstrom) for the parsed ``CELL|`` block.

    CP2K prints the vectors with three decimals but the lengths and angles with
    six. When the vectors are exactly the standard orientation of those lengths
    and angles (to within the vectors' rounding), the more precise rebuilt
    vectors are used. Otherwise (a rotated cell, say) the printed vectors are
    kept, since the atoms are in that frame.
    """
    printed = np.array(cell["vectors"], dtype=float)
    if "angles" in cell:
        rebuilt = standard_orientation_vectors(cell["lengths"], cell["angles"])
        if rebuilt is not None and np.max(np.abs(rebuilt - printed)) <= CELL_CONSISTENCY_TOLERANCE:
            return rebuilt
    return printed


def _atomic_numbers(atoms, warn):
    numbers = []
    for index, atom in enumerate(atoms, start=1):
        number = atom["atomic number"]
        if number == 0:
            # Ghost/dummy atoms print Z = 0; fall back on the element symbol.
            number = atomic_number_from_label(atom["element"])
            if number == 0:
                warn(f"atom {index} ({atom['element']}) has no atomic number; stored as a dummy atom")
        numbers.append(number)
    return numbers


def _electronic_spectrum(blocks, warn):
    """energies (eV) and oscillator strengths of the block that is the spectrum.

    Restricted runs: the singlet block only (restricted triplets are
    spin-forbidden and print f = 0). Unrestricted runs: the U block.
    """
    chosen = None
    for block in blocks:
        if block["spin"] == "R" and block["multiplicity"] == 1:
            chosen = block
            break
        if block["spin"] == "U":
            chosen = block
            break
    if chosen is None or not chosen["states"]:
        return None
    return {
        "energies": [s["energy"] for s in chosen["states"]],
        "intensities": [s["oscillator strength"] for s in chosen["states"]],
    }


def read_cp2k(content, filename=None, warn=_default_warn, debug=False):
    """Parse the text of a CP2K main output into a CJSON dictionary.

    ``filename`` (the real path of the output, when known) lets GEO_OPT, CELL_OPT
    and MD runs pick up the sibling ``<project>-pos-1.xyz`` trajectory.
    """
    parsed = parse_blocks(content)
    run = parsed.get("run info", {})

    block = parsed.get("atomic coordinates")
    if not block or not block["atoms"]:
        raise ValueError(
            "no 'ATOMIC COORDINATES' table found; this does not look like a CP2K main output, "
            "or the run stopped before the Quickstep setup finished"
        )
    atoms = block["atoms"]

    expected = parsed.get("atom count")
    if expected is not None and expected != len(atoms):
        raise ValueError(
            f"the coordinate table has {len(atoms)} atoms but the output reports {expected}; "
            "refusing to build a partial molecule"
        )

    numbers = _atomic_numbers(atoms, warn)
    positions = [c for atom in atoms for c in atom["position"]]
    natoms = len(atoms)

    cjson = {"chemicalJson": 1}
    if run.get("project"):
        cjson["name"] = run["project"]

    coords = {}
    properties = {}
    energy = parsed.get("total energy", {}).get("hartree")

    # ---- geometry: single structure, or the sibling trajectory ---------------
    run_type = (run.get("run type") or "").upper()
    frames = None
    if run_type in TRAJECTORY_RUN_TYPES:
        path = sibling_path(filename, run.get("project"))
        frames = load_trajectory(path, numbers, warn)
        if frames:
            positions = frames[-1]["positions"]
            if len(frames) > 1:
                coords["3dSets"] = [f["positions"] for f in frames]
                coords["3dSetsActive"] = len(frames) - 1  # show the final structure
                frame_energies = [f["energy"] for f in frames]
                if all(e is not None for e in frame_energies):
                    properties["energies"] = frame_energies
                    properties["energyUnit"] = "Hartree"
            last_energy = frames[-1]["energy"]
            if last_energy is not None and energy is not None and abs(last_energy - energy) > ENERGY_MATCH_TOLERANCE:
                warn(
                    f"the last trajectory frame's energy ({last_energy:.8f} Eh) differs from the "
                    f"output's final energy ({energy:.8f} Eh); {path} may belong to a different run"
                )
        else:
            where = path if path else "(input arrived on stdin, so no path is known)"
            warn(
                f"run type {run_type}: no usable {run.get('project', '<project>')}-pos-1.xyz was found "
                f"next to the output {where}. The geometry shown is the INITIAL structure from the "
                "output file, while the energy, charges and spectrum are from the FINAL step."
            )
    coords["3d"] = positions
    cjson["atoms"] = {"elements": {"number": numbers}, "coords": coords}

    # ---- cell ----------------------------------------------------------------
    cell = parsed.get("cell")
    if cell and (cell["periodicity"] or "").upper() != "NONE":
        vectors = cell_vectors(cell)
        cjson["unitCell"] = {"cellVectors": [float(x) for x in vectors.flatten()]}
        if run_type in CELL_RUN_TYPES:
            warn(
                "CELL_OPT: the unit cell is the last CELL| block in the output, which is normally "
                "the INITIAL cell; the optimized cell is in <project>-1.cell and is not read"
            )

    # ---- energy and run properties --------------------------------------------
    if energy is not None:
        properties["totalEnergy"] = energy
        properties["energyUnit"] = "Hartree"
    if "charge" in run:
        properties["totalCharge"] = run["charge"]
    if "multiplicity" in run:
        properties["totalSpinMultiplicity"] = run["multiplicity"]
    if properties:
        cjson["properties"] = properties

    # ---- partial charges ------------------------------------------------------
    partial = {}
    for key, name in (("mulliken population analysis", "Mulliken"), ("hirshfeld charges", "Hirshfeld")):
        entry = parsed.get(key)
        if not entry:
            continue
        charges = [a["charge"] for a in entry["per atom"]]
        if len(charges) != natoms:
            warn(f"{name} analysis lists {len(charges)} atoms but the molecule has {natoms}; ignored")
            continue
        partial[name] = charges
    if partial:
        cjson["atoms"]["partialCharges"] = partial

    # ---- excited states -------------------------------------------------------
    tddfpt = parsed.get("tddfpt")
    if tddfpt:
        for tblock in tddfpt["blocks"]:
            check_oscillator_strengths(
                tblock,
                lambda message, tblock=tblock: warn(
                    f"{tblock['spin']}-TDDFPT multiplicity {tblock['multiplicity']}: {message}"
                ),
            )
        spectrum = _electronic_spectrum(tddfpt["blocks"], warn)
        if spectrum:
            cjson["spectra"] = {"electronic": spectrum}

    if debug:
        print(
            f"avogadro-cp2k: version={run.get('version')!r} run_type={run_type or None!r} "
            f"atoms={natoms} cell={'unitCell' in cjson} charges={sorted(partial)} "
            f"states={len(cjson.get('spectra', {}).get('electronic', {}).get('energies', []))} "
            f"frames={len(frames) if frames else 1}",
            file=sys.stderr,
        )

    cjson = sanitize(cjson)
    if "coords" not in cjson.get("atoms", {}) or "3d" not in cjson["atoms"]["coords"]:
        raise ValueError("the atomic coordinates contain non-finite values")
    return cjson
