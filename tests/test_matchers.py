#  This source file is part of the Avogadro project.
#  This source code is released under the 3-Clause BSD License, (see "LICENSE").

"""Per-matcher tests on excerpts cut from the real PBE-SVP.out."""

import glob
import os

import pytest
from conftest import DATA, snippet

from avogadro_cp2k.blocks.cell import match_cell
from avogadro_cp2k.blocks.energy import match_total_energy
from avogadro_cp2k.blocks.geometry import match_atom_count, match_atomic_coordinates
from avogadro_cp2k.blocks.hirshfeld import match_hirshfeld_charges
from avogadro_cp2k.blocks.run_info import match_run_info
from avogadro_cp2k.blocks.tddfpt import (
    check_oscillator_strengths,
    match_tddfpt,
    normalize_formulation,
)
from avogadro_cp2k.vendored.cp2k_output_tools.blocks.mulliken import (
    match_mulliken_population_analysis,
)


def test_snippets_are_verbatim_lines_of_the_real_file(real_text):
    """Trimming only: every snippet line must occur in PBE-SVP.out."""
    real_lines = {line.rstrip() for line in real_text.splitlines()}
    paths = glob.glob(os.path.join(DATA, "snippets", "*.txt"))
    assert paths
    for path in paths:
        with open(path, encoding="utf-8") as handle:
            for line in handle.read().splitlines():
                assert line.rstrip() in real_lines, f"{os.path.basename(path)}: {line!r}"


# --- run info ---------------------------------------------------------------


def test_run_info():
    data = match_run_info(snippet("run_info.txt")).data["run info"]
    assert data["version"] == "CP2K version 2026.2"
    assert data["run type"] == "ENERGY"
    assert data["project"] == "h2o"
    assert data["spin treatment"] == "RKS"
    assert data["multiplicity"] == 1
    assert data["charge"] == 0


def test_run_info_absent():
    assert match_run_info("nothing to see here\n") is None


# --- cell -------------------------------------------------------------------


def test_cell_uses_cell_lines_only():
    match = match_cell(snippet("cell.txt"))
    cell = match.data["cell"]
    assert cell["vectors"] == [[9.215, 0.0, 0.0], [4.046, 9.845, 0.0], [-4.130, -1.613, 6.337]]
    assert cell["lengths"] == [9.215, 10.644, 7.734]
    assert cell["angles"] == [113.32, 122.28, 67.66]
    assert cell["periodicity"] == "XYZ"
    # the span covers the CELL| run, not the CELL_REF| run that follows
    start, end = match.spans
    assert "CELL_REF" not in snippet("cell.txt")[start:end]


def test_cell_ref_alone_is_not_a_cell():
    ref_only = "\n".join(line for line in snippet("cell.txt").splitlines() if "CELL_REF|" in line)
    assert match_cell(ref_only) is None


# --- geometry ---------------------------------------------------------------


def test_atomic_coordinates():
    block = match_atomic_coordinates(snippet("coordinates.txt")).data["atomic coordinates"]
    assert block["unit"] == "angstrom"
    atoms = block["atoms"]
    assert len(atoms) == 7
    assert atoms[0]["element"] == "S"
    assert atoms[0]["atomic number"] == 16
    assert atoms[0]["kind"] == 1
    assert atoms[0]["position"] == [1.086966, 1.304941, 1.378422]
    assert atoms[0]["z_eff"] == 6.0
    assert atoms[0]["mass"] == 32.065
    assert atoms[-1]["element"] == "H" and atoms[-1]["atomic number"] == 1
    assert atoms[-1]["position"] == [1.299369, -0.916197, 4.296283]


def test_atomic_coordinates_bohr_are_converted():
    text = snippet("coordinates.txt").replace("IN ANGSTROM", "IN BOHR")
    atoms = match_atomic_coordinates(text).data["atomic coordinates"]["atoms"]
    assert atoms[0]["position"][0] == pytest.approx(1.086966 * 0.529177210903)


def test_atomic_coordinates_unknown_unit_is_refused():
    text = snippet("coordinates.txt").replace("IN ANGSTROM", "IN FURLONGS")
    assert match_atomic_coordinates(text) is None


def test_atom_count():
    assert match_atom_count(snippet("atom_count.txt")).data["atom count"] == 46


# --- energy -----------------------------------------------------------------


def test_total_energy():
    data = match_total_energy(snippet("energy.txt")).data["total energy"]
    assert data["method"] == "QS"
    assert data["hartree"] == -233.72255572397043


# --- Mulliken (vendored matcher) ---------------------------------------------


def test_mulliken_restricted():
    block = match_mulliken_population_analysis(snippet("mulliken.txt")).data["mulliken population analysis"]
    per_atom = block["per atom"]
    assert len(per_atom) == 8  # trimmed excerpt: atoms 1-5 and 44-46
    assert per_atom[0] == {"element": "S", "kind": 1, "population": 5.54616, "charge": 0.45384}
    assert per_atom[-1]["charge"] == 0.204724
    assert block["total"] == {"population": 160.0, "charge": -0.0}


# --- Hirshfeld ---------------------------------------------------------------


def test_hirshfeld_restricted():
    block = match_hirshfeld_charges(snippet("hirshfeld.txt")).data["hirshfeld charges"]
    per_atom = block["per atom"]
    assert len(per_atom) == 7  # trimmed excerpt: atoms 1-5 and 45-46
    assert per_atom[0] == {
        "element": "S",
        "kind": 1,
        "reference charge": 6.0,
        "population": 6.11,
        "charge": -0.11,
    }
    assert per_atom[-1]["charge"] == 0.424
    assert block["total charge"] == 0.0  # printed as -0.000


# --- TDDFPT ------------------------------------------------------------------

EXPECTED_STATES = [
    (1, 1.51072, [1.1066e-01, 3.6585e-02, -7.3801e-03], 1.63771e-01),
    (2, 2.40916, [3.6983e-02, 2.8620e-02, -5.4534e-02], 3.88607e-02),
    (3, 2.43248, [-2.0598e-02, -1.6588e-02, -8.3171e-02], 5.68049e-02),
    (4, 2.54017, [-1.4214e-02, -9.3684e-03, -2.5656e-02], 6.77036e-03),
    (5, 2.64301, [-1.9917e-02, -5.2981e-02, -1.1993e-02], 2.29763e-02),
]


def test_tddfpt_block():
    blocks = match_tddfpt(snippet("tddfpt.txt")).data["tddfpt"]["blocks"]
    assert len(blocks) == 1
    block = blocks[0]
    assert block["spin"] == "R"
    assert block["multiplicity"] == 1
    assert block["formulation"] == "velocity"
    assert len(block["states"]) == 5
    for state, (number, energy, dipole, strength) in zip(block["states"], EXPECTED_STATES):
        assert state["number"] == number
        assert state["energy"] == energy
        assert state["transition dipole"] == dipole
        assert state["oscillator strength"] == strength


def test_tddfpt_excitation_analysis_is_attached():
    states = match_tddfpt(snippet("tddfpt.txt")).data["tddfpt"]["blocks"][0]["states"]
    first = states[0]["transitions"]
    assert first[0] == {"occupied": 80, "virtual": 81, "amplitude": 0.950274}
    assert len(first) == 8
    assert len(states[1]["transitions"]) == 11
    assert "transitions" in states[2]  # the trimmed tail still parses
    assert "transitions" not in states[3]  # not in the excerpt


def test_tddfpt_without_formulation_line():
    text = "\n".join(
        line for line in snippet("tddfpt.txt").splitlines() if "formulation" not in line
    )
    block = match_tddfpt(text).data["tddfpt"]["blocks"][0]
    assert block["formulation"] is None
    assert len(block["states"]) == 5
    assert check_oscillator_strengths(block, lambda message: None) is None  # skipped


def test_normalize_formulation():
    assert normalize_formulation("velocity") == "velocity"
    assert normalize_formulation("length") == "length"
    assert normalize_formulation("Berry phase") == "berry phase"
    assert normalize_formulation("something new") is None
    assert normalize_formulation(None) is None


def test_oscillator_strengths_validate_in_velocity_form():
    block = match_tddfpt(snippet("tddfpt.txt")).data["tddfpt"]["blocks"][0]
    messages = []
    assert check_oscillator_strengths(block, messages.append) == 0
    assert messages == []


def test_oscillator_check_flags_a_wrong_value_without_altering_it():
    block = match_tddfpt(snippet("tddfpt.txt")).data["tddfpt"]["blocks"][0]
    block["states"][0]["oscillator strength"] *= 1.2
    messages = []
    assert check_oscillator_strengths(block, messages.append) == 1
    assert len(messages) == 1 and "state 1" in messages[0]
    assert block["states"][0]["oscillator strength"] == pytest.approx(1.63771e-01 * 1.2)


def test_oscillator_check_distinguishes_forms():
    """The velocity data must NOT pass as length-form data (the two differ)."""
    block = match_tddfpt(snippet("tddfpt.txt")).data["tddfpt"]["blocks"][0]
    block["formulation"] = "length"
    assert check_oscillator_strengths(block, lambda message: None) > 0


# --- last occurrence ---------------------------------------------------------


@pytest.mark.parametrize(
    "matcher, name",
    [
        (match_run_info, "run_info.txt"),
        (match_cell, "cell.txt"),
        (match_atom_count, "atom_count.txt"),
        (match_atomic_coordinates, "coordinates.txt"),
        (match_total_energy, "energy.txt"),
        (match_mulliken_population_analysis, "mulliken.txt"),
        (match_hirshfeld_charges, "hirshfeld.txt"),
        (match_tddfpt, "tddfpt.txt"),
    ],
)
def test_matchers_select_the_last_occurrence(matcher, name):
    """Two copies back to back: the reported span must lie in the second copy.

    (This checks the selection logic only; it cannot check that GEO_OPT step
    values are read correctly. See the skipped tests for that.)
    """
    text = snippet(name)
    twice = text + "\n" + text
    spans = matcher(twice).spans
    starts = [s if isinstance(s, int) else s[0] for s in spans]
    assert min(starts) >= len(text)


# --- needs real files (Geoff) --------------------------------------------------

NEED_FILE = "needs a real CP2K output from Geoff; will not be fabricated"


@pytest.mark.skip(reason=f"unrestricted (UKS) Hirshfeld: {NEED_FILE}")
def test_hirshfeld_unrestricted():
    pass


@pytest.mark.skip(reason=f"unrestricted (UKS) Mulliken, alpha/beta columns: {NEED_FILE}")
def test_mulliken_unrestricted():
    pass


@pytest.mark.skip(reason=f"R-TDDFPT triplet block (multiplicity 3, f = 0): {NEED_FILE}")
def test_tddfpt_triplet_block():
    pass


@pytest.mark.skip(reason=f"U-TDDFPT block layout: {NEED_FILE}")
def test_tddfpt_unrestricted():
    pass


@pytest.mark.skip(reason=f"GEO_OPT output (last-step charges/energy/TDDFPT): {NEED_FILE}")
def test_geo_opt_last_match():
    pass
