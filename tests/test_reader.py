#  This source file is part of the Avogadro project.
#  This source code is released under the 3-Clause BSD License, (see "LICENSE").

"""End-to-end tests on the complete PBE-SVP.out (CP2K 2026.2, RKS, ENERGY)."""

import bz2
import json
import lzma
import math
import os
import subprocess
import sys

import pytest
from conftest import HERE, REAL_GZ, snippet

from avogadro_cp2k.reader import read_cp2k, standard_orientation_vectors
from avogadro_cp2k.trajectory import parse_xyz_trajectory
from avogadro_cp2k.utils import atomic_number_from_label, decode, decompress, sanitize

SRC = os.path.join(HERE, "..", "src")


@pytest.fixture(scope="module")
def result(real_text):
    warnings = []
    cjson = read_cp2k(real_text, warn=warnings.append)
    return cjson, warnings


def test_no_warnings(result):
    assert result[1] == []  # includes the oscillator-strength validation


def test_chemical_json_header(result):
    cjson = result[0]
    assert cjson["chemicalJson"] == 1
    assert cjson["name"] == "h2o"


def test_atoms_and_elements(result):
    atoms = result[0]["atoms"]
    numbers = atoms["elements"]["number"]
    assert len(numbers) == 46
    assert numbers.count(16) == 4  # S
    assert numbers.count(7) == 4  # N
    assert numbers.count(6) == 26  # C
    assert numbers.count(1) == 12  # H
    assert numbers[:4] == [16] * 4 and numbers[4:8] == [7] * 4 and numbers[-12:] == [1] * 12


def test_coordinates(result):
    coords = result[0]["atoms"]["coords"]["3d"]
    assert len(coords) == 3 * 46
    assert coords[:3] == [1.086966, 1.304941, 1.378422]
    assert coords[-3:] == [5.170507, 1.120944, 1.901010]
    # atom 11 lies outside the cell and is left exactly as printed
    assert coords[30:33] == [-1.985033, 0.665964, 5.452731]
    assert "3dSets" not in result[0]["atoms"]["coords"]


def test_unit_cell(result):
    vectors = result[0]["unitCell"]["cellVectors"]
    printed = [9.215, 0.0, 0.0, 4.046, 9.845, 0.0, -4.130, -1.613, 6.337]
    assert len(vectors) == 9
    for got, want in zip(vectors, printed):
        assert got == pytest.approx(want, abs=1e-3)  # printed with 3 decimals
    # the lengths and angles are printed to six decimals and the rebuilt
    # vectors honour them exactly
    a, b, c = (vectors[0:3], vectors[3:6], vectors[6:9])
    norm = lambda v: math.sqrt(sum(x * x for x in v))  # noqa: E731
    assert norm(a) == pytest.approx(9.215, abs=1e-9)
    assert norm(b) == pytest.approx(10.644, abs=1e-9)
    assert norm(c) == pytest.approx(7.734, abs=1e-9)
    dot = lambda u, v: sum(x * y for x, y in zip(u, v))  # noqa: E731
    gamma = math.degrees(math.acos(dot(a, b) / (norm(a) * norm(b))))
    assert gamma == pytest.approx(67.66, abs=1e-9)
    # volume check against the printed 574.882486 A^3
    volume = abs(
        a[0] * (b[1] * c[2] - b[2] * c[1])
        - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
    )
    assert volume == pytest.approx(574.882486, rel=1e-5)


def test_partial_charges(result):
    charges = result[0]["atoms"]["partialCharges"]
    assert set(charges) == {"Hirshfeld", "Mulliken"}
    assert len(charges["Hirshfeld"]) == 46
    assert len(charges["Mulliken"]) == 46
    assert charges["Hirshfeld"][0] == -0.110 and charges["Hirshfeld"][-1] == 0.424
    assert charges["Mulliken"][0] == 0.453840 and charges["Mulliken"][-1] == 0.204724
    assert sum(charges["Mulliken"]) == pytest.approx(0.0, abs=1e-5)  # neutral
    assert sum(charges["Hirshfeld"]) == pytest.approx(0.0, abs=2e-2)  # 3-decimal rounding


def test_electronic_spectrum(result):
    spectrum = result[0]["spectra"]["electronic"]
    assert spectrum["energies"] == [1.51072, 2.40916, 2.43248, 2.54017, 2.64301]
    assert spectrum["intensities"] == [
        1.63771e-01,
        3.88607e-02,
        5.68049e-02,
        6.77036e-03,
        2.29763e-02,
    ]


def test_energy_and_properties(result):
    props = result[0]["properties"]
    assert props["totalEnergy"] == -233.72255572397043
    assert props["energyUnit"] == "Hartree"
    assert props["totalCharge"] == 0
    assert props["totalSpinMultiplicity"] == 1


def test_valid_json_without_nan(result):
    text = json.dumps(result[0], allow_nan=False)  # raises on NaN/Infinity
    assert "NaN" not in text and "Infinity" not in text


# --- the command-line entry point ----------------------------------------------


def run_cli(stdin_bytes, *extra):
    code = "import sys; from avogadro_cp2k import main; main()"
    env = dict(os.environ, PYTHONPATH=SRC)
    return subprocess.run(
        [sys.executable, "-c", code, "cp2k", "--read", *extra],
        input=stdin_bytes,
        capture_output=True,
        env=env,
        check=False,
    )


def check_cli_output(proc):
    assert proc.returncode == 0, proc.stderr.decode()
    assert proc.stderr.decode() == ""
    cjson = json.loads(proc.stdout.decode())  # bare CJSON: no envelope or banner
    assert len(cjson["atoms"]["elements"]["number"]) == 46
    assert len(cjson["atoms"]["partialCharges"]["Hirshfeld"]) == 46
    assert len(cjson["spectra"]["electronic"]["energies"]) == 5
    assert cjson["properties"]["totalEnergy"] == -233.72255572397043


def test_cli_file_mode_gzip():
    payload = json.dumps({"operation": "read", "filename": REAL_GZ}).encode()
    check_cli_output(run_cli(payload))


def test_cli_file_mode_plain_text(tmp_path, real_text):
    path = tmp_path / "PBE-SVP.out"
    path.write_bytes(real_text.encode("utf-8"))
    payload = json.dumps({"operation": "read", "filename": str(path)}).encode()
    check_cli_output(run_cli(payload))


def test_cli_string_mode_raw_and_compressed(real_text):
    raw = real_text.encode("utf-8")
    check_cli_output(run_cli(raw))
    with open(REAL_GZ, "rb") as handle:
        check_cli_output(run_cli(handle.read()))
    check_cli_output(run_cli(bz2.compress(raw)))
    check_cli_output(run_cli(lzma.compress(raw)))


def test_cli_rejects_non_cp2k_input():
    proc = run_cli(b"this is not a CP2K output\n")
    assert proc.returncode == 1
    assert proc.stdout == b""
    assert b"ATOMIC COORDINATES" in proc.stderr


def test_cli_missing_file_is_an_error():
    proc = run_cli(json.dumps({"operation": "read", "filename": "/no/such/file.out"}).encode())
    assert proc.returncode == 1 and proc.stdout == b""


# --- guards ----------------------------------------------------------------------


def test_atom_count_mismatch_refuses_partial_molecule():
    text = snippet("coordinates.txt") + snippet("atom_count.txt")  # 7 rows vs "46"
    with pytest.raises(ValueError, match="46"):
        read_cp2k(text)


def test_periodicity_none_emits_no_cell(real_text):
    text = real_text.replace("CELL| Periodicity                                                           XYZ",
                             "CELL| Periodicity                                                          NONE")
    assert text != real_text
    assert "unitCell" not in read_cp2k(text, warn=lambda m: None)


def test_oscillator_disagreement_warns_but_keeps_value(real_text):
    text = real_text.replace("1.63771E-01", "2.63771E-01")
    warnings = []
    cjson = read_cp2k(text, warn=warnings.append)
    assert any("state 1" in w for w in warnings)
    assert cjson["spectra"]["electronic"]["intensities"][0] == 2.63771e-01


def test_sanitize_drops_keys_with_non_finite_values():
    cleaned = sanitize(
        {"a": [1.0, float("nan")], "b": {"c": float("inf"), "d": 2.0}, "e": [[1.0], [2.0]]}
    )
    assert cleaned == {"b": {"d": 2.0}, "e": [[1.0], [2.0]]}


def test_decompress_and_decode():
    data = "café".encode("utf-8")
    assert decompress(bz2.compress(data)) == data
    assert decompress(lzma.compress(data)) == data
    assert decompress(data) == data
    assert decode(b"\xef\xbb\xbfabc") == "abc"  # BOM stripped
    assert decode(b"ab\xffcd") == "ab�cd"  # never raises


def test_atomic_number_from_label():
    assert atomic_number_from_label("S") == 16
    assert atomic_number_from_label("C1") == 6
    assert atomic_number_from_label("Cl") == 17
    assert atomic_number_from_label("X") == 0


def test_standard_orientation_matches_cp2k_example():
    vectors = standard_orientation_vectors([9.215, 10.644, 7.734], [113.32, 122.28, 67.66])
    assert vectors[0].tolist() == [9.215, 0.0, 0.0]
    assert vectors[1][2] == 0.0
    assert vectors[2][2] > 0


# --- sibling-trajectory code paths ------------------------------------------------
# These exercise the logic with a SYNTHETIC trajectory derived from the real
# coordinates (the real file is an ENERGY run, so its run type is edited). They
# check the plumbing, not CP2K's real <project>-pos-1.xyz layout.


def _synthetic_trajectory(real_text, frames=3):
    from avogadro_cp2k.blocks.geometry import match_atomic_coordinates

    atoms = match_atomic_coordinates(real_text).data["atomic coordinates"]["atoms"]
    out = []
    for step in range(frames):
        out.append(str(len(atoms)))
        out.append(f" i = {step:8d}, time = {0.5 * step:12.3f}, E = {-233.0 - 0.1 * step:20.10f}")
        for atom in atoms:
            x, y, z = (c + 0.01 * step for c in atom["position"])
            out.append(f" {atom['element']:2s} {x:16.10f} {y:16.10f} {z:16.10f}")
    return "\n".join(out) + "\n"


def test_parse_xyz_trajectory_drops_a_truncated_last_frame(real_text):
    text = _synthetic_trajectory(real_text, 2)
    messages = []
    frames = parse_xyz_trajectory(text + "46\n i = 2\n S 0 0 0\n", messages.append)
    assert len(frames) == 2
    assert frames[1]["energy"] == pytest.approx(-233.1)
    assert frames[1]["step"] == 1
    assert messages and "partial frame" in messages[0]


def test_geo_opt_uses_sibling_trajectory(tmp_path, real_text):
    out = tmp_path / "run.out"
    out.write_text(real_text.replace("ENERGY\n", "GEO_OPT\n", 1))
    (tmp_path / "h2o-pos-1.xyz").write_text(_synthetic_trajectory(real_text, 3))
    warnings = []
    cjson = read_cp2k(out.read_text(), filename=str(out), warn=warnings.append)
    coords = cjson["atoms"]["coords"]
    assert len(coords["3dSets"]) == 3
    assert coords["3dSetsActive"] == 2
    assert coords["3d"] == coords["3dSets"][-1]
    assert coords["3d"][0] == pytest.approx(1.086966 + 0.02)
    assert cjson["properties"]["energies"] == pytest.approx([-233.0, -233.1, -233.2])
    # the synthetic last-frame energy deliberately differs from the .out energy
    assert any("different run" in w for w in warnings)


def test_geo_opt_without_sibling_warns_loudly(tmp_path, real_text):
    out = tmp_path / "run.out"
    out.write_text(real_text.replace("ENERGY\n", "GEO_OPT\n", 1))
    warnings = []
    cjson = read_cp2k(out.read_text(), filename=str(out), warn=warnings.append)
    assert len(warnings) == 1 and "INITIAL" in warnings[0] and "FINAL" in warnings[0]
    assert "3dSets" not in cjson["atoms"]["coords"]


def test_geo_opt_with_mismatched_sibling_is_refused(tmp_path, real_text):
    out = tmp_path / "run.out"
    out.write_text(real_text.replace("ENERGY\n", "GEO_OPT\n", 1))
    (tmp_path / "h2o-pos-1.xyz").write_text("1\n i = 0, E = -1.0\n He 0.0 0.0 0.0\n")
    warnings = []
    cjson = read_cp2k(out.read_text(), filename=str(out), warn=warnings.append)
    assert any("does not match" in w for w in warnings)
    assert any("INITIAL" in w for w in warnings)
    assert "3dSets" not in cjson["atoms"]["coords"]


NEED_FILE = "needs a real CP2K output from Geoff; will not be fabricated"


@pytest.mark.skip(reason=f"real GEO_OPT .out plus its -pos-1.xyz (format, last frame): {NEED_FILE}")
def test_geo_opt_real_sibling_pos_xyz():
    pass
