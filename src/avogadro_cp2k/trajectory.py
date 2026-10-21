#  This source file is part of the Avogadro project.
#  This source code is released under the 3-Clause BSD License, (see "LICENSE").

"""Reader for the ``<project>-pos-1.xyz`` trajectory CP2K writes beside the
main output of GEO_OPT, CELL_OPT and MD runs."""

import math
import os

import regex as re

from .utils import atomic_number_from_label, decode

# "i =        1, time =        0.000, E =      -233.7225557240"
_ENERGY_RE = re.compile(r"(?:^|[\s,])E\s*=\s*([+\-]?(?:\d+\.?\d*|\.\d+)(?:[Ee][+\-]?\d+)?)(?:\s|,|$)")
_STEP_RE = re.compile(r"(?:^|[\s,])i\s*=\s*(\d+)")


def sibling_path(out_filename, project):
    """Where CP2K would have written the position trajectory."""
    if not out_filename or not project:
        return None
    return os.path.join(os.path.dirname(os.path.abspath(out_filename)), f"{project}-pos-1.xyz")


def parse_xyz_trajectory(text, warn=lambda message: None):
    """Parse a multi-frame XYZ file.

    Returns a list of frames ``{"labels", "numbers", "positions", "energy",
    "step"}``; positions are flat lists in the file's units (angstrom by
    default in CP2K) and energy is hartree or None. A truncated final frame is
    dropped with a warning.
    """
    lines = text.splitlines()
    frames = []
    i = 0
    n = len(lines)
    while i < n:
        if not lines[i].strip():
            i += 1
            continue
        try:
            count = int(lines[i].split()[0])
        except (ValueError, IndexError):
            warn(f"unparseable atom-count line {i + 1} in trajectory; stopping there")
            break
        if i + 2 + count > n:
            warn("the trajectory ends in the middle of a frame; ignoring the partial frame")
            break
        comment = lines[i + 1]
        labels, numbers, positions = [], [], []
        ok = True
        for line in lines[i + 2 : i + 2 + count]:
            parts = line.split()
            if len(parts) < 4:
                ok = False
                break
            try:
                xyz = [float(parts[1]), float(parts[2]), float(parts[3])]
            except ValueError:
                ok = False
                break
            if not all(math.isfinite(c) for c in xyz):
                ok = False
                break
            labels.append(parts[0])
            numbers.append(atomic_number_from_label(parts[0]))
            positions.extend(xyz)
        if not ok:
            warn(f"malformed atom line in the frame starting at line {i + 1}; ignoring the rest")
            break
        energy = None
        m = _ENERGY_RE.search(comment)
        if m:
            value = float(m.group(1))
            if math.isfinite(value):
                energy = value
        step = _STEP_RE.search(comment)
        frames.append(
            {
                "labels": labels,
                "numbers": numbers,
                "positions": positions,
                "energy": energy,
                "step": int(step.group(1)) if step else None,
            }
        )
        i += 2 + count
    return frames


def load_trajectory(path, expected_numbers, warn):
    """Frames from ``path`` whose atoms match ``expected_numbers``, else None.

    Matching is by atom count and by element sequence, so a trajectory that
    belongs to a different system (or a stale file from an earlier run in the
    same directory) is refused instead of silently attached.
    """
    if not path or not os.path.isfile(path):
        return None
    try:
        with open(path, "rb") as handle:
            text = decode(handle.read())
    except OSError as exc:
        warn(f"could not read {path}: {exc}")
        return None
    frames = parse_xyz_trajectory(text, warn)
    if not frames:
        warn(f"{path} holds no complete frames")
        return None
    for frame in frames:
        if frame["numbers"] != list(expected_numbers):
            warn(
                f"{path} does not match the output file's atoms "
                f"({len(frame['numbers'])} atoms vs {len(expected_numbers)}, or different elements)"
            )
            return None
    return frames
