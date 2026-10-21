#  This source file is part of the Avogadro project.
#  This source code is released under the 3-Clause BSD License, (see "LICENSE").

"""Small helpers: decompression, element symbols, non-finite sanitising."""

import bz2
import gzip
import lzma
import math
import sys

_GZIP_MAGIC = b"\x1f\x8b"
_BZIP2_MAGIC = b"BZh"
_XZ_MAGIC = b"\xfd7zXZ\x00"

# CODATA 2018. CP2K prints its own (CODATA 2006) constants in the output
# header; the difference is ~1e-7 relative, far below the printed precision.
BOHR_TO_ANGSTROM = 0.529177210903
HARTREE_TO_EV = 27.211386245988

_SYMBOLS = (
    "H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni "
    "Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I "
    "Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt "
    "Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr "
    "Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og"
).split()
ATOMIC_NUMBER = {symbol: number for number, symbol in enumerate(_SYMBOLS, start=1)}


def atomic_number_from_label(label):
    """Atomic number for an element symbol or a CP2K kind-style label.

    CP2K kind names are user-chosen ("C1", "O_water"), so trailing digits and
    anything after an underscore, dash or digit is dropped before lookup.
    Returns 0 if the label is not an element.
    """
    letters = ""
    for char in label:
        if char.isalpha():
            letters += char
        else:
            break
    if not letters:
        return 0
    for candidate in (letters, letters[:2], letters[:1]):
        number = ATOMIC_NUMBER.get(candidate.capitalize())
        if number:
            return number
    return 0


def decompress(raw):
    """Decompress gzip/bzip2/xz content, sniffed by magic bytes."""
    if raw.startswith(_GZIP_MAGIC):
        return gzip.decompress(raw)
    if raw.startswith(_BZIP2_MAGIC):
        return bz2.decompress(raw)
    if raw.startswith(_XZ_MAGIC):
        return lzma.decompress(raw)
    return raw


def decode(raw):
    """Decode bytes as UTF-8 (BOM tolerated); never consult the locale."""
    return raw.decode("utf-8-sig", errors="replace")


def _sanitize_value(value):
    """Return (native value, saw_non_finite) for a leaf or nested list."""
    if isinstance(value, (list, tuple)):
        converted = []
        dirty = False
        for item in value:
            item_value, item_dirty = _sanitize_value(item)
            dirty = dirty or item_dirty
            converted.append(item_value)
        return converted, dirty
    if hasattr(value, "item"):  # numpy scalar
        value = value.item()
    if isinstance(value, float) and not math.isfinite(value):
        return value, True
    return value, False


def sanitize(obj, path=""):
    """Make ``obj`` safe for ``json.dumps``.

    A single bare NaN/Infinity in a CJSON document makes Avogadro discard the
    whole document. Any key whose value (however deeply nested) holds a
    non-finite float is therefore DROPPED, with a note on stderr; zeros or any
    other placeholder are never substituted.
    """
    sanitized = {}
    for key, value in obj.items():
        if isinstance(value, dict):
            sanitized[key] = sanitize(value, path + key + ".")
            continue
        converted, dirty = _sanitize_value(value)
        if dirty:
            print(
                f"avogadro-cp2k: dropping non-finite values found in '{path}{key}'",
                file=sys.stderr,
            )
            continue
        sanitized[key] = converted
    return sanitized
