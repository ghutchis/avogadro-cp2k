#  This source file is part of the Avogadro project.
#  This source code is released under the 3-Clause BSD License, (see "LICENSE").

import gzip
import os
import sys

import pytest

HERE = os.path.dirname(__file__)
sys.path.insert(0, os.path.join(HERE, "..", "src"))

DATA = os.path.join(HERE, "data")
REAL_GZ = os.path.join(DATA, "PBE-SVP.out.gz")


def snippet(name):
    """A hand-trimmed excerpt of PBE-SVP.out (lines only dropped, never edited)."""
    with open(os.path.join(DATA, "snippets", name), encoding="utf-8") as handle:
        return handle.read()


@pytest.fixture(scope="session")
def real_text():
    with gzip.open(REAL_GZ, "rb") as handle:
        return handle.read().decode("utf-8", errors="replace")
