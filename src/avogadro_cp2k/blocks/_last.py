#  This source file is part of the Avogadro project.
#  This source code is released under the 3-Clause BSD License, (see "LICENSE").


def last_match(pattern, content):
    """The last match of a compiled regex in ``content`` (or None).

    CP2K repeats most blocks once per step; upstream's ``.search()`` returns
    step 0 of a GEO_OPT or MD run, whereas the final step is what we want.
    """
    found = None
    for found in pattern.finditer(content):
        pass
    return found
