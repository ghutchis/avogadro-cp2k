# Vendored code: cp2k-output-tools

Upstream: https://github.com/cp2k/cp2k-output-tools (branch `develop`)
Commit:   `403cc70966d7284d9f71f56b8cd2ae08f969fa62`
License:  MIT, Copyright (c) 2022 Tiziano Muller (see `LICENSE` in this directory)

Only the files below were copied. Upstream's other block modules and the
`cli.py` / `levelparser.py` / `trajectories` / `scripts` code were not, because
they pull in `pint`, `click`, `rich` and friends, which this plugin does not
depend on.

| File | Status |
|------|--------|
| `LICENSE` | verbatim |
| `parser.py` | verbatim (imports only `.blocks`, `.blocks.common`) |
| `blocks/common.py` | verbatim |
| `blocks/mulliken.py` | modified, see below |
| `__init__.py` | replaced by a comment-only file (upstream imports `.blocks` and `.parser` and exports `__version__`) |
| `blocks/__init__.py` | replaced: exports only `builtin_matchers = [match_mulliken_population_analysis]` (upstream imports every block module, and through them `pint`) |

## Local changes

Each is marked in the source with a `# avogadro-cp2k:` comment.

1. `blocks/mulliken.py`, `match_mulliken_population_analysis`: replaced
   `MULLIKEN_POPULATION_ANALYSIS_RE.search(content)` (first occurrence, which
   is step 0 of a GEO_OPT/MD run) with a `finditer` loop that keeps the LAST
   complete match.

No other upstream line was changed. The regular expression itself is untouched.
