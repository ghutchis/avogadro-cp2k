# avogadro-cp2k

Avogadro plugin that reads [CP2K](https://www.cp2k.org/) main output files
(`.out` / `.log` / `.cp2k`, also when gzip-, bzip2- or xz-compressed) and
hands Avogadro Chemical JSON (CJSON). It was developed against CP2K 2026.2.

## What is read

| CP2K output | Avogadro (CJSON key) |
|-------------|----------------------|
| `MODULE QUICKSTEP: ATOMIC COORDINATES` table | `atoms.elements.number`, `atoms.coords.3d` (Angstrom; Bohr is converted) |
| `CELL\|` block, last occurrence | `unitCell.cellVectors` (omitted when `Periodicity` is `NONE`) |
| `ENERGY\| Total FORCE_EVAL ( QS ) energy`, last occurrence | `properties.totalEnergy` (Hartree) and `properties.energyUnit = "Hartree"` |
| `DFT\| Charge`, `DFT\| Multiplicity` | `properties.totalCharge`, `properties.totalSpinMultiplicity` |
| Mulliken population analysis, last occurrence | `atoms.partialCharges.Mulliken` (net charge, e) |
| Hirshfeld charges, last occurrence | `atoms.partialCharges.Hirshfeld` (net charge, e) |
| TDDFPT states (singlet block for restricted runs, the U block for unrestricted) | `spectra.electronic.energies` (eV) and `.intensities` (oscillator strengths) |
| `GLOBAL\| Project name` | `name` |

The element of each atom comes from the atomic-number column of the coordinate
table, so user-chosen kind names (`C1`, `O_water`, ...) are no problem. Bonds are
perceived by Avogadro after reading.

Not stored in CJSON, because Avogadro has no place for it: transition dipoles,
the per-state occupied -> virtual amplitudes of the "Excitation analysis"
table (these are parsed, and available to code that imports
`avogadro_cp2k.blocks.tddfpt`), spin moments and alpha/beta populations of
open-shell runs, and "Natural Orbital analysis".

## Oscillator strengths: velocity versus length form

CP2K prints a "transition dipole" and an oscillator strength per state, and the
meaning of the dipole columns depends on the formulation named on the line
`Transition dipoles calculated using ... formulation`:

* length form (and Berry phase): the columns are the dipole <mu>, and
  f = (2/3) dE |mu|^2;
* velocity form: the columns are the momentum <p>, and f = 2 |p|^2 / (3 dE),

with dE in Hartree (eV / 27.211386). The plugin recomputes f this way for every
state with f > 1e-3 and prints a warning on stderr if it differs from the printed
value by more than 5 percent. The printed value is always what is stored; the
check never alters or rejects anything. If the formulation line is missing or
not recognised, the check is skipped. (On the example file the velocity formula
reproduces the printed values to five significant figures.)

Restricted-spin triplet blocks (`R-TDDFPT states of multiplicity 3`) are parsed
but not put in the spectrum, since those transitions are spin-forbidden and
print f = 0.

## GEO_OPT, CELL_OPT and MD

CP2K does not reprint the geometry at every step of these runs; the main output
only holds the initial coordinates, and the trajectory goes to
`<project name>-pos-1.xyz` next to it. When the run type is `GEO_OPT`,
`CELL_OPT` or `MD` and Avogadro hands the plugin a real file path, the plugin
looks for that file:

* found, same atom count and same element sequence: its last frame is the
  geometry, and every frame becomes a coordinate set (`atoms.coords.3dSets`, the
  last one active). If every frame comment has an `E = ...` value (Hartree), these
  become `properties.energies`, which drives Avogadro's conformer energy plot.
  Coordinates are assumed to be in Angstrom, CP2K's default for this file.
* not found, unreadable or mismatched: the coordinates from the main output
  (the INITIAL structure) are used, and a warning on stderr says that the
  energy, charges and spectrum are from the FINAL step while the geometry is not.

If the last frame's energy differs from the output's final energy by more than
1e-5 Hartree, a warning notes that the trajectory may belong to another run.

## Known limitations

* **Charges, spectrum and energy always come from the last step** in the main
  output; earlier steps are not exposed.
* **CELL_OPT:** the unit cell is the last `CELL|` block in the main output,
  which is normally the initial cell. The relaxed cell is in `<project>-1.cell`
  and is not read, so the atoms and the cell can belong to different steps. A
  warning is printed.
* **Cell precision:** CP2K prints the cell vectors with three decimals but the
  lengths and angles with six. When the vectors are the standard orientation
  (a along x, b in the xy plane) of those lengths and angles, the plugin rebuilds
  them from the six-decimal values; otherwise it keeps the printed vectors.
* **Partial periodicity** (`XY`, `Z`, ...) is still given a full 3D cell.
* **QM/MM and multiple force evaluations:** the last coordinate table found is
  used and it may describe only the QM region.
* **Spin-polarised output formats** (UKS Mulliken/Hirshfeld, U-TDDFPT) are
  handled according to the layout described in CP2K's source and in the
  upstream parser, but have not yet been tested against a real file.
  Atoms whose coordinate table is cut off, or whose count disagrees with the
  `- Atoms:` total, make the read fail rather than produce a partial molecule.
* **Very long MD trajectories** are loaded in full; there is no striding.
* Writing CP2K files is not supported.

## Vendored code

`src/avogadro_cp2k/vendored/cp2k_output_tools/` contains parts of
[cp2k-output-tools](https://github.com/cp2k/cp2k-output-tools) (MIT licence,
Copyright (c) 2022 Tiziano Mueller), commit
`403cc70966d7284d9f71f56b8cd2ae08f969fa62`. The licence text and a list of
every local change are in that directory (`LICENSE`, `VENDORED.md`). The other
matchers in `src/avogadro_cp2k/blocks/` follow the same style but are original to
this plugin.

## File routing

`.out` and `.log` are written by many programs, so Avogadro cannot choose a
reader from the extension alone. The `patterns` key of this plugin's
`[[tool.avogadro.file-formats]]` entry in `pyproject.toml` (`"CP2K| version
string"`, a banner every CP2K run prints) lets Avogadro's generic-output
reader scan the file and hand CP2K output to this plugin. Older Avogadro
versions ignore the key, so there a CP2K `.out` file has to be renamed to
`.cp2k`, which is matched by extension. The plugin deliberately does not
register `out` or `log`: that would send every `.out`/`.log` file, whatever
program wrote it, to this plugin.

## Installing for testing

Use Extensions => Manage Plugins to install from a local directory. Avogadro
creates the pixi environment on first use, so the first read pauses while
`numpy` and `regex` are fetched.

## Development

```bash
pixi run -e test test        # pytest
```

The tests use excerpts of a real CP2K 2026.2 output (`tests/data/snippets`, which
are checked to be verbatim lines of the full file) and the full file itself,
`tests/data/PBE-SVP.out.gz`. Tests that need output nobody has supplied yet
(unrestricted runs, triplets, a real GEO_OPT with its trajectory) are marked as
skipped with a reason.

To exercise the plugin the way Avogadro does:

```bash
echo '{"operation":"read","filename":"tests/data/PBE-SVP.out.gz"}' \
  | avogadro-cp2k cp2k --read --debug
```
