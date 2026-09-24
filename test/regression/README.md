# Maxent regression suite

This suite records what Maxent computes **today** and flags any change.
It is the safety net for the modernization and the move to ALPS 3.0
(see `ALPS3_MIGRATION_PLAN.md`, step 2.0). The references record behavior,
not correctness: a case that runs through a known bug records the buggy
result, and is flagged in [`MANIFEST.md`](MANIFEST.md).

## Layout

| Path | Contents |
|---|---|
| `cases.py` | Definition of every case (name, set, inputs, parameter file, arguments, flags) |
| `make_inputs.py` | Generates the synthetic inputs in `inputs/` from known spectra (run once; output committed) |
| `inputs/<case>/` | Inputs of the targeted and CLI cases |
| `generate.py` | Runs the cases with a `maxent` binary and packs the results, one HDF5 file per case |
| `compare.py` | Compares results against `reference/` using `tolerances.json` |
| `components/` | `dump_components` (a target of the main build): kernel matrices, grids, default models and singular values at small size |
| `CMakeLists.txt` | Registers the suite with CTest |
| `reference/<case>.h5` | The references |
| `tolerances.json` | Comparison tolerances per dataset pattern |
| [`MANIFEST.md`](MANIFEST.md) | Every case: what it covers, runtime, known-issue flags; tolerance summary |
| [`PROVENANCE.md`](PROVENANCE.md) | How the references were produced |
| `baseline.patch` | The build-only patch applied to the baseline commit |
| [`REFERENCE_CHANGES.md`](REFERENCE_CHANGES.md) | Log of every regeneration |

## Sets

| Set | Cases | Runtime | Purpose |
|---|---|---|---|
| `fast` | the 11 example runs with `--NFREQ=200 --N_ALPHA=20` (Legendre `--NFREQ=500`) | seconds | default |
| `targeted` | 29 small runs on synthetic inputs | seconds | inputs, kernels, default models, grids, `MODEL_RUNS` |
| `cli` | `--help`, `--help.models`, `--help.grids`, missing `BETA` | instant | command-line behavior |
| `full` | the 11 example runs as shipped | about 1 min (Legendre dominates) | opt-in |
| `kk` | the `kk` utility on 3 inputs (Im G and Im Σ from the references, an analytic Re G) | seconds | registered when `MAXENT_BUILD_UTILITIES=ON` |
| components | `components.h5` from `dump_components` | instant | building blocks |

## Running

From a build tree (`MAXENT_BUILD_TESTS=ON`, the default):

```bash
ctest --test-dir build -L regression-fast   # fast, targeted, cli, components (seconds)
ctest --test-dir build -L regression-full   # needs -DMAXENT_REGRESSION_FULL=ON (about 1 min)
ctest --test-dir build -L regression        # every registered regression test
```

Each set is two CTest tests: `regression_<set>_generate` runs the cases into
`build/test/regression/results-<set>/`, and `regression_<set>_compare` checks
them against `reference/`. By hand:

```bash
python3 test/regression/generate.py --maxent build/maxent \
    --components build/test/regression/dump_components --out /tmp/maxent-results \
    [--kk build/kk/kk] [--sets fast,targeted,cli,full,kk] [--cases REGEX]
python3 test/regression/compare.py test/regression/reference /tmp/maxent-results \
    [--sets fast,targeted,cli,components] [--cases REGEX]
```

`compare.py` exits with status 1 if any case fails. `--report` prints the
difference of every dataset instead of judging it.

Requirements: Python 3 with numpy and h5py. Without them, CMake prints a
warning and does not register the regression tests.

## What is compared

Each result file holds:

* `/files/<output>`: every output file of the run. HDF5 outputs (`*.out.h5`) are
  copied dataset by dataset at full precision. Text tables (`*.dat`) are loaded
  as float64, but only carry the 10 significant digits Maxent writes.
  `*.spex.dat` and `*.fits.dat` (spectrum and fit at every α; large, and
  covered by `chi2.dat`, the α probabilities and the `*_back.dat` files) and
  `*.booterr.dat` (random, B10) are not packed.
* `/log/scalars/`: values parsed from the logs: minimal χ², posterior
  probability of the default model, Ng, number of singular values kept, and
  the number of α values where the minimizer hit `MAX_IT`.
* `/cli/`: stdout, stderr and exit code of command-line cases.
* `/log/stdout`, `/log/stderr`: kept for debugging, not compared.

Numeric datasets pass if `max|a - b| <= atol + rtol * max|reference|`, text
must match exactly, and both sides must contain the same datasets.

## Regenerating references

Regenerate only on purpose (a deliberate behavior change, or a bug fix), never
to make a failing comparison pass without understanding it:

1. Run `generate.py` into a scratch directory and `compare.py --report` against
   the current references; understand every difference.
2. Replace the affected files in `reference/`.
3. Add an entry to [`REFERENCE_CHANGES.md`](REFERENCE_CHANGES.md): date, commit,
   reason, affected cases, largest difference.
4. Update [`PROVENANCE.md`](PROVENANCE.md) if the build environment changed.
