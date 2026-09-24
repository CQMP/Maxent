# Maxent → ALPS 3.0 migration plan

Status: **planning** · Written 2026-09-23 · Baseline commit `a399253` (master)

This document holds the working notes and the plan for moving Maxent from
ALPSCore 2.x to the modernized ALPS package (3.0, <https://github.com/ALPSim/ALPS>).
It has five steps:

| Step | Goal | Section |
|---|---|---|
| 1 | Analyze the code, find the tests, record the current state | [§1](#1-current-state-step-1--done) |
| 2 | Clean up: modern compilers, no errors or warnings, modern CMake | [§2](#2-step-2--cleanup-and-modernization) |
| 3 | Find what is in ALPSCore but not in ALPS 3.0, and where the incompatibilities are | [§3](#3-step-3--alpscore-vs-alps-30-incompatibility-analysis) |
| 4 | Port Maxent so it builds and runs against ALPS 3.0 | [§4](#4-step-4--port-maxent-to-alps-30-standalone) |
| 5 | Integrate Maxent into the ALPS source tree | [§5](#5-step-5--integrate-into-the-alps-repository) |

Decisions that need a human are collected in [§6](#6-open-decisions).

Reference trees used for this analysis:

* ALPS 3.0.0: `~/Projects/ALPS` (commit `f28d428`, `ALPS_VERSION.txt` = 3.0.0)
* ALPSCore 2.2.0: `~/Projects/ALPSCore` (fork `egull/ALPSCore`, commit `35e03ae`)

---

## 1. Current state (Step 1, done)

*This section is a snapshot of the baseline (`a399253`, 2026-09-23). Later
changes are recorded in the work log (Appendix B); removed components are
marked as such.*

### 1.1 Repository inventory

| Component | Path | LOC | Built by default | ALPSCore use | Other deps |
|---|---|---|---|---|---|
| `libmaxent` core library | `src/*.cpp` (7 files) | ~1,900 | yes | params, hdf5, utilities | Eigen3, Boost (math, random, algorithm, lexical_cast, shared_ptr, exception), GSL (one call), optional LAPACK |
| `maxent` executable | `src/maxent.cpp` | 170 | yes | params (`help_requested`, `get_origin_name`, …), `fs::remove_extensions` | Boost exception |
| Unit tests | `test/*Test.cpp` (6 files) | ~1,640 | yes (`Testing=ON`) | params, hdf5, `temporary_filename` | bundled Google Test (`gtest-all.cc`, 21k-line `gtest.h`) |
| `kk` (Kramers–Kronig) | `kk/kk.cpp` | 209 | yes | **none** | GSL (spline), Boost.ProgramOptions, OpenMP |
| `legendre_convert` | `legendre_convert/` | 579 | yes | **none** (links `ALPSCore_LIBRARIES` for no reason) | Boost (program_options, random, math) |
| ~~`pade`~~ | ~~`pade/pade_arbitrary_degree/`~~ | — | — | — | **Removed** (D7, 2026-09-23). At baseline: ~1,300 LOC, off by default, subclassed `alps::params`, needed GMP. |

Pade (removed, D7) had five unbuilt alternative sources at baseline
(`main.cpp`, `main_gmp.cpp`, `main_bary.cpp`, `lu.cpp`, `pade_interpolator_old.cpp`)
and a CMake file that hard-coded `link_directories("/opt/local/lib")`.

Files that exist only in the local checkout (ignored, but they clutter the work
tree): `build/` (CMake 3.14 cache), `submission/`, `submission.zip` (3.5 MB),
`theory/`, `doc/`, `examples/SpM/`, Eclipse `.project`/`.cproject`/`.settings`,
and `.DS_Store` files.

CI is `.travis.yml` (Ubuntu trusty, Boost 1.58, Travis). **It is dead.** No CI
currently runs.

### 1.2 Architecture (core library)

```
ContiParameters            reads data (text / HDF5 / param file), builds grid
 └─ MaxEntParameters       default model, kernel, covariance rotation, SVD of kernel
     └─ MaxEntHelper       transforms (u <-> A), chi2, entropy, Q, log-prob, backcontinuation, bootstrap errors
         └─ MaxEntSimulation   alpha loop + Levenberg–Marquardt, evaluate() writes text + HDF5 output
kernel        (maxent_kernel.*)   13 kernel types: time / frequency / Legendre x fermionic / bosonic / anomalous / T=0
grid          (maxent_grid.*)     real-frequency grids: lorentzian, half lorentzian, quadratic, log, linear
DefaultModel  (default_model.*)   flat, (shifted/double/two/general) gaussian, lorentzians, exp-decay, tabulated
Backcont      (maxent_backcont.*) A(omega) -> G(X)
eigen_hdf5.hpp                    ALPSCore hdf5 traits for Eigen::VectorXd
eigen_lapack.hpp                  optional LAPACK dgesvd path (USE_LAPACK)
```

The classes pass `alps::params&` through almost every constructor. **The
parameter object is the main coupling to ALPSCore.** HDF5 is used in exactly
two places: output in `MaxEntSimulation::evaluate()`, and input in
`ContiParameters::read_data_from_hdf5_file()`.

### 1.3 ALPSCore API surface actually used

| API | Where | Count |
|---|---|---|
| `alps::params` (ctor from `argc/argv`, from file, default; `operator[]`, `.as<T>()`, implicit conversions, `==`/`!=` against string/bool) | everywhere | 86 references |
| `params::define<T>(name[, default], descr)` | `maxent_simulation.cpp`, `maxent.cpp`, `maxent_params.cpp`, `maxent_kernel.cpp` | ~50 |
| `params::exists`, `defined`, `defaulted`, `help_requested(std::ostream&)`, `description`, `get_origin_name` (deprecated) | `maxent.cpp`, `maxent_params.cpp`, `default_model.hpp`, `maxent_kernel.cpp` | ~15 |
| **Runtime `define` after parsing** (`X_i`, `SIGMA_i`, `TAU_i`, `RUN_i` are defined on the fly) | `maxent_params.cpp:186`, `maxent_kernel.cpp` (TAU), `maxent.cpp` (RUN_) | these rely on ALPSCore-specific semantics |
| `alps::hdf5::archive` (`"r"`/`"w"`), `alps::make_pvp`, `ar["/path"] << x`, `alps/hdf5/vector.hpp` | `maxent_simulation.cpp`, `maxent_params.cpp`, tests | 14 `make_pvp` |
| `alps::hdf5::is_continuous`, `detail::get_extent`, `archive::write`, `archive_error`, `ALPS_STACKTRACE` | `eigen_hdf5.hpp` | 1 file |
| params ↔ HDF5 (`oar["/parameters"] << p`, `iar["/parameters"] >> p`) | `paramsTest.cpp` | 2 tests |
| `alps::temporary_filename` | tests | 12 |
| `alps::fs::remove_extensions` | `maxent.cpp:88` | 1 |
| `alps::cast` | `eigen_hdf5.hpp` (`using` only, never called) | 0 real uses |
| `<alps/config.hpp>` ("needed to set up correct bindings") | 5 files | obsolete ublas-bindings leftover |
| `alps::numeric::matrix` | Pade, **unbuilt files only** (Pade removed, D7) | — |

Boost usage outside ALPSCore: `shared_ptr` (31), `lexical_cast` (15),
`throw_exception` (11), `to_lower` (5), `math::isnan` (4), `random`
(`mt19937`, `normal_distribution`, `variate_generator`),
`math::legendre_p`/`factorial`/`sph_bessel`, `diagnostic_information`, and
`program_options` in the utilities. Some `using namespace boost::numeric;`
lines are leftovers from the old ublas code.

### 1.4 Build results on a modern toolchain (measured)

Toolchain: macOS 26.6 arm64, AppleClang 21.0.0, CMake 4.2.1, MacPorts Boost
1.81/1.88, Eigen 3.4.1, HDF5 2.1.1, OpenBLAS. GSL is **not** installed.

The unmodified tree does **not** configure. It stops at each of these in turn:

1. `cmake_minimum_required(VERSION 3.1)`: CMake 4 has removed compatibility with CMake < 3.5. The same problem is in `test/` (2.8.12), `kk/`, `legendre_convert/` and `pade/` (2.8).
2. The installed ALPSCore rejects the build on a compiler patch-level mismatch (`21000101` vs `21000334`). The workaround is `-DALPS_FORCE_NO_COMPILER_CHECK=ON`.
3. `find_package(GSL REQUIRED)` fails because GSL is missing.
4. The ALPSCore imported targets name `Boost::headers` in their link interface, but Maxent never calls `find_package(Boost)`, so generation fails.
5. At runtime, the local ALPSCore install is stale: it links `libhdf5.310.dylib` (HDF5 1.14), which MacPorts replaced with 2.1.1. Every test aborts in `dyld`.

Compiling each translation unit directly with `-std=c++17 -Wall -Wextra`
against ALPSCore gives:

| File | Errors | Maxent-owned warnings |
|---|---|---|
| `maxent_params.cpp` | 2: `using namespace boost::numeric` (namespace no longer declared) | sign-compare ×4, unused `threshold`, misleading indentation (l. 394) |
| `maxent_helper.cpp` | 1: same | unused `beta`, sign-compare ×2 |
| `maxent_simulation.cpp` | 1: same | ~30 sign-compare, `-Wreorder-ctor` (`qvec`/`nfreq`) |
| `maxent_kernel.cpp` | 1: `gsl/gsl_integration.h` not found | — |
| `maxent.cpp` | 0 | deprecated `get_origin_name`, sign-compare |
| `eigen_hdf5.hpp` | 0 | sign-compare |
| everything else | 0 | 0 |

Most of the ~300 warnings in the logs come from ALPSCore headers
(`dict_value_impl.hpp`, `hdf5/archive.hpp`, `-Wunused-parameter`).

**Baseline test run.** I made a scratch copy (the repo was not changed) and
applied four workarounds:

* removed the `using namespace boost::numeric;` lines
* replaced the single `gsl_integration_qag` call with `boost::math::quadrature::gauss_kronrod<double,61>`
* dropped `kk` from the build and added `find_package(Boost CONFIG)`
* rebuilt ALPSCore 2.2 against Boost 1.88 and HDF5 2.1.1, without MPI

With those, everything builds and **all 6 test executables (35 gtest cases)
pass in about 9 s**. All shipped examples run to completion:

| Example | Wall time |
|---|---|
| U0 (freq), U0 (tau), U2 freq/time/self, Bosonic, Bosonic PH, Self Energy U1/U10 | 1–4 s each |
| Legendre | **111 s** (with Boost quadrature in place of GSL; still needs a timing check against GSL) |

### 1.5 Test inventory

| Executable | Cases | What it covers |
|---|---|---|
| `gridTest` | 10 | all 5 grids, odd and even NFREQ: endpoints and monotonicity |
| `default_modelTest` | 6 | TabFunction parsing/junk/out-of-range; TwoGaussians ≡ Gaussian / DoubleGaussian |
| `paramsTest` | 9 | data from param file / text file / HDF5, covariance (text, HDF5), T=0 kernel, high-frequency check |
| `paramFailureTest` | 2 | missing `X_i`/`SIGMA_i`; NDAT larger than file |
| `backcontTest` | 4 | backcontinuation for frequency PH / non-PH / bosonic / tau |
| `simulationTest` | 4 | full runs (frequency fermionic, frequency bosonic, tau, Legendre): only sanity bounds (endpoint ≈ 0, norm ≈ 1, output sizes) |

**Gaps.**

* No golden-value regression tests. The simulation tests would not catch a 10% change in a spectrum.
* Nothing tests the `maxent` executable: CLI parsing, `--help`, `MODEL_RUNS`, or the output files and HDF5 layout.
* No kernel unit tests. `kernelTest` is mentioned in `maxent_kernel.cpp:39` but does not exist.
* No tests for `kk` or `legendre_convert` (nor for `pade`, since removed, D7).
* The `examples/` directory is not exercised by CI.
* The LAPACK SVD path (`USE_LAPACK`) is untested.

### 1.6 Code issues found during the read-through

These are candidates for step 2. "Verify" means I suspect a bug but have not
confirmed it numerically.

| # | Location | Issue |
|---|---|---|
| B1 | `src/maxent.cpp:163-168` | `main` catches every exception, prints it and **returns 0**, so scripts cannot detect failure. |
| B2 | `CMakeLists.txt:41` | `-O2 -DNDEBUG` is forced into `CMAKE_CXX_FLAGS` for every build type. Debug builds are not debug builds, and `#ifndef NDEBUG` checks never run. |
| B3 | `CMakeLists.txt:45` | `-msse2` whenever the compiler is GCC. **This breaks GCC on aarch64** (Linux ARM, Apple Silicon with GCC). |
| B4 | `src/maxent_backcont.cpp:36` | The size check is under `#ifdef NDEBUG`, which is inverted. It runs only in release builds. |
| B5 | `src/maxent_kernel.cpp:202,207` | `if(dtype_==legendre_dataspace)` sits inside the `dtype_==time_dataspace` branch, so `time_*_legendre_kernel` can never be selected (dead enum values). |
| B6 | `src/maxent_kernel.cpp:100` | A second `else if(ktype_==time_fermionic_kernel)` branch (which sets K=-1) can never be reached. |
| B7 | `src/maxent_kernel.cpp:81` | Time-bosonic: `K_(i,0) = T_` is immediately overwritten by the `j=0` loop iteration, so the ω→0 limit is lost. The result is NaN if a grid point sits at ω=0. Verify. **Confirmed:** a fix exists on the unmerged branch `xi_bose_bugfix` (`7d1ab5d` starts the loop at `j=1`; `c270d7f` maps overflow NaN to 0). **Fixed 2026-09-23** by merging `xi_bose_bugfix` into `modernize/step2`. |
| B8 | `src/eigen_hdf5.hpp:44`, `src/eigen_lapack.hpp:18` | Non-inline, non-template functions are defined in headers. Including either header in two translation units is an ODR/link error. |
| B9 | `src/maxent_params.cpp:53,73` | `while (datstream)` reads past EOF, and `expectedDatIn -= 1` hacks around the extra iteration. This is fragile with trailing blank lines. Verify. |
| B10 | `src/maxent_helper.cpp:361` | The bootstrap RNG is seeded with `time(0)`, so `GENERATE_ERR` output cannot be reproduced. |
| B11 | `src/maxent_params.cpp:304-306` | The `threshold` variable is unused. `JacobiSVD` is O(n³) and slow for large kernels; `BDCSVD` is the modern choice. |
| B12 | `src/maxent_simulation.cpp:25` | Member initialization order (`qvec` vs `nfreq`) gives a warning. It is harmless now but fragile. |
| B13 | `src/maxent.cpp:88` | Uses the deprecated `params::get_origin_name()`. |
| B14 | throughout | `std::size_t` vs `Eigen::Index` sign-compare in about 40 loops. |
| B15 | `src/maxent_grid.cpp` (log grid) | Intermediate `float` casts reduce the precision of the log grid. |
| B16 | `src/maxent.cpp:67` vs `src/maxent_grid.cpp:25` | `--help.grids` advertises `half-lorentzian` with `CUT=0.1`, but the grid code only accepts `half lorentzian` (with a space) and the defined default is `CUT=0.01`. A user following the help text gets "No valid frequency grid specified". Fix: accept both spellings and print the real default. |
| B7b | `src/maxent_kernel.cpp` (time bosonic) | The B7 fix sets **column 0** to the ω→0 limit T, but column 0 is the lowest grid frequency, not ω=0. With `OMEGA_MIN=0` the error is small (0.2 vs 0.2065 at ω=0.125, β=5); with a symmetric grid column 0 is ω=−OMEGA_MAX and T is badly wrong. Confirmed with the component dump. Fix: evaluate the formula everywhere, using the limit only for \|ω\| below a small threshold (or a series expansion). |
| B17 | `src/maxent_kernel.cpp` (`setup_legendre_kernel`) | The Legendre **bosonic** kernel uses the fermionic integrand (1+e^{−βω} denominator); it is bit-identical to the fermionic kernel (confirmed with the component dump). |
| B18 | usability | The default Lorentzian grid is centered at (OMEGA_MIN+OMEGA_MAX)/2. With `OMEGA_MIN=0` (T=0, bosonic) it has almost no points near ω=0 (lowest points 0.82 and 2.06 for NFREQ=200, OMEGA_MAX=10), so spectra with weight at low frequency cannot be fitted. Consider a different default grid when `OMEGA_MIN=0`, or a warning. |
| B19 | `src/maxent_simulation.cpp` (`levenberg_marquardt`) | The minimizer can diverge. Reproducer: `test/regression/inputs/t_model_quadratic_rise_exp_decay` with `--LAMBDA=1`: from the second α on, Q ≈ 1e26 and norm ≈ 1e8, every α hits `MAX_IT` (258 s); with `--MAX_IT=100` it stops with `Q=NaN, something went wrong`. Also diverges with the quadratic grid and with `OMEGA_MIN=0.2`, so it is not caused by the default model vanishing at ω=0. `LAMBDA=2` converges. Needs step-size control / a trust region. |
| B20 | `pade/pade_arbitrary_degree/` | Pade does not compile: missing `#include <iostream>`, ALPSCore changed `params::help_requested()`, and `std::complex<mpf_class>` is not supported by libc++ (the standard only allows `std::complex` of floating-point types). Left off (D7). **Resolved 2026-09-23 by removing Pade (D7).** |
| B21 | `src/` (logging) | Progress and diagnostic messages go to `std::cerr` (25 places in `src/`) with no convention and no verbosity control (issue #46, open since 2018; the old PR #47 `cerr -> cout` was closed as outdated on 2026-09-23). Decide a convention (results/progress to `cout`, warnings/errors to `cerr`), add a verbosity setting, and update the CLI regression references deliberately (they record stdout/stderr). |
| B22 | `src/maxent_helper.cpp` (`log_prob`, `chi_scale_factor`) | Performance: for large NFREQ the run time is dominated by dense NFREQ x NFREQ work per alpha (`K^T K`, a Cholesky factorization of an NFREQ x NFREQ matrix in `log_prob`); the Legendre example (NFREQ=5000) takes ~57 s. These determinants/eigenvalues could be computed in the singular space (ns x ns). Not a correctness issue; candidate for a later optimization with the regression suite as the safety net. |

---

## 2. Step 2: cleanup and modernization

**Goal:** the code, still on ALPSCore, builds warning-free with
`-Wall -Wextra -Wpedantic` on current GCC and Clang, uses modern
target-based CMake, and has a regression safety net. **Numerical output stays
unchanged** unless we fix a bug on purpose, and each such fix is documented.

### 2.0 Safety net first (before touching any numerics)

Design agreed 2026-09-23. References come from the **original numerics**:
commit `a399253` plus a minimal build-only patch (`boost::numeric` lines
removed, Boost found before ALPSCore, `CMAKE_POLICY_VERSION_MINIMUM=3.5`),
built with **GSL 2.8** against ALPSCore `3606edfb`.

Measured facts that the design relies on:

* Runs are deterministic. Text outputs are byte-identical across runs, and HDF5 outputs differ only in file metadata (`h5diff`: 0 differences).
* HDF5 output is float64 at full precision. Text output has 10 significant digits.
* HDF5 output is about 25 KB per run; `spex.dat` is about 800 KB per run.

**What is captured**

| Set | Cases | Quantities |
|---|---|---|
| Full examples (CTest label `regression-full`, opt-in) | 11 runs: U0 `in`/`in_tau`; U2 `frequency`/`time`/`self`; Bosonic `bosonic`/`bosonic_ph`; Legendre `in`; Self Energy `U1/in`, `U10/in`, `U10/green` | HDF5: `/alpha/{values,probability}`, `/spectrum/{omega,average,maximum,chi,variance}`, `bosonic`/`anomalous` groups. Text-only: χ²(α), `*_back.dat`, `fits.dat`, `*_self.dat`, bosonic spectra. Log scalars: minimal χ², number of singular values, posterior probability of the default model, Ng, χ scale factor. Not captured: `spex.dat`. |
| Fast set (default `ctest`) | The same 11 inputs with command-line overrides `--NFREQ=200 --N_ALPHA=20` (Legendre `--NFREQ=500`); `NDAT` unchanged (truncating τ data changes the physics). About 1–2 s each. | same |
| Targeted cases | Covariance (text, HDF5); `DATA_IN_HDF5`; `X_i`/`SIGMA_i` input; explicit `TAU_i`; T=0, time-bosonic, anomalous (PH and non-PH) kernels; every default model incl. tabulated; every grid; `MODEL_RUNS` (varspec) | same, plus a component harness that dumps the kernel matrix, its singular values, the grid and the discretized default model for small sizes |
| CLI snapshots | `--help`, `--help.models`, `--help.grids`, missing-parameter error | exact text (exit codes once B1 is fixed) |

Not captured:

* `GENERATE_ERR`: seeded from the clock (B10). It gets a reference after the seed parameter is added.
* The LAPACK SVD path: compared against the Eigen results within tolerance, with no separate reference.

Cases that run through suspected bugs record current behavior and are flagged
in the manifest. The baseline is `modernize/step2`, so it already contains the
B7 fix (merged from `xi_bose_bugfix`). Grid cases use the spelling the code
accepts (`half lorentzian`, see B16).

**Builds**

* Reference build: the baseline in a `git worktree`, plus `baseline.patch`, with GSL 2.8 and ALPSCore `3606edfb`.
* Variant builds, used only to measure tolerances:
  * AppleClang Debug with `-O0`
  * `clang++-mp-22`
  * `g++-mp-15` (needs a second ALPSCore built with GCC, because GCC's libstdc++ can't link against the libc++ ALPSCore; `-msse2` removed, B3)
  * `USE_LAPACK=1`

**Layout and documentation** (`test/regression/`)

| File | Contents |
|---|---|
| `README.md` | Purpose, how to run, how to regenerate, tolerance policy |
| `cases/<name>/` | Parameter file and input data (or a pointer into `examples/`) |
| `reference/<name>.h5` | All captured quantities as float64, log scalars as attributes, provenance as root attributes |
| `MANIFEST.md` | Per case: source, coverage (kernel, data space, grid, model, features), runtime, compared quantities, tolerances, bug flags |
| `PROVENANCE.md` + `baseline.patch` | Commits, exact patch, compiler, flags, library versions (Boost 1.88, Eigen 3.4.1, GSL 2.8, HDF5 2.1.1), OS and architecture, date |
| `REFERENCE_CHANGES.md` | Log of every regeneration: reason, affected cases, largest difference |
| `generate.py`, `compare.py` | numpy/h5py tools. `compare.py` is registered with CTest and reports the worst element per quantity. |

**Tolerances**

The check is |a−b| ≤ atol + rtol·max|ref| per dataset. The values are set **by
measurement**, not guessed: generate each reference with several build
variants (AppleClang -O0/-O2, a second compiler, Eigen vs LAPACK SVD) and use
the observed spread times a safety margin. Record the chosen values in
`MANIFEST.md`.

### 2.1 CMake modernization

Decisions (2026-09-23): project version **2.0.0**; default build type
**Release**; old option names (`Testing`, `PADE`, `USE_LAPACK`) are **dropped**,
not aliased; GoogleTest via `FetchContent`, pinned to the **newest release
(1.18.0)**. Pade was kept off in 2.1 because it did not compile (B20); it has
since been **removed entirely** (D7, resolved).

**Done 2026-09-23.** Results are **bit-identical** to the references
(Release `-O3 -std=c++17`: 1157 fast/targeted/CLI/component datasets and 237
full-size datasets, 0 differences). CTest: 35 unit tests (listed individually)
plus the regression tests, all passing with AppleClang (Release, Debug,
ASan+UBSan with no sanitizer reports), clang 22 and GCC 15. Notes:

* ALPSCore's installed config never passes its Boost hint on: it sets
  `alps_BOOST_DIR_` but reads `alps_Boost_DIR_` (case mismatch). Maxent
  therefore requires `Boost_DIR` (or `CMAKE_PREFIX_PATH`) and checks that the
  version equals `ALPSCore_BOOST_VERSION`. **Fixed upstream** by ALPSCore#666
  (merged 2026-09-23, `8d2ed3a9`): against ALPSCore `master` at or after that
  commit, Maxent configures without `Boost_DIR`.
* The ALPSCore install ships gtest 1.16 headers in its include directory; they
  shadowed the fetched 1.18 headers (link errors). Tests now link GoogleTest
  first.
* GCC 15 cannot compile GoogleTest against the macOS 26 SDK (`mach/message.h`);
  it works with `-DCMAKE_OSX_SYSROOT=.../MacOSX27.sdk`. GCC also needs
  `-include algorithm` until ALPSCore#667 is merged, and
  `MAXENT_BUILD_UTILITIES=OFF` because MacPorts' Boost.ProgramOptions is built
  against libc++.
* Warnings (clean build, `-Wall -Wextra -Wpedantic`): 51 with Clang, all in
  Maxent sources (37 sign-compare, 10 unused-variable, 1 each unused-parameter,
  reorder, misleading-indentation, deprecated). GCC additionally reports 46
  `-Wmaybe-uninitialized` inside libstdc++ headers (known false positives at
  `-O3`). These are the targets for 2.3.

* One top-level `cmake_minimum_required(VERSION 3.22)`, matching ALPS 3.0. Subdirectories drop their own `cmake_minimum_required`/`project()` or become proper subprojects.
* `set(CMAKE_CXX_STANDARD 17)`, `CMAKE_CXX_STANDARD_REQUIRED ON`, `CMAKE_CXX_EXTENSIONS OFF`.
* **No global flag injection.** Remove the forced `-O2 -DNDEBUG` (B2) and `-msse2` (B3). Warnings go on through a `maxent_warnings` INTERFACE target, with an optional `MAXENT_WERROR`.
* Targets and imported dependencies:
  * `find_package(Boost 1.76 CONFIG REQUIRED COMPONENTS program_options)`, called **before** ALPSCore (fixes build blocker 4).
  * `find_package(Eigen3 3.3 CONFIG REQUIRED)` → `Eigen3::Eigen`. Delete `cmake/FindEigen3.cmake` and the `ALPSCore_HAS_EIGEN_VERSION` branch.
  * `find_package(LAPACK)` → `LAPACK::LAPACK`, behind `option(MAXENT_USE_LAPACK OFF)`.
  * Library target `maxent::core` (rename from `libmaxent` and drop the `PREFIX ""` hack), with `target_include_directories(... PUBLIC $<BUILD_INTERFACE:...>)` and `target_compile_definitions` in place of the generated config header, or keep `configure_file` but attach it to the target.
* Options are namespaced: `MAXENT_BUILD_TESTS`, `MAXENT_BUILD_UTILITIES` (instead of `Testing`, `USE_LAPACK`; `MAXENT_BUILD_PADE`/`PADE` went away with Pade, D7).
* Google Test: delete the bundled 2013-era `gtest-all.cc`/`gtest.h`. Use `find_package(GTest)` and fall back to `FetchContent`, then `gtest_discover_tests()`. Replace `cmake/EnableGtests.cmake`.
* Utilities: `kk`/`legendre_convert` link only what they use (`Boost::program_options`, GSL or its replacement, `OpenMP::OpenMP_CXX`).
* Add `install(EXPORT)` and a `maxentConfig.cmake` only if we want downstream consumers. This is probably unnecessary given step 5.
* Add a `CMakePresets.json` (dev, asan, release), mirroring ALPS 3.0, which ships one.

### 2.2 Dependency reduction

**Done 2026-09-24** (branch `modernize/step2.2`). Decisions: the Legendre kernel
uses its **closed form**; `kk` got regression cases first; the `e.what()` error
text, the `SEED` parameter and the old ublas code moved to 2.3.

* `kk` regression cases added (3 cases, references from the GSL build); `kk`
  writes 17 significant digits.
* `kk`: GSL spline replaced by a natural cubic spline with the same algorithm;
  references reproduced **bit for bit** (serial and OpenMP), same speed.
* Legendre kernel in closed form, `K = -sqrt(2l+1) beta (-1)^l i_l(a) / (2 cosh a)`
  with `a = beta omega / 2`, from scaled modified spherical Bessel functions. It
  reproduces the GSL kernels to **3.3e-16**; no reference changed. GSL removed
  from the build.
* Boost utilities replaced by the standard library (`shared_ptr`,
  `lexical_cast`, `throw_exception`, `to_lower`, `isnan`, Boost.Random);
  bit-identical results.
* Remaining Boost: in the core library only `boost::diagnostic_information`
  (goes in 2.3), so the core will then use no Boost directly; Boost.Math is no
  longer used there. The utilities keep Boost: `kk` uses `program_options`, and
  `legendre_convert` uses `program_options`, Boost.Random and Boost.Math
  (`legendre_p`, `factorial`, `sph_bessel`). Decision (2026-09-24): keep
  Boost.Math; it needs C++14 (Boost >= 1.82) and Maxent uses C++17, and Boost is
  a dependency anyway through ALPSCore and ALPS.
* Warnings: 51 → 44. The Legendre example runtime is unchanged (~57-60 s):
  the kernel was never the bottleneck (B22).

| Item | Action | Notes |
|---|---|---|
| `cmake/FindGSL.cmake` | Delete; it uses the deprecated `EXEC_PROGRAM` (CMake dev warnings). If GSL stays anywhere, use CMake's built-in `FindGSL` → `GSL::gsl`. | |
| GSL (core, one `gsl_integration_qag` call) | Replace with `boost::math::quadrature::gauss_kronrod` (header-only; ALPS ships Boost anyway) | **Validated in 2.0:** `boost::math::quadrature::gauss_kronrod<double,61>` reproduces the GSL Legendre kernels to 2e-16 (machine precision). **Performance:** the Legendre example takes 55 s with GSL and 111 s with the naive Boost swap; the replacement must at least match GSL. Profile it, cache `legendre_p` via recurrence, use an adaptive depth, or parallelize over (l, j). |
| GSL (`kk`, cubic spline) | Replace with a small natural-cubic-spline implementation, or `boost::math::interpolators::cardinal_cubic_b_spline` (needs a uniform grid; check), or keep GSL optional for `kk` only | Decide per [§6](#6-open-decisions) D6. |
| `boost::shared_ptr` | `std::shared_ptr` / `std::unique_ptr` (the default model is owned uniquely) | mechanical |
| `boost::lexical_cast<std::string>(int)` | `std::to_string` | mechanical |
| `boost::math::isnan` | `std::isnan` | mechanical |
| `boost::throw_exception` | `throw` | mechanical |
| `boost::to_lower` | a local `to_lower` helper | trivial |
| `boost::random` | `<random>` (`std::mt19937`, `std::normal_distribution`), with a `SEED` parameter (fixes B10) | changes the bootstrap output stream, which is acceptable |
| `boost::diagnostic_information` | `e.what()` | |
| `boost::math::{legendre_p, factorial, sph_bessel}` | **keep** | header-only; C++17 `std::legendre`/`std::sph_bessel` are missing from libc++ |
| `boost::program_options` (utilities) | keep for now | ALPS 3.0 builds `program_options` from its bundled Boost |

After this step, the core library's only Boost dependency is header-only Boost.Math.

### 2.3 Code fixes

* Fix B1 through B20 (§1.6), including B7b. B21 (logging convention and verbosity) can go here or into a later step; it changes the CLI regression references on purpose.
* Moved here from 2.2 (they change the CLI references on purpose): replace
  `boost::diagnostic_information` by `e.what()` (with B1), add a `SEED`
  parameter for the bootstrap (B10), and delete the commented-out ublas and
  LAPACK-bindings code.
* `legendre_convert` (decided 2026-09-24), in this order:
  1. Add regression cases for `legendre_convert` first (it has no tests), with
     references from the current Boost build, like the `kk` cases in 2.2.
  2. Replace Boost.Random by `<random>`. The `mt19937` engine gives identical
     numbers; the normal variates differ (different algorithm), which is
     harmless because the error estimate is seeded from the clock.
  3. Replace `boost::math::factorial` by a product (only small arguments
     occur).
  4. Replace `boost::math::legendre_p` by the standard three-term recurrence.
     `std::legendre` is not an option: libc++ (AppleClang) does not implement
     the C++17 special math functions.
  Keep `boost::math::sph_bessel` (also missing in libc++; our own version
  would need careful checking at large l and argument) and `program_options`
  (no standard equivalent). B5, B6 and B7 change reachable behavior, so each gets its own commit with a before/after test.
* Remove all `using namespace boost::numeric;` lines, `#include <alps/config.hpp>`, dead commented-out ublas and lapack-bindings code, and the unused `alps::cast`.
* Make `eigen_hdf5.hpp`/`eigen_lapack.hpp` functions `inline`, or move them into `.cpp` files (B8).
* Use `Eigen::Index` for loop indices (B14). Consider `BDCSVD` in place of `JacobiSVD` (B11); that one is a numerics change and needs checking against the references.
* Keep the `MaxEntSimulation` public getters stable, because the tests use them.

### 2.4 Hygiene and CI

* Delete `.travis.yml` and add GitHub Actions with this matrix:
  * Ubuntu: GCC 11/13/15, Clang 15/19
  * macOS arm64: AppleClang
  * one ASan/UBSan job
  * one job with `USE_LAPACK=ON`

  Build ALPSCore from source in CI and cache it.
* **Replace the allowlist `.gitignore`.** Done 2026-09-23: replaced with a minimal ignore list (build dirs, OS/editor files, Python caches). Local-only material (D10) is hidden per clone via `.git/info/exclude`, not in the repository.
* ~~Remove the dead Pade sources~~ Done: Pade was removed entirely (D7).
* Update the README: build instructions, dependency list (no GSL), and remove the Travis badge.

### 2.5 Isolate ALPSCore behind two seams (bridge into step 4)

This is the last task of step 2, still on ALPSCore:

1. **Parameters seam:** introduce `maxent::params`, Maxent's own class (§4.1). Every `alps::params` use in `src/` goes through it. At the end of step 2 it can still wrap `alps::params`.
2. **HDF5 seam:** move all archive access into `maxent_io_hdf5.{hpp,cpp}`, with functions such as `write_vector(path, Eigen::VectorXd)`, `read_vector(path)` and `write_params(...)`. Nothing else includes `<alps/hdf5...>`.
3. Tests use `std::filesystem::temp_directory_path()` in place of `alps::temporary_filename`.

**Step 2 acceptance:**

* Configures with CMake 3.22 through 4.x.
* Zero warnings with `-Wall -Wextra -Wpedantic` in Maxent-owned code on the CI matrix.
* All unit, regression and CLI tests pass.
* No GSL in the core.
* `grep -r "alps::" src/` matches only the two seam files.

---

## 3. Step 3: ALPSCore vs ALPS 3.0 incompatibility analysis

### 3.1 What ALPS 3.0 is

ALPS 3.0 is the revived "legacy" ALPS: the ALPS 2 code base modernized to
C++17, Boost ≥ 1.76 bundled from source, and CMake ≥ 3.22, with a Python
wheel (`pyalps`). **It is not a superset of ALPSCore.** ALPSCore was a
separate rewrite of a few ALPS pieces (params, hdf5, accumulators, mc, gf).
ALPS 3.0 only has the ALPS-2-era ancestors of those pieces.

### 3.2 API mapping

| ALPSCore feature used by Maxent | ALPS 3.0 counterpart | Compatibility | Consequence |
|---|---|---|---|
| `#include <alps/params.hpp>`, `alps::params` | `<alps/ngs/params.hpp>`, `alps::params` (NGS, 2010–2012 design) | **Incompatible: same name, different class.** NGS `params` has no `define<T>()`, no defaults or descriptions, no `exists`/`defaulted`, no `help_requested`/`description`, no `argc/argv` constructor, no `get_origin_name`. Defaults come from `p["X"] \| default`. | **The largest item.** See §4.1. |
| ALPSCore INI parser: `KEY=value`, `#` comments, quoted strings, dotted keys (`help.models`) | NGS `params(boost::filesystem::path)` parses with legacy `alps::Parameters` (Spirit; comments are `//` and `/* */`) | **Incompatible file format.** Every shipped `examples/*.param` uses `#` comments (`BETA=8  #inverse temperature`). | Keep Maxent's own parser so user files stay valid. |
| `--KEY=value` command-line overrides | none; ALPS apps use `alps::mcoptions` (`input_file`, `output_file`, time limit) | missing | Maxent-owned CLI parsing. |
| Runtime `define` of `X_i`/`SIGMA_i`/`TAU_i`/`RUN_i` | NGS params hold everything as untyped `paramvalue` strings; nothing needs defining | semantics differ | This simplifies things once we own the params class. |
| params ↔ HDF5 (`ar["/parameters"] << p`) | NGS `params::save/load(hdf5::archive&)` exist, **different on-disk layout** | incompatible | Only tests rely on this. Maxent reads user data from `/Data`, `/Error`, `/Covariance` (plain vectors), which stays compatible. |
| `alps::hdf5::archive(fname, "r"/"w")` | `archive(boost::filesystem::path, std::string mode)`; also `(string, int)` | compatible (same lineage) | Implicit string→path conversion works. |
| `alps::make_pvp`, `ar["/p"] << x`, `<alps/hdf5/vector.hpp>` | present (`hdf5/archive.hpp:419-448`, `archive::operator[]`) | compatible | |
| `hdf5::is_continuous`, `detail::get_extent`, `archive::write(path, ptr, size, chunk, offset)`, `is_group`, `delete_group` | present in `hdf5/archive.hpp` | **probably compatible; must compile-test** | Better: drop the Eigen trait specialization and write `std::vector<double>` or raw `write(path, data, extent)` through the I/O seam. |
| `archive_error`, `ALPS_STACKTRACE` | `alps/hdf5/errors.hpp`, `alps/ngs/stacktrace.hpp` | present, different header | Drop them from Maxent. |
| `<alps/config.hpp>` | `<alps/config.h>` (+ `alps/ngs/config.hpp`) | different header name | Maxent doesn't need it. Delete in step 2. |
| `alps/utilities/temporary_filename.hpp` | `alps/utility/temporary_filename.hpp` (**utility**, singular), same signature | path differs | Replace with `std::filesystem` in step 2. |
| `alps::fs::remove_extensions` | missing | missing | `std::filesystem::path` logic. Check whether ALPSCore strips *all* extensions or the last one, and keep the `BASENAME` default identical. |
| `alps::numeric::matrix` (Pade dead files only) | `alps/numeric/matrix.hpp` exists | irrelevant | No longer used: Pade removed (D7). |
| Eigen (ALPSCore optionally bundles it, `ALPSCore_HAS_EIGEN_VERSION`) | **ALPS 3.0 does not use Eigen anywhere** (it uses ublas plus LAPACK bindings) | missing | Eigen becomes a *new* dependency for ALPS. See D3. |
| GSL | not used by ALPS | missing | Removed in step 2. |
| Imported CMake targets (`alps::alps-params`, …) via `ALPSCoreConfig.cmake` | `ALPSConfig.cmake` sets **variables** (`ALPS_INCLUDE_DIRS`, `ALPS_LIBRARY_DIRS`, …). `UseALPS.cmake` **force-sets the compiler and `CMAKE_CXX_FLAGS` in the cache** unless `PREVENT_ALPS_COMPILERS` is set. No imported targets. | different model; UseALPS is intrusive | A standalone build needs a small `FindALPS3`/wrapper that creates an INTERFACE target, sets `PREVENT_ALPS_COMPILERS`, and avoids `UseALPS.cmake`. |
| ALPSCore compiled with C++11, strict | ALPS 3.0 compiles with C++17 **plus `-fpermissive`** and global `-DBOOST_NO_AUTO_PTR= -DBOOST_ALLOW_DEPRECATED_HEADERS …` | flags differ | `-fpermissive` can hide real errors. Keep a standalone strict build in CI after integration. |
| Boost from the system (`Boost::headers`) | ALPS builds Boost **from a bundled source tree** (FetchContent `boost_src`) as a single `boost` library; `program_options` is included | different provisioning | Inside ALPS, Maxent must use ALPS's Boost; never mix versions. |
| gtest | ALPS: `add_alps_test()` (runs an executable and diffs against `<name>.output`), optional Boost.Test, `pytest` | different framework | See D5. |

### 3.3 Conflicts that are not about APIs

1. **Name and functionality collision.** ALPS 3.0 already ships a Maxent:
   * `tool/maxent*.cpp` (~1,550 LOC; Fuchs, Pruschke, Troyer 2010, Gull 2012). It is the **direct ancestor of this code** (same `ContiParameters`/`MaxEntHelper`/`MaxEntSimulation` class names), still on ublas, NGS `mcbase` and `mcoptions`, and it reads its parameters from **HDF5/XML input files**. It installs as `bin/maxent`.
   * The Python module `pyalps.maxent_c.AnalyticContinuation(dict)`, which `lib/pyalps/maxent.py` wraps.
   * The test `maxent_linear_grid_numeric` (`tool/CMakeLists.txt:91-94`).

   Our binary name, class names (global namespace) and features overlap. **Integration must replace or rename, not add.** Existing pyalps users call `AnalyticContinuation(dict)` with the legacy parameter names.
2. **Licensing (resolved 2026-09-23).** All authors agreed to relicense. Maxent is now MIT and uses the ALPS header convention (`ALPS Project:` + `SPDX-License-Identifier: MIT`), so it is license-compatible with ALPS 3.0. `scripts/check_license_headers.py` enforces the headers.
3. **Output compatibility.** The legacy ALPS maxent writes a different `.out.h5` layout than ours, and pyalps plotting helpers may assume the legacy layout. This needs checking in step 5.
4. **Global namespace.** Maxent classes (`grid`, `kernel`, `Backcont`, `DefaultModel`, `Model`, `Gaussian`, …) and typedefs (`matrix_type`, `vector_type`) live in the global namespace. The legacy ALPS maxent uses the same names. Before the two can share a build, everything goes into `namespace alps::maxent` (or `maxent`).

### 3.4 Risk ranking

| Risk | Severity | Likelihood | Mitigation |
|---|---|---|---|
| Parameter semantics and file format differ | high | certain | Maxent-owned params (§4.1) with a regression suite over all example files |
| Collision with legacy ALPS `maxent` / pyalps users | medium | certain | Coordinate with ALPS maintainers; offer a compatibility shim (§5.3) |
| Eigen as a new ALPS dependency | medium | certain | Optional component, found or fetched (D3) |
| HDF5 archive subtle differences (Eigen traits, chunking, HDF5 2.x) | low–medium | possible | HDF5 seam plus the layout test (§2.0/5) |
| `-fpermissive` hides errors in the ALPS build | low | likely | Strict standalone CI job |
| Numerics drift from the GSL→Boost, Jacobi→BDC SVD, or RNG changes | medium | likely | Golden references and documented tolerances |
| Legendre-kernel performance regression | medium | observed (111 s) | Profile and optimize in step 2 |

---

## 4. Step 4: port Maxent to ALPS 3.0 (standalone)

**Goal:** Maxent, still in its own repository, builds against an *installed*
ALPS 3.0 (and no ALPSCore), and every step-2 test passes with identical or
tolerance-equal results.

### 4.1 Own the parameter layer (the largest work item)

Recommendation: a small Maxent-owned `maxent::params`, about 400–600 LOC with tests. The other options are worse:

* (a) Rewriting on top of NGS `alps::params` would lose `define`/help/defaults and break the `#`-comment file format.
* (b) Vendoring ALPSCore's `params` module (~2.8k LOC plus iniparser, MIT) would duplicate a class named `alps::params` inside ALPS, which is guaranteed confusion.

Requirements for `maxent::params`:

* `define<T>(name, default, description)` and `define<T>(name, description)` (required), for `T ∈ {bool, int, double, std::string}`.
* Parse `KEY=value` files with `#` comments and quoted strings exactly as ALPSCore does today (golden-test against every `examples/**/*.param`), then apply `--KEY=value` CLI overrides.
* `operator[]` returns a proxy with `.as<T>()`, implicit conversion, and `==`/`!=`. Keep this, so the ~200 call sites stay unchanged.
* `exists`, `defined`, `defaulted`, lazy definition of indexed keys (`X_i`, `SIGMA_i`, `TAU_i`, `RUN_i`), `help_requested(os)` with the same `--help` text, and `origin_name()`.
* Constructors from `alps::params` (NGS) and from a Python dict, for the pyalps path in step 5.

Develop it during step 2.5, while ALPSCore is still available for side-by-side comparison. That removes all risk from the switch.

### 4.2 HDF5 through ALPS 3.0

* Re-implement `maxent_io_hdf5.cpp` against `<alps/hdf5.hpp>` from ALPS 3.0.
* Prefer writing `std::vector<double>` (or `archive::write` with an explicit extent) over specializing `is_continuous<Eigen::VectorXd>`.
* The layout test from §2.0 item 5 must pass unchanged.
* Rewrite `paramsTest` HDF5 cases to write the input HDF5 files directly (`/Data`, `/Error`, `/Covariance`) rather than through params serialization.

### 4.3 CMake against ALPS 3.0

* `cmake/FindALPS.cmake` (or a thin wrapper around `ALPSConfig.cmake`): read `ALPS_ROOT_DIR`, `ALPS_INCLUDE_DIRS` and `ALPS_LIBRARY_DIRS`, create `ALPS::alps` as an IMPORTED/INTERFACE target, and **do not** include `UseALPS.cmake` (or set `PREVENT_ALPS_COMPILERS`).
* Use Boost from `ALPS_Boost_INCLUDE_DIR` so both sides compile against the same Boost headers.
* Keep an `MAXENT_ALPS_BACKEND={ALPSCore,ALPS3}` switch for one transition release, only if that is cheap. Since the seams isolate ALPS to one I/O file, it is.
* CI: build ALPS 3.0 from source with `-DALPS_BUILD_APPLICATIONS=OFF -DALPS_BUILD_PYTHON=OFF -DALPS_BUILD_TESTS=OFF` (library only), cache the install, then build and test Maxent. Include one job with ALPS's own flags (`-fpermissive …`) to find problems early.

### 4.4 Namespacing

Move everything into `namespace alps::maxent` (or `maxent`). Rename the generic global names (`grid`, `kernel`, `Model`) as needed. Do this here, before step 5, so the move into ALPS does not collide with `tool/maxent*`.

**Step 4 acceptance:**

* Builds against a clean ALPS 3.0 install on the CI matrix, with no ALPSCore present.
* All unit, regression and CLI tests pass.
* Every `examples/**/*.param` parses identically, with the same effective parameter set as under ALPSCore.
* HDF5 output layout is unchanged.

---

## 5. Step 5: integrate into the ALPS repository

**Goal:** Maxent lives in `ALPSim/ALPS`, builds as part of the ALPS CMake
build and the Python wheel, and replaces the legacy `tool/maxent` or
coexists with it in a controlled way.

### 5.1 Preconditions

* Relicensing to MIT is done (D1, resolved), and all Maxent files carry the ALPS SPDX header. Keep `scripts/check_license_headers.py` passing.
* Agreement with the ALPS maintainers on placement, naming, Eigen, and what happens to the legacy tool (D2–D4). Follow `CONTRIBUTING.md` and the repo's `CLAUDE.md`/`AGENTS.md`.

### 5.2 Layout (proposal)

```
applications/maxent/            (or tool/maxent/ replacing the legacy files)
  CMakeLists.txt                option ALPS_BUILD_MAXENT (default ON if Eigen found)
  src/  include/alps/maxent/    library alps_maxent + executable
  utilities/{kk,legendre_convert}
  test/                         unit tests + add_alps_test regression cases with .output files
tutorials/maxent-*/             the examples/ inputs (+ PDFs moved to documentation)
lib/pyalps/maxent.py            updated wrapper (see 5.3)
```

### 5.3 Handling the legacy ALPS maxent

Recommended path:

1. The new code installs as `maxent`. The legacy binary is renamed `maxent_legacy` for one release, marked deprecated, and then removed.
2. `pyalps.maxent_c.AnalyticContinuation(dict)` is re-implemented on the new library through `maxent::params(dict)`. Legacy parameter names get a translation table where the meaning matches, and a clear error where it does not.
3. The legacy test `maxent_linear_grid_numeric` either gets ported (same physics check: finite spectrum, sum rule) or retired with the legacy tool.
4. Check the pyalps helpers (`load.py`, `plot*.py`) for dependencies on the legacy `.out.h5` layout.

### 5.4 Build integration details

* Eigen: `find_package(Eigen3 CONFIG)`, falling back to `FetchContent` (header-only, MPL2, compatible with MIT). The wheel build (`pyproject.toml`/scikit-build-core) needs the same.
* Compile under ALPS's global flags (C++17, `-fpermissive`, Boost defines), and also keep a strict build job.
* Tests: unit tests through gtest (fetched) or converted to plain executables registered with `add_alps_test`. Regression cases through `add_alps_test` with `.output` references (see D5).
* Utilities: `kk`/`legendre_convert` link ALPS's bundled Boost (`program_options` is built).
* Docs: README content goes into the ALPS docs, and the PDFs go to `tutorials/` or the documentation site. Update `ACKNOWLEDGE.TXT`/`CITATION.md` (Levy, LeBlanc, Gull, CPC 215 (2017)).
* CI: add Maxent tests to the ALPS GitHub Actions matrix (GCC 10–15, Clang 13–22, macOS).

### 5.5 Afterwards

* Archive `CQMP/Maxent`, or keep it as a read-only mirror, with a README pointing to ALPS.
* Zenodo/DOI: keep the Maxent citation (arXiv:1606.00368, CPC 2017) in ALPS's citation metadata.

**Step 5 acceptance:**

* `cmake --build` and `ctest` in ALPS build and pass the Maxent tests on the ALPS CI matrix.
* The `pyalps` wheel contains the new `maxent_c`.
* The legacy tool is deprecated or removed as agreed.

---

## 6. Open decisions

| ID | Decision | Needed by | Recommendation |
|---|---|---|---|
| D1 | Relicense Maxent GPL v2 → MIT | step 5 | **Resolved 2026-09-23:** all authors agreed; Maxent is MIT in the ALPS format (see §3.3 item 2). The GPL v2 license and the unmerged GPLv3 attempt are superseded. |
| D2 | Replace the legacy ALPS `tool/maxent`, or coexist | step 5 | Replace, with a one-release `maxent_legacy` deprecation window |
| D3 | Eigen as a new ALPS dependency, or port numerics to ublas/LAPACK | step 4/5 | Keep Eigen (header-only; find, else FetchContent). A port would be large and risky. |
| D4 | Location in ALPS: `applications/maxent` vs `tool/maxent` | step 5 | `applications/maxent` (it is a full application with utilities and tests) |
| D5 | Test framework inside ALPS: keep gtest vs `add_alps_test`/Boost.Test | step 2 (framework choice), step 5 | Keep gtest for unit tests and use `add_alps_test` for example regressions. Ask ALPS maintainers whether a fetched gtest is acceptable. |
| D6 | GSL in `kk`: replace or keep optional | step 2 | **Resolved 2026-09-24:** replaced by a natural cubic spline (same algorithm, bit-identical results). GSL is no longer used anywhere. |
| D7 | Pade: keep (GMP dependency), fix up, or drop | step 2 | **Resolved 2026-09-23: removed.** Pade (all 13 files, including the five unbuilt ones) is deleted, together with `MAXENT_BUILD_PADE` and its README section. |
| D8 | Install GSL locally once to produce the original Legendre reference outputs | step 2.0 | **Resolved 2026-09-23:** GSL 2.8 installed via MacPorts (`/opt/local`) and detected by `cmake/FindGSL.cmake`. |
| D9 | Rebuild the stale local ALPSCore install | step 2 | **Resolved 2026-09-23:** `~/Projects/ALPSCore` switched to `master` (`3606edfb`, = v2.3.3 + merge; CMake package still reports 2.2.0), built in `build-master/` (C++17, MPI on, Boost 1.88, HDF5 2.1.1, RelWithDebInfo), 135/135 ALPSCore tests pass, clean install to `~/Projects/ALPSCore/install`. The install now also ships gtest/gmock 1.16. |
| D10 | What to do with untracked local directories | step 2 | **Resolved 2026-09-23:** deleted `doc/` (generated Doxygen output), the paper-submission snapshot (`submission/`, `submission.zip`), `theory/`, `examples/SpM/`, the Eclipse files in `src/`, and `pade/pade_arbitrary_degree.zip` (byte-identical copy of the tracked Pade sources on `master`). No local-only material remains. |

---

## Appendix A: Reproducing the baseline build (scratch, repo unchanged)

```bash
# 1. ALPSCore 2.2 rebuilt against current Boost/HDF5 (no MPI)
cmake ~/Projects/ALPSCore -DCMAKE_INSTALL_PREFIX=$S/alpscore-install \
  -DBoost_DIR=/opt/local/libexec/boost/1.88/lib/cmake/Boost-1.88.0 \
  -DCMAKE_PREFIX_PATH=/opt/local -DTesting=OFF -DENABLE_MPI=OFF \
  -DEIGEN3_INCLUDE_DIR=/opt/local/include/eigen3 -DCMAKE_BUILD_TYPE=Release
make -j10 install

# 2. Maxent copy with workarounds:
#    - delete `using namespace boost::numeric;` lines
#    - swap gsl_integration_qag -> boost::math::quadrature::gauss_kronrod<double,61>
#    - remove find_package(GSL) and add_subdirectory(kk)
#    - add find_package(Boost CONFIG REQUIRED) before find_package(ALPSCore)
cmake <copy> -DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
  -DALPSCore_DIR=$S/alpscore-install/share/ALPSCore \
  -DBoost_DIR=/opt/local/libexec/boost/1.88/lib/cmake/Boost-1.88.0 \
  -DCMAKE_PREFIX_PATH=/opt/local -DEIGEN3_INCLUDE_DIR=/opt/local/include/eigen3
make -j8 && ctest     # 6/6 executables, 35 cases pass, ~9 s
```

## Appendix B: Work log

| Date | Step | Note |
|---|---|---|
| 2026-09-23 | 1, 3 | Initial analysis; baseline build/test on AppleClang 21 / CMake 4.2; this plan written. |
| 2026-09-23 | 2.0 | GSL 2.8 installed (MacPorts); `FindGSL.cmake` finds it. Configure now gets past GSL and stops at the known Boost-before-ALPSCore issue. |
| 2026-09-23 | 2.0 | ALPSCore reinstalled from master (see D9). Maxent (only the `boost::numeric` fix, Boost found first, policy flag, **original GSL code**) builds incl. `kk` and passes 6/6 test executables against it. |
| 2026-09-23 | 0 | Branch audit: `new_alps` and `external_to_ALPS_branch` are fully merged into master (no hidden migration work). `xi_bose_bugfix` holds the B7 fix; `GPLv3` holds an unmerged relicensing attempt. Work branch `modernize/step2` created from master; this plan committed. |
| 2026-09-23 | 2.3 | Merged `xi_bose_bugfix` into `modernize/step2` (B7 fixed); 6/6 test executables pass. The regression references will therefore include the B7 fix. Local cleanup: deleted local `xi_bose_bugfix`, `GPLv3`, `HiroshiMethod`, `NNLS` (all fully on origin) and `tmp` (its unique commit `7967ffa` kept as tag `archive/tmp`); pruned stale remote-tracking refs. Six old stashes kept (2017 SpM/ADMM work in `stash@{1}`, `stash@{3}`). |
| 2026-09-23 | lic | Relicensed to MIT in the ALPS format: ALPS "Applications" header with the existing copyright line kept verbatim plus `ALPS Project:` and `SPDX-License-Identifier: MIT` lines (39 C++ files; third-party gtest/FindEigen3 notices untouched); `LICENSE.TXT` (GPL v2) → `LICENSE.txt` (ALPS MIT text, `Copyright 1998-2026 ALPS Collaboration`); `ACKNOWLEDGE.TXT` → `CITATION.md` (ALPS format); README MIT badge and license section; `.zenodo.json` license `MIT`; old Python 2 header scripts and `HEADER.TXT` replaced by `scripts/check_license_headers.py`. Build and 6/6 tests pass. |
| 2026-09-23 | 2.4 | Allowlist `.gitignore` replaced by a minimal one. In this clone, `doc/`, `submission/`, `submission.zip`, `theory/`, `examples/SpM/`, `pade/pade_arbitrary_degree.zip` and the Eclipse files in `src/` are listed in `.git/info/exclude` (local only) until D10 is decided. |
| 2026-09-23 | 2.4 | Deleted local-only `doc/`, `submission/`, `submission.zip`, `theory/`, `examples/SpM/`, the Eclipse files in `src/` and `pade/pade_arbitrary_degree.zip` (none were tracked). D10 resolved; `.git/info/exclude` back to its default. |
| 2026-09-23 | 2.0 | Dry run of the reference generation. Confirmed that CLI overrides (`--NFREQ=…`) work with parameter files. Found 11 (not 10) example runs, a filename-case bug in `examples/Legendre/in.param` (fixed), and B16. Fast-set sizes agreed. |
| 2026-09-23 | 2.0 | Reference generation started: baseline worktree at `6ba7250` + `baseline.patch`, GSL build; 56 cases (11 full, 11 fast, 29 targeted, 4 CLI, components) deterministic (two runs: 1209 datasets, 0 differences). Found B7b, B17, B18. Confirmed B3 (GCC on arm64: `unrecognized command-line option '-msse2'`). ALPSCore does not build with GCC 15 (missing `<algorithm>` in `params_impl.hpp`); fix and upstream PR delegated. Variant builds (-O0, LAPACK, clang 22, GCC 15) for tolerance measurement. |
| 2026-09-23 | 2.0 | **Step 2.0 done.** `test/regression/`: 56 cases (11 full, 11 fast, 29 targeted, 4 CLI, components), 1165 datasets, 4.7 MB of references, deterministic. Tolerances measured with four variant builds (-O0, LAPACK, clang 22, GCC 15); all pass. Detection check: the pre-B7 kernel is caught. GSL→Boost Legendre quadrature validated to 2e-16. Found B19 (minimizer divergence). ALPSCore GCC fix submitted as ALPSCore/ALPSCore#667. Not yet in CTest (step 2.1). |
| 2026-09-23 | 2.1 | **Step 2.1 done.** Build-only source fix committed; CMake rewritten (3.22...4.2, project version 2.0.0, C++17, default Release, `maxent::core` target, `MAXENT_*` options, GoogleTest 1.18.0 via FetchContent, regression suite in CTest, presets, GNUInstallDirs). Bit-identical to the references. See §2.1 for notes. Merged as PR #52 (`ccc3be4`) after two Copilot reviews (4 findings: 3 fixed, 1 kept with explanation; second review: approval recommended). |
| 2026-09-23 | 2.2 prep | Ready for step 2.2: `modernize/step2` at `ccc3be4` builds against ALPSCore `master` `8d2ed3a9` (C++17, Boost 1.88, found without `Boost_DIR`); 39/39 tests, regression bit-identical (920 + 237 datasets). Note: `~/Projects/ALPSCore/install` currently holds a C++11/Boost 1.81 build that Maxent rejects; use an install of current `master`. |
| 2026-09-23 | repo | Closed two outdated PRs: #35 (SpM notes, superseded by #45) and #47 (`cerr -> cout`, 2018). Issue #46 (logging to `cerr`) stays open; recorded as B21. |
| 2026-09-23 | 2.x | Pade removed (D7), on branch `modernize/remove-pade`: `pade/` (13 files: 6 built, 5 unbuilt alternatives, header, CMake), `MAXENT_BUILD_PADE`, and the README section. It required GMP (`mpf_class`, 256-bit default precision) and did not compile (B20). |
| 2026-09-24 | 2.2 | Step 2.2 done on `modernize/step2.2`: kk regression cases, kk spline, closed-form Legendre kernel (GSL removed), Boost utilities → std. No reference changed; 41/41 tests. ALPSCore#668 merged (GoogleTest 1.18, not installed, CMake 3.16). Found B22 (performance). |
