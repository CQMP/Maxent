# Reference changes

Every regeneration of files in `reference/` is logged here: date, commit,
reason, affected cases, and the largest difference to the previous
references (from `compare.py --report`). Newest entry first.

## 2026-09-24: time-bosonic zero-frequency limit corrected (B7b)

* **Source:** `modernize/step2.3c` after commit `46b5610`.
* **Reason:** the time-bosonic kernel previously assigned the ω→0 limit to
  column 0 regardless of that column's frequency. It now uses the limit only
  when ω is exactly zero and evaluates the kernel formula everywhere else.
* **Cases:** `t_kernel_time_bosonic` and `components`; added a symmetric-grid
  time-bosonic component that contains both negative and exactly zero ω.
* **Largest differences:** 2.202e-3 relative in the kernel component; 5.373e-4
  relative in χ² and 1.190e-2 relative in the small spectrum variance.

## 2026-09-24: reproducible bootstrap errors (B10)

* **Source:** `modernize/step2.3c` after commit `379b2f0`.
* **Reason:** bootstrap error generation now uses the user-selectable `SEED`
  parameter, which defaults to 0 instead of the wall clock.
* **Cases:** added `t_generate_err`, including `booterr.dat`; updated
  `cli_help` for the new parameter. No existing numerical reference changed.
* **Verification:** two independent runs with the default seed were identical.

## 2026-09-24: grid spelling and help corrected (B16)

* **Source:** `modernize/step2.3c`, based on merge commit `454f7ee`.
* **Reason:** expected help-text correction: `--help.grids` now reports the
  actual `CUT=0.01` default. The documented `half-lorentzian` spelling is now
  accepted as an alias for `half lorentzian`.
* **Cases:** `cli_help_grids` text only. The targeted half-Lorentzian case now
  uses the hyphenated alias and remains within its existing numerical
  tolerance; no numerical reference was regenerated.

## 2026-09-24: command-line failures return nonzero (B1)

* **Source:** `modernize/step2.3c`, based on merge commit `454f7ee`.
* **Reason:** expected behavior change: Maxent now returns status 1 on errors
  and reports standard exceptions with `e.what()` instead of Boost diagnostic
  formatting.
* **Cases:** `cli_missing_beta` only. Exit status changed from 0 to 1; stderr
  changed from the multi-line Boost diagnostic to
  `Caught Exception: Critical parameters not defined`.

## 2026-09-24: legendre_convert reference added

* **Source:** `modernize/step2` merge commit `a57f150`; the utility is still
  the original Boost implementation.
* **Reason:** `legendre_convert` had no tests; a reference was added before
  replacing Boost.Random, `factorial`, and `legendre_p` (plan step 2.3B).
* **Cases:** 1 new (`legendre_convert_transform`). No existing reference changed.
* **Coverage:** transform, tail enforcement, back-continuation, and Matsubara
  convergence. Zero input errors make the clock-seeded bootstrap deterministic.

## 2026-09-24: kk references added

* **Commit:** `811b3c6` (see [`PROVENANCE.md`](PROVENANCE.md), kk section).
* **Reason:** `kk` had no tests; references added before its GSL spline is
  replaced (plan step 2.2).
* **Cases:** 3 new (`kk_imag_to_real_green`, `kk_imag_to_real_self`,
  `kk_real_to_imag`). No existing reference changed.

## 2026-09-23: initial references

* **Commit:** baseline `6ba7250` plus `baseline.patch` (see [`PROVENANCE.md`](PROVENANCE.md)).
* **Reason:** first generation, before any modernization of the code (plan step 2.0).
* **Cases:** all 56 (11 full, 11 fast, 29 targeted, 4 CLI, components).
* **Largest difference:** none (no previous references). Two independent
  generations were bit-identical.
* **Notes:**
  * The numerics include the B7 fix (`xi_bose_bugfix`), merged before
    generation. The time-bosonic kernel still has B7b (see the plan).
  * `t_model_quadratic_rise_exp_decay` uses `LAMBDA=2`; with `LAMBDA=1` the
    minimizer diverges (B19), which would make the reference meaningless.
  * `*.fits.dat` is not packed (27 of 30 MB; see README.md). The references
    are 4.7 MB.
