# Reference changes

Every regeneration of files in `reference/` is logged here: date, commit,
reason, affected cases, and the largest difference to the previous
references (from `compare.py --report`). Newest entry first.

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
