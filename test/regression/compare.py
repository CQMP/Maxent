#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Compare regression results against the references.

    python3 test/regression/compare.py REFERENCE_DIR RESULT_DIR [--sets SETS] [--cases REGEX] [--report]

For every reference file <case>.h5 the result file of the same name must exist.
Every dataset under /files, /log/scalars and /cli is compared:

  numeric   same shape, same NaN pattern, and
            max|a - b| <= atol + rtol * max|reference|
            with (atol, rtol) from tolerances.json (first matching pattern wins)
  text      exactly equal

Datasets present on only one side fail.  /log/stdout and /log/stderr are
informational and not compared.  The run status must equal the expected status.

--report prints the difference of every numeric dataset (for measuring
tolerances) and never fails.
"""

import argparse
import json
import re
import sys
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
COMPARED = re.compile(r"^(files|log/scalars|cli)/")
BOOTSTRAP_DATASET = "files/case.out.booterr.dat"


def load_tolerances(path):
    rules = json.loads(path.read_text())["rules"]
    return [(re.compile(r["pattern"]), float(r["atol"]), float(r["rtol"])) for r in rules]


def tolerance(rules, case, dataset):
    key = f"{case}:{dataset}"
    for rx, atol, rtol in rules:
        if rx.search(key):
            return atol, rtol
    raise KeyError(f"no tolerance rule matches {key}")


def datasets(h5):
    out = {}

    def visit(name, obj):
        if isinstance(obj, h5py.Dataset) and COMPARED.match(name):
            out[name] = obj[()]

    h5.visititems(visit)
    return out


def as_text(x):
    return x.decode() if isinstance(x, bytes) else x


def validate_bootstrap(reference, result):
    """Validate stochastic bootstrap output without fixing a library's PRNG transform."""
    a, b = np.asarray(reference, dtype=float), np.asarray(result, dtype=float)
    problems = []
    if a.shape != b.shape or b.ndim != 2 or b.shape[1] != 4:
        return [f"booterr shape {b.shape}, expected {a.shape} with four columns"]
    if not np.all(np.isfinite(b)):
        problems.append("booterr contains non-finite values")
        return problems
    if not np.array_equal(a[:, :2], b[:, :2]):
        problems.append("booterr frequency or spectrum columns differ")
    if np.any(b[:, 2:] < 0):
        problems.append("booterr mean or error estimate is negative")

    spectrum_scale = float(np.max(np.abs(b[:, 1]), initial=0.0))
    if spectrum_scale == 0:
        problems.append("booterr spectrum is identically zero")
        return problems
    mean_deviation = float(np.max(np.abs(b[:, 2] - b[:, 1]), initial=0.0))
    if mean_deviation > 0.25 * spectrum_scale:
        problems.append("booterr bootstrap mean is inconsistent with the spectrum")
    max_error = float(np.max(b[:, 3], initial=0.0))
    if max_error == 0 or max_error > 5.0 * spectrum_scale:
        problems.append("booterr uncertainty estimate is zero or implausibly large")
    return problems


def compare_case(ref_path, res_path, rules, report):
    case = ref_path.stem
    problems, rows = [], []
    if not res_path.exists():
        return [f"{case}: result file missing"], rows
    with h5py.File(ref_path, "r") as fr, h5py.File(res_path, "r") as fc:
        status, expect = as_text(fc.attrs["status"]), as_text(fr.attrs.get("expect", "ok"))
        if status != expect:
            problems.append(f"{case}: status {status!r}, expected {expect!r}")
        ref, res = datasets(fr), datasets(fc)
    for name in sorted(set(ref) | set(res)):
        if name not in res:
            problems.append(f"{case}: {name} missing in result")
            continue
        if name not in ref:
            problems.append(f"{case}: {name} not in reference")
            continue
        a, b = ref[name], res[name]
        if isinstance(a, (bytes, str)) or getattr(a, "dtype", None) is not None and a.dtype.kind in "SOU":
            if as_text(a) != as_text(b):
                problems.append(f"{case}: {name} text differs")
            continue
        a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
        if name == BOOTSTRAP_DATASET:
            problems.extend(f"{case}: {p}" for p in validate_bootstrap(a, b))
            continue
        if a.shape != b.shape:
            problems.append(f"{case}: {name} shape {b.shape}, reference {a.shape}")
            continue
        if not np.array_equal(np.isnan(a), np.isnan(b)):
            problems.append(f"{case}: {name} NaN pattern differs")
            continue
        finite = ~np.isnan(a)
        diff = float(np.max(np.abs(a - b)[finite], initial=0.0))
        scale = float(np.max(np.abs(a)[finite], initial=0.0))
        rel = diff / scale if scale > 0 else (0.0 if diff == 0 else np.inf)
        rows.append((case, name, diff, scale, rel))
        if not report:
            atol, rtol = tolerance(rules, case, name)
            if diff > atol + rtol * scale:
                problems.append(f"{case}: {name} max|diff| {diff:.3e} (relative {rel:.3e}) "
                                f"exceeds atol {atol:g} + rtol {rtol:g} * {scale:.3e}")
    return problems, rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("reference", type=Path)
    ap.add_argument("result", type=Path)
    ap.add_argument("--cases", default=".*", help="regular expression on case names")
    ap.add_argument("--sets", help="comma-separated sets to compare "
                                   "(full,fast,targeted,cli,kk,legendre_convert,components); "
                                   "default: all")
    ap.add_argument("--tolerances", type=Path, default=HERE / "tolerances.json")
    ap.add_argument("--report", action="store_true", help="print all differences, never fail")
    args = ap.parse_args()

    rules = [] if args.report else load_tolerances(args.tolerances)
    refs = sorted(p for p in args.reference.glob("*.h5") if re.search(args.cases, p.stem))
    if args.sets:
        wanted = set(args.sets.split(","))

        def ref_set(path):
            with h5py.File(path, "r") as h5:
                return as_text(h5.attrs["set"])

        refs = [p for p in refs if ref_set(p) in wanted]
    if not refs:
        print("no reference files selected")
        return 1
    failed = 0
    for ref in refs:
        problems, rows = compare_case(ref, args.result / ref.name, rules, args.report)
        if args.report:
            for case, name, diff, scale, rel in rows:
                print(f"{case}\t{name}\t{diff:.3e}\t{scale:.3e}\t{rel:.3e}")
            for p in problems:
                print(f"# {p}")
            continue
        print(f"{'FAIL' if problems else 'ok  '} {ref.stem}")
        for p in problems:
            print(f"     {p}")
        failed += bool(problems)
    bootstrap_cases = {"t_generate_err", "t_generate_err_seed"}
    selected_cases = {ref.stem for ref in refs}
    bootstrap_results = [args.result / f"{case}.h5" for case in sorted(bootstrap_cases)]
    if (not args.report and bootstrap_cases <= selected_cases
            and all(path.exists() for path in bootstrap_results)):
        with (h5py.File(bootstrap_results[0], "r") as first,
              h5py.File(bootstrap_results[1], "r") as second):
            if BOOTSTRAP_DATASET in first and BOOTSTRAP_DATASET in second:
                default = first[BOOTSTRAP_DATASET][()][:, 2:]
                explicit = second[BOOTSTRAP_DATASET][()][:, 2:]
                if np.array_equal(default, explicit):
                    print("FAIL bootstrap seeds: default and explicit seeds produced identical estimates")
                    failed += 1
    if not args.report:
        print(f"{len(refs)} cases, {failed} failed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
