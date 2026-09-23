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
    ap.add_argument("--sets", help="comma-separated sets to compare (full,fast,targeted,cli,components); "
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
    if not args.report:
        print(f"{len(refs)} cases, {failed} failed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
