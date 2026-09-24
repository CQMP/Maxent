#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Run the regression cases with a given maxent binary and pack the results.

    python3 test/regression/generate.py --maxent build/maxent --out /tmp/results \
        [--sets fast,targeted,cli] [--cases REGEX] [--components build/dump_components] \
        [--kk build/kk/kk] [--legendre-convert build/legendre_convert/legendre_convert] \
        [--provenance prov.json] [--jobs 8]

For every case one HDF5 file <out>/<case>.h5 is written:
  /files/<output file>          every output of the run: datasets of *.out.h5 files are
                                copied as-is; text tables (*.dat) are loaded as float64
                                (they carry 10 significant digits)
  /log/{stdout,stderr}          run logs (informational; not compared)
  /log/scalars/<name>           values parsed from the logs (see SCALARS)
  /cli/{stdout,stderr,exit_code} for command-line cases
Root attributes hold the case metadata, the command, the status and the provenance.

Not packed: *.spex.dat and *.fits.dat (spectrum and fit for every alpha: they
grow like N_ALPHA x NFREQ or N_ALPHA x NDAT, and are covered by chi2.dat, the
alpha probabilities and the *_back.dat files).
"""

import argparse
import json
import re
import shutil
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
sys.path.insert(0, str(HERE))
from cases import all_cases  # noqa: E402

SKIP_OUTPUTS = re.compile(r"\.(spex|fits)\.dat$")
SCALARS = {
    "minimal_chi2": re.compile(r"^minimal chi2: (\S+)", re.M),
    "posterior_probability": re.compile(r"^posterior probability of the default model: (\S+)", re.M),
    "ng": re.compile(r"^Ng: (\S+)", re.M),
}
SINGULAR = re.compile(r"^# (\d+)\t(\S+)$", re.M)
TIMEOUT = 1800


def run(cmd, cwd, timeout=TIMEOUT):
    """Run cmd; a timeout is reported as returncode None instead of raising."""
    t0 = time.monotonic()
    try:
        p = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired as e:
        out = e.stdout.decode(errors="replace") if isinstance(e.stdout, bytes) else (e.stdout or "")
        err = e.stderr.decode(errors="replace") if isinstance(e.stderr, bytes) else (e.stderr or "")
        p = subprocess.CompletedProcess(cmd, None, out, err)
    return p, time.monotonic() - t0


def pack_outputs(h5, workdir, before):
    files = h5.create_group("files")
    for path in sorted(workdir.rglob("*")):
        rel = path.relative_to(workdir).as_posix()
        if not path.is_file() or rel in before or SKIP_OUTPUTS.search(rel):
            continue
        if path.suffix == ".h5":
            grp = files.create_group(rel)

            def copy(name, obj, grp=grp):
                # visititems stops at the first callback that returns a value, so return None
                if isinstance(obj, h5py.Dataset):
                    grp.create_dataset(name, data=obj[()])

            with h5py.File(path, "r") as src:
                src.visititems(copy)
        else:
            try:
                data = np.loadtxt(path, comments="#", ndmin=2)
            except ValueError:
                files.create_dataset(rel, data=path.read_text())
                continue
            files.create_dataset(rel, data=data)


def pack_scalars(h5, log):
    grp = h5.create_group("log/scalars")
    for key, rx in SCALARS.items():
        grp.create_dataset(key, data=np.array([float(v) for v in rx.findall(log)], dtype=float))
    # singular values are printed once per kernel set-up: '# <s>\t<value>', s restarting at 0
    counts = []
    for s, _ in SINGULAR.findall(log):
        if int(s) == 0:
            counts.append(0)
        if counts:
            counts[-1] += 1
    grp.create_dataset("n_singular", data=np.array(counts, dtype=np.int64))
    # alpha values for which the Levenberg-Marquardt loop did not converge within MAX_IT
    grp.create_dataset("max_it_warnings", data=np.int64(log.count("reached max_it")))


def run_case(case, programs, outdir, provenance, timeout=TIMEOUT):
    out = outdir / f"{case['name']}.h5"
    with tempfile.TemporaryDirectory(prefix="maxent-reg-") as tmp:
        workdir = Path(tmp) / "work"
        if case["inputs"]:
            shutil.copytree(ROOT / case["inputs"], workdir)
        else:
            workdir.mkdir()
        before = {p.relative_to(workdir).as_posix() for p in workdir.rglob("*")}
        cmd = [str(programs[case["program"]])] + ([case["param"]] if case["param"] else []) + case["args"]
        p, seconds = run(cmd, workdir, timeout)
        # Maxent reports exceptions with a non-zero exit code and a diagnostic;
        # the other programs report failure through their exit code alone.
        status = "timeout" if p.returncode is None else \
            "exception" if "Caught Exception" in p.stderr else \
            "failed" if case["program"] != "maxent" and p.returncode != 0 else "ok"
        with h5py.File(out, "w") as h5:
            for key in ("name", "set", "covers"):
                h5.attrs[key] = case[key]
            h5.attrs["inputs"] = case["inputs"] or ""
            h5.attrs["flags"] = json.dumps(case["flags"])
            h5.attrs["command"] = json.dumps([case["program"]] + cmd[1:])
            h5.attrs["runtime_seconds"] = seconds
            h5.attrs["status"] = status
            h5.attrs["expect"] = case["expect"]
            for key, value in provenance.items():
                h5.attrs["provenance." + key] = value
            h5.create_dataset("log/stdout", data=p.stdout)
            h5.create_dataset("log/stderr", data=p.stderr)
            if case["set"] == "cli":
                h5.create_dataset("cli/stdout", data=p.stdout)
                h5.create_dataset("cli/stderr", data=p.stderr)
                h5.create_dataset("cli/exit_code", data=-1 if p.returncode is None else p.returncode)
            pack_outputs(h5, workdir, before)
            pack_scalars(h5, p.stdout + p.stderr)
    return case["name"], seconds, status, case["expect"]


def run_components(binary, outdir, provenance):
    out = outdir / "components.h5"
    with tempfile.TemporaryDirectory(prefix="maxent-comp-") as tmp:
        p, seconds = run([str(binary), tmp], tmp)
        with h5py.File(out, "w") as h5:
            h5.attrs.update(name="components", set="components",
                            covers="kernel matrices, grids, default models, singular values",
                            inputs="", flags="[]", command=json.dumps(["dump_components"]),
                            runtime_seconds=seconds, status="ok" if p.returncode == 0 else "failed",
                            expect="ok")
            for key, value in provenance.items():
                h5.attrs["provenance." + key] = value
            h5.create_dataset("log/stdout", data=p.stdout)
            h5.create_dataset("log/stderr", data=p.stderr)
            files = h5.create_group("files")
            for path in sorted(Path(tmp).glob("*.txt")):
                files.create_dataset(path.name, data=np.loadtxt(path, ndmin=2))
    return "components", seconds, "ok" if p.returncode == 0 else "failed", "ok"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--maxent", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--sets", default="fast,targeted,cli",
                    help="comma-separated: full,fast,targeted,cli,kk,legendre_convert")
    ap.add_argument("--cases", default=".*", help="regular expression on case names")
    ap.add_argument("--components", type=Path, help="dump_components binary (adds components.h5)")
    ap.add_argument("--kk", type=Path, help="kk binary (needed for the 'kk' set)")
    ap.add_argument("--legendre-convert", type=Path,
                    help="legendre_convert binary (needed for the 'legendre_convert' set)")
    ap.add_argument("--provenance", type=Path, help="JSON file with provenance key/values")
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--timeout", type=float, default=TIMEOUT, help="seconds per case (default %(default)s)")
    args = ap.parse_args()

    sets = set(args.sets.split(","))
    selected = [c for c in all_cases() if c["set"] in sets and re.search(args.cases, c["name"])]
    provenance = json.loads(args.provenance.read_text()) if args.provenance else {}
    args.out.mkdir(parents=True, exist_ok=True)
    programs = {"maxent": args.maxent.resolve()}
    if args.kk:
        programs["kk"] = args.kk.resolve()
    if args.legendre_convert:
        programs["legendre_convert"] = args.legendre_convert.resolve()
    missing = sorted({c["program"] for c in selected} - set(programs))
    if missing:
        ap.error(f"selected cases need --{missing[0].replace('_', '-')}")

    results = []
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = [pool.submit(run_case, c, programs, args.out, provenance, args.timeout)
                   for c in selected]
        if args.components:
            futures.append(pool.submit(run_components, args.components.resolve(), args.out, provenance))
        for f in futures:
            results.append(f.result())
    bad = 0
    for name, seconds, status, expect in sorted(results):
        note = "" if status == expect else f"  (expected {expect})"
        print(f"{name:40s} {seconds:8.2f} s  {status}{note}")
        bad += status != expect
    print(f"{len(results)} cases, {bad} unexpected")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
