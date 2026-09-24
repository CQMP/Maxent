#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Generate the synthetic inputs for the targeted regression cases.

The inputs are made from spectral functions A(omega) that are known exactly,
transformed to the data space with the kernel conventions of
src/maxent_kernel.cpp, and given seeded Gaussian noise.  The output is
committed; the references never depend on re-running this script.

    python3 test/regression/make_inputs.py   # writes test/regression/inputs/

Spectra (all normalized to 1 on the support Maxent sees):
  S1  Gaussian at 0, sigma 1                        (particle-hole symmetric)
  S2  0.6 N(-1.5, 0.5) + 0.4 N(2.5, 0.5)            (not particle-hole symmetric)
  S3  2 N(0, 1.5) restricted to omega >= 0          (bosonic, OMEGA_MIN=0)
  S4  N(2.5, 0.5) on omega >= 0 (T=0, OMEGA_MIN=0)

The kk inputs come from the maxent references (Im parts) or analytically from
a Gaussian spectrum (Re part, needs scipy for the Dawson function).

Cases with OMEGA_MIN=0 use FREQUENCY_GRID=linear: the default Lorentzian grid is
centered at (OMEGA_MIN+OMEGA_MAX)/2 and has almost no points near omega=0.
"""

from pathlib import Path

import h5py
import numpy as np
from numpy.polynomial.legendre import leggauss, legval

HERE = Path(__file__).resolve().parent
OUT = HERE / "inputs"
SEED = 20260923
NOISE = 1e-4
BETA = 5.0


# ---------------------------------------------------------------- spectra
def gauss(w, mu, s):
    return np.exp(-0.5 * ((w - mu) / s) ** 2) / (np.sqrt(2 * np.pi) * s)


def omega_quadrature(lo, hi, panels=400, order=20):
    x, wt = leggauss(order)
    edges = np.linspace(lo, hi, panels + 1)
    a, b = edges[:-1, None], edges[1:, None]
    nodes = 0.5 * (b - a) * x[None, :] + 0.5 * (a + b)
    weights = 0.5 * (b - a) * wt[None, :]
    return nodes.ravel(), weights.ravel()


def spectrum(name):
    """Return quadrature nodes, weights and A(nodes) for a named spectrum."""
    if name == "S1":
        w, q = omega_quadrature(-10, 10)
        a = gauss(w, 0.0, 1.0)
    elif name == "S2":
        w, q = omega_quadrature(-10, 10)
        a = 0.6 * gauss(w, -1.5, 0.5) + 0.4 * gauss(w, 2.5, 0.5)
    elif name == "S3":
        w, q = omega_quadrature(0, 10)
        a = 2 * gauss(w, 0.0, 1.5)
    elif name == "S4":
        w, q = omega_quadrature(0, 10)
        a = gauss(w, 2.5, 0.5)
    else:
        raise ValueError(name)
    return w, q, a / np.sum(q * a)


# ---------------------------------------------------------------- kernels
def matsubara_fermionic(n):
    return (2 * np.arange(n) + 1) * np.pi / BETA


def k_time_fermionic(tau, w):
    # -exp(-w tau) / (1 + exp(-beta w)), evaluated stably
    t, w = tau[:, None], w[None, :]
    return -np.exp(-w * t - np.logaddexp(0.0, -BETA * w))


def k_time_bosonic(tau, w):
    # 0.5 w (exp(-w tau) + exp(-w (beta - tau))) / (1 - exp(-beta w)); w -> 0 limit is 1/beta
    t, w = tau[:, None], w[None, :]
    with np.errstate(divide="ignore", invalid="ignore"):
        k = 0.5 * w * (np.exp(-w * t) + np.exp(-w * (BETA - t))) / (-np.expm1(-BETA * w))
    return np.where(w == 0, 1.0 / BETA, k)


def k_tzero(tau, w):
    return -np.exp(-tau[:, None] * w[None, :])


def k_legendre(lmax, w, tau_order=200):
    x, wt = leggauss(tau_order)
    tau = 0.5 * BETA * (x + 1)
    dtau = 0.5 * BETA * wt
    k = np.empty((lmax, w.size))
    kernel = np.exp(-tau[:, None] * w[None, :] - np.logaddexp(0.0, -BETA * w)[None, :])
    for l in range(lmax):
        c = np.zeros(l + 1)
        c[l] = 1
        pl = legval(2 * tau / BETA - 1, c)
        k[l] = -np.sqrt(2 * l + 1) * (dtau * pl) @ kernel
    return k


def transform(kernel, a, q):
    return kernel @ (q * a)


# ---------------------------------------------------------------- writers
def rng_for(name):
    return np.random.default_rng([SEED, sum(map(ord, name))])


def fmt(x):
    return f"{x:.17g}"


def write_param(case_dir, entries):
    lines = [f"{k}={v}" for k, v in entries]
    (case_dir / "case.param").write_text("\n".join(lines) + "\n")


def write_columns(path, *cols):
    with open(path, "w") as f:
        for row in zip(*cols):
            f.write(" ".join(fmt(v) for v in row) + "\n")


def noisy(name, y, sigma=NOISE):
    return y + sigma * rng_for(name).standard_normal(y.shape)


def correlated_covariance(n, sigma=NOISE, rho=0.5, n_singular=2):
    """sigma^2 rho^|i-j| with the n_singular smallest eigenvalues set to 1e-13."""
    idx = np.arange(n)
    c = sigma ** 2 * rho ** np.abs(idx[:, None] - idx[None, :])
    lam, v = np.linalg.eigh(c)
    lam[:n_singular] = 1e-13
    return (v * lam) @ v.T


def correlated_noise(name, y, cov):
    lam, v = np.linalg.eigh(cov)
    z = rng_for(name).standard_normal(y.shape)
    return y + v @ (np.sqrt(np.clip(lam, 0, None)) * z)


COMMON = [("BETA", BETA), ("NFREQ", 200), ("N_ALPHA", 20)]


def freq_ph_data(name, n=24):
    w, q, a = spectrum("S1")
    wn = matsubara_fermionic(n)
    k = -wn[:, None] / (wn[:, None] ** 2 + w[None, :] ** 2)
    return wn, noisy(name, transform(k, a, q))


def freq_nonph_data(name, n=16, kernel="fermionic"):
    w, q, a = spectrum("S2")
    iwn = 1j * matsubara_fermionic(n)
    if kernel == "fermionic":
        kc = 1.0 / (iwn[:, None] - w[None, :])
    else:  # anomalous
        kc = -w[None, :] / (iwn[:, None] - w[None, :])
    g = transform(kc, a, q)
    rng = rng_for(name)
    re = g.real + NOISE * rng.standard_normal(n)
    im = g.imag + NOISE * rng.standard_normal(n)
    return iwn.imag, re, im


def base_frequency_ph(case_dir, name, extra=()):
    wn, y = freq_ph_data(name)
    write_columns(case_dir / "data.dat", wn, y, np.full_like(y, NOISE))
    write_param(case_dir, COMMON + [("NDAT", wn.size), ("DATASPACE", "frequency"),
                                    ("KERNEL", "fermionic"), ("PARTICLE_HOLE_SYMMETRY", "true"),
                                    ("DATA", '"data.dat"')] + list(extra))


def base_tzero(case_dir, name, extra=()):
    w, q, a = spectrum("S4")
    tau = np.linspace(0, 10, 41)
    y = noisy(name, transform(k_tzero(tau, w), a, q))
    write_columns(case_dir / "data.dat", tau, y, np.full_like(y, NOISE))
    write_param(case_dir, COMMON + [("NDAT", tau.size), ("DATASPACE", "time"), ("KERNEL", "tzero"),
                                    ("OMEGA_MIN", 0), ("OMEGA_MAX", 10), ("FREQUENCY_GRID", "linear"),
                                    ("DATA", '"data.dat"')] + list(extra))


# ---------------------------------------------------------------- cases
def make(name, case_dir):
    if name == "t_cov_text":
        w, q, a = spectrum("S2")
        tau = np.linspace(0, BETA, 33)
        y0 = transform(k_time_fermionic(tau, w), a, q)
        cov = correlated_covariance(tau.size)
        y = correlated_noise(name, y0, cov)
        write_columns(case_dir / "data.dat", tau, y, np.sqrt(np.diag(cov)))
        with open(case_dir / "cov.dat", "w") as f:
            for i in range(tau.size):
                for j in range(tau.size):
                    f.write(f"{i} {j} {fmt(cov[i, j])}\n")
        write_param(case_dir, COMMON + [("NDAT", tau.size), ("DATASPACE", "time"),
                                        ("KERNEL", "fermionic"), ("DATA", '"data.dat"'),
                                        ("COVARIANCE_MATRIX", '"cov.dat"')])
    elif name == "t_cov_hdf5":
        w, q, a = spectrum("S1")
        wn = matsubara_fermionic(24)
        y0 = transform(-wn[:, None] / (wn[:, None] ** 2 + w[None, :] ** 2), a, q)
        cov = correlated_covariance(wn.size)
        y = correlated_noise(name, y0, cov)
        with h5py.File(case_dir / "data.h5", "w") as f:
            f["/Data"] = y
            f["/Covariance"] = cov.ravel()  # row-major, i*ndat+j
        write_param(case_dir, COMMON + [("NDAT", wn.size), ("DATASPACE", "frequency"),
                                        ("KERNEL", "fermionic"), ("PARTICLE_HOLE_SYMMETRY", "true"),
                                        ("DATA_IN_HDF5", "true"), ("DATA", '"data.h5"'),
                                        ("COVARIANCE_MATRIX", '"in-hdf5"')])
    elif name == "t_hdf5_errors":
        _, re, im = freq_nonph_data(name)
        data = np.column_stack([re, im]).ravel()  # re0 im0 re1 im1 ...
        with h5py.File(case_dir / "data.h5", "w") as f:
            f["/Data"] = data
            f["/Error"] = np.full_like(data, NOISE)
        write_param(case_dir, COMMON + [("NDAT", data.size), ("DATASPACE", "frequency"),
                                        ("KERNEL", "fermionic"), ("PARTICLE_HOLE_SYMMETRY", "false"),
                                        ("DATA_IN_HDF5", "true"), ("DATA", '"data.h5"')])
    elif name == "t_param_xi":
        wn, y = freq_ph_data(name, n=16)
        entries = COMMON + [("NDAT", wn.size), ("DATASPACE", "frequency"), ("KERNEL", "fermionic"),
                            ("PARTICLE_HOLE_SYMMETRY", "true")]
        entries += [(f"X_{i}", fmt(v)) for i, v in enumerate(y)]
        entries += [(f"SIGMA_{i}", fmt(NOISE)) for i in range(y.size)]
        write_param(case_dir, entries)
    elif name == "t_param_taui":
        w, q, a = spectrum("S2")
        tau = np.linspace(0, BETA, 21)
        y = noisy(name, transform(k_time_fermionic(tau, w), a, q))
        entries = COMMON + [("NDAT", tau.size), ("DATASPACE", "time"), ("KERNEL", "fermionic")]
        entries += [(f"X_{i}", fmt(v)) for i, v in enumerate(y)]
        entries += [(f"SIGMA_{i}", fmt(NOISE)) for i in range(y.size)]
        entries += [(f"TAU_{i}", fmt(t)) for i, t in enumerate(tau)]
        write_param(case_dir, entries)
    elif name == "t_kernel_tzero":
        base_tzero(case_dir, name)
    elif name == "t_kernel_time_bosonic":
        w, q, a = spectrum("S3")
        tau = np.linspace(0, BETA, 33)
        y = noisy(name, transform(k_time_bosonic(tau, w), a, q))
        write_columns(case_dir / "data.dat", tau, y, np.full_like(y, NOISE))
        write_param(case_dir, COMMON + [("NDAT", tau.size), ("DATASPACE", "time"), ("KERNEL", "bosonic"),
                                        ("OMEGA_MIN", 0), ("OMEGA_MAX", 10), ("FREQUENCY_GRID", "linear"),
                                        ("DATA", '"data.dat"')])
    elif name == "t_kernel_anomalous_ph":
        w, q, a = spectrum("S1")
        wn = matsubara_fermionic(24)
        y = noisy(name, transform(w[None, :] ** 2 / (wn[:, None] ** 2 + w[None, :] ** 2), a, q))
        write_columns(case_dir / "data.dat", wn, y, np.full_like(y, NOISE))
        write_param(case_dir, COMMON + [("NDAT", wn.size), ("DATASPACE", "frequency"),
                                        ("KERNEL", "anomalous"), ("PARTICLE_HOLE_SYMMETRY", "true"),
                                        ("DATA", '"data.dat"')])
    elif name == "t_kernel_anomalous_nonph":
        wn, re, im = freq_nonph_data(name, kernel="anomalous")
        err = np.full_like(re, NOISE)
        write_columns(case_dir / "data.dat", wn, re, err, im, err)
        write_param(case_dir, COMMON + [("NDAT", 2 * wn.size), ("DATASPACE", "frequency"),
                                        ("KERNEL", "anomalous"), ("PARTICLE_HOLE_SYMMETRY", "false"),
                                        ("DATA", '"data.dat"')])
    elif name == "t_kernel_legendre_bosonic":
        w, q, a = spectrum("S1")
        lmax = 12
        y = noisy(name, transform(k_legendre(lmax, w), a, q))
        write_columns(case_dir / "data.dat", np.arange(lmax), y, np.full_like(y, NOISE))
        write_param(case_dir, COMMON + [("NDAT", lmax), ("DATASPACE", "legendre"), ("KERNEL", "bosonic"),
                                        ("PARTICLE_HOLE_SYMMETRY", "true"), ("DATA", '"data.dat"')])
    elif name.startswith("t_model_"):
        model = name[len("t_model_"):]
        options = {
            "flat": [("DEFAULT_MODEL", "flat")],
            "gaussian": [("DEFAULT_MODEL", "gaussian"), ("SIGMA", 1.5)],
            "shifted_gaussian": [("DEFAULT_MODEL", '"shifted gaussian"'), ("SIGMA", 1.5), ("SHIFT", 0.5)],
            "double_gaussian": [("DEFAULT_MODEL", '"double gaussian"'), ("SIGMA", 0.7), ("SHIFT", 1)],
            "two_gaussians": [("DEFAULT_MODEL", '"two gaussians"'), ("SIGMA1", 0.5), ("SIGMA2", 1),
                              ("SHIFT1", -1), ("SHIFT2", 1), ("NORM1", 0.4)],
            "general_double_gaussian": [("DEFAULT_MODEL", '"general double gaussian"'), ("SIGMA", 1),
                                        ("SHIFT", 1), ("BOSE_NORM", 0.5)],
            "lorentzian": [("DEFAULT_MODEL", "lorentzian"), ("GAMMA", 1)],
            "shifted_lorentzian": [("DEFAULT_MODEL", '"shifted lorentzian"'), ("GAMMA", 1), ("SHIFT", 0.5)],
            "double_lorentzian": [("DEFAULT_MODEL", '"double lorentzian"'), ("GAMMA", 0.7), ("SHIFT", 1)],
            "two_lorentzians": [("DEFAULT_MODEL", '"two lorentzians"'), ("GAMMA1", 0.5), ("GAMMA2", 1),
                                ("SHIFT1", -1), ("SHIFT2", 1)],
            "linear_rise_exp_decay": [("DEFAULT_MODEL", '"linear rise exp decay"'), ("LAMBDA", 1)],
            # LAMBDA=2: with LAMBDA=1 the minimizer diverges on this data (B19)
            "quadratic_rise_exp_decay": [("DEFAULT_MODEL", '"quadratic rise exp decay"'), ("LAMBDA", 2)],
            "tabulated": [("DEFAULT_MODEL", '"model.dat"')],
            "runs": [("MODEL_RUNS", 2), ("RUN_0", "flat"), ("RUN_1", "gaussian"), ("SIGMA", 1.5)],
        }[model]
        if model == "tabulated":
            wt = np.linspace(-10, 10, 201)
            write_columns(case_dir / "model.dat", wt, gauss(wt, 0.0, 1.5))
        if model.endswith("rise_exp_decay"):
            base_tzero(case_dir, name, options)
        else:
            base_frequency_ph(case_dir, name, options)
    elif name.startswith("t_grid_"):
        grid = {"lorentzian": "lorentzian", "half_lorentzian": "half-lorentzian",
                "quadratic": "quadratic", "log": "log", "linear": "linear"}[name[len("t_grid_"):]]
        base_frequency_ph(case_dir, name, [("FREQUENCY_GRID", grid)])
    elif name in ("t_generate_err", "t_generate_err_seed"):
        # Keep the scientific input identical so the two cases isolate SEED.
        wn, y = freq_ph_data("t_generate_err")
        write_columns(case_dir / "data.dat", wn, y, np.full_like(y, NOISE))
        entries = [("BETA", BETA), ("NFREQ", 40), ("N_ALPHA", 8),
                   ("NDAT", wn.size), ("DATASPACE", "frequency"),
                   ("KERNEL", "fermionic"), ("PARTICLE_HOLE_SYMMETRY", "true"),
                   ("DATA", '"data.dat"'), ("GENERATE_ERR", "true")]
        if name == "t_generate_err_seed":
            entries.append(("SEED", 1234))
        write_param(case_dir, entries)
    elif name in ("kk_imag_to_real_green", "kk_imag_to_real_self"):
        source, dataset, scale = {
            "kk_imag_to_real_green": ("u0_frequency", "files/in.out.avspec.dat", -np.pi),
            "kk_imag_to_real_self": ("self_u1", "files/in.out.avspec_self.dat", 1.0),
        }[name]
        with h5py.File(HERE / "reference" / f"{source}.h5", "r") as f:
            spec = f[dataset][()]
        write_columns(case_dir / "input.dat", spec[:, 0], scale * spec[:, 1])
    elif name == "kk_real_to_imag":
        from scipy.special import dawsn
        w = np.linspace(-10, 10, 401)
        # Re G of a normalized Gaussian spectrum (sigma 1): sqrt(2)/sigma D(w/(sqrt(2) sigma))
        write_columns(case_dir / "input.dat", w, np.sqrt(2) * dawsn(w / np.sqrt(2)))
    elif name == "legendre_convert_transform":
        # Zero error bars make the clock-seeded bootstrap deterministic.  The
        # data remain nontrivial and exercise all Legendre orders in the test.
        w, q, a = spectrum("S1")
        tau = np.linspace(0, BETA, 33)
        gtau = transform(k_time_fermionic(tau, w), a, q)
        write_columns(case_dir / "input.dat", tau, gtau, np.zeros_like(gtau))
    elif name == "cli_missing_beta":
        write_param(case_dir, [("NDAT", 4), ("X_0", 0.1), ("X_1", 0.2), ("X_2", 0.3), ("X_3", 0.4),
                               ("SIGMA_0", 0.5), ("SIGMA_1", 0.5), ("SIGMA_2", 0.5), ("SIGMA_3", 0.5)])
    else:
        raise ValueError(f"no input recipe for {name}")


def main():
    import sys
    sys.path.insert(0, str(HERE))
    from cases import all_cases
    names = [c["name"] for c in all_cases()
             if c["inputs"] and c["inputs"].startswith("test/regression/inputs/")]
    for name in names:
        case_dir = OUT / name
        case_dir.mkdir(parents=True, exist_ok=True)
        for old in case_dir.iterdir():
            old.unlink()
        make(name, case_dir)
        print(f"wrote {case_dir.relative_to(HERE.parent.parent)}: "
              f"{', '.join(sorted(p.name for p in case_dir.iterdir()))}")


if __name__ == "__main__":
    main()
