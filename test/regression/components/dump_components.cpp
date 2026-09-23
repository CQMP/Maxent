/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1998-2026 ALPS Collaboration
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

// Dumps Maxent's building blocks at small sizes for the regression suite:
// kernel matrices for every reachable kernel type, real-frequency grids,
// default models, and the singular values of the scaled kernel.
//
//   dump_components <output directory>
//
// Every quantity is written to <name>.txt as whitespace-separated rows with
// 17 significant digits; generate.py packs them into components.h5.

#include "maxent.hpp"
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>

namespace {

std::string outdir;

void write_matrix(const std::string &name, const matrix_type &m) {
  std::ofstream os((outdir + "/" + name + ".txt").c_str());
  os << std::setprecision(17);
  for (int i = 0; i < m.rows(); ++i) {
    for (int j = 0; j < m.cols(); ++j)
      os << (j ? " " : "") << m(i, j);
    os << "\n";
  }
}

void write_columns(const std::string &name, const std::vector<std::vector<double> > &cols) {
  matrix_type m(cols[0].size(), cols.size());
  for (std::size_t j = 0; j < cols.size(); ++j)
    for (std::size_t i = 0; i < cols[j].size(); ++i)
      m(i, j) = cols[j][i];
  write_matrix(name, m);
}

alps::params base_params() {
  alps::params p;
  MaxEntSimulation::define_parameters(p);
  p["BETA"] = 5.0;
  p["TEXT_OUTPUT"] = false;
  return p;
}

// ---------------------------------------------------------------- kernels
void dump_kernel(const std::string &name, const std::string &dataspace, const std::string &kernel_name,
                 bool ph, int ndat, double omega_min, double omega_max, bool time_input) {
  alps::params p = base_params();
  const int nfreq = 24;
  p["NDAT"] = ndat;
  p["NFREQ"] = nfreq;
  p["DATASPACE"] = dataspace;
  p["KERNEL"] = kernel_name;
  p["PARTICLE_HOLE_SYMMETRY"] = ph;
  vector_type freq(nfreq);
  for (int j = 0; j < nfreq; ++j)
    freq(j) = omega_min + (omega_max - omega_min) * (j + 0.5) / nfreq;
  vector_type input_grid = vector_type::Zero(ndat);
  if (time_input)
    for (int i = 0; i < ndat; ++i)
      input_grid(i) = 5.0 * i / (ndat - 1);
  kernel k(p, freq, input_grid);
  write_matrix("kernel_" + name, k());
  matrix_type grid_out(ndat, 1);
  grid_out.col(0) = input_grid;
  write_matrix("kernel_" + name + "_inputgrid", grid_out);
}

// ---------------------------------------------------------------- grids
void dump_grid(const std::string &name, const std::string &grid_name, int nfreq) {
  alps::params p = base_params();
  p["NFREQ"] = nfreq;
  p["FREQUENCY_GRID"] = grid_name;
  grid g(p);
  std::vector<double> t(g.t_array());
  write_columns("grid_" + name + "_" + std::to_string(nfreq), std::vector<std::vector<double> >(1, t));
}

// ---------------------------------------------------------------- default models
void dump_model(const std::string &name, alps::params p, double omega_min, double omega_max) {
  p["OMEGA_MIN"] = omega_min;
  p["OMEGA_MAX"] = omega_max;
  boost::shared_ptr<DefaultModel> model = make_default_model(p, "DEFAULT_MODEL");
  std::vector<double> omega, d, x, omega_of_x;
  for (int i = 0; i < 41; ++i) {
    omega.push_back(omega_min + (omega_max - omega_min) * (i + 0.5) / 41);
    d.push_back(model->D(omega.back()));
    x.push_back((i + 0.5) / 41);
    omega_of_x.push_back(model->omega(x.back()));
  }
  std::vector<std::vector<double> > cols;
  cols.push_back(omega);
  cols.push_back(d);
  cols.push_back(x);
  cols.push_back(omega_of_x);
  write_columns("model_" + name, cols);
}

alps::params model_params(const std::string &model) {
  alps::params p = base_params();
  p["DEFAULT_MODEL"] = model;
  return p;
}

// ---------------------------------------------------------------- singular values
void dump_svd(const std::string &name, alps::params p) {
  MaxEntParameters mp(p);
  matrix_type s(mp.ns(), 1);
  for (int i = 0; i < mp.ns(); ++i)
    s(i, 0) = mp.Sigma()(i, i);
  write_matrix("svd_" + name, s);
}

}  // namespace

int main(int argc, char **argv) {
  if (argc != 2) {
    std::fprintf(stderr, "usage: %s <output directory>\n", argv[0]);
    return 2;
  }
  outdir = argv[1];

  // kernels: every kernel type reachable through set_kernel_type (the time-Legendre types are not, B5)
  dump_kernel("time_fermionic", "time", "fermionic", false, 9, -6, 6, true);
  dump_kernel("time_bosonic", "time", "bosonic", false, 9, 0, 6, true);
  dump_kernel("time_tzero", "time", "tzero", false, 9, 0, 6, true);
  dump_kernel("legendre_fermionic", "legendre", "fermionic", true, 6, -6, 6, false);
  dump_kernel("legendre_bosonic", "legendre", "bosonic", true, 6, -6, 6, false);
  dump_kernel("frequency_fermionic_ph", "frequency", "fermionic", true, 8, -6, 6, false);
  dump_kernel("frequency_bosonic_ph", "frequency", "bosonic", true, 8, -6, 6, false);
  dump_kernel("frequency_anomalous_ph", "frequency", "anomalous", true, 8, -6, 6, false);
  dump_kernel("frequency_fermionic", "frequency", "fermionic", false, 8, -6, 6, false);
  dump_kernel("frequency_bosonic", "frequency", "bosonic", false, 8, -6, 6, false);
  dump_kernel("frequency_anomalous", "frequency", "anomalous", false, 8, -6, 6, false);

  // grids: odd and even sizes
  const char *grids[][2] = {{"lorentzian", "lorentzian"}, {"half_lorentzian", "half lorentzian"},
                            {"quadratic", "quadratic"}, {"log", "log"}, {"linear", "linear"}};
  for (auto &g : grids) {
    dump_grid(g[0], g[1], 20);
    dump_grid(g[0], g[1], 21);
  }

  // default models (parameters as in the targeted regression cases)
  {
    std::ofstream tab("tab_model.dat");
    tab << std::setprecision(17);
    for (int i = 0; i <= 200; ++i) {
      double w = -10 + 0.1 * i;
      tab << w << " " << std::exp(-w * w / 4.5) / std::sqrt(2 * M_PI * 2.25) << "\n";
    }
  }
  alps::params p;
  dump_model("flat", model_params("flat"), -10, 10);
  p = model_params("gaussian"); p["SIGMA"] = 1.5; dump_model("gaussian", p, -10, 10);
  p = model_params("shifted gaussian"); p["SIGMA"] = 1.5; p["SHIFT"] = 0.5; dump_model("shifted_gaussian", p, -10, 10);
  p = model_params("double gaussian"); p["SIGMA"] = 0.7; p["SHIFT"] = 1.0; dump_model("double_gaussian", p, -10, 10);
  p = model_params("two gaussians"); p["SIGMA1"] = 0.5; p["SIGMA2"] = 1.0; p["SHIFT1"] = -1.0; p["SHIFT2"] = 1.0;
  p["NORM1"] = 0.4; dump_model("two_gaussians", p, -10, 10);
  p = model_params("general double gaussian"); p["SIGMA"] = 1.0; p["SHIFT"] = 1.0; p["BOSE_NORM"] = 0.5;
  dump_model("general_double_gaussian", p, -10, 10);
  p = model_params("lorentzian"); p["GAMMA"] = 1.0; dump_model("lorentzian", p, -10, 10);
  p = model_params("shifted lorentzian"); p["GAMMA"] = 1.0; p["SHIFT"] = 0.5; dump_model("shifted_lorentzian", p, -10, 10);
  p = model_params("double lorentzian"); p["GAMMA"] = 0.7; p["SHIFT"] = 1.0; dump_model("double_lorentzian", p, -10, 10);
  p = model_params("two lorentzians"); p["GAMMA1"] = 0.5; p["GAMMA2"] = 1.0; p["SHIFT1"] = -1.0; p["SHIFT2"] = 1.0;
  dump_model("two_lorentzians", p, -10, 10);
  p = model_params("linear rise exp decay"); p["LAMBDA"] = 1.0; dump_model("linear_rise_exp_decay", p, 0, 10);
  p = model_params("quadratic rise exp decay"); p["LAMBDA"] = 1.0; dump_model("quadratic_rise_exp_decay", p, 0, 10);
  p = model_params("tab_model.dat"); dump_model("tabulated", p, -10, 10);

  // singular values of the error-scaled kernel (A = two delta peaks at +-1, sigma = 1e-4)
  {
    alps::params q = base_params();
    const int ndat = 16;
    q["NDAT"] = ndat; q["NFREQ"] = 100; q["DATASPACE"] = "frequency"; q["KERNEL"] = "fermionic";
    q["PARTICLE_HOLE_SYMMETRY"] = true;
    for (int i = 0; i < ndat; ++i) {
      double wn = (2 * i + 1) * M_PI / 5.0;
      q["X_" + std::to_string(i)] = -wn / (wn * wn + 1.0);
      q["SIGMA_" + std::to_string(i)] = 1e-4;
    }
    dump_svd("frequency_fermionic_ph", q);
  }
  {
    alps::params q = base_params();
    const int ndat = 17;
    q["NDAT"] = ndat; q["NFREQ"] = 100; q["DATASPACE"] = "time"; q["KERNEL"] = "fermionic";
    for (int i = 0; i < ndat; ++i) {
      double tau = 5.0 * i / (ndat - 1);
      q["X_" + std::to_string(i)] = -0.5 * (std::exp(-tau) / (1 + std::exp(-5.0)) + std::exp(tau) / (1 + std::exp(5.0)));
      q["SIGMA_" + std::to_string(i)] = 1e-4;
      q["TAU_" + std::to_string(i)] = tau;
    }
    dump_svd("time_fermionic", q);
  }
  return 0;
}
