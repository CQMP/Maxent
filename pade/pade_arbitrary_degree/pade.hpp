/*****************************************************************************
 *
 * ALPS Project Applications
 *
 * Copyright (C) 2013 by Emanuel Gull <egull@umich.edu>
 *
 * This software is part of the ALPS Applications, published under the ALPS
 * Application License; you can use, redistribute it and/or modify it under
 * the terms of the license, either version 1 or (at your option) any later
 * version.
 *
 * You should have received a copy of the ALPS Application License along with
 * the ALPS Applications; see the file LICENSE.txt. If not, the license is also
 * available from http://alps.comp-phys.org/.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
 * SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
 * FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 *
 *****************************************************************************/

#ifndef ALPS_TOOL_PADE_HPP
#define ALPS_TOOL_PADE_HPP
#include<alps/ngs.hpp>
#include<alps/numeric/matrix.hpp>
#include <gmpxx.h>
typedef std::vector<double> vector_type;
typedef std::vector<std::complex<double> > complex_vector_type;

//typedefs:
//typedef double pade_real_type;
typedef mpf_class pade_real_type;
typedef std::complex<pade_real_type> pade_complex_type;
typedef std::vector<pade_real_type> pade_vector_type;
typedef std::vector<pade_complex_type> pade_complex_vector_type;
typedef alps::numeric::matrix<pade_complex_type> pade_complex_matrix_type;

enum frequency_grid_type{
  lorentzian_grid,
  half_lorentzian_grid,
  quadratic_grid,
  log_grid,
  linear_grid
};
enum imaginary_domain{
  tau,
  omegan
};

//this class contains the real frequency discretization. A typical case is a logarithmic grid with many points near zero and few at high frequencies
class grid{
public:
  //constructor
  grid(const alps::params &p);
  //grid frequency points
  const vector_type &freq() const{ return freq_; }
  //delta_freq contains the half the distance between points i-1 and i+1
  const vector_type &delta_freq() const{ return delta_freq_; }
  //mapping from [0,1] to [omega_min, omega_max]
  double omega_of_t(const double t) const { return omega_min_ + (omega_max_-omega_min_)*t; }
  //mapping from [omega_min, omega_max] to [0,1]
  double t_of_omega(const double omega) const { return (omega-omega_min_)/(omega_max_-omega_min_); }
private:
  //construct a grid
  void setup_grid(const alps::params &p);
  //number of frequencies
  int N_freq_;
  //type of the grid
  frequency_grid_type grid_type_;
  //the integrated relative weight
  vector_type t_array_;
  //the actual frequency grid
  vector_type freq_;
  //delta_freq contains the half the distance between points i-1 and i+1
  vector_type delta_freq_;               
  
  //the minimum and maximum value of the grid domain
  double omega_min_;
  double omega_max_;
};

class imag_domain_grid{
public:
  imag_domain_grid(const alps::params &p);
  const double &freq(int i) const{return freq_[i];}
  const vector_type &freq() const{return freq_;}
private:
  //number of frequencies
  int N_freq_;
  //values of frequencies
  vector_type freq_;
  //inverse temperature
  double  T_;
};

//this class contains the input raw data which is stored in Matsubara frequency or imaginary time space
class imaginary_domain_data{
public:
  //constructor
  imaginary_domain_data(const alps::params &p);
  //get the data of frequency/time i
  const std::complex<double> &operator()(int i) const{ return val_[i];}
  //number of data points
  int N_imag() const{return N_imag_;}
  //return normalization of data
  double norm() const{return norm_;}
  //vector of data values, mutable access
  complex_vector_type &val(){return val_;}
  //vector of data values, immutable access
  const complex_vector_type &val() const {return val_;}
  //vector of data x points, immutable access
  const vector_type &freq() const {return G_.freq();}
  //text i/o function for the input data
  void write(const std::string &filename) const;
private:
  //Matsubara grid
  imag_domain_grid G_;
  //Matsubara data
  complex_vector_type val_;
  //normalization of the spectral function: (target) integral over the real axis of the continued function.
  double norm_;
  //number of input Matsubara data points.
  int N_imag_;
};

class real_domain_data{
public:
  //constructor
  real_domain_data(const alps::params &p);
  //number of real frequencies
  int N_real() const{return N_real_;}
  //function values of the continued function, immutable access
  const complex_vector_type &val() const{return val_;}
  //function values of the continued function, mutable access
  complex_vector_type &val() {return val_;}
  //evaluation point of the continued function, immutable access
  const vector_type &freq() const {return G_.freq();}

  void write(const std::string &filename) const;
private:
  const class grid G_;
  complex_vector_type val_;  //spectral function
  //number of real frequency discretization points.
  int N_real_;
};

class pade_interpolator{
public:
  //constructor
  pade_interpolator(const alps::params &p);
  //interpolation routine
  void pade_interpolate(const imaginary_domain_data &data, real_domain_data &real) const;
  
private:
  
  //private functions
  void find_epsilon();
  //evaluate a rational function in its barycentric form
  pade_complex_type evaluate_bary_poly(const pade_complex_vector_type &q, const pade_complex_vector_type &f, const pade_complex_vector_type &x, const pade_complex_matrix_type &Vinv, int m, int n, const pade_complex_type &x0)const;
  //see Eq. 3 for omega_j
  pade_complex_type compute_omega_j(int j, const pade_complex_vector_type &s, const pade_complex_type &x0)const;
  pade_complex_type compute_omegaprime_kp1_of_xj(int k, int j, const pade_complex_vector_type &s) const;
  //see Eq. 5 for the inverse of V
  void fill_Vinv_matrix(pade_complex_matrix_type &Vinv, const pade_complex_vector_type &s, int m, int n)const;
  //see Eq. 9, use q0=1, and build a matrix system
  void assemble_matrix_system(pade_complex_matrix_type &Lambda, pade_complex_vector_type &rhs, const pade_complex_matrix_type &Vinv, const pade_complex_vector_type &f, int m, int n)const;

  std::complex<double> to_simple_precision(const std::complex<mpf_class> &x) const{ return std::complex<double>(x.real().get_d(),x.imag().get_d()); }
  std::complex<double> to_simple_precision(const std::complex<double> &x) const{ return x; }
  //private variables:
  //polynomial degree (numerator)
  int pade_mu_;
  //polynomial degree (denominator)
  int pade_nu_;
  //floating point epsilon
  pade_real_type epsilon_;
};

 class pade_solver{
public:
  pade_solver(){}
  void backsub_lower(const pade_complex_matrix_type &Linv, const pade_complex_vector_type &rhs, pade_complex_vector_type &res);
  void backsub_upper(const pade_complex_matrix_type &Linv, const pade_complex_vector_type &rhs, pade_complex_vector_type &res);
  void solve(const pade_complex_matrix_type &A, const pade_complex_vector_type &rhs, pade_complex_vector_type &res);
};

#endif
