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

#pragma once
#include<complex>
#include<vector>
#include<gmpxx.h>
#include<Eigen/Core>
#include "alps/params.hpp"
//#include<tgmath.h>

typedef std::vector<double> vector_type;
typedef std::vector<std::complex<double> > complex_vector_type;

//typedefs:
typedef mpf_class pade_real_type;
typedef std::complex<pade_real_type> pade_complex_type;
typedef std::vector<pade_real_type> pade_vector_type;
typedef std::vector<pade_complex_type> pade_complex_vector_type;

//missing arithmetics
inline bool isnan(pade_complex_type x){ return false;}
inline bool isinf(pade_complex_type x){ return false;}
#ifdef copysign
#undef copysign
#endif
inline pade_real_type copysign(const pade_real_type &a, const pade_real_type &b){
  return sgn(a)==sgn(b)?a:-1*a;
}
inline pade_complex_type operator/(const pade_complex_type &p, const pade_complex_type &q){
  pade_real_type a=p.real(), b=p.imag(), c=q.real(), d=q.imag();
  return std::complex<pade_real_type>((a*c+b*d)/(c*c+d*d), (b*c-a*d)/(c*c+d*d));
}
namespace Eigen {
template<> struct NumTraits<pade_complex_type>
 : NumTraits<double> // permits to get the epsilon, dummy_precision, lowest, highest functions
{
  typedef pade_real_type Real;
  typedef pade_complex_type NonInteger;
  typedef pade_complex_type Nested;
  enum {
    IsComplex = 1,
    IsInteger = 0,
    IsSigned = 1,
    RequireInitialization = 1,
    ReadCost = 1,
    AddCost = 3,
    MulCost = 3
  };
};
}

typedef Eigen::Matrix<pade_complex_type, Eigen::Dynamic, Eigen::Dynamic> pade_complex_matrix_type;

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

class PadeParams:public alps::params {
public:
  PadeParams(){
    define_parameters();
  }
  PadeParams(int argc, const char* argv[]):alps::params(argc, argv){
    define_parameters();
  }

private:
  void define_parameters(){
    define<int>("real.NFREQ", "Number of real frequency points");
    define<std::string>("real.FREQUENCY_GRID", "Type of real frequency grid: Lorentzian, half Lorentzian, quadratic, log, or linear");
    define<double>("real.CUT", 0.01, "Lorentzian cutoff parameter");
    define<double>("real.SPREAD", 4, "Quadratic grid spread parameter");
    define<double>("real.LOG_MIN", 1.0e-4, "Log grid minimum point parameter");
    define<double>("real.OMEGA_MIN", -25, "lowest frequency point");
    define<double>("real.OMEGA_MAX", 25, "highest frequency point");
    define<std::string>("real.OUTPUT", "output file format");
    
    define<int>("imag.NDAT", "number of input frequency points");
    define<double>("imag.BETA", "inverse temperature");
    define<std::string>("imag.STATISTICS", "Fermi or Bose statistics");
    define<std::string>("imag.DATA", "text input data file in the format \"freq real imag\"");
    define<bool>("imag.NEGATIVE_DATA", false, "set to true if data for both pos and neg frequencies");
    
    define<int>("pade.PADE_NUMERATOR_DEGREE", "Degree of pade numerator");
    define<int>("pade.PADE_DENOMINATOR_DEGREE", "Degree of pade numerator");
    define<int>("pade.FLOAT_PRECISION", 256, "Precision of floating point arithmetics");
   
    if (help_requested(std::cout)) {
      exit(0);
    }
  }
};

//this class contains the real frequency discretization. A typical case is a logarithmic grid with many points near zero and few at high frequencies
class grid{
public:
  //constructor
  grid(const PadeParams &p);
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
  void setup_grid(const PadeParams &p);
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
  imag_domain_grid(const PadeParams &p);
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
  imaginary_domain_data(const PadeParams &p);
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
  real_domain_data(const PadeParams &p);
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
  pade_interpolator(const PadeParams &p);
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
