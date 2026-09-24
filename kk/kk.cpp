/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1998-2016 ALPS Collaboration
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include<iostream>
#include<fstream>
#include<stdexcept>
#include<vector>
#include<complex>
#include<cmath>
#include<cstdlib>
//#include<alps/params.hpp>
#include<algorithm>
#include<sstream>
#include<string>
#include <boost/program_options.hpp>

/// Natural cubic spline, the same algorithm as GSL's gsl_interp_cspline (whose
/// results it reproduces): second derivatives vanish at both ends, and the
/// symmetric tridiagonal system is solved by the same LDL^T elimination.
class natural_cubic_spline {
public:
  natural_cubic_spline(const std::vector<double> &x, const std::vector<double> &y)
  : x_(x), y_(y), c_(x.size(), 0.) {
    const std::size_t n = x.size();
    if (n < 3 || y.size() != n) throw std::invalid_argument("spline needs at least 3 points and matching x and y");
    // like gsl_spline_init: the grid must be strictly increasing
    for (std::size_t i = 0; i + 1 < n; ++i)
      if (!(x[i + 1] > x[i])) throw std::invalid_argument("spline: x values must be strictly increasing");
    const std::size_t m = n - 2;  // interior points
    std::vector<double> diag(m), offdiag(m), rhs(m);
    for (std::size_t i = 0; i < m; ++i) {
      const double h_i = x[i + 1] - x[i], h_ip1 = x[i + 2] - x[i + 1];
      const double g_i = h_i != 0. ? 1. / h_i : 0., g_ip1 = h_ip1 != 0. ? 1. / h_ip1 : 0.;
      offdiag[i] = h_ip1;
      diag[i] = 2. * (h_ip1 + h_i);
      rhs[i] = 3. * ((y[i + 2] - y[i + 1]) * g_ip1 - (y[i + 1] - y[i]) * g_i);
    }
    if (m == 1) {
      c_[1] = rhs[0] / diag[0];
      return;
    }
    // LDL^T decomposition and solve (as gsl_linalg_solve_symm_tridiag)
    std::vector<double> gamma(m), alpha(m), z(m);
    alpha[0] = diag[0];
    gamma[0] = offdiag[0] / alpha[0];
    for (std::size_t i = 1; i < m - 1; ++i) {
      alpha[i] = diag[i] - offdiag[i - 1] * gamma[i - 1];
      gamma[i] = offdiag[i] / alpha[i];
    }
    alpha[m - 1] = diag[m - 1] - offdiag[m - 2] * gamma[m - 2];
    z[0] = rhs[0];
    for (std::size_t i = 1; i < m; ++i) z[i] = rhs[i] - gamma[i - 1] * z[i - 1];
    for (std::size_t i = 0; i < m; ++i) z[i] /= alpha[i];
    c_[m] = z[m - 1];
    for (std::size_t i = m - 1; i-- > 0;) c_[i + 1] = z[i] - gamma[i] * c_[i + 2];
  }

  /// value at x, for xmin <= x <= xmax
  double operator()(double x) const {
    if (x < x_.front() || x > x_.back()) throw std::domain_error("spline evaluated outside its data range");
    // interval i with x_[i] <= x < x_[i+1] (the last interval for x == xmax);
    // like GSL's accelerator, try the interval of the previous call first
    std::size_t i = last_;
    if (x < x_[i] || (x >= x_[i + 1] && i + 2 < x_.size())) {
      i = std::upper_bound(x_.begin(), x_.end(), x) - x_.begin();
      i = std::min(std::max<std::size_t>(i, 1), x_.size() - 1) - 1;
      last_ = i;
    }
    const double h = x_[i + 1] - x_[i], dx = x - x_[i];
    if (h == 0.) return y_[i];
    const double b = (y_[i + 1] - y_[i]) / h - h * (c_[i + 1] + 2. * c_[i]) / 3.;
    const double d = (c_[i + 1] - c_[i]) / (3. * h);
    return y_[i] + dx * (b + dx * (c_[i] + dx * d));
  }

private:
  std::vector<double> x_, y_, c_;
  mutable std::size_t last_ = 0;  // cached interval (each OpenMP thread has its own copy)
};

inline double fun(double omegaprime, double xmin_data, double xmax_data, const natural_cubic_spline &input_spline){
  if(omegaprime <= xmin_data || omegaprime >= xmax_data) return 0.;
  return -1./M_PI*input_spline(omegaprime);
}
double fun2(double omegaprime, const natural_cubic_spline &input_spline){
  return input_spline(omegaprime);
}
inline double fun3(double omegaprime, double omega, double chi2_omega, double xmin_data, double xmax_data, const natural_cubic_spline &input_spline){
  if(omegaprime==omega) return 0.;
  return (fun(omegaprime, xmin_data, xmax_data, input_spline)-chi2_omega)/(omegaprime-omega);
}


double integrate(double lower_limit, double upper_limit, double pole_location, double xmin_data, double xmax_data,const natural_cubic_spline &input_spline){
  //set the global vars
  double omega=pole_location;
  //chi2_omega = -1./M_PI*Sigma_2(omega)
  double chi2_omega=fun(omega, xmin_data, xmax_data, input_spline);
  
  //compute the offset
  double offset=(pole_location<=xmin_data || pole_location >= xmax_data)?0.:chi2_omega*(-log(std::abs((upper_limit-pole_location)/(-pole_location+lower_limit))));
  
  
  //return values
  double result;
  int N=100000;
  double dx=(upper_limit-lower_limit)/(double)N;
  
  double I=0.;
  int i;
  {
    I=0.;
    omega=pole_location;
    for(i=1;i<N-1;i++){
      double xi=i*dx+lower_limit;
      I+=(i%2==0?2.:4.)*fun3(xi, omega, chi2_omega, xmin_data, xmax_data, input_spline);
    }
  }
  I+=fun(upper_limit, xmin_data, xmax_data, input_spline)+fun(lower_limit, xmin_data, xmax_data, input_spline);
  I*=dx/3.;
  result=I;
  
  return result-offset;
}

double integrate_norm(double lower_limit, double upper_limit, const natural_cubic_spline &input_spline){
  int N=100000;
  double dx=(upper_limit-lower_limit)/(double)N;
  
  double I;
  int i;
  {
    I=0.;
    for(i=1;i<N-1;i++){
      double xi=i*dx+lower_limit;
      I+=(i%2==0?2.:4.)*fun2(xi, input_spline);
    }
  }
  I+=fun2(upper_limit, input_spline)+fun2(lower_limit, input_spline);
  I*=dx/3.;
  double result=I;
  return result;
}

//the conventions for this are the Wikipedia conventions.
//Sigma(omega) = Sigma_1(omega) + i Sigma_2(omega) and, given Sigma_2, we
//compute Sigma_1:
//Sigma_1(omega) =

enum direction_type{
  imag_to_real,
  real_to_imag
};

int main(int argc, char **argv){
  //define parameter values
  namespace po = boost::program_options;
  std::string input_file_name;
  std::string output_file_name;
  double direction_sign=-1.;
  direction_type direction;
  double xmin=-100.;
  double xmax=100.;
  
  //read in and parse command line options
  po::options_description desc("Allowed options");
  desc.add_options()
  ("help", "show this help")
  ("input_file", po::value<std::string>(&input_file_name), "Input file, e.g. Im(Sigma) selfenergy file out of continuation")
  ("output_file", po::value<std::string>(&output_file_name), "Output file, e.g. Re(Sigma) and Im(Sigma)")
  ("imag_to_real", "input is imaginary part, produce real part")
  ("real_to_imag", "input is real part, produce imaginary part")
  ;
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);
  
  if (vm.count("help")) {
    std::cout<<desc;
    return 1;
  }
  if(!vm.count("input_file")) throw std::runtime_error("you need to specify the Im(Sigma) file.");
  if(!vm.count("output_file")) throw std::runtime_error("you need to specify the Re(Sigma) file.");
  if(!vm.count("imag_to_real") && !vm.count("real_to_imag")) throw std::runtime_error("you need to specify the direction with --imag_to_real or --real_to_imag.");
  if(vm.count("real_to_imag")){
    direction_sign=1.;
    direction=real_to_imag;
  }else{
    direction_sign=-1;
    direction=imag_to_real;
  }

  //deal with input and output files
  std::ofstream output_file(output_file_name.c_str());
  std::ifstream input_file(input_file_name.c_str());
  if(!input_file.is_open()) throw std::runtime_error("input sigma file not open.");
  if(!output_file.is_open()) throw std::runtime_error("output sigma file not open.");
  output_file.precision(17);
 
 
  //read in the file for sigma imag:
  std::vector<double> xgrid;
  std::vector<double> input_data;
  std::vector<double> output_xgrid;
  std::vector<std::pair<double,double> > output_data;
  do{
    double x,y;
    std::string line;
    getline(input_file, line);
    if(line.length()!=0){
      std::stringstream sstream(line);
      sstream>>x>>y>>std::ws;
      xgrid.push_back(x);
      input_data.push_back(y);
    }
  }while(!input_file.eof());
  double xmin_data=xgrid[0];
  double xmax_data=xgrid.back();
 
  
  //interpolate the self energy
  const natural_cubic_spline input_spline(xgrid, input_data);
  
  //compute the self energy normalization
  double norm=1.;
  if(direction_sign==1.){
    norm=integrate_norm(xgrid[0]+1.e-12, xgrid.back()-1.e-12, input_spline);
    std::cout<<"integrated norm is: "<<norm<<std::endl;
  }
  
  
  //compute an output xgrid
  for(double x=xmin;x<xmax+1.e-5;x+=(x<-2.000001 || x>1.999999 ?1.e-1:1.e-2)){
    //for(double x=xmin;x<xmax+1.e-5;x+=(x<-2.000001 || x>1.999999 ?1.e-2:1.e-4)){
    output_xgrid.push_back(x);
  }
  output_data.resize(output_xgrid.size());
#pragma omp parallel default(none) firstprivate(output_xgrid, xmin_data, xmax_data, input_spline, xgrid, direction_sign) shared(output_data)
  {
#pragma omp for
    for(std::size_t i=0;i<output_xgrid.size();++i){
      double omega=output_xgrid[i];
      double kk_integral=direction_sign*integrate(xgrid[0]+1.e-12, xgrid.back()-1.e-12, omega, xmin_data, xmax_data, input_spline);
      double kk_sourceval=omega<=xgrid[0]|| omega>=xgrid.back()?0.:input_spline(omega);
      output_data[i]=std::make_pair(kk_integral, kk_sourceval);
    }
  }
  for(std::size_t i=0;i<output_xgrid.size();++i){
    if(direction==imag_to_real){
      output_file<<output_xgrid[i]<<" "<<output_data[i].first<<" "<<output_data[i].second<<std::endl;
    }else{
      output_file<<output_xgrid[i]<<" "<<output_data[i].second<<" "<<output_data[i].first<<std::endl;
    }
  }
}
