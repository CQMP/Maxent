#include<iostream>
#include<fstream>
#include<stdexcept>
#include<vector>
#include<complex>
#include<cmath>
#include<cstdlib>
#include<alps/parameter.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <boost/program_options.hpp>

inline double fun(double omegaprime, double xmin_data, double xmax_data, const gsl_spline *input_spline, gsl_interp_accel *acc){
  if(omegaprime <= xmin_data || omegaprime >= xmax_data) return 0.;
  return -1./M_PI*gsl_spline_eval (input_spline, omegaprime, acc);
}
double fun2(double omegaprime, const gsl_spline *input_spline, gsl_interp_accel *acc){
  return gsl_spline_eval (input_spline, omegaprime, acc);
}
inline double fun3(double omegaprime, double omega, double chi2_omega, double xmin_data, double xmax_data, const gsl_spline *input_spline, gsl_interp_accel *acc){
  if(omegaprime==omega) return 0.;
  return (fun(omegaprime, xmin_data, xmax_data, input_spline, acc)-chi2_omega)/(omegaprime-omega);
}


double integrate(double lower_limit, double upper_limit, double pole_location, double xmin_data, double xmax_data,const gsl_spline *input_spline, gsl_interp_accel *acc){
  //set the global vars
  double omega=pole_location;
  //chi2_omega = -1./M_PI*Sigma_2(omega)
  double chi2_omega=fun(omega, xmin_data, xmax_data, input_spline, acc);
  
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
      I+=(i%2==0?2.:4.)*fun3(xi, omega, chi2_omega, xmin_data, xmax_data, input_spline, acc);
    }
  }
  I+=fun(upper_limit, xmin_data, xmax_data, input_spline, acc)+fun(lower_limit, xmin_data, xmax_data, input_spline, acc);
  I*=dx/3.;
  result=I;
  
  //std::cout<<omega<<" "<<chi2_omega<<" "<<offset<<" "<<result<<" "<<result-offset<<std::endl;
  
  return result-offset;
}

double integrate_norm(double lower_limit, double upper_limit, const gsl_spline *input_spline, gsl_interp_accel *acc){
  int N=100000;
  double dx=(upper_limit-lower_limit)/(double)N;
  
  double I;
  int i;
  {
    I=0.;
    for(i=1;i<N-1;i++){
      double xi=i*dx+lower_limit;
      I+=(i%2==0?2.:4.)*fun2(xi, input_spline, acc);
    }
  }
  I+=fun2(upper_limit, input_spline, acc)+fun2(lower_limit, input_spline, acc);
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
  namespace po = boost::program_options;
  std::string input_file_name;
  std::string output_file_name;
  //double global_scale_factor=1.;
  double direction_sign=-1.;
  direction_type direction;
  double xmin=-100.;
  double xmax=100.;
  
  po::options_description desc("Allowed options");
  desc.add_options()
  ("help", "show this help")
  ("input_file", po::value<std::string>(&input_file_name), "Input file, e.g. Im(Sigma) selfenergy file out of continuation")
  ("output_file", po::value<std::string>(&output_file_name), "Output file, e.g. Re(Sigma) and Im(Sigma)")
  //("global_scale_factor", po::value<double>(&global_scale_factor), "scale factor for the integral (real and imag part) to account for pis from maxent")
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
  //if(!vm.count("new_norm")) throw std::runtime_error("you need to specify the new norm with --new_norm.");
  if(!vm.count("imag_to_real") && !vm.count("real_to_imag")) throw std::runtime_error("you need to specify the direction with --imag_to_real or --real_to_imag.");
  if(vm.count("real_to_imag")){
    direction_sign=1.;
    direction=real_to_imag;
  }else{
    direction_sign=-1;
    direction=imag_to_real;
  }
  std::ofstream output_file(output_file_name.c_str());
  std::ifstream input_file(input_file_name.c_str());
  if(!input_file.is_open()) throw std::runtime_error("input sigma file not open.");
  if(!output_file.is_open()) throw std::runtime_error("output sigma file not open.");
 
 
  //read in the file for sigma imag:
  std::vector<double> xgrid;
  std::vector<double> input_data;
  std::vector<double> output_xgrid;
  std::vector<std::pair<double,double> > output_data;
  do{
    double x,y;
    input_file>>x>>y>>std::ws;
    xgrid.push_back(x);
    input_data.push_back(y);
    //input_data.push_back(y*global_scale_factor);
  }while(!input_file.eof());
  int N=input_data.size();
  double xmin_data=xgrid[0];
  double xmax_data=xgrid.back();
 
  
  //interpolate the self energy
  gsl_spline *input_spline= gsl_spline_alloc (gsl_interp_cspline, N);
  gsl_interp_accel *acc= gsl_interp_accel_alloc ();
  gsl_spline_init (input_spline, &(xgrid[0]), &(input_data[0]), N);
  
  //compute the self energy normalization
  double norm=1.;
  if(direction_sign==1.){
    norm=integrate_norm(xgrid[0]+1.e-12, xgrid.back()-1.e-12, input_spline, acc);
    std::cout<<"integrated norm is: "<<norm<<std::endl;
  }
  gsl_interp_accel_free(acc);
  
  
  //compute an output xgrid
  for(double x=xmin;x<xmax+1.e-5;x+=(x<-2.000001 || x>1.999999 ?1.e-1:1.e-2)){
    //for(double x=xmin;x<xmax+1.e-5;x+=(x<-2.000001 || x>1.999999 ?1.e-2:1.e-4)){
    output_xgrid.push_back(x);
  }
  output_data.resize(output_xgrid.size());
#pragma omp parallel default(none) firstprivate(output_xgrid, xmin_data, xmax_data, input_spline, xgrid, direction_sign) shared(output_data)
  {
    gsl_interp_accel *acc= gsl_interp_accel_alloc ();
#pragma omp for
    for(std::size_t i=0;i<output_xgrid.size();++i){
      double omega=output_xgrid[i];
      double kk_integral=direction_sign*integrate(xgrid[0]+1.e-12, xgrid.back()-1.e-12, omega, xmin_data, xmax_data, input_spline, acc);
      double kk_sourceval=omega<=xgrid[0]|| omega>=xgrid.back()?0.:gsl_spline_eval (input_spline, omega, acc);
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
