/*
 * Copyright (C) 1998-2016 ALPS Collaboration
 * 
 *     This program is free software; you can redistribute it and/or modify it
 *     under the terms of the GNU General Public License as published by the Free
 *     Software Foundation; either version 2 of the License, or (at your option)
 *     any later version.
 * 
 *     This program is distributed in the hope that it will be useful, but WITHOUT
 *     ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *     FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for 
 *     more details.
 * 
 *     You should have received a copy of the GNU General Public License along
 *     with this program; if not, write to the Free Software Foundation, Inc., 59
 *     Temple Place, Suite 330, Boston, MA 02111-1307 USA
 * For use in publications, see ACKNOWLEDGE.TXT
 */
#include "../src/NNLS_method.hpp"
#include "gtest.h"
#include <iostream>
#include <cmath>
#include <gsl/gsl_integration.h>

double gaussian(double x){
  //return gaussian(x,0,1,1);
  return 1./(std::sqrt(2*M_PI))*std::exp(-x*x/2.);
}
double general_gaussian(double x, double mu, double width, double weight){
  double sigma=width/std::sqrt(2*log(2));
  return weight/(sigma*std::sqrt(2*M_PI))*std::exp(-0.5*((x-mu)*(x-mu)/(sigma*sigma)));
}
double triple_gaussian(double x){
  return general_gaussian(x, 0, 0.15, 0.2)+general_gaussian(x,1,0.8,0.4)+general_gaussian(x,-1,0.8,0.4);
}
  //return gaussian(x,0,1,1);
double imag_backcont_gaussian(double x, void *params){
  double wn=*((double*)(params));
  return (gaussian(x)/(std::complex<double>(0,wn)-x)).imag();
}
double imag_backcont_triple_gaussian(double x, void *params){
  double wn=*((double*)(params));
  return (triple_gaussian(x)/(std::complex<double>(0,wn)-x)).imag();
}
void matsubara_gf_to_param(alps::params &p, double beta, int nmax, double (*fun)(double, void*)){

  double omega_min=-15.;
  double omega_max=15.;
  gsl_integration_workspace *w = gsl_integration_workspace_alloc (1000);

  gsl_function F;
  F.function = fun;
  for(int j=0;j<nmax;j++){
     double omegan=(2.*j+1)*M_PI/beta;

     double epsabs=1.49e-12; //python default resolution
     double epsrel=1.49e-12;
     double result,err;
     size_t limit=1000;

     F.params = &omegan;

     gsl_integration_qag(&F, omega_min,omega_max,epsabs,epsrel, limit, GSL_INTEG_GAUSS61, w, &result, &err);
     std::stringstream pval; pval<<"X_"<<j;
     p[pval.str()] = result;
  }
  gsl_integration_workspace_free (w);
}
TEST(NNLS,KernelCanDealWithGaussian){
    alps::params p;
    NNLS_Simulation::define_parameters(p);

    double beta=100;
    int ndat=500;
    p["BETA"]=beta;
    p["NDAT"]=ndat;
    p["NO_ERRORS"]=true;
    p["DATASPACE"]="frequency";
    p["KERNEL"]="fermionic";
    p["PARTICLE_HOLE_SYMMETRY"]=true;
    p["TEXT_OUTPUT"]=false;
    p["VERBOSE"]=false;

    //grid
    p["OMEGA_MAX"]=10;
    p["OMEGA_MIN"]=-10;
    p["CUT"]=0.1;
    p["FREQUENCY_GRID"]="Lorentzian";
    p["NFREQ"]=10000;

  matsubara_gf_to_param(p, beta, ndat, &imag_backcont_gaussian);

  NNLS_Simulation C(p);

  Eigen::VectorXd omega_vals=C.omega_coord();
  Eigen::VectorXd delta_omega=C.delta_omega();
  Eigen::VectorXd input_grid=C.inputGrid();

  //fabricate the input data
  Eigen::VectorXd real_comparison_data=delta_omega;
  for(int i=0;i<real_comparison_data.size();++i){
    real_comparison_data[i]*=gaussian(omega_vals[i]);
  }

  //let the kernel multiply to imag frequency
  Eigen::VectorXd backcont=C.K()*real_comparison_data;

  for(int i=0;i<ndat;++i){
    std::stringstream pval; pval<<"X_"<<i;
    EXPECT_NEAR(backcont[i], p[pval.str()], 1.e-7);
  }

}
TEST(NNLS,KernelCanDealWithTripleGaussian){
    alps::params p;
    NNLS_Simulation::define_parameters(p);

    double beta=20;
    int ndat=500;
    p["BETA"]=beta;
    p["NDAT"]=ndat;
    p["NO_ERRORS"]=true;
    p["DATASPACE"]="frequency";
    p["KERNEL"]="fermionic";
    p["PARTICLE_HOLE_SYMMETRY"]=true;
    p["VERBOSE"]=false;

    //grid
    p["OMEGA_MAX"]=10;
    p["OMEGA_MIN"]=-10;
    p["CUT"]=0.1;
    p["FREQUENCY_GRID"]="Lorentzian";
    p["NFREQ"]=1000;

  matsubara_gf_to_param(p, beta, ndat, &imag_backcont_triple_gaussian);

  NNLS_Simulation C(p);

  Eigen::VectorXd omega_vals=C.omega_coord();
  Eigen::VectorXd delta_omega=C.delta_omega();
  Eigen::VectorXd input_grid=C.inputGrid();

  //fabricate the input data
  Eigen::VectorXd real_comparison_data=delta_omega;
  for(int i=0;i<real_comparison_data.size();++i){
    real_comparison_data[i]*=triple_gaussian(omega_vals[i]);
  }

  //let the kernel multiply to imag frequency
  Eigen::VectorXd backcont=C.K()*real_comparison_data;

  for(int i=0;i<ndat;++i){
    std::stringstream pval; pval<<"X_"<<i;
    EXPECT_NEAR(backcont[i], p[pval.str()], 1.e-5);
  }
}
TEST(NNLS,FirstDerivativeWorksAsAdvertised){
    alps::params p;
    NNLS_Simulation::define_parameters(p);

    double beta=20;
    int ndat=500;
    p["BETA"]=beta;
    p["NDAT"]=ndat;
    p["NO_ERRORS"]=true;
    p["DATASPACE"]="frequency";
    p["KERNEL"]="fermionic";
    p["PARTICLE_HOLE_SYMMETRY"]=true;
    p["VERBOSE"]=false;

    //grid
    p["OMEGA_MAX"]=10;
    p["OMEGA_MIN"]=-10;
    p["CUT"]=0.1;
    p["FREQUENCY_GRID"]="Lorentzian";
    p["NFREQ"]=1000;

  matsubara_gf_to_param(p, beta, ndat, &imag_backcont_triple_gaussian);

  NNLS_Simulation C(p);

  Eigen::VectorXd omega_vals=C.omega_coord();
  Eigen::VectorXd delta_omega=C.delta_omega();
  Eigen::VectorXd input_grid=C.inputGrid();

  //fabricate the input data
  Eigen::VectorXd real_comparison_data(omega_vals.size());
  for(int i=0;i<real_comparison_data.size();++i){
    real_comparison_data[i]=triple_gaussian(omega_vals[i]);
  }

  C.compute_first_derivative_matrix();
  Eigen::MatrixXd L1=C.L1();
  Eigen::VectorXd Deriv=C.L1()*real_comparison_data;

  std::cout<<"Deriv size: "<<Deriv.size()<<" omega vals size: "<<omega_vals.size()<<std::endl;
  for(int i=0;i<Deriv.size();++i){
    std::cout<<omega_vals[i+2]<<" "<<Deriv[i]<<std::endl;
  }
}
TEST(NNLS,SecondDerivativeWorksAsAdvertised){
    alps::params p;
    NNLS_Simulation::define_parameters(p);

    double beta=20;
    int ndat=500;
    p["BETA"]=beta;
    p["NDAT"]=ndat;
    p["NO_ERRORS"]=true;
    p["DATASPACE"]="frequency";
    p["KERNEL"]="fermionic";
    p["PARTICLE_HOLE_SYMMETRY"]=true;
    p["VERBOSE"]=false;

    //grid
    p["OMEGA_MAX"]=10;
    p["OMEGA_MIN"]=-10;
    p["CUT"]=0.1;
    p["FREQUENCY_GRID"]="Lorentzian";
    p["NFREQ"]=1000;

  matsubara_gf_to_param(p, beta, ndat, &imag_backcont_triple_gaussian);

  NNLS_Simulation C(p);

  Eigen::VectorXd omega_vals=C.omega_coord();
  Eigen::VectorXd delta_omega=C.delta_omega();
  Eigen::VectorXd input_grid=C.inputGrid();

  //fabricate the input data
  Eigen::VectorXd real_comparison_data(omega_vals.size());
  for(int i=0;i<real_comparison_data.size();++i){
    real_comparison_data[i]=triple_gaussian(omega_vals[i]);
  }

  C.compute_second_derivative_matrix();
  Eigen::MatrixXd L2=C.L2();
  Eigen::VectorXd Deriv2=C.L2()*real_comparison_data;

  for(int i=0;i<Deriv2.size();++i){
    std::cout<<omega_vals[i+2]<<" "<<Deriv2[i]<<std::endl;
  }
}

