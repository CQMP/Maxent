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
#include "../src/SpM_method.hpp"
#include "gtest.h"
#include <iostream>
#include <cmath>
#include <gsl/gsl_integration.h>
#include <boost/random/normal_distribution.hpp>
#include <boost/random.hpp>


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
void imag_time_gf_to_param(alps::params &p, int nmax, double beta){
    
    //todo: put in the actual imaginary time green's function
    for(int j=0;j<nmax;j++){
        
        std::stringstream pval; pval<<"X_"<<j;
        std::stringstream pval2; pval2<<"TAU_"<<j;
        p[pval.str()] = 0;
        p[pval2.str()] = beta/(nmax-1)*j;
    }
}

TEST(SpM,KernelCanDealWithGaussian){
    alps::params p;
    SVDContinuation::define_parameters(p);

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

  SVDContinuation C(p);

  Eigen::VectorXd omega_vals=C.omega_coord();
  Eigen::VectorXd delta_omega=C.delta_omega();
  Eigen::VectorXd input_grid=C.inputGrid();

  //fabricate the input data
  Eigen::VectorXd real_comparison_data=delta_omega;
  for(int i=0;i<real_comparison_data.size();++i){
    real_comparison_data[i]*=gaussian(omega_vals[i]);
  }

  //let the kernel multiply to imag frequency
  Eigen::VectorXd backcont=C.U()*C.Sigma()*C.Vt()*real_comparison_data;

  for(int i=0;i<ndat;++i){
    std::stringstream pval; pval<<"X_"<<i;
    EXPECT_NEAR(backcont[i], p[pval.str()], 1.e-7);
  }

}
TEST(SpM,KernelCanDealWithTripleGaussian){
    alps::params p;
    SVDContinuation::define_parameters(p);

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

  SVDContinuation C(p);

  Eigen::VectorXd omega_vals=C.omega_coord();
  Eigen::VectorXd delta_omega=C.delta_omega();
  Eigen::VectorXd input_grid=C.inputGrid();

  //fabricate the input data
  Eigen::VectorXd real_comparison_data=delta_omega;
  for(int i=0;i<real_comparison_data.size();++i){
    real_comparison_data[i]*=triple_gaussian(omega_vals[i]);
  }

  //let the kernel multiply to imag frequency
  Eigen::VectorXd backcont=C.U()*C.Sigma()*C.Vt()*real_comparison_data;

  for(int i=0;i<ndat;++i){
    std::stringstream pval; pval<<"X_"<<i;
    EXPECT_NEAR(backcont[i], p[pval.str()], 1.e-5);
  }

}
TEST(SpM,ADMMTestChiSquareNormForExactSolution){
    alps::params p;
    SVDContinuation::define_parameters(p);

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

  SVDContinuation C(p);

  double rho=1, rhoprime=1., lambda=1.;
  ADMM A(C.Vt(), C.U().transpose()*C.y(), C.delta_omega(), C.Sigma().diagonal(), rho, rhoprime, lambda);

  //find the correct spectral function
  Eigen::VectorXd real_comparison_data=C.delta_omega();
  for(int i=0;i<real_comparison_data.size();++i){
    real_comparison_data[i]*=triple_gaussian(C.omega_coord(i));
  }
  
  //plug in the correct spectral function
  A.externally_set_xprime(C.Vt()*real_comparison_data);
  std::cout<<"chi square is: "<< A.chisquare_term()<<std::endl;
  EXPECT_NEAR(0, A.chisquare_term(), 1.e-8);
  std::cout<<"integral is: "<< ((A.spectral_function()).cwiseProduct(C.delta_omega())).sum()<<std::endl;
  EXPECT_NEAR(0, A.constraint_violation_norm(), 1.e-4);
}
TEST(SpM,ADMMTestPositivityForExactSolution){
    alps::params p;
    SVDContinuation::define_parameters(p);

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

  SVDContinuation C(p);

  double rho=1, rhoprime=1., lambda=1.;
  ADMM A(C.Vt(), C.U().transpose()*C.y(), C.delta_omega(), C.Sigma().diagonal(), rho, rhoprime, lambda);

  //find the correct spectral function
  Eigen::VectorXd real_comparison_data=C.delta_omega();
  for(int i=0;i<real_comparison_data.size();++i){
    real_comparison_data[i]*=triple_gaussian(C.omega_coord(i));
  }

  std::cout<<"positivity constraint is: "<< A.constraint_violation_positivity()<<std::endl;
  EXPECT_NEAR(0, A.constraint_violation_positivity(), 1.e-7);
}
TEST(SpM,ADMMSpectralFunctionAccuracy){
    alps::params p;
    SVDContinuation::define_parameters(p);

    double beta=200;
    int ndat=5000;
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

  SVDContinuation C(p);

  double rho=1, rhoprime=1., lambda=1.;
  ADMM A(C.Vt(), C.U().transpose()*C.y(), C.delta_omega(), C.Sigma().diagonal(), rho, rhoprime, lambda);

  //find the correct spectral function
  Eigen::VectorXd real_comparison_data=C.delta_omega();
  for(int i=0;i<real_comparison_data.size();++i){
    real_comparison_data[i]*=triple_gaussian(C.omega_coord(i));
  }

  std::cout<<"difference between spectral function reconstructed from singular space and direct calc:"<<std::endl;
  Eigen::VectorXd diff=(real_comparison_data-A.spectral_function());
  std::cout<<"diff norm: "<<diff.norm()<<std::endl;
  std::cout<<" integral: "<<diff.cwiseProduct(C.delta_omega()).sum()<<std::endl;
  std::cout<<" integral abs: "<<(diff.cwiseProduct(C.delta_omega())).lpNorm<1>()<<std::endl;
  
  //plug in the correct spectral function
  A.externally_set_xprime(C.Vt()*real_comparison_data);
  Eigen::VectorXd spectral_function=A.spectral_function();
  for(int i=0;i<spectral_function.size();++i){
    std::cout<<C.omega_coord(i)<<" "<<triple_gaussian(C.omega_coord(i))<<" "<<spectral_function[i]<<std::endl;
    EXPECT_NEAR(triple_gaussian(C.omega_coord(i)), spectral_function[i], 1.e-6);
  }
}
TEST(SpM,DISABLED_ADMML1OptimizationReducesL1Norm){
    alps::params p;
    SVDContinuation::define_parameters(p);

    double beta=200;
    int ndat=5000;
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

  SVDContinuation C(p);

  //find the correct spectral function
  Eigen::VectorXd real_comparison_data=C.delta_omega();
  for(int i=0;i<real_comparison_data.size();++i){
    real_comparison_data[i]*=triple_gaussian(C.omega_coord(i));
  }

  double rho=1, rhoprime=1., lambda=1.;
  ADMM A(C.Vt(), C.U().transpose()*C.y(), C.delta_omega(), C.Sigma().diagonal(), rho, rhoprime, lambda);
  //plug in the correct spectral function
  A.externally_set_xprime(C.Vt()*real_comparison_data);
  A.update_zprime();
  Eigen::VectorXd xprime=A.xprime();
  Eigen::VectorXd zprime=A.zprime();
  for(int i=0;i<zprime.size();++i){
    //need to continue here and actually check the L1 norm  reduction
    std::cout<<i<<" "<<xprime[i]<<" "<<zprime[i]<<std::endl;
  }

}

TEST(SpM,KernelInImagTimeCanDealWithTripleGaussian){
    alps::params p;
    SVDContinuation::define_parameters(p);
    
    double beta=20;
    int ndat=500;
    p["BETA"]=beta;
    p["NDAT"]=ndat;
    p["NO_ERRORS"]=true;
    p["DATASPACE"]="time";
    p["KERNEL"]="fermionic";
    p["PARTICLE_HOLE_SYMMETRY"]=true;
    p["VERBOSE"]=false;
    
    //grid
    p["OMEGA_MAX"]=10;
    p["OMEGA_MIN"]=-10;
    p["CUT"]=0.1;
    p["FREQUENCY_GRID"]="Lorentzian";
    p["NFREQ"]=1000;
    
    imag_time_gf_to_param(p, ndat, beta);
    
    SVDContinuation C(p);
    
    Eigen::VectorXd omega_vals=C.omega_coord();
    Eigen::VectorXd delta_omega=C.delta_omega();
    Eigen::VectorXd input_grid=C.inputGrid();
    
    //fabricate the input data
    Eigen::VectorXd real_comparison_data=delta_omega;
    for(int i=0;i<real_comparison_data.size();++i){
        real_comparison_data[i]*=triple_gaussian(omega_vals[i]);
    }
    
    //let the kernel multiply to imag frequency
    Eigen::VectorXd backcont=C.U()*C.Sigma()*C.Vt()*real_comparison_data;
    
    for(int i=0;i<ndat;++i){
        std::cout<<i*beta/(ndat-1)<<" "<<backcont[i]<<std::endl;
    }
    
}

TEST(SpM,TestRandom){
  boost::random::mt19937 rng;
  boost::random::normal_distribution<double> norm(0., 0.00001);
  for(int i=0;i<1001;++i){
      std::cout<<i<<" "<<norm(rng)<<std::endl;;
  }
  
    std::string filename1="noise.dat";
    std::ofstream noise_file(filename1);
    for(int i=0;i<1001;++i){
        noise_file<<i<<" "<<norm(rng)<<std::endl;
    }

}
