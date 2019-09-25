/*
 * Copyright (C) 1998-2018 ALPS Collaboration.
 * All rights reserved. Use is subject to license terms. See LICENSE.TXT
 * For use in publications, see ACKNOWLEDGE.TXT
 */
#include "maxent_kernel.hpp"
#include <cmath>
#include <boost/algorithm/string.hpp>    
#include <boost/lexical_cast.hpp>
#include <boost/math/special_functions/legendre.hpp>    //for Legendre kernel
#include <boost/math/special_functions/factorials.hpp>  //for Legendre kernel
#include <gsl/gsl_integration.h>                        //for Legendre Kernel

namespace bmth = boost::math;

kernel::kernel(alps::params &p, const vector_type& freq, vector_type &inputGrid):
ndat_(p["NDAT"]),
nfreq_(p["NFREQ"]),
T_(1./static_cast<double>(p["BETA"])),
K_(ndat_,nfreq_)
{
  using namespace boost::numeric;
  //K_.clear();
  K_=matrix_type::Zero(ndat_,nfreq_);
  std::string dataspace_name = p["DATASPACE"];
  std::string kernel_name = p["KERNEL"];
  boost::to_lower(dataspace_name);
  boost::to_lower(kernel_name);
  bool ph_symmetry=p["PARTICLE_HOLE_SYMMETRY"];
  std::cout<<"using kernel "<<kernel_name<<" in domain "<<dataspace_name;
  if(ph_symmetry) std::cout<<" with ph symmetry"; else std::cout<<" without ph symmetry"; std::cout<<std::endl;

  set_kernel_type(dataspace_name,kernel_name, ph_symmetry);

  //determine tau points; allow passing from various sources
	if(dataspace_name=="time"){
		//Test for tau
    tau_points_.resize(ndat_);
    //programatically added tau points (see kernelTest)
    if(p.defined("TAU_1")){
        std::cout<<"Using param direct input tau points"<<std::endl;
        for(int i=0;i<ndat_;i++){
          tau_points_[i]=p["TAU_"+boost::lexical_cast<std::string>(i)];
        }
    }
    //legacy tau points in param file
    else if(p.exists("TAU_0")){
      std::cout<<"Using param input tau points"<<std::endl;
      tau_points_[0]=p["TAU_0"];
      for(int i=1;i<ndat_;i++)
        p.define<double>("TAU_"+boost::lexical_cast<std::string>(i),"");
      for(int i=1;i<ndat_;i++){
        tau_points_[i]=p["TAU_"+boost::lexical_cast<std::string>(i)];
      }
    }
    else{
      if(p.exists("X_2")){
        //using data file entry so there is no input grid
        throw std::runtime_error("Missing input tau points! Define them with TAU_0=....");
      }
      std::cout<<"Using data file tau points"<<std::endl;
      tau_points_= inputGrid;  
    }
    inputGrid = tau_points_;
  }
    
  if(ktype_==legendre_fermionic_kernel || ktype_==legendre_bosonic_kernel){
      setup_legendre_kernel(p,freq,ndat_);
  }else if(ktype_==time_fermionic_kernel){
      for (int i=0; i<ndat_; ++i) {
        double tau=tau_points_[i];
        for (int j=0; j<nfreq_; ++j) {
          double omega = freq[j]; //Default().omega_of_t(double(j)/(nfreq_-1));
          K_(i,j) =  -1. / (std::exp(omega*tau) + std::exp(-omega*(1./T_-tau)));
        }
      }
    }
    else if (ktype_==time_bosonic_kernel) {
      for (int i=0; i<ndat_; ++i) {
        double tau=tau_points_[i];
        K_(i,0) = T_;
        for (int j=1; j<nfreq_; ++j) {
          double omega = freq[j];
          K_(i,j) = 0.5*omega * (std::exp(-omega*tau) + std::exp(-omega*(1./T_-tau))) / (1 - std::exp(-omega/T_));
        }
      }
    }
    else if(ktype_== time_fermionic_legendre_kernel || ktype_==time_bosonic_legendre_kernel)
        setup_legendre_kernel(p,freq,ndat_);
    //for zero temperature, only positive frequency matters
    else if (ktype_ == time_boris_kernel) {
      for (int i=0; i<ndat_; ++i) {
        double tau=tau_points_[i];
        for (int j=0; j<nfreq_; ++j) {
          double omega = freq[j];
          K_(i,j) = -std::exp(-omega*tau);
        }
      }
    }
    else if(ktype_==time_fermionic_kernel){
        for (int i=0; i<ndat_; ++i) {
						double tau=tau_points_[i];
            for (int j=0; j<nfreq_; ++j) {
                double omega = freq[j];
                K_(i,j) =  -1.;
            }
        }
    }
    else if(ktype_==frequency_fermionic_ph_kernel) {
    for (int i=0; i<ndat_; ++i) {
      double omegan = (2*i+1)*M_PI*T_;
      inputGrid(i) = omegan;
      for (int j=0; j<nfreq_; ++j) {
        double omega = freq[j];
        K_(i,j) =  -omegan / (omegan*omegan + omega*omega);
      }
    }
  }
  else if (ktype_==frequency_bosonic_ph_kernel) {
    for (int i=0; i<ndat_; ++i) {
      double Omegan = (2*i)*M_PI*T_;
      inputGrid(i) = Omegan;
      for (int j=0; j<nfreq_; ++j) {
        double Omega = freq[j];
        if(Omega ==0) throw std::runtime_error("Bosonic kernel is singular at frequency zero. Please use grid w/o evaluation at zero.");
        K_(i,j) =  Omega*Omega / (Omegan*Omegan + Omega*Omega);
      }
    }
  }else if (ktype_==frequency_anomalous_ph_kernel) {
    for(int i=0;i<ndat_;++i){
      double omegan = (2*i+1)*M_PI*T_;
      inputGrid(i) = omegan;
      for (int j=0; j<nfreq_; ++j) {
        double omega = freq[j];
        K_(i,j) =  omega*omega / (omegan*omegan + omega*omega);
      }
    }
  }
  else{
    complex_matrix_type Kc(ndat_/2, nfreq_);
    if (ktype_==frequency_fermionic_kernel) {
      ///ndat/2 is defined as such below because ndat=number of points inputed
      ///if ph symmetry, then ndat=number of imag points
      ///otherwise ndat=total number of real+imag points = ndat/2 data points
      for (int i=0; i<ndat_/2; ++i) {
        std::complex<double> iomegan(0, (2*i+1)*M_PI*T_);
        inputGrid(i) = iomegan.imag();
        for (int j=0; j<nfreq_; ++j) {
          double omega = freq[j];
          Kc(i,j) =  1. / (iomegan - omega);
        }
      }
    }
    else if (ktype_==frequency_bosonic_kernel){
      for (int i=0; i<ndat_/2; ++i) {
        std::complex<double> iomegan(0, 2*i*M_PI*T_);
        inputGrid(i) = iomegan.imag();
        inputGrid(i+1) = iomegan.imag();
        for (int j=0; j<nfreq_; ++j) {
          double omega = freq[j];
          Kc(i,j) =  omega / (iomegan + omega);
        }
      }
    }
    else if (ktype_==frequency_anomalous_kernel){
      for (int i=0; i<ndat_/2; ++i) {
        std::complex<double> iomegan(0, (2*i+1)*M_PI*T_);
        inputGrid(i) = iomegan.imag();
        for (int j=0; j<nfreq_; ++j) {
          double omega = freq[j];
          Kc(i,j) =  -omega / (iomegan - omega);
        }
      }
    }else
      throw std::logic_error("unknown kernel type");

    for (int i=0; i<ndat_; i+=2) {
      //TODO: understand the j=1 here
      for (int j=0; j<nfreq_; ++j) {
        K_(i,j) = Kc(i/2,j).real();
        K_(i+1,j) = Kc(i/2,j).imag();
      }
    }
  }
}

void kernel::set_kernel_type(const std::string &dataspace_name, const std::string &kernel_name,
                             bool ph_symmetry){
  if(dataspace_name=="time"){
    dtype_=time_dataspace;
     
  }else if(dataspace_name=="frequency"){
    dtype_=frequency_dataspace;
  }else if(dataspace_name=="legendre"){
    dtype_=legendre_dataspace;
  }
  else
    throw std::invalid_argument("unknown dataspace name. it should be time, frequency, or legendre");

  if(dtype_==time_dataspace){
    if(kernel_name=="fermionic")
            if(dtype_==legendre_dataspace)
                ktype_=time_fermionic_legendre_kernel;
            else
                ktype_=time_fermionic_kernel;
    else if(kernel_name=="bosonic")
            if(dtype_==legendre_dataspace)
                ktype_=time_bosonic_legendre_kernel;
            else
                ktype_=time_bosonic_kernel;
    else if(kernel_name=="tzero")
      ktype_=time_boris_kernel;
    else throw std::invalid_argument("unknown kernel name. In the time domain it should be fermionic, bosonic, or tzero.");
  }else if(dtype_ == legendre_dataspace){
      if(kernel_name=="fermionic")
          ktype_=legendre_fermionic_kernel;
      else if(kernel_name=="bosonic"){
          ktype_=legendre_bosonic_kernel;
      }
      else throw std::invalid_argument("unknown kernel name. In the legendre domain it should be fermionic or bosonic");
  }
  else{
    if(ph_symmetry){
      if(kernel_name== "fermionic")
        ktype_=frequency_fermionic_ph_kernel;
      else if (kernel_name=="bosonic")
        ktype_=frequency_bosonic_ph_kernel;
      else if (kernel_name=="anomalous")
        ktype_=frequency_anomalous_ph_kernel;
      else throw std::invalid_argument("unknown kernel name. In the particle hole symmetric frequency domain it should be fermionic, bosonic, or anomalous.");
    }else{
      if(kernel_name== "fermionic")
        ktype_=frequency_fermionic_kernel;
      else if (kernel_name=="bosonic")
        ktype_=frequency_bosonic_kernel;
      else if (kernel_name=="anomalous")
        ktype_=frequency_anomalous_kernel;
      else throw std::invalid_argument("unknown kernel name. In the non-particle hole symmetric frequency domain it should be fermionic, bosonic, or anomalous.");
    }
  }

}
struct integrand_params {int l; double omega;double T_;};

///Integrand of the Legendre Kernel for GSL integration
double  legendre_kernel_integrand(double x, void * params){
    double tau = x;
    //parms = [l,omega,T_]
    integrand_params *p = (integrand_params *)params;
    int l = p->l;
    double omega = p->omega;
    double T_ = p->T_;
    //std::cout<< l<< std::endl;
    return bmth::legendre_p(l, 2*tau*T_-1)*std::exp(-tau*omega)/(1+std::exp(-omega/T_));
}
void kernel::setup_legendre_kernel(const alps::params &p, const vector_type& freq,const int lmax){
    if(lmax>boost::math::max_factorial<double>::value)
        throw std::runtime_error("lmax is greater than boost factorial precision");
    
    const double PI = std::acos(-1);
    const std::complex<double> CONE(0,1);
    //recall that ndat()=lmax

    int N = 20000/2;
    gsl_integration_workspace *w = gsl_integration_workspace_alloc (1000);
    
    for(int l=0;l<lmax;l++){
        for(int j=0;j<nfreq_;j++){
            double I=1;
            double I1=0;
            double omega =freq[j];
            double h = (1/T_-0)/(2*N);
            //int Pl(x(tau))*exp(-tau*omega)/(1\pm exp(-beta*omega))
            
            //Simpsons with
            //Simpson's method of integrations
            //eval endpoints
            /*I1 += bmth::legendre_p(l, 1.0)*std::exp(-0*omega)/(1+sign*std::exp(-omega/T_));
            I1 += bmth::legendre_p(l, -1.0)*std::exp(-omega/T_)/(1+sign*std::exp(-omega/T_));
            for(int i=1;i<N;i++){
                double tau = 0 + 2*i*h;
                I1+=2*bmth::legendre_p(l, 2*tau*T_-1)*std::exp(-tau*omega)/(1+std::exp(-omega/T_));
            }
            for(int i=1;i<N+1;i++){
                double tau= 0+ (2*i-1)*h;
                I1+=4*bmth::legendre_p(l, 2*tau*T_-1)*std::exp(-tau*omega)/(1+std::exp(-omega/T_));
            }
            I1*=h/3;//*/
            
            double a = 0;
            double b = 1/T_;
            double epsabs=1.49e-08; //python default resolution
            double epsrel=1.49e-08;
            double result,err;
            size_t nval,limit=500;
            
            gsl_function F;
            F.function = &legendre_kernel_integrand;

            integrand_params p = {l,omega,T_};
            F.params = &p;

            gsl_integration_qag(&F, a,b,epsabs,epsrel, limit, GSL_INTEG_GAUSS61, w, &result, &err);
            I1=result; 
            
            K_(l,j) = -sqrt(2*l+1)*I1;
        }
    }
    gsl_integration_workspace_free (w);
}

