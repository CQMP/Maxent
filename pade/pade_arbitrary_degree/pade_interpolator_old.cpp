/*
 * Copyright (C) 1998-2015 ALPS Collaboration. See COPYRIGHT.TXT
 * All rights reserved. Use is subject to license terms. See LICENSE.TXT
 * For use in publications, see ACKNOWLEDGE.TXT
 */

#include "pade.hpp"
#include <fstream>


pade_interpolator::pade_interpolator(const alps::params &p){
  pade_mu_=p["PADE_NUMERATOR_DEGREE"];
  pade_nu_=p["PADE_DENOMINATOR_DEGREE"];
  if(pade_mu_+pade_nu_ +1 != p["NDAT"]) throw std::runtime_error("Ill defined interpolation: please make sure that numerator degree + denominator degree + 1 = # of data points");
  mpf_set_default_prec(p["FLOAT_PRECISION"]|256);
  find_epsilon();
}

//compute 'alpha' = x - xs (distance of interpolation to grid point). If we find an interpolation point identical to the grid point, there is no need to interpolate and we return true.
std::pair<bool, int> pade_interpolator::init_alpha(int N, pade_complex_vector_type &alpha_sx, const pade_complex_vector_type &s, const pade_complex_type &x)const{
  for(int i=0;i<N;++i){
    alpha_sx[i]=x-s[i];
    if(alpha_sx[i]==pade_complex_type(0.)){
      return std::make_pair(true, i);
    }
  }
  return std::make_pair(false, 0.);
}
//a single step in the Neville scheme, see above Eq. 2.2.3.7
void pade_interpolator::neville_step(pade_complex_vector_type &phi_mu, const pade_complex_vector_type &phi_mu_minus_one, const pade_complex_vector_type &alpha_sx, int mu) const{
  int N=phi_mu.size();
  std::cout<<"neville step. mu is: "<<mu<<" iterating up to: "<<N-mu<<std::endl;
  for(int n=0;n<N-mu;n++){
    phi_mu[n]=(alpha_sx[n]*phi_mu_minus_one[n+1]-alpha_sx[n+mu]*phi_mu_minus_one[n])/(alpha_sx[n]-alpha_sx[n+mu]);
    std::cout<<phi_mu[n]<<" "<<(alpha_sx[n]-alpha_sx[n+mu])<<" ";
  }
  std::cout<<std::endl;
}
//a single step in the inverse Neville scheme, see Eq. 2.2.3.7
void pade_interpolator::inverse_neville_step(pade_complex_vector_type &phi_nu, /*const*/ pade_complex_vector_type &phi_nu_minus_one, const pade_complex_vector_type &alpha_sx, int nu) const{
  int N=phi_nu.size();
  //std::cout<<"nu is: "<<nu<<std::endl;
  for(int n=0;n<N-nu;n++){
    if(phi_nu_minus_one[n+1]==pade_complex_type(0.)){phi_nu_minus_one[n+1]=pade_complex_type(10*epsilon_);}
    if(phi_nu_minus_one[n  ]==pade_complex_type(0.)){phi_nu_minus_one[n  ]=pade_complex_type(10*epsilon_);}
    phi_nu[n]=(alpha_sx[n]-alpha_sx[n+nu])/(alpha_sx[n]/phi_nu_minus_one[n+1]-alpha_sx[n+nu]/phi_nu_minus_one[n]);
    std::cout<<"numerator: "<<(alpha_sx[n]-alpha_sx[n+nu])<<" denominator terms: "<<alpha_sx[n]/phi_nu_minus_one[n+1]<<" "<<alpha_sx[n+nu]/phi_nu_minus_one[n]<<std::endl;
    //std::cout<<phi_nu[n]<<" ";
  }
  //std::cout<<std::endl;
}

//this function implements the Neville scheme according to Stoer and Bulirsch, chapter 2.2.3 etc
pade_complex_type pade_interpolator::pade_scheme(pade_complex_vector_type &phi_mu_nu, const pade_complex_vector_type &s, const pade_complex_vector_type &fs, const pade_complex_vector_type &alpha_sx, const pade_complex_type &x) const{
  phi_mu_nu=fs;
  pade_complex_vector_type phi_mu_nu_minus_one(phi_mu_nu);
  int nu=0;
  int mu=0;
  
  if(pade_nu_+pade_mu_+1 != phi_mu_nu.size()) throw std::runtime_error("problem in interpolation: degree of numerator and denominator+1 has to match degree of data.");
  
  //start by doing a nu-mu inverse neville steps, until
  for(int interp_step=1;interp_step<=pade_nu_+pade_mu_;++interp_step){
    phi_mu_nu.swap(phi_mu_nu_minus_one);
    if(interp_step < pade_nu_-pade_mu_+1){
      //std::cout<<" beginning inverse step"<<std::endl;
      inverse_neville_step(phi_mu_nu, phi_mu_nu_minus_one, alpha_sx, interp_step);
      nu++;
    }else if(interp_step < pade_mu_-pade_nu_+1){
      std::cout<<" beginning direct step"<<std::endl;
      neville_step(phi_mu_nu, phi_mu_nu_minus_one, alpha_sx, interp_step);
      mu++;
    }else if (interp_step%2==0){
      std::cout<<" inverse step"<<std::endl;
      inverse_neville_step(phi_mu_nu, phi_mu_nu_minus_one, alpha_sx, interp_step);
      nu++;
    }else{
      std::cout<<" direct step"<<std::endl;
      neville_step(phi_mu_nu, phi_mu_nu_minus_one, alpha_sx, interp_step);
      mu++;
    }
  }
  return phi_mu_nu[0];
}

void pade_interpolator::pade_interpolate(const imaginary_domain_data &data, real_domain_data &real){
  
  int N_real=real.N_real();
  int N_imag=data.N_imag();
  
  //extract Matsubara frequencies, Matsubara data, and real frequencies
  pade_complex_vector_type real_frequencies(N_real);
  pade_complex_vector_type matsubara_frequencies(N_imag);
  pade_complex_vector_type matsubara_values(N_imag);
  for(int i=0;i<N_real;++i){
    real_frequencies[i]=pade_complex_type(real.freq()[i], 0);
  }
  for(int i=0;i<N_imag;++i){
    matsubara_frequencies[i]=pade_complex_type(0.,data.freq()[i]);
    matsubara_values     [i]=pade_complex_type(data.val()[i].real(), data.val()[i].imag());
  }
  
  pade_complex_vector_type alpha_sx(N_imag);
  pade_complex_vector_type phi_mu_nu(N_imag);
  double norm=data.norm();
  
  
  ////DEBUG
  int N=21;
  pade_complex_vector_type input_x(N), input_y(N);
  for(int i=0;i<N;++i){ input_x[i]=0.1*i; input_y[i]=pade_complex_type(sin(0.1*i)); }
  pade_complex_vector_type pmn(N);
  pade_mu_=20;
  pade_nu_=0;
  for(pade_complex_type x=pade_complex_type(-2.); x.real()<pade_complex_type(-1.999).real(); x=x+pade_complex_type(0.2)){
    //pade_complex_type x=3.001;
    pade_complex_vector_type asx(N); for(int i=0;i<N;++i) asx[i]=x-input_x[i];
    //for(int i=0;i<N;++i){ std::cout<<input_x[i]<<" "<<input_y[i]<<" "<<asx[i]<<std::endl;}
    pade_complex_type z=pade_scheme(pmn, input_x, input_y, asx, x);
    std::cout<<x.real()<<" "<<z.real()<<" "<<z.imag()<<std::endl;
  }
  exit(0);
  ////DEBUG
  
  
  
  //interpolate 'alpha', the distance from each point to x
#pragma omp parallel for default(none) shared(real,std::cerr) firstprivate(N_imag,N_real,real_frequencies,alpha_sx,matsubara_frequencies,matsubara_values,phi_mu_nu, norm)
  for(int i=0;i<N_real;++i){
    pade_complex_type x=real_frequencies[i];
    std::pair<bool, int> v=init_alpha(N_imag, alpha_sx,matsubara_frequencies, x);
    if(v.first){ //tried to evaluate on a point where we know the function
      real.val()[i]=to_simple_precision(matsubara_values[v.second]);
      continue;
    }
    
    //run neville scheme
#pragma omp critical
    std::cerr<<"computing pade for i: "<<i<<std::endl;
    real.val()[i]=to_simple_precision(pade_scheme(phi_mu_nu, matsubara_frequencies, matsubara_values, alpha_sx, x))*norm;
  }
}

void pade_interpolator::find_epsilon(){
  int digits = 0;
  pade_real_type x, y=1.;
  do {
    y/=2.;
    x=y+1.;
    if(x==1.) break;
    digits++;
  }while(1);
  epsilon_=y;
  std::cout<<"epsilon is: "<<epsilon_<<" digits are: "<<digits<<std::endl;
}

