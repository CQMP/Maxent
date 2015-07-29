/*
 * Copyright (C) 1998-2015 ALPS Collaboration. See COPYRIGHT.TXT
 * All rights reserved. Use is subject to license terms. See LICENSE.TXT
 * For use in publications, see ACKNOWLEDGE.TXT
 */

#include "pade.hpp"
#include <fstream>


pade_interpolator::pade_interpolator(const PadeParams &p){
  pade_mu_=p["pade.PADE_NUMERATOR_DEGREE"];
  pade_nu_=p["pade.PADE_DENOMINATOR_DEGREE"];
  if(pade_mu_+pade_nu_ +1 != p["imag.NDAT"].as<int>()) throw std::runtime_error("Ill defined interpolation: please make sure that numerator degree + denominator degree + 1 = # of data points");
  mpf_set_default_prec(p["pade.FLOAT_PRECISION"]);
  find_epsilon();
}


//An implementation of the pade interpolation according to X. Zhu and G. Zhu, J. Comp. and Appl. Math 148, 341 (2002)
void pade_interpolator::pade_interpolate(const imaginary_domain_data &data, real_domain_data &real)const{
  
  int N_real=real.N_real();
  int N_imag=data.N_imag();
  
  //extract Matsubara frequencies, Matsubara data, and real frequencies
  pade_complex_vector_type real_frequencies(N_real);
  //pade_complex_vector_type imag_frequencies(N_real);
  pade_complex_vector_type matsubara_frequencies(N_imag);
  pade_complex_vector_type matsubara_values(N_imag);
  for(int i=0;i<N_real;++i){
    real_frequencies[i]=pade_complex_type(real.freq()[i], 0.);
    //imag_frequencies[i]=pade_complex_type(0., real.freq()[i]);
  }
  for(int i=0;i<N_imag;++i){
    matsubara_frequencies[i]=pade_complex_type(0.,data.freq()[i]);
    matsubara_values     [i]=pade_complex_type(data.val()[i].real(), data.val()[i].imag());
  }

  std::cerr<<"filling."<<std::endl;
  pade_complex_matrix_type Vinv(pade_mu_+pade_nu_+1,pade_mu_+pade_nu_+1);
  fill_Vinv_matrix(Vinv, matsubara_frequencies, pade_mu_, pade_nu_);
  
  std::cerr<<"assembling."<<std::endl;
  pade_complex_matrix_type Lambda(pade_mu_+pade_nu_,pade_mu_+pade_nu_);
  pade_complex_vector_type rhs(pade_mu_+pade_nu_);
  assemble_matrix_system(Lambda, rhs, Vinv, matsubara_values, pade_mu_, pade_nu_);

  std::cerr<<"solving."<<std::endl;
  pade_complex_vector_type res(rhs.size(), pade_complex_type(0., 0.));
  pade_complex_vector_type q(pade_mu_+pade_nu_+1);
  pade_solver s;
  s.solve(Lambda, rhs,res);
  //coefficients q, q[0] is set to 1.
  q[0]=1;
  for(int i=0;i<res.size();++i){
    q[i+1]=res[i];
  }
  
  std::cerr<<"evaluating."<<std::endl;
#pragma omp parallel for default(none) shared(real) firstprivate(q,matsubara_values,matsubara_frequencies,Vinv,real_frequencies,N_real)
  for(int k=0;k<N_real;++k){
    //real.val()[k]=to_simple_precision(evaluate_bary_poly(q, matsubara_values, matsubara_frequencies, Vinv,  pade_mu_,  pade_nu_, imag_frequencies[k]));
    real.val()[k]=to_simple_precision(evaluate_bary_poly(q, matsubara_values, matsubara_frequencies, Vinv,  pade_mu_,  pade_nu_, real_frequencies[k]));
  }

}
pade_complex_type pade_interpolator::evaluate_bary_poly(const pade_complex_vector_type &q, const pade_complex_vector_type &f, const pade_complex_vector_type &x, const pade_complex_matrix_type &Vinv, int m, int n, const pade_complex_type &x0)const{
  //evaluate numerator
  pade_complex_type numerator=f[0]*q[0];
  for(int k=1;k<=m;++k){
    pade_complex_type dfq=pade_complex_type(0., 0.);
    for(int i=0;i<=k;++i){
      dfq+=Vinv(k, i)*q[i]*f[i];
    }
    numerator+=dfq*compute_omega_j(k, x, x0);
  }
  
  //evaluate denominator
  pade_complex_type denominator=q[0];
  for(int k=1;k<=n;++k){
    pade_complex_type dq=pade_complex_type(0., 0.);
    for(int i=0;i<=k;++i){
      dq+=Vinv(k, i)*q[i];
    }
    denominator+=dq*compute_omega_j(k, x, x0);
  }
  return numerator/denominator;
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
  std::cout<<"epsilon is: "<<epsilon_.get_d()<<" digits are: "<<digits<<std::endl;
}

pade_complex_type pade_interpolator::compute_omega_j(int j, const pade_complex_vector_type &s, const pade_complex_type &x0)const {
  if(j==0) return pade_complex_type(1., 0.);
  pade_complex_type omega=pade_complex_type(1., 0.);
  for(int i=0;i<j;++i){
    omega*=(x0-s[i]);
  }
  return omega;
}

void pade_interpolator::fill_Vinv_matrix(pade_complex_matrix_type &Vinv, const pade_complex_vector_type &s, int m, int n) const{
  for(int k=0;k<=m+n;++k){
    for(int j=0;j<k;++j){
      Vinv(k,j)=Vinv(k-1,j)/(s[j]-s[k]);
    }
    Vinv(k,k)=pade_complex_type(1., 0.)/compute_omegaprime_kp1_of_xj(k,k,s);
  }
}
void pade_interpolator::assemble_matrix_system(pade_complex_matrix_type &Lambda, pade_complex_vector_type &rhs, const pade_complex_matrix_type &Vinv, const pade_complex_vector_type &f, int m, int n) const{
  pade_complex_matrix_type M(n+m,n+m+1);
  //pade_complex_matrix_type M(n+m,n+m+1, pade_complex_type(0., 0.));
  //assemble upper half
  //std::cout<<"assembling upper half"<<std::endl;
  for(int k=0;k<m;++k){
    for(int j=0;j<=n+1+k;++j){
      M(k,j)=Vinv(k+n+1, j);
    }
  }
  //std::cout<<"assembling lower half"<<std::endl;
  for(int k=0;k<n;++k){
    for(int j=0;j<=m+1+k;++j){
      M(m+k,j)=Vinv(m+k+1, j)*f[j];
      //if(j==0) std::cout<<"part of M(m+k,0): "<<M(m+k,0)<<" "<<Vinv(m+k+1,0)<<" "<<f[0]<<std::endl;
    }
  }
  //std::cout<<M<<std::endl;
  //std::cout<<"decomposing into left and right"<<std::endl;
  for(int i=0;i<n+m;++i){
    for(int j=1;j<n+m+1;++j){
      Lambda(i,j-1)=M(i,j);
    }
    rhs[i]=-M(i,0); //this encodes the case q0=1, which is the standard case we're interested in.
  }
  /*for(int i=0;i<rhs.size();++i){
   std::cout<<i<<" rhs: "<<rhs[i]<<std::endl;
   }
   for(int i=0;i<f.size();++i){
   std::cout<<i<<" "<<Vinv(i, 0)<<" "<<f[i]<<std::endl;
   }*/
  
}
pade_complex_type pade_interpolator::compute_omegaprime_kp1_of_xj(int k, int j, const pade_complex_vector_type &s) const{
  pade_complex_type omega=pade_complex_type(1., 0.);
  for(int i=0;i<=k;++i){
    if(i==j) continue;
    omega*=(s[j]-s[i]);
  }
  //std::cout<<"returning omegaprime k+1 for k: "<<k<<" and j: "<<j<<" as: "<<omega<<std::endl;
  return omega;
}


