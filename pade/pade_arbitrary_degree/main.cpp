/*
 * Copyright (C) 1998-2015 ALPS Collaboration. See COPYRIGHT.TXT
 * All rights reserved. Use is subject to license terms. See LICENSE.TXT
 * For use in publications, see ACKNOWLEDGE.TXT
 */
#include<iostream>
#include<vector>
#include<cmath>
#include<stdexcept>

void init(std::vector<double> &s, std::vector<double> &f, const int &N){
  s.resize(N);
  f.resize(N);
  for(int i=0;i<N;++i){
    s[i]=0.1*i+0.05;
    f[i]=sin(s[i]);
  }
}

//compute 'alpha' = x - xs (distance of interpolation to grid point). If we find an interpolation point identical to the grid point, there is no need to interpolate and we return true.
bool init_alpha(int N, std::vector<double> &alpha_sx, const std::vector<double> &s, const std::vector<double> &fs, const double &x){
  for(int i=0;i<N;++i){
    alpha_sx[i]=x-s[i];
    if(alpha_sx[i]==0){
      std::cout<<x<<" "<<fs[i]<<std::endl;
      return true;
    }
  }
  return false;
}
//a single step in the Neville scheme, see above Eq. 2.2.3.7
void neville_step(std::vector<double> &phi_mu, const std::vector<double> &phi_mu_minus_one, const std::vector<double> &alpha_sx, int mu){
  int N=phi_mu.size();
  //std::cout<<"running : "<<N-mu<<" points for neville "<<std::endl;
  for(int n=0;n<N-mu;n++){
    phi_mu[n]=(alpha_sx[n]*phi_mu_minus_one[n+1]-alpha_sx[n+mu]*phi_mu_minus_one[n])/(alpha_sx[n]-alpha_sx[n+mu]);
    //std::cout<<" numerator: "<<(alpha_sx[n]*phi_mu_minus_one[n+1]-alpha_sx[n+mu]*phi_mu_minus_one[n])<<" denominator: "<<(alpha_sx[n]-alpha_sx[n+mu])<<std::endl;
  }
}
//a single step in the inverse Neville scheme, see Eq. 2.2.3.7
void inverse_neville_step(std::vector<double> &phi_nu, const std::vector<double> &phi_nu_minus_one, const std::vector<double> &alpha_sx, int nu){
  int N=phi_nu.size();
  //std::cout<<"running "<<N-nu<<" points for inverse neville."<<std::endl;
  for(int n=0;n<N-nu;n++){
    phi_nu[n]=(alpha_sx[n]-alpha_sx[n+nu])/(alpha_sx[n]/phi_nu_minus_one[n+1]-alpha_sx[n+nu]/phi_nu_minus_one[n]);
    //std::cout<<"numerator: "<<(alpha_sx[n]-alpha_sx[n+nu])<<" denominator: "<<(alpha_sx[n]/phi_nu_minus_one[n+1]-alpha_sx[n+nu]/phi_nu_minus_one[n])<<std::endl;
  }
}

//this function implements the Neville scheme for interpolating polynomials of degree mumax through point x
void poly_neville_scheme(std::vector<double> &phi_mu, const std::vector<double> &s, const std::vector<double> &fs, const std::vector<double> &alpha_sx, const double &x, const int &mumax){
  phi_mu=fs;
  std::vector<double> phi_mu_minus_one(phi_mu);
  for(int mu=1;mu<=mumax;++mu){
    phi_mu.swap(phi_mu_minus_one);
    neville_step(phi_mu, phi_mu_minus_one, alpha_sx, mu);
  }
  std::cout<<x<<" "<<phi_mu[0]<<std::endl;
}
//this function implements the Neville scheme for interpolating rational functions of degree mumax through point x
void poly_inverse_neville_scheme(std::vector<double> &phi_nu, const std::vector<double> &s, const std::vector<double> &fs, const std::vector<double> &alpha_sx, const double &x, const int &numax){
  phi_nu=fs;
  std::vector<double> phi_nu_minus_one(phi_nu);
  for(int nu=1;nu<=numax;++nu){
    phi_nu.swap(phi_nu_minus_one);
    inverse_neville_step(phi_nu, phi_nu_minus_one, alpha_sx, nu);
  }
  std::cout<<x<<" "<<phi_nu[0]<<std::endl;
}

void pade_scheme(std::vector<double> &phi_mu_nu, const std::vector<double> &s, const std::vector<double> &fs, const std::vector<double> &alpha_sx, const double &x, const int &mumax, const int &numax){
  phi_mu_nu=fs;
  std::vector<double> phi_mu_nu_minus_one(phi_mu_nu);
  int nu=0;
  int mu=0;
  
  if(numax+mumax+1 != phi_mu_nu.size()) throw std::runtime_error("problem in interpolation: degree of numerator and denominator+1 has to match degree of data.");
  
  //start by doing a nu-mu inverse neville steps, until 
  for(int interp_step=1;interp_step<=numax+mumax;++interp_step){
    phi_mu_nu.swap(phi_mu_nu_minus_one);
    if(interp_step < numax-mumax+1){
      inverse_neville_step(phi_mu_nu, phi_mu_nu_minus_one, alpha_sx, interp_step);
      //std::cout<<"# running inverse step (beginning)"<<mu<<" "<<nu<<std::endl;
      nu++;
    }else if(interp_step < mumax-numax+1){
      neville_step(phi_mu_nu, phi_mu_nu_minus_one, alpha_sx, interp_step);
      //std::cout<<"# running step (beginning)"<<mu<<" "<<nu<<std::endl;
      mu++;
    }else if (interp_step%2==1){
      //std::cout<<"# running inverse step (intermed)"<<mu<<" "<<nu<<std::endl;
      inverse_neville_step(phi_mu_nu, phi_mu_nu_minus_one, alpha_sx, interp_step);
      nu++;
    }else{
      //std::cout<<"# running step (intermed)"<<mu<<" "<<nu<<std::endl;
      neville_step(phi_mu_nu, phi_mu_nu_minus_one, alpha_sx, interp_step);
      mu++;
    }
  }
  std::cout<<x<<" "<<phi_mu_nu[0]<<std::endl;
}


int main(){
  
  //N points on which function is known
  int N=62;
  //function values are fs, points of support are s
  std::vector<double> fs, s;
  //initialize function f
  init(s,fs, N);
  
  for(int i=0;i<N;++i){
    std::cout<<s[i]<<" "<<fs[i]<<std::endl;
  }
  std::cout<<std::endl;

  for(int q=0;q<1000;++q){
    //evaluation point x
    double x=-0.5+(q/100.);
    //double x=3.15;
    //compute distance to point x
    std::vector<double> alpha_sx(N);
    if(init_alpha(N, alpha_sx, s, fs, x)) continue;
    
    std::vector<double> phi_mu(N);
    pade_scheme(phi_mu, s, fs, alpha_sx, x, 30,31);
  }
  std::cout<<std::endl;
  
}

