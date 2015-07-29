/*
 * Copyright (C) 1998-2015 ALPS Collaboration. See COPYRIGHT.TXT
 * All rights reserved. Use is subject to license terms. See LICENSE.TXT
 * For use in publications, see ACKNOWLEDGE.TXT
 */
#include<iostream>
#include<vector>
#include<cmath>
#include<stdexcept>

#include <alps/numeric/matrix.hpp>
#include<gmpxx.h>
typedef mpf_class real_t;
typedef alps::numeric::matrix<real_t> matrix_t;
typedef std::vector<real_t> vector_t;
extern "C" void dgesv_(const int *N, const int *NRHS, double *A, const int *LDA, int *IPIV, double *B, const int *LDB, int *INFO);
//void solve(const matrix_t &A, const vector_t &rhs, vector_t &res);
template<typename T> void solve(const alps::numeric::matrix<T> &A, const std::vector<T> &rhs, std::vector<T>  &res);
void init_example1(std::vector<real_t> &s, std::vector<real_t> &f, const int &N){
  if(N!=5) throw std::logic_error("this is done for the example in the paper which has N=5");
  s.resize(N);
  f.resize(N);
  
  s[0]=-2;
  s[1]=-1;
  s[2]= 0;
  s[3]= 1;
  s[4]= 2;
  
  f[0]=-2;
  f[1]=-1;
  f[2]=-1;
  f[3]= 0;
  f[4]= 1;
}
void init_sin(std::vector<real_t> &s, std::vector<real_t> &f, const int &N){
  s.resize(N);
  f.resize(N);
  
  for(int i=0;i<N;++i){
    s[i]=0.2*i;
    f[i]=cos(0.2*i);
  }
}
//see equation 3
real_t compute_omega_j(int j, const std::vector<real_t> &s, const real_t &x0){
  if(j==0) return 1;
  real_t omega=1;
  for(int i=0;i<j;++i){
    omega*=(x0-s[i]);
  }
  return omega;
}
//see equation 6
real_t compute_omegaprime_kp1_of_xj(int k, int j, const std::vector<real_t> &s){
  real_t omega=1;
  for(int i=0;i<=k;++i){
    if(i==j) continue;
    omega*=(s[j]-s[i]);
  }
  std::cout<<"returning omegaprime k+1 for k: "<<k<<" and j: "<<j<<" as: "<<omega<<std::endl;
  return omega;
}
void fill_Vinv_matrix( matrix_t &Vinv, const std::vector<real_t> &s, int m, int n){
  for(int k=0;k<=m+n;++k){
    for(int j=0;j<k;++j){
      Vinv(k,j)=Vinv(k-1,j)/(s[j]-s[k]);
    }
    Vinv(k,k)=1./compute_omegaprime_kp1_of_xj(k,k,s);
  }
  std::cout<<"Vinv is: "<<std::endl;
  std::cout<<Vinv<<std::endl;
}
void assemble_matrix_system(matrix_t &Lambda, std::vector<real_t> &rhs, const matrix_t &Vinv, const std::vector<real_t> &f, int m, int n){
  matrix_t M(n+m,n+m+1, 0.);
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
  std::cout<<M<<std::endl;
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
real_t evaluate_bary_poly(const std::vector<real_t> &q, const std::vector<real_t> &f, const std::vector<real_t> &x, const matrix_t &Vinv, int m, int n, const real_t &x0){
  //evaluate numerator
  real_t numerator=f[0]*q[0];
  for(int k=1;k<=m;++k){
    real_t dfq=0.;
    for(int i=0;i<=k;++i){
      dfq+=Vinv(k, i)*q[i]*f[i];
    }
    numerator+=dfq*compute_omega_j(k, x, x0);
  }
  
  //evaluate denominator
  real_t denominator=q[0];
  for(int k=1;k<=n;++k){
    real_t dq=0.;
    for(int i=0;i<=k;++i){
      dq+=Vinv(k, i)*q[i];
    }
    denominator+=dq*compute_omega_j(k, x, x0);
  }
  //std::cout<<"numerator: "<<numerator<<" denominator: "<<denominator<<std::endl;
  return numerator/denominator;
}
template<> void solve(const alps::numeric::matrix<double> &Lambda, const std::vector<double> &rhs, std::vector<double> &res){
  std::cout<<"double precision solver. "<<std::endl;
  int size=rhs.size();
  int one =1;
  int info;
  std::vector<int> ipiv(size);
  alps::numeric::matrix<double> Lambda2(Lambda);
  res=rhs;
  dgesv_(&size, &one, &(Lambda2(0,0)), &size, &(ipiv[0]), &(res[0]), &size, &info);
}
int main(){
  
  //N points on which function is known
  int N=512;
  int n=255;  //degree of denominator
  int m=256; //degree of numerator
  mpf_set_default_prec(2048);
  //function values are f, points of support are x
  std::vector<real_t> s, f;
  //initialize function f
  //init_example1(s,f, N);
  init_sin(s,f, N);
  
  //compute Vinv:
  matrix_t Vinv(m+n+1,m+n+1);
  std::cout<<"filling"<<std::endl;
  fill_Vinv_matrix(Vinv, s, m, n);
    
  matrix_t Lambda(m+n,m+n);
  std::vector<real_t> rhs(m+n);
  std::vector<real_t> q(m+n+1);
  std::cout<<"assembling"<<std::endl;
  assemble_matrix_system(Lambda, rhs, Vinv, f, m, n);
  matrix_t Lambda_bup(Lambda);
  std::vector<real_t> res(rhs.size(), 0.);

  std::cout<<"right hand side: "<<std::endl;
  for(int i=0;i<rhs.size();++i){
    std::cout<<i<<" "<<rhs[i]<<std::endl;
  }

  std::cout<<"solving: "<<std::endl;
  solve(Lambda, rhs,res);
  
  //coefficients q, q[0] is set to 1.
  q[0]=1;
  for(int i=0;i<res.size();++i){
    q[i+1]=res[i];
  }
  std::cout<<"qs are: "<<std::endl;
  for(int i=0;i<q.size();++i){
    std::cout<<i<<" "<<q[i]<<std::endl;
  }
  std::cout<<std::endl;

  for(int i=0;i<N;++i){
    std::cout<<s[i]<<" "<<f[i]<<std::endl;
  }
  std::cout<<std::endl;
  
  //assemble polynomial:
  for(real_t x0=-25;x0<120;x0+=0.1){
    //std::cout<<"evaluating"<<std::endl;
    std::cout<<x0<<" "<<evaluate_bary_poly(q, f, s, Vinv,  m,  n, x0)<<std::endl;
  }
  std::cout<<std::endl;

}

