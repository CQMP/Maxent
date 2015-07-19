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

#include "pade.hpp"
#include <fstream>


void pade_solver::backsub_lower(const pade_complex_matrix_type &Linv, const pade_complex_vector_type &rhs, pade_complex_vector_type &res){
  int N=Linv.rows();
  for(int i=0;i<N;++i){
    pade_complex_type tmp(0.,0.);
    for(int j=0;j<i;++j){
      tmp+=Linv(i,j)*res[j];
    }
    res[i]=(rhs[i]-tmp)/Linv(i,i);
  }
}
void pade_solver::backsub_upper(const pade_complex_matrix_type &Uinv, const pade_complex_vector_type &rhs, pade_complex_vector_type &res){
  int N=Uinv.rows();
  for(int i=N-1;i>=0;--i){
    pade_complex_type tmp(0.,0.);
    for(int j=i+1;j<N;++j){
      tmp+=Uinv(i,j)*res[j];
    }
    res[i]=(rhs[i]-tmp)/Uinv(i,i);
  }
}

//run an LU decomposition with pivoting, then back-substitute
void pade_solver::solve(const pade_complex_matrix_type &A, const pade_complex_vector_type &rhs, pade_complex_vector_type &res){
  std::cerr<<"multiprecision solver."<<std::endl;
  int N=A.rows();
  pade_complex_matrix_type P=pade_complex_matrix_type::Identity(N,N);
  //pade_complex_matrix_type P(pade_complex_matrix_type::identity_matrix(N));
  pade_complex_matrix_type L(P);
  pade_complex_matrix_type U(A);
  std::vector<int> pivot(N); for(int i=0;i<N;++i) pivot[i]=i;

  int kprime;
  for(int k=0;k<N;++k){
    pade_real_type p(0.);
    for(int i=k;i<N;++i){
      if(std::abs(U(i,k))>p){
        p=abs(U(i,k));
        kprime=i;
      }
    }
    if(p==0.) throw std::runtime_error("encountered singular matrix!");
    if(k != kprime){
      
      std::swap(pivot[k],pivot[kprime]);
      for(int i=k;i<N;++i){
        std::swap(U(k,i), U(kprime, i));
      }
      for(int i=0;i<k;++i){
        std::swap(L(k,i), L(kprime, i));
      }
      //std::cout<<"Pivot U is: "<<U<<std::endl;
    }
    for(int j=k+1;j<N;++j){
      L(j,k)=U(j,k)/U(k,k);
      for(int l=k;l<N;++l){
        U(j,l)-=L(j,k)*U(k,l);
      }
    }
  }
  pade_complex_vector_type tmp(N, pade_complex_type(0.,0.)), rhs_pivot(rhs);
  //pivot
  for(int i=0;i<N;++i){
    if(i<pivot[i]){
      std::swap(rhs_pivot[i], rhs_pivot[pivot[i]]);
    }
  }
  std::cerr<<"backsubstitution."<<std::endl;
  
  backsub_lower(L, rhs_pivot, tmp);
  backsub_upper(U, tmp, res);
  std::cerr<<"matrix equation solved."<<std::endl;
}
