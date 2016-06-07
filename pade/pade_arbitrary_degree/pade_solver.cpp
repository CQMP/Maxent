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
