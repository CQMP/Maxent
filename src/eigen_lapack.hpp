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
#pragma once

#ifdef HAVE_LAPACK
extern "C" void dgesvd_( const char* jobu, const char* jobvt,
                                              const int* m, const int* n, double* a, const int* lda,
                                              double* s, double* u, const int* ldu,
                                              double* vt, const int* ldvt,
                                              double* work, const int* lwork, int* info);

///performs a SVD on input matrix K
/// returns K =  U (S) V^T where S is a vector of the diagonal values of
/// matrix Sigma
void lapack_svd(matrix_type &K, vector_type &S, matrix_type &Vt, matrix_type &U){

   /*boost::numeric::bindings::lapack::gesvd('S', 'S', Kt, S, U_, Vt_);*/
  //use 'S' for thin U,Vt matrices
  char jobu = 'S';
  char jobvt = 'S';
  int m = K.rows();
  int n = K.cols();
  double *a = K.data();
  int lda = K.outerStride();
  int ldu = U.outerStride();
  int ldvt = Vt.outerStride();

  double *vt = Vt.data();
  double *u = U.data();
  double *s =S.data();

  double dummywork;
  int LWORK=-1;
  int INFO=0;
  //get best workspace
  dgesvd_(&jobu,&jobvt,&m,&n,a,&lda,s,u,&ldu,vt,&ldvt,&dummywork,&LWORK,&INFO);

  LWORK=int(dummywork)+32;
  vector_type WORK(LWORK);

  dgesvd_(&jobu,&jobvt,&m,&n,a,&lda,s,u,&ldu,vt,&ldvt,WORK.data(),&LWORK,&INFO);
  if(INFO!=0)
    std::cerr<<"Warning, danger! SVD failed to converge" <<std::endl;



}
#endif
