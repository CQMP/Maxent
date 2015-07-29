/*
 * Copyright (C) 1998-2015 ALPS Collaboration. See COPYRIGHT.TXT
 * All rights reserved. Use is subject to license terms. See LICENSE.TXT
 * For use in publications, see ACKNOWLEDGE.TXT
 */

#pragma once

#include <alps/config.hpp> // needed to set up correct bindings
#include "maxent_config.hpp"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/bindings/ublas.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/bindings/blas/level3/gemm.hpp>
#include <boost/numeric/bindings/blas/level2/gemv.hpp>



typedef boost::numeric::ublas::matrix<double, boost::numeric::ublas::column_major> matrix_type;
typedef boost::numeric::ublas::vector<double> vector_type;
typedef boost::numeric::ublas::vector<std::complex<double> > complex_vector_type;
typedef std::pair<vector_type, complex_vector_type> omega_complex_type;

///matrix-vector multiplication. TODO: this should be delegated to a matrix library
inline vector_type maxent_prec_prod(const matrix_type &p, const vector_type &q) {
#ifdef HAVE_LAPACK
  vector_type r(p.size1());
  double alpha=1.;
  double beta=0.;
  boost::numeric::bindings::blas::gemv( alpha, p,q,beta,r);
  return r;
#else
  return boost::numeric::ublas::prec_prod(p,q);
#endif
}
///matrix-vector multiplication of transpose of matrix. TODO: this should be delegated to a matrix library
inline vector_type maxent_prec_prod_trans(const matrix_type &p, const vector_type &q) {
#ifdef HAVE_LAPACK
  vector_type r(p.size2());
  double alpha=1.;
  double beta=0.;
  boost::numeric::bindings::blas::gemv( alpha, boost::numeric::ublas::trans(p),q,beta,r);
  return r;
#else
  return boost::numeric::ublas::prec_prod(boost::numeric::ublas::trans(p),q);
#endif
}
///matrix-matrix multiplication. TODO: this should be delegated to a matrix library
inline matrix_type maxent_prec_prod(const matrix_type &p, const matrix_type &q) {
#ifdef HAVE_LAPACK
  matrix_type r(p.size1(), q.size2());
  double alpha=1.;
  double beta=0.;
  boost::numeric::bindings::blas::gemm(alpha, p, q, beta, r);
  return r;
#else
  return boost::numeric::ublas::prec_prod(p,q);
#endif
}
///matrix-matrix multiplication of transpose(p) with q. TODO: this should be delegated to a matrix library
inline matrix_type maxent_prec_prod_trans(const matrix_type &p, const matrix_type &q) {
#ifdef HAVE_LAPACK
  matrix_type r(p.size2(), q.size2());
  double alpha=1.;
  double beta=0.;
  boost::numeric::bindings::blas::gemm(alpha, boost::numeric::ublas::trans(p), q, beta, r);
  return r;
#else
  return boost::numeric::ublas::prec_prod(boost::numeric::ublas::trans(p),q);
#endif
}
