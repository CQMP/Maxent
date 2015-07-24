/*****************************************************************************
 *
 * ALPS Project Applications
 *
 * Copyright (C) 2010 by Sebastian Fuchs <fuchs@comp-phys.org>
 *                       Thomas Pruschke <pruschke@comp-phys.org>
 *                       Matthias Troyer <troyer@comp-phys.org>
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

#pragma once

#include <alps/config.hpp> // needed to set up correct bindings
#include "maxent_config.hpp"
#include <Eigen/Core>



typedef Eigen::MatrixXd matrix_type;
typedef Eigen::MatrixXcd complex_matrix_type;
typedef Eigen::VectorXd vector_type;
typedef Eigen::MatrixXcd complex_vector_type;
typedef std::pair<vector_type, complex_vector_type> omega_complex_type;

///matrix-vector multiplication. TODO: this should be delegated to a matrix library
inline vector_type maxent_prec_prod(const matrix_type &p, const vector_type &q) {
	return p*q;
}
///matrix-vector multiplication of transpose of matrix. TODO: this should be delegated to a matrix library
inline vector_type maxent_prec_prod_trans(const matrix_type &p, const vector_type &q) {
	return p.transpose()*q;
}
///matrix-matrix multiplication. TODO: this should be delegated to a matrix library
inline matrix_type maxent_prec_prod(const matrix_type &p, const matrix_type &q) {
	return p*q;
}
///matrix-matrix multiplication of transpose(p) with q. TODO: this should be delegated to a matrix library
inline matrix_type maxent_prec_prod_trans(const matrix_type &p, const matrix_type &q) {
	return p.transpose()*q; 
}
