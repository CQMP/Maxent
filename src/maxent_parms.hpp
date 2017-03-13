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
#include "default_model.hpp"
#include "maxent_grid.hpp"
#include "maxent_kernel.hpp"
#include "maxent_matrix_def.hpp"

///This class has all the information about general, non-SVD specific analytic continuation things.
///It constructs kernels, grids, etc.
class KernelAndGrid {

public:
  
  ///constructs the kernel and grid from the parameters p. Also reads in the data.
  KernelAndGrid(alps::params& p);
  ///define parameter defaults
  static void define_parameters(alps::params &p);
  
  ///value of the Matsubara data at index i
  double y(const int i) const { return y_[i]; }
  ///value of the Matsubara data
  const vector_type& y() const { return y_; }
  ///value of the covariance matrix
  double cov(const int i,const int j) const { return cov_(i,j); }
  ///value of the error (if covariance not given)
  double sigma(const int i) const { return sigma_[i]; }
  ///value of the inverse temperature
  double T() const { return T_; }
  ///number of data points
  int ndat() const { return ndat_; }
  ///number of frequency points on the real axis
  int nfreq() const { return nfreq_; }
  ///value of the Kernel matrix K(i,j).
  double K(const int i, const int j) const {return K_(i,j);}
  ///returns the entire kernel matrix
  const matrix_type& K() const { return K_; }
  ///returns kernel type of K
  kernel_type getKernelType() const {return k_type; }

private:
  ///temperature
  const double T_;

  void read_data_from_text_file(const alps::params& p);
  void read_data_from_hdf5_file(const alps::params& p);
  void read_data_from_param_file(const alps::params& p);
  void read_covariance_matrix_from_text_file(const std::string& fname);

protected:

  void decompose_covariance_matrix(const alps::params& p);
  ///This function scales both the y (data) and the kernel with the errors
  void scale_data_with_error(const int ntab);

  ///number of fitting data points / size of y
  int ndat_; //TODO: find a better way of changing this. If we use Gl, we change the size of y
  ///number of real frequencies
  const int nfreq_;
  ///vector of Matsubara data
  vector_type y_;
  ///vector of errors on y
  vector_type sigma_;
  ///Kernel matrix
  matrix_type K_;
  ///covariance matrix
  matrix_type cov_;
  ///real frequency grid
  grid grid_;
  ///vector containing input matsubara or tau data
  vector_type inputGrid_;
  ///type of kernel used
  kernel_type k_type;
};



///This class contains all of the information about singular value decompositions.
///It inherits from the KernelAndGrid class that was responsible for building
///the kernel.
class SVDContinuation : public KernelAndGrid
{
public:
  ///constructs and singular value decomposes the kernel
  SVDContinuation(alps::params& p);
  ///define parameter defaults
  static void define_parameters(alps::params &p);
  
  ///SVD decomposes the kernel into U*S*V^T. This is U
  const matrix_type& U() const { return U_; }
  ///SVD decomposes the kernel into U*S*V^T. This is V^T
  const matrix_type& Vt() const { return Vt_; }
  ///SVD decomposes the kernel into U*S*V^T. This is the vector Sigma that forms the diagonal of S
  const matrix_type& Sigma() const { return Sigma_; }
  double omega_coord(const int i) const { return omega_coord_[i]; }
  const vector_type& omega_coord() const { return omega_coord_; }
  double delta_omega(const int i) const { return delta_omega_[i]; }
  ///getter function for the number of singular values
  int ns() const { return ns_; }
  ///getter function for the default model
  const DefaultModel& Default() const { return *Default_; }
  ///getter function for input data grid
  const vector_type& inputGrid() const { return inputGrid_; }
  ///getter function for input data grid values
  double inputGrid(const int i) const { return inputGrid_(i); }

private:
  ///The default model
  boost::shared_ptr<DefaultModel> Default_;
  matrix_type U_;
  matrix_type Vt_;
  matrix_type Sigma_;
  vector_type omega_coord_;
  vector_type delta_omega_;
  ///the number of singular values
  int ns_;
  ///compute the minimal chi2 that is possible given the SVD of the kernel
  void compute_minimal_chi2() const;
  ///reduce the matrices to the number of non-zero singular values
  void truncate_to_singular_space(const vector_type& S);
  ///take the kernel and compute its singular value decomposition
  void singular_value_decompose_kernel(bool verbose, vector_type& S);
  ///check G~-1/iw_{n} and default model back continues same limit
  void check_high_frequency_limit(const vector_type& y, const kernel_type kt);
};
