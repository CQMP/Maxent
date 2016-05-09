/*
 * Copyright (C) 1998-2016 ALPS Collaboration. See COPYRIGHT.TXT
 * All rights reserved. Use is subject to license terms. See LICENSE.TXT
 * For use in publications, see ACKNOWLEDGE.TXT
 */
#pragma once
#include<vector>
#include<alps/params.hpp>
#include"maxent_matrix_def.hpp"
///enum that enumerates if we're in time or in frequency
enum dataspace_type{
  time_dataspace,
  frequency_dataspace,
  legendre_dataspace,
  real_dataspace
};
///we have a range of different kernels that all need different input parameters. This enum enumerates them.
enum kernel_type{
  time_fermionic_kernel,
  time_bosonic_kernel,
  time_boris_kernel,
  time_fermionic_legendre_kernel,
  time_bosonic_legendre_kernel, //TODO determine boris->legendre
  legendre_fermionic_kernel,
  legendre_bosonic_kernel,
  frequency_fermionic_ph_kernel,
  frequency_bosonic_ph_kernel,
  frequency_anomalous_ph_kernel,
  frequency_fermionic_kernel,
  frequency_bosonic_kernel,
  frequency_anomalous_kernel,
  real_kernel
};
class kernel{
public:
  kernel(alps::params &p, const vector_type& freq, vector_type &inputGrid);

  ///getter function for the kernel matrix
  const matrix_type &operator()()const{return K_;}
  ///getter function for kernel type
  const kernel_type getKernelType() const{return ktype_;}
  ///getter function for dataspace type
  const dataspace_type getDataspaceType() const{return dtype_;}
	///getter function for tau_points
	const vector_type getTauPoints() const{return tau_points_;}
    
private:
  ///figure out which kernel is to be used
  void set_kernel_type(const std::string &dataspace_name, const std::string &kernel_name,
                       bool ph_symmetry);
  ///set up kernel with the legendre transform
  void setup_legendre_kernel(const alps::params &p, const vector_type& freq, const int lmax);
  ///number; of Matsubara points
  int ndat_;
  ///number of real frequency points
  int nfreq_;
  ///temperature
  double T_;
  ///type of the kernel (see enums above)
  kernel_type ktype_;
  ///type of the dataspace
  dataspace_type dtype_;
  ///kernel matrix
  matrix_type K_;
	///array of tau points
	vector_type tau_points_;
};




