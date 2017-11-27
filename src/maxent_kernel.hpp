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
#include<vector>
#include<alps/params.hpp>
#include<iostream>
#include"maxent_matrix_def.hpp"
///enum that enumerates if we're in time, frequency, or legendre space
enum dataspace_type{
  time_dataspace,
  frequency_dataspace,
  legendre_dataspace
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
  frequency_anomalous_kernel
};
class kernel{
public:
  kernel(alps::params &p, const vector_type& freq, vector_type &inputGrid);

  ///getter function for the kernel matrix
  const matrix_type &operator()()const{return K_;}
  ///getter function for kernel type
  kernel_type getKernelType() const{return ktype_;}
  ///getter function for dataspace type
  dataspace_type getDataspaceType() const{return dtype_;}
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




