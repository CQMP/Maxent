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
#include<vector>
#include<alps/ngs/params.hpp>
#include"maxent_blas.hpp"

enum dataspace_type{
  time_dataspace,
  frequency_dataspace
};
enum kernel_type{
  time_fermionic_kernel,
  time_bosonic_kernel,
  time_boris_kernel,
  frequency_fermionic_ph_kernel,
  frequency_bosonic_ph_kernel,
  frequency_anomalous_ph_kernel,
  frequency_fermionic_kernel,
  frequency_bosonic_kernel,
  frequency_anomalous_kernel
};
class kernel{
public:
  kernel(const alps::params &p, const vector_type& freq);

  ///getter function for the kernel matrix
  const matrix_type &operator()()const{return K_;}
private:
  ///figure out which kernel is to be used
  void set_kernel_type(const std::string &dataspace_name, const std::string &kernel_name, bool ph_symmetry);
  ///number of Matsubara points
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
};




