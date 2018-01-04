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

#include <fstream>
#include "maxent_matrix_def.hpp"
#include "maxent_parms.hpp"
#include <boost/random/mersenne_twister.hpp>

class NNLS_Simulation:protected KernelAndGridIO{
public:
  ///setup of parameters
  NNLS_Simulation(alps::params& parms);
  ///setup parameters
  static void define_parameters(alps::params &parms);
  ///the fitting using ADMM calculation
  void run();
  ///the evaluation and writing of files
  void evaluate();
  ///compute the first derivative in a finite difference approximation
  void compute_first_derivative_matrix();
  ///compute the second derivative in a finite difference approximation
  void compute_second_derivative_matrix();

  ///const getter functions
  const vector_type &omega_coord() const{return omega_coord_;}
  const vector_type &delta_omega() const{return delta_omega_;}
  const vector_type &inputGrid() const {return inputGrid_;}
  const matrix_type &K() const {return K_;}
  const matrix_type &L1() const {return L1_;}
  const matrix_type &L2() const {return L2_;}

private:
  ///overall normalization
  const double norm;
  ///overall simulation base name
  std::string name;
  ///String identifying the kernel type, have a look at --help for available kernels
  std::string Kernel_type;
  ///number of real frequency points
  const int nfreq_;

  ///the real frequency grid
  vector_type omega_coord_;
  ///differences of the real frequency grid
  vector_type delta_omega_;
  ///finite difference matrix of first derivative
  matrix_type L1_;
  ///finite difference matrix of second derivative
  matrix_type L2_;
};
