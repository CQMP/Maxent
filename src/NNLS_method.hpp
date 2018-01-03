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
private:
  ///overall normalization
  const double norm;
  ///overall simulation base name
  std::string name;
  ///String identifying the kernel type, have a look at --help for available kernels
  std::string Kernel_type;
  ///number of real frequency points
  const int nfreq;
};
