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
#include<iomanip>
#include<alps/params.hpp>
#include "maxent_matrix_def.hpp"

///this class is responsible for creating real frequency grids. TODO: replace some of this functionality with the ALPSCore grids and eliminate it.
class grid{
public:
  grid(const alps::params &p);
  ///define parameter defaults
  static void define_parameters(alps::params &p);
  ///output a help function for the parameters that initialize this grid.
  static void print_help();

  const std::vector<double> &t_array() const{return t_array_;}
  double operator()(int i)const{return t_array_[i];}
private:
  ///the number of (real) frequency point for this grid.
  int nfreq_;
  std::vector<double> t_array_;

  void initialize_linear_grid();
  void initialize_logarithmic_grid(double t_min);
  void initialize_quadratic_grid(double spread);
  void initialize_half_lorentzian_grid(double cut);
  void initialize_lorentzian_grid(double cut);
};




