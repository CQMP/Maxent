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
#include<alps/hdf5/vector.hpp>
#include "maxent_matrix_def.hpp"

///this class is responsible for creating mapping data for real frequency grids.
///TODO: replace some of this functionality with the ALPSCore grids and eliminate it.
///it allocates an array with nfreq_+1 integers. Those go from zero to 1, including zero and 1,
///and are evenly (in the case of uniform grids) or unevenly (in the case of general grids)
///spaced. OMEGA_MIN will be mapped to zero, OMEGA_MAX will be mapped to one.
///Note that in the actual code, when we evaluate points we use some sort of trapezoidal rule that
///evaluates nfreq_ grid points where point n of, say, a spectral function is located between point
///n and n+1 of this grid.
class grid{
public:
  grid(const alps::params &p);
  ///define parameter defaults
  static void define_parameters(alps::params &p);
  ///output a help function for the parameters that initialize this grid.
  static void print_help();
  ///raw data for the array mapping the interval [0,1] to a non-uniform grid.
  const std::vector<double> &t_array() const{return t_array_;}
  ///get an integer between 0 and (including) nfreq_. Get back a double between 0 and one corresponding to the grid.
  double operator()(int i)const{return t_array_[i];}

  ///equidistant mapping from [0,1] to [omega_min, omega_max]
  double omega_of_t(const double t) const { return omega_min_ + (omega_max_-omega_min_)*t; }
  ///equidistant mapping from [omega_min, omega_max] to [0,1]
  double t_of_omega(const double omega) const { return (omega-omega_min_)/(omega_max_-omega_min_); }
  ///get the minimum value of the grid
  double omega_min() const{return omega_min_;}
  ///get the maximum value fo the grid
  double omega_max() const{return omega_max_;}

private:
  ///the number of (real) frequency point for this grid.
  int nfreq_;
  std::vector<double> t_array_;
  ///maximum value of grid
  double omega_max_;
  ///minimum value fo grid
  double omega_min_;

  ///this is for an equidistantly spaced grid
  void initialize_linear_map();
  ///this is for a logarithmic grid
  void initialize_logarithmic_map(double t_min);
  ///this has a quadratic distribution
  void initialize_quadratic_map(double spread);
  ///the NRG-type logarithmic grid but only on the positive half axis
  void initialize_half_lorentzian_grid(double cut);
  ///this maps a lorentzian grid, i.e. has a very good resolution near zero and a bad hifreq resolution. Good choice for maxent.
  void initialize_lorentzian_grid(double cut);
};




