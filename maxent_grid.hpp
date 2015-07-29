/*
 * Copyright (C) 1998-2015 ALPS Collaboration. See COPYRIGHT.TXT
 * All rights reserved. Use is subject to license terms. See LICENSE.TXT
 * For use in publications, see ACKNOWLEDGE.TXT
 */
#pragma once
#include<vector>
#include<alps/params.hpp>
#include "maxent_matrix_def.hpp"

class grid{
public:
  grid(const alps::params &p);
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




