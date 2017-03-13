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
#include "maxent_grid.hpp"
#include<cmath>
#include <boost/algorithm/string.hpp>    

grid::grid(const alps::params &p):
nfreq_(p["NFREQ"]),
t_array_(nfreq_+1){
  std::string p_f_grid = p["FREQUENCY_GRID"];
  boost::to_lower(p_f_grid);
   double cut = p["CUT"];
  if (p_f_grid =="lorentzian") {
    initialize_lorentzian_grid(cut);
  }
  else if (p_f_grid=="half lorentzian") {
    initialize_half_lorentzian_grid(cut);
  }
  else if (p_f_grid=="quadratic") {
    double s = p["SPREAD"];
    initialize_quadratic_grid(s);
  }
  else if (p_f_grid=="log") {
    double t_min = p["LOG_MIN"];
    initialize_logarithmic_grid(t_min);
  }
  else if (p_f_grid=="linear") {
    initialize_linear_grid();
  }
  else
    boost::throw_exception(std::invalid_argument("No valid frequency grid specified"));
}
///define parameter defaults
void grid::define_parameters(alps::params &p){
  //---------------------------------
  //    Grid
  //---------------------------------
  p.define<double>("CUT",0.01,"cut for lorentzian grids");
  p.define<double>("SPREAD",4,"spread for quadratic grid");
  p.define<double>("LOG_MIN",1.0e-4,"log_min for log grid");
  p.define<std::string>("FREQUENCY_GRID","Lorentzian","Type of frequency grid");
  p.define<int>("NFREQ",1000,"Number of A(omega) real frequencies");
}

void grid::initialize_linear_grid() {
  for (int i = 0; i < nfreq_+1; ++i)
    t_array_[i] = double(i) / (nfreq_);
}

void grid::initialize_logarithmic_grid(double t_min) {
  double  t_max = 0.5;
  double scale = std::log(t_max / t_min) / ((float) ((nfreq_ / 2 - 1)));
  t_array_[nfreq_ / 2] = 0.5;
  for (int i = 0; i < nfreq_ / 2; ++i) {
    t_array_[nfreq_ / 2 + i + 1] = 0.5
        + t_min * std::exp(((float) (i)) * scale);
    t_array_[nfreq_ / 2 - i - 1] = 0.5
        - t_min * std::exp(((float) (i)) * scale);
  }
  //if we have an odd # of frequencies, this catches the last element
  if (nfreq_ % 2 != 0)
    t_array_[nfreq_ / 2 + nfreq_ / 2 + 1] = 0.5
    + t_min * std::exp(((float) (nfreq_) / 2) * scale);
}

void grid::initialize_quadratic_grid(double spread) {
  if (spread < 1)
    boost::throw_exception(
        std::invalid_argument("the parameter SPREAD must be greater than 1"));

  std::vector<double> temp(nfreq_);
  double t = 0;
  for (int i = 0; i < nfreq_; ++i) {
    double a = double(i) / (nfreq_ - 1);
    double factor = 4 * (spread - 1) * (a * a - a) + spread;
    factor /= double(nfreq_ - 1) / (3. * (nfreq_ - 2))
                  * ((nfreq_ - 1) * (2 + spread) - 4 + spread);
    double delta_t = factor;
    t += delta_t;
    temp[i] = t;
  }
  t_array_[0] = 0.;
  for (int i = 1; i < nfreq_+1; ++i)
    t_array_[i] = temp[i - 1] / temp[temp.size() - 1];
}

void grid::initialize_half_lorentzian_grid(double cut) {
  std::vector<double> temp(nfreq_ + 1);
  for (int i = 0; i < nfreq_ + 1; ++i)
    temp[i] = tan(
        M_PI
        * (double(i + nfreq_) / (2 * nfreq_ - 1) * (1. - 2 * cut) + cut
            - 0.5));
  for (int i = 0; i < nfreq_ + 1; ++i)
    t_array_[i] = (temp[i] - temp[0]) / (temp[temp.size() - 1] - temp[0]);
}

void grid::initialize_lorentzian_grid(double cut) {
  std::vector<double> temp(nfreq_ + 1);
  for (int i = 0; i < nfreq_ + 1; ++i)
    temp[i] = tan(M_PI * (double(i) / (nfreq_) * (1. - 2 * cut) + cut - 0.5));
  for (int i = 0; i < nfreq_ + 1; ++i)
    t_array_[i] = (temp[i] - temp[0]) / (temp[temp.size() - 1] - temp[0]);
}




