/*
 * Copyright (C) 1998-2015 ALPS Collaboration. See COPYRIGHT.TXT
 * All rights reserved. Use is subject to license terms. See LICENSE.TXT
 * For use in publications, see ACKNOWLEDGE.TXT
 */

#include "pade.hpp"
#include <fstream>

grid::grid(const PadeParams &p){
  N_freq_=p["real.NFREQ"];
  std::cout<<"using a grid size of : "<<N_freq_<<" real frequency points."<<std::endl;
  t_array_.resize(N_freq_+1);
  freq_.resize(N_freq_);
  delta_freq_.resize(N_freq_);
  
  //set the grid type
  std::string p_f_grid = p["real.FREQUENCY_GRID"];
  if(p_f_grid=="Lorentzian") grid_type_=lorentzian_grid;
  else if (p_f_grid=="half Lorentzian") grid_type_=half_lorentzian_grid;
  else if (p_f_grid=="quadratic") grid_type_=quadratic_grid;
  else if (p_f_grid=="log") grid_type_=log_grid;
  else if (p_f_grid=="linear") grid_type_=linear_grid;
  else throw std::invalid_argument("grid type parameter FREQUENCY_GRID should be one of: Lorentzian, half Lorentzian, quadratic, log, linear");
  setup_grid(p);
}
void grid::setup_grid(const PadeParams &p){
  std::vector<double> temp_(N_freq_+1);
  if(grid_type_==lorentzian_grid){
    double cut = p["real.CUT"];
    std::cout<<"using Lorentzian grid with CUT parameter: "<<cut<<std::endl;
    for (int i=0; i<N_freq_+1; ++i)
      temp_[i] = tan(M_PI * (double(i)/(N_freq_)*(1.-2*cut)+cut - 0.5));
    for (int i=0; i<N_freq_+1; ++i)
      t_array_[i] = (temp_[i] - temp_[0])/(temp_[N_freq_] - temp_[0]);
  }
  else if (grid_type_==half_lorentzian_grid){
    double cut = p["real.CUT"];
    std::cout<<"using half Lorentzian grid with CUT parameter: "<<cut<<std::endl;
    for (int i=0; i<N_freq_; ++i)
      temp_[i] = tan(M_PI * (double(i+N_freq_)/(2*N_freq_-1)*(1.-2*cut)+cut - 0.5));
    for (int i=0; i<N_freq_+1; ++i)
      t_array_[i] = (temp_[i] - temp_[0])/(temp_[N_freq_] - temp_[0]);\
  }
  else if (grid_type_==quadratic_grid) {
    double s = p["real.SPREAD"];
    if (s<1)
      boost::throw_exception(std::invalid_argument("the parameter SPREAD must be greater than 1"));
    std::cout<<"using quadratic grid with SPREAD parameter: "<<s<<std::endl;
    double t=0;
    for (int i=0; i<N_freq_-1; ++i) {
      double a = double(i)/(N_freq_-1);
      double factor = 4*(s-1)*(a*a-a)+s;
      factor /= double(N_freq_-1)/(3.*(N_freq_-2))*((N_freq_-1)*(2+s)-4+s);
      double delta_t = factor;
      t += delta_t;
      temp_[i] = t;
    }
    t_array_[0] = 0.;
    for (int i=1; i<N_freq_; ++i)
      t_array_[i]  = temp_[i-1]/temp_[N_freq_];
  }
  else if (grid_type_==log_grid) {
    double t_min = p["real.LOG_MIN"];
    double t_max=0.5;
    std::cout<<"using log grid with LOG_MIN parameter: "<<p["real.LOG_MIN"]<<std::endl;
    double scale=std::log(t_max/t_min)/((float) (N_freq_/2-1));
    t_array_[N_freq_/2] = 0.5;
    for (int i=0; i<N_freq_/2; ++i) {
      t_array_[N_freq_/2+i+1]  = 0.5+t_min*std::exp(((float) i)*scale);
      t_array_[N_freq_/2-i-1]  = 0.5-t_min*std::exp(((float) i)*scale);
    }
  }
  else if (grid_type_==linear_grid) {
    std::cout<<"using equidistantly spaced linear grid."<<std::endl;
    for (int i=0; i<N_freq_; ++i){
      t_array_[i] = double(i)/(N_freq_-1);
    }
  }
  else
    boost::throw_exception(std::invalid_argument("No valid frequency grid specified"));
  
  omega_max_=p["real.OMEGA_MAX"];
  omega_min_=p["real.OMEGA_MIN"];
  
  //store the frequencies and frequency differences (integration weights)
  for(int i=0;i<freq_.size();++i){
    freq_[i] = (omega_of_t(t_array_[i]) + omega_of_t(t_array_[i+1]))/2.;
    delta_freq_[i] = omega_of_t(t_array_[i+1]) - omega_of_t(t_array_[i]);
  }
}


imag_domain_grid::imag_domain_grid(const PadeParams &p){
  N_freq_= p["imag.NDAT"];
  freq_.resize(N_freq_);
  T_=1./p["imag.BETA"].as<double>();
  std::cout<<"statistics is: "<<p["imag.STATISTICS"]<<std::endl;
  if(p["imag.STATISTICS"]==std::string("Fermi")){
    for(int i=0;i<N_freq_;++i){
      freq_[i]=(2*i+1)*M_PI*T_;
    }
  }else if(p["imag.STATISTICS"]==std::string("Bose")){
    if(p["imag.NEGATIVE_DATA"]){
      if(N_freq_%2 ==0) throw std::invalid_argument("make sure your grid has ODD entries when doing bosons with neg freq");
      for(int i=0;i<N_freq_;++i){
        freq_[i]=(2*(i-N_freq_/2))*M_PI*T_;
        std::cout<<"freq: "<<freq_[i]<<std::endl;
      }
    }else{
      for(int i=0;i<N_freq_;++i){
        freq_[i]=(2*i)*M_PI*T_;
      }
    }
  }else
   throw std::runtime_error("statistics (Fermi/Bose) not understood");
}
