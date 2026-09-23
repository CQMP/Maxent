/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1998-2016 ALPS Collaboration
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include "pade.hpp"
#include <fstream>


real_domain_data::real_domain_data(const PadeParams &p):G_(p){
  N_real_=p["real.NFREQ"];
  val_.resize(N_real_);
}
void real_domain_data::write(const std::string &s) const{
  std::ofstream file(s.c_str());
  for(int i=0;i<N_real_;++i){
    file<<G_.freq()[i]<<" "<<val_[i].real()<<" "<<val_[i].imag()<<std::endl;
  }
}
