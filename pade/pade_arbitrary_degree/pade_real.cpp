/*
 * Copyright (C) 1998-2016 ALPS Collaboration. See COPYRIGHT.TXT
 * All rights reserved. Use is subject to license terms. See LICENSE.TXT
 * For use in publications, see ACKNOWLEDGE.TXT
 */

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
