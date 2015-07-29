/*
 * Copyright (C) 1998-2015 ALPS Collaboration. See COPYRIGHT.TXT
 * All rights reserved. Use is subject to license terms. See LICENSE.TXT
 * For use in publications, see ACKNOWLEDGE.TXT
 */

#include "pade.hpp"

int main(int argc, const char* argv[])
{
  //parse parameters etc
  PadeParams parms(argc, argv);
  
  //read in and set up real and imaginary domain
  imaginary_domain_data f_iomega(parms);
  real_domain_data f_omega(parms);

  f_iomega.write("input_data.dat");
  
  //set up pade class
  pade_interpolator P(parms);
  P.pade_interpolate(f_iomega, f_omega);
  
  f_omega.write(parms["real.OUTPUT"]);
  return 0;
}
