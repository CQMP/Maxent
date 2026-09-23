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
