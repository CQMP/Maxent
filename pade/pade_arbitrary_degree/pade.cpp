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
