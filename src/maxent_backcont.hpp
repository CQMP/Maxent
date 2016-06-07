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
#include "maxent.hpp"
#include <alps/config.hpp> // needed to set up correct bindings

///This class contains necessities to analytically continue a 
//real function to the imaginary axis
class Backcont{
  public:
    Backcont(const MaxEntParameters *param_in);
    ///Backcontinue a given function A
    vector_type backcontinue(const vector_type &A);
    ///determine the maximum difference between the two functions
    double max_error(const vector_type &y1, const vector_type &y2);
  private:
   //pointer to Parameters, set up through MaxEntSimulation->MaxEntHelper
   const MaxEntParameters *param;
   ///holds backcontinued Green's function
   vector_type G;
   kernel_type k_type;
};
