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
  private:
   //pointer to Parameters, set up through MaxEntSimulation->MaxEntHelper
   const MaxEntParameters *param;
   ///holds backcontinued Green's function
   vector_type G;
   kernel_type k_type;
};
