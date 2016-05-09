/*
 * Copyright (C) 1998-2016 ALPS Collaboration. See COPYRIGHT.TXT
 * All rights reserved. Use is subject to license terms. See LICENSE.TXT
 * For use in publications, see ACKNOWLEDGE.TXT
 */

#include "../src/maxent.hpp"
#include "gtest.h"
#include <alps/utilities/temporary_filename.hpp>
#include <iostream>
#include "write_test_files.hpp"

TEST(Parameters,CatchMissingDataInParamFile){
  std::string pf=alps::temporary_filename("param_file.dat");
  write_minimal_param_file(pf);

  //fake input
  alps::params p(pf);
  MaxEntSimulation::define_parameters(p);
  p["NDAT"] = 5;

  try{
    ContiParameters c(p);
    boost::filesystem::remove(pf);
    FAIL() << "expected error on missing X_4 datapoint";

  }
  catch(...){
    //woo
    SUCCEED();
  }
  
  p["X_4"]=0.5;
  
  try{
    ContiParameters c2(p);
    boost::filesystem::remove(pf);
    FAIL() << "expected error on missing SIGMA_4";
  }
  catch(...){
    SUCCEED();
    boost::filesystem::remove(pf);
  }
}

TEST(Parameters,CatchMissingDataInDataFile){
  std::string pf=alps::temporary_filename("in_file.dat");
  write_minimal_input_file(pf);

  //fake input
  alps::params p;
  MaxEntSimulation::define_parameters(p);
  p["BETA"]=3;
  p["DATA"]=pf;
  //test for off by 1 error on data input
  p["NDAT"] = 6;

  try{
    ContiParameters c(p);
    boost::filesystem::remove(pf);
    FAIL() << "expected error on missing 4th datapoint";

  }
  catch(...){
    //error caught, user warned
    SUCCEED();
    boost::filesystem::remove(pf);
  }
}
