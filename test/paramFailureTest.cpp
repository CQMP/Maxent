/*
 * Copyright (C) 1998-2016 ALPS Collaboration. See COPYRIGHT.TXT
 * All rights reserved. Use is subject to license terms. See LICENSE.TXT
 * For use in publications, see ACKNOWLEDGE.TXT
 */

#include "../src/maxent.hpp"
#include "gtest/gtest.h"
#include <alps/utilities/temporary_filename.hpp>
#include <iostream>

void write_minimal_param_file(const std::string &str){
    std::ofstream tmpfile(str.c_str());
    tmpfile<<"BETA=2" <<std::endl;
    tmpfile<<"X_0=.1" <<std::endl;
    tmpfile<<"X_1=.2" <<std::endl;
    tmpfile<<"X_2=.3" <<std::endl;
    tmpfile<<"X_3=.4" <<std::endl;
    tmpfile<<"SIGMA_0=.5" <<std::endl;
    tmpfile<<"SIGMA_1=.5" <<std::endl;
    tmpfile<<"SIGMA_2=.5" <<std::endl;
    tmpfile<<"SIGMA_3=.5" <<std::endl;
    tmpfile.close();
}

TEST(Parameters,CatchMissingData){
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
