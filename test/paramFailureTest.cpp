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
    KernelAndGrid c(p);
    boost::filesystem::remove(pf);
    FAIL() << "expected error on missing X_4 datapoint";

  }
  catch(...){
    //woo
    SUCCEED();
  }
  
  p["X_4"]=0.5;
  
  try{
    KernelAndGrid c2(p);
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
    KernelAndGrid c(p);
    boost::filesystem::remove(pf);
    FAIL() << "expected error on missing 4th datapoint";

  }
  catch(...){
    //error caught, user warned
    SUCCEED();
    boost::filesystem::remove(pf);
  }
}
