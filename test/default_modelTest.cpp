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


#include "../src/default_model.hpp"
#include <alps/utilities/temporary_filename.hpp>
#include"gtest.h"
#include <fstream>
#include "write_test_files.hpp"


///create a tabulated default model
TEST(TabFunction, TabFunctionConstruction){
  alps::params p; p["OMEGA_MAX"]=20; p["OMEGA_MIN"]=-20;
  std::string tabparamname="TABFILE";
  std::string tf=alps::temporary_filename("tab_file.dat");
  p[tabparamname]=tf;
  write_minimal_tab_file(tf);
  TabFunction T(p, tabparamname);
  std::remove(tf.c_str());
  EXPECT_EQ(T(-20.),1./40);
  EXPECT_EQ(T(20.),1./40);
  EXPECT_EQ(T(0.),1./40);
}
TEST(TabFunction, TabFunctionConstructionWithMoreLines){
  alps::params p; p["OMEGA_MAX"]=20; p["OMEGA_MIN"]=-20;
  std::string tabparamname="TABFILE";
  std::string tf=alps::temporary_filename("tab_file.dat");
  p[tabparamname]=tf;
  write_minimal_tab_file(tf);
  TabFunction T(p, tabparamname);
  std::remove(tf.c_str());
  EXPECT_NEAR(T(-20.),1./40, 1.e-12);
  EXPECT_NEAR(T(20.),1./40, 1.e-12);
  EXPECT_NEAR(T(0.),1./40, 1.e-12);
}
TEST(TabFunction, OutsideMinMaxIsZero){
  alps::params p; p["OMEGA_MAX"]=20; p["OMEGA_MIN"]=-20;
  std::string tabparamname="TABFILE";
  std::string tf=alps::temporary_filename("tab_file.dat");
  p[tabparamname]=tf;
  write_minimal_tab_file(tf);
  TabFunction T(p, tabparamname);
  std::remove(tf.c_str());

  EXPECT_EQ(T(-21.),0.);
  EXPECT_EQ(T(21.),0.);

}
TEST(TabFunction, FailOutsideRange){
  //tab file is [-20,20]
  //ask for range [-20,25] so it will fail
  alps::params p; p["OMEGA_MAX"]=25; p["OMEGA_MIN"]=-20;
  std::string tabparamname="TABFILE";
  std::string tf=alps::temporary_filename("tab_file.dat");
  p[tabparamname]=tf;
  write_minimal_tab_file(tf);
  try{
    TabFunction T(p, tabparamname);
    FAIL() << "Expected to fail parameters out of range";
  }
  catch(std::logic_error const & err){
    std::remove(tf.c_str());
    EXPECT_EQ(err.what(),std::string("Input range outside of default model"));
  }
  catch(...){
    std::remove(tf.c_str());
    FAIL() << "expected parameters out of range error";
  }
}
TEST(Gaussian, TwoGaussiansIsGaussian){
  alps::params p;
  p["SIGMA1"]=1;
  p["SIGMA2"]=1;
  p["SHIFT1"]=0;
  p["SHIFT2"]=-1;
  p["NORM1"]=1;

  TwoGaussians TG(p);
  alps::params q;
  q["SIGMA"]=1;

  Gaussian G(q);

  for(int i=-5;i<5;++i){
    EXPECT_NEAR(TG(i),G(i),1.e-10);
  }
}
TEST(Gaussian, TwoGaussiansIsDoubleGaussian){
  alps::params p;
  p["SIGMA1"]=1;
  p["SIGMA2"]=1;
  p["SHIFT1"]=1;
  p["SHIFT2"]=-1;
  p["NORM1"]=0.5;

  TwoGaussians TG(p);
  alps::params q;
  q["SIGMA"]=1;
  q["SHIFT"]=1;

  DoubleGaussian D(q);
  for(int i=-5;i<5;++i){
    EXPECT_NEAR(TG(i),D(i),1.e-10);
  }
}
