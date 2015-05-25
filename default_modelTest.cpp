/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2010 by Sebastian Fuchs <fuchs@comp-phys.org>
*                       Thomas Pruschke <pruschke@comp-phys.org>
*                       Matthias Troyer <troyer@comp-phys.org>
*
* This software is part of the ALPS Applications, published under the ALPS
* Application License; you can use, redistribute it and/or modify it under
* the terms of the license, either version 1 or (at your option) any later
* version.
* 
* You should have received a copy of the ALPS Application License along with
* the ALPS Applications; see the file LICENSE.txt. If not, the license is also
* available from http://alps.comp-phys.org/.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
* FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
* SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
* FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
* DEALINGS IN THE SOFTWARE.
*
*****************************************************************************/


#include "default_model.hpp"
#include "alps/utility/temporary_filename.hpp"
#include"gtest/gtest.h"
#include <fstream>

void write_minimal_tab_file(const std::string &str){
  {
    std::ofstream tmpfile(str.c_str());
    tmpfile<<-20<<" "<<1./40<<std::endl;
    tmpfile<<20<<" "<<1./40<<std::endl;
    tmpfile.close();
  }
}
void write_tab_file_with_junk(const std::string &str){
  {
    std::ofstream tmpfile(str.c_str());
    tmpfile<<"# this is a file with junk in it."<<std::endl;
    tmpfile<<-20<<" "<<1./40<<" default model data" <<std::endl;
    tmpfile<<20<<" "<<1./40<<std::endl;
    tmpfile.close();
  }
}

///create a tabulated default model
TEST(TabFunction, TabFunctionConstruction){
  alps::params p; p["OMEGA_MAX"]=20; p["OMEGA_MIN"]=-20;
  std::string tabparamname="TABFILE";
  std::string tf=alps::temporary_filename("tab_file.dat");
  p[tabparamname]=tf;
  write_minimal_tab_file(tf);
  TabFunction T(p, tabparamname);
  boost::filesystem::remove(tf);
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
  boost::filesystem::remove(tf);
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
  boost::filesystem::remove(tf);

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
    boost::filesystem::remove(tf);
    EXPECT_EQ(err.what(),std::string("Input range outside of default model"));
  }
  catch(...){
    boost::filesystem::remove(tf);
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
