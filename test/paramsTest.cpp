/*
 * Copyright (C) 1998-2018 ALPS Collaboration.
 * All rights reserved. Use is subject to license terms. See LICENSE.TXT
 * For use in publications, see ACKNOWLEDGE.TXT
 */

#include "../src/maxent.hpp"
#include "gtest.h"
#include <alps/hdf5.hpp>
#include <alps/utilities/temporary_filename.hpp>
#include <alps/hdf5/vector.hpp>
#include <iostream>
#include "write_test_files.hpp"

TEST(Parameters,ContiParams){
	//set up parameters
	alps::params p;
	MaxEntSimulation::define_parameters(p);
	p["N_ALPHA"] = 60;
	p["ALPHA_MIN"] = .01;
	p["ALPHA_MAX"] = 20.0;
	p["NORM"] = 1.0;
	p["OMEGA_MAX"] = 6.0;
	p["KERNEL"] = "fermionic"; 
	p["BETA"] = 2.0;
	p["NFREQ"] = 2001;
	p["NDAT"] = 5;
	p["FREQUENCY_GRID"] = "Lorentzian";
	p["DATASPACE"] = "frequency";
	p["MAX_IT"] = 200;
	p["DEFAULT_MODEL"] = "flat";
	p["PARTICLE_HOLE_SYMMETRY"] = true;
	p["TEXT_OUTPUT"] = false;
	p["X_0"]=.1;
	p["X_1"]=.2;
	p["X_2"]=.3;
	p["X_3"]=.4;
	p["X_4"]=.5;
	p["SIGMA_0"]=0.5;
	p["SIGMA_1"]=0.5;
	p["SIGMA_2"]=0.5;
	p["SIGMA_3"]=0.5;
	p["SIGMA_4"]=0.5;
	ContiParameters c(p);

	EXPECT_EQ(c.ndat(),5);
	EXPECT_EQ(c.T(),0.5);
	
	//check input values
    for(int i=0;i<c.ndat();i++){
        EXPECT_NEAR(c.y(i),(i+1)*0.1,1e-10);
        EXPECT_EQ(c.sigma(i),0.5);
    }
}
TEST(Parameters,DataInParam){
  std::string pf=alps::temporary_filename("param_file.dat");
  write_minimal_param_file(pf);

  //fake input
  alps::params p(pf);
  MaxEntSimulation::define_parameters(p);
  p["NDAT"] = 4;
  
  ContiParameters c(p);
  EXPECT_EQ(c.ndat(),4);
  EXPECT_EQ(c.T(),0.5);

  for(int i=0;i<c.ndat();i++){
    EXPECT_NEAR(c.y(i),(i+1)*0.1,1e-10);
    EXPECT_EQ(c.sigma(i),0.5);
  }

  std::remove(pf.c_str());
}

TEST(Parameters,DataInFile){
std::string pf=alps::temporary_filename("in_file.dat");
  write_minimal_input_file(pf);

  //fake input
  alps::params p;
  MaxEntSimulation::define_parameters(p);
  p["BETA"]=2;
  p["DATA"]=pf;
  p["NDAT"] = 5;

  ContiParameters c(p);
  EXPECT_EQ(c.ndat(),5);
  EXPECT_EQ(c.T(),0.5);

  for(int i=0;i<c.ndat();i++){
    EXPECT_NEAR(c.y(i),(i+1)*0.1,1e-10);
    EXPECT_EQ(c.sigma(i),0.5);
  }

  std::remove(pf.c_str());
}

TEST(Parameters,MaxentParams){
    //set up parameters
	alps::params p;
	MaxEntSimulation::define_parameters(p);
	p["N_ALPHA"] = 60;
	p["ALPHA_MIN"] = .01;
	p["ALPHA_MAX"] = 20.0;
	p["NORM"] = 1.0;
	p["OMEGA_MAX"] = 6.0;
	p["KERNEL"] = "fermionic";
	p["BETA"] = 2.0;
	p["NFREQ"] = 2001;
	p["NDAT"] = 5;
	p["FREQUENCY_GRID"] = "Lorentzian";
	p["DATASPACE"] = "frequency";
	p["MAX_IT"] = 200;
	p["DEFAULT_MODEL"] = "flat";
	p["PARTICLE_HOLE_SYMMETRY"] = true;
	p["TEXT_OUTPUT"] = false;
	p["X_0"]=.1;
	p["X_1"]=.2;
	p["X_2"]=.3;
	p["X_3"]=.4;
	p["X_4"]=.5;
	p["SIGMA_0"]=0.5;
	p["SIGMA_1"]=0.5;
	p["SIGMA_2"]=0.5;
	p["SIGMA_3"]=0.5;
	p["SIGMA_4"]=0.5;
	MaxEntParameters c(p);
    
	EXPECT_EQ(c.ndat(),5);
	EXPECT_EQ(c.T(),0.5);
	
	//check input values; should be scaled by sigma
    for(int i=0;i<c.ndat();i++){
        EXPECT_NEAR(c.y(i),(i+1)*.1/0.5,1e-10);
        EXPECT_EQ(c.sigma(i),0.5);
    }
    
	//test flat model is setup
    EXPECT_NEAR(c.Default().omega(0),-6,1e-10);
	EXPECT_NEAR(c.Default().omega(1.0),6,1e-10);
}
TEST(Parameters,HighFrequencyCheck){
    alps::params p;
    MaxEntSimulation::define_parameters(p);
	p["N_ALPHA"] = 60;
	p["ALPHA_MIN"] = .01;
	p["ALPHA_MAX"] = 20.0;
	p["NORM"] = 1.0;
	p["OMEGA_MAX"] = 6.0;
	p["KERNEL"] = "fermionic";
	p["BETA"] = 2.0;
	p["NFREQ"] = 2001;
	p["NDAT"] = 5;
	p["FREQUENCY_GRID"] = "Lorentzian";
	p["DATASPACE"] = "frequency";
	p["PARTICLE_HOLE_SYMMETRY"] = true;
	p["TEXT_OUTPUT"] = false;
	p["X_0"]=.1;
	p["X_1"]=.2;
	p["X_2"]=.3;
	p["X_3"]=.4;
	p["X_4"]=.5;
	p["SIGMA_0"]=0.5;
	p["SIGMA_1"]=0.5;
	p["SIGMA_2"]=0.5;
	p["SIGMA_3"]=0.5;
	p["SIGMA_4"]=0.5;
    const int numModels = 5;
    std::string models[numModels] = {"flat","Gaussian","lorentzian","double gaussian","double lorentzian"};
    
    //specific params for default models
    p["SIGMA"]=1.0;
    p["GAMMA"]=.1;
    p["LAMBDA"]=.5;
    p["SHIFT"]=2.0;
    
    std::complex<double> G;
    for(int i=0;i<numModels;i++){
        p["DEFAULT_MODEL"] = models[i];
        MaxEntParameters c(p);
        
        // G(iw_{n})=\sum_{m}K_{nm}A_{m}
        G=0;
        std::complex<double> iwn(0,(2*2024+1)*M_PI*c.T());
        for(int j=0;j<c.nfreq();j++){
            G+= 1.0/(iwn-c.omega_coord(j))*c.Default().D(c.omega_coord(j)) * c.delta_omega(j);
        }
        std::cout<<c.delta_omega(0)<<" " << c.delta_omega(c.nfreq()-1) << std::endl;
        double limit = G.imag()*iwn.imag();

        EXPECT_TRUE(std::abs(1+limit)<.1);
    }
}

TEST(Parameters,HDF5ContiParams){
  //here we make a temporary HDF5 param file
  //and make sure Maxent does all the necessary reading

  alps::params p1; 
  std::string tf=alps::temporary_filename("hdf5_input.h5");
  alps::hdf5::archive oar(tf, "w");
  MaxEntSimulation::define_parameters(p1);
  p1["N_ALPHA"] = 60;
	p1["ALPHA_MIN"] = .01;
	p1["ALPHA_MAX"] = 20.0;
	p1["NORM"] = 1.0;
	p1["OMEGA_MAX"] = 6.0;
	p1["KERNEL"] = "fermionic"; 
	p1["BETA"] = 2.0;
	p1["NFREQ"] = 2001;
	p1["NDAT"] = 5;
	p1["FREQUENCY_GRID"] = "Lorentzian";
	p1["DATASPACE"] = "frequency";
	p1["MAX_IT"] = 200;
	p1["DEFAULT_MODEL"] = "flat";
	p1["PARTICLE_HOLE_SYMMETRY"] = true;
	p1["TEXT_OUTPUT"] = false;
  p1["DATA_IN_HDF5"] = true;
  p1["DATA"]=tf;
  oar["/parameters"] << p1;
  std::vector<double> x,sigma;
  for(int i=1;i<6;i++){
    x.push_back(i/10.0);
    sigma.push_back(0.5);
  }
  oar["/Data"] << x;
  oar["/Error"] << sigma; 
  oar.close();
  //now read in same file
  alps::hdf5::archive iar(tf, "r");
  alps::params p;
  iar["/parameters"] >> p;
  iar.close();
  
  ContiParameters c(p);

	EXPECT_EQ(c.ndat(),5);
	EXPECT_EQ(c.T(),0.5);
	
	//check input values
  for(int i=0;i<c.ndat();i++){
    EXPECT_NEAR(c.y(i),(i+1)*0.1,1e-10);
    EXPECT_EQ(c.sigma(i),0.5);
  }
  std::remove(tf.c_str());
}
TEST(Parameters,TZero){
	//set up parameters
	alps::params p;
	MaxEntSimulation::define_parameters(p);
	p["OMEGA_MAX"] = 6.0;
  p["OMEGA_MIN"]=0;
	p["KERNEL"] = "TZero";
  p["DATASPACE"] = "time"; 
	p["BETA"] = 9999999999; //should allow anything
	p["NDAT"] = 5;
	p["TEXT_OUTPUT"] = false;
	p["X_0"]=-0.5;
  p["X_1"]=-0.31968831176465;
  p["X_2"]=-0.22989567294666;
  p["X_3"]=-0.18028809394476;
  p["X_4"]=-0.15025032491734;
	p["SIGMA_0"]=0.5;
	p["SIGMA_1"]=0.5;
	p["SIGMA_2"]=0.5;
	p["SIGMA_3"]=0.5;
	p["SIGMA_4"]=0.5;
  p["TAU_0"]=0;
  p["TAU_1"]=0.3203125;
  p["TAU_2"]=0.640625;
  p["TAU_3"]=0.9609375;
  p["TAU_4"]=1.28125;

	MaxEntParameters c(p);
  //this is really just to make sure
  //nothing weird happens at T=0

	EXPECT_EQ(c.ndat(),5);
	
}

TEST(Parameters,CovarianceDataInFile){
  std::string pf=alps::temporary_filename("in_file.dat");
  write_minimal_input_file(pf);
  //create a diagonal covariance matrix as input
  //this should match the same tests of DataInFile
  std::string cov=alps::temporary_filename("cov.dat");
  std::ofstream tempfile(cov.c_str());
  tempfile<<0<< " "<<0<< " " << 0.3<<std::endl;
  tempfile<<0<< " "<<1<< " " << 0  <<std::endl;
  tempfile<<0<< " "<<2<< " " << 0  <<std::endl;
  tempfile<<0<< " "<<3<< " " << 0  <<std::endl;
  tempfile<<0<< " "<<4<< " " << 0  <<std::endl;
  tempfile<<1<< " "<<0<< " " << 0  <<std::endl;
  tempfile<<1<< " "<<1<< " " << 0.3<<std::endl;
  tempfile<<1<< " "<<2<< " " << 0  <<std::endl;
  tempfile<<1<< " "<<3<< " " << 0  <<std::endl;
  tempfile<<1<< " "<<4<< " " << 0  <<std::endl;
  tempfile<<2<< " "<<0<< " " << 0  <<std::endl;
  tempfile<<2<< " "<<1<< " " << 0  <<std::endl;
  tempfile<<2<< " "<<2<< " " << 0.3<<std::endl;
  tempfile<<2<< " "<<3<< " " << 0  <<std::endl;
  tempfile<<2<< " "<<4<< " " << 0  <<std::endl;
  tempfile<<3<< " "<<0<< " " << 0  <<std::endl;
  tempfile<<3<< " "<<1<< " " << 0  <<std::endl;
  tempfile<<3<< " "<<2<< " " << 0  <<std::endl;
  tempfile<<3<< " "<<3<< " " << 0.3<<std::endl;
  tempfile<<3<< " "<<4<< " " << 0  <<std::endl;
  tempfile<<4<< " "<<0<< " " << 0  <<std::endl;
  tempfile<<4<< " "<<1<< " " << 0  <<std::endl;
  tempfile<<4<< " "<<2<< " " << 0  <<std::endl;
  tempfile<<4<< " "<<3<< " " << 0  <<std::endl;
  tempfile<<4<< " "<<4<< " " << 0.3<<std::endl;
  tempfile.close();

  //fake input
  alps::params p;
  MaxEntSimulation::define_parameters(p);
  p["BETA"]=2;
  p["DATA"]=pf;
  p["NDAT"] = 5;
  p["COVARIANCE_MATRIX"]=cov;

  //MaxEntParameters handles covariance scaling
  MaxEntParameters c(p);
  EXPECT_EQ(c.ndat(),5);
  EXPECT_EQ(c.T(),0.5);

  for(int i=0;i<c.ndat();i++){
    EXPECT_NEAR(c.y(i),(i+1)*.1/0.3,1e-10);
    EXPECT_EQ(c.sigma(i),0.3);
  }

  std::remove(pf.c_str());
  std::remove(cov.c_str());
}
TEST(Parameters,CovarianceHDF5Params){
  //here we make a temporary HDF5 param file
  //and make sure Maxent does all the necessary reading

  alps::params p1; 
  std::string tf=alps::temporary_filename("hdf5_input.h5");
  alps::hdf5::archive oar(tf, "w");
  MaxEntSimulation::define_parameters(p1);
  p1["N_ALPHA"] = 60;
	p1["ALPHA_MIN"] = .01;
	p1["ALPHA_MAX"] = 20.0;
	p1["NORM"] = 1.0;
	p1["OMEGA_MAX"] = 6.0;
	p1["KERNEL"] = "fermionic"; 
	p1["BETA"] = 2.0;
	p1["NFREQ"] = 2001;
	p1["NDAT"] = 5;
	p1["FREQUENCY_GRID"] = "Lorentzian";
	p1["DATASPACE"] = "frequency";
	p1["MAX_IT"] = 200;
	p1["DEFAULT_MODEL"] = "flat";
	p1["PARTICLE_HOLE_SYMMETRY"] = true;
	p1["TEXT_OUTPUT"] = false;
  p1["DATA_IN_HDF5"] = true;
  p1["DATA"]=tf;
  //make it non-empty so it gets read in
  p1["COVARIANCE_MATRIX"]="Garbage";
  oar["/parameters"] << p1;
  std::vector<double> x,sigma;
  for(int i=1;i<6;i++){
    x.push_back(i/10.0);
    sigma.push_back(0.5);
  }
  Eigen::MatrixXd cov = Eigen::MatrixXd::Identity(5,5);
  cov *= 0.3; //change error from sigma
  std::vector<double> cov_vec;
  for(int i=0;i<5;i++){
    for(int j=0;j<5;j++)
      cov_vec.push_back(cov(i,j));
  }
  oar["/Data"] << x;
  oar["/Error"] << sigma; 
  oar["/Covariance"] << cov_vec;
  oar.close();
  //now read in same file
  alps::hdf5::archive iar(tf, "r");
  alps::params p;
  iar["/parameters"] >> p;
  iar.close();
  
  MaxEntParameters c(p);

	EXPECT_EQ(c.ndat(),5);
	EXPECT_EQ(c.T(),0.5);
	
	//check input values
  for(int i=0;i<c.ndat();i++){
    EXPECT_NEAR(c.y(i),(i+1)*.1/0.3,1e-10);
    EXPECT_EQ(c.sigma(i),0.3);
  }
  std::remove(tf.c_str());
}

