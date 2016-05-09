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

  boost::filesystem::remove(pf);
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

  boost::filesystem::remove(pf);
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
  p1.save(oar);
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
  p.load(iar);
  iar.close();
  
  ContiParameters c(p);

	EXPECT_EQ(c.ndat(),5);
	EXPECT_EQ(c.T(),0.5);
	
	//check input values
  for(int i=0;i<c.ndat();i++){
    EXPECT_NEAR(c.y(i),(i+1)*0.1,1e-10);
    EXPECT_EQ(c.sigma(i),0.5);
  }
  boost::filesystem::remove(tf);
}
