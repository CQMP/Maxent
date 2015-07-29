/*
 * Copyright (C) 1998-2015 ALPS Collaboration. See COPYRIGHT.TXT
 * All rights reserved. Use is subject to license terms. See LICENSE.TXT
 * For use in publications, see ACKNOWLEDGE.TXT
 */

#include "maxent.hpp"
#include "gtest/gtest.h"
#include <iostream>
#include "maxent_parms_default.hpp"

TEST(Parameters,ContiParams){
	//set up parameters
	alps::params p;
	set_defaults(p);
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
	p["PARTICLE_HOLE_SYMMETRY"] = 1;
	p["TEXT_OUTPUT"] = 0;
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
TEST(Parameters,MaxentParams){
    //set up parameters
	alps::params p;
	set_defaults(p);
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
	p["PARTICLE_HOLE_SYMMETRY"] = 1;
	p["TEXT_OUTPUT"] = 0;
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
TEST(Paramaters,HighFrequencyCheck){
    alps::params p;
    set_defaults(p);
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
	p["PARTICLE_HOLE_SYMMETRY"] = 1;
	p["TEXT_OUTPUT"] = 0;
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
