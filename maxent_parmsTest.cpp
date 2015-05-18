//write copywrite header when included with ALPSCore

#include "maxent.hpp"
#include "gtest/gtest.h"
#include <iostream>

TEST(Parameters,ContiParams){
	//set up parameters
	alps::params p;
	p["N_ALPHA"] = 60;
	p["ALPHA_MAX"] = .01;
	p["ALPHA_MAX"] = 20;
	p["NORM"] = 1.0;
	p["OMEGA_MAX"] = 6;
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
	
	//check random input values
	EXPECT_EQ(c.y(1),0.2);
	EXPECT_EQ(c.sigma(3),0.5);

	//test flat model is setup
	EXPECT_EQ(c.Default().omega(1.0),6);
}