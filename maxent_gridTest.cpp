/*
 * Copyright (C) 1998-2016 ALPS Collaboration. See COPYRIGHT.TXT
 * All rights reserved. Use is subject to license terms. See LICENSE.TXT
 * For use in publications, see ACKNOWLEDGE.TXT
 */

#include "maxent.hpp"
#include "maxent_grid.hpp"
#include "gtest/gtest.h"
#include <fstream>

//these tests only make sure the grids are initiallized
//to the regions of [0,1] within the NFREQ+1 size array

TEST(Grid,LorentzianOdd){
    alps::params p;
    MaxEntSimulation::define_parameters(p);
    const int NFREQ=2001;
    p["FREQUENCY_GRID"] = "Lorentzian";
    p["NFREQ"] = NFREQ;
    
    grid g(p);
    EXPECT_EQ(g(0),0);
    EXPECT_EQ(g(NFREQ),1);
}
TEST(Grid,HalfLorentzianOdd){
    alps::params p;
    MaxEntSimulation::define_parameters(p);
    
    const int NFREQ=2001;
    p["FREQUENCY_GRID"] = "half lorentzian";
    p["NFREQ"] = NFREQ;
    
    grid g(p);
    EXPECT_EQ(g(0),0);
    EXPECT_EQ(g(NFREQ),1);
}
TEST(Grid,QuadraticOdd){
    alps::params p;
    MaxEntSimulation::define_parameters(p);

    const int NFREQ=2001;
    p["FREQUENCY_GRID"] = "Quadratic";
    p["NFREQ"] = NFREQ;
    
    grid g(p);
    EXPECT_EQ(g(0),0);
    EXPECT_EQ(g(NFREQ),1);
}
TEST(Grid,LogOdd){
    alps::params p;
    MaxEntSimulation::define_parameters(p);

    const int NFREQ=2001;
    p["FREQUENCY_GRID"] = "log";
    p["NFREQ"] = NFREQ;
    
    grid g(p);
    EXPECT_NEAR(g(0),0,1e-10);
    //log overshoots the default value by ~0.001
    EXPECT_NEAR(g(NFREQ),1,0.01);
}
TEST(Grid,LinearOdd){
    alps::params p;
    MaxEntSimulation::define_parameters(p);

    const int NFREQ=2001;
    p["FREQUENCY_GRID"] = "linear";
    p["NFREQ"] = NFREQ;
    
    grid g(p);
    EXPECT_EQ(g(0),0);
    EXPECT_EQ(g(NFREQ),1);
}
//the following tests have an even amount of points
//so we test the same endpoints (which shouldn't change)
//as well as the middle of symmetric grids
TEST(Grid,LorentzianEven){
    alps::params p;
    MaxEntSimulation::define_parameters(p);

    const int NFREQ=2000;
    p["FREQUENCY_GRID"] = "Lorentzian";
    p["NFREQ"] = NFREQ;
    
    grid g(p);
    EXPECT_EQ(g(0),0);
    EXPECT_EQ(g(NFREQ/2), 0.5);
    EXPECT_EQ(g(NFREQ),1);
}
TEST(Grid,HalfLorentzianEven){
    //this is not symmetric so there is no new test
    alps::params p;
    MaxEntSimulation::define_parameters(p);

    const int NFREQ=2000;
    p["FREQUENCY_GRID"] = "half lorentzian";
    p["NFREQ"] = NFREQ;
    
    grid g(p);
    EXPECT_EQ(g(0),0);
    EXPECT_EQ(g(NFREQ),1);
}
TEST(Grid,QuadraticEven){
    alps::params p;
    MaxEntSimulation::define_parameters(p);

    const int NFREQ=2000;
    p["FREQUENCY_GRID"] = "Quadratic";
    p["NFREQ"] = NFREQ;
    
    grid g(p);
    EXPECT_EQ(g(0),0);
    EXPECT_NEAR(g(NFREQ/2), 0.5,1e-10);
    EXPECT_EQ(g(NFREQ),1);
}
TEST(Grid,LogEven){
    alps::params p;
    MaxEntSimulation::define_parameters(p);
    
    const int NFREQ=2000;
    p["FREQUENCY_GRID"] = "log";
    p["NFREQ"] = NFREQ;
    
    grid g(p);
    EXPECT_NEAR(g(0),0,1e-10);
    EXPECT_EQ(g(NFREQ/2), 0.5);
    //log overshoots the default value by ~0.001
    EXPECT_NEAR(g(NFREQ),1,0.01);
}
TEST(Grid,LinearEven){
    alps::params p;
    MaxEntSimulation::define_parameters(p);
    const int NFREQ=2000;
    p["FREQUENCY_GRID"] = "linear";
    p["NFREQ"] = NFREQ;
    
    grid g(p);
    EXPECT_EQ(g(0),0);
    EXPECT_EQ(g(NFREQ/2), 0.5);
    EXPECT_EQ(g(NFREQ),1);
}



