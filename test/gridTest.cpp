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
#include "../src/maxent_grid.hpp"
#include "gtest.h"
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



