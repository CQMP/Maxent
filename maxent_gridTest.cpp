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


#include "maxent_grid.hpp"
#include"gtest/gtest.h"
#include <fstream>

//these tests only make sure the grids are initiallized
//to the regions of [0,1] within the NFREQ+1 size array

TEST(Grid,Lorentzian){
    alps::params p;
    const int NFREQ=2001;
    p["FREQUENCY_GRID"] = "Lorentzian";
    p["NFREQ"] = NFREQ;
    
    grid g(p);
    EXPECT_EQ(g(0),0);
    EXPECT_EQ(g(NFREQ),1);
}
TEST(Grid,HalfLorentzian){
    alps::params p;
    const int NFREQ=2001;
    p["FREQUENCY_GRID"] = "half lorentzian";
    p["NFREQ"] = NFREQ;
    
    grid g(p);
    EXPECT_EQ(g(0),0);
    EXPECT_EQ(g(NFREQ),1);
}
TEST(Grid,Quadratic){
    alps::params p;
    const int NFREQ=2001;
    p["FREQUENCY_GRID"] = "Quadratic";
    p["NFREQ"] = NFREQ;
    
    grid g(p);
    EXPECT_EQ(g(0),0);
    EXPECT_EQ(g(NFREQ),1);
}
TEST(Grid,Log){
    alps::params p;
    const int NFREQ=2001;
    p["FREQUENCY_GRID"] = "log";
    p["NFREQ"] = NFREQ;
    
    grid g(p);
    EXPECT_NEAR(g(0),0,1e-10);
    //log overshoots the default value by ~0.001
    EXPECT_NEAR(g(NFREQ),1,0.01);
}
TEST(Grid,Linear){
    alps::params p;
    const int NFREQ=2001;
    p["FREQUENCY_GRID"] = "linear";
    p["NFREQ"] = NFREQ;
    
    grid g(p);
    EXPECT_EQ(g(0),0);
    EXPECT_EQ(g(NFREQ),1);
}



