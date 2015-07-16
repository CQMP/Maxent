#include "maxent.hpp"
#include "gtest/gtest.h"
#include <iostream>
#include "maxent_parms_default.hpp"
double getNorm(vector_type &omega, vector_type &y){
    int size = omega.size();
    double norm = 0.0;
    for(int i=0;i<size-1;i++){
	norm+=(omega[i+1]-omega[i])*y[i];
    }
    norm+=(omega[size-1]-omega[size-2])*y[size-1];
    return norm;
}
TEST(Simulation,FrequencySimulation){
    alps::params p;
    set_defaults(p);
    p["BETA"]=8;
    p["NDAT"]=15;
    p["PARTICLE_HOLE_SYMMETRY"]=1;
    p["DATASPACE"]="frequency";
    p["TEXT_OUTPUT"]=0;


    //data from DMFT;U=0;beta=8
    p["X_0"]=-0.58900239090596;
    p["X_1"]=-0.40986302909581;
    p["X_2"]=-0.32440959736089;
    p["X_3"]=-0.26834511172698;
    p["X_4"]=-0.2278124578356;
    p["X_5"]=-0.19711608409901;
    p["X_6"]=-0.173170382827;
    p["X_7"]=-0.15406437082322;
    p["X_8"]=-0.13853122078645;
    p["X_9"]=-0.12569699334634;
    p["X_10"]=-0.1149417307568;
    p["X_11"]=-0.10581568778351;
    p["X_12"]=-0.097986218747872;
    p["X_13"]=-0.09120298645574;
    p["X_14"]=-0.085274608164208;
    p["X_15"]=-0.080052654686189;

    p["SIGMA_0"]=1e-6;
    p["SIGMA_1"]=1e-6;
    p["SIGMA_2"]=1e-6;
    p["SIGMA_3"]=1e-6;
    p["SIGMA_4"]=1e-6;
    p["SIGMA_5"]=1e-6;
    p["SIGMA_6"]=1e-6;
    p["SIGMA_7"]=1e-6;
    p["SIGMA_8"]=1e-6;
    p["SIGMA_9"]=1e-6;
    p["SIGMA_10"]=1e-6;
    p["SIGMA_11"]=1e-6;
    p["SIGMA_12"]=1e-6;
    p["SIGMA_13"]=1e-6;
    p["SIGMA_14"]=1e-6;
    p["SIGMA_15"]=1e-6;
    //do the real work
    
    MaxEntSimulation my_sim(p);
    my_sim.run();
    my_sim.evaluate();
    int gridsize = my_sim.omegaGrid.size();
    
    const double minZero = 1e-4;
    //check endpoints of grid
    EXPECT_NEAR(my_sim.omegaGrid(0),-10,1);
    EXPECT_NEAR(my_sim.omegaGrid(gridsize-1),10,1);
    //endpoints of A(omega) should be <<1
    EXPECT_EQ(my_sim.avspec[0]<minZero,true);
    EXPECT_EQ(my_sim.avspec[gridsize-1]<minZero,true);

    EXPECT_EQ(my_sim.maxspec[0]<minZero,true);
    EXPECT_EQ(my_sim.maxspec[gridsize-1]<minZero,true);

    //check norm
    double max_norm = getNorm(my_sim.omegaGrid,my_sim.maxspec);
    double av_norm = getNorm(my_sim.omegaGrid,my_sim.avspec);
    EXPECT_NEAR(max_norm,1,1e-2);
    EXPECT_NEAR(av_norm,1,1e-2);
    SUCCEED();
}
