#include "maxent.hpp"
#include "gtest/gtest.h"
#include <iostream>
#include "maxent_parms_default.hpp"
double getNorm(const vector_type &omega, const vector_type &y){
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
    p["KERNEL"]="fermionic";
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
    int gridsize = my_sim.getOmegaGrid().size();
    
		//legendre has a harder time with high freq
    const double minZero = 1e-3;
    //check endpoints of grid
    EXPECT_NEAR(my_sim.getOmegaGrid()(0),-10,1);
    EXPECT_NEAR(my_sim.getOmegaGrid()(gridsize-1),10,1);
    //endpoints of A(omega) should be <<1
    EXPECT_EQ(my_sim.getAvspec()[0]<minZero,true);
    EXPECT_EQ(my_sim.getAvspec()[gridsize-1]<minZero,true);

    EXPECT_EQ(my_sim.getMaxspec()[0]<minZero,true);
    EXPECT_EQ(my_sim.getMaxspec()[gridsize-1]<minZero,true);

		//flat default model is not very Lorentzian
		EXPECT_EQ(my_sim.getPostProb()<minZero,true);

		//expect a converged solution/good minimum found
		vector_type q = my_sim.getQvec();
		EXPECT_EQ(q[q.size()-1]<1,true);

    //check norm
    double max_norm = getNorm(my_sim.getOmegaGrid(),my_sim.getMaxspec());
    double av_norm = getNorm(my_sim.getOmegaGrid(),my_sim.getAvspec());
    EXPECT_NEAR(max_norm,1,1e-2);
    EXPECT_NEAR(av_norm,1,1e-2);
    SUCCEED();
}
TEST(Simulation,TauSimulation){
    alps::params p;
    set_defaults(p);
    p["BETA"]=8;
    p["NDAT"]=26;
    p["PARTICLE_HOLE_SYMMETRY"]=1;
    p["DATASPACE"]="time";
    p["KERNEL"]="fermionic";
    p["TEXT_OUTPUT"]=0;


    //data from DMFT;U=0;beta=8
    p["X_0"]=-0.5;
    p["X_1"]=-0.31968831176465;
    p["X_2"]=-0.22989567294666;
    p["X_3"]=-0.18028809394476;
    p["X_4"]=-0.15025032491734;
    p["X_5"]=-0.13069255038955;
    p["X_6"]=-0.11727844091991;
    p["X_7"]=-0.10777271410719;
    p["X_8"]=-0.10094024798648;
    p["X_9"]=-0.096059700592185;
    p["X_10"]=-0.09269474697556;
    p["X_11"]=-0.090580337980363;
    p["X_12"]=-0.08956408577659;
    p["X_13"]=-0.089576497236633;
    p["X_14"]=-0.090618406660374;
    p["X_15"]=-0.09276112281222;
    p["X_16"]=-0.096159341080863;
    p["X_17"]=-0.10108143619853;
    p["X_18"]=-0.10796902395697;
    p["X_19"]=-0.11755259463163;
    p["X_20"]=-0.13108435868405;
    p["X_21"]=-0.15083404013491;
    p["X_22"]=-0.18121280672505;
    p["X_23"]=-0.23148524272195;
    p["X_24"]=-0.32270252375015;
    p["X_25"]=-0.5;

    p["TAU_0"]=0;
    p["TAU_1"]=0.3203125;
    p["TAU_2"]=0.640625;
    p["TAU_3"]=0.9609375;
    p["TAU_4"]=1.28125;
    p["TAU_5"]=1.6015625;
    p["TAU_6"]=1.921875;
    p["TAU_7"]=2.2421875;
    p["TAU_8"]=2.5625;
    p["TAU_9"]=2.8828125;
    p["TAU_10"]=3.203125;
    p["TAU_11"]=3.5234375;
    p["TAU_12"]=3.84375;
    p["TAU_13"]=4.1640625;
    p["TAU_14"]=4.484375;
    p["TAU_15"]=4.8046875;
    p["TAU_16"]=5.125;
    p["TAU_17"]=5.4453125;
    p["TAU_18"]=5.765625;
    p["TAU_19"]=6.0859375;
    p["TAU_20"]=6.40625;
    p["TAU_21"]=6.7265625;
    p["TAU_22"]=7.046875;
    p["TAU_23"]=7.3671875;
    p["TAU_24"]=7.6875;
    p["TAU_25"]=8;

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
    p["SIGMA_16"]=1e-06;
    p["SIGMA_17"]=1e-06;
    p["SIGMA_18"]=1e-06;
    p["SIGMA_19"]=1e-06;
    p["SIGMA_20"]=1e-06;
    p["SIGMA_21"]=1e-06;
    p["SIGMA_22"]=1e-06;
    p["SIGMA_23"]=1e-06;
    p["SIGMA_24"]=1e-06;
    p["SIGMA_25"]=1e-06;

    //do the real work
    
    MaxEntSimulation my_sim(p);
    my_sim.run();
    my_sim.evaluate();
    int gridsize = my_sim.getOmegaGrid().size();
    
		//legendre has a harder time with high freq
    const double minZero = 1e-3;
    //check endpoints of grid
    EXPECT_NEAR(my_sim.getOmegaGrid()(0),-10,1);
    EXPECT_NEAR(my_sim.getOmegaGrid()(gridsize-1),10,1);
    //endpoints of A(omega) should be <<1
    EXPECT_EQ(my_sim.getAvspec()[0]<minZero,true);
    EXPECT_EQ(my_sim.getAvspec()[gridsize-1]<minZero,true);

    EXPECT_EQ(my_sim.getMaxspec()[0]<minZero,true);
    EXPECT_EQ(my_sim.getMaxspec()[gridsize-1]<minZero,true);

		//flat default model is not very Lorentzian
		EXPECT_EQ(my_sim.getPostProb()<minZero,true);

		//expect a converged solution/good minimum found
		vector_type q = my_sim.getQvec();
		EXPECT_EQ(q[q.size()-1]<1,true);

    //check norm
    double max_norm = getNorm(my_sim.getOmegaGrid(),my_sim.getMaxspec());
    double av_norm = getNorm(my_sim.getOmegaGrid(),my_sim.getAvspec());
    EXPECT_NEAR(max_norm,1,1e-2);
    EXPECT_NEAR(av_norm,1,1e-2);
    SUCCEED();
}


TEST(Simulation,LegendreSimulation){
    alps::params p;
    set_defaults(p);
    p["BETA"]=8;
    p["NDAT"]=13;
    p["PARTICLE_HOLE_SYMMETRY"]=1;
    p["DATASPACE"]="legendre";
    p["KERNEL"]="fermionic";
    p["TEXT_OUTPUT"]=0;



		p["X_0"]=-1.2385922211475;
		p["X_1"]=-7.0664305609387e-08;
		p["X_2"]=-0.66408457100986;
		p["X_3"]=-8.5145830727646e-08;
		p["X_4"]=-0.26534486436104;
		p["X_5"]=2.0234010147938e-07;
		p["X_6"]=-0.091681224519598;
		p["X_7"]=3.9397565117767e-07;
		p["X_8"]=-0.027019280247812;
		p["X_9"]=2.0785849309379e-07;
		p["X_10"]=-0.0067707643917139;
		p["X_11"]=-2.8643971926556e-07;
		p["X_12"]=-0.0014885447940878;


		p["SIGMA_0"]=2.5043167301794e-07;
		p["SIGMA_1"]=2.6594271923164e-07;
		p["SIGMA_2"]=3.0438230279004e-07;
		p["SIGMA_3"]=2.5006246215834e-07;
		p["SIGMA_4"]=2.815176991835e-07;
		p["SIGMA_5"]=2.6399861311466e-07;
		p["SIGMA_6"]=2.6165696141285e-07;
		p["SIGMA_7"]=2.8880049016989e-07;
		p["SIGMA_8"]=2.4924030790157e-07;
		p["SIGMA_9"]=2.5205277343449e-07;
		p["SIGMA_10"]=2.2929459388671e-07;
		p["SIGMA_11"]=2.4635379139221e-07;
		p["SIGMA_12"]=2.4614799010746e-07;

    //do the real work
    
    MaxEntSimulation my_sim(p);
    my_sim.run();
    my_sim.evaluate();
    int gridsize = my_sim.getOmegaGrid().size();
    
		//legendre has a harder time with high freq
    const double minZero = 1e-3;
    //check endpoints of grid
    EXPECT_NEAR(my_sim.getOmegaGrid()(0),-10,1);
    EXPECT_NEAR(my_sim.getOmegaGrid()(gridsize-1),10,1);
    //endpoints of A(omega) should be <<1
    EXPECT_EQ(my_sim.getAvspec()[0]<minZero,true);
    EXPECT_EQ(my_sim.getAvspec()[gridsize-1]<minZero,true);

    EXPECT_EQ(my_sim.getMaxspec()[0]<minZero,true);
    EXPECT_EQ(my_sim.getMaxspec()[gridsize-1]<minZero,true);

		//flat default model is not very Lorentzian
		EXPECT_EQ(my_sim.getPostProb()<minZero,true);

		//expect a converged solution/good minimum found
		vector_type q = my_sim.getQvec();
		EXPECT_EQ(q[q.size()-1]<1,true);

    //check norm
    double max_norm = getNorm(my_sim.getOmegaGrid(),my_sim.getMaxspec());
    double av_norm = getNorm(my_sim.getOmegaGrid(),my_sim.getAvspec());
    EXPECT_NEAR(max_norm,1,1e-2);
    EXPECT_NEAR(av_norm,1,1e-2);
    SUCCEED();
}

