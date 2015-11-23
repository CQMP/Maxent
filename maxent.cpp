/*
 * Copyright (C) 1998-2015 ALPS Collaboration. See COPYRIGHT.TXT
 * All rights reserved. Use is subject to license terms. See LICENSE.TXT
 * For use in publications, see ACKNOWLEDGE.TXT
 */

#include "maxent.hpp"
#include <alps/utilities/remove_extensions.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/exception/diagnostic_information.hpp> 

int main(int argc,char** argv)
{
  //gtest requires char** while alps params requires const char**
  if(argc==2 && std::string(argv[1])==std::string("--test")){
    ::testing::InitGoogleTest(&argc, argv);
    exit(RUN_ALL_TESTS());
  }

  alps::params parms(argc,const_cast<const char**>(argv)); 
  MaxEntSimulation::define_parameters(parms);
  if(!parms.exists("BETA")){
    std::cout<<"Please supply BETA"<<std::endl;
    parms["help"] = true;
  }
  if(!parms.exists("NDAT")){
    std::cout<<"Please supply NDAT"<<std::endl;
    parms["help"] = true;
  }
  if (parms.help_requested(std::cout)) {
    return 0;
  }
  std::string basename;
  if(parms.defaulted("BASENAME")){
    basename = alps::remove_extensions(parms.get_origin_name()) + ".out";
    parms["BASENAME"] = basename;
  }
  else
    basename=parms["BASENAME"].as<std::string>();
  
  try{
        //allow for multiple default model runs
        //set MODEL_RUNS = #runs
        //then place:
        //RUN_0 = "Model1"
        //RUN_1 = "Model2"
        //..etc
        if(parms.exists("MODEL_RUNS")){
            int nruns=parms["MODEL_RUNS"];
            std::cout<<"Performing " << nruns <<" runs" <<std::endl;
	          //ALPSCore requires all params are defined
 	          for(int i=0;i<nruns;i++)
		          parms.define<std::string>("RUN_" + boost::lexical_cast<std::string>(i),"Run");
            for(int i=0;i<nruns;i++){
                std::string currModel = parms["RUN_" + boost::lexical_cast<std::string>(i)];
                parms["DEFAULT_MODEL"]= currModel;
                
                //run a simulation with the new default model.
                //Change the basename to match
                parms["BASENAME"] = basename+'.'+currModel;
                MaxEntSimulation my_sim(parms);
                my_sim.run();
                my_sim.evaluate();
            }
        }
        else{
          MaxEntSimulation * my_sim;
          /*if(parms.exists("RT_TIME"))
            my_sim = new MaxEntSimulationRT (parms);
          else */
            my_sim = new MaxEntSimulation (parms);

          my_sim->run();
          my_sim->evaluate();
        }
  }
    catch(const std::exception &e){
        std::cerr << "Caught Exception " << boost::diagnostic_information(e);
    }
    catch(...){
        std::cerr << "Caught Exception" << boost::current_exception_diagnostic_information();
    }
}

