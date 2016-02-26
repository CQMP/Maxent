/*
 * Copyright (C) 1998-2016 ALPS Collaboration. See COPYRIGHT.TXT
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
  if (parms.help_requested(std::cout)) {
    return 0;
  }
  bool exitEarly = false;

  if(!parms.exists("BETA")){
    std::cout<<"Please supply BETA"<<std::endl;
    exitEarly = true;
  }
  if(!parms.exists("NDAT")){
    std::cout<<"Please supply NDAT"<<std::endl;
    exitEarly = true;
  }
  if(parms["DATA"].as<std::string>().empty() && !parms.exists("X_0")
      && parms["DATA_IN_HDF5"]!=true){
    std::cout<<"Please supply input data"<<std::endl;
    exitEarly = true;
  }

  std::string basename;
  if(parms.defaulted("BASENAME")){
    basename = alps::remove_extensions(parms.get_origin_name()) + ".out";
    parms["BASENAME"] = basename;
  }
  else
    basename=parms["BASENAME"].as<std::string>();
  
  try{
        if(exitEarly){
          throw std::runtime_error("Critical parameters not defined");
        }

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
          MaxEntSimulation my_sim(parms);
          my_sim.run();
          my_sim.evaluate();
        }
  }
    catch(const std::exception &e){
        std::cerr << "Caught Exception " << boost::diagnostic_information(e);
    }
    catch(...){
        std::cerr << "Caught Exception" << boost::current_exception_diagnostic_information();
    }
}

