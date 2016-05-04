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


int main(int argc,const char** argv)
{
  alps::params parms(argc,argv); 
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
            //vectors to hold spectra
            std::vector<vector_type> max_spectra(nruns);
            std::vector<vector_type> av_spectra(nruns);

            vector_type omega_grid;
            for(int i=0;i<nruns;i++){
                std::string currModel = parms["RUN_" + boost::lexical_cast<std::string>(i)];
                parms["DEFAULT_MODEL"]= currModel;
                
                //run a simulation with the new default model.
                //Change the basename to match
                parms["BASENAME"] = basename+'.'+currModel;
                MaxEntSimulation my_sim(parms);
                my_sim.run();
                my_sim.evaluate();

                //save spectra for later
                max_spectra[i] = my_sim.getMaxspec();
                av_spectra[i]  = my_sim.getAvspec();
                if(i==0){
                  omega_grid = my_sim.getOmegaGrid();
                }
            }
          
            //calculate variance
            const int nfreq = max_spectra[0].size();
            vector_type mean_max = vector_type::Zero(nfreq);
            vector_type stdev_max(nfreq);
            vector_type mean_av = vector_type::Zero(nfreq);
            vector_type stdev_av(nfreq);

            determineVariance(max_spectra,mean_max,stdev_max);
            determineVariance(av_spectra,mean_av,stdev_av);

            //save to file
            if(parms["TEXT_OUTPUT"]){
            ofstream_ spec_file;
            spec_file.open((basename+".varspec.dat").c_str());
            spec_file << "#omega mean_maxspec stdev_maxspec mean_avspec stdev_avspec" <<std::endl;
            for (std::size_t  i=0; i<omega_grid.size(); ++i)
              spec_file << omega_grid(i)
                << " " << mean_max(i) << " " << stdev_max(i)
                << " " << mean_av(i)  << " " << stdev_av(i) 
                << std::endl;
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

