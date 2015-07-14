/*****************************************************************************
 *
 * ALPS Project Applications
 *
 * Copyright (C) 2010 by Sebastian  Fuchs <fuchs@comp-phys.org>
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

#include "maxent.hpp"
//#include <alps/mc/mcoptions.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/exception/diagnostic_information.hpp> 
#include "maxent_parms_default.hpp"

inline void checkInput(alps::params &p){
    if(!p.exists("BETA")){
	std::cout<<"Please supply BETA"<<std::endl;
	p["help"] = true;
    }
    if(!p.exists("NDAT")){
	std::cout<<"Please supply NDAT"<<std::endl;
	p["help"] = true;
    }
}

int main(int argc,const char** argv)
{
  //gtest requires char** while alps params requires const char**
    char ** argv_;
    size_t strlen_sum = 0;
    for (int i = 0; i < argc; i++) strlen_sum += strlen(argv[i]) + 1;
    argv_=(char **) malloc(sizeof(char *) * (argc + 1) + strlen_sum);
    for(int i=0;i<argc;i++){
        int width = strlen(argv[i]) + 1;
        argv_[i] = (char *)malloc(width*sizeof(char));
        memcpy(argv_[i], argv[i], width);
    }
  
  if(argc==2 && std::string(argv[1])==std::string("--test")){
    ::testing::InitGoogleTest(&argc, argv_);
    exit(RUN_ALL_TESTS());
  }
  //alps::mcoptions options(argc, argv);
  /*boost::filesystem::path input_path();
  alps::params parms;
  if(input_path.extension() == ".h5")
    parms =alps::params(alps::hdf5::archive(options.input_file));
  else
    parms =alps::params(options.input_file);
  */
  alps::params parms(argc,argv); 
  set_defaults(parms);
  checkInput(parms);
    if (parms.help_requested(std::cout)) {
        return 0;
    }
  std::string basename;
  //if basename="" then make a better one
  if(parms.defaulted("BASENAME")){
    std::string input_file = argv[1];
    if (input_file.find(".in.h5") != std::string::npos)
	parms["BASENAME"] = input_file.substr(0,input_file.find_last_of(".in.h5")-5)+ ".out.h5";
    else if (input_file.find(".out.h5") != std::string::npos)
	 parms["BASENAME"] = input_file;
    else
         parms["BASENAME" ]= input_file.substr(0,input_file.find_last_of('.'))+ ".out";
    basename=parms["BASENAME"].as<std::string>();
  }
  else
    basename=boost::lexical_cast<std::string>(parms["BASENAME"]);
      // If requested, we print the help message, which is constructed from the
    // information we gave when defining the parameters.
  try{
        //allow for multiple default model runs
        //set MODEL_RUNS = #runs
        //then place:
        //RUN_0 = "Model1"
        //RUN_1 = "Model2"
        //..etc
        if(parms.defined("MODEL_RUNS")){
            int nruns=parms["MODEL_RUNS"];
            std::cout<<"Performing " << nruns <<" runs" <<std::endl;
            for(int i=0;i<nruns;i++){
                if(!parms.defined("RUN_" + boost::lexical_cast<std::string>(i))) {
                    throw std::runtime_error("parameter RUN_i missing!");
                }
                std::string currModel = boost::lexical_cast<std::string>(
                                             parms["RUN_" + boost::lexical_cast<std::string>(i)]);
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

