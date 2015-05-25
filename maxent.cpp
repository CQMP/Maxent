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
#include <alps/ngs/mcoptions.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem/path.hpp>

int main(int argc, char** argv)
{
  if(argc==2 && std::string(argv[1])==std::string("--test")){
    ::testing::InitGoogleTest(&argc, argv);
    exit(RUN_ALL_TESTS());
  }
  alps::mcoptions options(argc, argv);
  boost::filesystem::path input_path(options.input_file);
  alps::params parms;
  if(input_path.extension() == ".h5")
    parms =alps::params(alps::hdf5::archive(options.input_file));
  else
    parms =alps::params(options.input_file);
   
  std::string basename;
  if(!parms.defined("BASENAME")){
    parms["BASENAME"]=options.output_file;
    basename=options.output_file;
    }
  else
    basename=boost::lexical_cast<std::string>(parms["BASENAME"]);
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

