#include "SpM_method.hpp"
#include <alps/utilities/fs/remove_extensions.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/exception/diagnostic_information.hpp>


int main(int argc,const char** argv)
{
  alps::params parms(argc,argv);
  SpMSimulation::define_parameters(parms);
  //other help messages
  parms.define("help.grids","show help for grids");

  if (parms.help_requested(std::cout)) {
    return 0;
  }
  //currently ALPSCore doesn't have a good way
  //to show multilple error messages
  //This is a temp hack, so we can have
  //inline grid and default model
  //descriptions

  bool exitEarly = false;

  if(parms["help.grids"]){
    if(exitEarly) //display both messages
      std::cout << std::endl;
    grid::print_help();
    exitEarly = true;
  }

  if(exitEarly)
    return 0;

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
    basename = alps::fs::remove_extensions(parms.get_origin_name()) + ".out";
    parms["BASENAME"] = basename;
  }
  else
    basename=parms["BASENAME"].as<std::string>();

  try{
    if(exitEarly){
      throw std::runtime_error("Critical parameters not defined");
    }


    SpMSimulation my_sim(parms);
    my_sim.run();
    my_sim.evaluate();
  }
  catch(const std::exception &e){
    std::cerr << "Caught Exception " << boost::diagnostic_information(e);
  }
  catch(...){
    std::cerr << "Caught Exception" << boost::current_exception_diagnostic_information();
  }
}
