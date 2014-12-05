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


bool stop_callback(boost::posix_time::ptime const & end_time) {
  static alps::ngs::signal signal;
  return !signal.empty() || boost::posix_time::second_clock::local_time() > end_time;
}


#ifdef BUILD_PYTHON_MODULE
//compile it as a python module (requires boost::python library)
using namespace boost::python;

void run_it(boost::python::dict parms_){
  alps::parameters_type<MaxEntSimulation>::type parms(parms_);
  std::string out_file = boost::lexical_cast<std::string>(parms["BASENAME"]|"results")+std::string(".out.h5");

#else

  int main(int argc, char** argv)
  {
    {
      if(argc==2 && std::string(argv[1])==std::string("--test"))
        ::testing::InitGoogleTest(&argc, argv);
        exit(RUN_ALL_TESTS());
    }
    alps::mcoptions options(argc, argv);

    alps::parameters_type<MaxEntSimulation>::type parms(alps::hdf5::archive(options.input_file));

    std::string out_file(boost::lexical_cast<std::string>(parms["BASENAME"]|options.output_file)+std::string(".out.h5"));
#endif
    MaxEntSimulation my_sim(parms,out_file); // create a simulation
    my_sim.run(boost::bind(&stop_callback, boost::posix_time::second_clock::local_time() + boost::posix_time::seconds((int)(parms["MAX_TIME"]|60)))); // run the simulation
#ifdef BUILD_PYTHON_MODULE
    return;
#else
    return 0;
#endif
  }

#ifdef BUILD_PYTHON_MODULE
  BOOST_PYTHON_MODULE(maxent_c)
  {
    def("AnalyticContinuation",run_it);//define python-callable run method
  };
#endif

