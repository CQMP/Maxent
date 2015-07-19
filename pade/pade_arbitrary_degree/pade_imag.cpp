/*****************************************************************************
 *
 * ALPS Project Applications
 *
 * Copyright (C) 2013 by Emanuel Gull <egull@umich.edu>
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

#include "pade.hpp"
#include <fstream>


imaginary_domain_data::imaginary_domain_data(const alps::params &p):G_(p){
  //if(!p.defined("NORM")){ throw std::runtime_error("normalization parameter NORM missing!"); }
  norm_=1.;//p["NORM"]; PADE DOES NOT KNOW HOW TO DEAL WITH NORM!
  N_imag_=p["NDAT"];
  val_.resize(N_imag_);
  //std::cerr<<"Data normalized to: "<<norm_<<std::endl;
  //our data is defined in a data file.
  if (p.defined("DATA")) {
    std::string fname = p["DATA"];
    //the data file is a hdf5 data file
    if(p.defined("DATA_IN_HDF5") && p["DATA_IN_HDF5"]) {
      throw std::runtime_error("check and validate hdf5 data reading!");
    } else {
      std::ifstream datstream(fname.c_str());
      if (!datstream)
        boost::throw_exception(std::invalid_argument("could not open data file: "+fname));
      std::cerr<<"reading column text data as #freq real imag"<<std::endl;
      double index, X_i_real,X_i_imag;
      for(int i=0;i<N_imag_;++i){
        //if(!datstream.good()) throw std::runtime_error("problem with reading input data, not enough data or wrong format?");
        datstream >> index >> X_i_real >> X_i_imag >> std::ws;
        std::cout<<" read: "<<index<<" "<<X_i_real<<" "<<X_i_imag<<std::endl;
        val_[i] = std::complex<double>(X_i_real, X_i_imag)/norm_;
      }
    }
  } else {
    std::cerr<<"reading data from parameter file"<<std::endl;
    for (int i=0; i<N_imag_; ++i){
      if(!p.defined("X_Re_"+boost::lexical_cast<std::string>(i))){ throw std::runtime_error("parameter X_Re_i missing!"); }
      if(!p.defined("X_Im_"+boost::lexical_cast<std::string>(i))){ throw std::runtime_error("parameter X_Im_i missing!"); }
      val_[i] = std::complex<double>(p["X_Re_"+boost::lexical_cast<std::string>(i)],p["X_Im_"+boost::lexical_cast<std::string>(i)])/static_cast<double>(p["NORM"]);
    }
  }
}

void imaginary_domain_data::write(const std::string &s) const{
  std::ofstream file(s.c_str());
  for(int i=0;i<N_imag_;++i){
    file<<G_.freq(i)<<" "<<val_[i].real()<<" "<<val_[i].imag()<<std::endl;
  }
}

