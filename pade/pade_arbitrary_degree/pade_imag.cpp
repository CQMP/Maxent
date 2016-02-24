/*
 * Copyright (C) 1998-2016 ALPS Collaboration. See COPYRIGHT.TXT
 * All rights reserved. Use is subject to license terms. See LICENSE.TXT
 * For use in publications, see ACKNOWLEDGE.TXT
 */

#include "pade.hpp"
#include <fstream>


imaginary_domain_data::imaginary_domain_data(const PadeParams &p):G_(p){
  //if(!p.defined("NORM")){ throw std::runtime_error("normalization parameter NORM missing!"); }
  norm_=1.;//p["NORM"]; PADE DOES NOT KNOW HOW TO DEAL WITH NORM!
  N_imag_=p["imag.NDAT"];
  val_.resize(N_imag_);
  
  //our data is defined in a data file.
  std::string fname = p["imag.DATA"];
  std::ifstream datstream(fname.c_str());
  if (!datstream)
    boost::throw_exception(std::invalid_argument("could not open data file: "+fname));
  std::cout<<"reading column text data as #freq real imag"<<std::endl;
  double index, X_i_real,X_i_imag;
  for(int i=0;i<N_imag_;++i){
    datstream >> index >> X_i_real >> X_i_imag >> std::ws;
    std::cout<<" read: "<<index<<" "<<X_i_real<<" "<<X_i_imag<<std::endl;
    val_[i] = std::complex<double>(X_i_real, X_i_imag)/norm_;
    if(std::abs(G_.freq(i)-index)>1.e-5) throw std::invalid_argument("Grid mismatch. make sure the first entry contains the grid points!!");
  }
}

void imaginary_domain_data::write(const std::string &s) const{
  std::ofstream file(s.c_str());
  for(int i=0;i<N_imag_;++i){
    file<<G_.freq(i)<<" "<<val_[i].real()<<" "<<val_[i].imag()<<std::endl;
  }
}

