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

