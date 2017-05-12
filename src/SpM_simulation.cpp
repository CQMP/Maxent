/*
 * Copyright (C) 1998-2016 ALPS Collaboration
 *
 *     This program is free software; you can redistribute it and/or modify it
 *     under the terms of the GNU General Public License as published by the Free
 *     Software Foundation; either version 2 of the License, or (at your option)
 *     any later version.
 *
 *     This program is distributed in the hope that it will be useful, but WITHOUT
 *     ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *     FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 *     more details.
 *
 *     You should have received a copy of the GNU General Public License along
 *     with this program; if not, write to the Free Software Foundation, Inc., 59
 *     Temple Place, Suite 330, Boston, MA 02111-1307 USA
 * For use in publications, see ACKNOWLEDGE.TXT
 */

#include <iomanip>
#include "SpM_method.hpp"
#include <alps/hdf5/vector.hpp>

SpMSimulation::SpMSimulation(alps::params &parms):
 SVDContinuation(parms)
, norm(parms["NORM"])                                             //The integral is normalized to NORM (use e.g. for self-energies
, Kernel_type(parms["KERNEL"].as<std::string>())
, nfreq(parms["NFREQ"].as<int>())
, rho_(parms["RHO"])
, rhoprime_(parms["RHOPRIME"])
, lambda_(parms["LAMBDA"])
, admm_(Vt(), U().transpose()*y(), delta_omega(), Sigma().diagonal(), rho_, rhoprime_, lambda_)
{
  std::string bn=parms["BASENAME"]; name=bn+'.';

  if(norm != 1.) std::cout<<"Data is assumed to be normalized to: "<<norm<<std::endl;

}
///define parameter defaults
void SpMSimulation::define_parameters(alps::params &p){
  SVDContinuation::define_parameters(p);
  p.description("SpM - a utility for "
    "performing analytic continuation \n \t a la SpM. Experimental, do not use.\n");

  p.define<std::string>("BASENAME","","Specified output name \n(generated if not given)");
  p.define<double>("RHO",20,"ADMM Relaxation Parameter for L1 norm");
  p.define<double>("RHOPRIME",20,"ADMM Relaxation Parameter for positivity");
  p.define<double>("LAMBDA","L1 fit norm");
}

void SpMSimulation::run(){
  std::cout<<"This is an SpM run. "<<std::endl;
  std::cout<<"the Kernel has eigenvalues: "<<Sigma().diagonal()<<std::endl;

  admm_.print_info(std::cout); std::cout<<std::endl;
  for(int j=0;j<100;++j){
    for(int i=0;i<100;++i){
      admm_.iterate();
    }
    admm_.print_info(std::cout); std::cout<<std::endl;
    //std::cout<<"chi2 admm: "<<admm_.chisquare_term()<<" direct: "<<0.5*(y_-K_*(admm_.spectral_function().cwiseProduct(delta_omega()))).squaredNorm()
    //    <<" trans: "<<0.5*(U().transpose()*y_-Sigma()*Vt()*(admm_.spectral_function().cwiseProduct(delta_omega()))).squaredNorm();
  }
}
void SpMSimulation::evaluate(){
  Eigen::VectorXd omega=omega_coord();
  Eigen::VectorXd spectrum=admm_.spectral_function();
  std::string filename="spectrum.dat";
  std::ofstream spec_file(filename);
  for(int i=0;i<nfreq;++i){
    spec_file<<omega[i]<<" "<<spectrum[i]<<" "<<delta_omega(i)<<std::endl;
  }

  /*std::string filename="backcont.dat";
  for(int i=0;i<ndat_;++i){
    spec_file<<inputGrid_[i]<<" "<<y_[i]<<" "<<back(i)<<std::endl;
  }*/
}
