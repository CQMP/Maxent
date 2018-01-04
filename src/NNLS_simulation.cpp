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
#include "NNLS_method.hpp"
#include <alps/hdf5/vector.hpp>

NNLS_Simulation::NNLS_Simulation(alps::params &parms):
 KernelAndGridIO(parms)
, norm(parms["NORM"])                                             //The integral is normalized to NORM (use e.g. for self-energies
, Kernel_type(parms["KERNEL"].as<std::string>())
, nfreq_(parms["NFREQ"].as<int>())
, omega_coord_(nfreq_)
, delta_omega_(nfreq_)
{
  std::string bn=parms["BASENAME"]; name=bn+'.';

  if(norm != 1.) std::cout<<"Data is assumed to be normalized to: "<<norm<<std::endl;

  for (int i=0; i<nfreq_; ++i) {
    //set the omega_coord into the middle of two grid points
    omega_coord_[i] = (grid_.omega_of_t(grid_(i)) + grid_.omega_of_t(grid_(i+1)))/2.;
    //and set delta_omega to the distance between the two neighboring grid points
    delta_omega_[i] = grid_.omega_of_t(grid_(i+1)) - grid_.omega_of_t(grid_(i));
  }
  //build a kernel matrix
  kernel ker(parms,omega_coord_,inputGrid_);
  K_=ker();
  k_type = ker.getKernelType();
}
///define parameter defaults
void NNLS_Simulation::define_parameters(alps::params &p){
  KernelAndGridIO::define_parameters(p);
  p.define<bool>("DATA_IN_HDF5",false,"1 if data is in HDF5 format");
  p.define<std::string>("DATA","","data file input");
  p.define<double>("NORM",1.0,"NORM");
  p.define<bool>("VERBOSE",true,"1 to print verbose output");

  p.description("NNLS - a utility for "
    "performing analytic continuation using non-linear least squares.\n Experimental, do not use.\n");

  p.define<std::string>("BASENAME","","Specified output name \n(generated if not given)");
}

void NNLS_Simulation::run(){
  std::cout<<"This is an NNLS run. "<<std::endl;


  std::cout<<"the kernel is: "<<K_<<std::endl;
  std::cout<<"the kernel has dimension: "<<K_.rows()<<" "<<K_.cols()<<std::endl;
  std::cout<<"the rhs (Matsubara data) has dimension: "<<y_.size()<<std::endl;

  return;
}
void NNLS_Simulation::evaluate(){

  return;
}
void NNLS_Simulation::compute_first_derivative_matrix(){
  L1_=Eigen::MatrixXd::Zero(nfreq_-2,nfreq_);
  for(int i=2;i<nfreq_-1;++i){
    //-> this needs to be verified
    double d_im2=delta_omega_[i-1]+0.5*(delta_omega_[i-2]+delta_omega_[i]);
    double d_im1=0.5*(delta_omega_[i-1]+delta_omega_[i]);
    double d_ip1=0.5*(delta_omega_[i]+delta_omega_[i+1]);
    double alpha=d_im1*d_ip1/(d_im2*(d_im2+d_ip1)*(d_im2-d_im1));
    double beta=-d_im2*d_ip1/(d_im1*(d_im2-d_im1)*(d_im1+d_ip1));
    double gamma=d_im2*d_im1/(d_ip1*(d_im1+d_ip1)*(d_im2+d_ip1));
    L1_(i-2,i-2)=alpha;
    L1_(i-2,i-1)=beta;
    L1_(i-2,i  )=-(alpha+beta+gamma);
    L1_(i-2,i+1)=gamma;
  }
}
void NNLS_Simulation::compute_second_derivative_matrix(){
  L2_=Eigen::MatrixXd::Zero(nfreq_-2,nfreq_);
  for(int i=2;i<nfreq_-1;++i){
    //-> this needs to be verified

    double d_im2=delta_omega_[i-1]+0.5*(delta_omega_[i-2]+delta_omega_[i]);
    double d_im1=0.5*(delta_omega_[i-1]+delta_omega_[i]);
    double d_ip1=0.5*(delta_omega_[i]+delta_omega_[i+1]);
    double alpha=2*(d_ip1-d_im1)/(d_im2*(d_im2+d_ip1)*(d_im2-d_im1));
    double beta=2*(d_im2-d_ip1)/(d_im1*(d_im2-d_im1)*(d_im1+d_ip1));
    double gamma=2*(d_im2+d_im1)/(d_ip1*(d_im1+d_ip1)*(d_im2+d_ip1));
    L2_(i-2,i-2)=alpha;
    L2_(i-2,i-1)=beta;
    L2_(i-2,i  )=-(alpha+beta+gamma);
    L2_(i-2,i+1)=gamma;
  }
}
