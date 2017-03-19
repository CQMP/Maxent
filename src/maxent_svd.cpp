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

#include "maxent.hpp"
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/math/special_functions/legendre.hpp> //needed for Legendre transform
namespace bmth = boost::math;

  // We provide a file with data points and error bars, the latter are used only if
  // COVARIANCE_MATRIX is not set. The format is
  //
  // index data error
  // OR
  // index data error data error
  //
  // index is stored. If tau points are not specified in another location
  // they must be used as the index.
  // Note: for the second option ndat = 2*number of points
  //
  // In case we provide data in a HDF5 file (DATA_IN_HDF5 = 1) we use
  // the following convention:
  // The data are consecutively contained in directory /Data. If we have complex
  // data, we expect the sequence real1 imag1 real2 imag2 ...
  // If we do not provide the covariance matrix, error bars will be read from directory
  // /Error, otherwise the covariance matrix will be read from directory /Covariance
  // as a ndat*ndat field, i.e. it should have been stored according to i*ndat+j
  // As with the data, we adopt the convention that for complex data the sequence
  // will be real imag.



void SVDContinuation::compute_minimal_chi2()const {
  vector_type y2 = U_*(U_.transpose()*y_); //y2 has dimension ndat(), which is dimension of y
  double chi = (y_ - y2).norm(); //this measures the loss of precision when transforming to singular space and back.
  std::cout << "minimal chi2: " << chi * chi / y_.size() << std::endl;
}

void SVDContinuation::truncate_to_singular_space(const vector_type& S) {
  //(truncated) U has dimension ndat() * ns_; ndat() is # of input (matsubara frequency/imag time) points
  //(truncated) Sigma has dimension ns_*ns_ (number of singular eigenvalues)
  //(truncated) V^T has dimensions ns_* nfreq(); nfreq() is number of output (real) frequencies
  //ns_ is the dimension of the singular space.
  U_.conservativeResize(ndat(), ns_);
  Vt_.conservativeResize(ns_, nfreq());
  Sigma_= Eigen::MatrixXd::Zero(ns_,ns_);
  for (int s = 0; s < ns_; ++s) {
    std::cout << "# " << s << "\t" << S[s] << "\n";
    Sigma_(s, s) = S[s];
  }
}

void SVDContinuation::singular_value_decompose_kernel(bool verbose,
    vector_type& S) {
  const double threshold = std::sqrt(std::numeric_limits<double>::epsilon())
      * nfreq();
  Eigen::JacobiSVD<matrix_type> svd(K_,Eigen::ComputeThinU | Eigen::ComputeThinV);
  //svd.setThreshold(threshold);
  S=svd.singularValues();
  U_=svd.matrixU();
  Vt_=svd.matrixV().transpose(); 

  if (verbose)
    std::cout << "# Singular values of the Kernel:\n";

  const double prec = std::sqrt(std::numeric_limits<double>::epsilon())
      * nfreq() * S[0];
  if (verbose)
    std::cout << "# eps = " << sqrt(std::numeric_limits<double>::epsilon())
        << std::endl << "# prec = " << prec << std::endl;

  for (unsigned int s = 0; s < S.size(); ++s) {
    if (verbose)
      std::cout << "# " << s << "\t" << S[s] << "\n";

    ns_ = (S[s] >= prec) ? s + 1 : ns_;
  }
  if (ns() == 0)
    boost::throw_exception(
        std::logic_error("all singular values smaller than the precision"));
  
}
void SVDContinuation::check_high_frequency_limit(const vector_type& y,const kernel_type kt){
    //TODO: expand check into time/bosonic domain
    //we know that the limit of a green's function is ~ -1/iw_{n}
    //this checks that we have a good high frequency limit
    //within the tolerance of the error bar
    
    int n=0;
    if (kt==frequency_fermionic_kernel)
        n = (ndat()-2)/2; //this is the real n; due to weird vector structure
    if (kt==frequency_fermionic_ph_kernel)
        n=ndat()-1;
    
    if(n!=0){
        //limit = G(iwn)*iwn
        std::complex<double> iwn(0,(2*n+1)*M_PI*T());
        double limit = y(ndat()-1)*iwn.imag();
        if(std::abs(1+limit)>sigma(ndat()-1)){
            std::cerr<<"The high frequency limit is not 1!: " << std::abs(limit)
            <<" Check norm?"<< std:: endl;
        }
        //now backcontinue default model and check high frequency limit
        // G(iw_{n})=\sum_{m}K_{nm}A_{m}
        std::complex<double> G;
        iwn = std::complex<double>(iwn.real(), (2*2000+1)*M_PI*T());
        for(int j=0;j<nfreq();j++){
            G+= 1.0/(iwn-omega_coord(j))*Default().D(omega_coord(j)) * delta_omega(j);
        }
        limit = G.imag()*iwn.imag();
        if(std::abs(1+limit)>.01){
            std::cerr<<"The high frequency limit of the default model is not 1!: "
            << std::abs(limit) <<" Check norm?"<< std:: endl;
        }
    }
    //time space=> tail_1 = -G(0)-G(beta) = 1
    if(kt==time_fermionic_kernel){
        double err = sqrt(sigma(0)*sigma(0)+sigma(ndat()-1)*sigma(ndat()-1));
        double limit = -y(0)-y(ndat()-1);
        if(std::abs(limit-1)>err)
            std::cerr<<"The high frequency limit is not +1!: " << limit
            <<" Check norm?"<< std:: endl;
    }
}

SVDContinuation::SVDContinuation(alps::params& p) :
    KernelAndGrid(p),
    Default_(make_default_model(p, "DEFAULT_MODEL")),
    U_(ndat(), ndat()), Vt_(ndat(), nfreq()), Sigma_(ndat(), ndat()),
    omega_coord_(nfreq()), delta_omega_(nfreq()), ns_(0)
{
  for (int i=0; i<nfreq(); ++i) {
    //set the omega_coord into the middle of two map_to_zeroone_interval points
    omega_coord_[i] = (Default().omega_of_t(grid_(i)) + Default().omega_of_t(grid_(i+1)))/2.;
    //and set delta_omega to the distance between the two neighboring map_to_zeroone_interval points
    delta_omega_[i] = Default().omega_of_t(grid_(i+1)) - Default().omega_of_t(grid_(i));
  }
  //build a kernel matrix
  kernel ker(p,omega_coord_,inputGrid_);
  K_=ker();
  k_type = ker.getKernelType();

  //scale lhs and rhs according to errors, etc.
  if (p.defined("COVARIANCE_MATRIX"))
    decompose_covariance_matrix(p);
    
  check_high_frequency_limit(y(),k_type);

  //Look around Eq. D.5 in Sebastian's thesis. We have sigma_ = sqrt(eigenvalues of covariance matrix) or, in case of a diagonal covariance matrix, we have sigma_=SIGMA_X. The then define y := \bar{G}/sigma_ and K := (1/sigma_)\tilde{K}
  scale_data_with_error(nfreq());

  std::cerr << "Kernel is set up\n";

  vector_type S(ndat());
  bool verbose=p["VERBOSE"];
  //perform the SVD decomposition K = U Sigma V^T
  singular_value_decompose_kernel(verbose, S);

  //truncate all matrix to singular space
  truncate_to_singular_space(S);

  //compute Ut and 
  compute_minimal_chi2();
}

void SVDContinuation::define_parameters(alps::params &p){
  KernelAndGrid::define_parameters(p);
  DefaultModel::define_parameters(p);

  p.define<bool>("DATA_IN_HDF5",false,"1 if data is in HDF5 format");
  //TODO: revisit covariance matrix handling.
  //p.define<bool>("COVARIANCE_MATRIX",false,"1 if covariance matrix needs to be diagonalized");
  p.define<std::string>("DATA","","data file input");
  p.define<double>("NORM",1.0,"NORM");
}


