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
#include <alps/hdf5/vector.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/math/special_functions/legendre.hpp> //needed for Legendre transform
namespace bmth = boost::math;
#include <iomanip>

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

///Read data from a text file, with filename given by p["DATA"] in the parameters.
///The format should be index data error
void KernelAndGridIO::read_data_from_text_file(const alps::params& p) {
  std::string fname = p["DATA"];
  std::ifstream datstream(fname.c_str());
  if (!datstream){
    boost::throw_exception(
        std::invalid_argument("could not open data text file: " + fname+". data should be specified in parameter DATA"));
  }
  int datIn=0; //counts up to ndat
  double norm=p["NORM"];
  std::string dataspace = p["DATASPACE"].as<std::string>();
  boost::to_lower(dataspace);
  if(dataspace == "time" || dataspace == "legendre" || p["PARTICLE_HOLE_SYMMETRY"]==true){
    while (datstream) {
      double index, X_i, dX_i;
      datstream >> index >> X_i;
      if(!no_errors_) datstream >> dX_i;
      if (datIn < ndat()) {
        inputGrid_[datIn] = index;
        y_[datIn] = X_i / norm;
        if(!no_errors_)  sigma_[datIn] = dX_i / norm;
        datIn++;
      }
      std::string filename1="y.dat";
      std::ofstream y_file(filename1);
      y_file<<std::setprecision(14);
      for(int i=0;i<ndat();++i){
          y_file<<y_[i]<<std::endl;
      }
    }
    if(ndat()!=datIn){
      throw std::invalid_argument(std::string("The NDAT value ("+boost::lexical_cast<std::string>(ndat_)
          +") is not <= the elements in your input file ("
          +p["DATA"].as<std::string>()+")"));
    }
  }
  else{
    if(ndat()%2 != 0){
      std::cerr << "WARNING: frequency data without particle-hole symmetry"
          << " requires an even amount of input data. Fix parameter file"
          << std::endl;
      throw std::invalid_argument("Your NDAT is odd!");\
    }
    while (datstream) {
      double index, X_i_re, dX_i_re, X_i_im, dX_i_im;
      datstream >> index >> X_i_re;
      if(!no_errors_) datstream >> dX_i_re;
      datstream >> X_i_im;
      if(!no_errors_) datstream>> dX_i_im;
      if (datIn < ndat()) {
        inputGrid_(datIn) = index;
        inputGrid_(datIn+1) = index; 
        y_(datIn) = X_i_re / norm;
        y_(datIn+1) = X_i_im / norm;
        if(!no_errors_){
          sigma_(datIn) = dX_i_re / norm;
          sigma_(datIn+1) = dX_i_im / norm;
        }
        datIn+=2;
      }
    }
    if(ndat()!=datIn){
      throw std::invalid_argument(std::string("The NDAT value ("+boost::lexical_cast<std::string>(ndat_)
          +") is not <= the elements in your input file ("
          +p["DATA"].as<std::string>()+")"));
    }
  }
  if(p.defined("COVARIANCE_MATRIX")) {
    std::string fname = p["COVARIANCE_MATRIX"];
    read_covariance_matrix_from_text_file(fname);
  }

  //check for user error
  /*if(expectedDatIn<ndat()){
    throw std::runtime_error(
        std::string("The NDAT value ("+boost::lexical_cast<std::string>(ndat_) 
            +") is not <= the elements in your input file ("
            +p["DATA"].as<std::string>()+")"));
  }*/
}

///Read data from a hdf5 file, with filename given by p["DATA"] in the parameters.
///The data is stored at /Data.
///The error is stored at /Error
///if the parameter COVARIANCE_MATRIX is specified, then the covariance matrix is read in instead of the error.
///The covariance matrix is expected to be stored at /Covariance

void KernelAndGridIO::read_data_from_hdf5_file(const alps::params& p) {
  std::string fname = p["DATA"];
  //attempt to read from h5 archive
  alps::hdf5::archive ar(fname, "r");
  std::vector<double> tmp(ndat());
  std::stringstream path;
  path << "/Data";
  ar >> alps::make_pvp(path.str(), tmp);
  for (std::size_t i = 0; i < ndat(); i++)
    y_(i) = tmp[i] / p["NORM"].as<double>();
  path.str("");
  if(!no_errors_){
    if (!p.defined("COVARIANCE_MATRIX")) {
      path << "/Error";
      ar >> alps::make_pvp(path.str(), tmp);
      for (std::size_t i = 0; i < ndat(); i++){
        sigma_(i) = tmp[i] / p["NORM"].as<double>();
      }
    } else {
      path << "/Covariance";
      cov_.resize(ndat(), ndat());
      tmp.clear();
      tmp.resize(ndat() * ndat());
      ar >> alps::make_pvp(path.str(), tmp);
      for (std::size_t i = 0; i < ndat(); i++)
        for (std::size_t j = 0; j < ndat(); j++)
          cov_(i, j) = tmp[i * ndat() + j];
    }
  }
}

void KernelAndGridIO::read_data_from_param_file(const alps::params& p) {
  if (!p.defined("NORM")) {
    throw std::runtime_error("parameter NORM missing!");
  } else
    std::cerr << "Data normalized to: " << p["NORM"].as<double>()<< std::endl;

  for (int i = 0; i < ndat(); ++i) {
    if (!p.exists("X_" + boost::lexical_cast<std::string>(i))) {
      throw std::runtime_error("parameter X_"+ boost::lexical_cast<std::string>(i)+ " missing!");
    }
    y_(i) = p["X_" + boost::lexical_cast<std::string>(i)].as<double>()/ p["NORM"].as<double>();
    if(!no_errors_){
      if (!p.defined("COVARIANCE_MATRIX")) {
        if (!p.exists("SIGMA_" + boost::lexical_cast<std::string>(i))) {
          throw std::runtime_error(
              std::string("parameter SIGMA_"+boost::lexical_cast<std::string>(i)+ " missing! "));
        }
        sigma_(i) = p["SIGMA_"+ boost::lexical_cast<std::string>(i)].as<double>()/p["NORM"].as<double>();
      }else{
        throw std::runtime_error("parameter COVARIANCE_MATRIX is defined but there is no way for reading it from a parameter file. Please use a text file instead.");
      }
    }
  }
}

KernelAndGridIO::KernelAndGridIO(alps::params& p) :
    T_(1./p["BETA"].as<double>()),no_errors_(p["NO_ERRORS"]),ndat_(p["NDAT"]), nfreq_(p["NFREQ"]),
    y_(ndat_),sigma_(vector_type::Ones(ndat_)),K_(),grid_(p),inputGrid_(ndat_)
{
  if (ndat_<4) 
    boost::throw_exception(std::invalid_argument("NDAT too small"));

  if (p.defined("DATA") && p["DATA"].as<std::string>() != "") {
    if(p.defined("DATA_IN_HDF5") && (p["DATA_IN_HDF5"])) {
      //attempt to read from h5 archive
      read_data_from_hdf5_file(p);
    } else {
      read_data_from_text_file(p);
    }
  } else {
    //if using input file with X_i, need to define them first
    for (int i = 1; i < ndat(); ++i) {
      //first check for explictly assigned
      std::string x_str = "X_"+boost::lexical_cast<std::string>(i);
      if(!p.defined(x_str)){
        p.define<double>(x_str,"");
        p.define<double>("SIGMA_"+boost::lexical_cast<std::string>(i),"");
      }
    }
    read_data_from_param_file(p);
  }
}
void KernelAndGridIO::define_parameters(alps::params &p){
  p.define<bool>("NO_ERRORS",false,"Set =1 if no error estimates are provided, e.g. for continuation of semi-analytical data");
  kernel::define_parameters(p);
  grid::define_parameters(p);
}


void KernelAndGridIO::scale_data_with_error(const int ntab) {
  //Look around Eq. D.5 in Sebastian's thesis. We have sigma_ = sqrt(eigenvalues of covariance matrix) or, in case of a diagonal covariance matrix, we have sigma_=SIGMA_X. The then define y := \bar{G}/sigma_ and K := (1/sigma_)\tilde{K}
  for (int i = 0; i < ndat(); i++) {
    y_[i] /= sigma_[i];
    for (int j = 0; j < ntab; ++j) {
      K_(i, j) /= sigma_[i];
    }
  }
}

void KernelAndGridIO::read_covariance_matrix_from_text_file(
    const std::string& fname) {
  if(no_errors_) throw std::runtime_error("if you specify no-errors, you should not try to read a covariance matrix.");
  cov_.resize(ndat(), ndat());
  std::cerr << "Reading covariance matrix\n";
  std::ifstream covstream(fname.c_str());
  if (!covstream)
    boost::throw_exception(
        std::invalid_argument(
            "could not open covariance matrix file: " + fname));

  int i, j;
  double covariance;
  while (covstream) {
    covstream >> i >> j >> covariance;
    if (i < ndat() && j < ndat())
      cov_(i, j) = covariance;
  }
}

void KernelAndGridIO::decompose_covariance_matrix(const alps::params& p){
  vector_type var(ndat());
  //bindings::lapack::syev('V', bindings::upper(cov_) , var, bindings::lapack::optimal_workspace());
  //TODO: check if this truly implements lapack's expected overwrite of cov_
  Eigen::SelfAdjointEigenSolver<matrix_type> es(cov_);
  var=es.eigenvalues();
  cov_=es.eigenvectors();
  matrix_type K_loc = cov_.transpose()*K_;
  vector_type y_loc = cov_.transpose()*y_;
  if (p["VERBOSE"])
    std::cout << "# Eigenvalues of the covariance matrix:\n";
  // We drop eigenvalues of the covariance matrix which are smaller than 1e-10
  // as they represent bad data directions (usually there is a steep drop
  // below that value)
  int new_ndat_,old_ndat_=ndat();
  for (new_ndat_ =0;new_ndat_<ndat();new_ndat_++)
    if (var[new_ndat_]>1e-10) break;
  // This is the number of good data
  ndat_ = old_ndat_ - new_ndat_;
  std::cout << "# Ignoring singular eigenvalues (0-" << new_ndat_-1 << " out of " << old_ndat_ << ")\n";
  // Now resize kernel and data matrix and fill it with the values for the
  // good data directions
  K_.resize(ndat_,nfreq_);
  y_.resize(ndat_);
  sigma_.resize(ndat_);
  for (int i=0; i<ndat(); i++) {
    y_(i) = y_loc(new_ndat_+i);
    for (int j=0; j<nfreq_; j++) {
      K_(i,j)=K_loc(new_ndat_+i,j);
    }
  }
  for (int i=0; i<ndat(); ++i) {
    sigma_[i] = std::abs(var(new_ndat_+i))/static_cast<double>(p["NORM"]);
    if (p["VERBOSE"])
      std::cout << "# " << var(new_ndat_+i) << "\n";
  }
}

