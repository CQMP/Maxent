/*****************************************************************************
 *
 * ALPS Project Applications
 *
 * Copyright (C) 2010 by Sebastian Fuchs <fuchs@comp-phys.org>
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
#include "maxent_kernel.hpp"
#include <alps/config.h> // needed to set up correct bindings
#include <boost/numeric/bindings/ublas.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector_expression.hpp>
#include <boost/numeric/bindings/lapack/driver/syev.hpp>
#include <boost/numeric/bindings/lapack/driver/gesvd.hpp>
#include <boost/numeric/bindings/upper.hpp>
#include <alps/hdf5/vector.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

  // We provide a file with data points and error bars, the latter are used only if
  // COVARIANCE_MATRIX is not set. The format is
  //
  // index data error
  //
  // index is ignored, but MUST be an integer
  // If we wish to continue imaginary frequency data, the structure must be:
  // index_re data_re error_re index_im data_im error_im
  //
  // again, index_re and index_im MUST be integers and are ignored
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
void ContiParameters::read_data_from_text_file(const alps::params& p) {
  std::string fname = p["DATA"];
  std::ifstream datstream(fname.c_str());
  if (!datstream){
    boost::throw_exception(
        std::invalid_argument("could not open data text file: " + fname+"data should be specified in parameter DATA"));
  }
  while (datstream) {
    int i;
    double X_i, dX_i;
    datstream >> i >> X_i >> dX_i;
    if (i < ndat()) {
      y_(i) = X_i / static_cast<double>(p["NORM"]);
      sigma_(i) = dX_i / static_cast<double>(p["NORM"]);
    }
  }
  if(p.defined("COVARIANCE_MATRIX")) {
    std::string fname = p["COVARIANCE_MATRIX"];
    read_covariance_matrix_from_text_file(fname);
  }
}

///Read data from a hdf5 file, with filename given by p["DATA"] in the parameters.
///The data is stored at /Data.
///The error is stored at /Error
///if the parameter COVARIANCE_MATRIX is specified, then the covariance matrix is read in instead of the error.
///The covariance matrix is expected to be stored at /Covariance

void ContiParameters::read_data_from_hdf5_file(const alps::params& p) {
  std::string fname = p["DATA"];
  //attempt to read from h5 archive
  alps::hdf5::archive ar(fname, alps::hdf5::archive::READ);
  std::vector<double> tmp(ndat());
  std::stringstream path;
  path << "/Data";
  ar >> alps::make_pvp(path.str(), tmp);
  for (std::size_t i = 0; i < ndat(); i++)
    y_(i) = tmp[i] / static_cast<double>(p["NORM"]);
  path.str("");
  if (!p.defined("COVARIANCE_MATRIX")) {
    path << "/Error";
    ar >> alps::make_pvp(path.str(), tmp);
    for (std::size_t i = 0; i < ndat(); i++)
      sigma_(i) = tmp[i] / static_cast<double>(p["NORM"]);
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

void ContiParameters::read_data_from_param_file(const alps::params& p) {
  if (!p.defined("NORM")) {
    throw std::runtime_error("parameter NORM missing!");
  } else
    std::cerr << "Data normalized to: " << static_cast<double>(p["NORM"])
        << std::endl;

  for (int i = 0; i < ndat(); ++i) {
    if (!p.defined("X_" + boost::lexical_cast<std::string>(i))) {
      throw std::runtime_error("parameter X_i missing!");
    }
    y_(i) = static_cast<double>(p["X_" + boost::lexical_cast<std::string>(i)])
        / static_cast<double>(p["NORM"]);
    if (!p.defined("COVARIANCE_MATRIX")) {
      if (!p.defined("SIGMA_" + boost::lexical_cast<std::string>(i))) {
        throw std::runtime_error(
            std::string("parameter SIGMA_i missing!") + "SIGMA_"
                + boost::lexical_cast<std::string>(i));
      }
      sigma_(i) = static_cast<double>(p["SIGMA_"
          + boost::lexical_cast<std::string>(i)])
          / static_cast<double>(p["NORM"]);
    }
  }
}

ContiParameters::ContiParameters(const alps::params& p) :
T_(p["T"]|1./static_cast<double>(p["BETA"])),
ndat_(p["NDAT"]), nfreq_(p["NFREQ"]),
y_(ndat_),sigma_(ndat_),K_(),grid_(p)
{
  if (ndat_<4) 
    boost::throw_exception(std::invalid_argument("NDAT too small"));


  if (p.defined("DATA")) {
    if(p.defined("DATA_IN_HDF5") && (p["DATA_IN_HDF5"]|false)) {
      //attempt to read from h5 archive
      read_data_from_hdf5_file(p);
    } else {
      read_data_from_text_file(p);
    }
  } else {
    read_data_from_param_file(p);
  }
}

void ContiParameters::enforce_strict_normalization(double sigma_normalization,
    double norm, const int ntab) {
  std::cout << "enforcing strict normalization." << std::endl;
  double artificial_norm_enforcement_sigma = sigma_normalization / norm;
  for (int j = 0; j < ntab; ++j) {
    K_(ndat() - 1, j) = 1. / artificial_norm_enforcement_sigma;
  }
  y_[ndat() - 1] = 1. / artificial_norm_enforcement_sigma;
}

void ContiParameters::scale_data_with_error(const int ntab) {
  //Look around Eq. D.5 in Sebastian's thesis. We have sigma_ = sqrt(eigenvalues of covariance matrix) or, in case of a diagonal covariance matrix, we have sigma_=SIGMA_X. The then define y := \bar{G}/sigma_ and K := (1/sigma_)\tilde{K}
  for (int i = 0; i < ndat(); i++) {
    y_[i] /= sigma_[i];
    for (int j = 0; j < ntab; ++j) {
      K_(i, j) /= sigma_[i];
    }
  }
}

void ContiParameters::read_covariance_matrix_from_text_file(
    const std::string& fname) {
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

void ContiParameters::decompose_covariance_matrix(const alps::params& p){
  using namespace boost::numeric;
    vector_type var(ndat());
    bindings::lapack::syev('V', bindings::upper(cov_) , var, bindings::lapack::optimal_workspace());
    matrix_type cov_trans = ublas::trans(cov_);
    matrix_type K_loc = ublas::prec_prod(cov_trans, K_);
    vector_type y_loc = ublas::prec_prod(cov_trans, y_);
    if (p["VERBOSE"]|false)
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
      if (p["VERBOSE"]|false)
        std::cout << "# " << var(new_ndat_+i) << "\n";
    }
}

void MaxEntParameters::compute_minimal_chi2()const {
  using namespace boost::numeric;
  matrix_type Ut = ublas::trans(U_); //U^T has dimension ns_*ndat()
  vector_type t = ublas::prec_prod(Ut, y_); //t has dimension ns_
  vector_type y2 = ublas::prec_prod(U_, t); //y2 has dimension ndat(), which is dimension of y
  double chi = ublas::norm_2(y_ - y2); //this measures the loss of precision when transforming to singular space and back.
  std::cout << "minimal chi2: " << chi * chi / y_.size() << std::endl;
}

void MaxEntParameters::truncate_to_singular_space(const vector_type& S) {
  //(truncated) U has dimension ndat() * ns_; ndat() is # of input (matsubara frequency/imag time) points
  //(truncated) Sigma has dimension ns_*ns_ (number of singular eigenvalues)
  //(truncated) V^T has dimensions ns_* nfreq(); nfreq() is number of output (real) frequencies
  //ns_ is the dimension of the singular space.
  U_.resize(ndat(), ns_, true);
  Vt_.resize(ns_, nfreq(), true);
  Sigma_.resize(ns_, ns_);
  Sigma_.clear();
  for (int s = 0; s < ns_; ++s) {
    std::cout << "# " << s << "\t" << S[s] << "\n";
    Sigma_(s, s) = S[s];
  }
}

void MaxEntParameters::singular_value_decompose_kernel(bool verbose,
    vector_type& S) {
  matrix_type Kt = K_; // gesvd destroys K!
  boost::numeric::bindings::lapack::gesvd('S', 'S', Kt, S, U_, Vt_);
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

MaxEntParameters::MaxEntParameters(const alps::params& p) :
    ContiParameters(p),
    Default_(make_default_model(p, "DEFAULT_MODEL")),
    U_(ndat(), ndat()), Vt_(ndat(), nfreq()), Sigma_(ndat(), ndat()),
    omega_coord_(nfreq()), delta_omega_(nfreq()), ns_(0)
{

  for (int i=0; i<nfreq(); ++i) {
    omega_coord_[i] = (Default().omega_of_t(grid_(i)) + Default().omega_of_t(grid_(i+1)))/2.;
    delta_omega_[i] = Default().omega_of_t(grid_(i+1)) - Default().omega_of_t(grid_(i));
  }

  //build a kernel matrix
  kernel ker(p,omega_coord_);
  K_=ker();

  //scale lhs and rhs according to errors, etc.
  if (p.defined("COVARIANCE_MATRIX"))
    decompose_covariance_matrix(p);

  //Look around Eq. D.5 in Sebastian's thesis. We have sigma_ = sqrt(eigenvalues of covariance matrix) or, in case of a diagonal covariance matrix, we have sigma_=SIGMA_X. The then define y := \bar{G}/sigma_ and K := (1/sigma_)\tilde{K}
  scale_data_with_error(ndat());

  //this enforces a strict normalization if needed. not sure that this is done properly. recheck!
  if(p["ENFORCE_NORMALIZATION"]|false) {
    double sigma_normalization=p["SIGMA_NORMALIZATION"];
    double norm=p["NORM"];
    enforce_strict_normalization(sigma_normalization, norm, ndat());
  }
  std::cerr << "Kernel set up\n";

  vector_type S(ndat());
  bool verbose=p["VERBOSE"]|false;
  //perform the SVD decomposition K = U Sigma V^T
  singular_value_decompose_kernel(verbose, S);

  //truncate all matrix to singular space
  truncate_to_singular_space(S);

  //compute Ut and 
  compute_minimal_chi2();
}




