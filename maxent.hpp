/*****************************************************************************
 *
 * ALPS Project Applications
 *
 * Copyright (C) 2010 by Sebastian Fuchs <fuchs@comp-phys.org>
 *                       Thomas Pruschke <pruschke@comp-phys.org>
 *                       Matthias Troyer <troyer@comp-phys.org>
 *               2011 by Emanuel Gull <gull@phys.columbia.edu>
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

/* $Id: factory.h 1101 2004-08-18 18:11:43Z troyer $ */

#pragma once

#include <fstream>
#include "maxent_matrix_def.hpp"
#include "maxent_parms.hpp"

#include"gtest/gtest.h"

class MaxEntHelper : private MaxEntParameters
{
public : 

  MaxEntHelper(alps::params& p);

  double omega_coord(const int i) const { return MaxEntParameters::omega_coord(i); }

  ///getter function for the discretized default model
  double Default(const int i) const { return def_[i]; }  
  ///getter function for the discretized default model
  const vector_type& Default() const { return def_; }

  //compute the 'free energy' Q (eq. 6.8 in Sebastian's thesis) according to equation D.8 in Sebastian's thesis
  double Q(const vector_type& u, const double alpha) const {
    vector_type A=transform_into_real_space(u);
    return 0.5*chi2(A)-alpha*entropy(A);
  }

  ///number of points for the Matsubara data
  int ndat() const { return MaxEntParameters::ndat(); }

  ///A->u
  vector_type transform_into_singular_space(vector_type A) const;
  ///u-A
  vector_type transform_into_real_space(vector_type u) const;
  ///A=A_i/\sigma_i; this removes \sigma_i
  vector_type get_spectrum(const vector_type& u) const;
  vector_type PrincipalValue(const vector_type &w,const vector_type &a) const;
  /// \Sigma*(V^T*RealSpace(u)*V)*\Sigma
  matrix_type left_side(const vector_type& u) const;
  /// \Sigma*U^T*(K*RealSpace(u)-y)
  vector_type right_side(const vector_type& u) const;
  ///  \delta \dot (V^T*RealSpace(u)*V)
  double step_length(const vector_type& delta, const vector_type& u) const;
  double convergence(const vector_type& u, const double alpha) const;
  double log_prob(const vector_type& u, const double alpha) const;
  double chi_scale_factor(vector_type A, const double chi_sq, const double alpha) const;
  double chi2(const vector_type& u) const;
  void print_chi2(const vector_type& u, std::ostream &os) const;
  /// compute S=\int d\omega(A(\omega)-D(\omega)-A(\omega)\ln[A(\omega)/D(\omega))
  double entropy(const vector_type& u) const;

private:
  ///discretized and normalized version of the default model.
  vector_type def_;
  ///check that the default model is non-zero
  void checkDefaultModel(const vector_type &D) const;
};




class MaxEntSimulation : private MaxEntHelper
{

public:

  ///setup of parameters
  MaxEntSimulation(alps::params& parms);
  ///the maxent calculation
  void run();
  ///the evaluation and writing of files
  void evaluate();
  vector_type levenberg_marquardt(vector_type u, const double alpha) const;
  vector_type iteration(vector_type u, const double alpha, const double mu) const;

private:

  ///grid of alpha values
  vector_type alpha;
  const double norm;
  const int max_it;
  std::string name,Kernel_type;
  bool verbose,text_output,self;
  const int nfreq;

  vector_type lprob;
  vector_type chi_sq;
  std::vector<vector_type> spectra;
  vector_type u;
  ///averaged spectrum
  vector_type avspec;
  ///classic MaxEnt
  vector_type maxspec;
  ///grid of Omega points
  vector_type omegaGrid;
  //posterior probability of the default model
  double postprobdef;
  ///vector of calculated Q values for each alpha iteration
  //this is one method of checking if a minimum was found
  vector_type qvec;

public:
  ///getter for avspec, the averaged spectrum
  const vector_type getAvspec() const{return avspec;}
  ///getter for maxspec, the most probable spectrum
  const vector_type getMaxspec() const{return maxspec;}
  ///getter for the grid of omega points used in determining 
  ///the spectral function A(omega)
  const vector_type getOmegaGrid() const{return omegaGrid;}
  ///getter for the posterior probability of the default model
  const double getPostProb() const{return postprobdef;}
  ///getter for alpha values
  const vector_type getAlphaGrid() const{return alpha;} 
  ///getter for vector of Q value per alpha iteration
  const vector_type getQvec() const{return qvec;}
}; 

