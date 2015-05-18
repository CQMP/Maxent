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


#pragma once

#include <math.h>
#include <alps/ngs/params.hpp>
#include <boost/shared_ptr.hpp>


//Note the slightly crooked structure here:

//class DefaultModel
// -> FlatDefaultModel : public DefaultModel
// -> GeneralDefaultModel : public DefaultModel
//  and GeneralDefaultModel contains a 'Model' object.

//class Model
// -> Gaussian : public Model 
//    -> ShiftedGaussian : public Gaussian
//      -> DoubleGaussian : public ShiftedGaussian
//        -> GeneralDoubleGaussian : public ShiftedGaussian
// ->  TabFunction : public Model
// ->  LinearRiseExpDecay : public Model
// ->  QuadraticRiseExpDecay : public Model

class DefaultModel 
{
public:
  DefaultModel(const alps::params& p) :
    omega_max(p["OMEGA_MAX"]),
    omega_min(static_cast<double>(p["OMEGA_MIN"]|(-omega_max))){ //we had a 0 here in the bosonic case. That's not a good idea if you're continuing symmetric functions like chi(omega)/omega. Change omega_min to zero manually if you need it.
  }

  virtual ~DefaultModel(){}

  ///omega returns a frequency point, given x between 0 and 1.
  virtual double omega(const double x) const = 0;

  ///D returns the derivative of the integrated default model
  virtual double D(const double omega) const = 0;

  ///returns the integrated default model
  virtual double x(const double t=0) const = 0;

  ///equidistant mapping from [0,1] to [omega_min, omega_max]
  double omega_of_t(const double t) const { return omega_min + (omega_max-omega_min)*t; }

  ///equidistant mapping from [omega_min, omega_max] to [0,1]
  double t_of_omega(const double omega) const { return (omega-omega_min)/(omega_max-omega_min); }

protected:
  const double omega_max;
  const double omega_min;
};



class FlatDefaultModel : public DefaultModel
{
public:

  FlatDefaultModel(const alps::params& p) : DefaultModel(p) {}

  double omega(const double x) const {
    return x*(omega_max-omega_min) + omega_min;
  }

  double D(const double) const {
    return 1./(omega_max-omega_min);
  }

  double x(const double t) const {
    return t;
  }

};


///general class for a model function. Implements operator() which, given a frequency omega, will give back the value of the default model at that frequency.
class Model
{
public:
  virtual double operator()(const double omega)=0;
  virtual ~Model(){}
};


///a model function that implements a Gaussian, i.e.
/// \f$D(\omega)= \frac{1}{\sqrt{2\pi}\sigma}e^{-\frac{\omega^2}{2\sigma^2}} \f$.
///
/// \f$\sigma\f$ has to be specified as a parameter SIGMA.
/// This is one of the standard model functions you should always try.
class Gaussian : public Model 
{
public:
  Gaussian(const alps::params& p) : sigma_(static_cast<double>(p["SIGMA"])) {}

  virtual double operator()(const double omega) {
    return std::exp(-omega*omega/2./sigma_/sigma_)/sqrt(2*M_PI)/sigma_;
  }

private:
  const double sigma_;
};

///a model function that implements two Gaussians with arbitrary relative norms and arbitrary points around which they are centered.
///
/// \f$\sigma\f$ has to be specified as a parameter SIGMA1 and SIGMA2 for the two gaussians.
/// shift has to be defined as SHIFT1 and SHIFT2 for each Gaussian.
/// NORM1 has to be defined to specify the relative weight of the two spectral functions.
///
/// \f$D(\omega)=\frac{norm1}{\sqrt{2\pi}\sigma1}e^{-\frac{(\omega-shift1)^2}{2\sigma1^2}}+ \frac{(1-norm1)}{\sqrt{2\pi}\sigma2}e^{-\frac{(\omega-shift2)^2}{2\sigma2^2}}\f$

class TwoGaussians : public Model
{
public:
  TwoGaussians(const alps::params& p) : sigma1(static_cast<double>(p["SIGMA1"])),
  sigma2(static_cast<double>(p["SIGMA2"])),
  shift1(static_cast<double>(p["SHIFT1"]|0.0)),
  shift2(static_cast<double>(p["SHIFT2"])),
  norm1(static_cast<double>(p["NORM1"]|0.5)) {}

  virtual double operator()(const double omega) {
    return norm1*std::exp(-(omega-shift1)*(omega-shift1)/2./sigma1/sigma1)/sqrt(2*M_PI)/sigma1+(1.0-norm1)*std::exp(-(omega-shift2)*(omega-shift2)/2./sigma2/sigma2)/sqrt(2*M_PI)/sigma2;
  }

private:
  const double sigma1,sigma2,shift1,shift2,norm1;
};


///a model function that implements a Gaussian that is not centered at zero but at some other frequency, i.e.
/// \f$D(\omega)= \frac{1}{\sqrt{2\pi}\sigma}e^{-\frac{(\omega-shift)^2}{2\sigma^2}} \f$.
///
/// \f$\sigma\f$ has to be specified as a parameter SIGMA, and shift has to be specified as parameter SHIFT.
class ShiftedGaussian : public Gaussian
{
public:
  ShiftedGaussian(const alps::params& p) :
    Gaussian(p), shift_(static_cast<double>(p["SHIFT"])){}

  double operator()(const double omega) {
    return Gaussian::operator()(omega-shift_);
  }

protected:
  const double shift_;
};


///a model function that implements a sum of two Gaussians, each of them shifted by +/- shift.
/// \f$D(\omega)=\frac{1}{2} \left( \frac{1}{\sqrt{2\pi}\sigma}e^{-\frac{(\omega-shift)^2}{2\sigma^2}}+ \frac{1}{\sqrt{2\pi}\sigma}e^{-\frac{(\omega+shift)^2}{2\sigma^2}}\right)\f$
/// \f$\sigma\f$ has to be specified as a parameter SIGMA, and \f$shift_\f$ has to be specified as parameter SHIFT.
/// Try this to try to model a system with 'side peaks'.
class DoubleGaussian : public ShiftedGaussian
{
public:
  DoubleGaussian(const alps::params& p) :
    ShiftedGaussian(p){}

  double operator()(const double omega) {
    return 0.5*(Gaussian::operator()(omega-shift_) + Gaussian::operator()(omega+shift_));
  }
};

class LinearRiseExpDecay : public Model{
public:
  LinearRiseExpDecay(const alps::params &p): lambda_(p["LAMBDA"]){}
  double operator()(const double omega) {
    return lambda_*lambda_*omega*std::exp(-lambda_*omega);
  }

private:
  const double lambda_;
};

class QuadraticRiseExpDecay : public Model{
public:
  QuadraticRiseExpDecay(const alps::params &p): lambda_(p["LAMBDA"]){}
  double operator()(const double omega) {
    return (lambda_*lambda_*lambda_)/2.*(omega*omega)*std::exp(-lambda_*omega);
  }

private:
  const double lambda_;
};

class GeneralDoubleGaussian : public ShiftedGaussian
{
public:
  GeneralDoubleGaussian(const alps::params& p) :
    ShiftedGaussian(p), bnorm_(static_cast<double>(p["BOSE_NORM"])) {}

  double operator()(const double omega) {
    if (omega > 0)
      return Gaussian::operator()(omega);
    else
      return bnorm_*Gaussian::operator()(omega+shift_);
  }

private:
  const double bnorm_; 
};

///This class deals with tabulated model functions
class TabFunction : public Model
{
public:
  ///constructor will read in a default model from a file. File format is:
  /// First column: frequency
  /// Second column: value of default model
  /// anything after that: ignored.
  TabFunction(const alps::params& p, std::string const& name);

  ///return value of default model. If INSIDE interval we have data in: return linearly interpolated data. Otherwise: return zero.
  double operator()(const double omega);

private:
  ///private variable to store the frequency grid
  std::vector<double> Omega_;
  ///private variable to store the tabulated value of the default model at a frequency belonging to Omega_
  std::vector<double> Def_;
};



class GeneralDefaultModel : public DefaultModel
{
public:

  GeneralDefaultModel(const alps::params& p, boost::shared_ptr<Model> mod);

  ///given a number x between 0 and 1, find the frequency omega belonging to x.
  double omega(const double x) const;

  ///D returns the derivative of the integrated default model
  ///i.e. just the value of the model function at frequency omega
  double D(const double omega) const;

  ///returns the integrated default model
  double x(const double t) const;

private:
  boost::shared_ptr<Model> Mod;
  const int ntab;
  std::vector<double> xtab; //xtab has an equidistantly tabulated discretized model function

  double norm();
};



boost::shared_ptr<DefaultModel> make_default_model(const alps::params& parms, std::string const& name);
