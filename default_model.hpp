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
    if(p.defined("BLOW_UP"))
      throw std::logic_error("Previous versions supported a parameter \'blowup\'. I've removed this from the code, I don't think it should exist");
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



class Model
{
public:
  virtual double operator()(const double omega)=0;
  virtual ~Model(){}
};



class Gaussian : public Model 
{
public:
  Gaussian(const alps::params& p) : sigma(static_cast<double>(p["SIGMA"])) {}

  virtual double operator()(const double omega) {
    return std::exp(-omega*omega/2./sigma/sigma)/sqrt(2*M_PI)/sigma;
  }

private:
  const double sigma;
};

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


class ShiftedGaussian : public Gaussian
{
public:
  //  ShiftedGaussian(const alps::Parameters& p) :
  ShiftedGaussian(const alps::params& p) :
    Gaussian(p), shift(static_cast<double>(p["SHIFT"])){}

  double operator()(const double omega) {
    return Gaussian::operator()(omega-shift);
  }

protected:
  const double shift;
};


class DoubleGaussian : public ShiftedGaussian
{
public:
  DoubleGaussian(const alps::params& p) :
    ShiftedGaussian(p){}

  double operator()(const double omega) {
    return 0.5*(Gaussian::operator()(omega-shift) + Gaussian::operator()(omega+shift));
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
    ShiftedGaussian(p), bnorm(static_cast<double>(p["BOSE_NORM"])) {}

  double operator()(const double omega) {
    if (omega > 0)
      return Gaussian::operator()(omega);
    else
      return bnorm*Gaussian::operator()(omega+shift);
  }

private:
  const double bnorm; 
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

  //return value of default model. If INSIDE interval we have data in: return linearly interpolated data. Otherwise: return zero.
  double operator()(const double omega);

private:
  std::vector<double> Omega;
  std::vector<double> Def;
};



class GeneralDefaultModel : public DefaultModel
{
public:

  GeneralDefaultModel(const alps::params& p, boost::shared_ptr<Model> mod);

  ///given a number x between 0 and 1, find the frequency omega belonging to x.
  double omega(const double x) const;
  /// returns the value of the model function at frequency omega
  double D(const double omega) const;

  //I have no idea what this does.
  double x(const double t) const;

private:
  boost::shared_ptr<Model> Mod;
  const int ntab;
  std::vector<double> xtab; //xtab has an equidistantly tabulated discretized model function

  double norm();
};



boost::shared_ptr<DefaultModel> make_default_model(const alps::params& parms, std::string const& name);
