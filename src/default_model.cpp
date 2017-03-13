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

#include "default_model.hpp"
#include <fstream>
#include <boost/algorithm/string.hpp>



///This class deals with tabulated model functions
TabFunction::TabFunction(const alps::params& p, std::string const& name){
  std::string p_name = p[name].as<std::string>();
  std::ifstream defstream(p_name.c_str());
  if (!defstream)
    boost::throw_exception(std::invalid_argument("could not open default model file: "+p[name].as<std::string>()));
  double om, D;
  std::string line;
  while (getline(defstream, line)) {
    if(line.size()>0 && line[0]=='#') continue;
    std::stringstream linestream(line);
    linestream>>om>>D;
    Omega_.push_back(om);
    Def_.push_back(D);
  }
  double omega_max = p["OMEGA_MAX"];
  double omega_min = p.exists("OMEGA_MIN") ? p["OMEGA_MIN"] : -omega_max; //we had a 0 here in the bosonic case. That's not a good idea if you're continuing symmetric functions like chi(omega)/omega. Change omega_min to zero manually if you need it.
  if(Omega_[0]>omega_min || Omega_.back()<omega_max)
      boost::throw_exception(std::logic_error(std::logic_error("Input range outside of default model")));
    
  if (Omega_[0]!=omega_min || Omega_.back()!=omega_max){
    std::cout<<"Omega[ 0] "<<Omega_[0]<<" omega min: "<<omega_min<<std::endl;
    std::cout<<"Omega[-1] "<<Omega_.back()<<" omega max: "<<omega_max<<std::endl;
  }
}

//return value of default model. If INSIDE interval we have data in: return linearly interpolated data. Otherwise: return zero.
double TabFunction::operator()(const double omega) {
  //for values outside the grid point: return zero:
  if(omega < Omega_[0] || omega > Omega_.back())
    return 0.;
  //if we happen to exactly hit the first or last point
  if(omega == Omega_[0])
    return Def_[0];
  if(omega == Omega_.back())
    return Def_.back();

  //otherwise linear interpolation
  std::vector<double>::const_iterator ub = std::upper_bound(Omega_.begin(), Omega_.end(), omega);
  int index = ub - Omega_.begin();

  double om1 = Omega_[index-1];
  double om2 = Omega_[index];
  double D1 = Def_[index-1];
  double D2 = Def_[index];
  return -(D2-D1)/(om2-om1)*(om2-omega)+D2;
}


GeneralDefaultModel::GeneralDefaultModel(const alps::params& p, boost::shared_ptr<Model> mod)
: DefaultModel(p)
, Mod(mod)
, ntab(p["NTAB"])
, xtab(ntab) {
  tabulate_integral();
}

///given a number x between 0 and 1, find the frequency omega belonging to x.
double GeneralDefaultModel::omega(const double x) const {
  //range check for x
  if(x>1. || x<0.)
    throw std::logic_error("parameter x is out of bounds!");
  //binary search in ordered xtab
  std::vector<double>::const_iterator ub = std::upper_bound(xtab.begin(), xtab.end(), x);
  int omega_index = ub - xtab.begin();
  if (ub==xtab.end())
    omega_index = xtab.end() - xtab.begin() - 1;
  double om1 = omega_min + (omega_index-1)*(omega_max-omega_min)/(ntab-1);
  double om2 = omega_min + omega_index*(omega_max-omega_min)/(ntab-1);
  double x1 = xtab[omega_index-1];
  double x2 = xtab[omega_index];
  //linear interpolation
  return -(om2-om1)/(x2-x1)*(x2-x)+om2;
}

/// returns the value of the model function at frequency omega
double GeneralDefaultModel::D(const double omega) const {
  return (*Mod)(omega);
}

///x returns the integral of the default model, \f$x(t) = \int_0^t dt' D(\omega(t')) \f$
///This function performs a linear interpolation of the default model integral values given t between 0 and 1.
double GeneralDefaultModel::x(const double t) const {
  if(t>1. || t<0.) throw std::logic_error("parameter t is out of bounds!");
  int od = (int)(t*(ntab-1));
  if (od==(ntab-1))
    return 1.;
  double x1 = xtab[od];
  double x2 = xtab[od+1];
  return -(x2-x1)*(od+1-t*ntab)+x2;
}

void GeneralDefaultModel::tabulate_integral() {
  double sum = 0;
  xtab[0] = 0.;
  //this is an evaluation on an equidistant grid; sum integrated by trapezoidal rule
  double delta_omega = (omega_max - omega_min) / (ntab - 1);
  for (int o = 1; o < ntab; ++o) {
    double omega1 = omega_min + (o - 1) * delta_omega;
    double omega2 = omega_min + o * delta_omega;
    sum += ((*Mod)(omega1) + (*Mod)(omega2)) / 2. * delta_omega;
    xtab[o] = sum;
  }
  if(std::abs(sum -1)>1.e-7) std::cerr<<"careful: your default model does not integrate to one."<<std::endl;
  for (int o=0; o<ntab; ++o) {
    xtab[o] *= 1./sum;
  }
}

boost::shared_ptr<DefaultModel> make_default_model(const alps::params& parms, std::string const& name){
  std::string p_name = parms[name].as<std::string>();
  boost::to_lower(p_name);
  if (p_name == "flat") {
    std::cout << "Using flat default model" << std::endl;
    return boost::shared_ptr<DefaultModel>(new FlatDefaultModel(parms));
  }
  else if (p_name == "gaussian") {
    std::cout << "Using Gaussian default model" << std::endl;
    boost::shared_ptr<Model> Mod(new Gaussian(parms));
    return boost::shared_ptr<DefaultModel>(new GeneralDefaultModel(parms, Mod));
  }
  else if (p_name == "twogaussians" || p_name == "two gaussians") {
    std::cout << "Using sum of two Gaussians default model" << std::endl;
    boost::shared_ptr<Model> Mod(new TwoGaussians(parms));
    return boost::shared_ptr<DefaultModel>(new GeneralDefaultModel(parms, Mod));
  }
  else if (p_name == "shifted gaussian" || p_name == "shiftedgaussian") {
    std::cout << "Using shifted Gaussian default model" << std::endl;
    boost::shared_ptr<Model> Mod(new ShiftedGaussian(parms));
    return boost::shared_ptr<DefaultModel>(new GeneralDefaultModel(parms, Mod));
  }
  else if (p_name == "double gaussian" || p_name == "doublegaussian") {
    std::cout << "Using double Gaussian default model" << std::endl;
    boost::shared_ptr<Model> Mod(new DoubleGaussian(parms));
    return boost::shared_ptr<DefaultModel>(new GeneralDefaultModel(parms, Mod));
  }
  else if (p_name == "general double gaussian" || p_name == "generaldoublegaussian") {
    std::cout << "Using general double Gaussian default model" << std::endl;
    boost::shared_ptr<Model> Mod(new GeneralDoubleGaussian(parms));
    return boost::shared_ptr<DefaultModel>(new GeneralDefaultModel(parms, Mod));
  }
  else if (p_name == "linear rise exp decay" || p_name == "linearriseexpdecay") {
    std::cout << "Using linear rise exponential decay default model" << std::endl;
    boost::shared_ptr<Model> Mod(new LinearRiseExpDecay(parms));
    return boost::shared_ptr<DefaultModel>(new GeneralDefaultModel(parms, Mod));
  }
  else if (p_name == "quadratic rise exp decay" || p_name == "quadraticriseexpdecay") {
    std::cout << "Using quadratic rise exponential decay default model" << std::endl;
    boost::shared_ptr<Model> Mod(new QuadraticRiseExpDecay(parms));
    return boost::shared_ptr<DefaultModel>(new GeneralDefaultModel(parms, Mod));
  }
  else if (p_name == "lorentzian") {
    std::cout << "Using Lorentzian default model" << std::endl;
    boost::shared_ptr<Model> Mod(new Lorentzian(parms));
    return boost::shared_ptr<DefaultModel>(new GeneralDefaultModel(parms, Mod));
  }
  else if (p_name == "twolorentzians" || p_name == "two lorentzians") {
    std::cout << "Using sum of two Lorentzians default model" << std::endl;
    boost::shared_ptr<Model> Mod(new TwoLorentzians(parms));
    return boost::shared_ptr<DefaultModel>(new GeneralDefaultModel(parms, Mod));
  }
  else if (p_name == "shifted lorentzian" || p_name == "shiftedlorentzian") {
    std::cout << "Using shifted Lorentzian default model" << std::endl;
    boost::shared_ptr<Model> Mod(new ShiftedLorentzian(parms));
    return boost::shared_ptr<DefaultModel>(new GeneralDefaultModel(parms, Mod));
  }
  else if (p_name == "double lorentzian" || p_name == "doublelorentzian") {
    std::cout << "Using double Lorentzian default model" << std::endl;
    boost::shared_ptr<Model> Mod(new DoubleLorentzian(parms));
    return boost::shared_ptr<DefaultModel>(new GeneralDefaultModel(parms, Mod));
  }
  else { 
    std::cout << "Using tabulated default model" << std::endl;
    boost::shared_ptr<Model> Mod(new TabFunction(parms, name));
    return boost::shared_ptr<DefaultModel>(new GeneralDefaultModel(parms, Mod));
  }
}

///define parameter defaults
void DefaultModel::define_parameters(alps::params &p){
  //---------------------------------
  //      Default Model
  //---------------------------------
  p.define<double>("OMEGA_MAX",10,"Maximum frequency for A(omega) grid");
  p.define<double>("OMEGA_MIN","Minimum frequency, or =-OMEGA_MAX");
  p.define<std::string>("DEFAULT_MODEL","flat","Default model for entropy");
  p.define<double>("NORM1",0.5,"for Two Gaussians model");
  p.define<double>("SHIFT",0.0,"shift of a model");
  p.define<double>("SHIFT1",0.0,"for Two Gaussians model");
  p.define<double>("SHIFT2","for Two Gaussians model");
  p.define<double>("SIGMA1","for Two Gaussians model");
  p.define<double>("SIGMA2","for Two Gaussians model");
  p.define<double>("SIGMA","stddev - For Gaussian models");
  p.define<double>("GAMMA","width of Lorentzian model");
  p.define<double>("GAMMA1","for Two Lorentzian models");
  p.define<double>("GAMMA2","for Two Lorentzian models");

  p.define<double>("LAMBDA","for ___ExpDecay models");
  p.define<double>("BOSE_NORM","General Double Gaussian Norm");
  p.define<int>("NTAB",5001,"Discretization points of the Default model");

}
void DefaultModel::print_help(){
  std::cout << "Default model help - choices for generated default models"<< std::endl;
  std::cout << "For more information see examples/default_models.pdf\n"   << std::endl;
  std::cout << std::left << std::setw(23)<< "Model Name"           <<'\t' << "option=default" << std::endl;
  std::cout << std::left << std::setw(23)<< "=========="           <<'\t' << "==============" << std::endl;
  std::cout << std::left << std::setw(23)<< "flat"                 <<'\t' << "---" << "\n"
            << "------------------" << "\n"
            << std::left << std::setw(23)<< "gaussian"             <<'\t' << "SIGMA" << "\n"
            << std::left << std::setw(23)<< "double gaussian"      <<'\t' << "SIGMA;SHIFT=0" << "\n"
            << std::left << std::setw(23)<< "two gaussians"        <<'\t' << "SIGMA1;SIGMA2;SHIFT2;NORM1=0.5;SHIFT1=0" << "\n"
            << std::left << std::setw(23)<< "shifted gaussians"    <<'\t' << "SIGMA;SHIFT=0"  << "\n"
            << "------------------" << "\n"
            << std::left << std::setw(23)<< "lorentzian"           <<'\t' << "GAMMA" << "\n"
            << std::left << std::setw(23)<< "double lorentzian"    <<'\t' << "GAMMA;SHIFT=0" << "\n"
            << std::left << std::setw(23)<< "two lorentzians"      <<'\t' << "GAMMA1;GAMMA2;SHIFT2;SHIFT1=0" << "\n"
            << std::left << std::setw(23)<< "shifted lorentzian"   <<'\t' << "GAMMA;SHIFT=0" << "\n"
            << "------------------" << "\n"
            << std::left << std::setw(23)<< "LinearRiseExpDecay"   <<'\t' << "LAMBDA" << "\n"
            << std::left << std::setw(23)<< "QuadraticRiseExpDecay"<<'\t' << "LAMBDA" << std::endl;

}
