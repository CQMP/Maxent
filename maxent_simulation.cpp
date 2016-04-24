/*
 * Copyright (C) 1998-2016 ALPS Collaboration. See COPYRIGHT.TXT
 * All rights reserved. Use is subject to license terms. See LICENSE.TXT
 * For use in publications, see ACKNOWLEDGE.TXT
 */

#include "maxent.hpp"
#include <alps/config.hpp> // needed to set up correct bindings
#include <alps/hdf5/vector.hpp>
#include <boost/math/special_functions/fpclassify.hpp> //needed for boost::math::isnan
#include <Eigen/LU>
#include "eigen_hdf5.hpp"
#include <iomanip>

MaxEntSimulation::MaxEntSimulation(alps::params &parms)
: MaxEntHelper(parms)
, alpha((int)parms["N_ALPHA"])              //This is the # of \alpha parameters that should be tried.
, norm(parms["NORM"])                                             //The integral is normalized to NORM (use e.g. for self-energies
, max_it(parms["MAX_IT"])                                       //The number of iterations done in the root finding procedure
, Kernel_type(parms["KERNEL"].as<std::string>())
, verbose(parms["VERBOSE"])
, text_output(parms["TEXT_OUTPUT"])
, self(parms["SELF"])
, make_back(parms["BACKCONTINUE"])
, gen_err(parms["GENERATE_ERR"])
, qvec((int)parms["N_ALPHA"])
, nfreq(parms["NFREQ"].as<int>())
{
  std::string bn=parms["BASENAME"]; name=bn+'.';

  if(norm != 1.) std::cerr<<"WARNING: Redefinition of parameter NORM: Input (and output) data are assumed to be normalized to NORM."<<std::endl;
  const double alpha_min = parms["ALPHA_MIN"];                                          //Smallest value of \alpha that is tried
  const double alpha_max = parms["ALPHA_MAX"];                                          //Largest  value of \alpha that is tried
  alpha[0] = alpha_max;
  for (std::size_t a=1; a<alpha.size(); ++a)                                            //These are all the alpa values on a log grid
    alpha[a] =  alpha[a-1] * std::pow(alpha_min/alpha_max, 1./double(alpha.size()-1));
}
///define parameter defaults
void MaxEntSimulation::define_parameters(alps::params &p){
  p.description("Maxent - a utility for " 
    "performing analytical continuation \n \t using the method of Maximum Entropy\n");
  //---------------------------------
  //    General
  //---------------------------------
  p.define<bool>("DATA_IN_HDF5",false,"1 if data is in HDF5 format");
  p.define<bool>("TEXT_OUTPUT",true,"1 if results should be output to text files");
  p.define<bool>("ENFORCE_NORMALIZATION",false,"1 to renormalize last data point");
  p.define<bool>("VERBOSE",false,"1 to print verbose output");
  p.define<bool>("SELF",false,"input is a self energy");
  p.define<int>("MAX_IT",1000,"Maximum Iterations for the fitting routine");
  p.define<int>("N_ALPHA",60,"Number of alpha samples");
  p.define<double>("ALPHA_MIN",0.01,"Minimum alpha");
  p.define<double>("ALPHA_MAX",20,"Maximum alpha");
  p.define<double>("NORM",1.0,"NORM");
  p.define<bool>("BACKCONTINUE",true,"Output A(omega) back to imaginary axis");
  p.define<bool>("GENERATE_ERR",false,"Generate a bootstrap approximation for error bars");
  //*********************************
  p.define<double>("BETA","beta, inverse temperature");
  p.define<int>("NDAT","# of input points");
  p.define<std::string>("DATA","","data file input");
  p.define<std::string>("BASENAME","","Specified output name \n(generated if not given)");
  p.define<int>("MODEL_RUNS","How many default model runs");
  p.define<double>("X_0","G input for param file entry");
  p.define<double>("TAU_0","Used for input tau points");


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
  p.define<double>("SIGMA","stddev - For Gaussian models");
  p.define<double>("GAMMA","width of Lorentzian model"); 
  p.define<double>("GAMMA2","for Two Lorentzian models");

  p.define<double>("LAMBDA","for ___ExpDecay models");
  p.define<double>("BOSE_NORM","General Double Gaussian Norm");
  

  //---------------------------------
  //    Grid
  //---------------------------------
  p.define<double>("CUT",0.01,"cut for lorentzian grids");
  p.define<double>("SPREAD",4,"spread for quadratic grid");
  p.define<double>("LOG_MIN",1.0e-4,"log_min for log grid");
  p.define<int>("NFREQ",1000,"Number of A(omega) frequencies");
  p.define<std::string>("FREQUENCY_GRID","Lorentzian","Type of frequency grid");

  //---------------------------------
  //    Kernel
  //---------------------------------
  p.define<std::string>("DATASPACE","time","Time or Frequency space");
  p.define<std::string>("KERNEL","fermionic","Type of kernel: \nFermionic,Bosonic,Boris,Legendre"); 
  p.define<bool>("PARTICLE_HOLE_SYMMETRY",false,"If particle hole symmetric"); 
}
void MaxEntSimulation::run()
{
  lprob.resize(alpha.size());
  chi_sq.resize(alpha.size());
  spectra.resize(alpha.size());
  u = transform_into_singular_space(Default());

  ofstream_ spectral_function_file;
  ofstream_ fits_file;

  if (text_output) {
    spectral_function_file.open((name+"spex.dat").c_str());
    fits_file.open((name+"fits.dat").c_str());
  }
  //this loop is the 'core' of the maxent program: iterate over all alphas, compute the spectra, normalization, and probabilities
  //loop over all alpha values
  for (std::size_t a=0; a<alpha.size(); ++a) {
    std::cerr << "alpha it: " << a << "\t";
    //fitting procedure for 'u'
    u = levenberg_marquardt(u, alpha[a]);
    //computation of spectral function out of 'u'
    vector_type A = get_spectrum(u);
    //computation of normalization
    std::cerr << "norm: " << transform_into_real_space(u).sum() << "\t";
    if (text_output) {
      spectral_function_file<<"# alpha: "<<alpha[a]<<std::endl;
      for (std::size_t i=0; i<A.size(); ++i)
        spectral_function_file << omega_coord(i) << " " << A[i] << "\n";
      spectral_function_file << "\n";
    }
    //computation of probability
    lprob[a] = log_prob(u, alpha[a]);
    spectra[a] = A;
    //computation of chi2
    double chi_squared = chi2(transform_into_real_space(u));
    chi_sq[a] = chi_squared;
    if (verbose) std::cerr << "0.5*chi2  : " << 0.5*chi_squared;
    std::cerr << std::endl;
    if (text_output) print_chi2(transform_into_real_space(u), fits_file);
    qvec(a)=Q(u,alpha[a]); 
  }
    omegaGrid.resize(nfreq);
    for(std::size_t i=0;i<nfreq;i++)
	    omegaGrid(i)=omega_coord(i);
}
  //everything from here on down is evaluation.
void MaxEntSimulation::evaluate(){
  if (text_output) {
    ofstream_ chi_squared_file;
    chi_squared_file.open((name+"chi2.dat").c_str());
    for (std::size_t a=0; a<chi_sq.size(); ++a){
      chi_squared_file << alpha[a] << " " << chi_sq[a] << std::endl;
    }
  }
  int a_chi = 0;
  double diff = std::abs(chi_sq[0]-ndat());
  for (std::size_t a=1; a<chi_sq.size(); ++a) {
    double diff_new = std::abs(chi_sq[a]-ndat());
    if (diff_new < diff) {
      diff = diff_new;
      a_chi = a;
    }
  }

  vector_type def = get_spectrum(transform_into_singular_space(Default()));
  if (text_output){
    ofstream_ chispec_file;
    chispec_file.open((name+"chispec.dat").c_str());
    for (std::size_t i=0; i<spectra[0].size(); ++i){
      chispec_file << omega_coord(i) << " " << spectra[a_chi][i]*norm << " " << def[i]*norm << std::endl;
    }
  }
  //boost::numeric::ublas::vector<double>::const_iterator max_lprob = std::max_element(lprob.begin(), lprob.end());  
  //const int max_a = max_lprob-lprob.begin();
  int max_a,nothing; double max_lprob;
  max_lprob=lprob.maxCoeff(&max_a,&nothing);
  const double factor = chi_scale_factor(spectra[max_a], chi_sq[max_a], alpha[max_a]);
  if (verbose) std::cerr << "chi scale factor: " << factor << std::endl;

	alps::hdf5::archive ar(name+"out.h5", alps::hdf5::archive::WRITE);
	ar << alps::make_pvp("/alpha/values",alpha);

  vector_type om(spectra[0].size());
  for (int i=0;i<om.size();i++) om[i] = omega_coord(i);        
  ar<<alps::make_pvp("/spectrum/omega",om);

  //output 'maximum' spectral function (classical maxent metod)
  if (text_output){
    ofstream_ maxspec_file;
    maxspec_file.open((name+"maxspec.dat").c_str());
    for (std::size_t i=0; i<spectra[0].size(); ++i)
      maxspec_file << omega_coord(i) << " " << spectra[max_a][i]*norm << " " << def[i]*norm << std::endl;
  }
	
  maxspec = spectra[max_a]*norm;
  vector_type specchi = spectra[a_chi]*norm;
	ar << alps::make_pvp("/spectrum/chi",specchi);
  ar << alps::make_pvp("/spectrum/maximum",maxspec);
  
	vector_type prob(lprob.size());
  for (std::size_t a=0; a<prob.size(); ++a) 
    prob[a] = exp(lprob[a]-max_lprob);
  double probnorm = 0;
  for (std::size_t a=0; a<prob.size()-1; ++a) 
    probnorm += 0.5*(prob[a]+prob[a+1])*(alpha[a]-alpha[a+1]);
  prob /= probnorm;
  ar << alps::make_pvp("/alpha/probability",prob);
  if (text_output){
    ofstream_ prob_str;
    prob_str.open((name+"prob.dat").c_str());
    for (std::size_t a=0; a<prob.size(); ++a) {
      prob_str << alpha[a] << "\t" << prob[a] << "\n";
    }
  }
  postprobdef = 0;
  for (std::size_t a=0; a<lprob.size()-1; ++a) 
    postprobdef += 0.5*(exp(lprob[a])+exp(lprob[a+1]))*(alpha[a]-alpha[a+1]);
  std::cout << "posterior probability of the default model: " << postprobdef << std::endl;

  //compute 'average' spectral function (Brian's method)
  avspec.resize(spectra[0].size());
  for (std::size_t i=0; i<avspec.size(); ++i) {
    avspec[i] = 0.;
    for (std::size_t a=0; a<prob.size()-1; ++a) 
      avspec[i] += 0.5*(prob[a]*spectra[a][i] +prob[a+1]*spectra[a+1][i])*(alpha[a]-alpha[a+1]);
  }
  //Estimate the variance for the spectrum
  vector_type varspec(spectra[0].size());
  for (std::size_t i=0; i<varspec.size(); ++i) {
    varspec[i] = 0.;
    for (std::size_t a=0; a<prob.size()-1; ++a)
      varspec[i] += 0.5*(prob[a]*(spectra[a][i]-avspec[i])*(spectra[a][i]-avspec[i]) + prob[a+1]*(spectra[a+1][i]-avspec[i])*(spectra[a+1][i]-avspec[i]))*(alpha[a]-alpha[a+1]);
  }
  avspec *= norm;
  varspec *= norm*norm;

  if (text_output){
    ofstream_ avspec_file;
    avspec_file.open((name+"avspec.dat").c_str());
    for (std::size_t  i=0; i<avspec.size(); ++i)
      avspec_file << omega_coord(i) << " " << avspec[i] << " " << def[i]*norm << std::endl;
  }
  ar << alps::make_pvp("/spectrum/average",avspec);
  ar << alps::make_pvp("/spectrum/variance",varspec);
  
  if(gen_err){
    //Bootstrap errors
    ofstream_ boot_file;
    boot_file.open((name+"booterr.dat").c_str());
    generateCovariantErr(maxspec,alpha[max_a],boot_file);
  }
  

  if(Kernel_type=="anomalous"){ //for the anomalous function: use A(omega)=Im Sigma(omega)/(pi omega).
    ofstream_ maxspec_anom_str;maxspec_anom_str.open((name+"maxspec_anom.dat").c_str());
    ofstream_ avspec_anom_str; avspec_anom_str.open((name+"avspec_anom.dat").c_str());
    vector_type spec(avspec.size());
    for (std::size_t  i=0; i<avspec.size(); ++i){ 
      //if(omega_coord(i)>=0.)
      spec[i] = avspec[i]*omega_coord(i)*M_PI;
      avspec_anom_str << omega_coord(i) << " " << avspec[i]*omega_coord(i)*M_PI<<std::endl;
    }
    ar << alps::make_pvp("/spectrum/anomalous/average",spec);
    for (std::size_t i=0; i<spectra[0].size(); ++i){
      //if(omega_coord(i)>=0.)
      spec[i] = spectra[max_a][i]*norm*omega_coord(i)*M_PI;
      maxspec_anom_str << omega_coord(i) << " " << spectra[max_a][i]*norm*omega_coord(i)*M_PI << std::endl;
    }
    ar << alps::make_pvp("/spectrum/anomalous/maximum",spec);
  }
  if(Kernel_type=="bosonic"){ //for the anomalous function: use A(Omega_)=Im chi(Omega_)/(pi Omega_) (as for anomalous)
    vector_type spec(avspec.size());
    for (std::size_t  i=0; i<avspec.size(); ++i){
      spec[i] = avspec[i]*omega_coord(i);
    }
    if (text_output) {
      ofstream_ avspec_anom_str;avspec_anom_str.open((name+"maxspec_bose.dat").c_str());
      for (std::size_t  i=0; i<avspec.size(); ++i){
        //if(omega_coord(i)>=0.)
        avspec_anom_str << omega_coord(i) << " " << spec[i]<<std::endl;
      }
    }
    ar << alps::make_pvp("/spectrum/bosonic/average",spec);
    for (std::size_t i=0; i<spectra[0].size(); ++i){
      //if(omega_coord(i)>=0.)
      spec[i] = spectra[max_a][i]*norm*omega_coord(i);
    }
    if (text_output) {
      ofstream_ maxspec_anom_str;maxspec_anom_str.open((name+"avspec_bose.dat").c_str());
      for (std::size_t i=0; i<spectra[0].size(); ++i){
        maxspec_anom_str << omega_coord(i) << " " << spectra[max_a][i]*norm*omega_coord(i) << std::endl;
      }
    }
    ar << alps::make_pvp("/spectrum/bosonic/maximum",spec);
  }

  //don't understand why this was commented out...
  if(self){
    // A quick word about normalization here. Usually we have G(iomega_n) = -1/pi \int_{-\infty}^\infty Im G(omega)/(omega_n - omega).
    // However, we are not interested in Im G but instead in A. In the case of the self-energy we have, analogously,
    // Sigma(i\omega_n) = -1/pi \int_{-\infty}^\infty Im \Sigma(omega)/(omega_n - omega); and we define A_\Sigma(omega) = -1/pi Sigma(omega). This makes
    // A_\Sigma be always positive, whereas Im Sigma(omega) is always negative.
    // here we compute Im Sigma out of A:
    //
    // for the self energy: use Im Sigma(omega)=-A(omega)*pi
    ofstream_ maxspec_self_str;maxspec_self_str.open((name+"maxspec_self.dat").c_str());
    ofstream_ avspec_self_str; avspec_self_str.open((name+"avspec_self.dat").c_str());
    for (std::size_t  i=0; i<avspec.size(); ++i){ 
      avspec_self_str << omega_coord(i) << " " << -avspec[i]*M_PI<< " " << -def[i]*norm*M_PI<<std::endl;
    }
    for (std::size_t i=0; i<spectra[0].size(); ++i){
      maxspec_self_str << omega_coord(i) << " " << -spectra[max_a][i]*norm*M_PI<< " " << -def[i]*norm*M_PI << std::endl;
    }
    //for public facing variables
    avspec*=-M_PI;
    maxspec*=-M_PI;
  }

  if(text_output && make_back){
    ofstream_ avspec_back_file,maxspec_back_file,chispec_back_file;
    avspec_back_file.open((name+"avspec_back.dat").c_str());
    maxspec_back_file.open((name+"maxspec_back.dat").c_str());
    chispec_back_file.open((name+"chispec_back.dat").c_str());

    const std::string sp = "    ";
    double norm_fix = norm;
    //fix self-energy normalization 
    if(self){
      norm_fix *=-M_PI;
    }

    std::cerr << "spectra"<<sp<< " max backcont diff" <<sp<<  "chi^2 value " <<std::endl;
    std::cerr << "======="<< sp<<" ================="<< sp<< "=========== " <<std::endl;
    backcontinue(chispec_back_file,specchi,norm,"chispec");
    backcontinue(avspec_back_file,avspec,norm_fix,"avspec ");
    backcontinue(maxspec_back_file,maxspec,norm_fix,"maxspec");
  }
  ar.close();
}



//this is the levenberg marquardt fitting procedure. It minimizes the quantity Q = 1/2 chi^2 - \alpha S
// 
vector_type MaxEntSimulation::levenberg_marquardt(vector_type u, const double alpha) const
{
  using namespace boost::numeric;
  double mu = 1e-18;
  const double nu = 1.3;
  double Q1=0.;
  int it = 0;
  int it2 = 0;
  for (; it<max_it; it++) {
    vector_type delta;
    if(boost::math::isnan(Q1))
        throw std::logic_error("Q=NaN, something went wrong");
    for (it2=0; it2<max_it; ++it2) {
      //compute change vector delta to u
      delta = iteration(u, alpha, mu);
      /*std::cout<<"delta is: "<<delta<<std::endl;
      vector_type z=transform_into_real_space(delta);
      for(int i=0;i<z.size();++i){
        std::cout<<omega_coord(i)<<" "<<z(i)<<std::endl;
      }*/
      //compute Q = 1/2 chi^2 - \alpha S
      Q1 = Q(u+delta, alpha);
      if (step_length(delta, u)<=0.02) {
        break;
      }
      else if (mu<1e20) {
        mu *= nu;
      }

    } 
    u += delta;
    if (convergence(u, alpha)<=1e-4)
      break;
  }
  if (it == max_it) std::cerr<<"WARNING: iteration reached max_it without converging, your minimizer is having problems. Please be careful!"<<std::endl;
  if (verbose) std::cerr <<"Iterations: " << it+1 << "\t";
  std::cerr << "Q = 0.5chi^2-\\alpha*entropy: " << Q1 << "\t";
  if (verbose) std::cerr << "entropy: "<<entropy(transform_into_real_space(u))<<"\talpha*entropy: "<<alpha*entropy(transform_into_real_space(u))<<"\t ";
  return u;
}



//this function computes the change delta to the vector 'u' 
//to be used in the Levenberg Marquardt fitting procedure
vector_type MaxEntSimulation::iteration(vector_type u, const double alpha, const double mu) const
{
  matrix_type M = left_side(u);
  for (std::size_t i=0; i<M.rows(); ++i) 
    M(i,i) += alpha + mu;
  vector_type b = right_side(u) + alpha*u;
  matrix_type B(b.size(),1);
  for (std::size_t i=0; i<M.rows(); ++i) 
    B(i,0) = -b[i];
   //bindings::lapack::gesv(M, ipiv, B);
  //NOTE: gesv uses LU decomp, but we can switch to a safe QR routine as well
  matrix_type Bp = M.lu().solve(B);
  //may need a transposeInPlace();
  return Bp;
}
