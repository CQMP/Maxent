#pragma once

#include <alps/params.hpp>

/// write all default values and descriptions 
//  for a maxent param file
inline void set_defaults(alps::params& p){
	//---------------------------------
	//		General
	//---------------------------------
	p.define<int>("DATA_IN_HDF5",false,"1 if data is in HDF5 format");
	p.define<int>("TEXT_OUTPUT",true,"1 if results should be output to text files");
    	p.define<int>("ENFORCE_NORMALIZATION",false,"1 to renormalize last data point");
	p.define<int>("VERBOSE",false,"1 to print verbose output");
	p.define<int>("SELF",false,"input is a self energy");
	p.define<int>("MAX_IT",1000,"Maximum Iterations for the fitting routine");
	p.define<int>("N_ALPHA",60,"Number of alpha samples");
	p.define<double>("ALPHA_MIN",0.01,"Minimum alpha");
	p.define<double>("ALPHA_MAX",20,"Maximum alpha");
	p.define<double>("NORM",1.0,"NORM");
	//*********************************
	p.define<double>("BETA","beta, inverse temperature");
	p.define<int>("NDAT","# of input points");
	p.define<std::string>("DATA","","data file input");
	p.define<std::string>("BASENAME","","Specified output name \n(generated if not given)");
	p.define<int>("MODEL_RUNS","How many default model runs");
  p.define<double>("TAU_0","Used for input tau points");


	//---------------------------------
	//	    Default Model
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
	//		Grid
	//---------------------------------
	p.define<double>("CUT",0.01,"cut for lorentzian grids");
	p.define<double>("SPREAD",4,"spread for quadratic grid");
	p.define<double>("LOG_MIN",1.0e-4,"log_min for log grid");
	p.define<int>("NFREQ",1000,"Number of A(omega) frequencies");
	p.define<std::string>("FREQUENCY_GRID","Lorentzian","Type of frequency grid");

	//---------------------------------
	//		Kernel
	//---------------------------------
	p.define<std::string>("DATASPACE","time","Time or Frequency space");
	p.define<std::string>("KERNEL","fermionic","Type of kernel: \nFermionic,Bosonic,Boris,Legendre");	
	p.define<int>("PARTICLE_HOLE_SYMMETRY",false,"If particle hole symmetric"); 
	//*********************************

	//---------------------------------
	//		Legendre
	//---------------------------------
	p.define<int>("LEGENDRE",0,"LEGENDRE");
	p.define<int>("MAXL","Maximum L cutoff");
}
