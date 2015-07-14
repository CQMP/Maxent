#pragma once

#include <alps/params.hpp>

/// write all default values and descriptions 
//  for a maxent param file
inline void set_defaults(alps::params& p){
	//---------------------------------
	//		General
	//---------------------------------
	p.define<int>("DATA_IN_HDF5",false,"DATA_IN_HDF5");
	p.define<int>("TEXT_OUTPUT",true,"TEXT_OUTPUT");
    	p.define<int>("ENFORCE_NORMALIZATION",false,"ENFORCE_NORMALIZATION");
	p.define<int>("VERBOSE",false,"VERBOSE");
	p.define<int>("SELF",false,"input is a self energy");
	p.define<int>("MAX_IT",1000,"Maximum Iterations for the fitting routine");
	p.define<int>("N_ALPHA",60,"N_ALPHA");
	p.define<double>("ALPHA_MIN",0.01,"ALPHA_MIN");
	p.define<double>("ALPHA_MAX",20,"ALPHA_MAX");
	p.define<double>("NORM",1.0,"NORM");
	//*********************************
	p.define<double>("BETA","BETA");
	p.define<int>("NDAT","# of input points");
	p.define<std::string>("DATA","","data file input");
	p.define<std::string>("BASENAME","","Specified output name");


	//---------------------------------
	//	    Default Model
	//---------------------------------
	p.define<double>("OMEGA_MAX",10,"OMEGA_MAX");
	p.define<double>("OMEGA_MIN","OMEGA_MIN");
	p.define<std::string>("DEFAULT_MODEL","flat","DEFAULT_MODEL");
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
	p.define<double>("CUT",0.01,"CUT");
	p.define<double>("SPREAD",4,"SPREAD");
	p.define<double>("LOG_MIN",1.0e-4,"LOG_MIN");
	p.define<int>("NFREQ",1000,"Number of A(omega) frequencies");
	p.define<std::string>("FREQUENCY_GRID","Lorentzian","Freq_Grid");

	//---------------------------------
	//		Kernel
	//---------------------------------
	p.define<std::string>("DATASPACE","time","DATASPACE");
	p.define<std::string>("KERNEL","ferminonic","KERNEL");	
	p.define<int>("PARTICLE_HOLE_SYMMETRY",false,"PH SYM"); 
	//*********************************

	//---------------------------------
	//		Legendre
	//---------------------------------
	p.define<int>("LEGENDRE",0,"LEGENDRE");
	p.define<int>("MAXL","Maximum L cutoff");
}
