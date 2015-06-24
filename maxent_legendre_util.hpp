#pragma once
#include<alps/ngs/params.hpp>
#include"maxent_blas.hpp"

class legendre_util{
public:
	legendre_util(const double T,const int ndat, int maxl, const vector_type y);
	///convert G(tau) data into G(l)
	void convertTauToGl(const alps::params &p);
	///put Gl data into the appropriate external vector
	void reassignData(vector_type &y, int &ndat,int &lmax, bool verbose);
	///calculated lmax getter
	int lmax() const{return lmax_;}
protected:
	///set up vector that contains corresponding tau points to the G(tau) data points 
	void constructTauPoints(const alps::params &p);
	///vector of size ndat+1 that holds legendre coefficients
	vector_type Gl;
	///error of Gl based on bootstrap
	vector_type err;
	///legendre coefficient cutoff; also size of true Gl values
	int lmax_=-1;
private:
	///temperature = 1/beta
	const double T_;
	///number of input data points
	const int ndat_;
	///maximum coefficient cutoff supplied by user
	int maxl_;
	///vector of G(tau) points
	const vector_type y_;
	///corresponding tau points to the G(tau) data points 
	vector_type tau_points;
};