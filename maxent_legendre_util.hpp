#pragma once
#include<alps/params.hpp>
#include"maxent_blas.hpp"
#include <boost/random/mersenne_twister.hpp>
class Legendre_util{
public:
	Legendre_util(const double T,const int ndat, int maxl, const vector_type y, const vector_type sigma);
	///convert G(tau) data into G(l)
	void convertTauToGl(const alps::params &p);
	///put Gl data into the appropriate external vector
	void reassignData(vector_type &y, vector_type &sigma, int &ndat,int &lmax, bool verbose);
	///calculated lmax getter
	int lmax() const{return lmax_;}
	//struct return_type {doublex data; vector_type err;};
    typedef std::pair<double,double> return_type;
	///boostrap a function f to determine error on said function
	return_type bootstrap(double (Legendre_util::*f)(vector_type,void*),const vector_type data, const vector_type err,
            				void *args, const int maxit);
	///check the high frequency (tail) behaivor on a Gl vector
	double checkTail(vector_type gl_in,int order);
protected:
	///set up vector that contains corresponding tau points to the G(tau) data points 
	void constructTauPoints(const alps::params &p);
	///determine Gl coefficient based on input
	double generateGl(vector_type gtau,int l);
	///bootstrap wrapper for generateGl
	double generateGlBoot(vector_type gtau, void* arg);
	///enforce high frequency behaivor
	double enforceTails(vector_type gl_in,int l, std::pair<double,double> tails);
	///bootstrap wrapper for enforceTails
	double enforceTailsBoot(vector_type gl_in, void* arg);
	///vector of size ndat+1 that holds legendre coefficients
	vector_type Gl;
	///error of Gl based on bootstrap
	vector_type err_;
	///legendre coefficient cutoff; also size of true Gl values
	int lmax_;
	
private:
	///temperature = 1/beta
	const double T_;
	///number of input data points
	const int ndat_;
	///maximum coefficient cutoff supplied by user
	int maxl_;
	///vector of G(tau) points
	const vector_type y_;
	///vector of error of G(tau) points
	const vector_type sigma_;
	///corresponding tau points to the G(tau) data points 
	vector_type tau_points;
	/// generate normal noise about each data[i] with mean=err[i]
	vector_type generateGaussNoise(vector_type data, vector_type err,boost::mt19937 &rng);
	///high frequency tail coefficients in legendre space;coefficient l, order(tail) p
	double tl(int l, int p);


};
