
#include "maxent_legendre_util.hpp"
#include <boost/math/special_functions/legendre.hpp> //needed for Legendre transform
namespace bmth = boost::math;
legendre_util::legendre_util(const double T,const int ndat, int maxl, const vector_type y):T_(T),ndat_(ndat),maxl_(maxl),y_(y) {}

void legendre_util::constructTauPoints(const alps::params &p){
    tau_points.resize(ndat_);
    if(p.defined("TAU_1")) //hack to see if we have imported tau points
        for(int j=0;j<ndat_;j++)
            tau_points[j]=p["TAU_"+boost::lexical_cast<std::string>(j)];
    else
        for(int j=0;j<ndat_;j++)
            tau_points[j] = j / ((ndat_)* T_); //TODO: standardize tau grid

}

void legendre_util::convertTauToGl(const alps::params &p){

    Gl.resize(maxl_+1);
    constructTauPoints(p);
    
    double I,tau;
    double G0_lmax=0; //use a point to backcontinue to find cutoff
    double G0_prev=0;
    
    
    while(lmax_<maxl_){
        lmax_++;
        I=0;
        //int [0,beta] P_l(x(tau))*G(tau)
        for(int i=0;i<ndat_-1;i++){
            I+= bmth::legendre_p(lmax_, 2*tau_points[i]*T_-1)
            *y_[i]*(tau_points[i+1]-tau_points[i]);
        }
        I+= bmth::legendre_p(lmax_, 2*tau_points[ndat_-1]*T_-1)
        *y_[ndat_-1]*(tau_points[ndat_-1]-tau_points[ndat_-2]);
        Gl[lmax_] = sqrt(2*lmax_+1)*I;
        
        //after some l, subsequent points only contribute noise
        //this is essentially a convergence check
        G0_lmax=0;
        G0_prev=0;
        for(int l=0;l<lmax_;l++){
            G0_prev+= sqrt(2*l+1)*T_*bmth::legendre_p(l, 2*tau_points[ndat_/2]*T_-1)
            *Gl[l];
        }
        G0_lmax=G0_prev + sqrt(2*lmax_+1)*T_*bmth::legendre_p(lmax_, 2*tau_points[ndat_/2]*T_-1)
        *Gl[lmax_];
        //if(std::abs(y()[ndat()/2]-G0_lmax)>std::abs(y()[ndat()/2]-G0_prev) && 0.0001>std::abs(y()[ndat()/2]-G0_prev))
        //    break;
    }
    
    lmax_--;
    std::cout<<"Using " << lmax_ << " Legendre points" << std::endl;
    if(p["VERBOSE"]|false)
        std::cout << "With an error of:" <<std::abs(y_[0]-G0_prev) << std::endl;
    //return Gl;
}

void legendre_util::reassignData(vector_type &y, int &ndat,int &lmax, bool verbose){
    //switch data
    y.resize(lmax_);
    for(int i=0;i<lmax_;i++)
        y[i]=Gl[i];
    if(verbose)
        std::cout<<"Gl points:"<<std::endl<<y<<std::endl;
    ndat=lmax_;
    lmax=lmax_;
}