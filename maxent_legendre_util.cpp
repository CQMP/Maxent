
#include "maxent_legendre_util.hpp"
#include <boost/math/special_functions/legendre.hpp> //needed for Legendre transform
#include <numeric>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random.hpp>

namespace bmth = boost::math;
Legendre_util::Legendre_util(const double T,const int ndat, int maxl, const vector_type y, const vector_type sigma):lmax_(-1),T_(T),ndat_(ndat),maxl_(maxl),y_(y),sigma_(sigma) {}

void Legendre_util::constructTauPoints(const alps::params &p){
    tau_points.resize(ndat_);
    if(p.defined("TAU_1")) //hack to see if we have imported tau points
        for(int j=0;j<ndat_;j++)
            tau_points[j]=p["TAU_"+boost::lexical_cast<std::string>(j)];
    else
        for(int j=0;j<ndat_;j++)
            tau_points[j] = j / ((ndat_)* T_); //TODO: standardize tau grid

}
vector_type Legendre_util::generateGaussNoise(vector_type data, vector_type err,boost::mt19937 &rng){
    
    typedef boost::variate_generator<boost::mt19937&,boost::normal_distribution<> > ran_gen;
    //notice the & in the first template argument and function rng argument.
    //If we omit this, it will compile and run
    //however, the numbers will be less(/not) random b/c it will copy the generator
    //each time, outputting the mean with some noise, rather than truly random
    
    const int N = data.size();
    vector_type data_noise(N);
    for(int i=0;i<N;i++){
        boost::normal_distribution<> s(data[i],err[i]);
        data_noise[i] = ran_gen(rng,s)();
    }
    return data_noise;
}
Legendre_util::return_type Legendre_util::bootstrap(double (Legendre_util::*f)(vector_type,void*),
                                                    const vector_type data,const vector_type err,
                                                    void *args, const int maxit){
    //bootstrap consists of generating noisy data within error bars
    //and determinging the variation on the output
    std::vector<double> newData(maxit);
    std::cout << std::setprecision(14);
    boost::mt19937 rng;
    rng.seed(static_cast<unsigned int>(std::time(0)));
    for(int i=0;i<maxit;i++){
         vector_type temp_data= generateGaussNoise(data, err,rng);
        //std::cout<<temp_data[0]<<" " << &temp_data<<std::endl;
        double point = (this->*f)(temp_data,args);
        newData[i]=point;
    }
    //calculate statistics
    std::vector<double>::iterator it=newData.begin();
    double mean=0.0;
    while(it!=newData.end()){
        mean+= *it;
        it++;
    }
    mean/=newData.size();

    it=newData.begin();
    double stdev = 0.0;
    while(it!=newData.end()){
        stdev+= (*it-mean)*(*it-mean);
        it++;
    }
    stdev/=(newData.size()-1);
    stdev=std::sqrt(stdev);
    double std_err = stdev/std::sqrt(newData.size());

    //return (data,err)
    return_type r;
    r.first = mean;
    r.second = std_err;
    return r;
}
///Gl=sqrt(2l+1) \int d\tau P_l(x(\tau))*G(\tau)
double Legendre_util::generateGl(vector_type gtau, int l){
    double gsum= 0.0;
    //TODO: implement higher order integration
    //      given good number of data
    
    for(int i=0;i<ndat_-1;i++){
        double dtau = tau_points[i+1]-tau_points[i];
        gsum+=dtau*(bmth::legendre_p(l, 2*tau_points[i]*T_-1)*gtau[i]);
    }
    //endpoint:
    double dtau=tau_points[ndat_-1]-tau_points[ndat_-2];
    gsum+=dtau*(bmth::legendre_p(l, 2*tau_points[ndat_-1]*T_-1)*gtau[ndat_-1]);
    
    return std::sqrt(2*l+1)*gsum;
}
double Legendre_util::generateGlBoot(vector_type gtau, void* arg){
    int l =  *(intptr_t *) arg;
    return Legendre_util::generateGl(gtau,l);
}

void Legendre_util::convertTauToGl(const alps::params &p){

    Gl.resize(maxl_+1);
    err_.resize(maxl_+1);
    constructTauPoints(p);
    
    while(lmax_<maxl_){
        lmax_++;

        //Gl[lmax_] = generateGl(y_,lmax_);//sqrt(2*lmax_+1)*I;
        double (Legendre_util::*f)(vector_type,void*);
        f = &Legendre_util::generateGlBoot;
        return_type Gl_tmp = bootstrap(f, y_,sigma_,&lmax_,500);
        std::cout  << Gl_tmp.first << " " << Gl_tmp.second << std::endl;
        Gl[lmax_]=Gl_tmp.first;
        err_[lmax_]=Gl_tmp.second;
        //TODO: check to see if Gl is within error bars
    }
    
    lmax_--;
    std::cout<<"Using " << lmax_ << " Legendre points" << std::endl;
    if(p["VERBOSE"]|false)
        std::cout << "With an error of:" <<std::abs(1) << std::endl;
    //return Gl;
}

void Legendre_util::reassignData(vector_type &y,vector_type &sigma, int &ndat,int &lmax, bool verbose){
    //switch data
    y.resize(lmax_);
    for(int i=0;i<lmax_;i++){
        y[i]=Gl[i];
        sigma[i]=err_[i];
    }
    if(verbose)
        std::cout<<"Gl points:"<<std::endl<<y<<std::endl;
    ndat=lmax_;
    lmax=lmax_;
}