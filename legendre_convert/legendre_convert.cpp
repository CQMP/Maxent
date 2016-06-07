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

#include <iostream>
#include <fstream>
#include <boost/program_options.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/legendre.hpp> //needed for Legendre transform
#include <boost/math/special_functions/bessel.hpp> 

typedef std::vector<double> vector_type;
typedef std::pair<double,double> return_type;
typedef std::complex<double> Complex ;
namespace bmth = boost::math;


///E8 of Boehnke, et al
///high frequency expansion of Gl (see Tl)
double tl(int l, int p){
    //NOTE: p is tail number (order), starting at 1
    if((l+p)%2==0)
        return 0;
    else{
        
    //(l+p-1)!/(l-p-1)! => product(l+q for q in range(-p+2,p))
        double qsum =1;
        for(int q=-p+2;q<p;q++)
            qsum*=l+q;
        return std::pow(-1.0,p)*2*std::sqrt(2*l+1)*qsum/bmth::factorial<double>(p-1);
    }
}
///G(i\omega)=\sumT_{n\ell}G_\ell
//i.e. the unitary transformation
Complex Tl(int n, int l){
  //(5) of Boehnke, et al
  const Complex CONE(0,1);
  return std::pow(-1.0,n)*std::pow(CONE,l+1)*std::sqrt(2*l+1)*bmth::sph_bessel(l,(2*n+1)*M_PI/2);
}
///fix a set of Legendre coefficients to 2 high frequency tails
double enforceTails(double beta, const vector_type &gl_in,int l, return_type tails){
    //TODO: generalize pair->vector for higher order tails?
    //(10) of Boehnke, et al
    //Gl=>Gl+correction=(beta^n*c_n-sum(tl(n)*gl))*tl(n)/sum(tl(n)^2)
    
    double lmax = gl_in.size();
    //tails = [c1,c2]
    double c1 = tails.first;
    double c2 = tails.second;
    double correction=0.0;
    
    double proj1sum = 0.0,proj2sum=0.0;
    double tailsum1=0.0,tailsum2=0.0;
    for(int lp=0;lp<lmax;lp++){
        proj1sum+=tl(lp,1)*gl_in[lp];
        proj2sum+=tl(lp,2)*gl_in[lp];
        
        tailsum1+=tl(lp,1)*tl(lp,1);
        tailsum2+=tl(lp,2)*tl(lp,2);
    }
    correction+=(c1*beta-proj1sum)*tl(l,1)/tailsum1;
    correction+=(c2*beta*beta-proj2sum)*tl(l,2)/tailsum2;
    return correction;
}
///Given G(tau) and errors, generate G(l)
double generateGl(const vector_type &gtau,const vector_type &tau_points, int l,double beta){
	int ndat_ = gtau.size();
    double gsum= 0.0;
    //TODO: non-equispaced grid
    if((ndat_-1)%2==0){
    	//even number so we can use Simpson's rule

    	--ndat_;
    	double dtau = tau_points[1]-tau_points[0];
    	gsum+=gtau[0]*bmth::legendre_p(l,2*tau_points[0]/beta-1.0);
	    for(int i=1;i<ndat_/2;i++)
	    	gsum+=2*gtau[2*i]*bmth::legendre_p(l,2*tau_points[2*i]/beta-1);
	    for(int i=1;i<ndat_/2+1;i++)
	    	gsum+=4*gtau[2*i-1]*bmth::legendre_p(l,2*tau_points[2*i-1]/beta-1);

	    gsum+=gtau[ndat_]*bmth::legendre_p(l,2*tau_points[ndat_]/beta-1.0);
	    gsum*= dtau/3;
    }
    else{
	    
		//Riemann sums
	    /*for(int i=0;i<ndat_-1;i++){
	        double dtau = tau_points[i+1]-tau_points[i];
	        gsum+=dtau*(bmth::legendre_p(l, 2*tau_points[i]/beta-1)*gtau[i]);
	    }
	    //endpoint:
	    double dtau=tau_points[ndat_-1]-tau_points[ndat_-2];
	    gsum+=dtau*(bmth::legendre_p(l, 2*tau_points[ndat_-1]/beta-1)*gtau[ndat_-1]);*/

	    //trapezoidal rule
    /*    double dtau = tau_points[1]-tau_points[0];
	    gsum+=dtau*(bmth::legendre_p(l, 2*tau_points[0]/beta-1)*gtau[0]);
        
        for(int i=1;i<ndat_-1;i++){
        	 double dtau = tau_points[i+1]-tau_points[i];
	        gsum+=2*dtau*(bmth::legendre_p(l, 2*tau_points[i]/beta-1)*gtau[i]);
        }
     	dtau = tau_points[ndat_-1]-tau_points[ndat_-2];
        gsum+=dtau*(bmth::legendre_p(l, 2*tau_points[ndat_-1]/beta-1)*gtau[ndat_-1]);
        gsum/=2; */
        //try Simpson's rule within ndat-1 region, then trapezoidal for endpoint

        //if we have an odd number of points, we want simpsons to average
        //going between  1,N-1 and 0,N-2 (because [N] is undefined)

        double gsum1=0.0,gsum2=0.0;
        ndat_-=2;
        double dtau = tau_points[1]-tau_points[0];
    	gsum1+=gtau[0]*bmth::legendre_p(l,2*tau_points[0]/beta-1.0);
	    for(int i=1;i<ndat_/2;i++)
	    	gsum1+=2*gtau[2*i]*bmth::legendre_p(l,2*tau_points[2*i]/beta-1);
	    for(int i=1;i<ndat_/2+1;i++)
	    	gsum1+=4*gtau[2*i-1]*bmth::legendre_p(l,2*tau_points[2*i-1]/beta-1);

	    gsum1+=gtau[ndat_]*bmth::legendre_p(l,2*tau_points[ndat_]/beta-1.0);
	    gsum1*= dtau/3;
	    
      //-----------------------
      ndat_++;
      dtau = tau_points[1]-tau_points[0];
      gsum2+=gtau[0]*bmth::legendre_p(l,2*tau_points[0]/beta-1.0);
      for(int i=2;i<ndat_/2;i++)
        gsum2+=2*gtau[2*i]*bmth::legendre_p(l,2*tau_points[2*i]/beta-1);
      for(int i=2;i<ndat_/2+1;i++)
        gsum2+=4*gtau[2*i-1]*bmth::legendre_p(l,2*tau_points[2*i-1]/beta-1);

      gsum2+=gtau[ndat_]*bmth::legendre_p(l,2*tau_points[ndat_]/beta-1.0);
      gsum1*= dtau/3;
      //----------------------
      gsum = (gsum1+gsum2)/2;

    }
    return std::sqrt(2*l+1)*gsum;
}
///given Legendre coefficients, generate a G(tau) point
double Gt(const vector_type &gl_in, double tau, double beta ){
	//sum l>=0 sqrt(2l+1)/beta*P_l(x(\tau))G_l

    double gsum = 0.0;
    int N = gl_in.size();
    for(int l=0;l<N;l++){
    	gsum+=std::sqrt(2*l+1)/beta*bmth::legendre_p(l,2.0*tau/beta-1.0)*gl_in[l];
    }
    return gsum;
}
///checking the high frequency tails
double checkTail(const vector_type &gl_in, int order,double beta){
    //(8) of Boehnke, et al
    // cp=1/beta^p * sum(tl(p)*gl)
    double lmax=gl_in.size();
    double tailsum=0.0;
    for(int i=0;i<lmax;i++){
        tailsum+=tl(i,order)*gl_in[i];
    }
    return tailsum/std::pow(beta,order);
}
vector_type plot_convergence(const vector_type &gtau,const vector_type &tau_points, const double beta,const double maxl,
						return_type tails,int taupoint, const bool tailEnforced){
	vector_type errpoints;
    for(int l=5;l<maxl;l++){
    	vector_type gl_temp;
    	gl_temp.clear();
    	for(int i=0;i<l;i++)
    		gl_temp.push_back(generateGl(gtau,tau_points,i,beta));
    	//************************
    	if(tailEnforced){
    		vector_type correction(l);
	    	for(int lp=0;lp<l;lp++)
	  			correction[lp]= enforceTails(beta, gl_temp,lp,tails);
	  		for(int lp=0;lp<l;lp++)
	  			gl_temp[lp]+=correction[lp];
    	}
    	//************************
    	errpoints.push_back(gtau[taupoint]-Gt(gl_temp,tau_points[taupoint],beta));
    	//errpoints.push_back(Gt(gl_temp,tau_points[taupoint],beta));
    }
    return errpoints;
}
///plot the convergence of matsubara points as a function lmax cutoff
vector_type plot_convergence_matsubara(const vector_type &gtau,const vector_type &tau_points, const double beta,const double maxl,
            return_type tails,int n, const bool tailEnforced){
  vector_type errpoints;
    for(int l=5;l<maxl;l++){
      vector_type gl_temp;
      gl_temp.clear();
      for(int i=0;i<l;i++)
        gl_temp.push_back(generateGl(gtau,tau_points,i,beta));
      //************************
      if(tailEnforced){
        vector_type correction(l);
        for(int lp=0;lp<l;lp++)
          correction[lp]= enforceTails(beta, gl_temp,lp,tails);
        for(int lp=0;lp<l;lp++)
          gl_temp[lp]+=correction[lp];
      }
      //************************
      Complex Gm = 0.0;
      for(int lp=0;lp<l;lp++){
        Gm+=Tl(n,lp)*gl_temp[lp];
      }
      errpoints.push_back(std::imag(Gm));
      //errpoints.push_back(Gt(gl_temp,tau_points[taupoint],beta));
    }
    return errpoints;
}

///iterate through Gl points to double check that we're not outside error bars. 
///resize Gl to a new cutoff that is within error bars
void checkErrorBars(vector_type &gl,vector_type &error,int &lmax,int step,bool enforceErrCheck){
      for(int l=0;l<lmax;l+=step){
      if(std::abs(error[l]/gl[l])>=1){
        std::cout<< "Warning! l="<< l<<" has G(l) on the order of errorbars!" << std::endl;
        if(enforceErrCheck){
          std::cout<<"resizing to lmax=" <<l-1<<std::endl;
          lmax=l-1;
          gl.resize(lmax);
          error.resize(lmax);
        }
      }
    }
}
//************************
/// struct for bootstrap routine
struct argstruct{
        double beta;int l;double tau; 
        std::pair<double,double> tails; vector_type *tau_points;
 };
///genereate a normally distrubted noisy vector. 
//i.e. output[i] = normally dist number with mean=data[i] and stddev=err[i]
vector_type generateGaussNoise(vector_type data, vector_type err,boost::mt19937 &rng){
    
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
///Generalized bootstrap routine. Requires the non-linear function to be
// f(vector_type v,void *arg) where *arg is most easily a struct
return_type bootstrap(double (*f)(vector_type,void*),
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
        double point = (*f)(temp_data,args);
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
///wrapper for generateGl in a bootstrap routine
double generateGlBoot(vector_type gtau,void* arg){
	argstruct *a = (argstruct *) arg;
    //int l =  *(intptr_t *) arg;
    vector_type tau_points = *(a->tau_points);
    return generateGl(gtau,tau_points,a->l,a->beta);
}
///wrapper for enforceTails in a bootstrap routine
double enforceTailsBoot(vector_type gl_in, void *arg){

    argstruct *a = (argstruct *) arg;
    //double beta, const vector_type &gl_in,int l, return_type tails
    return enforceTails(a->beta,gl_in,a->l,a->tails);
}
///wrapper for Gt in a bootstrap routine 
double GtBoot(vector_type gl_in, void *arg){
	//Gt(const vector_type &gl_in, double tau, double beta )
	argstruct *a = (argstruct *) arg;
	return Gt(gl_in,a->tau,a->beta);
}
//************************
int main(int argc, char**argv){
	namespace po = boost::program_options;
	double beta,c1,c2,c3;
	int lmax,numConvergence,maxit,maxn;
	bool backcontinue = false,enforceTailOff = false,enforceErrCheck=true,continueMatsubara=false;
	std::string input_gtau_filename,gl_filename;
	po::options_description desc("Allowed options");
	desc.add_options()
	("help", "show this help")
	("beta", po::value<double>(&beta), "inverse temperature")
	("input_gtau_file", po::value<std::string>(&input_gtau_filename)->default_value("Gtau.dat"), "input G(tau). Format: tau G(tau) err")
	("output_gl_file",po::value<std::string>(&gl_filename)->default_value("Gl.dat"),"Output G(l) file. Format: l Gl err")
	("tail1",po::value<double>(&c1)->default_value(1),"c1/iwn high frequency tail")
	("tail2",po::value<double>(&c2)->default_value(0),"c2/(iwn)^2 high frequency tail")
	//("tail3",po::value<double>(&c3),"c3/(iwn)^3 high frequnecy tail. Optional input NOT IMPLEMENTED")
	("lmax",po::value<int>(&lmax)->default_value(40),"maximum Legendre coefficient cutoff")
	("maxit",po::value<int>(&maxit)->default_value(100),"maximum bootstrap iterations")
	("backcontinue","Output a G(tau) conversion of the computed Gl with error bars. Output file Gtau_back.dat")
	("notail","don't enforce tails. Works for both plot_convergence and main GL output")
	("noerr","only warn when exceeding error bars; outputs full Gl from 0 to lmax")
	("plot_convergence",po::value<int>(&numConvergence),
		"output lmax vs G(0,beta,beta/2,beta/8). \narg = maximum l cutoff; if negative then don't enforce tails. note: w/o bootstrap")
	 ("maxn",po::value<int>(&maxn),"if used with plot_convergence, will plot maxn matsubara point convergences")
	;

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);
	if (vm.count("help")) {
	std::cout<<desc;
	return 1;
	}
	if(vm.count("backcontinue"))
		backcontinue = true;
	if(!vm.count("plot_convergence"))
		numConvergence =0;
  if(vm.count("maxn"))
    continueMatsubara=true;
	if(vm.count("notail"))
		enforceTailOff = true;
  if(vm.count("noerr"))
    enforceErrCheck = false;
	std::ifstream orig_file(input_gtau_filename.c_str());
  	if(!orig_file.good()) throw std::invalid_argument("Gtau file: "+input_gtau_filename+" could not be opened. specify with --input_gtau_file");
  	if(!vm.count("beta")) throw std::runtime_error("you need to specify the inverse temperature with --beta.");

  	vector_type tarr,Gtau,sigma;
  	do{
  		double tau, Gt,sig;
  		orig_file>>tau>>Gt>>sig>>std::ws;
  		tarr.push_back(tau);
  		Gtau.push_back(Gt);
  		sigma.push_back(sig);
  	} while(orig_file.good());
  	//************************
    //basic error checks
    if(tarr.back()>beta|| tarr.back()*1.5<beta){ 
        std::cout<< "Beta != tau[-1] beta: "<<beta<<" tau[-1]:"<<tarr.back() <<std::endl;
        throw std::logic_error("Beta != tau[-1] beta");
      }
    if(tarr.back() != beta && tarr[0] != 0){
      std::cout << "Please   try to supply one endpoint, tau=0 or tau=beta" << std::endl;
      throw std::logic_error("Neither endpoint was supplied");
    }
    //************************
  	return_type tails;
  	tails.first = c1; tails.second=c2;
  	int N = tarr.size(); 	//ndat aka size of data
  	std::cout<<"Found " << N << " G(tau) points" <<std::endl;
  	vector_type Gl(lmax), err(lmax);
  	std::cout << std::setprecision(5); //for formatting
  	//************************
  	  	//see if we can add a point to make Simpon's rule work
  	if ((N-1)%2==1){//((N-1)%2==1)
  		if(tarr.back()!=beta && tarr[0] == 0){
  			std::cout<<"Adding endpoint G(beta) for integration"<<std::endl;
  			tarr.push_back(beta);
  			Gtau.push_back(-c1-Gtau[0]);
  			sigma.push_back(sigma.back());
  			std::cout << std::setw(5)<< std::left<< "tau" << " " <<std::setw(12)<< std::left<<"G(tau=beta)" << " " << "sigma" << std::endl;
  			std::cout<< std::setw(5)<< std::left << tarr.back() << " "<<std::setw(12)<<std::left<<Gtau.back() << " "<<sigma.back() << std::endl;
  			N++;
  		}
  		else if(tarr[0] != 0 && tarr.back() == beta){
  			std::cout<<"Adding endpoint G(0) for integration"<<std::endl;
  			tarr.insert(tarr.begin(),0);
  			Gtau.insert(Gtau.begin(),-c1-Gtau[0]);
  			sigma.insert(sigma.begin(),sigma[0]);
  			std::cout << std::setw(5)<< std::left<< "tau" << " " <<std::setw(12)<< std::left<<"G(tau=0)" << " " << "sigma" << std::endl;
  			std::cout<< std::setw(5)<< std::left << tarr[0] << " "<<std::setw(12)<<std::left<<Gtau[0] << " "<<sigma[0] << std::endl;
  			N++;
  		}
    }
  	//************************
  	//check endpoints
  	double dtau_check = tarr[N-1]-tarr[N-2];
  	if(std::abs(beta-tarr[N-1])>dtau_check)
  		std::cout<<"Warning! End point G(tau=beta)="<<tarr[N-1]<<" is not within dtau("<<dtau_check<<") of beta"<< std::endl;
  	std::cout<<"Beta = "<<beta<<std::endl;
  	//double check tail input
  	std::cout << std::setw(12)<< std::left<< "Tail check:"<< std::setw(7)<< std::left<<"input "<< std::setw(5)<< std::left<<"measured" << std::endl;
  	double c1_test = -Gtau[0]-Gtau.back();
  	std::cout <<std::setw(12)<< std::left<<""<<std::setw(7)<< std::left<<c1<<std::setw(5)<< std::left<<c1_test<<std::endl;
  	//use forward difference derivative 
  	double Gp0 = (Gtau[1]-Gtau[0])/(tarr[1]-tarr[0]);
  	double Gpbeta = (Gtau[N-1]-Gtau[N-2])/(tarr[N-1]-tarr[N-2]);
  	std::cout <<std::setw(12)<< std::left<<""<<std::setw(7)<< std::left<<c2<<std::setw(5)<< std::left<<Gp0+Gpbeta<<std::endl;
  	//************************
  	std::cout << "Generating Gl values" << std::endl;
  	std::cout <<"Using a bootstrap with " << maxit << " iterations"<<std::endl;
    argstruct a; a.beta=beta;a.tails=tails;a.tau_points = &tarr;
  	for(int l=0;l<lmax;l++){
  		a.l=l;
  		return_type Gl_tmp = bootstrap(&generateGlBoot, Gtau,sigma,&a,maxit);
        Gl[l]=Gl_tmp.first;
       /*if(l%2==1 && c1==1 && c2==0 ){
        Gl[l]=0;
      }*/
        err[l]=Gl_tmp.second;
  	//Gl[l]=generateGl(Gtau,tarr,l,beta);
  	}
   	 //check if data is within errorbars!
    int step; //if PH then odd ls are =0 within error bars
    (std::abs(c2)<1e-6) ? step=2 : step=1;
    checkErrorBars(Gl,err,lmax,step,enforceErrCheck);
   	if(!enforceTailOff){
	   	//enforce tails
	   	//it is important to keep the correction seperate from the Gl vals while calculating
	   	std::cout <<"Enforcing Tails" << std::endl;
	  	vector_type correction(lmax), corr_err(lmax);
	  	for(int l=0;l<lmax;l++){
	  		a.l=l;
	  		return_type corr_tmp = bootstrap(&enforceTailsBoot,Gl,err,&a,maxit);
	  		correction[l]=corr_tmp.first;
	  		corr_err[l] = corr_tmp.second;
	  		//correction[l]=enforceTails(beta, Gl,l,tails);
	  	}
	  	for(int l=0;l<lmax;l++){
	  		Gl[l]+=correction[l];
	  		err[l]=std::sqrt(corr_err[l]*corr_err[l]+err[l]*err[l]);
	  	}
	}
    //************************
    //check if data is within errorbars!
    //TODO:not sure if we should do this twice or not
    checkErrorBars(Gl,err,lmax,step,enforceErrCheck);
    //************************
    //write to file
    std::ofstream gl_file(gl_filename.c_str());
    gl_file.precision(14); //keep output from truncation
    for(int l=0;l<lmax;l++){
	    gl_file << l<<" "<<Gl[l]<<" " <<err[l]<<std::endl;
    }

    std::cout << std::setw(15)<< std::left<< "Gl Tail check: "<< std::setw(7)<< std::left
	      <<"input "<< std::setw(5)<< std::left<<"measured" << std::endl;
    std::cout <<std::setw(15)<< std::left<<""<<std::setw(7)<< std::left<<c1<<std::setw(5)<< std::left<<checkTail(Gl,1,beta)<<std::endl;
    std::cout <<std::setw(15)<< std::left<<""<<std::setw(7)<< std::left<<c2<<std::setw(5)<< std::left<<checkTail(Gl,2,beta)<<std::endl;
    //************************
    if(backcontinue){
	std::cout << "Continuing back to time axis" <<std::endl;
	vector_type gtau_back(N);
	vector_type gtau_back_err(N);
	std::ofstream gt_back_file("Gtau_back.dat");
	for(int t=0;t<N;t++){
		//gtau_back[t] = Gt(Gl,tarr[t],beta);
		a.tau = tarr[t];
		return_type gtau_tmp = bootstrap(GtBoot,Gl,err,&a,maxit);
		gtau_back[t] = gtau_tmp.first;
		gtau_back_err[t] = gtau_tmp.second;
		gt_back_file << tarr[t] <<" "<<gtau_back[t]<<" " << gtau_back_err[t]<<std::endl;
	}

	//************************
	vector_type back_err(N),back_errbars(N); //a convenient way to visualize if the error bars are within the value
	for(int t=0;t<N;t++){
		back_err[t] = Gtau[t]-gtau_back[t];
		back_errbars[t] = std::sqrt(sigma[t]*sigma[t]+gtau_back_err[t]*gtau_back_err[t]);
	}
	std::ofstream gt_err_file("Gtau_back_errors.dat");
	gt_err_file << std::setprecision(14);
	gt_err_file << "#errors of the backcontinuation" << std::endl;
	int count=0;
	for(int i=0;i<N;i++)
		if(std::abs(back_err[i])<back_errbars[i])
			count++;
	int percent = (double)count/N *1000;
	std::cout<<"Points within errorbars: " << count<< " or "<< (double) percent/10<< "%" << std::endl;
	gt_err_file << "#points within errorbars: " << count << " or "<< (double)count/N<< "%" <<std::endl;
	gt_err_file << "#tau gtau_input-gtau_back err" << std::endl;
	for(int i=0;i<N;i++)
		gt_err_file<< tarr[i] << " " << back_err[i] << " "<< back_errbars[i] << std::endl;
    }
    //************************

    if(numConvergence !=0){
	std::ofstream convergence_file("Gl_convergence.dat");
	convergence_file <<std::setprecision(14);
	convergence_file << "#convergence of points. lmax vs err of point"<<std::endl;
	std::vector<vector_type> points;
	int t_points[4] = {0,N/2,N/8,N-1};
	if(numConvergence>0 && enforceTailOff == false && !continueMatsubara){
		std::cout << "Creating convergence file with enforcing tails" <<std::endl;
		convergence_file << "#enforcing tails: " << tails.first << " "<< tails.second<<std::endl;
		for(int i=0;i<4;i++) 
			points.push_back(plot_convergence(Gtau,tarr,beta,numConvergence,
							 tails,t_points[i],true));

	}
	else if(!continueMatsubara){
		std::cout << "Creating convergence file without enforcing tails" <<std::endl;
		convergence_file << "#w/0 enforcing tails" <<std::endl;
		numConvergence*=-1;
		for(int i=0;i<4;i++)
			points.push_back(plot_convergence(Gtau,tarr,beta,numConvergence,
							 tails,t_points[i],false));
	}
      else if(continueMatsubara){
        std::cout << "Creating convergence file with " << maxn << " matsubara points" << std::endl;
        numConvergence= std::abs(numConvergence);
        for(int n=0;n<maxn;n++)
          points.push_back(plot_convergence_matsubara(Gtau,tarr,beta,numConvergence,tails,n,true));
        points.push_back(plot_convergence_matsubara(Gtau,tarr,beta,numConvergence,tails,1000,true));
      }
      //************************
      //Look at the convergence of Matsubara points rather than tau points
      if(!continueMatsubara){
        convergence_file << "#lmax G(0) G(beta/2) G(beta/8) G(beta)" <<std::endl;
    		for(int i=0;i<numConvergence-5;i++){
    			convergence_file << i+5<< " ";
    			for(int p=0;p<4;p++)
    				convergence_file<< points[p][i] << " ";
    			convergence_file<<std::endl;
    		}
      }
      else{
        convergence_file << "#matsubara points. Last point is n=1000 for tail check" <<std::endl;
        for(int i=0;i<numConvergence-5;i++){
	    convergence_file << i+5<< " ";
	    for(int p=0;p<maxn+1;p++)
		convergence_file<< points[p][i] << " ";
	convergence_file<<std::endl;
        }
      }
    }//numConvergence end
    //************************
 }//main end
