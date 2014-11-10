#include<iostream>
#include<vector>

extern "C" void dqawc_(double (*f)(double *),double &lower_limit ,double &upper_limit,double &pole,double &abstol,double &reltol, double &result, double &err,int &neval,int &ier,
     int &limit,int &lenw,int &last,int *iwork, double *work);
double f(double*x){ return 1.; }
int main(int argc, char **argv){
  //Integral and domain characterization
  double lower_limit=-1.;
  double upper_limit=2.;
  double pole=0.;
  //error estimate
  double abstol=1.e-12;
  double reltol=1.e-12;

  //return values
  double result;
  double err;
  int neval;
  int ier;

  //dimensioning parameters
  int limit=500;
  int lenw=limit*8;
  int last; //output: number of subdivisions
  std::vector<int> iwork(limit);
  std::vector<double> work(lenw);
/*
c***purpose  the routine calculates an approximation result to a
c            cauchy principal value i = integral of f*w over (a,b)
c            (w(x) = 1/((x-c), c.ne.a, c.ne.b), hopefully satisfying
c            following claim for accuracy
c            abs(i-result).le.max(epsabe,epsrel*abs(i)).
*/
  dqawc_(f,lower_limit, upper_limit,pole,abstol,reltol,result,err,neval,ier,limit,lenw,last,&(iwork[0]),&(work[0]));

  //analysis of output values:
  std::cout<<"result: "<<result<<std::endl;
  std::cout<<"error: "<<err<<std::endl;
  std::cout<<"neval: "<<neval<<std::endl;
  std::cout<<"ier: "<<ier<<std::endl;
}
