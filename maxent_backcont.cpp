#include "maxent_backcont.hpp"

Backcont::Backcont(const MaxEntParameters *param_in) :
 param(param_in),k_type(param->getKernelType()){
} 

vector_type Backcont::backcontinue(const vector_type &A){
  const int ndat = param->ndat();
  const int nfreq = param->nfreq();

  G = vector_type(ndat);
  // G(X) =\int K(X,\omega)A(\omega)d\omega
  for(int n=0;n<ndat;n++){
    //recall that in MaxEntParameters we scaled the kernel by the 1/error
    double err = param->sigma(n);
    G(n)=0;
    for(int i=1;i<nfreq-1;i++){
      //G(n) += -1./M_PI*(param->K(n,i))*A(i)*(param->delta_omega(i)); 
      double delta = (param->omega_coord(i+1)) - (param->omega_coord(i-1));
      G(n) += err*(param->K(n,i))*A(i)*delta/2; 
    }

    G(n) += err*(param->K(n,0))*A(0)*(param->omega_coord(1)-param->omega_coord(0))/2;
    G(n) += err*(param->K(n,nfreq-1))*A(nfreq-1)*(param->omega_coord(nfreq-1)-param->omega_coord(nfreq-2))/2;
  } 

  return G;
}
