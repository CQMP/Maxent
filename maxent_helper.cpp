/*
 * Copyright (C) 1998-2016 ALPS Collaboration. See COPYRIGHT.TXT
 * All rights reserved. Use is subject to license terms. See LICENSE.TXT
 * For use in publications, see ACKNOWLEDGE.TXT
 */

#include "maxent.hpp"
#include <alps/config.hpp> // needed to set up correct bindings
#include <boost/math/special_functions/fpclassify.hpp> //needed for boost::math::isnan
#include <Eigen/Eigenvalues>
#include <Eigen/Cholesky>
#include "maxent_backcont.hpp"
//NOTE: size1= rows; size2=columns

MaxEntHelper::MaxEntHelper(alps::params& p) :
MaxEntParameters(p) , def_(nfreq())
{
    for (int i=0; i<nfreq(); ++i)
        def_[i] = MaxEntParameters::Default().D(omega_coord(i)) * delta_omega(i);
    //normalizing the default model
    //def_ /= sum(def_);
    def_ /= def_.sum();
    checkDefaultModel(def_);
}

/// check that the default model is non-zero
/// this is needed for transform_into_singular_space
/// to work safely
void MaxEntHelper::checkDefaultModel(const vector_type &D) const{
    for(int i=0;i<D.size();i++){
        double Di=D(i);
        if(Di==0 || boost::math::isnan(Di))
          throw std::logic_error("dude, your D is zero");
    }
}

//the opposite of 'transform_into_real_space'; takes a vector 'A' and makes a vector 'u' out of it:
// u = V^T* log(A/Default)
vector_type MaxEntHelper::transform_into_singular_space(vector_type A) const
{
  double D;
  for (unsigned int i=0; i<A.size(); ++i) {
    D=Default(i);
    A[i] /= D;
    A[i] = A[i]==0. ? 0. : log(A[i]);
  }
  return maxent_prec_prod(Vt(), A);
}

//returns exp(V^T*u)*Default(i). This quantity is then usually called 'A'
vector_type MaxEntHelper::transform_into_real_space(vector_type u) const
{
  u = maxent_prec_prod(Vt().transpose(), u);
  for (unsigned int i=0; i<u.size(); ++i) {
    u[i] = exp(u[i]);
    u[i] *= Default(i);
  }
  return u; 
}


vector_type MaxEntHelper::get_spectrum(const vector_type& u) const
{
  vector_type A = transform_into_real_space(u);
  for (unsigned int i=0; i<A.size(); ++i) 
    A[i] /= delta_omega(i);
  return A;
}


//'left side' is defined as Sigma*(V^T*RealSpace(u)*V)*Sigma
//see Bryan's paper near Eq. 11 
matrix_type MaxEntHelper::left_side(const vector_type& u) const
{
  vector_type A = transform_into_real_space(u);
  matrix_type M = Vt().transpose();
  for (unsigned int i=0; i<M.rows(); ++i) 
    for (unsigned int j=0; j<M.cols(); ++j) 
      M(i,j) *= A[i];
  M = maxent_prec_prod(Vt(), M);
  M = maxent_prec_prod(Sigma() ,M);
  M = maxent_prec_prod(Sigma(), M);
  M *= 2./ndat();
  return M;
}


//this function computes
//Sigma*U^T*(K*RealSpace(u)-y)
//up to a factor of 2./ndat(). Compare this to Eq. D.12 in Sebastian's thesis
vector_type MaxEntHelper::right_side(const vector_type& u) const
{
  vector_type b = 2./ndat()*(maxent_prec_prod(K(), transform_into_real_space(u)) - y());
  b = maxent_prec_prod(U().transpose(), b);
  b = maxent_prec_prod(Sigma(), b);
  return b;
}

//this function constructs delta \dot (V^T*RealSpace(u)*V)
double MaxEntHelper::step_length(const vector_type& delta, const vector_type& u) const 
{
  vector_type A = transform_into_real_space(u);
  matrix_type L = Vt().transpose();
  for (unsigned int i=0; i<L.rows(); ++i) 
    for (unsigned int j=0; j<L.cols(); ++j) 
      L(i,j) *= A[i];
  L = maxent_prec_prod(Vt(), L);
  return delta.dot(maxent_prec_prod(L, delta));
}

//Bryan's paper section 2.3 (or after eq 22)
double MaxEntHelper::convergence(const vector_type& u, const double alpha) const 
{
  //using namespace boost::numeric::ublas;
  vector_type A = transform_into_real_space(u);
  matrix_type L = Vt().transpose();
  for (unsigned int i=0; i<L.rows(); ++i) 
    for (unsigned int j=0; j<L.cols(); ++j) 
      L(i,j) *= A[i];
  L = maxent_prec_prod(Vt(), L);
  vector_type alpha_dSdu = -alpha*maxent_prec_prod(L, u);
  vector_type dLdu = maxent_prec_prod(L, right_side(u));
  vector_type diff = alpha_dSdu - dLdu;
  double denom = alpha_dSdu.norm() + dLdu.norm();
  denom = denom*denom;
  return 2*diff.dot(diff)/denom;
}

double MaxEntHelper::log_prob(const vector_type& u, const double alpha) const
{
  matrix_type L = maxent_prec_prod_trans(K(), K());
  const vector_type A = transform_into_real_space(u);
  for (unsigned int i=0; i<L.rows(); ++i)
    for (unsigned int j=0; j<L.cols(); ++j)
      L(i,j) *= sqrt(A[i])*sqrt(A[j]);
  for (unsigned int i=0; i<L.rows(); ++i)
    L(i,i) += alpha;
  //boost::numeric::bindings::lapack::potrf(boost::numeric::bindings::lower(L));
  
  //LAPACK does potrf in place, while Eigen does not
  Eigen::LLT<matrix_type,Eigen::Lower> lltofL(L);
  L=lltofL.matrixL();

  double log_det = 0.;
  for (unsigned int i=0; i<L.rows(); ++i) 
    log_det  += log(L(i,i)*L(i,i));
  return 0.5*( (nfreq())*log(alpha) - log_det ) - Q(u, alpha);
}

double MaxEntHelper::chi_scale_factor(vector_type A, const double chi_sq, const double alpha) const
{
  for (unsigned int i=0; i<A.size(); ++i) 
    A[i] *= delta_omega(i);
  using namespace boost::numeric;
  matrix_type L = maxent_prec_prod_trans(K(), K());
  for (unsigned int i=0; i<L.rows(); ++i)
    for (unsigned int j=0; j<L.cols(); ++j)
      L(i,j) *= sqrt(A[i])*sqrt(A[j]);
  vector_type lambda(L.rows());
  //bindings::lapack::syev('N', bindings::upper(L) , lambda, bindings::lapack::optimal_workspace());
  Eigen::SelfAdjointEigenSolver<matrix_type> es(L);
  lambda = es.eigenvalues();
  double Ng = 0.;
  for (unsigned int i=0; i<lambda.size(); ++i) {
    if (lambda[i]>=0) 
      Ng += lambda[i]/(lambda[i]+alpha);
  }
  std::cerr << "Ng: " << Ng << std::endl;
  std::cerr << "chi2 max: " << chi_sq << std::endl;
  return sqrt(chi_sq/(ndat()-Ng));
}

//This function computes chi^2 as in equation D.6 in Sebastian's thesis
double MaxEntHelper::chi2(const vector_type& A) const 
{
  vector_type del_G = maxent_prec_prod(K(), A) - y();
  
  /*std::cout<<"in computation of chi2:"<<std::endl;
   for(int i=0;i<y().size();++i){
   std::cout<<i<<" "<<maxent_prec_prod(K(), A)[i]<<" "<<y()[i]<<std::endl;
   }*/
  
  double c = 0;
  for (unsigned int i=0; i<del_G.size(); ++i) 
    c += del_G[i]*del_G[i];
  return c;
}

void MaxEntHelper::print_chi2(const vector_type& A, std::ostream &os) const 
{
  vector_type backcont=maxent_prec_prod(K(), A);
  vector_type defaultm=maxent_prec_prod(K(), Default());
  os<<"#first column: index (Matsubara frequency). second column: fitted function. third: input data. fourth: default model."<<std::endl;
  for(int i=0;i<y().size();++i){
    os<<i<<" "<<backcont[i]<<" "<<y()[i]<<" "<<defaultm[i]<<std::endl;
  }
  os<<std::endl;
}

//This function computes the entropy as in Eq. D.7 in Sebastian's thesis
//Really eq 3.17 from Jarrell and Gubernatis
double MaxEntHelper::entropy(const vector_type& A) const 
{
  double S = 0;
  for (unsigned int i=0; i<A.size(); ++i) {
    double lg = A[i]==0. ? 0. : log(A[i]/Default(i));
    S += A[i] - Default(i) - A[i]*lg;
  }
  return S;
}

// Generic principal-value integration. Assumes that spectral function does not vary too
// strongly within [x_{i+1},x_i]
vector_type MaxEntHelper::PrincipalValue(const vector_type &w,const vector_type &a) const
{
    int N = w.size();
    vector_type r(N);
    for (int i=2;i<N-2;i++) {
      double scr = -a[i]*std::log(std::abs((w[i+1]-w[i])/(w[i]-w[i-1])));
      for (int j=0;j<i-1;j++) scr -= 0.5*(a[j]+a[j+1])*std::log(std::abs((w[i]-w[j+1])/(w[i]-w[j])));
      for (int j=i+1;j<N-2;j++) scr -= 0.5*(a[j]+a[j+1])*std::log(std::abs((w[j+1]-w[i])/(w[j]-w[i])));
      r[i] = scr;
    }
    r[0] = w[2]*r[2]/w[0];
    r[1] = w[2]*r[2]/w[1];
    r[N-2] = w[N-3]*r[N-3]/w[N-2];
    r[N-1] = w[N-3]*r[N-3]/w[N-1];
    return r;
}

void MaxEntHelper::backcontinue(ofstream_ &os, const vector_type &A) const
{
    const MaxEntParameters *pp = this;
    Backcont bc(pp);
    vector_type G = bc.backcontinue(A);
    kernel_type k_type = pp->getKernelType();
    double beta = 1/(pp->T());
    for(int n=0; n<G.size();n++)
      os << pp->inputGrid(n) << " " << G(n) << std::endl;
}
