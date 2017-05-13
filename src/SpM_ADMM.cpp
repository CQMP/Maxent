#include"SpM_ADMM.hpp"
#include<Eigen/Cholesky>
#include<iostream>


void ADMM::soft_threshold(double alpha, const Eigen::VectorXd &x, Eigen::VectorXd &z) const{
  for(int i=0;i<ns_;++i){
    if(x[i]>=alpha)
      z[i]=x[i]-alpha;
    else if(std::abs(x[i])<alpha)
      z[i]=0;
    else
      z[i]=x[i]+alpha;
  }
}
void ADMM::positive_projection(const Eigen::VectorXd &v, Eigen::VectorXd &z) const{
  for(int i=0;i<v.size();++i){
    z[i]=std::max(0., v[i]);
    //std::cout<<i<<" "<<v[i]<<" "<<z[i]<<std::endl;
  }
  //exit(0);
}

void ADMM::compute_A(){
  for(int i=0;i<ns_;++i){
    A_(i,i)=1.;//this is the 'x' part of the constraint x'-z'=0. All other elements are zero
  }
  for(int i=0;i<ns_;++i){
    //this is the integral over the spectral function that 'c' will enforce to be 1.
    for(int k=0;k<nw_;++k)
      A_(ns_,i)+=Vt_(i,k);
  }
}
void ADMM::compute_B(){
  for(int i=0;i<ns_;++i){
    B_(i,i)=-1.; //this is the '-' of the constraint x'-z'=0. All other elements are zero
  }
}
void ADMM::compute_c(){
  c_(ns_)=1.; //this is the norm. All other elements are zero
}
void ADMM::decompose_P(){
  Eigen::MatrixXd Sprod=S_.cwiseProduct(S_).asDiagonal();
  P_.compute(
      Sprod
      +rho_*(A_.transpose()*A_)
      );
}
void ADMM::compute_q(){
  q_=-yprime_.cwiseProduct(S_);
      +rho_*(A_.transpose()*(B_*zprime_-c_+uprime_));
  if(enforce_positivity_){
    q_-=rhoprime_*(wprime_-sprime_);
  }
}

void ADMM::update_xprime() {
  compute_q();
  xprime_=P_.solve(q_);
}
void ADMM::update_zprime(){
  Eigen::VectorXd w=A_*xprime_+uprime_-c_;
  Eigen::VectorXd v=-B_.transpose()*w;
  soft_threshold(lambda_/rho_, v, zprime_);
}
void ADMM::update_uprime(){
  uprime_+=A_*xprime_+B_*zprime_-c_;
}

void ADMM::update_sprime(){
  //tmp variable, s'=Vt*s
  static Eigen::VectorXd s(nw_);
  positive_projection(Vt_.transpose()*(xprime_+wprime_), s);
  sprime_=Vt_*s;
}
void ADMM::update_wprime(){
  wprime_+=xprime_-sprime_;
}

void ADMM::iterate(){
  update_xprime();
  update_zprime();
  update_uprime();
  if(enforce_positivity_){
    update_sprime();
    update_wprime();
  }
}

void ADMM::print_info(std::ostream &os) const{
  os<<"constraint violations: z'-x': "
      <<constraint_violation_zprime_xprime()
      <<" s'-x: "
      <<constraint_violation_sprime_xprime()<<std::endl;
  os<<" norm: "
      <<constraint_violation_norm()
      <<" positivity: "
      <<constraint_violation_positivity()<<std::endl;
  os<<" chi2 term: "
      <<chisquare_term()
      <<" l1 of x': "
      <<l1_of_xprime()<<std::endl;
  os<<" objective: "
      <<global_objective_functional()<<std::endl;
  std::cout<<"moving on."<<std::endl;
}
double ADMM::constraint_violation_zprime_xprime() const{
  return (zprime_-xprime_).norm();
}
double ADMM::constraint_violation_sprime_xprime() const{
  return (sprime_-xprime_).norm();
}
double ADMM::constraint_violation_norm() const{
  return std::abs((Vt_.transpose()*xprime_).sum()-1);
}
double ADMM::constraint_violation_positivity() const{
  Eigen::VectorXd x=spectral_function().cwiseProduct(domega_);
  for(int i=0;i<x.size();++i) if(x[i]>0) x[i]=0;
  return x.norm();
}
double ADMM::global_objective_functional() const{
  return chisquare_term()+lambda_*l1_of_xprime();
}

double ADMM::l1_of_xprime() const{
  return xprime_.lpNorm<1>();
}
double ADMM::chisquare_term() const{
  return 0.5*(yprime_-S_.cwiseProduct(xprime_)).squaredNorm();
}

Eigen::VectorXd ADMM::spectral_function() const{
  return (Vt_.transpose()*xprime_).cwiseQuotient(domega_);
}
void ADMM::check_regular(const Eigen::VectorXd &v, const std::string &msg) const{
  for(int i=0;i<v.size();++i) if(std::isnan(v[i])||std::isinf(v[i])) throw std::runtime_error(msg+"encountered nan or inf.");
}
