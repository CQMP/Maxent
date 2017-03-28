#include"SpM_ADMM.hpp"
#include<Eigen/Cholesky>
#include<iostream>

void ADMM::compute_denominator(){
  x1_update_denominator_.resize(ns_);
  for(int i=0;i<ns_;++i){
    x1_update_denominator_[i]=1./(S_[i]*S_[i]/lambda_+(muprime_+mu_));
  }
}

void ADMM::solve_for_xi1(Eigen::VectorXd &xi1) const{
  Eigen::VectorXd rhs=1./lambda_*S_.cwiseProduct(yprime_)
      +muprime_*(zprime_-uprime_)+
      mu_*Vt_scaled_*(z_-u_);
  xi1=x1_update_denominator_.cwiseProduct(rhs);
}
void ADMM::solve_for_xi2(Eigen::VectorXd &xi2) const{
  Eigen::VectorXd rhs=Vt_scaled_*Eigen::VectorXd::Ones(nw_);
  xi2=x1_update_denominator_.cwiseProduct(rhs);
}
void ADMM::update_xprime(Eigen::VectorXd &xprime)const {
   xprime=xi1_+nu_*xi2_;
}
void ADMM::update_zprime(Eigen::VectorXd &zprime)const{
  soft_threshold(1./muprime_, xprime_+uprime_, zprime);
}
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
void ADMM::update_uprime(Eigen::VectorXd &uprime) const{
  uprime+=xprime_-zprime_;
}
void ADMM::update_z(Eigen::VectorXd &z) const{
  positive_projection(spectral_function()+u_, z);
}
void ADMM::positive_projection(const Eigen::VectorXd &v, Eigen::VectorXd &z) const{
  for(int i=0;i<v.size();++i){
    z[i]=std::max(0., v[i]);
  }
}
void ADMM::update_u(Eigen::VectorXd &u) const{
  u+=spectral_function()-z_;
}
double ADMM::update_nu() const{
  return (1.-(V_scaled_*xi1_).sum())/(V_scaled_*xi2_).sum();
}

void ADMM::iterate(){
  solve_for_xi1(xi1_);
  solve_for_xi2(xi2_);
  nu_=update_nu();
  update_xprime(xprime_);
  update_zprime(zprime_);
  update_uprime(uprime_);
  update_z(z_);
  update_u(u_);
}

void ADMM::print_info(std::ostream &os) const{
  os<<"constraint violations: z'-x': "
      <<constraint_violation_zprime_xprime()
      <<" z-Vx': "
      <<constraint_violation_z_Vxprime()
      <<" norm: "
      <<constraint_violation_norm()
      <<" positivity: "
      <<constraint_violation_positivity()
      <<" chi2 term: "
      <<chisquare_term()
      <<" l1 of x': "
      <<l1_of_xprime()
      <<" objective: "
      <<objective_functional();
}
double ADMM::constraint_violation_zprime_xprime() const{
  return (zprime_-xprime_).norm();
}
double ADMM::constraint_violation_z_Vxprime() const{
  return (z_-spectral_function()).norm();
}
double ADMM::constraint_violation_norm() const{
  return std::abs(spectral_function().sum()-1);
}
double ADMM::constraint_violation_positivity() const{
  Eigen::VectorXd x=spectral_function();
  for(int i=0;i<x.size();++i) if(x[i]>0) x[i]=0;
  return x.norm();
}
double ADMM::objective_functional() const{
  return chisquare_term()+lambda_*l1_of_xprime();
}
double ADMM::l1_of_xprime() const{
  return xprime_.lpNorm<1>();
}
double ADMM::chisquare_term() const{
  return 0.5*(yprime_-S_.cwiseProduct(xprime_)).squaredNorm();
}

Eigen::VectorXd ADMM::spectral_function() const{
  return (V_scaled_*xprime_);//.cwiseQuotient(domega_);
}
void ADMM::check_regular(const Eigen::VectorXd &v, const std::string &msg) const{
  for(int i=0;i<v.size();++i) if(std::isnan(v[i])||std::isinf(v[i])) throw std::runtime_error(msg+"encountered nan or inf.");
}
