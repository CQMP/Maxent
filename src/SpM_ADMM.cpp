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
  Eigen::VectorXd rhs=1./lambda_*S_.cwiseProduct(yprime_)+muprime_*(zprime_-uprime_)+mu_*Vt_*(z_-u_);
  std::cout<<"mu Vt (z-u) norm: "<<(mu_*Vt_*(z_-u_)).norm()<<" cwise prod: "<<(x1_update_denominator_.cwiseProduct(mu_*Vt_*(z_-u_))).norm()<<std::endl;
  std::cout<<"x1 denom: "<<x1_update_denominator_<<std::endl;
  xi1=x1_update_denominator_.cwiseProduct(rhs);
}
void ADMM::solve_for_xi2(Eigen::VectorXd &xi2) const{
  Eigen::VectorXd rhs=Vt_*Eigen::VectorXd::Ones(nw_);
  xi2=x1_update_denominator_.cwiseProduct(rhs);
}
void ADMM::update_xprime(Eigen::VectorXd &xprime)const {
   xprime=xi1_+nu_*xi2_;
}
void ADMM::update_zprime(Eigen::VectorXd &zprime)const{
  //std::cout<<"st with alpha: "<<1./muprime_<<" sum: "<<xprime_+uprime_<<std::endl;
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
  positive_projection(Vt_.transpose()*xprime_+u_, z);
}
void ADMM::positive_projection(const Eigen::VectorXd &v, Eigen::VectorXd &z) const{
  for(int i=0;i<v.size();++i){
    z[i]=std::max(0., v[i]);
  }
}
void ADMM::update_u(Eigen::VectorXd &u) const{
  u+=Vt_.transpose()*xprime_-z_;
}
double ADMM::update_nu() const{
  return (1.-(Vt_.transpose()*xi1_).sum())/(Vt_.transpose()*xi2_).sum();
}

void ADMM::iterate(){
  solve_for_xi1(xi1_);
  solve_for_xi2(xi2_);
  nu_=update_nu();
  //std::cout<<"old xprime: "<<xprime_<<std::endl;
  update_xprime(xprime_);
  //std::cout<<" new xprime: "<<xprime_<<std::endl;
  //std::cout<<" old zprime: "<<zprime_<<std::endl;
  update_zprime(zprime_);
  //std::cout<<" new zprime: "<<zprime_<<std::endl;
  //std::cout<<" old uprime: "<<uprime_<<std::endl;
  update_uprime(uprime_);
  //std::cout<<" new uprime: "<<uprime_<<std::endl;
  update_z(z_);
  update_u(u_);
  std::cout<<" xprime: "<<xprime_.norm()<<" zprime: "<<zprime_.norm()<<" uprime: "<<uprime_.norm()<<" z: "<<z_.norm()<<" u: "<<u_.norm()<<std::endl;
}

void ADMM::print_info(std::ostream &os) const{
  //os<<Vt_.transpose()*xprime_<<std::endl;
  os<<"constraint violations: z'-x': "
      <<constraint_violation_zprime_xprime()
      <<" z-Vx': "
      <<constraint_violation_z_Vxprime()
      <<" norm: "
      <<constraint_violation_norm()
      <<" positivity: "
      <<constraint_violation_positivity()
      <<" objective: "
      <<objective_functional();

}
double ADMM::constraint_violation_zprime_xprime() const{
  return (zprime_-xprime_).norm();
}
double ADMM::constraint_violation_z_Vxprime() const{
  return (z_-Vt_.transpose()*xprime_).norm();
}
double ADMM::constraint_violation_norm() const{
  return std::abs((Vt_.transpose()*xprime_).sum()-1);
}
double ADMM::constraint_violation_positivity() const{
  Eigen::VectorXd x=Vt_.transpose()*xprime_;
  for(int i=0;i<x.size();++i) if(x[i]>0) x[i]=0;
  return x.norm();
}
double ADMM::objective_functional() const{
  return 0.5*(yprime_-S_.cwiseProduct(xprime_)).squaredNorm()+lambda_*xprime_.lpNorm<1>();
}
