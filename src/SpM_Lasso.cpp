#include"SpM_Lasso.hpp"
#include <Eigen/Cholesky>
#include<iostream>
double Lasso::find_lambda_max() const{
/*lambda_max = norm( A'*b, 'inf' );
 */
  return ATb_.lpNorm<Eigen::Infinity>();
}
void Lasso::factor(){
/*
 *     [m, n] = size(A);
    if ( m >= n )    % if skinny
       L = chol( A'*A + rho*speye(n), 'lower' );
    else            % if fat
       L = chol( speye(m) + 1/rho*(A*A'), 'lower' );
    end

    % force matlab to recognize the upper / lower triangular structure
    L = sparse(L);
    U = sparse(L');
 *
 */
  if(m_>=n_){
    llt_.compute(A_.transpose()*A_+rho_*Eigen::MatrixXd::Identity(n_,n_)); // compute the Cholesky decomposition of A, skinny case
  }else{
    llt_.compute(Eigen::MatrixXd::Identity(m_,m_)+1./rho_*(A_*A_.transpose())); // compute the Cholesky decomposition of A, fat case
  }
}
void Lasso::shrinkage(const Eigen::VectorXd &x, double kappa, Eigen::VectorXd &z){
  //z = max( 0, x - kappa ) - max( 0, -x - kappa );
  for(int i=0;i<z.size();++i){
    z[i]=std::max(0., x[i]-kappa)-std::max(0., -x[i]-kappa);
  }
}
void Lasso::lasso_solve(double lambda, double alpha, Eigen::VectorXd &x){
  std::cout<<"LASSO solver for lambda: "<<lambda<<" alpha: "<<alpha<<std::endl;
  printf("%3s\t%10s\t%10s\t%10s\t%10s\t%10s\n", "iter", "r norm", "eps pri", "s norm", "eps dual", "objective");

  Eigen::VectorXd q=Eigen::VectorXd::Zero(n_);
  z_=Eigen::VectorXd::Zero(n_);
  u_=Eigen::VectorXd::Zero(n_);
  x=Eigen::VectorXd::Zero(n_);

  for(int k=0;k<lasso_max_iter_;++k){
    /*    % x-update
    q = Atb + rho*(z - u);    % temporary value
    if( m >= n )    % if skinny
       x = U \ (L \ q);
    else            % if fat
       x = q/rho - (A'*(U \ ( L \ (A*q) )))/rho^2;
    end
     */

    q=ATb_+rho_*(z_-u_);
    if(m_>=n_){ //skinny case
      x=llt_.solve(q);
    }
    else{ //fat case
      x=q/rho_- (A_.transpose()*(llt_.solve(A_*q) ))/(rho_*rho_);
    }

    /*z-update with relaxation
    zold = z;
    x_hat = alpha*x + (1 - alpha)*zold;
    z = shrinkage(x_hat + u, lambda/rho);*/
    Eigen::VectorXd zold=z_;
    Eigen::VectorXd x_hat=alpha*x+(1-alpha)*zold;

    shrinkage(x_hat+u_, lambda/rho_, z_);

    /*    % u-update
    u = u + (x_hat - z);
     */
    u_+=x_hat-z_;

    //..and we're done with the step. now compute objectives and norm
    double obj=objective(lambda, x);
    double norm_r=r_norm(x);
    double norm_s=z_norm(zold);
    double eps_pri=std::sqrt(n_)*lasso_abs_tol_ + lasso_rel_tol_*std::max(x.norm(), (-z_).norm()); //minus z seems silly but is in the orig code
    double eps_dual=std::sqrt(n_)*lasso_abs_tol_ + lasso_rel_tol_*((rho_*u_).norm());

    printf("%3d\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.2f\n", k, norm_r, eps_pri,norm_s, eps_dual, obj);

    if (norm_r < eps_pri && norm_s < eps_dual)
         break;
  }
}
double Lasso::objective(double lambda, Eigen::VectorXd &x) const{
  return 0.5*((A_*x - b_).squaredNorm()) + lambda*(z_.lpNorm<1>());
}
double Lasso::r_norm(const Eigen::VectorXd &x) const{
 return (x - z_).norm();
}
double Lasso::z_norm(const Eigen::VectorXd &zold) const{
  return (-rho_*(z_ - zold)).norm();
}
double Lasso::matrix_norm(const Eigen::MatrixXd &M){
  //there has to be a better way of computing a norm than doing an SVD.
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(M, Eigen::ComputeThinU | Eigen::ComputeThinV);
  return svd.singularValues()[0];
}
