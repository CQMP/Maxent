#include"maxent_matrix_def.hpp"
#include<iostream>
///Class for performing  ADMM fit as described in paper, Eq. S5

class ADMM{
public:
  ADMM(const Eigen::MatrixXd &Vt,const Eigen::VectorXd &yprime, const Eigen::VectorXd &domega, const Eigen::VectorXd &S, double rho, double rhoprime, double lambda):
  S_(S)
, Vt_(Vt)
, yprime_(yprime)
, domega_(domega)
, ns_(S.size())
, nw_(Vt.cols())
, rho_(rho)
, rhoprime_(rhoprime)
, lambda_(lambda)
//, nu_(0)
, xprime_(Eigen::VectorXd::Zero(ns_))
, zprime_(Eigen::VectorXd::Zero(ns_))
, uprime_(Eigen::VectorXd::Zero(ns_+1))
, sprime_(Eigen::VectorXd::Zero(ns_))
, wprime_(Eigen::VectorXd::Zero(ns_))
, A_(Eigen::MatrixXd::Zero(ns_+1, ns_))
, B_(Eigen::MatrixXd::Zero(ns_+1, ns_))
, c_(Eigen::VectorXd::Zero(ns_+1))
, q_(Eigen::VectorXd::Zero(ns_))
, enforce_positivity_(true)
{
    if(Vt_.rows()!=ns_) throw std::runtime_error("Vt should have ns rows");
    if(Vt_.cols()!=nw_) throw std::runtime_error("Vt should have nw rows");
    if(domega_.size()!=nw_) throw std::runtime_error("domega should have nw rows");
    if(yprime_.size()!=ns_) throw std::runtime_error("yprime_ should have size ns");
    compute_A();
    compute_B();
    compute_c();
    decompose_P();
  }
  ///update xprime
  void update_xprime();
  ///update zprime
  void update_zprime();
  ///update uprime
  void update_uprime();
  ///update sprime
  void update_sprime();
  ///update wprime
  void update_wprime();


  ///run one ADMM iteration
  void iterate();
  ///statistics output
  void print_info(std::ostream &os) const;
  ///at convergence, x'=z'. This gives the norm of the difference.
  double constraint_violation_zprime_xprime() const;
  ///at convergence, s'=z'. This gives the norm of the difference
  double constraint_violation_sprime_xprime() const;
  ///the resulting vector x should be normalized. this is the deviation from the norm.
  double constraint_violation_norm() const;
  ///this checks the positivity of the solution
  double constraint_violation_positivity() const;
  ///this returns the current objective (Eq. S1)
  double global_objective_functional() const;
  ///this returns the chi2 part of the objective
  double chisquare_term() const;
  ///this returns the chi2 part of the objective
  double l1_of_xprime() const;

  ///this transforms the continued spectral function back from singular to real space
  Eigen::VectorXd spectral_function() const;

  ///getter functions, mostly for testing
  const Eigen::VectorXd &zprime() const{return zprime_;}
  const Eigen::VectorXd &xprime() const{return xprime_;}

  void externally_set_xprime(const Eigen::VectorXd &xprime){xprime_=xprime;}
  void externally_set_zprime(const Eigen::VectorXd &zprime){zprime_=zprime;}

private:
  ///the constraint matrix Ax'+By'=c
  void compute_A();
  ///the constraint matrix Ax'+By'=c
  void compute_B();
  ///the constraint matrix Ax'+By'=c
  void compute_c();
  ///the q side of the matrix equation minimization x^T P x+q^T x+c=0 for the x' minimization
  void compute_q();
  ///the P side of the matrix equation minimization x^T P x+q^T x+c=0 for the x' minimization
  void decompose_P();

  ///soft threshold function S7
  void soft_threshold(double alpha, const Eigen::VectorXd &x, Eigen::VectorXd &z) const;
  ///take a vector v, copy to z, set all negative values to zero.
  void positive_projection(const Eigen::VectorXd &v, Eigen::VectorXd &sprime) const;

  ///helper debug function
  void check_regular(const Eigen::VectorXd &v, const std::string &msg) const;
  ///vector of singular values
  const Eigen::VectorXd S_;
  ///transformation matrix Vt, the transpose of V. 'real space' on the right side, singular space on the left side, x'=V^T\rho.
  Eigen::MatrixXd Vt_;
  ///vector of data values
  const Eigen::VectorXd yprime_;
  ///vector of frequency grid discretizations
  const Eigen::VectorXd domega_;

  ///number of singular values
  const int ns_;
  ///number of frequency values
  const int nw_;


  ///penalty parameter for the L1 norm
  const double rho_;
  ///penalty parameter for the positivity
  const double rhoprime_;
  ///L1 norm parameter
  double lambda_;

  ///temporary doubles
  //double nu_;

  ///singular space vector for x
  Eigen::VectorXd xprime_;
  ///dual and running residual for L1 optimization
  Eigen::VectorXd zprime_;
  Eigen::VectorXd uprime_;
  ///dual and running residual for positivity optimization
  Eigen::VectorXd sprime_;
  Eigen::VectorXd wprime_;

  ///Constraint matrix of xprime update, Ax'+Bz'=c
  Eigen::MatrixXd A_;
  ///Constraint matrix of xprime update, Ax'+Bz'=c
  Eigen::MatrixXd B_;
  ///Constraint vector of xprime update, Ax'+Bz'=c
  Eigen::VectorXd c_;
  ///decomposed matrix P for the inverse
  Eigen::LDLT<Eigen::MatrixXd > P_;
  ///righthand side q for the inverse, for optimizing x^TPx+qx^t+c=0
  Eigen::VectorXd q_;

  const bool enforce_positivity_;

};
