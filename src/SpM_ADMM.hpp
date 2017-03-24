#include"maxent_matrix_def.hpp"

///Class for performing  ADMM fit as described in paper, Eq. S5

class ADMM{
public:
  ADMM(const Eigen::MatrixXd &Vt, const Eigen::VectorXd &yprime, const Eigen::VectorXd &domega, const Eigen::VectorXd &S, double muprime, double mu, double lambda):
  S_(S)
, Vt_(Vt)
, yprime_(yprime)
, domega_(domega)
, ns_(S.size())
, nw_(Vt.cols())
, muprime_(muprime)
, mu_(mu)
, lambda_(lambda)
, nu_(0)
, xprime_(Eigen::VectorXd::Zero(ns_))
, zprime_(Eigen::VectorXd::Zero(ns_))
, uprime_(Eigen::VectorXd::Zero(ns_))
, xi1_(Eigen::VectorXd::Zero(ns_))
, xi2_(Eigen::VectorXd::Zero(ns_))
, z_(Eigen::VectorXd::Zero(nw_))
, u_(Eigen::VectorXd::Zero(nw_))
, x1_update_denominator_(Eigen::VectorXd::Zero(ns_))
{
    if(Vt_.rows()!=ns_) throw std::runtime_error("Vt should have ns rows");
    if(yprime_.size()!=ns_) throw std::runtime_error("yprime_ should have size ns");
    compute_denominator();
  }
  ///solve for the first part of Eq. S5a
  void solve_for_xi1(Eigen::VectorXd &xi1) const;
  ///solve for the second part of Eq. S5a
  void solve_for_xi2(Eigen::VectorXd &xi2) const;
  ///update xprime as in Eq. S5a
  void update_xprime(Eigen::VectorXd &xprime) const;
  ///update zprime as in Eq. S5b
  void update_zprime(Eigen::VectorXd &zprime) const;
  ///update uprime as in Eq. S5c
  void update_uprime(Eigen::VectorXd &uprime) const;
  ///update z as in Eq. S5d
  void update_z(Eigen::VectorXd &z) const;
  ///update u as in Eq. S5e
  void update_u(Eigen::VectorXd &u) const;

  ///run one ADMM iteration
  void iterate();
  ///statistics output
  void print_info(std::ostream &os) const;
  ///at convergence, xprime=zprime. This gives the norm of the difference.
  double constraint_violation_zprime_xprime() const;
  ///at convergence, z=V xprime. This gives the norm of the difference
  double constraint_violation_z_Vxprime() const;
  ///the resulting vector x should be normalized. this is the deviation from the norm.
  double constraint_violation_norm() const;
  ///this checks the positivity of the solution
  double constraint_violation_positivity() const;
  ///this returns the current objective (Eq. S1)
  double objective_functional() const;

  ///this transforms the continued spectral function back from singular to real space
  Eigen::VectorXd spectral_function() const;
private:
  ///the first line of Eq. S5a
  void compute_denominator();
  ///soft threshold function S7
  void soft_threshold(double alpha, const Eigen::VectorXd &x, Eigen::VectorXd &z) const;
  ///take a vector v, copy to z, set all negative values to zero.
  void positive_projection(const Eigen::VectorXd &v, Eigen::VectorXd &z) const;
  ///compute the value of nu, Eq. S6
  double update_nu() const;
  ///vector of singular values
  const Eigen::VectorXd S_;
  ///transformation matrix Vt
  const Eigen::MatrixXd Vt_;
  ///vector of data values
  const Eigen::VectorXd yprime_;
  ///vector of frequency grid discretizations
  const Eigen::VectorXd domega_;

  ///number of singular values
  const int ns_;
  ///number of frequency values
  const int nw_;


  ///penalty parameter for z'-u'
  const double muprime_;
  ///penalty parameter for Vt(z-u)
  const double mu_;
  ///L1 norm parameter
  double lambda_;

  ///temporary doubles
  double nu_;

  ///temporary vectors
  Eigen::VectorXd xprime_;
  Eigen::VectorXd zprime_;
  Eigen::VectorXd uprime_;
  Eigen::VectorXd xi1_;
  Eigen::VectorXd xi2_;
  Eigen::VectorXd z_;
  Eigen::VectorXd u_;

  ///denominator vector of xprime update
  Eigen::VectorXd x1_update_denominator_;

};
