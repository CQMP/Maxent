#include"maxent_matrix_def.hpp"

///Class for performing LASSO fitting using ADMM
///minimize \f[1/2*|| Ax - b ||_2^2 + \lambda || x ||_1 \f]
/// rho is the augmented Lagrangian parameter.

class Lasso{
public:
  Lasso(const Eigen::MatrixXd &A, const Eigen::VectorXd &b, double rho):
     z_(A.cols()),u_(A.cols()),
  lasso_max_iter_(200),
  lasso_abs_tol_(1.e-4),
  lasso_rel_tol_(1.e-2),
  A_(A),
  b_(b),
  ATb_(A_.transpose()*b_),
  n_(A.cols()),
  m_(A.rows()),
  rho_(rho){
    if(b_.rows() !=m_) throw std::runtime_error("b should have size m");
    lambda_max_=find_lambda_max();
    factor();
  }
  ///there is no solution for lambdas bigger than \f$\lambda_\max \f$
  double find_lambda_max() const;

  ///This is the main solver routine. The solution is stored in the vector x.
  /// alpha is the over-relaxation parameter (typical values for alpha are
  /// between 1.0 and 1.8).
  /// More information can be found in the paper linked at:
  /// http://www.stanford.edu/~boyd/papers/distr_opt_stat_learning_admm.html
  void lasso_solve(double lambda, double alpha, Eigen::VectorXd &x);
  ///compute the objective
  double objective(double lambda, Eigen::VectorXd &x) const;
  ///convergence norm between x and z
  double r_norm(const Eigen::VectorXd &x) const;
  ///convergence norm between present and old z
  double z_norm(const Eigen::VectorXd &zold) const;
  ///matrix norm the same way matlab has it
  static double matrix_norm(const Eigen::MatrixXd &M);
private:
  ///cholesky-decompose the matrix A with rho
  void factor();
  ///normalization needed in update of z
  void shrinkage(const Eigen::VectorXd &x, double kappa, Eigen::VectorXd &z);

  ///maximum lambda
  double lambda_max_;

  ///ADMM values
  Eigen::VectorXd z_, u_;

  ///Maximum number of LASSO iterations
  const int lasso_max_iter_;
  ///absolute convergence tolerance
  const double lasso_abs_tol_;
  ///relative convergence tolerance
  const double lasso_rel_tol_;

  ///reference to the matrix A
  const Eigen::MatrixXd &A_;
  ///reference to the vector b
  const Eigen::VectorXd &b_;
  ///A^T times b premultiplied
  const Eigen::VectorXd ATb_;

  ///number of rows of A
  const int n_;
  ///number of cols of A
  const int m_;

  ///LASSO augmented Lagrangian parameter rho
  const double rho_;

  ///Decomposition of matrix A
  Eigen::LLT<Eigen::MatrixXd> llt_;

};
