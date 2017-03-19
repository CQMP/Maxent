#include"maxent_matrix_def.hpp"

///Class for performing LASSO fitting using ADMM
///minimize \f[1/2*|| Ax - b ||_2^2 + \lambda || x ||_1 \f]

class Lasso{
public:
  Lasso(const Eigen::MatrixXd &A, const Eigen::VectorXd &b):
  lasso_max_iter_(1000),
  lasso_abs_tol_(1.e-4),
  lasso_rel_tol_(1.e-2),
  A_(A),
  b_(b),
  n_(A.cols()),
  m_(A.rows()){
    if(b_.rows() !=m_) throw std::runtime_error("b should have size m");
  }
private:
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

  ///number of rows of A
  const int n_;
  ///number of cols of A
  const int m_;
};
