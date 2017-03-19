/*
 * Copyright (C) 1998-2016 ALPS Collaboration
 * 
 *     This program is free software; you can redistribute it and/or modify it
 *     under the terms of the GNU General Public License as published by the Free
 *     Software Foundation; either version 2 of the License, or (at your option)
 *     any later version.
 * 
 *     This program is distributed in the hope that it will be useful, but WITHOUT
 *     ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *     FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for 
 *     more details.
 * 
 *     You should have received a copy of the GNU General Public License along
 *     with this program; if not, write to the Free Software Foundation, Inc., 59
 *     Temple Place, Suite 330, Boston, MA 02111-1307 USA
 * For use in publications, see ACKNOWLEDGE.TXT
 */

#include "../src/SpM_Lasso.hpp"
#include "gtest.h"
#include <fstream>

class lassoTest:public ::testing::Test{
public:
protected:
  virtual void SetUp() {
    int cols=40; //number of columns of A
    int rows=20; //number of rows of A
    
    A.resize(rows, cols);
    b.resize(rows);
  }
  virtual void TearDown() {
  }
  lassoTest(){}
  Eigen::MatrixXd A;
  Eigen::VectorXd b;

  void load_A(std::string &filename, int rows, int cols);
  void load_b(std::string &filename, int cold);
};

void lassoTest::load_A(std::string &filename, int rows, int cols){
  std::ifstream A_file(filename);
  if(!A_file.good()) throw std::runtime_error("problem opening file: "+filename);
  A.resize(rows, cols);
  for(int i=0;i<rows;++i){
    for(int j=0;j<cols;++j){
      A_file>>A(i,j);
    }
  }
  std::cout<<A<<std::endl;
}
void lassoTest::load_b(std::string &filename, int rows){
  std::ifstream b_file(filename);
  if(!b_file.good()) throw std::runtime_error("problem opening file: "+filename);
  b.resize(rows);
    for(int j=0;j<rows;++j){
      b_file>>b(j);
  }
}

TEST_F(lassoTest, Init) {
  double rho=1.;
  Lasso L(A, b, rho);
}
TEST_F(lassoTest, FileIO) {
  int rows=3;
  int cols=4;
  std::string filename_A="../test/TestFiles/AFile_3_4.txt";
  std::string filename_b="../test/TestFiles/bFile_3.txt";
  load_A(filename_A, rows, cols);
  load_b(filename_b, rows);
  double rho=1.;
  double alpha=1.2;
  double lambda=1.;
  
  EXPECT_NEAR(1.69935438625108, Lasso::matrix_norm(A), 1.e-10);
  EXPECT_NEAR(1.70480451723583, b.norm(), 1.e-10);
}
/*TEST_F(lassoTest, LassoSolve_15_50) {
  int rows=15;
  int cols=50;
  std::string filename_A="../test/TestFiles/AFile_15_50.txt";
  std::string filename_b="../test/TestFiles/bFile_15.txt";
  load_A(filename_A, rows, cols);
  load_b(filename_b, rows);
  double rho=1.;
  double alpha=1.2;
  double lambda=1.;
  Lasso L(A, b, rho);

  Eigen::VectorXd x;
  L.lasso_solve(lambda, alpha, x);
}*/
TEST_F(lassoTest, LassoSolve_3_4) {
  int rows=3;
  int cols=4;
  std::string filename_A="../test/TestFiles/AFile_3_4.txt";
  std::string filename_b="../test/TestFiles/bFile_3.txt";
  load_A(filename_A, rows, cols);
  load_b(filename_b, rows);
  double rho=1.;
  double alpha=1;
  Lasso L(A, b, rho);
  EXPECT_NEAR(1.70451721668332, L.find_lambda_max(), 1.e-10);
  double lambda=L.find_lambda_max()/10.;

  Eigen::VectorXd x;
  L.lasso_solve(lambda, alpha, x);
}

