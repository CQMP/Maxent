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
#include "boost/random.hpp"

class lassoTest:public ::testing::Test{
public:
  typedef boost::variate_generator< boost::mt19937&, boost::normal_distribution<> > dist_n;
  typedef boost::variate_generator< boost::mt19937&, boost::uniform_real<> > dist_u;
protected:
  virtual void SetUp() {
    int n=40; //number of columns of A
    int m=20; //number of rows of A
    double p=1./n;
    
    //boost random crap
    boost::mt19937 rng;
    boost::normal_distribution<> nd(0.0, 1.0);
    boost::uniform_real<> ud(0.0, 1.0);
    dist_n var_nor(rng, nd);
    dist_u var_uni(rng, ud);
    rng.seed(0);
    
    //initialization as in LASSO sample
    initialize_A(var_nor, n, m);
    initialize_b(var_nor, var_uni, n, m, p);

    L=new Lasso(A, b);
  }
  virtual void TearDown() {
    delete L;
  }
  lassoTest(){}
  Eigen::MatrixXd A;
  Eigen::VectorXd b;
  Lasso *L;
private:
  boost::mt19937 rng;
  void initialize_A(dist_n &var_nor, int n, int m);
  void initialize_b(dist_n &var_nor, dist_u &var_uni, int n, int m, double p);
};

void lassoTest::initialize_A(dist_n &var_nor, int n, int m){
/*  A = randn(m,n);
    A = A*spdiags(1./sqrt(sum(A.^2))',0,n,n); % normalize columns
 */ 
  A.resize(m,n);
  for(int i=0;i<A.cols();++i){
    for(int j=0;j<A.rows();++j){
      A(j,i)=var_nor();
    }
  }
  //now normalize columns
  for(int j=0;j<A.cols();++j){
    double csum=0;
    for(int i=0;i<A.rows();++i){
      csum+=A(i,j)*A(i,j);
    }
    for(int i=0;i<A.rows();++i){
      A(i,j)=A(i,j)/std::sqrt(csum);
    }
  }
}
void lassoTest::initialize_b(dist_n &var_nor, dist_u &var_uni, int n, int m, double p){
  /*x0 = sprandn(n,1,p);
  b = A*x0 + sqrt(0.001)*randn(m,1);*/
  b.resize(m);
  Eigen::VectorXd x0(n);
  for(int i=0;i<n;++i){
    if(var_uni()<p){
      x0[i]=var_nor();
    }else{
      x0[i]=0.;
    }
  }
  b=A*x0;
  for(int i=0;i<m;++i){
    b[i]+=std::sqrt(0.001)*var_nor();
  }
}

TEST_F(lassoTest, Init) {
  ;//do nothing
}
