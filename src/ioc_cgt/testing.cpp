#include "ros/ros.h"
#include "std_msgs/String.h"
#include <eigen_matrix_utils/overloads.h>
#include <cppoptlib/problem.h>
#include <cppoptlib/solver/neldermeadsolver.h>
#include <Eigen/Dense>
#include <rosdyn_core/primitives.h>




template<typename T> class ObjFunc : public cppoptlib::Problem<T, 2> {
  public:
    using typename cppoptlib::Problem<T, 2>::TVector;
    
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    
    void setXi(Eigen::MatrixXd M) { xi = M;}
    void setU (Eigen::MatrixXd M) { u  = M;}
     
    T value(const TVector &x) {
      
      double res = 0;
      std::vector<double> lst;
      
      for (int j=0; j<xi.cols();j++)
      {
        lst.push_back(std::pow((x.transpose()*xi.col(j) + u.col(j)).norm(), 2.0));
      }
      
      for (int i=0;i<lst.size();i++)
        res += lst[i];
      
      return (res);
      
    }
    
  private:
    Eigen::MatrixXd xi;
    Eigen::MatrixXd u;
   
};


bool solveRiccati(const Eigen::MatrixXd &A,
                  const Eigen::MatrixXd &B,
                  const Eigen::MatrixXd &Q,
                  const Eigen::MatrixXd &R, Eigen::MatrixXd &P)
{

  const uint dim_x = A.rows();
  const uint dim_u = B.cols();

  Eigen::MatrixXd Ham = Eigen::MatrixXd::Zero(2 * dim_x, 2 * dim_x);
  Ham << A, -B * R.inverse() * B.transpose(), -Q, -A.transpose();

  Eigen::EigenSolver<Eigen::MatrixXd> Eigs(Ham);

  Eigen::MatrixXcd eigvec = Eigen::MatrixXcd::Zero(2 * dim_x, dim_x);
  int j = 0;
  for (int i = 0; i < 2 * dim_x; ++i) {
    if (Eigs.eigenvalues()[i].real() < 0.) {
      eigvec.col(j) = Eigs.eigenvectors().block(0, i, 2 * dim_x, 1);
      ++j;
    }
  }

  Eigen::MatrixXcd Vs_1, Vs_2;
  Vs_1 = eigvec.block(0, 0, dim_x, dim_x);
  Vs_2 = eigvec.block(dim_x, 0, dim_x, dim_x);
  P = (Vs_2 * Vs_1.inverse()).real();

  return true;
}

void fbn(const Eigen::MatrixXd &A,
         const Eigen::MatrixXd &B1,
         const Eigen::MatrixXd &B2,
         const Eigen::MatrixXd &Q1,
         const Eigen::MatrixXd &Q2,
         const Eigen::MatrixXd &R1,
         const Eigen::MatrixXd &R2, 
         Eigen::MatrixXd &K1,Eigen::MatrixXd &K2)
{
  Eigen::MatrixXd R12(R1.rows(),R1.cols());
  Eigen::MatrixXd R21(R2.rows(),R2.cols());
  
  R12.setZero();
  R21.setZero();
  
  Eigen::MatrixXd S1  = B1 * R1.inverse() * B1.transpose();
  Eigen::MatrixXd S2  = B2 * R2.inverse() * B2.transpose();
  Eigen::MatrixXd S12 = B1 * R1.inverse() * R21 * R1.inverse() * B1.transpose();
  Eigen::MatrixXd S21 = B2 * R2.inverse() * R12 * R2.inverse()* B2.transpose();

  solveRiccati(A,B1,Q1,R1,K1);
  solveRiccati(A,B2,Q2,R2,K2);
  
  Eigen::MatrixXd K1_prev = K1;
  Eigen::MatrixXd K2_prev = K2;
  double err_1 = 1;
  double err_2 = 1;
  double toll = 0.00001;

  while (err_1>toll && err_2>toll)
  {    
    Eigen::MatrixXd A1 = A - S2*K2;
    Eigen::MatrixXd A2 = A - S1*K1;
    
    Eigen::MatrixXd Q_1 = Q1 + K1*S21*K1;
    solveRiccati(A1,B1,Q_1,R1,K1);
    Eigen::MatrixXd Q_2 = Q2 + K2*S12*K2;
    solveRiccati(A2,B2,Q_2,R2,K2);
   
    err_1 = (K1-K1_prev).norm();
    err_2 = (K2-K2_prev).norm();
    
    K1_prev = K1;
    K2_prev = K2;
  }
  
  return;
}

Eigen::MatrixXd kron (const Eigen::MatrixXd &A, const Eigen::MatrixXd &B)
{
  Eigen::MatrixXd C(A.rows() * B.rows(), A.cols() * B.cols());
  C.setZero();

  for (int i = 0; i < A.rows(); i++)
  {
    for (int j = 0; j < A.cols(); j++)
    {
      C.block(i*B.rows(), j*B.cols(), B.rows(), B.cols()) = A(i, j) * B;
    }
  }
  return C;
}


bool kronsum( const Eigen::MatrixXd &A, Eigen::MatrixXd &B, Eigen::MatrixXd &C)
{
    if (A.cols() != A.rows())
      return false;
    if (B.cols() != B.rows())
      return false;
    
    C = kron(A, Eigen::MatrixXd::Identity(B.rows(), B.rows())) + kron(Eigen::MatrixXd::Identity(A.rows(), A.rows()), B);
    
    return true;
  }

int main(int argc, char **argv)
{
  
  typedef double T;
  typedef ObjFunc<T> ObjFunc;
  
  
  ros::init(argc, argv, "ioc_test");
  ros::NodeHandle n;
  
  double c = 2;
  double m = 6;
  
  Eigen::MatrixXd A(2,2);
  A << 0,1,
       0,-c/m;
  Eigen::MatrixXd Bh(2,1);
  Bh <<0,
       1/m;
  Eigen::MatrixXd Br(2,1);
  Br <<0,
       1/m;
  Eigen::MatrixXd Qhat(2,2);
  Qhat<<2,0,
        0,5;
  Eigen::MatrixXd Rhat(1,1);
  Rhat<<0.3;
  Eigen::MatrixXd Qr(2,2);
  Qr<<20,0,
        0,2;
  Eigen::MatrixXd Rr(1,1);
  Rr<<0.5;
  
  
  //FeedbackNashEquilibrium
  
  Eigen::MatrixXd Ph;
  Eigen::MatrixXd Pr;
  
  fbn(A,Bh,Br,Qhat,Qr,Rhat,Rr,Ph,Pr);
  
  ROS_FATAL_STREAM("P1\n"<<Ph);
  ROS_FATAL_STREAM("P2\n"<<Pr);
  
  
  // K_hat
  
  ROS_FATAL("yeyeye");

  ObjFunc f;

  ObjFunc::TVector x(2); x << 5, 5;

  Eigen::MatrixXd xi(2,10);
  xi <<  -0.1000   ,-0.0997   ,-0.0992   ,-0.0984   ,-0.0973   ,-0.0961   ,-0.0947   ,-0.0931   ,-0.0914   ,-0.0895,
          0.0057   , 0.0110    ,0.0159   , 0.0204    ,0.0246   , 0.0284    ,0.0318    ,0.0349    ,0.0377    ,0.0402;
  ROS_FATAL_STREAM("xi: \n"<<xi);
          
  Eigen::MatrixXd u(1,10);
  u <<  0.1014    ,0.0967    ,0.0920    ,0.0873    ,0.0827    ,0.0782    ,0.0738    ,0.0695    ,0.0653    ,0.0612;
  ROS_FATAL_STREAM("u: \n"<<u);
  
  f.setXi(xi);
  f.setU(u);
  cppoptlib::NelderMeadSolver<ObjFunc> solver;

  solver.minimize(f, x);

  std::cout << "   argmin: " << x.transpose() << std::endl;
  std::cout << "   f in argmin: " << f(x) << std::endl;
  
  
  Eigen::MatrixXd Khat = x.transpose();
  
  // IOC
  
  Eigen::MatrixXd Kr = Rr.inverse()*Br.transpose()*Pr;
  
  Eigen::MatrixXd F = A - Bh*Khat - Br*Kr;
  Eigen::MatrixXd F_kron;
  Eigen::MatrixXd Ft = F.transpose();
  
  if (!kronsum(Ft,Ft, F_kron))
    ROS_ERROR("kronsum failed");
  
  Eigen::MatrixXd Z1 = kron(Eigen::MatrixXd::Identity(Bh.rows(), Bh.rows()),Bh.transpose())*(F_kron.inverse());
  
  
  ROS_FATAL_STREAM("\n"<<F_kron);

  Eigen::MatrixXd K1_kron = kron(Khat.transpose(),Khat.transpose());
  Eigen::MatrixXd K2_kron = kron(Kr.transpose(),Kr.transpose());

  Eigen::MatrixXd M1_(Z1.rows(),6); //TODO:: dimensioni variabili
  M1_ << Z1, (Z1*K1_kron + kron(Khat.transpose(), Eigen::MatrixXd::Identity(Khat.rows(), Khat.rows()))), Z1*K2_kron;
  Eigen::MatrixXd M1(Z1.rows(),3);  //TODO dimensioni variabili 
  M1 << M1_.col(0),M1_.col(3),M1_.col(4) ;
  
  //   QP problem
  //   The problem is in the form:
  //   min 0.5 * x G x + g0 x
  //   s.t.
  //   CE^T x + ce0 = 0
  //   CI^T x + ci0 >= 0
 
  Eigen::MatrixXd H = (M1.transpose()*M1);
  Eigen::VectorXd l((H.cols()),1); l.setZero();

  Eigen::MatrixXd A_eq(H.rows(),H.rows()); A_eq.setZero(); A_eq(0,0) = 1;
  Eigen::VectorXd b_eq(H.rows());          b_eq.setZero(); b_eq(0) = -1;
  
  Eigen::MatrixXd A_ic = Eigen::MatrixXd::Identity(H.rows(),H.rows());
  Eigen::VectorXd b_ic(H.rows()); b_ic.setZero();
  
  Eigen::VectorXd sol(3);
  Eigen::solve_quadprog(H,l,A_eq,b_eq,A_ic,b_ic,sol );
  
  Eigen::Vector2d qq; qq<<sol(0), sol(1);
  
  Qhat = qq.asDiagonal();
  Rhat << sol(2);

  fbn(A,Bh,Br,Qhat,Qr,Rhat,Rr,Ph,Pr);
  Kr = Rr.inverse()*Br.transpose()*Pr;
  
  ROS_FATAL_STREAM("Kr\n"<<Kr);
  
  
  return 0;
}
