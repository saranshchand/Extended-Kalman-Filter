#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;
using std::cout;
using std::endl;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
   * TODO: Calculate the RMSE here.
   */
  VectorXd rmse(4);
  rmse << 0,0,0,0;

  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size
  if (estimations.size() != ground_truth.size()
      || estimations.size() == 0) {
    cout << "Invalid estimation or ground_truth data" << endl;
    return rmse;
  }

  // accumulate squared residuals
  for (unsigned int i=0; i < estimations.size(); ++i) {

    VectorXd residual = estimations[i] - ground_truth[i];

    // coefficient-wise multiplication
    residual = residual.array()*residual.array();
    rmse += residual;
  }

  // calculate the mean
  rmse = rmse/estimations.size();

  // calculate the squared root
  rmse = rmse.array().sqrt();

  // return the result
  return rmse;

}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
   * TODO:
   * Calculate a Jacobian here.
   */
  MatrixXd Hj(3, 4);
  
  double px, py, vx, vy;
  px = x_state(0);  
  py = x_state(1); 
  vx = x_state(2); 
  vy = x_state(3);   
  double px_2_plus_py_2 = pow(px, 2) + pow(py, 2);
  //checking division by zero
  if (px_2_plus_py_2 == 0)
  {
    cout << "Jacobian divide by zero error" << endl;
    return Hj;
  }
  
  //compute Jacobian
  else
  {
    Hj << px/sqrt(px_2_plus_py_2), py/sqrt(px_2_plus_py_2), 0, 0,
            -py/px_2_plus_py_2, px/px_2_plus_py_2, 0, 0,
            py*(vx*py - vy*px)/pow(sqrt(px_2_plus_py_2), 3), px*(vy*px - vx*py)/pow(sqrt(px_2_plus_py_2), 3), px/sqrt(px_2_plus_py_2), py/sqrt(px_2_plus_py_2);
  }
  
  return Hj;
}
