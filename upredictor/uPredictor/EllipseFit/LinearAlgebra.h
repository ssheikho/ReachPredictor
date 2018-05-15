#ifndef LINEAR_ALGEBRA_H
#define LINEAR_ALGEBRA_H

#include <Eigen/Dense>

#include <map>

std::pair<Eigen::MatrixXd, Eigen::MatrixXd> solveEigensystem(Eigen::MatrixXd a);
Eigen::MatrixXd buildCovMatPCA(Eigen::MatrixXd inPts3D);
#endif
