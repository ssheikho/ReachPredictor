#include "LinearAlgebra.h"

#include <iostream>

using namespace Eigen;
using namespace std;

std::pair<Eigen::MatrixXd, Eigen::MatrixXd>
	solveEigensystem(Eigen::MatrixXd a) {
	EigenSolver<MatrixXd> eigensolver(a);

	//cout << eigensolver.eigenvalues() << endl;
	//cout << eigensolver.eigenvectors() << endl;

	//SelfAdjointEigenSolver<MatrixXd> eigensolver(a);

	std::map<double, MatrixXd> retVal;

	for(size_t i = 0; i < a.rows(); i++)
		retVal.insert(
			pair<double, MatrixXd>(
				eigensolver.eigenvalues()(i).real()
				, eigensolver.eigenvectors().real().block(0, i, a.rows(), 1)));

	MatrixXd eigenValues(a.rows(), 1);
	MatrixXd eigenVectors(a.rows(), a.rows());
	size_t i = 0;
	for(std::map<double, Eigen::MatrixXd>::reverse_iterator iter =
			retVal.rbegin();
		iter != retVal.rend(); iter++) {
		eigenValues(i,0) = iter->first;
		eigenVectors.block(0,i,a.rows(),1) = iter->second;
		i++;
	}

	return std::pair<MatrixXd, MatrixXd>(eigenValues, eigenVectors);
}

//X is a 3 by n matrix in which  all the 3D coordinates of the 
//points are row-stacked.
// X_mean is the mean coordinate, i.e., the centroid of the points. 
MatrixXd buildCovMatPCA(MatrixXd inPts3D){
	size_t n = inPts3D.cols();
	MatrixXd covMat = MatrixXd::Zero(3,3);	
	covMat = (inPts3D.colwise() - inPts3D.rowwise().mean())
	* (inPts3D.colwise() - inPts3D.rowwise().mean()).transpose()
	/ n;

	std::pair<Eigen::MatrixXd, Eigen::MatrixXd> eigenSys 
	= solveEigensystem(covMat); 

	// diagEigMat: the diagonal matrix  with the 3  eigenvalues of covMat
	// OrthEigVMat:orthogonal matrix row-stacking the corresponding eigenvectors
	MatrixXd diagEigMat = MatrixXd::Identity(3,3);
	MatrixXd OrthEigVMat = MatrixXd::Identity(3,3);
	for (int i =0; i<3; i++)
		diagEigMat(i,i) = eigenSys.first(i,0);

	OrthEigVMat = eigenSys.second.transpose();

	MatrixXd diagCov = OrthEigVMat * diagEigMat * OrthEigVMat.transpose();

	//std::pair<Eigen::MatrixXd, Eigen::MatrixXd> eigSysDiagCov
	//= solveEigensystem(diagCov); 
	//cout << "diagEigMat:\n" << diagEigMat << endl;
	//cout << "OrthEigVMat:\n" << OrthEigVMat << endl;
	return diagCov;
}


