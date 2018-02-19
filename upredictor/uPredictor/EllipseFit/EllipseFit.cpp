#include "Ellipse3D.h"
#include "ParseMathematica.h"
#include <UBCUtil.h>

/*
#include "FitFunctions.h"
#include "LinearAlgebra.h"
#include "RigidTrans2D.h"
*/

#include <Eigen/Dense>
//#include <Eigen/SVD>

//#include <cmath>
#include <iostream>

using namespace Eigen;
using namespace std;

MatrixXd fixThetas(MatrixXd a) {
	MatrixXd retVal(a.rows(), a.cols());

	double curTheta = a(0,0);
	while(curTheta >= -M_PI) curTheta -= M_PI;
	while(curTheta <= M_PI) curTheta += M_PI;
	retVal(0,0) = curTheta;

	for(size_t i = 1; i < a.cols(); i++) {
		double curTheta = a(0,i);
		while(curTheta >= -M_PI) curTheta -= 2.0 * M_PI;
		while(curTheta <= M_PI) curTheta += 2.0 * M_PI;

		double testA = curTheta + 2.0 * M_PI;
		double testB = curTheta - 2.0 * M_PI;
		
		double distNull = curTheta - retVal(0, i - 1);
		double distA = testA - retVal(0, i - 1);
		double distB = testB - retVal(0, i - 1);

		retVal(0, i) = distNull < distA ? curTheta : testA;
		retVal(0, i) = distNull < distB ? curTheta : testB;
	}

	return retVal;
}

MatrixXd simpleGradientAcrossCols(MatrixXd a) {
	MatrixXd retVal(a.rows(), a.cols() - 1);

	for(size_t i = 1; i < a.cols(); i++) {
		for(size_t j = 0; j < a.rows(); j++) {
			retVal(j, i - 1) = a(j, i) - a(j, i - 1);
		}
	}

	return retVal;
}

int main(int argc, char **argv) {
	srand(time(NULL));
	MatrixXd inPtsAlongColsHom = parseMathematica(getWholeFile(argv[1]));

	//MatrixXd inPtsAlongCols = inPtsAlongRows.transpose();

	//printEigenMathematica(inPtsAlongCols, cout, "inPtsAlongCols");

	Ellipse3D<double> e3D(inPtsAlongColsHom);
	MatrixXd thetas = e3D.findThetas();
	double totalXYCost = e3D.getTotalXYCost();
	MatrixXd pts2d = e3D.getEllipse().getPointAtThetasH(thetas);
	MatrixXd pts3d = e3D.getPointAtThetasH(thetas);

	//Added by Sara
	Ellipse ellipse = e3D.getEllipse();
	double cX = ellipse.getRT2D().getCX();
	double cY = ellipse.getRT2D().getCY();
	double alpha = ellipse.getRT2D().getTheta();
	double a = ellipse.getA();
	double b = ellipse.getB();
	

	cout << "cX:" << cX << endl;
	cout << "cY:" << cY << endl;
	cout << "alpha:" << alpha << endl;
	cout << "a:" << a << endl;
	cout << "b:" << b << endl;
	cout << "totalXYCost:" << totalXYCost << endl;
	

	MatrixXd fixedThetas = fixThetas(thetas);
	MatrixXd speed = simpleGradientAcrossCols(fixedThetas);
	MatrixXd acceleration = simpleGradientAcrossCols(speed);
	MatrixXd jerk = simpleGradientAcrossCols(acceleration);

	printEigenMathematica(e3D.getEIF().getPtsXYAlongCols(), cout, "inPts");
	printEigenMathematica(pts2d, cout, "outPts2d");
	printEigenMathematica(pts3d, cout, "outPts3d");
	printEigenMathematica(fixedThetas, cout, "fixedThetas");
	printEigenMathematica(speed, cout, "speed");
	printEigenMathematica(acceleration, cout, "acceleration");
	printEigenMathematica(jerk, cout, "jerk");

	return 0;
}
