
#include "ParseMathematica.h"
#include "ParseCSV.h"
#include "UBCUtil.h"
#include "FitFunctions.h"
#include "Ellipse3D.h"


#include "LinearAlgebra.cpp"
#include "RigidTrans2D.h"

#include "BaysFitFunctions.h"

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

int main(int argc, char **argv) {
	srand(time(NULL));

	//MatrixXd inPtsAlongColsHom = parseMathematica(getWholeFile(argv[1]));
	//cout << "inPtsAlongColsHom :\n" << inPtsAlongColsHom.rows()<< endl;

	int colNo = 0;
	int nRows = countRowsFColCSV(argv[1]);
	// 3xn
	MatrixXd inWPtsAlongCols(3,nRows);
	fillMatWholeCSV(colNo, argv[1], inWPtsAlongCols);
	//cout << "inWPtsAlongCols :\n" << inWPtsAlongCols<< endl;

	/*********** Plane Fitting ***********/
	//cout << "/*********** Plane Fitting ***********/\n"
	//vector<std::pair<Eigen::MatrixXd, Eigen::MatrixXd>>* eigenSysS;
  MatrixXd covMat = buildCovMatPCA(inWPtsAlongCols);
	std::pair<Eigen::MatrixXd, Eigen::MatrixXd> eigenSys 
		= solveEigensystem(covMat); 		
	//cout << eigenSys.first << "\n" << endl;
	//cout << "EigenVals for " << argv[1] <<":\n" << eigenSys.first<< endl;
	/*********** End - Plane Fitting ***********/

	/*********** Ellipse Fitting ***********/
	Ellipse3D e3D(cartToHom(inWPtsAlongCols));
	MatrixXd thetas = e3D.findThetas();
	double totalXYCost = e3D.getTotalXYCost();
	MatrixXd fitCosts = e3D.getXYCostMat();
	MatrixXd pts2d = e3D.getEllipse().getPointAtThetasH(thetas);
	MatrixXd pts3d = e3D.getPointAtThetasH(thetas);
	
	// 3xn
	MatrixXd inPts2dRotatedAlongCols  = e3D.getEIF().getPtsXYAlongCols();

	Ellipse ellipse = e3D.getEllipse();
	double cX = ellipse.getRT2D().getCX();
	double cY = ellipse.getRT2D().getCY();
	double alpha = ellipse.getRT2D().getTheta();
	double a = ellipse.getA();
	double b = ellipse.getB();
	double EfAxisRatio = b/a;

	//printEigenMathematica(e3D.getEIF().getPtsXYAlongCols(), cout, "inPtsRotated");
	//cout << "Ellipse-Fit-Costs:\n" << inPtsRotated << endl;
	
	//cout << "cX:" << cX << endl;
	//cout << "cY:" << cY << endl;
	//cout << "alpha:" << alpha << endl;
	//cout << "a:" << a << endl;
	//cout << "b:" << b << endl;	
	if(totalXYCost >= 0.0) {
		//cout << eigenSys.first(2,0)<< endl;
		//cout << /*"Axis-Ratio:" <<*/ EfAxisRatio << endl;	
		//cout << /*"Ellipse-Total-Fit-Cost:" <<*/ totalXYCost << endl;
		//cout << "Ellipse-Fit-Costs:\n" << fitCosts << endl;
	
	

	/*********** Bayes Fitting ***********/

		int n = inPts2dRotatedAlongCols.cols();
		MatrixXd xVec = (inPts2dRotatedAlongCols.block(0,0,1,n)).transpose();
		MatrixXd yVec = (inPts2dRotatedAlongCols.block(1,0,1,n)).transpose();	
		int maxPolyOrder = 5;
		MatrixXd evidenceMat(maxPolyOrder+1,1);		
		MatrixXd AMetaLi = MatrixXd::Zero(maxPolyOrder+1,maxPolyOrder+1);	
		MatrixXd AMetaEv = MatrixXd::Zero(maxPolyOrder+1,maxPolyOrder+1);
		MatrixXd likelihoodMat(maxPolyOrder+1,1);
		Eigen::MatrixXd yVecHatMetaLi	= 
			MatrixXd::Zero(maxPolyOrder+1,n);
		Eigen::MatrixXd yVecHatMetaEv	= 
			MatrixXd::Zero(maxPolyOrder+1,n);
size_t k = 4;
		//for(size_t k = 0; k < maxPolyOrder; k++) {
			Eigen::MatrixXd Xmat = MatrixXd::Zero(k+1,n);
			Xmat = buildXMatIn(xVec, k);
			
			vector<double> AvecLi = findALi(Xmat, yVec.transpose());	
		
			for(size_t j = 0; j < AvecLi.size(); j++) 
				AMetaLi(k,j) = AvecLi[j];

		
		Eigen::MatrixXd yVecHatLi = 
				predictY(Xmat, AMetaLi.block(k,0,1,k+1));
		yVecHatMetaLi.block(k,0,1,n) = yVecHatLi;

		
		evidenceMat(k,0) = 
				computePy_xvalpha(Xmat, yVec.transpose());
		likelihoodMat(k,0) = 
				computePy_xav(Xmat, yVec.transpose());
	//}

	
	//cout << evidenceMat(k,0) << endl;	
	cout << /*"likelihood for " << argv[1] <<":\n" <<*/ likelihoodMat(k,0) << endl;

	/*********** End - Bayes Fitting ***********/

/*
	MatrixXd fixedThetas = fixThetas(thetas);
	MatrixXd speed = simpleGradientAcrossCols(fixedThetas);
	MatrixXd acceleration = simpleGradientAcrossCols(speed);
	MatrixXd jerk = simpleGradientAcrossCols(acceleration);*/


	/*printEigenMathematica(inWPtsAlongCols, cout, "inPts");
	printEigenMathematica(e3D.getEIF().getPtsXYAlongCols(), cout, "inPtsRotated");
	printEigenMathematica(pts2d, cout, "outPts2d");
	printEigenMathematica(pts3d, cout, "outPts3d");
	printEigenMathematica(fixedThetas, cout, "Thetas");
	printEigenMathematica(speed, cout, "speed");
	printEigenMathematica(acceleration, cout, "acceleration");
	printEigenMathematica(jerk, cout, "jerk");*/

	/*********** End - Ellipse Fitting ***********/

	
	return 0;
}
	return 0;
}
