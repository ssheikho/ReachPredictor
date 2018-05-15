#ifndef ELLIPSE_H
#define ELLIPSE_H


#include "FitFunctions.h"
#include "LinearAlgebra.h"
#include "EllipseImplicitFit.h"

#include <Eigen/Dense>
#include <cmath>

#include "RigidTrans2D.h"

using namespace Eigen;

//template<typename T>
class Ellipse {
public:
	Ellipse(double a, double b, double alpha, double cX, double cY);

	Ellipse();
	Ellipse(EllipseImplicitFit eif);

	~Ellipse();
	
	void setAB(double a, double b);

	double getA();

	double getB();

	RigidTrans2D getRT2D();

	MatrixXd getPointAtThetaH(double theta);

	double getPointXAtThetaH(double theta);

	double getPointYAtThetaH(double theta);

	double getPointWAtThetaH(double theta);


	MatrixXd getPointAtThetasH(MatrixXd thetas);

protected:
	double _a, _b;
	RigidTrans2D _rt2D;

};
#endif
