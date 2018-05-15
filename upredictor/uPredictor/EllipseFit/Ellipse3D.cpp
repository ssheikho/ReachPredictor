#include "Ellipse3D.h"
#include "Ellipse.h"
#include "EllipseImplicitFit.h"
#include "FitFunctions.h"

Ellipse3D::Ellipse3D(MatrixXd inPtsAlongColsHom) :
	_eif(inPtsAlongColsHom)
	, _ellipse(_eif) 
	, _totalXYCost(0.0)
	, _XYCosts(Eigen::MatrixXd::Zero
			(_eif.getPtsXYAlongCols().cols(),1)){}
 
Ellipse3D::~Ellipse3D() {}

EllipseImplicitFit Ellipse3D::getEIF() {
	return _eif;
}

Ellipse Ellipse3D::getEllipse() {
	return _ellipse;
}


MatrixXd Ellipse3D::getXYCostMat() {
	return _XYCosts;
}

double Ellipse3D::getTotalXYCost() {
	return _totalXYCost;
}

MatrixXd Ellipse3D::getPointAtThetaH(double theta) {
	return	cart3DRotMatToHomRT(_eif.getRotToXY().transpose()) *
			twoDHomToThreeDHom(_ellipse.getPointAtThetaH(theta)
				, -_eif.getPlane().getPlaneN()(3));
}

MatrixXd Ellipse3D::getPointAtThetasH(MatrixXd thetas) {
	// IMPORTANT: Transformation back from 2D->3D
	return	cart3DRotMatToHomRT(_eif.getRotToXY().transpose()) *
			twoDHomToThreeDHom(_ellipse.getPointAtThetasH(thetas)
				, -_eif.getPlane().getPlaneN()(3));
}

pair<double, double> Ellipse3D::findOneTheta(double startTheta, double x, double y) {
	pair<double, double> retVal;
	retVal.second = startTheta;
	ceres::CostFunction *cf =
		EllipseThetaMin::Create(
		_ellipse.getRT2D().getCX()
		, _ellipse.getRT2D().getCY()
		, _ellipse.getRT2D().getTheta()
		, _ellipse.getA()
		, _ellipse.getB()
		, x, y);

	ceres::Problem problem;
	problem.AddResidualBlock(cf,
		new ceres::CauchyLoss(0.5) /* squared loss */,
		&retVal.second);

	ceres::Solver::Options options;
	ceres::Solver::Summary summary;

	ceres::Solve(options, &problem, &summary);
	retVal.first = summary.final_cost;
	return retVal;
}

pair<double, double> Ellipse3D::findOneThetaQuad(
	double &startTheta, double x, double y) {
	double incr = M_PI / 2.0;

	pair<double, double> a = findOneTheta(startTheta, x, y);
	pair<double, double> b = findOneTheta(startTheta + incr, x, y);
	pair<double, double> c = findOneTheta(startTheta + incr * 2.0, x, y);
	pair<double, double> d = findOneTheta(startTheta + incr * 3.0, x, y);

	a = a.first < b.first ? a : b;
	c = c.first < d.first ? c : d;

	return a.first < c.first ? a : c;
}

MatrixXd Ellipse3D::findThetas() {
	//getPtsXYAlongCols in cartesian 
	MatrixXd pts =_eif.getPtsXYAlongCols();
	size_t cols = pts.cols();
	MatrixXd thetas(1,cols);
	thetas(0,0) = (double)rand() / (double)RAND_MAX * M_PI * 2.0 - M_PI;


	pair<double, double> a =
		findOneThetaQuad(thetas(0,0), pts(0,0), pts(1,0));
	thetas(0,0) = a.second;
	for(size_t i = 1; i < cols; i++) {
		double startTheta = thetas(0,i-1);
		a = findOneThetaQuad(startTheta, pts(0,i), pts(1,i));
		thetas(0,i) = a.second;
		_XYCosts(i,0) = a.first;
//cout<<XYCosts(0,i)<<endl;
	}
	_totalXYCost = _XYCosts.sum();
	return thetas;
}
