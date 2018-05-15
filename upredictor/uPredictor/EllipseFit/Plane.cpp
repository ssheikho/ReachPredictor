#include "Plane.h"
#include "UBCUtil.h"

Plane::Plane(MatrixXd inPts) :
	_plane(svdSolve((inPts).transpose()))
	, _planeNormalVect(_plane.block(0,0,3,1))
	, _norm(_planeNormalVect.norm())
	, _planeN(_plane / _norm)
	, _planeNormalVectN(_planeN.block(0,0,3,1)) {}

Plane::~Plane() {}

Vector3d Plane::getPlaneCentroid(MatrixXd inPts) {
	size_t cols = inPts.cols();

	Vector3d planeCentroid(Vector3d::Zero(3,1));
	
	for(size_t i = 0; i < cols; i++) {
		planeCentroid(0,0) +=  inPts(0,i);
		planeCentroid(1,0) +=  inPts(1,i);
		planeCentroid(2,0) +=  inPts(2,i);
	}
	//planeCentroid(0,0) = inPts.block(0, 0, 1, inPts.cols()).sum();
	//planeCentroid(1,0) = inPts.block(1, 0, 1, inPts.cols()).sum();
	//planeCentroid(2,0) = inPts.block(2, 0, 1, inPts.cols()).sum();

	return (planeCentroid / double(cols));
}

MatrixXd Plane::centerInPts(MatrixXd inPts){
	size_t cols = inPts.cols();
	MatrixXd inPtsC(MatrixXd::Zero(inPts.rows(),cols));
	inPtsC.block(inPts.rows(), 0, 1, inPts.cols()) =
		Eigen::MatrixXd::Constant(1, inPts.cols(), 1.0);

	Vector3d planeCentroid = getPlaneCentroid(inPts);

	for(size_t i = 0; i < cols; i++) 
		inPtsC.block(0,i,3,1) = inPts.block(0,i,3,1) - planeCentroid;

	return inPtsC;
}

Vector4d Plane::getPlane() {
	return _plane;
}

Vector3d Plane::getPlaneNormalVect() {
	return _planeNormalVect;
}

double Plane::getNorm() {
	return _norm;
}

Vector4d Plane::getPlaneN() {
	return _planeN;
}

Vector3d Plane::getPlaneNormalVectN() {
	return _planeNormalVectN;
}





