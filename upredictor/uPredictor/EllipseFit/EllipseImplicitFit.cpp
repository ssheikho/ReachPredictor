#include "EllipseImplicitFit.h"

EllipseImplicitFit::EllipseImplicitFit(MatrixXd inPtsAlongColsHom) :
	_inPtsAlongCols(inPtsAlongColsHom)
	, _plane(_inPtsAlongCols)
{
	Vector3d zAxis(3,1);
	zAxis << 0.0, 0.0, 1.0;
	_dotProduct = _plane.getPlaneNormalVectN().dot(zAxis);
	_crossProduct = _plane.getPlaneNormalVectN().cross(zAxis);

	_rotToXY =
		AngleAxisd(atan2(_crossProduct.norm(), _dotProduct)
			, _crossProduct.normalized()).toRotationMatrix();
	//assuming the last element of the hom coordinate is 1!
	_ptsXYAlongCols = _rotToXY * _inPtsAlongCols.block(0,0,3,
				_inPtsAlongCols.cols());
	_conicConstraintMat = buildConicConstraintMat(_ptsXYAlongCols);
	_conicVect = svdSolve(_conicConstraintMat);
	_conicImplicit = arrangeConicVectIntoImplicit(_conicVect);
}

EllipseImplicitFit::~EllipseImplicitFit() {}


Plane EllipseImplicitFit::getPlane() {
	return _plane;
}

MatrixXd EllipseImplicitFit::getRotToXY() {
	return _rotToXY;
}

MatrixXd EllipseImplicitFit::getPtsXYAlongCols() {
	return _ptsXYAlongCols;
}

MatrixXd EllipseImplicitFit::getConicConstraintMat() {
	return _conicConstraintMat;
}

MatrixXd EllipseImplicitFit::getConicVect() {
	return _conicVect;
}

MatrixXd EllipseImplicitFit::getConicImplicit() {
	return _conicImplicit;
}
