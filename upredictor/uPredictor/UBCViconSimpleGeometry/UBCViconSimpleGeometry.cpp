
#include "csv.h"
#include "ParseCSV.h"

#include "BallCentroidsOnRack.h"
#include <UBCUtil.h>

#include <Eigen/Dense>

#include <iostream>

using namespace Eigen;
using namespace std;

int main(int argc, char **argv) {

	string inLine = "2016_10_26_participant1_1_test.csv";
	BallCentroidsOnRack bcr(inLine);
	
	MatrixXd ballsCentroids = bcr.getBallCentroids();
	cout <<ballsCentroids<<endl;
	
	return 0;
}

