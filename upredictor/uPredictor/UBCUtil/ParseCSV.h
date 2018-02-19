#ifndef PARSE_CSV_H
#define PARSE_CSV_H

#include <Eigen/Core>
#include "csv.h"
#include <string>

using namespace Eigen;
using namespace std;

string getWholeFile(string fileName);
int countCols(string inLine);
int countRows(string inMat, string inCol);
void fillMat(int inColNo, string inLine, string inX,
	 string inY, string inZ, MatrixXd &mat);
MatrixXd parseCSV(string inMat);

#endif
