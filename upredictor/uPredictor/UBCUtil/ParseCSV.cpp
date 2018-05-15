#include "ParseCSV.h"

#include <fstream>

string getWholeFile(string fileName) {
	ifstream infile(fileName);
	string s, out = "";
	while (!infile.eof()) {	
		getline(infile, s);
		out += s;
	}
	infile.close();
	return out;
}

int countCols(string inLine) {
	inLine = inLine.substr(1);
	int commaDelim = inLine.find(',');
	int endBracketDelim = inLine.find('}');
	int ctr = 0;
	while((commaDelim >= 0) && (commaDelim < endBracketDelim)) {
		//cout << inLine.substr(0, commaDelim) << endl;		
		inLine = inLine.substr(commaDelim + 1);
		commaDelim = inLine.find(',');
		endBracketDelim = inLine.find('}');
		ctr++;
	}
	ctr++;
	//cout << inLine.substr(0, endBracketDelim) << endl;
	return ctr;
}

int countRowsCSV(string inMat, string inCol) {

	io::CSVReader<1> in(inMat);
	in.read_header(io::ignore_extra_column, inCol);
	int ctr = 0;
	while(in.read_row(inCol)){
		ctr++;

	}
	return ctr - 1;

}

int countRowsFColCSV(string inMat) {

	io::CSVReader<3> in(inMat);
	//in.read_header(io::ignore_extra_column, inCol);
	int ctr = 0 ;
	double x, y, z;
	while(in.read_row(x, y, z)) ctr++;
	return ctr;
}

void fillMatCSV(int inColNo, string inLine, string inX, 
	string inY, string inZ, MatrixXd &mat) {

	int nRows = countRowsCSV(inLine, inX);
	int nCols = 3;

	io::CSVReader<3> in(inLine);
	in.read_header(io::ignore_extra_column, inX, inY, inZ);
	double inX_, inY_, inZ_;

	int rowNo = 0;
	while(in.read_row(inX_, inY_, inZ_) && rowNo < mat.rows()){	
		// inX_, inY_, inZ_: contain the value from the file
		mat(rowNo, inColNo) = inX_;
		mat(rowNo, inColNo+1) = inY_;
		mat(rowNo, inColNo+2) = inZ_;

		rowNo++;
		inColNo + 2;
	}
}

void fillMatWholeCSV(int inColNo, string inLine, MatrixXd &mat) {

	int nRows = countRowsFColCSV(inLine);
	int nCols = 3;

 	io::CSVReader<3/*, trim_chars<' '>, double_quote_escape<',','\"'>*/ >in(inLine);
	//in.read_header(io::ignore_extra_column, inX, inY, inZ);
	double inX_, inY_, inZ_;

	int rowNo = 0;
	while(in.read_row(inX_, inY_, inZ_) && rowNo < mat.cols()){	
		mat(inColNo, rowNo) = inX_;
		mat(inColNo+1, rowNo) = inY_;
		mat(inColNo+2, rowNo) = inZ_;

		rowNo++;
		inColNo + 2;
	}
}



