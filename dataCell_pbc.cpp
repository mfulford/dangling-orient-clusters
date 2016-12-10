#include <cstdlib>      
#include <iostream>     
#include <algorithm>    
#include <numeric>		// accumulate (for sum of vector elements)
#include <math.h>       // sqrt, pow, round
#include <vector>
#include "dataCell_pbc.hpp"

using namespace std;

float pbc_x;     //Define it once for global use. Cell dimensions from cell data
float pbc_y;
float pbc_z;

// Open cell data file
dataCell::dataCell(const char cellfilename[]) {
	cout << "opening cell data file" << endl;
	cellf.open(cellfilename);
	cellf >> pbc_x >> pbc_y >> pbc_z;  	//read in cell data	
	cellf.close();
}

pbc_check::pbc_check() {
}

// Return true (1) if phase is prism, otherwise false (0)
bool pbc_check::IsPhasePrism() {
	if (pbc_x > pbc_y) {
		cout << "Prism Ice! " << endl;
		return true;
	} else if (pbc_y > pbc_x) {
		cout << "Basal Ice! " << endl;
		return false;
	}
}

// Determine box center
float pbc_check::box_centre(int NATOM, const float* coord_norm) {

	box_middle = *min_element(coord_norm, coord_norm + NATOM)
			+ (*max_element(coord_norm, coord_norm + NATOM)
					- *min_element(coord_norm, coord_norm + NATOM)) / 2.0;
	cout << "Box centre: " << box_middle << endl;
	return box_middle;
}

// Average theta (orientation) and q3 of dangling molecules:   
void Stats::average(const vector<vector<int> > &DanglingTop,
		const vector<vector<int> > &DanglingBot,
		const vector<vector<float> > &orient, const vector<float> &q3or) {

	avOrient_Top = 0.0;
	avOrient_Bot = 0.0;
	avQ3_Top = 0.0;
	avQ3_Bot = 0.0;

	// Determine number dangling OH on each surface:
	int numTop = DanglingTop[0].size() + DanglingTop[1].size();
	int numBot = DanglingBot[0].size() + DanglingBot[1].size();

	// Determine the average orientation and q3 of dangling OH on each surface:
	for (int H = 0; H < 2; H++) {

		for (int i = 0; i < DanglingTop[H].size(); i++) {
			avOrient_Top = orient[0][DanglingTop[H][i]] + avOrient_Top;
			avQ3_Top = q3or[DanglingTop[H][i]] + avQ3_Top;
		}
		for (int i = 0; i < DanglingBot[H].size(); i++) {
			avOrient_Bot = orient[0][DanglingBot[H][i]] + avOrient_Bot;
			avQ3_Bot = q3or[DanglingBot[H][i]] + avQ3_Bot;
		}
	}

	// Average q3 and  orientation of dangling molecules:
	avOrient_Top = avOrient_Top / numTop;
	avOrient_Bot = avOrient_Bot / numBot;
	avQ3_Top = avQ3_Top / numTop;
	avQ3_Bot = avQ3_Bot / numBot;

	cout << "#ORIENTDANG Top Layer Average: " << avOrient_Top
			<< " Bottom Layer: " << avOrient_Bot << endl;
	cout << "#Q3DANG Top Layer Average: " << avQ3_Top << " Bottom Layer: "
			<< avQ3_Bot << endl;

}

// Average q3 and orientation of QLL molecules:   
void Stats::average(const vector<int> &QLLTop, const vector<int> &QLLBot,
		const vector<vector<float> > &orient, const vector<float> &q3or) {

	avOrient_Top = 0.0;
	avOrient_Bot = 0.0;
	avQ3_Top = 0.0;
	avQ3_Bot = 0.0;

	// Determine number dangling OH on each surface:
	int numTop = QLLTop.size();
	int numBot = QLLBot.size();

	// Determine the average theta of QLL OH on each surface:
	for (int i = 0; i < QLLTop.size(); i++) {
		avOrient_Top = orient[0][QLLTop[i]] + avOrient_Top;
	}
	for (int i = 0; i < QLLBot.size(); i++) {
		avOrient_Bot = orient[0][QLLBot[i]] + avOrient_Bot;
	}

	for (int i = 0; i < QLLTop.size(); i++) {
		avQ3_Top = q3or[QLLTop[i]] + avQ3_Top;
	}
	for (int i = 0; i < QLLBot.size(); i++) {
		avQ3_Bot = q3or[QLLBot[i]] + avQ3_Bot;
	}

	// average q3 and orientation of dangling molecules:
	avOrient_Top = avOrient_Top / numTop;
	avOrient_Bot = avOrient_Bot / numBot;
	avQ3_Top = avQ3_Top / numTop;
	avQ3_Bot = avQ3_Bot / numBot;

	cout << "#ORIENTQLL Top Layer Average: " << avOrient_Top
			<< " Bottom Layer: " << avOrient_Bot << endl;
	cout << "#Q3QLL Top Layer Average: " << avQ3_Top << " Bottom Layer: "
			<< avQ3_Bot << endl;

}

float pbc_check::pbc_xcoord(const float xx) {
	xcoord = xx - (pbc_x * round(xx / pbc_x));
	distx_sqr = pow(xcoord, 2.0);
}

float pbc_check::pbc_ycoord(const float yy) {
	ycoord = yy - (pbc_y * round(yy / pbc_y));
	disty_sqr = pow(ycoord, 2.0);
}

float pbc_check::pbc_zcoord(const float zz) {
	zcoord = zz - (pbc_z * round(zz / pbc_z));
	distz_sqr = pow(zcoord, 2.0);
}

float pbc_check::getdistx() const {
	return distx_sqr;
}
float pbc_check::getdisty() const {
	return disty_sqr;
}
float pbc_check::getdistz() const {
	return distz_sqr;
}

float pbc_check::pbc_xyz(const float xx, const float yy, const float zz) {
	x_bound = xx - (pbc_x * round(xx / pbc_x)); // round() gives nearest integral 3.8 --> 4.0.
	y_bound = yy - (pbc_y * round(yy / pbc_y));
	z_bound = zz - (pbc_z * round(zz / pbc_z));
	dist = sqrt(pow(x_bound, 2.0) + pow(y_bound, 2.0) + pow(z_bound, 2.0)); // pow squares
}

vector<float> pbc_check::h_midpoint(const float x1, const float y1,
		const float z1, const float x2, const float y2, const float z2) {

	vector<float> coord_mid(3);

	float x_dist = abs(x1 - x2);
	float y_dist = abs(y1 - y2);
	float z_dist = abs(z1 - z2);

	float x_mid = (x1 + x2) / 2.0;
	float y_mid = (y1 + y2) / 2.0;
	float z_mid = (z1 + z2) / 2.0;

	if (x_dist >= pbc_x / 2) {
		x_mid = x_mid + (pbc_x / 2.0);
	}  
	if (y_dist >= pbc_y / 2) {
		y_mid = y_mid + (pbc_y / 2.0);
	}
	if (z_dist >= pbc_z / 2) {
		z_mid = z_mid + (pbc_z / 2.0);
	}

	coord_mid = {x_mid, y_mid, z_mid};
	return coord_mid;
}

float pbc_check::getdist_pbc() const {
	return dist;
}

float pbc_check::getxpbc() const {
	return x_bound;
}

float pbc_check::getypbc() const {
	return y_bound;
}

float pbc_check::getzpbc() const {
	return z_bound;
}

void pbc_check::clear_data() {

}