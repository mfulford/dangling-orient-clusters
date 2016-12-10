#include <iostream>     	
#include <vector>
#include <ctime> 
#include <algorithm>	 	// max vector size
#include <math.h>       	// isnan, sqrt 
#include <sstream>			// istringstream
#include "q3_r.hpp"
#include "dataCell_pbc.hpp"

// In this class the q3 file is read. Liquid/ice like molecules identified and indexed from q3. 

using namespace std;

q3_r::q3_r(const char q3filename[], const int nox, const float q3thresh) {
	cout << "opening q3 data file" << endl;
	num_ox = nox;
	threshold = q3thresh;
	q3f.open(q3filename);
}

void q3_r::q3_frame() {

	q3or.clear();

	string q3_str;
	float q3_float; 

	//here we read in q3 frame by frame. 
	for (int i = 0; i < num_ox; i++) {
		q3f >> q3_str;

		if (q3_str == "NaN") {
			cout << "#### q3 is NaN. i: " << i << endl;
			cout << "#### Instead q3_tmp is now =0.000" << endl;
			q3_float = 0.000;
		} else {
			q3_float = stof(q3_str);
		}

		q3or.push_back(q3_float);

		if (q3f.eof()) {
			cout << "### End of q3 file " << i << endl;
		}

		if (q3or[i] < threshold) {		// check if liquid-like
			liquid.push_back(i);		// liquid vector contains id of liquid-like Ox
		} else {
			crystal.push_back(i);
		}
	}


	cout << "#Q3FRAME Number Ice-like Ox: " << crystal.size()
			<< " Number Liquid-like Ox: " << liquid.size() << endl;
	cout << "#Q3FRAME Threshold Value used: " << threshold << endl;

} //end function

// in this function: determine all oxygens that are crystal-like AND within the QLL
void q3_r::qlloxygens(const vector<int> &crystal, const float* coord_norm,
		const float topInterface, const float botInterface) {

	//Determine if molecules are in top or bottom QLL
	for (int i = 0; i < num_ox; i++) {
		if (coord_norm[i * 3] > topInterface) {
			qllTop.push_back(i);
		} else if (coord_norm[i * 3] < botInterface) {
			qllBot.push_back(i);
		}
	}

	// loop through ice-like Ox:
	for (int i = 0; i < crystal.size(); i++) {// crystal[i] gives Ox id of ith "crystal Ox"
		// determine if crystal Ox i is in QLL
		if (coord_norm[crystal[i] * 3] > topInterface) {
			qllCrystalTop.push_back(crystal[i]);
			qllCrystalOx.push_back(crystal[i]); 	// qllCrystal[i] == Ox id	
		} else if (coord_norm[crystal[i] * 3] < botInterface) {
			qllCrystalBot.push_back(crystal[i]);
			qllCrystalOx.push_back(crystal[i]);
		}
	}

	int liq_count = 0;
	int liqTop_count = 0;
	for (int i = 0; i < liquid.size(); i++) {
		if (coord_norm[liquid[i] * 3] > topInterface
				|| coord_norm[liquid[i] * 3] < botInterface) {
			liq_count = liq_count + 1;
		}
		if (coord_norm[liquid[i] * 3] > topInterface) {
			liqTop_count = liqTop_count + 1;
		}
	}

	// loop through liquid-like Ox:
	for (int i = 0; i < liquid.size(); i++) {  		// liquid[i] gives Ox id of ith "liquid" Ox"
		// determine if liquid Ox i is in QLL
		if (coord_norm[liquid[i] * 3] > topInterface) {
			qllLiquidTop.push_back(liquid[i]);
			qllLiquidOx.push_back(liquid[i]);     	// qllLiquidOx[i] == Ox id
		} else if (coord_norm[liquid[i] * 3] < botInterface) {
			qllLiquidBot.push_back(liquid[i]);
			qllLiquidOx.push_back(liquid[i]);
		}
	}

	cout << "#QLLOXYGENS number top layer: " << qllTop.size()
			<< " bottom layer: " << qllBot.size() << endl;

	cout << "#QLLOXYGENS Number Ice-like Oxygens: " << crystal.size()
			<< " Number QLL Ice-like Oxygens:   " << qllCrystalOx.size()
			<< " Number TOP QLL Ice-like Oxygens: " << qllCrystalTop.size()
			<< endl;
	cout << "#QLLOXYGENS Number liquid-like Oxygens:  " << liquid.size()
			<< " Number QLL liquid-like Oxygens: " << liq_count
			<< " Number TOP QLL liquid-like Oxygens: " << liqTop_count << endl;

}

void q3_r::clear_data() {
	liquid.clear();
	crystal.clear();
	LiquidConnectVec.clear();
	qllCrystalOx.clear();
	qllCrystalConnectVec.clear();
	q3or.clear();
	qllCrystalTop.clear();
	qllCrystalBot.clear();

	qllCrystalTop.clear();
	qllCrystalBot.clear();
	qllLiquidTop.clear();
	qllLiquidBot.clear();
	qllLiquidOx.clear();
	qllTop.clear();
	qllBot.clear();

}

q3_r::~q3_r() {
	q3f.close();
}
