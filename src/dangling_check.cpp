#include <cstdlib> 	
#include <iostream> 	
#include <vector>
#include <ctime> 
#include <numeric>			// std::inner_product
#include <math.h>       	// sqrt, pow, round
#include <string.h>			// memset
#include "dcd_r.hpp"
#include "dangling_check.hpp"
#include "dataCell_pbc.hpp"

using namespace std;

// Determine which OH bonds are dangling. Build neighbour list


dangling::dangling(int nox) {
	num_ox = nox;
}

void dangling::danglingOx(const int frame, const float* x, const float* y,
		const float* z, bool IsPrism, const float box_midpoint) {

	pbc_check cellf;

	// vector of vector iniliased with num_ox rows and 0 columns (each with a value 0). 
	// Row num is Ox id. Number columns will correspond to number of neighbours. 
	Odist_Array.resize(num_ox, vector<int>(0, 0));

	vector<int> neighlist(num_ox, 0); // vector with num_ox rows with ever element =0

	vector<int> dang1;
	vector<int> dang2;

	float xx, yy, zz;
	float dist_sqr;
	float cut_sqr = pow(3.5, 2.0);
	float cut_ang = sqrt(3) / 2;

	// loop over oxygens:
	for (int i = 0; i < num_ox; i++) {

		// loop over oxygens ii>i:
		for (int ii = i + 1; ii < num_ox; ii++) {
			// Check if O-O are neighbours. 
			// To speed up code, skip ii if any of x, y or z are separated by more than cut-off
			// Check axis normal to surface first as this speeds up code (more entries to skip)		

			dist_sqr = 0.0;

			if (IsPrism) {
				xx = x[i * 3] - x[ii * 3];
				cellf.pbc_xcoord(xx);
				dist_sqr = cellf.getdistx() + dist_sqr;
				if (dist_sqr >= cut_sqr) {
					continue;
				}

				yy = y[i * 3] - y[ii * 3];
				cellf.pbc_ycoord(yy);
				dist_sqr = dist_sqr + cellf.getdisty();
				if (dist_sqr >= cut_sqr) {
					continue;
				}
			} else {

				yy = y[i * 3] - y[ii * 3];
				cellf.pbc_ycoord(yy);
				dist_sqr = dist_sqr + cellf.getdisty();
				if (dist_sqr >= cut_sqr) {
					continue;
				}

				xx = x[i * 3] - x[ii * 3];
				cellf.pbc_xcoord(xx);
				dist_sqr = cellf.getdistx() + dist_sqr;
				if (dist_sqr >= cut_sqr) {
					continue;
				}
			}

			zz = z[i * 3] - z[ii * 3];
			cellf.pbc_zcoord(zz);
			dist_sqr = dist_sqr + cellf.getdistz();
			if (dist_sqr >= cut_sqr) {
				continue;
			}
			// i and ii are neighbours if reached here
			// Array populates row i (Ox id i) with integer value ii (neighbour id number ii):
			Odist_Array.at(i).push_back(ii);
			Odist_Array.at(ii).push_back(i);

		} // for loop for neighbour list

		// if i has no neigbours then it is dangling:		 	
		if (Odist_Array[i].empty()) {
			cout << "#DANGLING Dangling(No Neighbours), i: " << i
					<< " Number Neighbours: " << Odist_Array[i].size() << endl;
			dang1.push_back(i);
			dang2.push_back(i);
		}
		// if not, it still may be dangling depending on angles:
		else {
			int count_ang[] = { 0, 0 };
			// vector HO--->O(neigh):
			vector<float> OO(3);
			// vector O-->H				
			vector<float> OH(3);
			float mag_OO, mag_OH;
			float xxOO, yyOO, zzOO;
			float xxOH, yyOH, zzOH;
			float dot_OO_OH, cosAng;
			int neigh_tmp;

			// loop over neighbours of Ox i:
			for (int n = 0; n < Odist_Array[i].size(); n++) {

				// neigh_tmp is the id num of neighbour to Ox i
				neigh_tmp = Odist_Array[i][n];
				xxOO = x[neigh_tmp * 3] - x[i * 3];
				yyOO = y[neigh_tmp * 3] - y[i * 3];
				zzOO = z[neigh_tmp * 3] - z[i * 3];

				// pbc; return updated x,y,z separations
				cellf.pbc_xyz(xxOO, yyOO, zzOO);
				xxOO = cellf.getxpbc();
				yyOO = cellf.getypbc();
				zzOO = cellf.getzpbc();

				// vector and mag of HO--->O
				OO = {xxOO, yyOO, zzOO};
				mag_OO = sqrt(
						inner_product(OO.begin(), OO.end(), OO.begin(), 0.0));

				// loop over both hydrogen of Ox i:
				for (int j = 0; j < 2; j++) {

					xxOH = x[(i * 3) + j + 1] - x[i * 3];
					yyOH = y[(i * 3) + j + 1] - y[i * 3];
					zzOH = z[(i * 3) + j + 1] - z[i * 3];

					cellf.pbc_xyz(xxOH, yyOH, zzOH);
					xxOH = cellf.getxpbc();
					yyOH = cellf.getypbc();
					zzOH = cellf.getzpbc();

					// vec O-->H (NB this syntax only works with c++11):
					OH = {xxOH, yyOH, zzOH};
					mag_OH = sqrt(
							inner_product(OH.begin(), OH.end(), OH.begin(),
									0.0));

					// dot product of OH and OO vectors:
					dot_OO_OH = inner_product(OH.begin(), OH.end(), OO.begin(),
							0.0);
					cosAng = dot_OO_OH / (mag_OO * mag_OH);

					// check if angle greater than 30:
					if (cosAng < cut_ang) {
						count_ang[j] = count_ang[j] + 1;
					}

				} 	// end loop over H 

			}       // end loop over neighbours n

			// OH1 i is dangling if all neighbours have angle > 30
			// dangling if num neighbours with angle > 30 = total number neighbours	
			if (count_ang[0] == Odist_Array[i].size()) {
				dang1.push_back(i);
			}

			// check if OH2 is dangling:
			if (count_ang[1] == Odist_Array[i].size()) {
				dang2.push_back(i);
			}

		}	// end loop over neighbours ii of Ox i

	} 		//end for loop over all Oxygens i

	danglinglist.push_back(dang1);
	danglinglist.push_back(dang2);
	dang1.clear();
	dang2.clear();

	vector<int> dangtopOH;
	vector<int> dangbotOH;

	// loop over OH1 and OH2. DanglingTop/Bot are each 2d vectors for each H in a molecule.
	// Must distinguish between both OH bonds for when we calculate orientation 
	for (int H = 0; H < 2; H++) {

		// loop over dangling OH(#H) and determine which QLL belongs to using midpoint of system: 
		for (int i = 0; i < danglinglist[H].size(); i++) {
			if (IsPrism) {
				if (x[danglinglist[H][i] * 3] >= box_midpoint) {
					dangtopOH.push_back(danglinglist[H][i]);
				} else {
					dangbotOH.push_back(danglinglist[H][i]);
				}
			}
			else {
				if (y[danglinglist[H][i] * 3] >= box_midpoint) {
					dangtopOH.push_back(danglinglist[H][i]);
				} else {
					dangbotOH.push_back(danglinglist[H][i]);
				}
			}
		}

		DanglingTop.push_back(dangtopOH);
		DanglingBot.push_back(dangbotOH);
		dangtopOH.clear();
		dangbotOH.clear();

	}
	cout << "Complete Dangling! " << endl;

}	

vector<vector<int> > dangling::getDangVec() const {
	return danglinglist;
}

vector<vector<int> > dangling::getOxArrayVec() const {
	return Odist_Array;
}

void dangling::clear_data() {
	danglinglist.clear();
	Odist_Array.clear();
	DanglingTop.clear();
	DanglingBot.clear();
}

dangling::~dangling() {
}

pbc_check::~pbc_check() {
}

