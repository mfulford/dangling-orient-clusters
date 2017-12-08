#include <iostream>
#include <vector>
#include <numeric>
#include <math.h>
#include <cmath>
#include "orientation.hpp"
#include "dataCell_pbc.hpp"

using namespace std;

/* In the following: determine the orientation of molecules and calulcate Euler angles. See wikipedia for more info
 * For EACH FRAME:
 * For each molecule calculate cos(theta) between dipole and z axis unit vector
 * z-axis unit vector points to surface. Therefore 2 options:
 * 		A) calculate cos(theta) for all molecules with same z unit vector (so that on each surface will get opposite result
 * 		as it only points to one of the 2 surfaces)
 * 		B) in the middle of the box/bulk swap direction of unit vector so that it always points to surface.
 *
 * GOING FOR OPTION A).
 * Pass dangling statistics into this class for stats.
 *
 * Normal axes: z points up, y points right, x points towards you.
 * dipole = (xdipole,ydipole,zdipole). mag = magnitude. Dipole vector determined using midpoint of H (pbc applied) and Oxygen atom.
 * pbc applied for final vector
 * xd = xdipole/mag,  yd = ydipole/mag etc
 * cos(theta) = zd (component of dipole along z - normalised)
 * Second angle (gamma) calculated from dot product of line of node N and H-H vector
 *
 * Details of basal and prism molecular frames
 * Basal Surface: rotate about x 90 degrees( pi/2)
 * Basal Surface:  x--> x | y--> -z | z--> y    (alternative x --> z | y --> x | z--> y. Shift of pi/2 relative to other basis)
 * Top/bottom half: rotate about x (after initial rotation): x-->x | y-->z | z --> -y	(alt x--> -z | y --> x | z --> -y)
 *
 * Prism Surface: rotate about y 90 degrees
 * Prism Surface: x--> -z | y--> y | z--> x    (alternative x--> y | y--> z | z--> x. Shift pi/2)
 * Top/bottom half: rotate about y (after initial): x--> z | y --> y | z--> -x. 	(alt: x-->y | y-->-z | z-->-x)
 *
 * Orientation of molecule obtained as follows:
 * 1) Rotation about "z" axis by gamma. Gamma=0 H-H in line with N or x axis (angle of N and H-H = gamma)
 * 2) Followed by rotation of theta about "x" axis. (Angle of dipole and z= theta) (in wiki, theta is beta)
 *
 */

orientation::orientation(){
	cout << "check orientation of OH bonds" << endl;
}

// CALCULATE ORIENTATION OF MOLECULES. Main Function in Class.
void orientation::orientationOH(const int num_ox, const int frame, const float* x, const float* y, const float* z, vector<bool> IcePhase, const float box_midpoint) {

	// Intialise for periodic boundary conditions
	pbc_check cellf;

	vector< vector<float> > orient_tmp(2, vector <float>(num_ox,0.0));

	vector<float> OH_vect(3);      
	vector<float> H_vect(3);        
	vector<float> N_vect(3);

	float xx, yy, zz, xx_H, yy_H, zz_H, Nx, Ny, Nz;
	float mag_OH, mag_H, mag_N;


	//loop over all oxygens:	
	for (int i=0;i<num_ox;i++){

		// Mid point between hydrogens
		vector<float> h_mid;		 
		h_mid = cellf.h_midpoint(x[i*3+1], y[i*3+1], z[i*3+1], x[i*3+2], y[i*3+2], z[i*3+2]); // mid point between hydrogens following pbc

		// H_mid --> O vector (dipole)	
		xx = x[i*3] - h_mid[0]; 
		yy = y[i*3] - h_mid[1];
		zz = z[i*3] - h_mid[2];
		cellf.pbc_xyz(xx,yy,zz);
		xx=cellf.getxpbc(); yy=cellf.getypbc(); zz=cellf.getzpbc();
		OH_vect = { xx, yy, zz };
		mag_OH = sqrt( inner_product( OH_vect.begin(), OH_vect.end(), OH_vect.begin(),0.0));
		xx = xx/mag_OH; yy = yy/mag_OH; zz = zz/mag_OH;

		// H-H vector:
		xx_H = x[i*3+2] - x[i*3+1]; 	// doesn't matter if H1 or H2. Will use magnitude
		yy_H = y[i*3+2] - y[i*3+1];	
		zz_H = z[i*3+2] - z[i*3+1];	
		cellf.pbc_xyz(xx_H,yy_H,zz_H);
		xx_H=cellf.getxpbc(); yy_H=cellf.getypbc(); zz_H=cellf.getzpbc();
		H_vect = {xx_H, yy_H, zz_H};
		mag_H  = sqrt( inner_product( H_vect.begin(), H_vect.end(), H_vect.begin(),0.0));
		xx_H = xx_H/mag_H; yy_H = yy_H/mag_H; zz_H = zz_H/mag_H;	


		// Determine if Prism of Basal and top or bottom half:
		if (IcePhase[0] == true){ // X-surface: Secondary Prism
			orient_tmp[0][i] = xx;
			Nx = 0.0;
			Ny=zz;
			Nz=-yy;
			if (x[i*3] < box_midpoint){ 
				orient_tmp[0][i] = -xx;
			}
		}
		else if (IcePhase[1] == true) { // Y-surface: Prism
			orient_tmp[0][i] = yy;
			Nx = -zz;
			Ny = 0.0;
			Nz = xx;
			if (y[i*3] < box_midpoint) {
				orient_tmp[0][i] = -yy;
			}
		}
		else if (IcePhase[2] == true) { // Z-surface: Basal
			orient_tmp[0][i] = zz;
			Nx = yy;
			Ny = -xx;
			Nz = 0.0;
			if (y[i*3] < box_midpoint) {
				orient_tmp[0][i] = -zz;
			}
		}

		//N = cross product of dipole and "z" N = dipole x z
		//Basal "z"=y=(0,1,0). Prism "z"=x=(1,0,0)
		//Gamma is angle of H-H with N
		N_vect = {Nx, Ny, Nz}; // perpendicular to dipole and "z" axis (y for basal)			
		mag_N = sqrt( inner_product( N_vect.begin(), N_vect.end(), N_vect.begin(),0.0));
		Nx = Nx/mag_N; Ny = Ny/mag_N; Nz = Nz/mag_N;
		float gamma = (Nx*xx_H) + (Ny*yy_H) + (Nz*zz_H); // H-H

		if ((gamma > 1.00) | (gamma < -1.00)) {
			if (gamma > 1) {
				gamma = 1;
			}
			else if (gamma < -1) {
				gamma = -1;
			}
		}
		orient_tmp[1][i] = acos(gamma)/3.14159265;

		// Legacy - if wanted alpha would do this.
		//float alpha = Nx * alpha_x[0] + Ny * alpha_x[1] + Nz * alpha_x[2];
		//orient_tmp[1][i] = acos(alpha) / 3.1419265;


		OH_vect.clear();
		H_vect.clear();
		N_vect.clear();

	}
	// COMPLETED Loop Over Oxygens

	// Copy 2d orient_tmp vector into orient and clear vector
	copy(orient_tmp.begin(), orient_tmp.end(), back_inserter(orient)) ;	 

	orient_tmp.clear();

	cout << "Completed H2O orientation! " << endl;

}

void orientation::clear_data(){
	orient.clear();
}

orientation::~orientation() {
}

// Calculate average theta (orientation) of dangling molecules:
void orientation::average(const vector< vector<int> > &DanglingTop, const vector< vector<int> > &DanglingBot) {

	avTop=0.0;
	avBot=0.0;

	// Determine number dangling OH on each surface:
	int numTop = DanglingTop[0].size() + DanglingTop[1].size();
	int numBot = DanglingBot[0].size() + DanglingBot[1].size();

	// Determine the average orientation of dangling OH on each surface:
	for (int H=0;H<2;H++) {

		for (int i=0;i<DanglingTop[H].size();i++) {
			avTop=orient[0][DanglingTop[H][i]] + avTop;
		}
		for (int i=0;i<DanglingBot[H].size();i++) {
			avBot=orient[0][DanglingBot[H][i]] + avBot;
		}
	}

	// average orientation of dangling molecules:
	avTop=avTop/numTop;
	avBot=avBot/numBot;

	cout << "#11ORIENTDANG Top Layer Average: " << avTop << " Bottom Layer: " << avBot << endl;

}


// Calculate average dipole theta of QLL molecules:
void orientation::average(const vector<int> &QLLTop, const vector<int> &QLLBot) {

	avTop=0.0;
	avBot=0.0;

	// Determine number dangling OH on each surface:
	int numTop = QLLTop.size();
	int numBot = QLLBot.size();

	// Determine the average theta of dangling OH on each surface:

	for (int i=0;i<QLLTop.size();i++) {
		avTop=orient[0][QLLTop[i]] + avTop;
	}
	for (int i=0;i<QLLBot.size();i++) {
		avBot=orient[0][QLLBot[i]] + avBot;
	}

	avTop=avTop/numTop;
	avBot=avBot/numBot;

	cout << "#ORIENTQLL Top Layer Average: " << avTop << " Bottom Layer: " << avBot << endl;

}

