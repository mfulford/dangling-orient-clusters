#include <cstdlib>      // need for "exit" 
#include <iostream>     // input/output streams library (cout) 
#include <algorithm>    // max vector size
#include <cstring> 		// memset
#include <algorithm>    // std::partial_sort
#include <math.h>       /* fabs */
#include "bin.hpp"

using namespace std;

/* In this class: 
 * 1) System sliced into bins to obtain density/q3/orient profiles along axis normal to surfaces
 * 2) QLLs split into grid (in plane) and average density/q3/orient obtained. This allows heat maps to be plotted
 */

// Slice Up box into bins along axis normal to surface for averaging: 
float bin::Initialise(bool FirstFrame, bool IsPrism, int NOX,
		const float* coord_norm) {
	float ox_coord[NOX];
	int binID, binID_short, binID_q3;

	// Oxygen coordinates: 
	for (int iox = 0; iox < NOX; iox++) {
		ox_coord[iox] = coord_norm[iox * 3];
	}

	// Initialise bin params. if 1st frame
	if (FirstFrame) {
		nbins_short = 50; 	// nbins_short must be a factor of nbins.  
		nbins = 1000;	
		nbins_q3 = 25;
		ratio_bins = (float) (nbins / nbins_short);
		ratio_q3 = (float) (nbins / nbins_q3);
	}

	if (FirstFrame) {
		// determine max and minimum ox_coord coodinates to help position bins
		// Taking into account distance between outermost and 4th to avoid bias from "vapour" molecules
		vector<float> ox_tmp(ox_coord, ox_coord + NOX);

		partial_sort(ox_tmp.begin(), ox_tmp.begin() + 20, ox_tmp.end(),
				greater<float>());
		int tmp_id = 0;
		while (ox_tmp[tmp_id] - ox_tmp[tmp_id + 3] > 2.5) {
			tmp_id = tmp_id + 1;
		}
		maxbin = ox_tmp[tmp_id];

		partial_sort(ox_tmp.begin(), ox_tmp.begin() + 20, ox_tmp.end(),
				less<float>());
		tmp_id = 0;
		while (fabs(ox_tmp[tmp_id] - ox_tmp[tmp_id + 3]) > 2.5) {
			tmp_id = tmp_id + 1;
		}
		minbin = ox_tmp[tmp_id];

		maxbin = maxbin + 1.5;
		minbin = minbin - 1.5;
		range = maxbin - minbin;

		binwidth = range / nbins;
		float binwidth_short = range / nbins_short;
		float binwidth_q3 = range / nbins_q3;
		binDataArray[0] = maxbin;
		binDataArray[1] = minbin;
		binDataArray[2] = binwidth;
		binDataArray[3] = binwidth_short;
		binDataArray[4] = binwidth_q3;

		cout << "#BinInfo. maxbin: " << maxbin << " minbin: " << minbin
				<< " range: " << range << " bindiwth: " << binwidth << endl;
	}

	NumElementsVec.resize(nbins, 0);
	NumElementsVec_short.resize(nbins_short, 0);
	NumElementsVec_q3.resize(nbins_q3, 0);

	// vector with number of rows = number bins. For each row/bin, each column gives the id of Ox atoms that are within the bin: 
	binVec.resize(nbins, vector<int>(0, 0));
	binVec_short.resize(nbins_short, vector<int>(0, 0)); 	
	binVec_q3.resize(nbins_q3, vector<int>(0, 0));

	// determine which bin each oxygen belongs to:
	for (int i = 0; i < NOX; i++) {
		if ((ox_coord[i] < minbin) | (ox_coord[i] > maxbin)) {
			cout << "##BinError: atom " << i
					<< "  outside orignal range. Cycle " << endl;
		} else {
			// floor rounds down. First bin is 0 (at top of box). Last is 23 (at bottom minbin) 
			binID = floor((maxbin - ox_coord[i]) / binwidth);
			if (binID == nbins) { 	// this occurs if coordinate = maxbin
				binID = nbins - 1;
			}
			binVec.at(binID).push_back(i);
			binID_short = floor(binID / ratio_bins);
			binID_q3 = floor(binID / ratio_q3);
			binVec_short.at(binID_short).push_back(i);
			binVec_q3.at(binID_q3).push_back(i);
		}
	}

	// density profile along coord_norm axis
	for (int i = 0; i < nbins; i++) {
		NumElementsVec[i] = binVec[i].size();
	}
	for (int i = 0; i < nbins_short; i++) {
		NumElementsVec_short[i] = binVec_short[i].size();
	}
	for (int i = 0; i < nbins_q3; i++) {
		NumElementsVec_q3[i] = binVec_q3[i].size();
	}
	for (int i = 0; i < nbins; i++) {
		BinAvVec.push_back((float) binVec[i].size());
	}

}

// In the following q3 is summed up so that the bin average can later be calculated:
void bin::Average(const std::vector<float> &q3or,
		const vector<vector<int> > &binVec) {

	int idOx;
	float binAv;

	// loop over number bins: 
	for (int i = 0; i < binVec.size(); i++) {
		binAv = 0.0;
		for (int j = 0; j < binVec[i].size(); j++) {
			idOx = binVec[i][j];
			binAv = q3or[idOx] + binAv;
		}
		BinAvVec.push_back(binAv);
	}
}

// In the following the sum of theta (orientation) in each bin is performed. Later the average will be calculated for each bin:
void bin::Average(const std::vector<std::vector<float> > &orient,
		const vector<vector<int> > &binVec) {

	int idOx;
	float binAv;

	for (int i = 0; i < binVec.size(); i++) {
		binAv = 0.0;
		for (int j = 0; j < binVec[i].size(); j++) {
			idOx = binVec[i][j];
			binAv = orient[0][idOx] + binAv;
		}
		BinAvVec.push_back(binAv);

	}
}

// Average Profile of binAvVec over time if QLL thickness is constant:
void bin::TimeAverage(bool FirstFrame, const vector<int> &NumElementsVec) {

	int nbins = BinAvVec.size();

	if (FirstFrame) {
		TimeBinAvVec = BinAvVec;
		tot_NumElementsVec = NumElementsVec;
	} else {
		for (int i = 0; i < nbins; i++) {
			TimeBinAvVec[i] = TimeBinAvVec[i] + BinAvVec[i];
			tot_NumElementsVec[i] = tot_NumElementsVec[i] + NumElementsVec[i];
		}
	}
}

// Average Profile of binAvVec over time if QLL isn't constant. 
// Added dimensionality corresponding to different nqllTop
void bin::TimeAverage(int nqllTop, bool FirstFrame, bool IsPrism,
		const vector<int> &NumElementsVec) {

	width_qll; 				 	// width of nQLL range  
	if (IsPrism) {
		width_qll = 288;
	} else {
		width_qll = 240;
	}

	int nbins = BinAvVec.size();
	int minQLL = 0;
	int maxQLL = 3000;
	int qllbinID;				// nqllbins bin ID nqll belongs to
	int QLLrange = maxQLL - minQLL;

	// Initialise values in the first frame:
	if (FirstFrame) {
		nqllbins = ceil(QLLrange / width_qll) + 1;					// +1 as count starts from 0	
		TimeAv_nQLL.resize(nqllbins);
		count_nQLL = new int[nqllbins];
		memset(count_nQLL, 0, sizeof(count_nQLL[0]) * nqllbins);	// set every element =0;	
		tot_NumElementsVec_nqll.resize(nqllbins);
	}

	// Determine what nqll bin this frame belongs to. Floor rounds down so eg n=77 is bin 0 and 300 is bin 1	
	qllbinID = floor(nqllTop / width_qll);

	// Check if particular nqllbin has been visited. If not: initilaise
	if (count_nQLL[qllbinID] == 0) {
		TimeAv_nQLL[qllbinID] = BinAvVec;
		count_nQLL[qllbinID] = 1;
		tot_NumElementsVec_nqll[qllbinID] = NumElementsVec;
	} else {
		count_nQLL[qllbinID] = 1 + count_nQLL[qllbinID];
		for (int i = 0; i < nbins; i++) {
			TimeAv_nQLL[qllbinID][i] = TimeAv_nQLL[qllbinID][i] + BinAvVec[i];
			tot_NumElementsVec_nqll[qllbinID][i] =
					tot_NumElementsVec_nqll[qllbinID][i] + NumElementsVec[i];
		}
	}

}

// 1-D orientation distribution for each slice
void bin::orient_slice(int nqllTop, const vector<vector<int> > &slices,
		const vector<vector<float> > &orient, const float* coord_inplane,
		const float* z, bool IsConstQLL, bool FirstFrame) {
	if (FirstFrame) {
		n_slice_bins = 100;
	}
	float or0;
	if (IsConstQLL) {
		if (FirstFrame) {
			OrientSlice.resize(num_slices);
			for (int s = 0; s < num_slices; s++) {
				OrientSlice[s].resize(n_slice_bins, 0.0);
			}
		}
		for (int s = 0; s < num_slices; s++) {
			for (int i = 0; i < slices[s].size(); i++) {
				or0 = floor(
						(((orient[0][slices[s][i]]) * n_slice_bins)
								+ n_slice_bins) / 2); 		// range -1 to 1 --> bin 0 to nbins-1 (0 to 99))
				if (or0 == n_slice_bins) {
					or0 = or0 - 1;
				}
				OrientSlice[s][or0] = OrientSlice[s][or0] + 1;
			}
		}
	} else {
		int qllbinID = floor(nqllTop / width_qll);
		if (FirstFrame) {
			OrientSlice_nqll.resize(nqllbins);
			for (int dim1 = 0; dim1 < nqllbins; dim1++) {
				OrientSlice_nqll[dim1].resize(num_slices);
				for (int s = 0; s < num_slices; s++) {
					OrientSlice_nqll[dim1][s].resize(n_slice_bins, 0.0);
				}
			}
		}

		for (int s = 0; s < num_slices; s++) {
			for (int i = 0; i < slices[s].size(); i++) {
				or0 = floor(
						(((orient[0][slices[s][i]]) * n_slice_bins)
								+ n_slice_bins) / 2); 		// range -1 to 1 --> bin 0 to nbins-1 (0 to 99))
				if (or0 == n_slice_bins) {
					or0 = or0 - 1;
				}
				OrientSlice_nqll[qllbinID][s][or0] =
						OrientSlice_nqll[qllbinID][s][or0] + 1;
			}
		}
	}

}
// Normalise slice 1-d orientation distribution after cycled through all frames
void bin::orient_sliceSUM(bool IsConstQLL) {

	if (IsConstQLL) {
		for (int s = 0; s < num_slices; s++) {
			for (int b = 0; b < n_slice_bins; b++) {
				OrientSlice[s][b] = OrientSlice[s][b] / countSlice[s];
			}
		}
	} else {
		for (int n = 0; n < nqllbins; n++) {
			if (count_nQLL[n] > 0) {
				for (int s = 0; s < num_slices; s++) {
					for (int b = 0; b < n_slice_bins; b++) {
						OrientSlice_nqll[n][s][b] = OrientSlice_nqll[n][s][b]
								/ countSlice_nqll[n][s];
					}
				}
			}
		}

	}
}

// Return q3 and orientation profile for writing output if QLL CONSTANT
const vector<float> bin::getTimeBinAvVec(bool IsDensity) {

	int nbins = TimeBinAvVec.size();

	if (IsDensity) {
		int totalnum_elements = 0;
		for (int i = 0; i < nbins; i++) {
			totalnum_elements = tot_NumElementsVec[i] + totalnum_elements;
		}
		for (int i = 0; i < nbins; i++) {
			TimeBinAvVec[i] = TimeBinAvVec[i] / totalnum_elements;
		}
	} else {
		for (int i = 0; i < nbins; i++) {
			if (tot_NumElementsVec[i] > 0) {
				TimeBinAvVec[i] = TimeBinAvVec[i] / tot_NumElementsVec[i];
			}
		}
	}
	return TimeBinAvVec;
}

// Return q3 and orientation profile for writing output if QLL DYNAMIC
const vector<vector<float> > bin::getTimeAv_nQLLVec(bool IsDensity) {
	int count_tmp = 0;

	// loop over nqllbins:
	for (int n = 0; n < nqllbins; n++) {
		// if nQLL has been visited, loop over "z" nbins
		if (count_nQLL[n] > 0) {
			int nbins = TimeAv_nQLL[n].size();
			if (IsDensity) {
				int totalnum_elements = 0;
				for (int i = 0; i < nbins; i++) {
					totalnum_elements = tot_NumElementsVec_nqll[n][i]
							+ totalnum_elements;
				}
				for (int i = 0; i < nbins; i++) {
					TimeAv_nQLL[n][i] = TimeAv_nQLL[n][i] / totalnum_elements;
				}
			} else {
				for (int i = 0; i < nbins; i++) {
					if (tot_NumElementsVec_nqll[n][i] > 0) {
						TimeAv_nQLL[n][i] = TimeAv_nQLL[n][i]
								/ tot_NumElementsVec_nqll[n][i];
					}
				}
			}
		}
	}

	return TimeAv_nQLL;
}

// ---------------------------------------------------------- In Plane Bin ----------------------------------------------
// ---------------------------------------------------------- HEATMAPS --------------------------------------------------


/* Calculate 2-D in-plane distributions/heatmaps. 
 * CONSTANT QLL
 * Following distributions calculated: 
 * A) in-plane coordinates with colour = q3
 * B) in-plane coordinates with colour = cos(angle)
 * C) in-plane coordinates with colour = density
 * E) cos(angle) v gamma   with colour = frequency
 */ 
void bin::InPlaneDensity(int NOX, const vector<vector<int> > &slices,
		const vector<vector<float> > &orient, const vector<float> &q3or,
		bool IsPrism, const float* coord_inplane, const float* z) {

	int n_zbins = 45;
	int n_planebins = 45;
	int n_orientbins = 100;

	num_slices = slices.size(); 

	// orientPlaneVec : num_slices * 100 (cos(theta)) *100 (rho). Each slice has a heatmap of costheta and rho
	orientPlaneVec.resize(num_slices);
	for (int dim1 = 0; dim1 < num_slices; dim1++) {
		orientPlaneVec[dim1].resize(n_orientbins);
		for (int dim2 = 0; dim2 < n_orientbins; dim2++) {
			orientPlaneVec[dim1][dim2].resize(n_orientbins, 0.0);
		}
	}

	PlaneVec.resize(num_slices); 						// Vector dimension: 2*3*n_zbins*n_planes. (top/bot, dens/q3/orient, grid)
	for (int dim1 = 0; dim1 < num_slices; dim1++) {		
		PlaneVec[dim1].resize(3);
		for (int dim2 = 0; dim2 < 3; dim2++) {
			PlaneVec[dim1][dim2].resize(n_zbins);
			for (int dim3 = 0; dim3 < n_zbins; dim3++) {
				PlaneVec[dim1][dim2][dim3].resize(n_planebins, 0.0);
			}
		}
	}

	float z_ox_tmp[NOX];
	float inplane_ox_tmp[NOX];
	// Oxygen coordinates:
	for (int iox = 0; iox < NOX; iox++) {
		z_ox_tmp[iox] = z[iox * 3];
		inplane_ox_tmp[iox] = coord_inplane[iox * 3];
	}

	// Min and max coordinates of z and plane
	float minZ, maxZ, minPlane, maxPlane;
	minZ = *min_element(z_ox_tmp, z_ox_tmp + (NOX));
	maxZ = *max_element(z_ox_tmp, z_ox_tmp + (NOX));
	minPlane = *min_element(inplane_ox_tmp, inplane_ox_tmp + (NOX));
	maxPlane = *max_element(inplane_ox_tmp, inplane_ox_tmp + (NOX));

	float rangeZ = maxZ - minZ;
	float rangePlane = maxPlane - minPlane;
	float binwidthZ = rangeZ / (n_zbins);
	float binwidthPlane = rangePlane / (n_planebins);

	binPlane_data[0] = binwidthZ;
	binPlane_data[1] = binwidthPlane;

	for (int s = 0; s < num_slices; s++) {
		for (int i = 0; i < slices[s].size(); i++) {
			float theta = orient[0][slices[s][i]];
			float rho = orient[1][slices[s][i]];
			int bin_theta = floor(n_orientbins * (theta + 1) / 2); //(- mintheta == - -1 == +1) 2 is range (-1to1). -1 --> bin 0. 1--> bin 99 
			int bin_rho = floor(n_orientbins * rho); 				// range is 0 to +1. Units pi.  
			if (bin_theta == n_orientbins) {
				bin_theta = bin_theta - 1;
			} // occurs if rho or theta = 1
			if (bin_rho == n_orientbins) {
				bin_rho = bin_rho - 1;
			}
			if (bin_rho < 0) {
				cout << "ERROR_RO_BIN " << bin_rho << " " << rho << endl;
			}
			float tmp_z = z_ox_tmp[slices[s][i]];
			float tmp_plane = inplane_ox_tmp[slices[s][i]];
			int ID_zbin = floor(n_zbins * (tmp_z - minZ) / rangeZ); //zbin 0 is at zmin. 
			int ID_planebin = floor(
					n_planebins * (tmp_plane - minPlane) / rangePlane);
			if (ID_zbin == n_zbins) {
				ID_zbin = n_zbins - 1;
			}
			if (ID_planebin == n_planebins) {
				ID_planebin = n_planebins - 1;
			}

			if ((ID_zbin >= n_zbins) || (ID_planebin >= n_zbins)) {
				cout << "ERRORR " << ID_zbin << " " << ID_planebin << endl;
			}

			orientPlaneVec[s][bin_theta][bin_rho] =
					orientPlaneVec[s][bin_theta][bin_rho] + 1.00;
			PlaneVec[s][0][ID_zbin][ID_planebin] =
					PlaneVec[s][0][ID_zbin][ID_planebin] + 1.00;
			PlaneVec[s][1][ID_zbin][ID_planebin] =
					PlaneVec[s][1][ID_zbin][ID_planebin] + q3or[slices[s][i]]; 	// <q3> slice
			PlaneVec[s][2][ID_zbin][ID_planebin] =
					PlaneVec[s][2][ID_zbin][ID_planebin]
							+ orient[0][slices[s][i]]; 							// cos(theta). 
		}
	}

}


// Build up average heatmaps. Import heatmap of current frame: 
void bin::InPlaneTimeAv(int nqllTop, bool IsConstQLL, bool FirstFrame,
		const vector<vector<int> > &slices) {

	if (IsConstQLL) {
		if (FirstFrame) {
			avOrientPlaneVec = orientPlaneVec;
			avPlaneVec = PlaneVec; // vector dimension: 2*3*n_zbins*n_planes. (top/bot, dens/q3/orient, grid)

			countSlice.resize(num_slices);
			for (int s = 0; s < num_slices; s++) {
				countSlice[s] = slices[s].size();
			}
		} else {
			// determine number of columns and rows in grid (45*45)
			int n_z = avPlaneVec[0][0].size();
			int n_plane = avPlaneVec[0][0][0].size();
			for (int layer = 0; layer < num_slices; layer++) {
				countSlice[layer] = countSlice[layer] + slices[layer].size();

				for (int i = 0; i < 100; i++) { // cos(theta)
					for (int j = 0; j < 100; j++) { //rho
						avOrientPlaneVec[layer][i][j] =
								avOrientPlaneVec[layer][i][j]
										+ orientPlaneVec[layer][i][j];
					}
				}

				for (int param = 0; param < 3; param++) {
					//loop over grid (n_z * n_plane). Data in vector of vector (n_z *n_plane) 
					for (int i = 0; i < n_z; i++) {
						for (int j = 0; j < n_plane; j++) {
							avPlaneVec[layer][param][i][j] =
									avPlaneVec[layer][param][i][j]
											+ PlaneVec[layer][param][i][j]; //Layer 0=top, 1=bot; param 1,2,3 = dens,q3,orient

						}
					}
				}
			}
		}
	} else {
		int n_z = PlaneVec[0][0].size();
		int n_plane = PlaneVec[0][0][0].size();
		int qllbinID = floor(nqllTop / width_qll);
		if (FirstFrame) {
			AvPlane_nQLL.resize(nqllbins);
			OrientAvPlane_nQLL.resize(nqllbins);
			Unvisited_nQLL.resize(nqllbins, true);

			countSlice_nqll.resize(nqllbins);
			for (int dim1 = 0; dim1 < nqllbins; dim1++) {
				countSlice_nqll[dim1].resize(num_slices);
			}
		}
		// if qll size unvisited initiliase
		if (Unvisited_nQLL[qllbinID]) {

			OrientAvPlane_nQLL[qllbinID] = orientPlaneVec;
			AvPlane_nQLL[qllbinID] = PlaneVec;		// AvPlane [nqll][layer][param][n_z][n_plane]
			Unvisited_nQLL[qllbinID] = false;

			for (int s = 0; s < num_slices; s++) {
				countSlice_nqll[qllbinID][s] = slices[s].size();
			}

		} else {
			for (int layer = 0; layer < num_slices; layer++) {	// was 2
				countSlice_nqll[qllbinID][layer] =
						countSlice_nqll[qllbinID][layer] + slices[layer].size();

				for (int i = 0; i < 100; i++) { // cos(theta)
					for (int j = 0; j < 100; j++) { //rho
						OrientAvPlane_nQLL[qllbinID][layer][i][j] =
								OrientAvPlane_nQLL[qllbinID][layer][i][j]
										+ orientPlaneVec[layer][i][j];
					}
				}

				for (int param = 0; param < 3; param++) {
					for (int z_grid = 0; z_grid < n_z; z_grid++) {
						for (int pl_grid = 0; pl_grid < n_plane; pl_grid++) {
							AvPlane_nQLL[qllbinID][layer][param][z_grid][pl_grid] =
									PlaneVec[layer][param][z_grid][pl_grid]
											+ AvPlane_nQLL[qllbinID][layer][param][z_grid][pl_grid];
						}
					}
				}
			} // layers loop               

		}

	} // end if QLLconst

}

// Normalise heatmaps after loop over all frames: 	
void bin::NormaliseHeatmaps() {
// avPLaneVec   (vector 2*3*n_zbins*n_planes. (top/bot, dens/q3/orient, grid)
	int n_z = avPlaneVec[0][0].size();
	int n_plane = avPlaneVec[0][0][0].size();

	for (int layer = 0; layer < num_slices; layer++) { 
		for (int param = 1; param < 3; param++) { 		// 0 = dens, 1,2 q3 or
			for (int i = 0; i < n_z; i++) {
				for (int j = 0; j < n_plane; j++) {
					if (avPlaneVec[layer][0][i][j] > 0) {
						avPlaneVec[layer][param][i][j] =
								avPlaneVec[layer][param][i][j]
										/ avPlaneVec[layer][0][i][j];
					}
					if (param == 2) { 					// can normalise density as don't need it anymore for q3 or orient
						avPlaneVec[layer][0][i][j] = avPlaneVec[layer][0][i][j]
								/ countSlice[layer];
					}
				}
			}
		}
	}

	for (int layer = 0; layer < num_slices; layer++) {
		for (int i = 0; i < 100; i++) { 				//costheta
			for (int j = 0; j < 100; j++) { 			// gamma
				avOrientPlaneVec[layer][i][j] = avOrientPlaneVec[layer][i][j]
						/ countSlice[layer];
			}
		}
	}
}

const vector<vector<vector<float> > > bin::getAvOrientPlaneVec() {
	return avOrientPlaneVec;
}

const vector<vector<vector<vector<float> > > > bin::getAvPlaneVec() {
	return avPlaneVec;
}

void bin::NormaliseHeatmaps_nQLL() {
	for (int n = 0; n < nqllbins; n++) {
		if (count_nQLL[n] > 0) {
			int n_z = AvPlane_nQLL[n][0][0].size();
			int n_plane = AvPlane_nQLL[n][0][0][0].size();
			cout << "Check Vec " << count_nQLL[n] << " " << n_z << " "
					<< n_plane << endl;
			for (int layer = 0; layer < num_slices; layer++) { 
				for (int param = 1; param < 3; param++) {
					for (int i = 0; i < n_z; i++) {
						for (int j = 0; j < n_plane; j++) {
							if (AvPlane_nQLL[n][layer][0][i][j] > 0) {
								AvPlane_nQLL[n][layer][param][i][j] =
										AvPlane_nQLL[n][layer][param][i][j]
												/ AvPlane_nQLL[n][layer][0][i][j];
							}
							if (param == 2) {
								AvPlane_nQLL[n][layer][0][i][j] =
										AvPlane_nQLL[n][layer][0][i][j]
												/ countSlice_nqll[n][layer];
							}

						}
					}
				}

				for (int i = 0; i < 100; i++) { 				//costheta
					for (int j = 0; j < 100; j++) { 			// rho
						OrientAvPlane_nQLL[n][layer][i][j] =
								OrientAvPlane_nQLL[n][layer][i][j]
										/ countSlice_nqll[n][layer];
					}
				}

			} // end layers
		} // end if
	} // end n (qllbins)
}

// Return heatmap constant QLL
const vector<vector<vector<vector<vector<float> > > > > bin::getAvPlane_nQLL() {
	return AvPlane_nQLL;
}
// Return heatmap dynamic QLL
const vector<vector<vector<vector<float> > > > bin::getOrientAvPlane_nQLL() {
	return OrientAvPlane_nQLL;
}

void bin::clear_vectors() {
	NumElementsVec.clear();
	NumElementsVec_short.clear();
	NumElementsVec_q3.clear();
	BinAvVec.clear();
	binVec.clear();
	binVec_short.clear();
	binVec_q3.clear();
	PlaneVec.clear();
	orientPlaneVec.clear();
}

