#include <cstdlib>      // need for "exit" 
#include <iostream>     // input/output streams library (cout) 
#include "output.hpp"

using namespace std;

// In this class the data is written to files with appropriate format for any post-processing and visualisation 

Output::Output(const char* dangfile, const char* orientfile,
		const char* clusterfile, const char* noxygenfile,
		const char* interfacefile, const char* clusVisfile,
		const char* q3Visfile, const char* qllVisfile, const int natoms,
		const int num_frames, bool IsPrism) {
	num_atoms = natoms;
	int width;                               // width of nQLL range  
	if (IsPrism) {
		width = 288;
	} else {
		width = 240;
	}
	width_qll = width;

	dangf.open(dangfile);
	dangf
			<< "#Frame; ##  Num Dang Top; ## Num Dang Bot; ## Orient Dang Top; ## Orient Dang Bot; ## Q3 Dang Top; ## Q3 Dang Bot "
			<< endl;

	orientf.open(orientfile);
	orientf << "Number Atoms: " << num_atoms << " Total number frames: "
			<< num_frames << endl;

	clusterf.open(clusterfile);
	clusterf
			<< "#1: Frame, #2: Water-Cluster Top, #3: Water-Cluster Bot, #4 QLL-Ice-Cluster Top, #5 QLL-Ice-Cluster Bot"
			<< endl;

	noxf.open(noxygenfile);
	noxf
			<< "#Frame; #Total nliquid; #Total ncrystal; #Total QLL Top; #Total QLL Bot; # Orient QLL Top; ## Orient QLL Bot; ## Q3 QLL Top; ## Q3 QLL Bot "
			<< endl;

	interfacef.open(interfacefile);
	interfacef
			<< "#Frame; top QLL/Ice interface position; bot QLL/Ice interface position; top QLL width; bot QLL width "
			<< endl;

	ClusVisf.open(clusVisfile);
	q3Visf.open(q3Visfile);
	qllVisf.open(qllVisfile);

}

void Output::output_noxygens(const int frame, const int ntotLiquid,
		const int ntotIce, const int nqll_top, const int nqll_bot,
		const float OrientTopQLL, const float OrientBotQLL,
		const float Q3TopQLL, const float Q3BotQLL) {
	noxf << frame << " " << ntotLiquid << " " << ntotIce << " " << nqll_top
			<< " " << nqll_bot << " " << OrientTopQLL << " " << OrientBotQLL
			<< " " << Q3TopQLL << " " << Q3BotQLL << endl;
}

void Output::output_dang(const int frame, const int num_dangTop,
		const int num_dangBot, const float OrientTopDang,
		const float OrientBotDang, const float Q3TopDang,
		const float Q3BotDang) {
	dangf << frame << " " << num_dangTop << " " << num_dangBot << " "
			<< OrientTopDang << " " << OrientBotDang << " " << Q3TopDang << " "
			<< Q3BotDang << endl;
}

void Output::output_orient(const int frame,
		const std::vector<std::vector<float> > orientvec) {
	orientf << "#Frame: " << frame << endl;
	for (int j = 0; j < num_atoms / 3; j++) {
		orientf << orientvec[0][j] << "  " << orientvec[1][j] << endl;
	}
}

void Output::interface(int frame, float top, float bot, float top_size,
		float bot_size) {
	interfacef << frame << " " << top << " " << bot << " " << top_size << " "
			<< bot_size << endl;
}

void Output::output_cluster(int frame, int WC_top, int WC_bot, int QLLIC_top,
		int QLLIC_bot) {
	clusterf << frame << " " << WC_top << " " << WC_bot << " " << QLLIC_top
			<< " " << QLLIC_bot << endl;
}

void Output::cluster_visual(vector<vector<int> > bulkIce_clusterId,
		vector<vector<int> > water_clusterId, int idGlobalIC,
		int idGlobalWC_top, int idGlobalWC_bot) {

	vector<int> clus_arr(num_atoms, 0);
	for (int i = 0; i < bulkIce_clusterId[idGlobalIC].size(); i++) {
		clus_arr[bulkIce_clusterId[idGlobalIC][i] * 3] = 1;
	}
	for (int i = 0; i < water_clusterId[idGlobalWC_top].size(); i++) {
		clus_arr[water_clusterId[idGlobalWC_top][i] * 3] = 2;
	}
	for (int i = 0; i < water_clusterId[idGlobalWC_bot].size(); i++) {
		clus_arr[water_clusterId[idGlobalWC_bot][i] * 3] = 3;
	}
	for (int i = 0; i < num_atoms; i++) {
		ClusVisf << clus_arr[i] << " ";
	}
	ClusVisf << endl;
	clus_arr.clear();
}

void Output::q3_visual(vector<int> id_liquid, vector<int> id_ice) {

	vector<int> q3_arr(num_atoms, 0);
	int num_liq = id_liquid.size();
	int num_ice = id_ice.size();
	for (int i = 0; i < num_liq; i++) {
		q3_arr[id_liquid[i] * 3] = 1;
	}
	for (int i = 0; i < num_ice; i++) {
		q3_arr[id_ice[i] * 3] = 2;
	}
	for (int i = 0; i < num_atoms; i++) {
		q3Visf << q3_arr[i] << " ";
	}
	q3Visf << endl;
	q3_arr.clear();
}

void Output::qll_visual(vector<int> id_qllTop, vector<int> id_qllBot) {
	vector<int> qll_arr(num_atoms, 0);
	int num_top = id_qllTop.size();
	int num_bot = id_qllBot.size();
	for (int i = 0; i < num_top; i++) {
		qll_arr[id_qllTop[i] * 3] = 1;
	}
	for (int i = 0; i < num_bot; i++) {
		qll_arr[id_qllBot[i] * 3] = 2;
	}
	for (int i = 0; i < num_atoms; i++) {
		qllVisf << qll_arr[i] << " ";
	}
	qllVisf << endl;
	qll_arr.clear();
}

void Output::InPlane(const char* planefile_angles, const char* planefile_dens,
		const char* planefile_q3, const char* planefile_orient,
		const vector<vector<vector<float> > > &avOrientPlaneVec,
		const vector<vector<vector<vector<float> > > > &avPlaneVec,
		const float binwidthZ, const float binwidthPlane) {
	planef_dens.open(planefile_dens);
	planef_q3.open(planefile_q3);
	planef_orient.open(planefile_orient);
	planef_angles.open(planefile_angles); // theta and gamma

	int n_z = avPlaneVec[0][0].size();
	int n_plane = avPlaneVec[0][0][0].size();
	int n_slices = avPlaneVec.size();
	int n_angles = avOrientPlaneVec[0].size();

	planef_dens << "#Z grid, #plane grid, # Density" << endl;
	planef_q3 << "#Z grid, #plane grid, #q3" << endl;
	planef_orient << "#Z grid, #plane grid, #cos(theta)" << endl;
	planef_angles << "#cos(theta), #rho, #frequency heatmap" << endl;

	planef_dens << "#Number of Slices " << n_slices
			<< " #Range is 0 (min z coord - bottom surface) to " << n_slices - 1
			<< " (max z coord - top surface)" << endl;
	planef_q3 << "#Number of Slices " << n_slices
			<< " #Range is 0 (min z coord - bottom surface) to " << n_slices - 1
			<< " (max z coord - top surface)" << endl;
	planef_orient << "#Number of Slices " << n_slices
			<< " #Range is 0 (min z coord - bottom surface) to " << n_slices - 1
			<< " (max z coord - top surface)" << endl;
	planef_angles << "#Number of Slices " << n_slices
			<< " #Range is 0 (min z coord - bottom surface) to " << n_slices - 1
			<< " (max z coord - top surface)" << endl;

	float zbin, planebin;
	int new_slice;

	for (int i = 0; i < n_z; i++) {
		zbin = i * binwidthZ + binwidthZ / 2.0;
		for (int j = 0; j < n_plane; j++) {
			planebin = j * binwidthPlane + binwidthPlane / 2.0;
			planef_dens << zbin << " " << planebin << " ";
			planef_q3 << zbin << " " << planebin << " ";
			planef_orient << zbin << " " << planebin << " ";

			for (int slice = 0; slice < n_slices; slice++) {
				new_slice = (slice * -1) + (n_slices - 1); // Now slice 0 --> min z coord. Order in file is from min to max z coordinate as column number increases
				planef_dens << avPlaneVec[new_slice][0][i][j] << " ";
				planef_q3 << avPlaneVec[new_slice][1][i][j] << " ";
				planef_orient << avPlaneVec[new_slice][2][i][j] << " ";
			}
			planef_dens << " " << endl;
			planef_q3 << " " << endl;
			planef_orient << " " << endl;
		}
		planef_dens << " " << endl;
		planef_q3 << " " << endl;
		planef_orient << " " << endl;
	}

	float rhobin, thetabin;
	float binwidth_angle = 2.0 / n_angles; // IF range isn't -1 to 1 then need to update this
	float binwidth_gamma = 1.0 / n_angles;
	for (int i = 0; i < n_angles; i++) { // theta
		thetabin = i * binwidth_angle + binwidth_angle / 2.0 - 1.0;
		for (int j = 0; j < n_angles; j++) {
			rhobin = j * binwidth_gamma + binwidth_gamma / 2.0;
			planef_angles << thetabin << " " << rhobin << " ";
			for (int slice = 0; slice < n_slices; slice++) {
				new_slice = (slice * -1) + (n_slices - 1); // Now slice 0 --> min z coord. Order in file is from min to max z coordinate as column number increases
				planef_angles << avOrientPlaneVec[new_slice][i][j] << " ";
			}
			planef_angles << " " << endl;
		}
		planef_angles << " " << endl;
	}

}

// non const QLL
void Output::InPlane(const char* planefile_angles, const char* planefile_dens,
		const char* planefile_q3, const char* planefile_orient,
		const vector<vector<vector<vector<float> > > > &OrientAvPlane_nQLL,
		const vector<vector<vector<vector<vector<float> > > > > &AvPlane_nQLL,
		const int* count_nQLL, const float binwidthZ,
		const float binwidthPlane) {

	//AvPlane_nQLL [nqll][layer][param][n_z][n_plane]

	planef_dens.open(planefile_dens);
	planef_q3.open(planefile_q3);
	planef_orient.open(planefile_orient);
	planef_angles.open(planefile_angles); //theta and gamma

	int nqllbins = AvPlane_nQLL.size();
	int temp_slice;
	// Determine grid size---------------------------------:
	int n_z, n_plane, n_angles;
	for (int i = 0; i < nqllbins; i++) {
		if (count_nQLL[i] > 0) {
			n_z = AvPlane_nQLL[i][0][0].size();
			n_plane = AvPlane_nQLL[i][0][0][0].size();
			n_angles = OrientAvPlane_nQLL[i][0].size();
			temp_slice = AvPlane_nQLL[i].size();
			break;
		}
	}
	// ------------------------------------------------------

	//Determine which QLL thickness are present-------------: 
	vector<int> nqllrow_vec;
	for (int i = 0; i < nqllbins; i++) {
		if (count_nQLL[i] > 0) {
			nqllrow_vec.push_back(i);
		}
	}
	int nQLL = nqllrow_vec.size();
	// --------------------------------------------------------

	float zbin, planebin;
	int n_slices = AvPlane_nQLL[nqllrow_vec[0]].size(); // first dimension is for QLL size
	int new_slice;

	cout << "AAAAAAAOOOOOO " << n_slices << " # " << AvPlane_nQLL[0].size()
			<< " # " << AvPlane_nQLL[1].size() << " # " << nqllrow_vec[0]
			<< " # " << nQLL << " # " << nqllbins << endl;
	cout << "BBBBBBBBBBBBB " << temp_slice << endl;

	planef_dens << "#Z grid, #plane grid, # Density" << endl;
	planef_q3 << "#Z grid, #plane grid, #q3" << endl;
	planef_orient << "#Z grid, #plane grid, #orientation" << endl;
	planef_angles << "#cos(theta), #rhp, #frequency" << endl;

	planef_dens << "#Number of Slices " << n_slices
			<< " #Range is 0 (min z coord - bottom surface) to " << n_slices - 1
			<< " (max z coord - top surface)" << endl;
	planef_q3 << "#Number of Slices " << n_slices
			<< " #Range is 0 (min z coord - bottom surface) to " << n_slices - 1
			<< " (max z coord - top surface)" << endl;
	planef_orient << "#Number of Slices " << n_slices
			<< " #Range is 0 (min z coord - bottom surface) to " << n_slices - 1
			<< " (max z coord - top surface)" << endl;
	planef_angles << "#Number of Slices " << n_slices
			<< " #Range is 0 (min z coord - bottom surface) to " << n_slices - 1
			<< " (max z coord - top surface)" << endl;

	planef_dens
			<< "#Columns ordered by number of QLL molecules in top surface. Mid point of QLL bins sampled listed below (+- "
			<< width_qll / 2 << ") molecules in bin: " << endl;
	planef_q3
			<< "#Columns ordered by number of QLL molecules in top surface. Mid point of QLL bins sampled listed below (+- "
			<< width_qll / 2 << ") molecules in bin: " << endl;
	planef_orient
			<< "#Columns ordered by number of QLL molecules in top surface. Mid point of QLL bins sampled listed below (+- "
			<< width_qll / 2 << ") molecules in bin: " << endl;
	planef_angles
			<< "#Columns ordered by number of QLL molecules in top surface. Mid point of QLL bins sampled listed below (+- "
			<< width_qll / 2 << ") molecules in bin: " << endl;

	planef_dens << "#QLL size: ";
	planef_q3 << "#QLL size: ";
	planef_orient << "#QLL size: ";
	planef_angles << "#QLL size: ";
	for (int i = 0; i < nQLL; i++) {
		planef_dens << "#" << (nqllrow_vec[i] * width_qll + width_qll / 2.0)
				<< "(slice columns {min z to max z}: " << (i * n_slices) + 3
				<< " to " << (i * n_slices) + 2 + n_slices << ")       ";
		planef_q3 << "#" << (nqllrow_vec[i] * width_qll + width_qll / 2.0)
				<< "(slice columns {min z to max z}: " << (i * n_slices) + 3
				<< " to " << (i * n_slices) + 2 + n_slices << ")       ";
		planef_orient << "#" << (nqllrow_vec[i] * width_qll + width_qll / 2.0)
				<< "(slice columns {min z to max z}: " << (i * n_slices) + 3
				<< " to " << (i * n_slices) + 2 + n_slices << ")       ";
		planef_angles << "#" << (nqllrow_vec[i] * width_qll + width_qll / 2.0)
				<< "(slice columns {min z to max z}: " << (i * n_slices) + 3
				<< " to " << (i * n_slices) + 2 + n_slices << ")       ";
	}
	planef_dens << endl;
	planef_q3 << endl;
	planef_orient << endl;
	planef_angles << endl;

	for (int i = 0; i < n_z; i++) {
		zbin = i * binwidthZ + binwidthZ / 2.0;
		for (int j = 0; j < n_plane; j++) {
			planebin = j * binwidthPlane + binwidthPlane / 2.0;
			planef_dens << zbin << " " << planebin << " ";
			planef_q3 << zbin << " " << planebin << " ";
			planef_orient << zbin << " " << planebin << " ";
			for (int k = 0; k < nQLL; k++) {
				int n = nqllrow_vec[k];
				for (int slice = 0; slice < n_slices; slice++) {
					new_slice = (slice * -1) + (n_slices - 1); // Now slice 0 --> min z coord. Order in file is from min to max z coordinate as column number increases
					planef_dens << AvPlane_nQLL[n][new_slice][0][i][j] << " ";
					planef_q3 << AvPlane_nQLL[n][new_slice][1][i][j] << " ";
					planef_orient << AvPlane_nQLL[n][new_slice][2][i][j] << " ";
				}
			}

			planef_dens << " " << endl;
			planef_q3 << " " << endl;
			planef_orient << " " << endl;

		}
		// required for format for plotting in gnuplot (need a line between "y" values)
		planef_dens << " " << endl;
		planef_q3 << " " << endl;
		planef_orient << " " << endl;
	}

	float rhobin, thetabin;
	float binwidth_angle = 2.0 / n_angles; // IF range isn't -1 to 1 then need to update this
	for (int i = 0; i < n_angles; i++) { // theta
		thetabin = i * binwidth_angle + binwidth_angle / 2.0;
		for (int j = 0; j < n_angles; j++) {
			rhobin = j * binwidth_angle + binwidth_angle / 2.0;
			planef_angles << thetabin << " " << rhobin << " ";
			for (int k = 0; k < nQLL; k++) {
				int n = nqllrow_vec[k];
				for (int slice = 0; slice < n_slices; slice++) {
					new_slice = (slice * -1) + (n_slices - 1); // Now slice 0 --> min z coord. Order in file is from min to max z coordinate as column number increases
					planef_angles << OrientAvPlane_nQLL[n][new_slice][i][j]
							<< " ";
				}
			}
			planef_angles << " " << endl;
		}
		planef_angles << " " << endl;
	}

}

// Av orientation profile for constant QLL
void Output::AverageTimeProfile(const char* binfile,
		const vector<float> TimeBinAvVec, const float binwidth) {

	binAvf.open(binfile);

	int nbins = TimeBinAvVec.size();

	float bin_pos;
	int new_bin;

	binAvf
			<< "#1st Column is 'z'-coordinate. 2nd is profile averaged over entire simulation"
			<< endl;
	binAvf << "#max z-coordinate is top of system" << endl;

	for (int i = 0; i < nbins; i++) {
		// invert bin pos. 
		// Before: bin=0 --> max "z". 
		// Now bin=0 --> new_bin=nbins-1 --> max z
		new_bin = (-1 * i) + (nbins - 1);
		bin_pos = new_bin * binwidth + binwidth / 2.0;
		binAvf << bin_pos << " " << TimeBinAvVec[i] << endl;
	}

	binAvf.close();

}

// Av orientation profile for dynamic QLL size
void Output::AverageTimeProfile(const char* binfile,
		const vector<vector<float> > TimeAv_nQLL, const int* count_nQLL,
		const float binwidth) {
	/* TimeAv_nQLL has n rows (where n=range of QLL size/binwdith) 
	 * For each QLL size sampled the correct row is populated with the average profile. 
	 * These rows will have nbins columns depending on how many bins the system was divided into
	 * Rows which have 0 columns correspond to QLL sizes which were not sampled during the simulation
	 */

	binAvf.open(binfile);

	vector<int> nqllrow_vec;
	int n;

	// Number of bins system is split into:
	int nqllbins = TimeAv_nQLL.size();
	int nbins;
	for (int i = 0; i < nqllbins; i++) {
		if (count_nQLL[i] > 0) {
			nbins = TimeAv_nQLL[i].size();
			break;
		}
	}

	for (int i = 0; i < nqllbins; i++) {
		if (count_nQLL[i] > 0) {
			nqllrow_vec.push_back(i);
		}
	}

	// Number of QLL sizes were sampled:
	int nQLL = nqllrow_vec.size();

	binAvf
			<< "#1st Column is 'z'-coordinate. Remaining Columns are profile for the different number of QLL molecules  sampled. Mid point of QLL bin listed below."
			<< endl;
	binAvf << "#max z-coordinate is top of system " << endl;
	binAvf << "#Number QLL Molecules (+-" << width_qll / 2.0
			<< ") molecules in bin: ";
	for (int i = 0; i < nQLL; i++) {
		binAvf << "#" << (nqllrow_vec[i] * width_qll + width_qll / 2.0)
				<< " (column number: " << i + 1 << ")     ";
	}
	binAvf << endl;

	float bin_pos;
	int new_bin;
	// loop over number rows (bins)
	// For each bin, loop over number columns (that are populated)	
	for (int bin = 0; bin < nbins; bin++) {
		// invert bin pos. 
		// Before: bin=0 --> max "z". 
		// Now bin=0 --> new_bin=nbins-1 --> max z
		new_bin = (-1 * bin) + (nbins - 1);
		bin_pos = new_bin * binwidth + binwidth / 2.0;
		binAvf << bin_pos << " ";
		for (int col = 0; col < nQLL; col++) {
			n = nqllrow_vec[col];
			binAvf << TimeAv_nQLL[n][bin] << " ";
		}
		binAvf << endl;
	}

	nqllrow_vec.clear();

	binAvf.close();
}

void Output::orientSlice(const char* orientSlicefile,
		const vector<vector<float> > &OrientSlice) {

	orientSlicef.open(orientSlicefile); // file: col1=bin, col2--n prob dist orient for slice 0 to ..
	int n_orient_bins = OrientSlice[0].size();
	int n_slices = OrientSlice.size();
	int new_slice;
	float binwidth_orient = 2.00 / n_orient_bins;
	float b_pos;

	cout << "OrientSlice Check: " << binwidth_orient << " " << n_orient_bins
			<< " " << endl;

	orientSlicef
			<< "#orientation probability distribution in each slice. Averaged over entire simulation (with slice position constant)"
			<< endl;
	orientSlicef
			<< "#Col1 is orientation bin mid-point value. #Col 2 (min z coord slice) to "
			<< n_slices + 1
			<< " are probability distribution for different slices" << endl;
	orientSlicef << "#Number of Slices " << n_slices
			<< " #Range is 0 (min z coord - bottom surface) to " << n_slices - 1
			<< " (max z coord - top surface)" << endl;

	for (int b = 0; b < n_orient_bins; b++) {
		b_pos = (((float(b) / n_orient_bins) * 2) - 1)
				+ (binwidth_orient / 2.0); // range from 0 to nbins-1 to -1 to +1 (with value at midpoint) 
		orientSlicef << b_pos << " ";
		for (int s = 0; s < n_slices; s++) {
			new_slice = (s * -1) + (n_slices - 1); // Now slice 0 --> min z coord. Order in file is from min to max z coordinate as column number increases
			orientSlicef << OrientSlice[new_slice][b] << " ";
		}
		orientSlicef << endl;
	}
	cout << endl;
}

void Output::orientSlice(const char* orientSlicefile,
		const vector<vector<vector<float> > > &OrientSlice_nqll,
		const int* count_nQLL) {
	orientSlicef.open(orientSlicefile);

	int nqllbins = OrientSlice_nqll.size();
	int n_slices, n_orient_bins;
	for (int i = 0; i < nqllbins; i++) {
		if (count_nQLL[i] > 0) {
			n_slices = OrientSlice_nqll[i].size();
			n_orient_bins = OrientSlice_nqll[i][0].size();
			break;
		}
	}

	vector<int> nqllrow_vec;
	for (int i = 0; i < nqllbins; i++) {
		if (count_nQLL[i] > 0) {
			nqllrow_vec.push_back(i);
		}
	}
	// Number of QLL sizes were sampled:
	int nQLL = nqllrow_vec.size();
	float b_pos;
	float binwidth_orient = 2.00 / n_orient_bins;

	orientSlicef
			<< "#orientation probability distribution in each slice. Non Const QLL size: slices split according to QLL size of top surface (with slice positions remaining constant)"
			<< endl;
	orientSlicef
			<< "#Col1 is orientation bin mid-point value. Remaining columns are orientation prob. distribution for all QLL size and slices"
			<< endl;
	orientSlicef << "#Number of Slices " << n_slices
			<< " #Range is 0 (min z coord - bottom surface) to " << n_slices - 1
			<< " (max z coord - top surface)" << endl;
	orientSlicef
			<< "#Columns ordered by number of QLL molecules in top surface. Number of sampled QLL sizes: "
			<< nQLL << endl;
	orientSlicef << "#Number QLL Molecules (+-" << width_qll / 2
			<< ") molecules in bin: ";
	for (int i = 0; i < nQLL; i++) {
		orientSlicef << "#" << (nqllrow_vec[i] * width_qll + width_qll / 2.0)
				<< "(slice columns {min z to max z}: " << (i * n_slices) + 3
				<< " to " << (i * n_slices) + 2 + n_slices << ")       ";
	}

	for (int b = 0; b < n_orient_bins; b++) {
		b_pos = (((float(b) / n_orient_bins) * 2) - 1)
				+ (binwidth_orient / 2.0); //range from 0 to nbins-1 to -1 to +1 (with value at midpoint) 
		orientSlicef << b_pos << " ";
		for (int q = 0; q < nQLL; q++) {
			int n = nqllrow_vec[q];
			for (int s = 0; s < n_slices; s++) {
				orientSlicef << OrientSlice_nqll[n][s][b] << " ";
			}
		}
		orientSlicef << endl;
	}

}

Output::~Output() {
	orientf.close();
	dangf.close();
}
