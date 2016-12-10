#include <cstdlib> 		 
#include <iostream> 	 
#include <ctime> 		
#include <vector>
#include <algorithm>    
#include "dcd_r.hpp"
#include "dangling_check.hpp" 
#include "dataCell_pbc.hpp" 
#include "orientation.hpp" 
#include "q3_r.hpp"
#include "cluster.hpp"
#include "bin.hpp"
#include "output.hpp"

// g++ -std=c++11 *.cpp -o Executable

/* Description: Code to analyze lammps binary trajectories. Read in DCD file, cell data, q3 bond order data
 * Code will:
 * 	 A) build a neighbour list of oxygen atoms and determine the dangling molecules. 
 * 	 B) calculate orientation of molecules from the Euler angles. Ref. frame is xyz. 
 * 	 C) determine cluster sizes of water-like and ice-like molecules used the Depth-first search algorithm. 
 * 	 D) determine the QLL/Ice interface based distribution of coordinates of largest ice-like and water-like clusters. 
 *   E) calculate heatmaps of q3, orientation, density and euler angles. Heatmaps calculated for slices along z-axis
 *   F) produce density, q3, orientation profiles. 
 *   G) calculate collective statistics such as q3/orientation of dangling mols, ice qll molecules
 *   H) determine number of ice-like and liquid-like molecules 
 *   
 *   q3 refers to the 3rd order Steinhard parameter and must be calculated beforehand and inputed
 */

using namespace std;

int main(int argc, char* argv[]) {
	if (argc < 4) {
		cout << "USAGE: executable.x [dcd file] [cell data] [q3file] " << endl;
		exit(1);
	}

	// Find out about simulation and conditions. Answers can be inputed as text file from command line. 
	char answer_qll;
	bool IsConstQLL;
	cout << "Type Y if this simulation is constant QLL thickness (ie normal MD)"
			<< endl;
	cout
			<< "Otherwise type N if this simulation is metadynamics or MD with dynamic QLL Size"
			<< endl;
	cout
			<< "NB: If option \"N\", code will only work if 1 surface QLL size is dynamic (not both)"
			<< endl;
	cin >> answer_qll;
	if ((answer_qll == 'Y') || (answer_qll == 'y')) {
		IsConstQLL = true;
	} else if ((answer_qll == 'N') || (answer_qll == 'n')) {
		IsConstQLL = false;
	} else {
		cout << "ERROR with input. Please enter N or Y " << endl;
		exit(1);
	}

	char answer_q3;
	float q3threshold;
	float q3or = 0.670;
	float q3m = 0.280;
	float q3tot = 0.665;

	cout << "What q3 are you using?" << endl;
	cout << "Type A if q3 original. Threshold = " << q3or << endl;
	cout << "Type B if <q3m>. Threshold = " << q3m << endl;
	cout << "Type C if q3tot. Threshold = " << q3tot << endl;
	cout << "Type D if you would like to enter a different threshold" << endl;
	cin >> answer_q3;
	if ((answer_q3 == 'A') || (answer_q3 == 'a')) {
		q3threshold = q3or;
	} else if ((answer_q3 == 'B') || (answer_q3 == 'b')) {
		q3threshold = q3m;
	} else if ((answer_q3 == 'C') || (answer_q3 == 'c')) {
		q3threshold = q3tot;
	} else if ((answer_q3 == 'D') || (answer_q3 == 'D')) {
		cout << "Enter your chosen threshold: " << endl;
		cin >> q3threshold;
	} else {
		cout << "Error with input. Please enter A, B, C or D " << endl;
		exit(1);
	}
	cout << "q3 threhsold used: " << q3threshold << endl;

	/* instance of a new object DCD_R attached to a dcd file and cell data file: */
	DCD_R dcdf(argv[1]);
	dataCell cellf(argv[2]);

	/* read the header and print it: */
	dcdf.read_header();
	dcdf.printHeader();

	/* determine Ice phase  [if prism true (1); if basal false (0)]: */
	pbc_check phase_check;
	bool IfPrismPhase = phase_check.IsPhasePrism();

	/* Istance of q3 file and open: */
	q3_r q3f(argv[3], dcdf.getNATOM() / 3, q3threshold);
	/* instance of data output */
	Output outf("dang_time.txt", "orient_frame.txt", "cluster_time.txt",
			"noxygens_time.txt", "interface_time.txt", "cluster_Bulk_tcl.txt",
			"q3_tcl.txt", "qll_tcl.txt", dcdf.getNATOM(), dcdf.getNFILE(),
			IfPrismPhase);
	/* Instance of dangling and orient. Use to calculate stats every frame */
	dangling dangframe(dcdf.getNATOM() / 3);
	orientation orient_frame;

	const float *x, *y, *z;
	float box_midpoint;
	cout << " " << endl;
	bool IsFirstFrame = true;

	//average Q3 and orientation profile
	bin Q3binAv;
	bin OrientbinAv;
	bin Density;
	// Instance of class which will calculate average q3 and cos(theta) of groups of molecules
	Stats Averages_q3Orient;

	//In plane average bin width:
	float binwidthPlane = 0.0;
	float binwidthZ = 0.0;
	float count_planebin = 0.0;

	//--------------------------------------------------------------------------------
	//-------------------- BEGIN MAIN LOOP ------------------------------------------
	// in this loop the coordinates are read in frame by frame 
	for (int i = 0; i < dcdf.getNFILE(); i++) {
		cout << "############################ FRAME " << i << endl;

		//timing
		clock_t begin = clock();

		// Read In One Frame and save xyz coordinates 
		dcdf.read_oneFrame();
		x = dcdf.getX();
		y = dcdf.getY();
		z = dcdf.getZ();

		/* Set coord_norm=coordinate normal to surface x or y */
		static const float* coord_norm;
		static const float* coord_inplane;
		if (IfPrismPhase) {
			coord_norm = x;
			coord_inplane = y;
		} else {
			coord_norm = y;
			coord_inplane = x;
		}

		/* Work out centre of box along axis perpendicular to surface. Do it only for frame 0 
		 * and assume it doesn't drift during the simulation. This avoids worrying about molecules 
		 * which may "evaporate" during metaD */
		if (i == 0) {
			box_midpoint = phase_check.box_centre(dcdf.getNATOM(), coord_norm);
		}

		clock_t t_readData = clock();

		// --------------------------------------------------------------------------------------------------------	
		// ----------------- CALCULATE DANGLING and ORIENTATION stats: --------------------------------------------
		dangframe.danglingOx(i, x, y, z, IfPrismPhase, box_midpoint);
		int numDangTop = dangframe.getDangTopVec()[0].size()
				+ dangframe.getDangTopVec()[1].size();
		int numDangBot = dangframe.getDangBotVec()[0].size()
				+ dangframe.getDangBotVec()[1].size();
		orient_frame.orientationOH(dcdf.getNATOM() / 3, i, x, y, z,
				IfPrismPhase, box_midpoint);
		// ----------------- Calculated dangling + orientation -----------------------------------------------------
		
		clock_t t1_statsDangOrient = clock();

		// IMPORT q3 data for current frame and determien liquid and ice molecules:
		q3f.q3_frame();
		
		// ----------------------------------------------------------------------------------------------------------
		// ----------------- CLUSTERS ANALYSIS: ---------------------------------------------------------------------
		// 1) Water-like cluster size
		// 2) (Bulk)Ice-like cluster size
		// 3) QLL/Ice interface from coordinates of 1) and 2)
		// 4) QLL Ice-like cluster size from QLL/ice Interface (3)
		// Input Array where row is ox id and cols gives id of neighbours. Determine connectivity matrix:

		// 1) Water-like cluster ------------------------------------------------------------------------------------
		// Input id of liquid-like molecules (q3f.liquid) and neighbours of molecules
		cluster WaterCluster(dangframe.getOxArrayVec(), q3f.liquid,
				box_midpoint);
		// connectivity matrix of liquid Ox
		WaterCluster.connectivity();    
		WaterCluster.size(coord_norm, IfPrismPhase);
		// ----------------------------------------------------------------------------------------------------------     
		clock_t t2_WaterClust = clock();

		
		/* In the following we will determine the largest Ice-like cluster. This will correspond to the bulk ice. 
		 * This allows us to assume that all remaining ice-like molecules belong to the QLL. We can then calculate 
		 * the QLL-ice-like cluster size without explicitly defining the QLL/Ice interface. 
		 * UPDATE: BulkIce cluster "infiltrates" QLL (particularly the thin layer). QLL ice-cluster binded to bulk ice
		 * THEREFORE: need to define the interface to determine the QLL ice-cluster size.
		 * 2nd and 3rd largest ice-like cluster will belong to the QLL. (Determine largest on each surface by double checking min/max "z" coordinate)
		 * QLL/Ice interface can be estimate by taking 0.5*(max coord (bulk)ice-cluster + max coord water-cluster). Average over many atoms
		 * This is a more accurate estimate than slicing up the system
		 */

		// 2) (Bulk)Ice-like cluster --------------------------------------------------------------------------------
		//	Input id of ice-like mols (q3f.crystal)
		cluster BulkIceCluster(dangframe.getOxArrayVec(), q3f.crystal,
				box_midpoint);
		BulkIceCluster.connectivity();
		BulkIceCluster.size(coord_norm, IfPrismPhase);
		// ---------------------------------------------------------------------------------------------------------- 
		clock_t t3_IceClust = clock(); 

		// 1) QLL/Ice interface -------------------------------------------------------------------------------------
		// Determine the QLL interface based on the cluster coordinates. FUTURE: Create Class
		int idMaxWC_top = WaterCluster.TopCluster_info[0][0];
		int idGlobalWC_top = WaterCluster.TopCluster_info[0][2];

		int idMaxWC_bot = WaterCluster.BotCluster_info[0][0];
		int idGlobalWC_bot = WaterCluster.BotCluster_info[0][2];

		int idMaxIC = BulkIceCluster.TopCluster_info[0][0];
		int idGlobalIC = BulkIceCluster.TopCluster_info[0][2];

		float minCoordWC_top = WaterCluster.TopCluster_coordInfo[idMaxWC_top][1];
		float maxCoordWC_top = WaterCluster.TopCluster_coordInfo[idMaxWC_top][0];
		float minCoordWC_bot = WaterCluster.BotCluster_coordInfo[idMaxWC_bot][1];
		float maxCoordWC_bot = WaterCluster.BotCluster_coordInfo[idMaxWC_bot][0];
		float maxCoordIC = BulkIceCluster.TopCluster_coordInfo[idMaxIC][0];
		float minCoordIC = BulkIceCluster.TopCluster_coordInfo[idMaxIC][1];

		float topInterface = maxCoordIC + 3.0;
		float botInterface = minCoordIC - 3.0;
		float topSurfaceSize = fabs(maxCoordWC_top - topInterface);
		float botSurfaceSize = fabs(botInterface - minCoordWC_bot);
		// ----------------------------------------------------------------------------------------------------------

		cout << "####QLL-INTERFACE (cluster) TOP: " << topInterface << " BOT: "
				<< botInterface << " Bulk: " << maxCoordIC << " " << minCoordIC
				<< " Water " << minCoordWC_top << " " << maxCoordWC_bot << endl;

		// Determine QLL molecules from QLL interface position-------------------------------------------------------
		q3f.qlloxygens(q3f.crystal, coord_norm, topInterface, botInterface);

		// 4) QLL Ice cluster ---------------------------------------------------------------------------------------
		cluster IceCluster(dangframe.getOxArrayVec(), q3f.qllCrystalOx,
				box_midpoint);
		IceCluster.connectivity();
		IceCluster.size(coord_norm, IfPrismPhase);
		// ----------------------------------------------------------------------------------------------------------
		clock_t t4_QLLClust = clock();

		// Cluster sizes: -------------------------------------------------------------------------------------------

		int sizeWC_top = WaterCluster.TopCluster_info[0][1];
		int sizeWC_bot = WaterCluster.BotCluster_info[0][1];
		int sizeQLLIC_top = IceCluster.TopCluster_info[0][1];
		int sizeQLLIC_bot = IceCluster.BotCluster_info[0][1];

		cout << "#MAIN Water-like Cluster Size. Top Surface: " << sizeWC_top
				<< " Bottom Surface: " << sizeWC_bot << endl;
		cout << "#MAIN QLL-Ice-like Cluster Size. Top Surface: "
				<< sizeQLLIC_top << " Bottom Surface: " << sizeQLLIC_bot
				<< endl;

		// ----------------------------------------------------------------------------------------------------------
		// ----------------- END CLUSTERS ANALYSIS:------------------------------------------------------------------
		// ----------------------------------------------------------------------------------------------------------

		
		
		// ----------------------------------------------------------------------------------------------------------
		// ----------------- DENSITY/Q3/ORIENTATION PROFILE: --------------------------------------------------------
		
		Density.Initialise(IsFirstFrame, IfPrismPhase, dcdf.getNATOM() / 3,
				coord_norm);

		Q3binAv.Average(q3f.q3or, Density.binVec_q3); 

		// Check if constant or dynamic QLL and act accordingly:
		if (IsConstQLL) {
			Q3binAv.TimeAverage(IsFirstFrame, Density.NumElementsVec_q3); 
			Density.TimeAverage(IsFirstFrame, Density.NumElementsVec);
		} else {
			Q3binAv.TimeAverage(q3f.getqllTopVec().size(), IsFirstFrame,
					IfPrismPhase, Density.NumElementsVec_q3);
			Density.TimeAverage(q3f.getqllTopVec().size(), IsFirstFrame,
					IfPrismPhase, Density.NumElementsVec);
		}
		// ----------------------------------------------------------------------------------------------------------
		clock_t t5_Q3DensProfile = clock(); 

		
		// ----------------------------------------------------------------------------------------------------------
		// ----------------- Dangling/QLL MOL STATS: ----------------------------------------------------------------
		
		// Average Orientation and q3 of Dangling Molecules:
		Averages_q3Orient.average(dangframe.getDangTopVec(),
				dangframe.getDangBotVec(), orient_frame.getOrient(), q3f.q3or);
		float OrientDangTop = Averages_q3Orient.avOrient_Top;
		float OrientDangBot = Averages_q3Orient.avOrient_Bot;
		float Q3DangTop = Averages_q3Orient.avQ3_Top;
		float Q3DangBot = Averages_q3Orient.avQ3_Bot;

		// Average Orientation and q3 of QLL molecules
		Averages_q3Orient.average(q3f.getqllTopVec(), q3f.getqllBotVec(),
				orient_frame.getOrient(), q3f.q3or);
		float OrientQLLTop = Averages_q3Orient.avOrient_Top;
		float OrientQLLBot = Averages_q3Orient.avOrient_Bot;
		float Q3QLLTop = Averages_q3Orient.avQ3_Top;
		float Q3QLLBot = Averages_q3Orient.avQ3_Bot;
		// ----------------------------------------------------------------------------------------------------------
		clock_t t6_AvOrientQ3 = clock(); 

		// ----------------------------------------------------------------------------------------------------------
		// ----------------- ORIENTATION PROFILE: -------------------------------------------------------------------
		OrientbinAv.Average(orient_frame.getOrient(), Density.binVec);
		if (IsConstQLL) {
			OrientbinAv.TimeAverage(IsFirstFrame, Density.NumElementsVec);
		} else {
			OrientbinAv.TimeAverage(q3f.getqllTopVec().size(), IsFirstFrame,
					IfPrismPhase, Density.NumElementsVec);
		}
		// ----------------------------------------------------------------------------------------------------------
		clock_t t7_OrienProfile = clock();

		// ----------------------------------------------------------------------------------------------------------
		// ----------------- IN PLANE STATS: ------------------------------------------------------------------------	
		// Determine density/q3/orientation distribution in-plane (ie in slices along z-axis)
		// Split plane into 45*45 grid and split System into slices
		// In each grid and slice calculate average orientation and average q3
		// Non Const QLL: Average over frames with similar qll size. 
		// Centre grid based on min and max coordinates in plane for each frame (safety measure in case system drifting laterally)
		// PBC conditions applied so don't have to worry about using a fixed position grid

		Density.InPlaneDensity(dcdf.getNATOM() / 3, Density.binVec_short,
				orient_frame.getOrient(), q3f.q3or, IfPrismPhase, coord_inplane,
				z);
		Density.InPlaneTimeAv(q3f.getqllTopVec().size(), IsConstQLL,
				IsFirstFrame, Density.binVec_short);

		// Calculate probability distribution of orientation in each slice:
		Density.orient_slice(q3f.getqllTopVec().size(), Density.binVec_short,
				orient_frame.getOrient(), coord_inplane, z, IsConstQLL,
				IsFirstFrame);

		binwidthZ = Density.binPlane_data[0] + binwidthZ;
		binwidthPlane = Density.binPlane_data[1] + binwidthPlane;
		count_planebin = count_planebin + 1.0;
		// ----------------------------------------------------------------------------------------------------------
		// ----------------- END IN PLANE STATS ---------------------------------------------------------------------
		clock_t t8_inPlane = clock(); 

		
		
		// ----------------------------------------------------------------------------------------------------------
		// ----------------- WRITE OUTPUT: --------------------------------------------------------------------------
		cout << "Write Output: " << endl;
		outf.output_dang(i + 1, numDangTop, numDangBot, OrientDangTop,
				OrientDangBot, Q3DangTop, Q3DangBot);
		outf.output_orient(i + 1, orient_frame.getOrient());
		outf.output_cluster(i + 1, sizeWC_top, sizeWC_bot, sizeQLLIC_top,
				sizeQLLIC_bot);
		outf.output_noxygens(i + 1, q3f.liquid.size(), q3f.crystal.size(),
				q3f.getqllTopVec().size(), q3f.getqllBotVec().size(),
				OrientQLLTop, OrientQLLBot, Q3QLLTop, Q3QLLBot);
		outf.interface(i + 1, topInterface, botInterface, topSurfaceSize,
				botSurfaceSize);
		outf.cluster_visual(BulkIceCluster.clusterId, WaterCluster.clusterId,
				idGlobalIC, idGlobalWC_top, idGlobalWC_bot);
		outf.q3_visual(q3f.liquid, q3f.crystal);
		outf.qll_visual(q3f.getqllTopVec(), q3f.getqllBotVec());
		// ----------------------------------------------------------------------------------------------------------
		// ----------------- COMPLETED OUTPUT: ----------------------------------------------------------------------

		// ----------------------------------------------------------------------------------------------------------
		// ----------------- CLEAR VECTORS and MEMORY: --------------------------------------------------------------
		Density.clear_vectors();
		OrientbinAv.clear_vectors();
		Q3binAv.clear_vectors();
		dangframe.clear_data();
		orient_frame.clear_data();
		q3f.clear_data();
		phase_check.clear_data();
		IceCluster.clear_data();
		WaterCluster.clear_data();
		// ----------------------------------------------------------------------------------------------------------
		// ----------------- COMPLETED CLEAR VECTORS ----------------------------------------------------------------

		// ----------------------------------------------------------------------------------------------------------
		// ----------------- TIMINGS: -------------------------------------------------------------------------------			
		clock_t end = clock();
		double ti = double(t_readData - begin) / CLOCKS_PER_SEC;
		double t1 = double(t1_statsDangOrient - t_readData) / CLOCKS_PER_SEC;
		double t2 = double(t2_WaterClust - t1_statsDangOrient) / CLOCKS_PER_SEC;
		double t3 = double(t3_IceClust - t2_WaterClust) / CLOCKS_PER_SEC;
		double t4 = double(t4_QLLClust - t3_IceClust) / CLOCKS_PER_SEC;
		double t5 = double(t5_Q3DensProfile - t4_QLLClust) / CLOCKS_PER_SEC;
		double t6 = double(t6_AvOrientQ3 - t5_Q3DensProfile) / CLOCKS_PER_SEC;
		double t7 = double(t7_OrienProfile - t6_AvOrientQ3) / CLOCKS_PER_SEC;
		double t8 = double(t8_inPlane - t7_OrienProfile) / CLOCKS_PER_SEC;
		double t9 = double(end - t8_inPlane) / CLOCKS_PER_SEC;
		double total = double(end - begin) / CLOCKS_PER_SEC;

		cout << "####TIMINGS: Read Data: " << ti << endl;
		cout << "####TIMINGS: 1) Dangling Orientation Stats: " << t1
				<< "    2) WaterCluster: " << t2 << "    3) IceCluster: " << t3
				<< endl;
		cout << "####TIMINGS: 4) QLL Cluster: " << t4
				<< "    5) Q3 & Density Profile: " << t5
				<< "    6) Average Orient & Q3 (QLL/Dang): " << t6 << endl;
		cout << "####TIMINGS: 7) Orient Profile: " << t7 << "    8) InPlane : "
				<< t8 << "   9) Output: " << t9 << "  TOTAL: " << total << endl;
		cout << " " << endl;
		// ----------------------------------------------------------------------------------------------------------
		// ----------------- COMPLETED TIMINGS ----------------------------------------------------------------------

		// ----------------------------------------------------------------------------------------------------------
		// ----------------- COMPLETED ANALYSIS OF FRAME ------------------------------------------------------------
		// ----------------- LOOP TO NEXT FRAME ---------------------------------------------------------------------
		// ----------------------------------------------------------------------------------------------------------
		IsFirstFrame = false;

	} 
	
	// --------------------------------------------------------------------------------------------------------------
	// ----------------- COMPLETED LOOP OVER ALL FRAMES -------------------------------------------------------------
	// --------------------------------------------------------------------------------------------------------------

	// --------------------------------------------------------------------------------------------------------------
	// ----------------- WRITE OUTPUT -------------------------------------------------------------------------------
	// --------------------------------------------------------------------------------------------------------------

	binwidthZ = binwidthZ / count_planebin;
	binwidthPlane = binwidthPlane / count_planebin;

	float b_width = Density.binDataArray[2];
	float b_width_q3 = Density.binDataArray[4];
	Density.orient_sliceSUM(IsConstQLL);

	if (IsConstQLL) {
		Density.NormaliseHeatmaps();
		bool IsDensity = false;
		outf.AverageTimeProfile("ProfileOrient.txt",
				OrientbinAv.getTimeBinAvVec(IsDensity), b_width);
		outf.AverageTimeProfile("ProfileQ3.txt",
				Q3binAv.getTimeBinAvVec(IsDensity), b_width_q3);
		IsDensity = true;
		outf.AverageTimeProfile("ProfileDensity.txt",
				Density.getTimeBinAvVec(IsDensity), b_width);

		outf.InPlane("Heatmap_OrientationThetaGamma.txt", "Heatmap_Density.txt",
				"Heatmap_Q3.txt", "Heatmap_DipoleTheta.txt",
				Density.getAvOrientPlaneVec(), Density.getAvPlaneVec(),
				binwidthZ, binwidthPlane);
		outf.orientSlice("ProfileSlice_dipoleTheta.txt", Density.OrientSlice);

	} else {
		Density.NormaliseHeatmaps_nQLL();
		bool IsDensity = false;
		outf.AverageTimeProfile("ProfileOrient.txt",
				OrientbinAv.getTimeAv_nQLLVec(IsDensity),
				OrientbinAv.getTimeAv_countnQLL(), b_width);
		outf.AverageTimeProfile("ProfileQ3.txt",
				Q3binAv.getTimeAv_nQLLVec(IsDensity),
				Q3binAv.getTimeAv_countnQLL(), b_width_q3);
		IsDensity = true;
		outf.AverageTimeProfile("ProfileDensity.txt",
				Density.getTimeAv_nQLLVec(IsDensity),
				Density.getTimeAv_countnQLL(), b_width);

		outf.InPlane("Heatmap_OrientationThetaGamma.txt", "Heatmap_Density.txt",
				"Heatmap_Q3.txt", "Heatmap_DipoleTheta.txt",
				Density.getOrientAvPlane_nQLL(), Density.getAvPlane_nQLL(),
				Density.getTimeAv_countnQLL(), binwidthZ, binwidthPlane);
		outf.orientSlice("ProfileSlice_dipoleTheta.txt",
				Density.OrientSlice_nqll, Density.getTimeAv_countnQLL());

	}

	return EXIT_SUCCESS;
}
