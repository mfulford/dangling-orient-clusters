#include <iostream>     
#include <vector>
#include <algorithm>    
#include <numeric> 
#include <math.h>  
#include "cluster.hpp"

using namespace std;

/* Determine clusters of groups of molecules.
 * Use depth-first search algorithm to efficiently determine clusters
 */

// Determine connectivity matrix
void cluster::connectivity() {

	int phase_size = OxygenPhase.size();

	PhaseConnectVec.resize(phase_size, vector<int>(0, 0)); // vector has number of rows = number of [phase] Ox. Number of columns=0. 

	// loop over all liquid Ox or crystal qll Ox:
	for (int i = 0; i < phase_size; i++) {

		if (GlobalConnectVec[OxygenPhase[i]].size() == 0) {
			cout << " ############### Num neighbours =0. i: " << i << endl;
		}

		// loop over neighbours of Ox i:
		for (int j = 0; j < GlobalConnectVec[OxygenPhase[i]].size(); j++) {

			int GlobConnect_id = GlobalConnectVec[OxygenPhase[i]][j];
			// loop over OxygenPhase to see if neighbours j of i are liquid/crystal-qll:
			for (int k = 0; k < phase_size; k++) {

				// check if neighbour k is [phase]:
				if (GlobConnect_id == OxygenPhase[k]) {
					if (OxygenPhase[i] == OxygenPhase[k]) {
						cout
								<< "####CLUSTERconnect ERROR in Connectivity: i and k "
								<< i << " " << k << endl;
					}
					PhaseConnectVec.at(i).push_back(k); //IMPORTANT push back k, not liquid[k]
					break;
				}

			} // end loop checking if neighbours j are [phase]

		} // end loop over neighbours j of Ox i

	} // end loop over Ox i

} // end function

// OxygenId = liquid 
// OxygenConnectVec= LiquidConnectVec
// Determine cluster and size:
void cluster::size(const float* coord_norm, bool IsPrism) { // read in LiquidConnectVec (matrix connectivity)  and liquid (which tells Ox id of row number in LiquidConnectVec) 

	state.clear();
	state.resize(OxygenPhase.size(), Unvisited);
	vector<int> size;
	size.clear(); 
	int sum = 0;

	cout << "#CLUSTERsize Cluster count: ";

	for (int v = 0; v < OxygenPhase.size(); v++) {           // for all vertices

		if (state[v] == Unvisited) {

			count = 0;
			dfsfromVertex(v);
			size.push_back(count);
			clusterId.push_back(vertexid);
			vertexid.clear();
			sum = sum + count;
			// output size of all clusters greater than 5
			if (count > 5) {
				cout << count << " ";
			}
		}

	}
	cout << endl;
	partial_sort(size.begin(), size.begin() + 2, size.end(), greater<int>());
	Max = size[0];
	SecMax = size[1];
	cout << "#CLUSTERsize Largest cluster size: " << Max << " Second Largest: "
			<< SecMax << endl;
	cout << "#CLUSTERsize Sum of All cluster sizes: "
			<< accumulate(size.begin(), size.end(), 0) << endl;

	vector<float> ClusterCoordVec;
	int count_clusterTop = 0;
	int count_clusterBot = 0;

	float minCoord, maxCoord;

	for (int i = 0; i < clusterId.size(); i++) {
		int ClusterSize = clusterId[i].size();

	
		/*---------------------------------------------------------
		* Populate ClusterCoordVec with "z" coordinates
		* Determine min/max position of each cluster. Use this later to determine interface position
		* Depending on size of cluster, determine min/max position by averaging of num_atoms number of positions
		* Shape taken into account by averaging over a cutoff. If min position - max postion > cut-off, 
		* then discard outermost atom and retry
		*/
		
		for (int j = 0; j < ClusterSize; j++) {
			ClusterCoordVec.push_back(coord_norm[clusterId[i][j] * 3]); /*push back y coordinate of Oxygen memebers of Max cluster */
		}

		int num_atoms;
		float cutoff;
		if (ClusterSize < 2000) {
			minCoord = *min_element(ClusterCoordVec.begin(),
					ClusterCoordVec.end());
			maxCoord = *max_element(ClusterCoordVec.begin(),
					ClusterCoordVec.end());
		} else {

			//cutoff=3.0;
			//if ( (ClusterSize >=30) && (ClusterSize < 70)) {num_atoms=5;}
			//else if ( (ClusterSize >=70) && (ClusterSize < 100)) {num_atoms=10;}
			//else if ( (ClusterSize >=100) && (ClusterSize < 150)) {num_atoms=15;}
			//else if ( (ClusterSize >= 150) && (ClusterSize < 550)) {num_atoms = 20; cutoff = 1.6;}                                                          //num_atoms=20; cutoff =1.0;}
			//else if ( (ClusterSize >= 550) && (ClusterSize < 1000)) {num_atoms = 100; cutoff = 2.0;}                                //num_atoms=50; cutoff=1.6; }
			//else if ( (ClusterSize >= 1000) && (ClusterSize < 2000)) {num_atoms=40; cutoff=1.6;} //num_atoms=75;}
			if (ClusterSize >= 2000) {
				if (IsPrism) {
					num_atoms = 200;
					cutoff = 4.0;
				} else {
					num_atoms = 200;
					cutoff = 4.0;
				} 	// 110; 1.5
			}
			maxCoord = 0.0;
			minCoord = 0.0;

			int i_top = 0;
			int i_bot = num_atoms;

			sort(ClusterCoordVec.begin(), ClusterCoordVec.end(),
					greater<float>());
			bool topComplete = false;
			while (topComplete == false) {
				float max_tmp = ClusterCoordVec[i_top];
				float min_tmp = ClusterCoordVec[i_bot];
				float diff = fabs(max_tmp - min_tmp);
				if (diff > cutoff) {
					i_top = i_top + 1;
					i_bot = i_bot + 1;
				} else {
					topComplete = true;
				}
			}

			if (ClusterSize >= 2000) {
				if (IsPrism) {
					num_atoms = 110;
					cutoff = 1.8;
				} else {
					num_atoms = 130;
					cutoff = 1.8;
				}
				i_bot = num_atoms + i_top;
				topComplete = false;
				while (topComplete == false) {
					float max_tmp = ClusterCoordVec[i_top];
					float min_tmp = ClusterCoordVec[i_bot];
					float diff = fabs(max_tmp - min_tmp);
					if (diff > cutoff) {
						i_top = i_top + 1;
						i_bot = i_bot + 1;
					} else {
						topComplete = true;
					}
				}
			}
			for (int ii = i_top; ii < i_bot; ii++) {
				maxCoord = ClusterCoordVec[ii] + maxCoord;
			}    
			maxCoord = maxCoord / num_atoms;

			if (ClusterSize >= 2000) {
				if (IsPrism) {
					num_atoms = 200;
					cutoff = 4.0;
				} else {
					num_atoms = 200;
					cutoff = 4.0;
				}
			}

			i_top = 0;
			i_bot = num_atoms;
			bool botComplete = false;
			sort(ClusterCoordVec.begin(), ClusterCoordVec.end(), less<float>());

			while (botComplete == false) {
				float max_tmp = ClusterCoordVec[i_bot];
				float min_tmp = ClusterCoordVec[i_top];
				float diff = fabs(max_tmp - min_tmp);	// float absolute
				if (diff > cutoff) {
					i_top = i_top + 1;
					i_bot = i_bot + 1;
				} else {
					botComplete = true;
				}
			}

			if (ClusterSize >= 2000) {
				if (IsPrism) {
					num_atoms = 110;
					cutoff = 1.8;
				} else {
					num_atoms = 130;
					cutoff = 1.8;
				}
				i_bot = i_top + num_atoms;
				botComplete = false;
				while (botComplete == false) {
					float max_tmp = ClusterCoordVec[i_bot];
					float min_tmp = ClusterCoordVec[i_top];
					float diff = fabs(max_tmp - min_tmp);
					if (diff > cutoff) {
						i_top = i_top + 1;
						i_bot = i_bot + 1;
					} else {
						botComplete = true;
					}
				}
			}

			for (int ii = i_top; ii < i_bot; ii++) {
				minCoord = ClusterCoordVec[ii] + minCoord;
			}
			minCoord = minCoord / num_atoms;
		}

		//-----------------------------------------------------------
		vector<int> clus_inf;
		vector<float> clus_coorInf;
		clus_coorInf.push_back(maxCoord);
		clus_coorInf.push_back(minCoord);

		// Determine if cluster is in top or bottom half
		// NB largest bulk Ice-like cluster will be in top half based on this definition
		// NB largest top qll ice-like cluster is 2nd largest in top half
		// NB largest bot qll ice-like cluster is 1st largest in bottom half 
		if (maxCoord > BoxCentre) {
			clus_inf.push_back(count_clusterTop);
			clus_inf.push_back(ClusterSize);
			clus_inf.push_back(count_clusterTop + count_clusterBot); // this gives ID of cluster (in clusterId)
			TopCluster_info.push_back(clus_inf); /* push back cluster id number and size of cluster */
			TopCluster_coordInfo.push_back(clus_coorInf);
			count_clusterTop = count_clusterTop + 1;
		} else {
			clus_inf.push_back(count_clusterBot);
			clus_inf.push_back(ClusterSize);
			clus_inf.push_back(count_clusterTop + count_clusterBot);
			BotCluster_info.push_back(clus_inf);
			BotCluster_coordInfo.push_back(clus_coorInf);
			count_clusterBot = count_clusterBot + 1;
		}

		clus_inf.clear();
		clus_coorInf.clear();

		ClusterCoordVec.clear();
	}

	sort(TopCluster_info.begin(), TopCluster_info.end(),
			[](const vector< int >& a, const vector< int >& b) {return a[1] > b[1];});
	sort(BotCluster_info.begin(), BotCluster_info.end(),
			[](const vector< int >& a, const vector< int >& b) {return a[1] > b[1];});

	if (TopCluster_info.size() == 0) {
		vector<int> tmp_clus_inf;
		tmp_clus_inf.push_back(10000);
		tmp_clus_inf.push_back(0);
		tmp_clus_inf.push_back(10000);
		TopCluster_info.push_back(tmp_clus_inf);
	}
	if (BotCluster_info.size() == 0) {
		vector<int> tmp_clus_inf;
		tmp_clus_inf.push_back(10000);
		tmp_clus_inf.push_back(0);
		tmp_clus_inf.push_back(10000);
		BotCluster_info.push_back(tmp_clus_inf);
	}

	// Top/BotCluster_info 1st col is id of cluster. 2nd col is size. Row is largest to smallest (Largest is [0][1] Third largest is [2][1]) 

}

// Depth-first search
void cluster::dfsfromVertex(int v) {
	int neigh_id;

	state[v] = Visiting;
	count = count + 1;
	vertexid.push_back(OxygenPhase[v]);	// Ox id = OxygenPhase[v].   

	for (int w = 0; w < PhaseConnectVec[v].size(); w++) { // for all neighbours of v

		neigh_id = PhaseConnectVec[v][w]; // id of Ox ---> actual Ox is liquid[id]

		if (state[neigh_id] == Unvisited) {
			dfsfromVertex(neigh_id);
		}
	}
	state[v] = Visited;
}

void cluster::clear_data() {
	vertexid.clear();
	clusterId.clear();
	GlobalConnectVec.clear();
	PhaseConnectVec.clear();
	state.clear();
	TopCluster_info.clear();
	BotCluster_info.clear();
	TopCluster_coordInfo.clear();
	BotCluster_coordInfo.clear();
}

