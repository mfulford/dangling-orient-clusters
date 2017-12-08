#include <vector>

#ifndef CLUSTER_INCLUDED 
#define CLUSTER_INCLUDED

class cluster{
	enum VertexState { Unvisited, Visiting, Visited};
	std::vector<VertexState> state;
	int count;

	const std::vector<int> &OxygenPhase;
	std::vector< std::vector<int> > PhaseConnectVec;
	std::vector <std::vector<int> > GlobalConnectVec; // making a copy of Odist_Array

	std::vector<int> vertexid;

	float BoxCentre;
public:
	std::vector< std::vector<int> > clusterId;

	int Max,SecMax;      
	int Size1Bot, Size2Bot, Size1Top, Size2Top; 

	std::vector< std::vector <int> > TopCluster_info;
	std::vector< std::vector <float> > TopCluster_coordInfo;
	std::vector< std::vector <int> > BotCluster_info;
	std::vector< std::vector <float> > BotCluster_coordInfo;


	cluster(const std::vector <std::vector <int> > &GlobalConnect, std::vector<int>  &Phase, const float boxcentre) : GlobalConnectVec(GlobalConnect), OxygenPhase(Phase), BoxCentre(boxcentre)  {}
	void connectivity();

	void size(const float*, std::vector<bool>);
	void dfsfromVertex(int);
	void clear_data();
};



#endif



