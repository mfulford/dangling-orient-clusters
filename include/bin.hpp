#include <vector>

#ifndef BIN_CHECK_INCLUDED 
#define BIN_CHECK_INCLUDED

class bin{
private:
	std::vector <float> densityVec;

	std::vector <std::vector <std::vector <std::vector<float> > > > PlaneVec;
	std::vector <std::vector <std::vector<float> > > orientPlaneVec;
	std::vector <std::vector <std::vector <std::vector<float> > > > avPlaneVec;
	std::vector <std::vector <std::vector<float> > > avOrientPlaneVec;
	std::vector <std::vector <std::vector <std::vector <std::vector<float> > > > > AvPlane_nQLL;
	std::vector <std::vector <std::vector <std::vector<float> > > > OrientAvPlane_nQLL;
	std::vector<bool> Unvisited_nQLL;	

	std::vector <float> TimeBinAvVec;
	std::vector <std::vector<float> > TimeAv_nQLL;	

	int nqllbins;
	int nbins_short;
	int nbins;
	int nbins_q3;
	float maxbin, minbin, range, binwidth;
	int width_qll;
	int num_slices;
	float ratio_bins, ratio_q3;	
	int *count_nQLL;
	int n_slice_bins;

	std::vector<int> countSlice;
	std::vector<std::vector<int> > countSlice_nqll;


public:
	std::vector <float> BinAvVec;
	std::vector< std::vector <int> > binVec;
	std::vector< std::vector <int> > binVec_short,  binVec_q3;
	std::vector<int> tot_NumElementsVec;
	std::vector<std::vector<int> > tot_NumElementsVec_nqll;
	std::vector<int> NumElementsVec, NumElementsVec_short, NumElementsVec_q3;

	std::vector<std::vector<float> > OrientSlice;
	std::vector<std::vector<std::vector<float> > > OrientSlice_nqll;


	float  binDataArray[5];
	float binPlane_data[2];
	float Initialise(bool, int, const float*);
	void Average(const std::vector<float> &, const std::vector<std::vector <int> > &);
	void Average(const std::vector< std::vector<float> > &, const std::vector<std::vector <int> > &);

	void TimeAverage(bool, const std::vector<int> &);
	void TimeAverage(int, bool, std::vector<bool>, const std::vector<int> &);
	void InPlaneTimeAv(int, bool, bool, const std::vector<std::vector <int> > &);
	void InPlaneDensity(int,const std::vector< std::vector<int> > &, const std::vector< std::vector<float> > &, const std::vector<float> &, const float*, const float*);

	void NormaliseHeatmaps();
	void NormaliseHeatmaps_nQLL();

	const std::vector<float> getTimeBinAvVec(bool);
	const std::vector<std::vector<float> > getTimeAv_nQLLVec(bool);
	const int* getTimeAv_countnQLL() const  {return count_nQLL;};
	const std::vector<std::vector<std::vector<std::vector<float> > > > getAvPlaneVec();
	const std::vector<std::vector<std::vector<float> > > getAvOrientPlaneVec(); 
	const std::vector<std::vector<std::vector<std::vector<std::vector<float> > > > > getAvPlane_nQLL();
	const std::vector<std::vector<std::vector<std::vector<float> > > > getOrientAvPlane_nQLL();

	std::vector<float> const getDensityVec() {return densityVec;};	

	void orient_slice(int, const std::vector<std::vector<int> > &, const std::vector<std::vector<float> > &, const float*, const float*, bool,bool);
	void orient_sliceSUM(bool);


	void clear_vectors();

};

#endif
