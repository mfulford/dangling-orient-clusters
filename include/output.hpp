#include <fstream>
#include <vector>

#ifndef OUTPUT_CHECK_INCLUDED 
#define OUTPUT_CHECK_INCLUDED

class Output {
protected:
	const char* dangfile;
	const char* orientfile;
	const char* clusterfile;
	const char* binfile;
	const char* noxygenfile;
	const char* interfacefile;
	const char* clusVisfile;
	const char* q3Visfile;
	const char* qllVisfile;
	const char* orientSlicefile;

	int num_atoms;
	int width_qll;

	std::ofstream dangf; // file with number of dangling OH and average orientation for each timestep	
	std::ofstream orientf; // file with orientation of each OH bond in every frame
	std::ofstream clusterf;
	std::ofstream binAvf;
	std::ofstream noxf;
	std::ofstream planef_Top;
	std::ofstream planef_Bot;
	std::ofstream planef_dens;
	std::ofstream planef_q3;
	std::ofstream planef_orient;
	std::ofstream planef_angles; //theta rho
	std::ofstream interfacef;
	std::ofstream ClusVisf;
	std::ofstream q3Visf;
	std::ofstream qllVisf;
	std::ofstream orientSlicef;

public:
	Output(const char* dangfile, const char* orientfile,
			const char* clusterfile, const char* noxygenfile,
			const char* interfacefile, const char* clusVisfile,
			const char* q3Visfile, const char* qllVisfile, const int natoms,
			const int num_frames, bool);
	void output_dang(const int frame, const int num_dangTop,
			const int num_dangBot, const float OrientTopDang,
			const float OrientBotDang, const float Q3TopDang,
			const float Q3BotDang);
	void output_orient(const int frame,
			const std::vector<std::vector<float> > orientvec);
	void output_cluster(int, int, int, int, int);
	void interface(int, float, float, float, float);
	void cluster_visual(std::vector<std::vector<int> >,
			std::vector<std::vector<int> >, int, int, int);
	void q3_visual(std::vector<int>, std::vector<int>);
	void qll_visual(std::vector<int>, std::vector<int>);

	void AverageTimeProfile(const char* binfile, const std::vector<float>,
			const float);
	void AverageTimeProfile(const char* binfile,
			const std::vector<std::vector<float> >, const int *, const float);
	void output_noxygens(const int frame, const int ntotLiquid,
			const int ntotIce, const int nqll_top, const int nqll_bot,
			const float OrientTopQLL, const float OrientBotQLL,
			const float Q3TopQLL, const float Q3BotQLL);

	void InPlane(const char*, const char*, const char*, const char*,
			const std::vector<std::vector<std::vector<float> > > &,
			const std::vector<std::vector<std::vector<std::vector<float> > > > &,
			const float, const float); // constant QLL
	void InPlane(const char*, const char*, const char*, const char*,
			const std::vector<std::vector<std::vector<std::vector<float> > > > &,
			const std::vector<
					std::vector<std::vector<std::vector<std::vector<float> > > > > &,
			const int*, const float, const float); // non const QLL

	void orientSlice(const char*, const std::vector<std::vector<float> > &);
	void orientSlice(const char*,
			const std::vector<std::vector<std::vector<float> > > &, const int*);

	~Output();

};

#endif
