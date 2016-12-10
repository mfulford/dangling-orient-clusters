#include <fstream>
#include <vector>

#ifndef Q3_R_INCLUDED 
#define Q3_R_INCLUDED

class q3_r {
protected:
	std::fstream q3f; 
	int num_ox;
	
	std::vector<int> qllCrystalTop;
	std::vector<int> qllCrystalBot;
	std::vector<int> qllLiquidTop;
	std::vector<int> qllLiquidBot;
	std::vector<int> qllLiquidOx;
	std::vector<int> qllTop;
	std::vector<int> qllBot;

	float threshold = 0.670; 

public:

	float topFace, botFace;

	std::vector<float> q3or;

	std::vector<int> liquid;    	// vector with id of Ox that are liquid-like     
	std::vector<int> crystal;     	// vector with id of Ox that are ice-like
	std::vector<int> qllCrystalOx;
	std::vector<std::vector<int> > LiquidConnectVec;
	std::vector<std::vector<int> > qllCrystalConnectVec;
	std::vector<std::vector<int> > PhaseConnectVec;

	std::vector<int> getqllCrystalTopVec() const {
		return qllCrystalTop;
	}
	;
	std::vector<int> getqllCrystalBotVec() const {
		return qllCrystalBot;
	}
	;
	std::vector<int> getqllLiquidTopVec() const {
		return qllLiquidTop;
	}
	;
	std::vector<int> getqllLiquidBotVec() const {
		return qllLiquidBot;
	}
	;

	std::vector<int> getqllTopVec() const {
		return qllTop;
	}
	;
	std::vector<int> getqllBotVec() const {
		return qllBot;
	}
	;

	std::vector<int> bulk;

	q3_r(const char q3filename[], const int, const float);
	void q3_frame();

	void qlloxygens(const std::vector<int> &, const float*, const float,
			const float);

	~q3_r();
	void clear_data();
};

#endif
