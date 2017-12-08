#include <vector>
#include "dataCell_pbc.hpp"

#ifndef ORIENTATION_CHECK_INCLUDED
#define ORIENTATION_CHECK_INCLUDED


class orientation : public pbc_check{
protected:
	std::vector< std::vector<float> > orient;  // 2 * num_ox vector. 2 for each H.

public:
	float avTop, avBot;

	orientation();
	void orientationOH(const int, const int, const float*,const float*, const float*, std::vector<bool>, const float);
	std::vector< std::vector<float> > getOrient() const {return orient;};

	void average(const std::vector< std::vector<int> > &, const std::vector< std::vector<int> > &);
	void average(const std::vector<int> &, const std::vector<int> &);


	void clear_data();
	~orientation();	
};


#endif

