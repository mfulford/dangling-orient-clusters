#include "dataCell_pbc.hpp"
#include <fstream>
#include <vector>


#ifndef DANGLING_CHECK_INCLUDED 
#define DANGLING_CHECK_INCLUDED

class dangling : public pbc_check{
private:

	std::vector< std::vector<int> > danglinglist;                    // 2D vector -- 2 rows: each H -- n columns
	std::vector< std::vector<int> > DanglingTop;                     // 2D vector -- 2 rows: each H -- n columns
	std::vector< std::vector<int> > DanglingBot;                     // 2D vector -- 2 rows: each H -- n columns

	int num_ox;
	std::vector< std::vector<int> > Odist_Array;

public:
	dangling(int);   

	void danglingOx(const int, const float*,const float*, const float*, std::vector<bool>, const float);
	std::vector< std::vector<int> > getDangVec() const;	// returns danglinglist vector
	std::vector< std::vector<int> > getDangTopVec() const {return DanglingTop;};     	// returns DanglingTop vector        
	std::vector< std::vector<int> > getDangBotVec() const {return DanglingBot;};     	// returns DanglingBot vector

	std::vector< std::vector<int> > getOxArrayVec() const;	// returns Odist_Array vector
	void clear_data();
	~dangling();
};

#endif
