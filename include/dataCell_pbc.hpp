#include <fstream>
#include <vector>

#ifndef DATACELL_CHECK_INCLUDED 
#define DATACELL_CHECK_INCLUDED

//global float. cell dimensions from cell data
extern float pbc_x;
extern float pbc_y;
extern float pbc_z;


class dataCell{

protected:
	std::fstream cellf;
public:
	dataCell(const char cellfilename[]);
};

class pbc_check{ 
private:
	float x_bound, y_bound, z_bound, dist;
	float box_middle;       	

	float xcoord, ycoord, zcoord, distx_sqr, disty_sqr, distz_sqr;

public:
	pbc_check();

	float pbc_xcoord(const float);
	float pbc_ycoord(const float);
	float pbc_zcoord(const float); 
	float getdistx() const;
	float getdisty() const;
	float getdistz() const;
	std::vector<float> h_midpoint(const float, const float, const float, const float, const float, const float);


	float pbc_xyz(const float, const float, const float);
	/* const after a function declaration means the function is not allowed to change any class members
	 * (except ones that are marked mutable). */
	float getdist_pbc() const;
	float getxpbc() const;
	float getypbc() const;
	float getzpbc() const;
	std::vector<bool> Phase();

	float box_centre(int,const float*);

	void clear_data(); 
	virtual ~pbc_check();

};

class Stats{
public:
	float avOrient_Top, avOrient_Bot, avQ3_Top, avQ3_Bot;
	void average(const std::vector< std::vector<int> > &, const std::vector< std::vector<int> > &,
			     const std::vector< std::vector<float> > &, const std::vector<float> &);
	void average(const std::vector<int> &, const std::vector<int> &,  const std::vector< std::vector<float> > &,
			     const std::vector<float> &);
};


#endif

