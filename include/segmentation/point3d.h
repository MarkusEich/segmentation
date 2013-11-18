#ifndef __point3d_h__
#define __point3d_h__

#include <math.h>
#include <iostream>
#include <stdlib.h>
#include <stdint.h>
#include <boost/shared_ptr.hpp>

namespace Pointcloud
{

inline double sq(double a) {
    return a*a;
}

class Point_3d
{

public:

	Point_3d(double x=0,double y=0, double z=0, uint8_t alpha=0, uint8_t r=0, uint8_t g=0, uint8_t b=0):
		x(x),
		y(y),
		z(z),
		color_alpha(alpha),
		color_red(r),
		color_green(g),
		color_blue(b)

    	{
    	}

	double x;
	double y;
	double z;

	uint8_t color_alpha;
	uint8_t color_red;
	uint8_t color_green;
	uint8_t color_blue;

	friend std::ostream& operator<<(std::ostream& os,const Point_3d& p)
	{
		os << p.x<<" "<<p.y<<" "<<p.z;
		return os;
	}

	/**
	 * Inline (only for the optimization), this method will be copied if possible by the 
	 * compiler, where it is called.
	*/
	inline Point_3d& operator=(const Point_3d& p2)
	{
		x=p2.x;
		y=p2.y;
		z=p2.z;
		color_alpha = p2.color_alpha;
		color_blue = p2.color_blue;
		color_green = p2.color_green;
		color_red = p2.color_red;
		return *this;
	}

	inline bool operator<(const Point_3d &p2) const
	{
		double dist1 = fabs(p2.x) + fabs(p2.y) + fabs(p2.z);
		double dist2 = fabs(x) + fabs(y) + fabs(z);
		if (dist1>dist2)
			return true;
		return false;
	}


	inline bool operator==(const Point_3d &p2) const
	{
	  double d = sqrt (sq(p2.x-x) + sq(p2.y-y) + sq(p2.z-z));
		return d < 0.001;
	}

	inline bool operator!=(const Point_3d &p2) const
	{
		return x!=p2.x || y!=p2.y || z!=p2.z;
	}

	double distance(Point_3d const& p) const
	{
		double dist = 0;
		dist += (x-p.x)*(x-p.x);
		dist += (y-p.y)*(y-p.y);
		dist += (z-p.z)*(z-p.z);
		return sqrt(dist);
	}

	inline double length() {
		return sqrt(pow(x,2)+pow(y,2)+pow(z,2));
	}

};
	   
inline double elementValue(Point_3d *t, size_t k ) { 
	if(t == 0){
		std::cerr<<"Cannot handle NULL pointer %s\n" <<__FUNCTION__<<std::endl;
		exit(-100);
	}
	if(k==0) return t->x; 
	if(k==1) return t->y; 
	if(k==2) return t->z; 
	return -1;
}

class Face{
public:
	Point_3d p1;
	Point_3d p2;
	Point_3d p3;

};


}

typedef boost::shared_ptr<Pointcloud::Point_3d>   Shared_Point_3d; 

#endif //__pointcloud_h__
