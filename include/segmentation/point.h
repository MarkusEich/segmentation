#ifndef _rg_point_h__
#define _rg_point_h__

#include <math.h>
#include <iostream>
#include <stdlib.h>
#include <list>
#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>

namespace segmentation{

namespace {
	/**
	 * Square of a number.
	 */
	inline double sq(double a) {
		return a*a;
	}
}


class Point
{
public:

	/**
	 * Mini-class to store the nearest neighbors of a point.
	 */
	class Neighbour {
	public:
		/**
		 * Saves point and distance.
		 */
		Neighbour(boost::shared_ptr<Point> p, double d)
			: point(p), distance(d) {}
		~Neighbour(){
			point.reset();
		}

		boost::weak_ptr<Point> point;
		double distance;
		/**
		 * Comparison of two neighbors based of their distance from the source point.
		 */
		bool operator<(const Neighbour& other) const {
			return distance < other.distance;
		}
	};

	/**
	 * Generates point, take values from the constructor.
	 */
	Point(double x=0,double y=0, double z=0) :
		x(p[0]),	//add references to the array elements, must not rewrite any code, witch directly accesses on x, y, z 
            	y(p[1]),
            	z(p[2])
    	{

		p[0] = x;
		p[1] = y;
		p[2] = z;
    	}

    	~Point(){
    		nearestNeighbours.clear();
    	}

	/**
	 * List of the nearst neighbours.
	 */
	std::list<Neighbour> nearestNeighbours;

	/**
	 * Coordinates of the point.
	 */
	double p[3];

	/**
	 * Angle of the point.
	 */
	double angle;

	/**
	 * Other names for the above coordinates.
	 */
	double &x;
	double &y;
	double &z;
	size_t index;
    
	inline Point& operator=(const Point& p2)
    	{
		p[0] = p2.p[0];
		p[1] = p2.p[1];
		p[2] = p2.p[2];
		x=p[0];
		y=p[1];
		z=p[2];
		angle=p2.angle;
		index=p2.index;
        	return *this;
    	}
    
    	inline bool operator==(const Point &p2) const
    	{ 

        return p[0]==p2.p[0] && p[1]==p2.p[1] && p[2]==p2.p[2];

    	}

	/**
	 * Distance of the current point to the given point.
	 */
    	inline double distance(Point const& p) const
    	{
        	double dist = 0;
        	dist += (x-p.x)*(x-p.x);
        	dist += (y-p.y)*(y-p.y);
        	dist += (z-p.z)*(z-p.z);
        	return sqrt(dist);
    	}


    	/**
	 * Compares 2 angle, auxiliary function for priority_queue.
	 */
    	bool operator() (const Point *p, const Point *pp) const {
        	return p->angle<pp->angle;
    	}


};

typedef boost::shared_ptr<Point>  Shared_Point;
}
#endif //__pointcloud_h__
