/*****************************************
 ** Markus Eich (2013)
 ** Malgorzata Goldhoorn (2010)
 ** DFKI GmbH
 **
 **
*/

#ifndef REGIONGROWING_H_
#define REGIONGROWING_H_

#include <vector>
#include <set>
#include <queue>
#include <functional>
#include <string.h>
#include <boost/thread/recursive_mutex.hpp>
#include <Eigen/Geometry>
#include <Eigen/Core>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <segmentation/point3d.h>
#include <segmentation/point.h>

#ifndef TRANSPARENCY
#define TRANSPARENCY 255
#endif


namespace segmentation{

class Point;
class Plane;


class RegionGrowing{

public:
	/**
	 * Creates a region growing object.
	 */
	RegionGrowing();
	virtual ~RegionGrowing();

    	static void calculateNeighbours(std::vector<Shared_Point> &points, double delta);

	/**
	 * Searchs regions, takes time.
	 */
	bool regionFind(std::list<boost::shared_ptr<Plane> > &result);
	
	/**
	 * Sets the callback function, which is called everytime plane is found.
	 */
	void setCallbackFunction(void(*callback)(boost::shared_ptr<Plane>)){
		this->funcPtr=callback;
	}

	//same as above but using boost bind
	void setBoostCallback(boost::function<void(boost::shared_ptr<Plane>)> const &cb)
	  {
	     callback_ = cb;
	  }

	void setPointcloud(Pointcloud::Point_3d* points, const int count);

	void setPointcloud(std::vector<Shared_Point>& points);

	void setGamma(double d);

	void setTheta(int d);

	void setEpsilon(double d);

	void setDelta(double d);

	void setNeighbours(int v);

 protected:

	/**
	 * Maximum distance between each point and its neighbors.
	 */
	double delta; 

	/**
	 * Maximum distance between points and the plane.
	 */
 	double gamma; 

	/**
	 * Maximum MSE for planes.
	 */
	double epsilon; 

	/**
	 * Minimum number of points in the plane.
	 */
	int theta;

	/**
	 * Number of nearest neighbors, which are considered.
	 */
	static int neighbours_count;
 private:

	/**
	 * Pre-calculating of neighbours_count nearest points with maximal distance delta.
	 */
	void precomputeNearestNeighbours();

	/**
	 * Number of already segmented points.
	 */
	size_t n_p;

	/**
	 * Entire point set. Is slowly re-sorted to select suitable starting point to simplify the point selection.
	 */
	std::vector<Shared_Point> points;

	/**
	 * Searchs of starting points and put them in the plane.
	 */
	void selectP1P2(boost::shared_ptr<Plane> pl);

	/**
	 * Resortings the set of points: at the end points are already proven, in the beginning still to be verified.
	 */
	void moveToEnd(size_t i);
	void movePtrToEnd(Shared_Point p);


	//Callbackfunktion für nebenläufiges Zurückgeben der Ebene.
	/**
	 * Callback function for returning the already processed plane.
	 */
	void (*funcPtr)(boost::shared_ptr<Plane>);
	boost::function<void(boost::shared_ptr<Plane>)> callback_;

	bool neighboursCalculated;

	boost::recursive_mutex mutex;

};
}
#endif /* REGIONGROWING_H_ */
