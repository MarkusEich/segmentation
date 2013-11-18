/*****************************************
 ** Markus Eich (2013)
 ** Malgorzata Goldhoorn (2010)
 ** DFKI GmbH
 **
 **
*/


#ifndef PLANE_H_
#define PLANE_H_

#include <segmentation/spatialObject.h>
#include <segmentation/rotatedRectangle.h>
#include <segmentation/point.h>
#include <boost/shared_ptr.hpp>
#include <Eigen/Core>
#include <Eigen/Geometry>

#include <queue>
#include <set>



namespace cv{
class Mat;
}

namespace segmentation{

class Plane : public SpatialObject{


public:
    Plane();
    Plane(Eigen::Vector3d pos, Eigen::Vector3d rotation, Eigen::Vector2d size, double alpha,double intervall=0.03, bool minimal=false);
    Plane(Eigen::Vector3d cog, Eigen::Vector3d ineriaAxisX, Eigen::Vector3d ineriaAxisY,Eigen::Vector3d ineriaAxisZ, Eigen::Vector3d S, Eigen::Matrix<double,9,1> C, Eigen::Matrix<double,9,1> productMatrix, double bias,std::set<Shared_Point> points,double alpha);
	~Plane();

	/**
	 * Calculates a 2D hull by using alpha shapes and returns a list of pairs.
	 */
	std::list<std::list<std::pair<Shared_Point, Shared_Point> > > getHull();
	
	/**
	 * Draws the alpha shape with the help of edges and return the filled area as a bitmap.
	 */	
    //Replaced QImage
    boost::shared_ptr<cv::Mat> getHullCVImages(double &pixSize);

    //Get hull only
    boost::shared_ptr<cv::Mat> getHullCVImagesNofill(double &pixSize);

	/**
	 * Returns the translated points of the plane.
	 */
	std::list<std::list<std::pair<Shared_Point,Shared_Point> > > getTransPoints();

	/**
	 * Calculates the size of the plane.
	 */
	double getSizeOfPlane();

    //set alpha value for hull detection
    void setAlpha(double alpha){alpha_=alpha;}

	/**
	 * Points in the plane, should never be modified directly,
	 * but can be used to read this attribute
	 * change only with addPoints
	 */
	std::set<Shared_Point> points;

	/**
	 * Returns the euler angles of the plane in  space.
	 */
	double getRoll();
	double getPitch();
	double getYaw();
	
	Eigen::Quaterniond getQuaternion();
	
	Eigen::Matrix<double,3,3> getRotationMatrix();


	/**
	 * Returns rotated square.
	 */
	RotatedRectangle getMaxExpansion();
	
	/**
	 * Returns the position of plane.
	 */
	Eigen::Vector3d getPosition();

	/**
	 * Returns the angle of plane as an quaternion.
	 */
	Eigen::Quaterniond getPlaneAngle();

	/**
	 * Returns percentage squareness.
	 */
	double getRectangularness();

	/**
	 * Adds a point. 
	 */
	void addPoint(Shared_Point p);

	/**
	 * Calculates the mean square error from internal state.
	 */
	double mse();

	/**
	 * Removes point from the plane.
	 */
	void removePoint(Shared_Point p);

	/**
	 * Checks if there are neighbors known.
	 */
	bool knowNeighbours();

	/**
	 * Return the number of known neighbors.
	 */
	int neighbourCount();

	/**
	 * Merges the given point to the list of neighbours.
	 */
	void mergeNearestNeighbours(Shared_Point p);

	/**
	 * Returns nearest neighbor. 
	 */
	Shared_Point findNearestNeighbour();

	/**
	 * Save internal state.
	 */
	void saveState();

	/**
	 * Restore state.
	 */
	void restoreState();

	/**
	 * Distance of the point of the hypothetical plane.
	 */
	double planeDist(Shared_Point p);

    /**
     * get the expansion of the plane in world cooridinates
     */
    Eigen::Vector3d getSpatialExpansion();



private:

	//alpha value of the alpha shapes
	double alpha_;

	/**
	 * Inertia axis.
	 */
//	Eigen::Vector3d inertiaAxisX;
//	Eigen::Vector3d inertiaAxisY;


	/**
	 * Calculates normal vector.
	 */
//	Eigen::Vector3d	normalVector;

	/**
	 * Center of mass of the points in the plane.
	 */
//	double	centerOfMass[3];
	
	RotatedRectangle rectangle;

	bool rectangleCalculated;
	
	/**
	 * Recalculates the normal vector of internal state. 
	 */
	void updateNormal();

	double calculateAngle(Point &p1, Point &p2, Point &p3);

	int matrixIndex(int row, int col) const;

	/**
	 * List of known nearst neighbours of points in the plane.
	 */
	std::priority_queue<Point::Neighbour> neighbour_queue;
	
	/**
	 * Backup to restore the list.
	 */
	std::priority_queue<Point::Neighbour> neighbour_queue_backup;


	/**
	 * Sum of the coordinates of all points.
	 */
	double	S[3];

	/**
	 * Covariance matrix of all points.
	 */
	double	C[9];

	/**
	 * Distance of the plane from the zero point in the direction of the normal vector.
	 */
	double	bias;

	/**
	 * Will be required for incremental computing the MSE.
	 */
	double	productMatrix[9];

	/**
	 * Saved internal state, is required for single-step "backtracking". 
	 */
	double savedS[3];
	double savedC[9];
	Eigen::Vector3d savedNormalVector;
	double savedBias;
	double savedProductMatrix[9];
	double savedCenterOfMass[3];

	/**
	 * Saves translated points.
	 */
	std::list<std::list<std::pair<Shared_Point,Shared_Point> > > transPoints;

	friend std::istream& operator >>(std::istream &is, Plane &obj);
	friend std::ostream& operator <<(std::ostream &os, const Plane &obj);
};
typedef boost::shared_ptr<Plane> Shared_Plane;
}
#endif /* PLANE_H_ */
