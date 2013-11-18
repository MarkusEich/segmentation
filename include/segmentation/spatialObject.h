#ifndef SPATIAL_OBJECT_H
#define SPATIAL_OBJECT_H

#include <boost/shared_ptr.hpp>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <segmentation/rotatedRectangle.h>
#include <iostream>
#include <math.h>
#include <fstream>
#include <stdlib.h>
#include <list>
#include <vector>
#include <string.h>


namespace segmentation{

struct OAL;
struct FOP;


class SpatialObject{

public:
	SpatialObject();
	~SpatialObject();

	void getCoG(double& x, double& y, double& z)const;
	Eigen::Vector3d getCoG()const;
	virtual double getSizeOfPlane()=0;
	virtual double getRectangularness()=0;
	virtual RotatedRectangle getMaxExpansion()=0;
	const Eigen::Vector3d getInertiaAxisX()const;
	const Eigen::Vector3d getInertiaAxisY()const;
	const Eigen::Vector3d getNormalVector()const;

	double onTheObject(SpatialObject &obj);
	bool above(const SpatialObject &obj);
	bool below(const SpatialObject &obj);
	double sameHeight(const SpatialObject &obj);
	static FOP heighest(std::vector<FOP> &objects);
	static OAL heighest(std::vector<OAL> &objects);
	boost::shared_ptr<SpatialObject> heighest(std::list<boost::shared_ptr<SpatialObject> > &objects);
	static OAL lowest(std::vector<OAL> &objects);
	boost::shared_ptr<SpatialObject>lowest(std::list<boost::shared_ptr<SpatialObject> > &objects);
	static FOP lowest(std::vector<FOP> &objects);
	bool lowest(std::list<FOP> &objects);
	static FOP heighest(std::list<FOP> &objects);
	bool isHeighest(std::list<FOP> &objects);
	double orthogonal(const SpatialObject &obj);
	double parallel(const SpatialObject &obj);
	double sameDepth(const SpatialObject &obj);
	bool mostRight(std::list<boost::shared_ptr<SpatialObject> > &objects);
	bool mostLeft(std::list<boost::shared_ptr<SpatialObject> > &objects);
	boost::shared_ptr<SpatialObject> mostRearward(std::list<boost::shared_ptr<SpatialObject> > &objects);
	static FOP mostForward(std::list<FOP> &objects);
	static FOP mostForward(std::vector<FOP> &objects);
	static OAL mostForward(std::vector<OAL> &objects);
	boost::shared_ptr<SpatialObject> mostForward(std::list<boost::shared_ptr<SpatialObject> > &objects);

protected:
	double	centerOfMass[3];
	Eigen::Vector3d inertiaAxisX;
	Eigen::Vector3d inertiaAxisY;
    Eigen::Vector3d normalVector;

};


/**
 * Struct for type and probablility of a object(Type-And-Probability).
 */
struct TAP{
	std::string type;
	double probability;
};

/**
 * Struct for representation of the feature with type and probability(Feature-Object-Probablilty).
 */
struct FOP{
	boost::shared_ptr<SpatialObject> obj;
	boost::shared_ptr<TAP> tap;
};
/**
 * Struct for representation of the object with the list of types and probabilities(Object-And-List).
 */
struct OAL{
	boost::shared_ptr<SpatialObject> obj;
	std::vector<boost::shared_ptr<TAP> > tap;
};
}
#endif

