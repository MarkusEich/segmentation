#ifndef ROTATEDRECTANGLE_H_
#define ROTATEDRECTANGLE_H_

#include <Eigen//Core>
#include <Eigen/Geometry>

namespace segmentation{

class RotatedRectangle{

public:
	RotatedRectangle(double width=0, double height=0, double angle=0, double pos_x=0, double pos_y=0):
	width(width), height(height), angle(angle), pos_x(pos_x), pos_y(pos_y){
	}

	~RotatedRectangle(){
	}

	double width;
	double height;
	double angle;
	double pos_x;
	double pos_y;
	Eigen::Vector3d edges[4];
};

}
#endif





