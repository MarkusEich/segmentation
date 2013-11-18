#include <segmentation/spatialObject.h>
#include <segmentation/plane.h>
#include <unistd.h>
#include <iostream>

#ifdef USE_DEBUG
#include <glwidget.h>
#include <point.h>
#endif

using namespace segmentation;

SpatialObject::SpatialObject(){
}

SpatialObject::~SpatialObject(){
}

void SpatialObject::getCoG(double& x, double& y, double& z)const{
	x=centerOfMass[0];
	y=centerOfMass[1];
	z=centerOfMass[2];
}

Eigen::Vector3d SpatialObject::getCoG()const{
	return Eigen::Vector3d(centerOfMass[0], centerOfMass[1], centerOfMass[2]);
}

const Eigen::Vector3d SpatialObject::getInertiaAxisX()const{
	return inertiaAxisX;
}

const Eigen::Vector3d SpatialObject::getInertiaAxisY()const{
	return inertiaAxisY;
}

const Eigen::Vector3d SpatialObject::getNormalVector()const{
	return normalVector;
}

double SpatialObject::onTheObject(SpatialObject &obj){
	double x,y,z, x1,y1,z1;
	this->getCoG(x,y,z);
	obj.getCoG(x1,y1,z1);
	Eigen::Vector3d vertices[4];

	SpatialObject *lower = (z<z1)?this:&obj;

	RotatedRectangle rect = (z<z1)?obj.getMaxExpansion():this->getMaxExpansion();
	RotatedRectangle rect_lower = (z1<z)?obj.getMaxExpansion():this->getMaxExpansion();
	Eigen::Vector3d centerOfMass = (z<z1)?Eigen::Vector3d(x,y,z):Eigen::Vector3d(x1,y1,z1);
	Eigen::Vector3d inertiaAxisX = (z<z1)?this->getInertiaAxisX():obj.getInertiaAxisX();
	Eigen::Vector3d inertiaAxisY = (z<z1)?this->getInertiaAxisY():obj.getInertiaAxisY();
	Eigen::Vector3d n_vector = (z<z1)?this->getNormalVector():obj.getNormalVector();
	for(int i=0; i<4; i++){
		vertices[i]=rect.edges[i];
	}
	//Lotfusspunktverfahren
	Eigen::Vector3d dist1 = (vertices[0]-centerOfMass).array() * (n_vector).array();
	double lenght1= dist1.norm();
	Eigen::Vector3d dist2 = (vertices[1]-centerOfMass).array() * (n_vector).array();
	double lenght2= dist2.norm();
	Eigen::Vector3d dist3 = (vertices[2]-centerOfMass).array() * (n_vector).array();
	double lenght3= dist3.norm();
	Eigen::Vector3d dist4 = (vertices[3]-centerOfMass).array() * (n_vector).array();
	double lenght4= dist4.norm();
	
	//std::cout << "Dist 1: " << dist1.transpose() << std::endl;
	//std::cout << "Dist 2: " << dist2.transpose() << std::endl;
	//std::cout << "Dist 3: " << dist3.transpose() << std::endl;
	//std::cout << "Dist 4: " << dist4.transpose() << std::endl;
	Eigen::Vector3d edge;
	Eigen::Vector3d vertex;
	if(lenght1<lenght2 && lenght1<lenght3 && lenght1<lenght4){
		edge = dist1;
		vertex = vertices[0];
	}else if(lenght2<lenght3 && lenght2<lenght4){
		edge = dist2;
		vertex = vertices[1];
	}else if(lenght3<lenght4){
		edge = dist3;
		vertex = vertices[2];
	}else{
		edge = dist4;
		vertex = vertices[3];
	}
	Eigen::Vector3d p = vertex - edge;
	
	double alpha, betha;
	
	Eigen::AngleAxisd rot(-rect_lower.angle,Eigen::Vector3d::UnitZ());
	Eigen::Vector3d axisX = rot * inertiaAxisX;
	Eigen::Vector3d axisY = rot * inertiaAxisY;

	Eigen::Vector3d vx = (p-centerOfMass).array() * axisX.array();
	Eigen::Vector3d vy = (p-centerOfMass).array() * axisY.array();
//	Eigen::Vector3d vx = (p).cwise() * inertiaAxisX;
//	Eigen::Vector3d vy = (p).cwise() * inertiaAxisY;

	

	alpha = vx.norm();
	betha = vy.norm();
	
	vx.normalize();
	vy.normalize();


	//alpha = -(((centerOfMass[1]-centerOfMass[0]-p[1]+p[0])*inertiaAxisY[2]+(p[2]-centerOfMass[2])*inertiaAxisY[1]+(centerOfMass[2]-p[2])*inertiaAxisY[0])/
	//	((inertiaAxisX[1]-inertiaAxisX[0])*inertiaAxisY[2]-inertiaAxisX[2]*inertiaAxisY[1]+inertiaAxisX[2]*inertiaAxisY[0]));
	//betha = -((alpha*inertiaAxisX[2]+centerOfMass[2]-p[2])/(inertiaAxisY[2]));
	//alpha=fabs(alpha);
	//betha=fabs(betha);
	
	//std::cout << "Point 1: " << vertices[0].transpose() << std::endl;
	//std::cout << "Point 2: " << vertices[1].transpose() << std::endl;
	//std::cout << "Point 3: " << vertices[2].transpose() << std::endl;
	//std::cout << "Point 4: " << vertices[3].transpose() << std::endl;
	//std::cout << "p: " << p.transpose() << std::endl;
	//std::cout << "Vertex: " << vertex.transpose() << std::endl;
	//std::cout << "Edge: " << edge.transpose() << std::endl;
	//std::cout << "Normalvector: " << n_vector.transpose() << std::endl;

	//printf("Lengths: %f,%f,%f,%f\n",lenght1,lenght2,lenght3,lenght4);
	//std::cout << "COG: " << centerOfMass << std::endl;
	//std::cout << "Point: " << p << std::endl;
	//std::cout << "Intertia Axis: " << inertiaAxisX << " Axis 2: " << inertiaAxisY << std::endl;
	if(alpha<(rect_lower.width*1.0) && betha<(rect_lower.height*1.0)){
	

		//Double check for plane, if point is near egnouth	
		Plane *plane = dynamic_cast<Plane*>(lower);
		if(plane != 0){
			//printf("HÄÄÄÄÄÄÄÄÄÄÄÄÄÄHHHHHHHHHHH: %i\n",(long int)plane);
			bool found=false;
			Point ref(p[0],p[1],p[2]);
			for(std::set<Shared_Point>::iterator it = plane->points.begin(); it != plane->points.end(); it++){
				if(ref.distance(*it->get()) < 0.3){
					found=true;
					break;
				}
			}
			if(!found){
				return -1;
			}
		}

#ifdef USE_DEBUG		
		if(edge.norm() < 0.5){
			printf("alpha: %f, Beta: %f sizes: %f %f\n",alpha,betha,rect_lower.width,rect_lower.height);
			printf("Point on plane: %f,%f,%f\n",p[0],p[1],p[2]);
			printf("Analge is: %f\n",rect_lower.angle);
			//GLWidget::debug_points.clear();
			for(double x=0.01;x<betha;x+=0.01){
					Eigen::Vector3d point = (x*vx) +centerOfMass;
					Shared_Point_3d _p(new Point_3d(point[0],point[1],point[2])); 
					GLWidget::debug_points.push_back(_p);
			}
			for(double y=0.01;y<alpha;y+=0.01){
					Eigen::Vector3d point = (y*vy) +centerOfMass;
					Shared_Point_3d _p(new Point_3d(point[0],point[1],point[2])); 
					GLWidget::debug_points.push_back(_p);
			}

		
			Shared_Point_3d _p(new Point_3d(p[0],p[1],p[2])); 
			GLWidget::debug_points.push_back(_p);
		}

		return edge.norm();
#endif
	}
	return -1;
}

bool SpatialObject::above(const SpatialObject &obj){
	return centerOfMass[2]>obj.getCoG()[2];
}

bool SpatialObject::below(const SpatialObject &obj){
	return centerOfMass[2]<obj.getCoG()[2];
}

double SpatialObject::orthogonal(const SpatialObject &obj){
	double orthogonality = this->getNormalVector().dot(obj.getNormalVector());
	return 1.0-orthogonality;
}

double SpatialObject::parallel(const SpatialObject &obj){
	Eigen::Vector3d parallelism = this->getNormalVector().cross(obj.getNormalVector());
	double lenght = parallelism.norm();
	return 1.0-lenght;
}

double SpatialObject::sameHeight(const SpatialObject &obj){
	return (fabs(centerOfMass[2]-obj.getCoG()[2])<0.15)?1.0:0.0;
}
	
double SpatialObject::sameDepth(const SpatialObject &obj){
	return (fabs(centerOfMass[0]-obj.getCoG()[0])<0.10)?1.0:0.0;
}
	
FOP SpatialObject::heighest(std::vector<FOP> &objects){
	std::vector<FOP>::iterator it;
	FOP candidate = *objects.begin();
	for(it=objects.begin(); it!=objects.end(); it++){
		if((*it).obj->above(*candidate.obj.get())){
			candidate = (*it);
		}
	}
	return candidate;
}

FOP SpatialObject::heighest(std::list<FOP> &objects){
	std::list<FOP>::iterator it;
	FOP candidate = *objects.begin();
	for(it=objects.begin(); it!=objects.end(); it++){
		if((*it).obj->above(*candidate.obj.get())){
			candidate = (*it);
		}
	}
	return candidate;
}

OAL SpatialObject::heighest(std::vector<OAL> &objects){
	std::vector<OAL>::iterator it;
	OAL candidate = *objects.begin();
	for(it=objects.begin(); it!=objects.end(); it++){
		if((*it).obj->above(*candidate.obj.get())){
			candidate = (*it);
		}
	}
	return candidate;
}

boost::shared_ptr<SpatialObject> SpatialObject::heighest(std::list<boost::shared_ptr<SpatialObject> > &objects){
	std::list<boost::shared_ptr<SpatialObject> > ::iterator it;
	boost::shared_ptr<SpatialObject> candidate = *objects.begin();
	for(it=objects.begin(); it!=objects.end(); it++){
		if((*it)->above(*candidate.get())){
			candidate = (*it);
		}
	}
	return candidate;
}
	
OAL SpatialObject::lowest(std::vector<OAL> &objects){
	std::vector<OAL>::iterator it;
	OAL candidate = *objects.begin();
	for(it=objects.begin(); it!=objects.end(); it++){
		if((*it).obj->below(*candidate.obj.get())){
			candidate = (*it);
		}
	}
	return candidate;

}

bool SpatialObject::isHeighest(std::list<FOP> &objects){
	std::list<FOP>::iterator it;
	for(it=objects.begin(); it!=objects.end(); it++){
		if(this != (*it).obj.get()){
			if(!this->above(*(*it).obj.get())){
				return false;
			}
		}
	}
	return true;

}

bool SpatialObject::lowest(std::list<FOP> &objects){
	std::list<FOP>::iterator it;
	for(it=objects.begin(); it!=objects.end(); it++){
		if(this != (*it).obj.get()){
			if(!this->below(*(*it).obj.get())){
				return false;
			}
		}
	}
	return true;

}

FOP SpatialObject::lowest(std::vector<FOP> &objects){
	std::vector<FOP>::iterator it;
	FOP candidate = *objects.begin();
	for(it=objects.begin(); it!=objects.end(); it++){
		if((*it).obj->below(*candidate.obj.get())){
			candidate = (*it);
		}
	}
	return candidate;

}

boost::shared_ptr<SpatialObject> SpatialObject::lowest(std::list<boost::shared_ptr<SpatialObject> > &objects){
	std::list<boost::shared_ptr<SpatialObject> > ::iterator it;
	boost::shared_ptr<SpatialObject> candidate = *objects.begin();
	for(it=objects.begin(); it!=objects.end(); it++){
		if((*it)->below(*candidate.get())){
			candidate = (*it);
		}
	}
	return candidate;
}


bool SpatialObject::mostRight(std::list<boost::shared_ptr<SpatialObject> > &objects){
	std::list<boost::shared_ptr<SpatialObject> > ::iterator it;
	for(it=objects.begin(); it!=objects.end(); it++){
		if(this->centerOfMass[1]<((*it)->getCoG()[1])){
			return false;
		}
	}
	return true;
}

bool SpatialObject::mostLeft(std::list<boost::shared_ptr<SpatialObject> > &objects){
	std::list<boost::shared_ptr<SpatialObject> > ::iterator it;
	for(it=objects.begin(); it!=objects.end(); it++){
		if(this->centerOfMass[1]>((*it)->getCoG()[1])){
			return false;
		}
	}
	return true;
}

boost::shared_ptr<SpatialObject> SpatialObject::mostRearward(std::list<boost::shared_ptr<SpatialObject> > &objects){
	std::list<boost::shared_ptr<SpatialObject> > ::iterator it;
	boost::shared_ptr<SpatialObject> candidate = *objects.begin();
	for(it=objects.begin(); it!=objects.end(); it++){
		if(candidate->centerOfMass[0]>((*it)->getCoG()[0])){
			candidate = (*it);
		}
	}
	return candidate;
}

boost::shared_ptr<SpatialObject> SpatialObject::mostForward(std::list<boost::shared_ptr<SpatialObject> > &objects){
	std::list<boost::shared_ptr<SpatialObject> > ::iterator it;
	boost::shared_ptr<SpatialObject> candidate = *objects.begin();
	for(it=objects.begin(); it!=objects.end(); it++){
		if(candidate->centerOfMass[0]<((*it)->getCoG()[0])){
			candidate = (*it);
		}
	}
	return candidate;
}
	
FOP SpatialObject::mostForward(std::vector<FOP> &objects){
	std::vector<FOP>::iterator it;
	FOP candidate = *objects.begin();
	for(it=objects.begin(); it!=objects.end(); it++){
		if(candidate.obj->centerOfMass[0]<((*it).obj->getCoG()[0])){
			candidate = (*it);
		}
	}
	return candidate;

}
FOP SpatialObject::mostForward(std::list<FOP> &objects){
	std::list<FOP>::iterator it;
	FOP candidate = *objects.begin();
	for(it=objects.begin(); it!=objects.end(); it++){
		if(candidate.obj->centerOfMass[0]<((*it).obj->getCoG()[0])){
			candidate = (*it);
		}
	}
	return candidate;

}
OAL SpatialObject::mostForward(std::vector<OAL> &objects){
	std::vector<OAL>::iterator it;
	OAL candidate = *objects.begin();
	for(it=objects.begin(); it!=objects.end(); it++){
		if(candidate.obj->centerOfMass[0]<((*it).obj->getCoG()[0])){
			candidate = (*it);
		}
	}
	return candidate;

}







