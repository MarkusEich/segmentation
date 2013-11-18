/*****************************************
 ** Markus Eich (2013)
 ** Malgorzata Goldhoorn (2010)
 ** DFKI GmbH
 **
 **
*/


#include <segmentation/regionGrowing.h>
#include <segmentation/plane.h>
#include <segmentation/hull.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <vector>
#include <unistd.h>
#include <opencv2/opencv.hpp>
#include <boost/foreach.hpp>


using namespace segmentation;

/**
 * Creates plane and initialize internal state.
 */
Plane::Plane() {
	bias = 0.0;
	savedBias = 0.0;
	normalVector.setZero();
    alpha_=0.01;
	for(int i = 0; i < 3; i++) {
		S[i] = 0.0;
		centerOfMass[i] = 0.0;
		savedS[i] = 0.0;
		savedNormalVector[i] = 0.0;
		savedCenterOfMass[i] = 0.0;
	}
	for(int i = 0; i < 9; i++) {
		C[i] = 0.0;
		productMatrix[i] = 0.0;
		savedC[i] = 0.0;
		savedProductMatrix[i] = 0.0;
	}
	rectangleCalculated = false;
}

Plane::~Plane() {
	points.clear();
}

Plane::Plane(Eigen::Vector3d cog, Eigen::Vector3d ineriaAxisX, Eigen::Vector3d ineriaAxisY,Eigen::Vector3d ineriaAxisZ, Eigen::Vector3d S, Eigen::Matrix<double,9,1> C, Eigen::Matrix<double,9,1> productMatrix, double bias,std::set<Shared_Point> points,double alpha):
	points(points)	

{
	for(int i=0;i<3;i++){
	centerOfMass[i] = cog[i];
	this->inertiaAxisX[i] = ineriaAxisX[i];
	this->inertiaAxisY[i]=ineriaAxisY[i];
	this->normalVector[i]=ineriaAxisZ[i];
	this->S[i]=S[i];
	}
	this->bias=bias;
	for(int i=0;i<9;i++){
	this->C[i]=C(i,0);
	this->productMatrix[i]=productMatrix(i,0);
	}
	
	rectangleCalculated = false;

    alpha_=alpha;
}


Plane::Plane(Eigen::Vector3d pos, Eigen::Vector3d rotation, Eigen::Vector2d size, double alpha, double intervall, bool minimal){
	Eigen::Quaterniond q(1,0,0,0);
	Eigen::AngleAxisd rot_x(rotation[0],Eigen::Vector3d::UnitX());
	Eigen::AngleAxisd rot_y(rotation[1],Eigen::Vector3d::UnitY());
	Eigen::AngleAxisd rot_z(rotation[2],Eigen::Vector3d::UnitZ());
	q =  q* rot_x * rot_y * rot_z;
	inertiaAxisX = q.toRotationMatrix().col(0);
	inertiaAxisY = q.toRotationMatrix().col(1);
	normalVector = q.toRotationMatrix().col(2);
	
	Eigen::Vector3d S;
	Eigen::Matrix<double,9,1> C;
	Eigen::Matrix<double,9,1> productMatrix;
	S.setZero();
	C.setZero();
	productMatrix.setZero();

	bias = 0;	

	if(!minimal){
		for(double x = -size[0]/2.0; x < +size[0]/2.0; x+=intervall){
			for(double y = -size[1]/2.0; y <+size[1]/2.0; y+=intervall){
				Eigen::Vector3d p(x,y,0);
				p = (q*p)+pos;
				points.insert(Shared_Point(new Point(p[0],p[1],p[2])));
			}	
		}
	}else{
		Eigen::Vector3d p0(-size[0]/2.0,-size[1]/2.0,0);
		Eigen::Vector3d p1(-size[0]/2.0,size[1]/2.0,0);
		Eigen::Vector3d p2(size[0]/2.0,-size[1]/2.0,0);
		Eigen::Vector3d p3(size[0]/2.0,size[1]/2.0,0);
		p0 = (q*p0)+pos;
		p1 = (q*p1)+pos;
		p2 = (q*p2)+pos;
		p3 = (q*p3)+pos;
		points.insert(Shared_Point(new Point(p0[0],p0[1],p0[2])));
		points.insert(Shared_Point(new Point(p1[0],p1[1],p1[2])));
		points.insert(Shared_Point(new Point(p2[0],p2[1],p2[2])));
		points.insert(Shared_Point(new Point(p3[0],p3[1],p3[2])));
	}
	 
	for(int i=0;i<3;i++){
		centerOfMass[i] = pos[i];
		this->S[i]=S[i];
	}
	for(int i=0;i<9;i++){
		this->C[i]=C(i,0);
		this->productMatrix[i]=productMatrix(i,0);
	}
	rectangleCalculated = false;

    alpha_=alpha;
}


int Plane::matrixIndex(int row, int col) const {
	return row*3+col;
}

/**
 * Adds point to the plane.
 */
void Plane::addPoint(Shared_Point p) {
	points.insert(p);

	//part of the state is required in the old version for the update formulas, so cached
	double oldSum[3];
	double oldCoM[3];
	memcpy(oldSum, S, sizeof(S));
	memcpy(oldCoM, centerOfMass, sizeof(centerOfMass));

	//updating sum vector and center of mass
	for( int i = 0; i < 3; i++ ) {
		S[i] += p->p[i];
		//CoM ist the summ of all points coordinates divied point size
		centerOfMass[i] = S[i] / points.size();
	}

	//updating the covariance- and product matrix
	for(int r = 0; r < 3; r++)                                      //1 = x, 2 = y, 3 = z
		for(int c = 0; c < 3; c++) {
			C[matrixIndex(r,c)] += p->p[r]*p->p[c]		//r_(n+1) * c_(n+1)
			                       - centerOfMass[r]*S[c]	//-m_r(n+1) * S_c(n+1)
			                       + oldCoM[r]*oldSum[c]; 	//+ m_r(n)*S_c(n)

			productMatrix[matrixIndex(r,c)] += p->p[r]*p->p[c];
		}

	//calculate normal vector
	updateNormal();
}

/**
 * MSE very fast (constant time) calculated from collected incrementally data.
 */
double Plane::mse() {
	//with less than 4 points the calculation is irrelevant, because three points are always in a plane (Thus MES = 0)
	//less than 3 points, there is no plane
	if(points.size() < 4)
		return 0.0;

	//see paper for formula
	double m = 0.0;
	double n = 0.0;
	for( int i = 0; i < 3; i++) {
		n += normalVector[i]*S[i];
		for( int j = 0; j < 3; j++) {
			m += normalVector[i]*normalVector[j]*productMatrix[matrixIndex(i,j)];
		}
	}
	return fabs((m + 2.0*bias*n)/points.size() + bias*bias);
}

/**
 * Removes point of the plane. 
 */
void Plane::removePoint(Shared_Point p) {
	points.erase(p);
}

/**
 * Checks if there are more neighbors known.
 */
bool Plane::knowNeighbours() {
	return !neighbour_queue.empty();
}

/**
 * Counts number of neighbors.
 */
int Plane::neighbourCount(){
	return neighbour_queue.size();
}

/**
 * Merges the given point to the list of neighbours.
 */
void Plane::mergeNearestNeighbours(Shared_Point p) {	
	std::list<Point::Neighbour>::iterator it;
	for(it = p->nearestNeighbours.begin(); it != p->nearestNeighbours.end(); it++){
		neighbour_queue.push(*it);
	}
}

/**
 * Find next neighbour, who still not belongs to the plane.
 */
Shared_Point Plane::findNearestNeighbour() {
	while(neighbour_queue.size()) {
		Shared_Point p = neighbour_queue.top().point.lock();
		neighbour_queue.pop();
		if(points.find(p) == points.end())
			return p;
	}
	return Shared_Point();
}

/**
 * Calculates the normal vector.
 */
void Plane::updateNormal() {
	//less than 3 points are not unique normal vector, because they do not form a plane
	if(points.size() < 3)
		return;

	//calculated eigen values with the GNU Scientific Library, somewhat complicated, because complex eigen values may occur (will not happen in practice)
	//therefore, after processing, we have only real values 
	double CBackup[9];
	for(int i=0;i<9;i++)
		CBackup[i] = C[i];

	 Eigen::Vector3d eigenVector;
	 gsl_matrix_view m = gsl_matrix_view_array (CBackup, 3, 3);
	 gsl_vector_complex *eval = gsl_vector_complex_alloc (3);
	 gsl_matrix_complex *evec = gsl_matrix_complex_alloc (3, 3);
	 gsl_eigen_nonsymmv_workspace * w = gsl_eigen_nonsymmv_alloc (3);
	 gsl_eigen_nonsymmv (&m.matrix, eval, evec, w);
	 gsl_eigen_nonsymmv_free (w);

	 gsl_complex eval_0 = gsl_vector_complex_get (eval, 0);
	 gsl_complex eval_1 = gsl_vector_complex_get (eval, 1);
	 gsl_complex eval_2 = gsl_vector_complex_get (eval, 2);

	 gsl_vector_complex_view evec_0 = gsl_matrix_complex_column (evec, 0);
	 gsl_vector_complex_view evec_1 = gsl_matrix_complex_column (evec, 1);
	 gsl_vector_complex_view evec_2 = gsl_matrix_complex_column (evec, 2);
	//search for the eigen vector, whose eigen value is the smallest
	 if(GSL_IMAG(eval_0) == 0 && GSL_IMAG(eval_1) == 0 && GSL_IMAG(eval_2) == 0 && GSL_REAL(eval_0) <= GSL_REAL(eval_1) && GSL_REAL(eval_0) <= GSL_REAL(eval_2)){
		//eigenValue = eval_0;
		 eigenVector[0] = GSL_REAL(gsl_vector_complex_get(&evec_0.vector,0));
		 eigenVector[1] = GSL_REAL(gsl_vector_complex_get(&evec_0.vector,1));
		 eigenVector[2] = GSL_REAL(gsl_vector_complex_get(&evec_0.vector,2));
		if(GSL_REAL(eval_1) <= GSL_REAL(eval_2)){
			inertiaAxisX[0] = GSL_REAL(gsl_vector_complex_get(&evec_1.vector,0));
			inertiaAxisX[1] = GSL_REAL(gsl_vector_complex_get(&evec_1.vector,1));
			inertiaAxisX[2] = GSL_REAL(gsl_vector_complex_get(&evec_1.vector,2));

			inertiaAxisY[0] = GSL_REAL(gsl_vector_complex_get(&evec_2.vector,0));
			inertiaAxisY[1] = GSL_REAL(gsl_vector_complex_get(&evec_2.vector,1));
			inertiaAxisY[2] = GSL_REAL(gsl_vector_complex_get(&evec_2.vector,2));
		}else{
			inertiaAxisX[0] = GSL_REAL(gsl_vector_complex_get(&evec_2.vector,0));
			inertiaAxisX[1] = GSL_REAL(gsl_vector_complex_get(&evec_2.vector,1));
			inertiaAxisX[2] = GSL_REAL(gsl_vector_complex_get(&evec_2.vector,2));

			inertiaAxisY[0] = GSL_REAL(gsl_vector_complex_get(&evec_1.vector,0));
			inertiaAxisY[1] = GSL_REAL(gsl_vector_complex_get(&evec_1.vector,1));
			inertiaAxisY[2] = GSL_REAL(gsl_vector_complex_get(&evec_1.vector,2));
		}
	 }else if(GSL_IMAG(eval_1) == 0 && GSL_IMAG(eval_0) == 0 && GSL_IMAG(eval_2) == 0 && GSL_REAL(eval_1) <= GSL_REAL(eval_2)){
		//eigenValue = eval_1;
		 eigenVector[0] = GSL_REAL(gsl_vector_complex_get(&evec_1.vector,0));
		 eigenVector[1] = GSL_REAL(gsl_vector_complex_get(&evec_1.vector,1));
		 eigenVector[2] = GSL_REAL(gsl_vector_complex_get(&evec_1.vector,2));

		if(GSL_REAL(eval_0) <= GSL_REAL(eval_2)){
			inertiaAxisX[0] = GSL_REAL(gsl_vector_complex_get(&evec_0.vector,0));
			inertiaAxisX[1] = GSL_REAL(gsl_vector_complex_get(&evec_0.vector,1));
			inertiaAxisX[2] = GSL_REAL(gsl_vector_complex_get(&evec_0.vector,2));

			inertiaAxisY[0] = GSL_REAL(gsl_vector_complex_get(&evec_2.vector,0));
			inertiaAxisY[1] = GSL_REAL(gsl_vector_complex_get(&evec_2.vector,1));
			inertiaAxisY[2] = GSL_REAL(gsl_vector_complex_get(&evec_2.vector,2));
		}else{
			inertiaAxisX[0] = GSL_REAL(gsl_vector_complex_get(&evec_2.vector,0));
			inertiaAxisX[1] = GSL_REAL(gsl_vector_complex_get(&evec_2.vector,1));
			inertiaAxisX[2] = GSL_REAL(gsl_vector_complex_get(&evec_2.vector,2));

			inertiaAxisY[0] = GSL_REAL(gsl_vector_complex_get(&evec_0.vector,0));
			inertiaAxisY[1] = GSL_REAL(gsl_vector_complex_get(&evec_0.vector,1));
			inertiaAxisY[2] = GSL_REAL(gsl_vector_complex_get(&evec_0.vector,2));
		}
	 }else if(GSL_IMAG(eval_2) == 0 && GSL_IMAG(eval_0) == 0 && GSL_IMAG(eval_1) == 0) {
		 //eigenValue = eval_2;
		 eigenVector[0] = GSL_REAL(gsl_vector_complex_get(&evec_2.vector,0));
		 eigenVector[1] = GSL_REAL(gsl_vector_complex_get(&evec_2.vector,1));
		 eigenVector[2] = GSL_REAL(gsl_vector_complex_get(&evec_2.vector,2));;
		if(GSL_REAL(eval_0) <= GSL_REAL(eval_1)){
			inertiaAxisX[0] = GSL_REAL(gsl_vector_complex_get(&evec_0.vector,0));
			inertiaAxisX[1] = GSL_REAL(gsl_vector_complex_get(&evec_0.vector,1));
			inertiaAxisX[2] = GSL_REAL(gsl_vector_complex_get(&evec_0.vector,2));

			inertiaAxisY[0] = GSL_REAL(gsl_vector_complex_get(&evec_1.vector,0));
			inertiaAxisY[1] = GSL_REAL(gsl_vector_complex_get(&evec_1.vector,1));
			inertiaAxisY[2] = GSL_REAL(gsl_vector_complex_get(&evec_1.vector,2));
		}else{
			inertiaAxisX[0] = GSL_REAL(gsl_vector_complex_get(&evec_1.vector,0));
			inertiaAxisX[1] = GSL_REAL(gsl_vector_complex_get(&evec_1.vector,1));
			inertiaAxisX[2] = GSL_REAL(gsl_vector_complex_get(&evec_1.vector,2));

			inertiaAxisY[0] = GSL_REAL(gsl_vector_complex_get(&evec_0.vector,0));
			inertiaAxisY[1] = GSL_REAL(gsl_vector_complex_get(&evec_0.vector,1));
			inertiaAxisY[2] = GSL_REAL(gsl_vector_complex_get(&evec_0.vector,2));
		}
	 } else {
	 	inertiaAxisX.setIdentity();
		inertiaAxisY.setIdentity();
	 	eigenVector[0] = 1.0;
	 	eigenVector[1] = 0.0;
	 	eigenVector[2] = 0.0;
	//std::cout << "no real eigenvalues found" << std::endl;
	 }
	 gsl_vector_complex_free (eval);
	 gsl_matrix_complex_free (evec);

	 //calculation of normal vector
	 //the eigen vector with the smallest eigen value corresponds to the normal vector of the plane, 
	 //these normalized for comfortable distance calculation
	 normalVector = eigenVector.normalized();

	 //calculation of bias
	 bias = -(normalVector[0]*centerOfMass[0] + normalVector[1]*centerOfMass[1] + normalVector[2]*centerOfMass[2]);
}

/**
 * Distance point-by-plane projection method.
 */
double Plane::planeDist(Shared_Point p) {
	double diffVector[3];
	diffVector[0]=(p->x-centerOfMass[0]);
	diffVector[1]=(p->y-centerOfMass[1]);
	diffVector[2]=(p->z-centerOfMass[2]);

	double dist = fabs(diffVector[0]*normalVector[0]+
			   diffVector[1]*normalVector[1]+
			   diffVector[2]*normalVector[2]);

	return dist;
}

/**
 * Stores the internal state.
 */
void Plane::saveState() {
	memcpy(savedS, S, sizeof(S));
	memcpy(savedC, C, sizeof(C));
	memcpy(savedProductMatrix, productMatrix, sizeof(productMatrix));
	memcpy(savedCenterOfMass, centerOfMass, sizeof(centerOfMass));
	savedNormalVector = normalVector;
	savedBias = bias;
	neighbour_queue_backup =  neighbour_queue;
}

/**
 * Restores the stored internal state.
 */
void Plane::restoreState() {
	memcpy(S, savedS, sizeof(S));
	memcpy(C, savedC, sizeof(C));
	normalVector = savedNormalVector;
	memcpy(centerOfMass, savedCenterOfMass, sizeof(centerOfMass));
	memcpy(productMatrix, savedProductMatrix, sizeof(productMatrix));
	bias = savedBias;
	neighbour_queue =  neighbour_queue_backup;
}

/**
 * Returns the translated elements of the plane.
 */
std::list<std::list<std::pair<Shared_Point,Shared_Point> > > Plane::getTransPoints(){
	if(this->transPoints.empty()){
               //printf("Calculating hull of %i points\n",(int)points.size());
		std::list<Point> trans_points;
		double angle, angle2;
		Eigen::Vector3d rotationAxis;
		Eigen::Vector3d centerOfMass(this->centerOfMass[0], this->centerOfMass[1], this->centerOfMass[2]);
		rotationAxis = normalVector.cross(Eigen::Vector3d::UnitY());
		rotationAxis.normalize();
		angle = acos(normalVector[1]);
		Eigen::AngleAxisd m = Eigen::AngleAxisd(angle, rotationAxis);
		Eigen::Vector3d inertiaAxisXNew = m*inertiaAxisX;
		inertiaAxisXNew.normalize();
		angle2 = acos(inertiaAxisXNew.dot(Eigen::Vector3d::UnitX()));	
		Eigen::AngleAxisd m2; 
		if(angle2 != 0.0){ //It crashed if the angle is 0.
			Eigen::Vector3d rotationAxis2 = inertiaAxisXNew.cross(Eigen::Vector3d::UnitX());
			rotationAxis2.normalize(); 
			m2 = Eigen::AngleAxisd(angle2, rotationAxis2); 
		}
		
		for(std::set<Shared_Point>::iterator it1=points.begin(); it1!=points.end(); it1++){
			Eigen::Vector3d p((*it1)->p[0], (*it1)->p[1], (*it1)->p[2]);
			p = p-centerOfMass;	 
			p = m * p;
			p = m2 * p;
			p[1] = 0;
			trans_points.push_back(Point(p[0], p[1], p[2]));
		}


#ifdef EXPORT_SINGLE_PLANES
		float cr = float(rand()) / RAND_MAX;
		float cg = float(rand()) / RAND_MAX;
		float cb = float(rand()) / RAND_MAX;

		std::list<Shared_Point_3d> coloured_points;
		for(std::list<Point>::iterator xit = trans_points.begin(); xit != trans_points.end(); xit++) {
			Shared_Point_3d p( new Pointcloud::Point_3d((*xit).p[0], (*xit).p[1], (*xit).p[2], TRANSPARENCY, 255*cr, 255*cg, 255*cb));
			coloured_points.push_back(p);
			
		}
		Ply ply;
		char fileName[500];
		sprintf(fileName,"Plane_Points-%i.ply",coloured_points.size());
		ply.exportPly(coloured_points, fileName);
#endif

        transPoints = Hull::sortedPair(trans_points,alpha_);

	}else{
		//printf("Returning already computed plane-hull\n");
	}

	std::list<std::list<std::pair<Shared_Point,Shared_Point> > > erg;
	for(std::list<std::list<std::pair<Shared_Point,Shared_Point> > >::iterator it=transPoints.begin(); it != transPoints.end(); it++){
		std::list<std::pair<Shared_Point,Shared_Point> > erg1;
		for(std::list<std::pair<Shared_Point,Shared_Point> >::iterator it1 = (*it).begin(); it1 != (*it).end(); it1++){
			Shared_Point p1(new Point((*it1).first->p[0],(*it1).first->p[1],(*it1).first->p[2]));
			Shared_Point p2(new Point((*it1).second->p[0],(*it1).second->p[1],(*it1).second->p[2]));
			erg1.push_back(std::pair<Shared_Point,Shared_Point>(p1,p2));
		}
		erg.push_back(erg1);
	}
	return erg; 
}


/**
 * Calculates a 2D hull by using alpha shapes and returns a list of pairs.
 */
std::list<std::list<std::pair<Shared_Point,Shared_Point> > > Plane::getHull(){
	std::list<std::list<std::pair<Shared_Point,Shared_Point> > > edgesHull;
	std::set<Shared_Point>::iterator it1;
	double angle, angle2;
	Eigen::Vector3d rotationAxis;
	Eigen::Vector3d centerOfMass(this->centerOfMass[0], this->centerOfMass[1], this->centerOfMass[2]);
	rotationAxis = normalVector.cross(Eigen::Vector3d::UnitY());
	rotationAxis.normalize();
	angle = acos(normalVector[1]);
	Eigen::AngleAxisd m = Eigen::AngleAxisd(angle, rotationAxis);

	Eigen::Vector3d inertiaAxisXNew = m*inertiaAxisX;
	inertiaAxisXNew.normalize();
	angle2 = acos(inertiaAxisXNew.dot(Eigen::Vector3d::UnitX()));	
	Eigen::AngleAxisd m2; 
	if(angle2 != 0.0){ //It crashed if the angle is 0.
		Eigen::Vector3d rotationAxis2 = inertiaAxisXNew.cross(Eigen::Vector3d::UnitX());
		rotationAxis2.normalize(); 
		m2 = Eigen::AngleAxisd(angle2, rotationAxis2); 
	}
	edgesHull = getTransPoints(); 
	if(edgesHull.size()==0)
		return edgesHull;
	
	Shared_Point p = edgesHull.begin()->begin()->first;
	Shared_Point p2= edgesHull.begin()->begin()->first;
	std::list< std::list<std::pair<Shared_Point, Shared_Point> > >::iterator it;
	std::list< std::pair<Shared_Point, Shared_Point> >::iterator it_list;
	for(it=edgesHull.begin(); it!=edgesHull.end(); it++){
		for(it_list=(*it).begin(); it_list!=(*it).end(); it_list++){
			if(fabs((*it_list).first->x)>fabs(p->x))
				p=(*it_list).first;
			if(fabs((*it_list).first->z)>fabs(p2->z))
				p2=(*it_list).first;

			if(fabs((*it_list).second->x)>fabs(p->x))
				p=(*it_list).second;
			if(fabs((*it_list).second->z)>fabs(p2->z))
				p2=(*it_list).second;
		}
	}

	std::list< std::list< std::pair<Shared_Point, Shared_Point> > >::iterator it2;
	std::list< std::pair<Shared_Point, Shared_Point> >::iterator it3;
	for(it2=edgesHull.begin(); it2!=edgesHull.end(); it2++){
		for(it3=(*it2).begin(); it3!=(*it2).end(); it3++){
			Eigen::Vector3d p((*it3).first->p[0], (*it3).first->p[1], (*it3).first->p[2]);
			Eigen::Vector3d p1((*it3).second->p[0], (*it3).second->p[1], (*it3).second->p[2]);
			p = m2.inverse() * p;
			p = m.inverse() * p;
			p = p + centerOfMass;
			(*it3).first->p[0]=p[0];
			(*it3).first->p[1]=p[1];
			(*it3).first->p[2]=p[2];
			
			p1 = m2.inverse() * p1;
			p1 = m.inverse() * p1;
			p1 = p1 + centerOfMass;
			(*it3).second->p[0]=p1[0];
			(*it3).second->p[1]=p1[1];
			(*it3).second->p[2]=p1[2];
		}
	}
	return edgesHull;
}

Eigen::Matrix<double,3,3> Plane::getRotationMatrix(){
	Eigen::Matrix<double,3,3> rotMat;
	rotMat.col(0) = inertiaAxisX.normalized();
	rotMat.col(1) = inertiaAxisY.normalized();
	rotMat.col(2) = Eigen::Vector3d(normalVector[0], normalVector[1], normalVector[2]).normalized();
	return rotMat;
}

Eigen::Quaterniond Plane::getQuaternion(){
	Eigen::Matrix<double,3,3> rotMat = getRotationMatrix();
	return Eigen::Quaterniond(rotMat);
}

double Plane::getRoll(){
	Eigen::Matrix<double,3,3> rotMat = getRotationMatrix();
	Eigen::Vector3d vector = rotMat.eulerAngles(0,1,2);
	return vector[1];
}

double Plane::getYaw(){
	Eigen::Matrix<double,3,3> rotMat = getRotationMatrix();
	Eigen::Vector3d vector = rotMat.eulerAngles(0,1,2);
	if(vector[2]<-(M_PI*0.5)){
		vector[2]+=M_PI;
	}else if(vector[2]>(M_PI*0.5)){
		vector[2]-=M_PI;
	}
	 return vector[2];
}

double Plane::getPitch(){
	 Eigen::Matrix<double,3,3> rotMat = getRotationMatrix();
	 Eigen::Vector3d vector = rotMat.eulerAngles(0,1,2);
	 return vector[0];
}

/**
 * Returns the bounding box, which includes the plane.
 */
RotatedRectangle Plane::getMaxExpansion(){
	if(!rectangleCalculated){
	std::list<std::list<std::pair<Shared_Point,Shared_Point> > > edgesHull = getTransPoints();
	std::vector<cv::Point2f> _points;
	std::list< std::list<std::pair<Shared_Point, Shared_Point> > >::iterator it;
	std::list< std::pair<Shared_Point, Shared_Point> >::iterator it_list;
	for(it=edgesHull.begin(); it!=edgesHull.end(); it++){
		for(it_list=(*it).begin(); it_list!=(*it).end(); it_list++){
			cv::Point2f p((*it_list).first->x, (*it_list).first->z);
			_points.push_back(p);
			//printf("Point: %f,%f\n",p.x,p.y);
		} 
	}
	if(_points.size()==0){
		fprintf(stderr,"Point size is zero so not calculating maximum expansion, edgesHull Size is: %i\n", (int)edgesHull.size());
		return RotatedRectangle();
	}
    
 
	cv::RotatedRect rect = cv::minAreaRect(cv::Mat(_points));

	double width = (rect.size.width>rect.size.height)?rect.size.width:rect.size.height;
	double height = (rect.size.width<rect.size.height)?rect.size.width:rect.size.height;
	double angle = rect.angle;

	//Berchnung von Eckpunkten mit Hile von Trygonometrie und dann die Punkte in die Welt rotieren
	Eigen::Vector3d centerOfMass(this->centerOfMass[0], this->centerOfMass[1], this->centerOfMass[2]);
	Eigen::Vector3d rotationAxis = normalVector.cross(Eigen::Vector3d::UnitY());
	rotationAxis.normalize();
	Eigen::AngleAxisd m = Eigen::AngleAxisd(acos(normalVector[1]), rotationAxis);
	Eigen::Vector3d inertiaAxisXNew = m*inertiaAxisX;
	inertiaAxisXNew.normalize();
	double angle2 = acos(inertiaAxisXNew.dot(Eigen::Vector3d::UnitX()));	
	Eigen::AngleAxisd m2; 
	if(angle2 != 0.0){ 
		Eigen::Vector3d rotationAxis2 = inertiaAxisXNew.cross(Eigen::Vector3d::UnitX());
		rotationAxis2.normalize(); 
		m2 = Eigen::AngleAxisd(angle2, rotationAxis2); 
	}
	if(angle>20){
		angle-=90;
	}
	angle = angle/180.0*M_PI;
	Eigen::Vector3d point[4];
	for(int i=0; i<4; i++){
		point[i].setZero();
	}
	/*
	point[0][0] = centerOfMass[0] + (width/2.0) * cos(angle) + (height/2.0) * sin(angle);
	point[0][2] = centerOfMass[2] + (width/2.0) * sin(angle) + (height/2.0) * cos(angle);

	point[1][0] = centerOfMass[0] + (-width/2.0) * cos(angle) + (height/2.0) * sin(angle);
	point[1][2] = centerOfMass[2] + (-width/2.0) * sin(angle) + (height/2.0) * cos(angle);

	point[2][0] = centerOfMass[0] + (width/2.0) * cos(angle) + (-height/2.0) * sin(angle);
	point[2][2] = centerOfMass[2] + (width/2.0) * sin(angle) + (-height/2.0) * cos(angle);

	point[3][0] = centerOfMass[0] + (-width/2.0) * cos(angle) + (-height/2.0) * sin(angle);
	point[3][2] = centerOfMass[2] + (-width/2.0) * sin(angle) + (-height/2.0) * cos(angle);
	*/
	double angle_tmp = angle;
	angle=-angle;
	double w = height;
	double h = width;

	point[0][0] = (w/2.0) * cos(angle) - (h/2.0) * sin(angle);
	point[0][2] = (w/2.0) * sin(angle) + (h/2.0) * cos(angle);

	point[1][0] = (-w/2.0) * cos(angle) - (h/2.0) * sin(angle);
	point[1][2] = (-w/2.0) * sin(angle) + (h/2.0) * cos(angle);

	point[2][0] = (w/2.0) * cos(angle) - (-h/2.0) * sin(angle);
	point[2][2] = (w/2.0) * sin(angle) + (-h/2.0) * cos(angle);

	point[3][0] = (-w/2.0) * cos(angle) - (-h/2.0) * sin(angle);
	point[3][2] = (-w/2.0) * sin(angle) + (-h/2.0) * cos(angle);

	angle = angle_tmp;

	RotatedRectangle rot_rect(width, height, angle, rect.center.x, rect.center.y);
	for(int i=0; i<4; i++){
		point[i] += Eigen::Vector3d(rect.center.x,0,rect.center.y);
		point[i] = m2.inverse() * point[i];
		point[i] = m.inverse() * point[i];
		point[i] = point[i] + centerOfMass;
		rot_rect.edges[i]=point[i];
	}
	rectangle = rot_rect;
	}
	return rectangle;
}

/**
 * Returns position of the surface.
 */
Eigen::Vector3d Plane::getPosition(){
	RotatedRectangle rect = getMaxExpansion();
	Eigen::Vector3d vec(-rect.pos_x, rect.pos_y, 0);
	vec = getQuaternion()*vec;
	Eigen::Vector3d result(centerOfMass[0], centerOfMass[1], centerOfMass[2]);
	result = result-vec;
	return result;
	
}

/**
 * Returns angle of the surface.
 */
Eigen::Quaterniond Plane::getPlaneAngle(){
	RotatedRectangle rect = getMaxExpansion();
	double angle = rect.angle;
	Eigen::AngleAxisd d = Eigen::AngleAxisd(angle, Eigen::Vector3d::Unit(2));
	Eigen::Quaterniond quat = getQuaternion();
	Eigen::Quaterniond result = quat*d;
	return result;	
}

/**
 * Angle calculation between points.
 */
double Plane::calculateAngle(Point &p0, Point &p1, Point &p2){

	return ((p1.x-p0.x)*(p2.z-p0.z))-((p2.x-p0.x)*(p1.z-p0.z));
	
}

/**
 * Draws the alpha shape with the help of edges and return the filled area as a bitmap.
 * Filled the space with Odd-Even algorithm. Using OpenCV Mat to get rid of QT
 */
boost::shared_ptr<cv::Mat> Plane::getHullCVImagesNofill(double &pixSize){

    boost::shared_ptr<cv::Mat> img(new cv::Mat(480,640,CV_8U));
    //Set image background to white
    uchar bkg_color(255);
    *img=bkg_color;

    std::list<std::list<std::pair<Shared_Point,Shared_Point> > >::iterator it;
    RotatedRectangle rot_rect = this->getMaxExpansion();    
    double scaleX = ((double)img->cols/rot_rect.width);
    double scaleY = ((double)img->rows/rot_rect.height);

    double scale = (scaleX<scaleY)?scaleX:scaleY;
    scale*=0.4;
    pixSize = (1.0)/scale*(1.0)/scale;

    std::list<std::pair<Shared_Point, Shared_Point> >::iterator it_list;
    std::list<std::list<std::pair<Shared_Point, Shared_Point> > > resultingEdges = getTransPoints();

    //draws only the unfilled hulls
   for(it=resultingEdges.begin(); it!=resultingEdges.end(); it++){
        for(it_list=(*it).begin(); it_list!=(*it).end(); it_list++){
        double x1 = (*it_list).first->x*scale+img->cols/2;
        double y1 = (*it_list).first->z*scale+img->rows/2;
        double x2 = (*it_list).second->x*scale+img->cols/2;
        double y2 = (*it_list).second->z*scale+img->rows/2;
        cv::line(*img,cv::Point(x1,y1),cv::Point(x2,y2),cv::Scalar(0));
        }
    }

   return img;
}


/**
 * Draws the alpha shape with the help of edges and return the filled area as a bitmap.
 * Filled the space with Odd-Even algorithm. Using OpenCV Mat to get rid of QT
 */
boost::shared_ptr<cv::Mat> Plane::getHullCVImages(double &pixSize){

    //the cv image for drawing the alpha shapes
    boost::shared_ptr<cv::Mat> img(new cv::Mat(480,640,CV_8U));
    //Sequence of contours for fillPoly
    std::vector<std::vector<cv::Point>> seq;
    //one single contour for fillPoly
    std::vector<cv::Point> contour;

    //Set image background to white
    uchar bkg_color(255);
    *img=bkg_color;

    std::list<std::list<std::pair<Shared_Point,Shared_Point> > >::iterator it;
    RotatedRectangle rot_rect = this->getMaxExpansion();    
    double scaleX = ((double)img->cols/rot_rect.width);
    double scaleY = ((double)img->rows/rot_rect.height);

    double scale = (scaleX<scaleY)?scaleX:scaleY;
    scale*=0.4;
    pixSize = (1.0)/scale*(1.0)/scale;

    std::list<std::pair<Shared_Point, Shared_Point> >::iterator it_list;
    std::list<std::list<std::pair<Shared_Point, Shared_Point> > > resultingEdges = getTransPoints();

    //draws filled hulls using openCV filled poly
    for(it=resultingEdges.begin(); it!=resultingEdges.end(); it++){
        contour.clear();
        for(it_list=(*it).begin(); it_list!=(*it).end(); it_list++){
            if (it_list==(*it).begin()){ //first pair? Take first and second point of the pair
                contour.push_back(cv::Point((*it_list).first->x*scale+img->cols/2,(*it_list).first->z*scale+img->rows/2));
                contour.push_back(cv::Point((*it_list).second->x*scale+img->cols/2,(*it_list).second->z*scale+img->rows/2));

            }else{
                //consequitive points contain only the second of the pair;
                contour.push_back(cv::Point((*it_list).second->x*scale+img->cols/2,(*it_list).second->z*scale+img->rows/2));
            }
        }
        seq.push_back(contour);
    }    


   /****test block
    contour.push_back(cv::Point(10,10));
    contour.push_back(cv::Point(10,100));
    contour.push_back(cv::Point(100,100));
    contour.push_back(cv::Point(100,10));

    seq.push_back(contour);
    contour.clear();
    contour.push_back(cv::Point(20,20));
    contour.push_back(cv::Point(20,50));
    contour.push_back(cv::Point(50,50));
    contour.push_back(cv::Point(50,20));

    seq.push_back(contour);
    end test block*****/


    uchar pen_color(0); //fill contour with black color
    cv::fillPoly(*img, seq, pen_color);


    return img;
}

/**
 * Calculates the size of the plane.
 */
double Plane::getSizeOfPlane(){
	double scale;
    boost::shared_ptr<cv::Mat> image = getHullCVImages(scale);
	//boost::shared_ptr<QImage> image = getHullImages(scale);
	double size;
	int value = 0;
	for(int row=0; row<image->rows; row++){
		for(int col=0; col<image->cols; col++){
            //QColor(image->pixel(x,y)) == Qt::white){
            //check if color is white at given pixel
            if(image->at<uchar>(row,col)==0){
				value++;
			}
		}
	}
	size = value*scale;
	return size;
}

/**
 * Returns the percentage squareness.
 */
double Plane::getRectangularness(){
	 RotatedRectangle rot_rect = getMaxExpansion();
	 return getSizeOfPlane()/(rot_rect.width*rot_rect.height);
}

/**
 * get the expansion of the plane in world cooridinates
 */
Eigen::Vector3d Plane::getSpatialExpansion()
{

    //hacked code. Only works for rectangular objects.
    //update if there will be time

    double minX=0;
    double maxX=0;
    double minY=0;
    double maxY=0;
    double minZ=0;
    double maxZ=0;

    Eigen::Vector3d diagonal;

    std::set<Shared_Point>::iterator iter=points.begin();

    minX=(*iter)->x;
    maxX=(*iter)->x;

    minY=(*iter)->y;
    maxY=(*iter)->y;

    minZ=(*iter)->z;
    maxZ=(*iter)->z;

    //get points have the largest distance in between (i.e. marking the plane diagonal)
    BOOST_FOREACH(Shared_Point sharedPoint,points){

        if (sharedPoint->x>maxX) maxX=sharedPoint->x;
        if (sharedPoint->y>maxY) maxY=sharedPoint->y;
        if (sharedPoint->z>maxZ) maxZ=sharedPoint->z;

        if (sharedPoint->x<minX) minX=sharedPoint->x;
        if (sharedPoint->y<minY) minY=sharedPoint->y;
        if (sharedPoint->z<minZ) minZ=sharedPoint->z;
    }


    diagonal=Eigen::Vector3d(maxX,maxY,maxZ)-Eigen::Vector3d(minX,minY,minZ);

    //std::cout<<diagonal.x()<<std::endl;
    //std::cout<<diagonal.y()<<std::endl;
    //std::cout<<diagonal.z()<<std::endl;

    return diagonal;
}


/*
std::istream& operator >>(std::istream &is, Plane &obj){
	for(int i=0; i<3; i++){
		is>>obj.centerOfMass[i];
	}
	for(int i=0; i<3; i++){
		is>>obj.inertiaAxisX[i];
	}
	for(int i=0; i<3; i++){
		is>>obj.inertiaAxisY[i];
	}
	for(int i=0; i<3; i++){
		is>>obj.normalVector[i];
	}
	for(int i=0; i<3; i++){
		is>>obj.S[i];
	}
	for(int i=0; i<9; i++){
		is>>obj.C[i];
	}
	is>>obj.bias;
	for(int i=0; i<9; i++){
		is>>obj.productMatrix[i];
	}
	int size;
	is >> size;

	printf("Size is: %i\n",size);

	for(int i=0; i<size; i++){
		Shared_Point p = Shared_Point(new Point);
		is>>p->x >> p->y >> p->z;
		obj.points.insert(p);
	}
	obj.rectangleCalculated = false;
	
	return is;
	
}

std::ostream& operator <<(std::ostream &os, const Plane &obj){
	for(int i=0; i<3; i++){
		os<<obj.centerOfMass[i] << " ";
	}
	for(int i=0; i<3; i++){
		os<<obj.inertiaAxisX[i] << " ";
	}
	for(int i=0; i<3; i++){
		os<<obj.inertiaAxisY[i] << " ";
	}
	for(int i=0; i<3; i++){
		os<<obj.normalVector[i] << " ";
	}
	for(int i=0; i<3; i++){
		os<<obj.S[i] << " ";
	}
	for(int i=0; i<9; i++){
		os<<obj.C[i] << " ";
	}
	os<<obj.bias << " ";
	for(int i=0; i<9; i++){
		os<<obj.productMatrix[i] << " ";
	}
	os<<obj.points.size()<< " ";
	std::set<Shared_Point>::const_iterator it;
	for(it=obj.points.begin(); it!=obj.points.end(); it++){
		os<<(*it)->x<< " " << (*it)->y << " " <<  (*it)->z << " " ;
	}
	return os;

}
*/
