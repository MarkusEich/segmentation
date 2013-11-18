/*****************************************
 ** Markus Eich (2013)
 ** Malgorzata Goldhoorn (2010)
 ** DFKI GmbH
 **
 **
*/


#include <iostream>
#include <math.h>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <segmentation/kdtree.hpp>
#include <segmentation/plane.h>
#include <segmentation/regionGrowing.h>
#include <segmentation/hull.h>

#ifdef USE_OPENMP
#include <omp.h>
#endif

using namespace Gamera::Kdtree;
using namespace std;
using namespace segmentation;

int RegionGrowing::neighbours_count;

RegionGrowing::RegionGrowing() : delta(0), gamma(0), epsilon(0), theta(0), n_p (0),  funcPtr(0), neighboursCalculated(false)
{


}

RegionGrowing::~RegionGrowing(){ 
	points.clear();
}

/**
 * Chooses starting points.
 */
void RegionGrowing::selectP1P2(Shared_Plane pl) {
	//only if there are two points existing
	while(n_p<=points.size()-2) {
		//add point to the plane, take the nearest neighbors of the point, then move the point to the end of the point set,
		//that future planes not from this point sample
		pl->mergeNearestNeighbours(points[0]);
		pl->addPoint(points[0]);
		moveToEnd(0);

		//search nearest neighbors to the point, which still is not in the plane
		Shared_Point p2 = pl->findNearestNeighbour();
		
		if(p2) {
			movePtrToEnd(p2);
			pl->addPoint(p2);
			pl->mergeNearestNeighbours(p2);
			return;
		}
	}
}

/**
 * Sets the pointcloud, which should be processed.
 */
void RegionGrowing::setPointcloud(std::vector<Shared_Point>& points){
	mutex.lock();
	this->points=points;
	mutex.unlock();
}

/**
 * Sets the pointcloud as Point_3d array, which should be processed.
 */
void RegionGrowing::setPointcloud(Pointcloud::Point_3d* points, const int count){
	mutex.lock();
	neighboursCalculated=false;

    for(size_t i=0; i<this->points.size(); i++){
		this->points[i].reset();
	}
	this->points.clear();

	for(int i=0; i<count; i++){
		Shared_Point p(new Point(points[i].x, points[i].y, points[i].z));
		this->points.push_back(p);
	}
	mutex.unlock();
}

/**
 * Sets the callback function.
 */


/**
 * Moves a point to the end of the point array, move the local point to the old position of the point.
 */
void RegionGrowing::movePtrToEnd(Shared_Point p) {
	moveToEnd(p->index);

}

void RegionGrowing::moveToEnd(size_t i) {
	size_t tmpi = points.size()-n_p-1;

	Shared_Point tmpp = points[i];
	points[i] = points[tmpi];
	points[tmpi] = tmpp;
	points[tmpi]->index = tmpi;
	points[i]->index = i;
	
	n_p++;
}

	
void RegionGrowing::calculateNeighbours(std::vector<Shared_Point> &points, double delta){
	//all points entered in the tree, they need to be converted to the internal messaging format gamera
	KdnodeVector nodes;
	nodes.reserve(points.size());
	std::vector<Shared_Point>::iterator it;
	size_t ind = 0;

	for(it = points.begin(); it != points.end(); it++) {
		//register index of the point in the point array, this will make a mapping of the points on their array position possible
		//(necessary for the reordering of the points)
		(*it)->index = ind++;
		CoordPoint p;
		p.reserve(3);
		for(int i = 0; i < 3; i++)
			p.push_back((*it)->p[i]);
		//in the node element can also store a pointer to the point in our format,
		//we may find that after the calculation quickly
		nodes.push_back(kdnode(p, &(*it)));
	}

	//create tree of the above assembled nodes
	kdtree tree(&nodes);

	//iterating over all the points and find its nearest neighbors
	//the first one is (probably) always be the point itself
	for(it = points.begin(); it != points.end(); it++) {
		KdnodeVector neighbours;
		CoordPoint ref;
		ref.reserve(3);
		for(int i = 0; i < 3; i++)
			ref.push_back((*it)->p[i]);
		tree.k_nearest_neighbors(ref, neighbours_count, &neighbours);
		KdnodeVector::iterator rit;
		//going throught neighbors to consider, whether they really are not the reference point and whether the distance is actually 
		//smaller than delta
		//if so, it will be added to the list of nearby points of the point
		for(rit = neighbours.begin(); rit != neighbours.end(); rit++) {
			Shared_Point p2 = *reinterpret_cast<Shared_Point*>(rit->data);
			if(p2 != *it && (*it)->distance(*p2) < delta) {
				(*it)->nearestNeighbours.push_back(Point::Neighbour(p2, (*it)->distance(*p2)));
			}

		}
	}
}

/**
 * Precomputation of the nearest neighbors by k-d tree.
 */
void RegionGrowing::precomputeNearestNeighbours() {
	if(neighboursCalculated){
		return;
	}
	RegionGrowing::calculateNeighbours(points, delta);
	neighboursCalculated = true;
	//std::cout << "finished precomputing nearest neighbours" << std::endl;
}

/**
 * Actual search function, gets a reference to a list handed over by planes.
 */
bool RegionGrowing::regionFind(std::list<Shared_Plane>& result) {
	if(!mutex.try_lock()){ //Only sucsessful if no setters are active, so it prevent race-conditidions on the pointlists.
		printf("Cannot lock lists, so skipping RegionGrowing\n");
		return false;
	}

	std::cout<< "Using following paramters: delta , epsilon , theta , gamma , KNN "<<std::endl;
	std::cout<<  delta << " " << epsilon << " "<< epsilon<< " "<< theta << " "<< gamma <<" "<< neighbours_count<<std::endl;

	if(points.empty()){
		mutex.unlock();
		printf("Points are Empty, cannot start\n");
		return false;
	}

	n_p = 0;
	//precalculate neighbors of the points 
	
	std::cout << "Start Calculation of N.N."<< std::endl;
	precomputeNearestNeighbours();
	std::cout << "Start computing of " <<  points.size() << " points " << std::endl;

	//as long as enough points are available

	#ifdef USE_OPENMP
	#pragma omp parallel
	{
	#pragma omp signle private(n_p,maxiter,)
	{
	#endif

	while(n_p <= points.size()-2) {
		int maxiter = points.size();
		Shared_Plane cur_plane(new Plane);
		//first point from the vector select, this is by no means in a pre-existing plane
		//(due to reordering of the points in the array!)
		
		selectP1P2(cur_plane);
		//as long as the plane does not nearest neighbors and not all points were processed
		while(cur_plane->knowNeighbours() && n_p < points.size() && maxiter--) {
			//read nearest neighbors (in constant time) and check whether this is already contained in the plane
			//if such a point exists, it is added to the plane
			Shared_Point pprime = cur_plane->findNearestNeighbour();
			if(pprime.get()!=0 && pprime->index < points.size()-n_p) {
			
				//has the point maximum gamma distance to the plane?
				if(cur_plane->planeDist(pprime) < gamma) {
					//save the state of the plane and add point, because of the stored state, we can afterward go back to the old state,
					//if the point not fit to the plane (MSE is too large)
					//simple backtracking
					cur_plane->saveState();
					cur_plane->mergeNearestNeighbours(pprime);
					cur_plane->addPoint(pprime);

					//check with the added data, if MSE is too large, if so, remove the point again and saved the state
					//and restore plane
					if(cur_plane->mse() >= epsilon) {
						cur_plane->removePoint(pprime);
						cur_plane->restoreState();
					} else {
						//if not, the point will be moved to the end of the point list and its neighbors in the queue of the plane
						movePtrToEnd(pprime);
					}
				}	
			}
		}

		//calculation of a plane is completed, check that it contains enough points,
		//if so, add to the result, if not, discard
        if(static_cast<int>(cur_plane->points.size()) > theta) {
			//std::cout << "Found plane with " << cur_plane->points.size() << " points." << std::endl;
			result.push_back(cur_plane);
			if(funcPtr!=0)
				(*funcPtr)(cur_plane);
			//same as above but with proper boost
			if(callback_!=0)
				callback_(cur_plane);

		}
	}
	#ifdef USE_OPENMP
	}
	}
	#endif

	//printf("Calculation ended\n");
	mutex.unlock();
	return true;
}

void RegionGrowing::setDelta(double v){
	neighboursCalculated=false;
	delta = v;
	fprintf(stdout,"New Delta %f\n",v);
}

void RegionGrowing::setGamma(double v){
	gamma = v;
	fprintf(stdout,"New Gamma %f\n",v);
}

void RegionGrowing::setEpsilon(double v){
	epsilon = v;
	fprintf(stdout,"New Epsilon %f\n",v);
}

void RegionGrowing::setTheta(int v){
	theta = v;
	fprintf(stdout,"New Theta %i\n",v);
}

void RegionGrowing::setNeighbours(int v){
	neighboursCalculated=false;
	RegionGrowing::neighbours_count = v;
	fprintf(stdout,"New Neighnour Count %i\n",v);
}
