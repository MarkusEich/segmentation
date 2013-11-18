#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_hierarchy_2.h>
#include <CGAL/Alpha_shape_2.h>
#include <CGAL/Polygon_2.h>
#include <iostream>
#include <fstream>
#include <list>
#include <cassert>
#include <segmentation/hull.h>
#include <segmentation/regionGrowing.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polygon_2<K> Polygon_2;

//typedef for  AlphaShape 2D
typedef CGAL::Alpha_shape_vertex_base_2<K>                      Avb;
typedef CGAL::Triangulation_hierarchy_vertex_base_2<Avb>        Av; 
typedef CGAL::Triangulation_face_base_2<K>                      Tf;
typedef CGAL::Alpha_shape_face_base_2<K,Tf>                     Af;

//typedef CGAL::Triangulation_data_structure_2<Avb,Af>            Tds;
typedef CGAL::Triangulation_default_data_structure_2<K,Av,Af>   Tds;
typedef CGAL::Delaunay_triangulation_2<K,Tds>                   Dt;
typedef CGAL::Triangulation_hierarchy_2<Dt>                     Ht;
typedef CGAL::Alpha_shape_2<Ht>                                 Alpha_shape_2;
typedef CGAL::Quotient<CGAL::MP_Float> 							NT;

//typedef K::Point_2                                  Point;
typedef Alpha_shape_2::Face_handle          		Face_handle;
typedef Alpha_shape_2::Face                 		Face;
typedef Alpha_shape_2::Edge                 		Edge;
typedef Alpha_shape_2::Alpha_iterator               Alpha_iterator;
//typedef Alpha_shape_2::NT                           NT;


using namespace segmentation;

/**
 * Calculates the alpha shape in 2D.
 */

std::vector<Shared_Point> Hull::alphaShape2D(const std::list<Point> points,double alpha){
	std::vector<Shared_Point> hull;
	if(points.size() < 3) return hull;

	#ifdef NO_HULL_CALCULATION
	for(std::list<Point>::const_iterator it=points.begin();it!= points.end();it++){
		hull.push_back(Shared_Point(new Point(it->p[0], 0, it->p[2])));
		return hull;		
	}
	#endif

	//iHt dt;
	std::vector<K::Point_2> dt;
	dt.reserve(points.size());

	for(std::list<Point>::const_iterator it=points.begin(); it!=points.end(); it++){
		K::Point_2 p((*it).p[0],(*it).p[2]);
		dt.push_back(p);
	}
	//std::cout << "Delaunay computed." << std::endl;
	//const NT alpha_value = 1;

	//compute alpha shape
	//printf("Size of input %i\n",(int)dt.size());
	try{
        Alpha_shape_2 as(dt.begin(),dt.end(),alpha);//Alpha_shape_2::REGULARIZED);
		//std::cout << "Alpha shape computed in REGULARIZED mode by defaut."<< std::endl;

		for(Alpha_shape_2::Alpha_shape_vertices_iterator it = as.alpha_shape_vertices_begin(); it != as.alpha_shape_vertices_end(); it++){
			Alpha_shape_2::Vertex_handle vh = *it; 
			K::Point_2& p = vh->point();
			Shared_Point p2(new Point(p.x(),0,p.y()));
			hull.push_back(p2); 
		}
		hull = sortAlphaShapePoints(hull);
	}catch( ...){
		printf("Catched Unknown Cgal Exception adding plane without generatig shape\n");
		for(std::list<Point>::const_iterator it = points.begin(); it != points.end() ; it++){
			hull.push_back(Shared_Point(new Point(it->p[0],it->p[1],it->p[2])));
		}
	}
	//printf("Hull ready \n");
	return hull;

}

/**
 * Returns an alpha shape but as lines instead of only points.
 */
std::list<std::pair<Shared_Point,Shared_Point> > Hull::alphaHull2D(const std::list<Point> points,double alpha){
	//iHt dt;
	std::list<std::pair<Shared_Point,Shared_Point> > resultingEdges;
	if(points.size() < 3) return resultingEdges;

	#ifdef NO_HULL_CALCULATION
	if(points.size() == 4){
		std::list<Point>::const_iterator it = points.begin();
		Shared_Point p0(new Point(it->p[0], 0, it->p[2]));
		it++;
		Shared_Point p1(new Point(it->p[0], 0, it->p[2]));
		it++;
		Shared_Point p2(new Point(it->p[0], 0, it->p[2]));
		it++;
		Shared_Point p3(new Point(it->p[0], 0, it->p[2]));

		resultingEdges.push_back(std::pair<Shared_Point,Shared_Point>(p0,p1));
		resultingEdges.push_back(std::pair<Shared_Point,Shared_Point>(p1,p2));
		resultingEdges.push_back(std::pair<Shared_Point,Shared_Point>(p2,p3));
		resultingEdges.push_back(std::pair<Shared_Point,Shared_Point>(p3,p0));
		
		return resultingEdges;
	}else{
			fprintf(stderr,"Cannot Compute no hull\n");
	}
	#endif

	std::vector<K::Point_2> dt;
	dt.reserve(points.size());

	for(std::list<Point>::const_iterator it=points.begin(); it!=points.end(); it++){
		K::Point_2 p((*it).p[0],(*it).p[2]);
		dt.push_back(p);
	}

	try{
        Alpha_shape_2 as(dt.begin(),dt.end(),alpha);//Alpha_shape_2::REGULARIZED);
		//Alpha_shape_2 as(dt.begin(),dt.end(),Alpha_shape_2::REGULARIZED);
		//std::cout << "Alpha shape computed in REGULARIZED mode by defaut."<< std::endl;

		//Alpha_iterator opt = as.find_optimal_alpha(1);
		//as.set_alpha(*opt);

		for(Alpha_shape_2::Alpha_shape_edges_iterator it = as.alpha_shape_edges_begin(); it != as.alpha_shape_edges_end(); it++){
			Edge e = (*it);
			Face_handle f = e.first;
			int i = e.second;

			Alpha_shape_2::Vertex_handle vh1 = f->vertex(f->cw(i));
			Alpha_shape_2::Vertex_handle vh2 = f->vertex(f->ccw(i));
			K::Point_2& p1 = vh1->point();
			K::Point_2& p2 = vh2->point();
			
			std::pair<Shared_Point,Shared_Point> pair(Shared_Point(new Point(p1.x(),0,p1.y())),Shared_Point(new Point(p2.x(),0,p2.y())));
			resultingEdges.push_back(pair);
		}

	}catch( ...){
		printf("Catched Unknown Cgal Exception adding plane without generatig shape\n");
	}
	return resultingEdges;

}

/**
 * Sorts all points on the angle (like graham scan).
 */
std::vector<Shared_Point> Hull::sortFromAngle(const std::vector<Shared_Point> hull_points){
	std::list<Shared_Point> sorted;
	std::vector<Shared_Point> result;
	if(hull_points.size()<3){
		return result;
	}

    for(size_t i=0; i<hull_points.size(); i++){
		sorted.push_back(hull_points[i]);
	}

	sorted.sort(angleSorter);
	for(std::list<Shared_Point>::const_iterator it=sorted.begin(); it!=sorted.end(); it++){
		result.push_back(*it);
	}
	return result;


}

/**
 * Function, that can be calld by std::list.sort.
 */
bool Hull::angleSorter(const Shared_Point p1, const Shared_Point p2){
	double angle1 = atan2(1-p1->z,-p1->x);
	double angle2 = atan2(1-p2->z,-p2->x);
	return angle1<=angle2;
}

/**
 * Sorting function for 2D point sets, which should result in a hull.
 */
std::vector<Shared_Point> Hull::sortAlphaShapePoints(const std::vector<Shared_Point> hull_points){
	std::list<Shared_Point> unsorted;
	std::vector<Shared_Point> sorted;
	if(hull_points.size()<3){
		return sorted;
	}

    for(size_t i=0; i<hull_points.size(); i++){
		unsorted.push_back(hull_points[i]);
	}
	Shared_Point point = *unsorted.begin();
	sorted.push_back(point);
	unsorted.remove(point);
	std::list<Shared_Point>::const_iterator it;
	while(!unsorted.empty()){
		Shared_Point p = sorted[sorted.size()-1];
		double minDist = 999999;
		Shared_Point pMin;
		for(it=unsorted.begin(); it!=unsorted.end(); it++){
			double dist = sqrt((p->x-(*it)->x)*(p->x-(*it)->x)+(p->z-(*it)->z)*(p->z-(*it)->z));
			if(dist < minDist){
				minDist = dist;
				pMin = (*it);
			}
		}
		
		sorted.push_back(pMin);
		unsorted.remove(pMin);
	}
	return sorted;
}

/**
 * Sorts the pairs according to the connectivity.
 */
std::list<std::list<std::pair<Shared_Point,Shared_Point> > > Hull::sortedPair(const std::list<Point> points,double alpha){
	std::list<std::list< std::pair<Shared_Point, Shared_Point> > > sortedEdges;
	std::list<std::pair<Shared_Point, Shared_Point> >::const_iterator it;
    std::list<std::pair<Shared_Point, Shared_Point> > edges = alphaHull2D(points,alpha);
	
	if(edges.begin() == edges.end()){
		fprintf(stderr,"Warning sortedPair dont get any hull, so skipping calculation\n");
		return sortedEdges;
	}
	bool addElement = true;
	int count = 0;
	while(edges.size()!=0){
		std::list<std::pair<Shared_Point,Shared_Point> > current_list;
		current_list.push_back(*edges.begin());
		edges.remove(*edges.begin());
		//printf("Start Searching\n");
		while(addElement){
			addElement = false;
			for(it=edges.begin(); it!=edges.end(); it++){
				std::pair<Shared_Point, Shared_Point> last = (current_list.back());
				//printf("Last Element: %f,%f,%f\n",last.first->p[0],last.first->p[1],last.first->p[2]);
				if(*(*it).first==*(last.second)){ 
					current_list.push_back(std::pair<Shared_Point, Shared_Point>((*it).first, (*it).second));
					addElement = true;
					count++;
					edges.remove(*it);
					break;
				}else if(*(*it).second==*(last.second)){
					current_list.push_back(std::pair<Shared_Point, Shared_Point>((*it).second, (*it).first));
					addElement = true;
					edges.remove(*it);
					count++;
					break;
				}
			}
		}
		//this should be commented in, but is currently not fully functional
		//if(*current_list.front().first==*current_list.back().second || *current_list.back().first==*current_list.front().second){
			sortedEdges.push_back(current_list);
		//}
	}
	return sortedEdges;
}
