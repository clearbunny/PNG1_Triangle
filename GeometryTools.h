/***********************************************************

Filename: GeometryTools.h

Description:

Collection of computational geometry tools: check if a point is inside of a
triangle or a polygon, project a point to a plane, distance between point and
line/plane, span a plane from 3 points


***********************************************************/

#ifndef  GEOMETRYTOOLS_H
#define  GEOMETRYTOOLS_H

//#include "TriMesh.h"
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/IO/Options.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

typedef OpenMesh::Vec3d Point;

#define SQ_DIST3(a, b)          (((a)[0]-(b)[0])*((a)[0]-(b)[0]) +      \
	((a)[1]-(b)[1])*((a)[1]-(b)[1]) +      \
	((a)[2]-(b)[2])*((a)[2]-(b)[2]))

#define DOTPROD3(a, b)		 ((a)[0]*(b)[0] + (a)[1]*(b)[1] + (a)[2]*(b)[2])

#define CROSSPROD3(a,b,c)       {(a)[0]=(b)[1]*(c)[2]-(b)[2]*(c)[1]; \
	(a)[1]=(b)[2]*(c)[0]-(b)[0]*(c)[2]; \
	(a)[2]=(b)[0]*(c)[1]-(b)[1]*(c)[0];}

class GeometryTools
{
public:
	GeometryTools(void);
	~GeometryTools(void);

	static void normalize(Point& p);
	static double length2(Point p);
	static double length(Point p);
	static double dotProduct(Point p1, Point p2);
	static Point crossProduct(Point p1, Point p2);
	// Check if point is inside of a triangle: http://mathworld.wolfram.com/TriangleInterior.html
	// p0, p1 and p2: three points consisting of a triangle; p: point in the same plane
	// return true if it's inside
	static bool isInsideOfTriangle(Point p0, Point p1, Point p2, float p[]);
	static bool isInsideOfTriangle(Point p0, Point p1, Point p2, Point p);
	//point-point distance
	static double distPoint2Point(Point p1, Point p2);
	static double squareDistPoint2Point(Point p1, Point p2);
	// point-line distance: http://mathworld.wolfram.com/Point-LineDistance2-Dimensional.html
	// p1 and p2: two points consisting of a line; p: point in the space
	// return the distance between p and this line
	static double distPoint2Line(Point p1, Point p2, float p[]);
	static double distPoint2Line(Point p1, Point p2, Point p);
	// point-plane distance: http://mathworld.wolfram.com/Point-PlaneDistance.html
	// plane[4]: plane parameters, ax+by+cz+d=0, p: point in the space
	// return the SIGNED distance between p and plane
	static double distPoint2Plane(double plane[4], Point p);
	// Given three points, span the plan by cross product
	// http://local.wasp.uwa.edu.au/~pbourke/geometry/linefacet/
	// p1, p2 and p3: three points in space; result: abcd in ax+by+cz+d=0
	static void spanPlane(Point p1, Point p2, Point p3, double result[4]);
	// project a point in 3d space to a plane, which is defined as ax+by+cz+d=0
	// http://local.wasp.uwa.edu.au/~pbourke/geometry/linefacet/
	// plane[4]: ax+by+cz+d=0; p1: a point in space; p: project p1 to the plane
	static void projectPoint2Plane(double plane[4], Point p1, Point &p);
	static void project(Point &p, Point normal);
	// check if point is inside of polygon
	// http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/isInsideOfPolygon.html
	// vertx and verty: vertices on polygon; testx and testy: test polygon; nvert: number of vertices
	static int isInsideOfPolygon(int nvert, float *vertx, float *verty, float testx, float testy);

	//// radius of inscribed circle
	//// http://en.wikipedia.org/wiki/Incircle
	//static double radiusOfInscribedCircle(TriMesh *mesh, FaceHandle &fh);
	//// radius of circumscribed circle
	//// http://en.wikipedia.org/wiki/Circumcircle
	//static double radiusOfCircumscribedCircle(TriMesh *mesh, FaceHandle &fh);
	//// return minimum angle in the triangle represented by fh
	//static double minAngleInTriangle(TriMesh *mesh, FaceHandle &fh);

	// test if line segment(s1p1, s1p2) and (s2p1, s2p2) are intersect and if yes, the intersect point is returned as Point p
	static bool intersect_segment_to_segment(Point A, Point B, Point C, Point D, Point& p);
	// test if Point t1 and t2 are in the same side of line (p1, p2)
	static bool is_same_side(Point p1, Point p2, Point t1, Point t2);
	// compute barycentric coordinate of Point p in the Triangle ( p0, p1, p2) 
	static bool barycentrc_coord(Point curgeo, Point A, Point B, Point C, Point& coord);
};
#endif
