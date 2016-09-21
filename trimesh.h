/***********************************************************
Filename: TriMesh.h

Description:

***********************************************************/

#ifndef   TRIMESH_H
#define  TRIMESH_H

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/IO/Options.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Geometry/VectorT.hh>
#include <QGLViewer/qglviewer.h>
#include <QMessageBox>
#include <math.h>
#include <windows.h>
#include <GL.h>
#include <GL/GLU.h>
#include <GLAUX.H>
#include "GlobalConstant.h"
#include <vector>
#include <limits>
#include <algorithm>

struct TriMeshTraits : public OpenMesh::DefaultTraits
{
	typedef OpenMesh::Vec3d Point;
	typedef OpenMesh::Vec3d Normal;

	VertexAttributes  ( OpenMesh::Attributes::Normal );
	FaceAttributes      ( OpenMesh::Attributes::Normal );
};

class TriMesh: public OpenMesh::TriMesh_ArrayKernelT<TriMeshTraits>
{
public:
	TriMesh(void);
	virtual ~TriMesh(void);

	double getAverageEdgeLength();
	void  needBoundingBox();
	Point getSceneCenter() { return (bbox_min + bbox_max) / 2.0; };
	double getSceneRadius() { return (bbox_max-bbox_min).norm() / 2.0; };
	TriMesh::VertexHandle approximatePointToVertex(Point center);
	bool selectInitialPoint(Point center);
	void init_boundary();
	void set_vertex_color(int R, int G, int B);
	bool barycentric_coord(FaceHandle fh, Point the, double& alpha, double& beta, double& gamma);
	//Cubic Bezier patch
	void EvaluateBezierPoint(double coeff[][3], double normCoeff[][3], double b1, double b2, double b3, Point& destVertices, Normal& destNormals);
	void Cubic_Bezier_Parameters(FaceHandle fh, double coeff[ ][3], double normCoeff[ ][3]);
	//evaluate surface Point and Normal from Parameters
	void para_to_point(FaceHandle fh, double alpha,double beta,double gamma, Point &p, Normal &n);
	void para_to_point_PN_triangle(FaceHandle fh, double u,double v,double w, Point &p, Normal &n);

	//get mesh statistic
	void getMeshStat(std::vector<double> &meshStat, TriMesh* meshCopy);
	static TriMesh* read(const char* filename, OpenMesh::IO::Options* opt = NULL);
	static bool save(const char* filename, TriMesh* mesh = NULL, OpenMesh::IO::Options* opt = NULL);
	void write();


public:
	double eps_;
	double average_edge_length;
	std::vector<unsigned int> vertex_group0; // All vertices for drawing
	std::vector<unsigned int> vertex_group1; //Vertices with texture
	std::vector<unsigned int> vertex_group2; //Vertices without texture
	std::vector<int> boundary;
private:
	Point bbox_min, bbox_max;

};

#endif // TRIMESH_H
