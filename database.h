#pragma once

#include "TriMesh.h"
#include "GeometryTools.h"

class database
{
protected:
	static database* instance;
	static int refCount;

public:
	database(void);
	~database(void);

	static database* getInstance();
	static void Release();
	static void Destroy();

	void PN_TRIANGLE(int level_num);
	void exact_pn_triangle(int level_num);
	void exact_torus(int level_num);
	void initalize_cubic_bezier();
	void clear();
	void save();
	void read();

	QImage *textureimage; //texture image
	TriMesh *themesh;
	TriMesh *cubic_bezier;
};

