#include "TriMesh.h"
#include "GeometryTools.h"

typedef OpenMesh::Vec3d Point;
typedef OpenMesh::Vec3d Vec;
typedef OpenMesh::Vec3d Normal;

double startu,startv;

TriMesh::TriMesh(void)
{
	eps_ = 1e-12;
	average_edge_length=0.0;
}

TriMesh::~TriMesh(void)
{
	if(&boundary)
		boundary.clear();

	if (&vertex_group0)
		vertex_group0.clear();

	if (&vertex_group1)
		vertex_group1.clear();

	if (&vertex_group2)
		vertex_group2.clear();
}

double TriMesh::getAverageEdgeLength(){
	TriMesh::HalfedgeIter he_it=halfedges_begin();
	Point from, to;
	double edge_length=0.0;
	for ( ; he_it!=halfedges_end(); ++he_it)
	{
		from=point(from_vertex_handle(he_it.handle()));
		to=point(to_vertex_handle(he_it.handle()));
		edge_length+=GeometryTools::distPoint2Point(from, to);
	}
	return edge_length/n_halfedges();
}

void TriMesh::needBoundingBox()
{
	TriMesh::ConstVertexIter vIt = vertices_begin();
	TriMesh::ConstVertexIter vEnd = vertices_end();

    bbox_min = bbox_max = point(vIt);

    for (; vIt != vEnd; ++vIt) {
		const Point& p = point(vIt);
		bbox_min.minimize(p);
		bbox_max.maximize(p);
    }
}

void TriMesh::set_vertex_color(int R, int G, int B){
	Color color(R, G, B);
	TriMesh::VertexIter v_it, v_end(vertices_end());
	for(v_it=vertices_begin();v_it!=v_end;++v_it)
		this->set_color(v_it,color);
	
}

bool TriMesh::barycentric_coord(FaceHandle fh, Point the, double& alpha, double& beta, double& gamma){
	FaceVertexIter fv_it=fv_iter(fh);
	Point A=point(fv_it.handle());
	++fv_it;
	Point B=point(fv_it.handle());
	++fv_it;
	Point C=point(fv_it.handle());
	if (!(A[0]==C[0]&&B[0]==C[0]))
	{
		if ((B[0]-C[0])*(A[1]-C[1]+A[2]-C[2])==(A[0]-C[0])*(B[1]-C[1]+B[2]-C[2]))
		{
			double ab=sqrt(SQ_DIST3(A,B));
			double bc=sqrt(SQ_DIST3(B,C));
			double ac=sqrt(SQ_DIST3(C,A));

			double a=sqrt(SQ_DIST3(B,C));
			double b=sqrt(SQ_DIST3(B, the));
			double c= sqrt(SQ_DIST3(C,the));

			double aa=sqrt(SQ_DIST3(A,C));
			double bb=sqrt(SQ_DIST3(A, the));
			double cc= sqrt(SQ_DIST3(C,the));
			if( (ab+bc>ac)&&(ab+ac>bc)&&(ac+bc>ab)&&(a+b>c)&&(b+c>a)&&(a+c>b)&&(aa+bb>cc)&&(bb+cc>aa)&&(cc+aa>bb))
			{
				alpha=sqrt((a+b+c)/2.0*((a+b+c)/2.0-a)*((a+b+c)/2.0-b)*((a+b+c)/2.0-c))/sqrt((ab+bc+ac)/2.0*((ab+bc+ac)/2.0-ab)*((ab+bc+ac)/2.0-bc)*((ab+bc+ac)/2.0-ac));
				beta=sqrt((aa+bb+cc)/2.0*((aa+bb+cc)/2.0-aa)*((aa+bb+cc)/2.0-bb)*((aa+bb+cc)/2.0-cc))/sqrt((ab+bc+ac)/2.0*((ab+bc+ac)/2.0-ab)*((ab+bc+ac)/2.0-bc)*((ab+bc+ac)/2.0-ac));
				gamma=1.0-alpha-beta;
				return true;
			}else{
				//not a triangle, there is two edge length sum equal to the third edge length
				//QMessageBox::information(NULL, "not a triangle in barycentric_coord", "three point in the same line");
				return false;
			}
		}else{
			if (((B[0]-C[0])*(C[1]-the[1]+C[2]-the[2])-(C[0]-the[0])*(B[1]-C[1]+B[2]-C[2]))==0)
				alpha=0.0;
			else
				alpha=((B[0]-C[0])*(C[1]-the[1]+C[2]-the[2])-(C[0]-the[0])*(B[1]-C[1]+B[2]-C[2]))
				/((A[0]-C[0])*(B[1]-C[1]+B[2]-C[2])-(B[0]-C[0])*(A[1]-C[1]+A[2]-C[2]));
			if (((A[0]-C[0])*(C[1]-the[1]+C[2]-the[2])-(C[0]-the[0])*(A[1]-C[1]+A[2]-C[2]))==0)
				beta=0.0;
			else
				beta=((A[0]-C[0])*(C[1]-the[1]+C[2]-the[2])-(C[0]-the[0])*(A[1]-C[1]+A[2]-C[2]))
				/((B[0]-C[0])*(A[1]-C[1]+A[2]-C[2])-(A[0]-C[0])*(B[1]-C[1]+B[2]-C[2]));
			gamma=1.0-alpha-beta;
			return true;
		}
	}else{
		Point ba,bc,bd,cb,ca,cd;
		ba=A-B;
		bc=C-B;
		bd=the-B;
		cb=B-C;
		ca=A-C;
		cd=the-C;
		double test0=SQ_DIST3(A,B)*SQ_DIST3(C,B)-DOTPROD3(ba,bc)*DOTPROD3(ba,bc);
		double test1=(SQ_DIST3(the,B)*SQ_DIST3(C,B)-DOTPROD3(bd,bc)*DOTPROD3(bd,bc));
		if ((SQ_DIST3(the,B)*SQ_DIST3(C,B)-DOTPROD3(bd,bc)*DOTPROD3(bd,bc))<eps_)
			alpha=0.0;
		else
			alpha=sqrt(SQ_DIST3(the,B)*SQ_DIST3(C,B)-DOTPROD3(bd,bc)*DOTPROD3(bd,bc))
			/sqrt(SQ_DIST3(A,B)*SQ_DIST3(C,B)-DOTPROD3(ba,bc)*DOTPROD3(ba,bc));

		if ((SQ_DIST3(the,C)*SQ_DIST3(C,A)-DOTPROD3(cd,ca)*DOTPROD3(cd,ca))<eps_ )
			beta=0.0;
		else
			beta=sqrt(SQ_DIST3(the,C)*SQ_DIST3(C,A)-DOTPROD3(cd,ca)*DOTPROD3(cd,ca))
			/sqrt(SQ_DIST3(C,B)*SQ_DIST3(C,A)-DOTPROD3(cb,ca)*DOTPROD3(cb,ca));
		gamma=1.0-alpha-beta;
		return true;
	}	
}


void TriMesh::EvaluateBezierPoint(double coeff[][3], double normCoeff[][3], double b1, double b2, double b3, Point &destVertices, Normal& destNormals)
{
	int i;
	int x, y;

	float buff[6][3];
	float tmpf;
	x = 0;
	y = 1;
	for (i=0; i<6; i++)	{
		buff[i][0] = coeff[x][0]*b2 + coeff[y][0]*b1 + coeff[(y+1)][0]*b3;
		buff[i][1] = coeff[x][1]*b2 + coeff[y][1]*b1 + coeff[(y+1)][1]*b3;
		buff[i][2] = coeff[x][2]*b2 + coeff[y][2]*b1 + coeff[(y+1)][2]*b3;
		x++;
		if ((i == 0) || (i == 2))
			y += 2;
		else
			y++;
	}
	x = 0;
	y = 1;
	for (i=0; i<3; i++)	{
		buff[i][0] = buff[x][0]*b2 + buff[y][0]*b1 + buff[(y+1)][0]*b3;
		buff[i][1] = buff[x][1]*b2 + buff[y][1]*b1 + buff[(y+1)][1]*b3;
		buff[i][2] = buff[x][2]*b2 + buff[y][2]*b1 + buff[(y+1)][2]*b3;
		x++;
		if (i == 0)
			y += 2;
		else
			y++;
	}

	destVertices[0] = buff[0][0]*b2 + buff[1][0]*b1 + buff[2][0]*b3;
	destVertices[1] = buff[0][1]*b2 + buff[1][1]*b1 + buff[2][1]*b3;
	destVertices[2] = buff[0][2]*b2 + buff[1][2]*b1 + buff[2][2]*b3;
	//destVertices->w = 1;

	/**********/
	/* Normal */
	/**********/
	x = 0;
	y = 1;
	for (i=0; i<3; i++)	{
		buff[i][0] = normCoeff[x][0]*b2 + normCoeff[y][0]*b1 + normCoeff[(y+1)][0]*b3;
		buff[i][1] = normCoeff[x][1]*b2 + normCoeff[y][1]*b1 + normCoeff[(y+1)][1]*b3;
		buff[i][2] = normCoeff[x][2]*b2 + normCoeff[y][2]*b1 + normCoeff[(y+1)][2]*b3;
		x++;
		if (i == 0)
			y += 2;
		else
			y++;
	}

	destNormals[0] = buff[0][0]*b2 + buff[1][0]*b1 + buff[2][0]*b3;
	destNormals[1] = buff[0][1]*b2 + buff[1][1]*b1 + buff[2][1]*b3;
	destNormals[2] = buff[0][2]*b2 + buff[1][2]*b1 + buff[2][2]*b3;

	tmpf = 1.0f / sqrt(destNormals[0]*destNormals[0] + destNormals[1]*destNormals[1] + destNormals[2]*destNormals[2]);

	destNormals[0] *= tmpf;
	destNormals[1] *= tmpf;
	destNormals[2] *= tmpf;

}

void TriMesh::Cubic_Bezier_Parameters(FaceHandle fh, double coeff[ ][3], double normCoeff[ ][3]){
	double tmpf;

	Point v1,v2,v3;
	Normal n1,n2,n3;

	FaceVertexIter fv_it=this->fv_iter(fh);

	v1=this->point(fv_it.handle());
	n1=this->normal(fv_it.handle());
	++fv_it;
	v2=this->point(fv_it.handle());
	n2=this->normal(fv_it.handle());
	++fv_it;
	v3=this->point(fv_it.handle());
	n3=this->normal(fv_it.handle());

	Point edge;
	Normal avgNorm;
	Point tangent1;
	Point tangent2;

	/***************************************************/
	/* Compute cubic Bezier control mesh for positions */
	/***************************************************/
	//    v2   0
	//         1 2
	//        3 4 5
	// v1  6 7 8 9  v3

	//3 corner vertices
	coeff[6][0] = v1[0];			coeff[6][1] = v1[1];			coeff[6][2] = v1[2];			
	coeff[0][0] = v2[0];			coeff[0][1] = v2[1];			coeff[0][2] = v2[2];			
	coeff[9][0] = v3[0];			coeff[9][1] = v3[1];			coeff[9][2] = v3[2];

	/***********/
	/* Edge #1 */
	/***********/
	edge = v2 - v1;	
	//E - (E.N)N
	tangent1 = edge;
	tangent1 -= n1 * DOTPROD3(tangent1,n1);
	coeff[3][0] = v1[0] + (tangent1[0] * 0.3333f);
	coeff[3][1] = v1[1] + (tangent1[1] * 0.3333f);
	coeff[3][2] = v1[2] + (tangent1[2] * 0.3333f);

	//E - (E.N)N
	tangent2 = edge;
	tangent2 -= n2 * DOTPROD3(tangent2, n2);
	coeff[1][0] = v2[0] - (tangent2[0] * 0.3333f);
	coeff[1][1] = v2[1] - (tangent2[1] * 0.3333f);
	coeff[1][2] = v2[2] - (tangent2[2] * 0.3333f);

	/***********/
	/* Edge #2 */
	/***********/
	edge = v3 - v2;
	//E - (E.N)N
	tangent1 = edge;
	tangent1 -= n2 * DOTPROD3(tangent1, n2);
	coeff[2][0] = v2[0] + (tangent1[0] * 0.3333f);
	coeff[2][1] = v2[1] + (tangent1[1] * 0.3333f);
	coeff[2][2] = v2[2] + (tangent1[2] * 0.3333f);
	//E - (E.N)N
	tangent2 = edge;
	tangent2 -= n3 * DOTPROD3(tangent2, n3);
	coeff[5][0] = v3[0] - (tangent2[0] * 0.3333f);
	coeff[5][1] = v3[1] - (tangent2[1] * 0.3333f);
	coeff[5][2] = v3[2] - (tangent2[2] * 0.3333f);

	/***********/
	/* Edge #3 */
	/***********/
	edge = v3 - v1;
	//E - (E.N)N
	tangent1 = edge;
	tangent1 -= n1 * DOTPROD3(tangent1, n1);

	coeff[7][0] = v1[0] + (tangent1[0] * 0.3333f);
	coeff[7][1] = v1[1] + (tangent1[1] * 0.3333f);
	coeff[7][2] = v1[2] + (tangent1[2] * 0.3333f);

	//E - (E.N)N
	tangent2 = edge;
	tangent2[0] -= n3[0] * DOTPROD3(tangent2, n3);

	coeff[8][0] = v3[0] - (tangent2[0] * 0.3333f);
	coeff[8][1] = v3[1] - (tangent2[1] * 0.3333f);
	coeff[8][2] = v3[2] - (tangent2[2] * 0.3333f);

	/****************/
	/* Middle Point */
	/****************/
	// Farin's Method
	//    v2  0
	//       1 2
	//      3 4 5
	// v1  6 7 8 9  v3
	coeff[4][0] = (coeff[1][0] + coeff[3][0] + coeff[2][0] + coeff[5][0] + coeff[7][0] + coeff[8][0]) / 4.0f;
	coeff[4][1] = (coeff[1][1] + coeff[3][1] + coeff[2][1] + coeff[5][1] + coeff[7][1] + coeff[8][1]) / 4.0f;
	coeff[4][2] = (coeff[1][2] + coeff[3][2] + coeff[2][2] + coeff[5][2] + coeff[7][2] + coeff[8][2]) / 4.0f;

	coeff[4][0] -= (coeff[0][0] + coeff[6][0] + coeff[9][0]) / 6.0f;
	coeff[4][1] -= (coeff[0][1] + coeff[6][1] + coeff[9][1]) / 6.0f;
	coeff[4][2] -= (coeff[0][2] + coeff[6][2] + coeff[9][2]) / 6.0f;

	/***********************************************************/
	/* Compute control mesh for quadratic normal interpolation */
	/***********************************************************/
	//    v2  0
	//       1 2
	// v1  3 4 5  v3

	//3 corner normals
	normCoeff[3][0] = n1[0];
	normCoeff[3][1] = n1[1];
	normCoeff[3][2] = n1[2];

	normCoeff[0][0] = n2[0];
	normCoeff[0][1] = n2[1];
	normCoeff[0][2] = n2[2];

	normCoeff[5][0] = n3[0];
	normCoeff[5][1] = n3[1];
	normCoeff[5][2] = n3[2];

	/*********************/
	/* Normal for edge 1 */
	/*********************/
	avgNorm = n1 + n2;
	//N - 2*(E.N)E
	edge = v2 - v1;
	tmpf = sqrt((edge[0]*edge[0]) + (edge[1]*edge[1]) + (edge[2]*edge[2]));
	if (tmpf > 0.0f)	{
		edge /= tmpf;

		tmpf = 2.0f * DOTPROD3(edge, avgNorm); //2*(E.N)

		avgNorm -= tmpf*edge;

		tmpf = 1.0f / sqrt((avgNorm[0]*avgNorm[0]) + (avgNorm[1]*avgNorm[1]) + (avgNorm[2]*avgNorm[2]));
		normCoeff[1][0] = avgNorm[0] * tmpf;
		normCoeff[1][1] = avgNorm[1] * tmpf;
		normCoeff[1][2] = avgNorm[2] * tmpf;
	}	else {//Don't reflect normal about plane
		tmpf = 1.0f / sqrt((avgNorm[0]*avgNorm[0]) + (avgNorm[1]*avgNorm[1]) + (avgNorm[2]*avgNorm[2]));
		normCoeff[1][0] = avgNorm[0] * tmpf;
		normCoeff[1][1] = avgNorm[1] * tmpf;
		normCoeff[1][2] = avgNorm[2] * tmpf;
	}

	/*********************/
	/* Normal for edge 2 */
	/*********************/
	avgNorm = n3 + n2;
	//N - 2*(E.N)E
	edge= v3 - v2;
	tmpf = sqrt((edge[0]*edge[0]) + (edge[1]*edge[1]) + (edge[2]*edge[2]));
	if (tmpf > 0.0f)	{
		edge /= tmpf;

		tmpf = 2.0f * DOTPROD3(edge, avgNorm); //2*(E.N)

		avgNorm -= tmpf*edge;

		tmpf = 1.0f / sqrt((avgNorm[0]*avgNorm[0]) + (avgNorm[1]*avgNorm[1]) + (avgNorm[2]*avgNorm[2]));
		normCoeff[2][0] = avgNorm[0] * tmpf;
		normCoeff[2][1] = avgNorm[1] * tmpf;
		normCoeff[2][2] = avgNorm[2] * tmpf;
	}	else {//Don't reflect normal about plane
		tmpf = 1.0f / sqrt((avgNorm[0]*avgNorm[0]) + (avgNorm[1]*avgNorm[1]) + (avgNorm[2]*avgNorm[2]));
		normCoeff[2][0] = avgNorm[0] * tmpf;
		normCoeff[2][1] = avgNorm[1] * tmpf;
		normCoeff[2][2] = avgNorm[2] * tmpf;
	}

	/*********************/
	/* Normal for edge 3 */
	/*********************/
	avgNorm= n1 + n3;
	//N - 2*(E.N)E
	edge = v1 - v3;
	tmpf = sqrt((edge[0]*edge[0]) + (edge[1]*edge[1]) + (edge[2]*edge[2]));
	if (tmpf > 0.0f)	{
		edge /= tmpf;

		tmpf = 2.0f * ((edge[0]*avgNorm[0]) + (edge[1]*avgNorm[1]) + (edge[2]*avgNorm[2])); //2*(E.N)

		avgNorm -= tmpf*edge;
		tmpf = 1.0f / sqrt((avgNorm[0]*avgNorm[0]) + (avgNorm[1]*avgNorm[1]) + (avgNorm[2]*avgNorm[2]));
		normCoeff[4][0] = avgNorm[0] * tmpf;
		normCoeff[4][1] = avgNorm[1] * tmpf;
		normCoeff[4][2] = avgNorm[2] * tmpf;
	}	else {//Don't reflect normal about plane
		tmpf = 1.0f / sqrt((avgNorm[0]*avgNorm[0]) + (avgNorm[1]*avgNorm[1]) + (avgNorm[2]*avgNorm[2]));
		normCoeff[4][0] = avgNorm[0] * tmpf;
		normCoeff[4][1] = avgNorm[1] * tmpf;
		normCoeff[4][2] = avgNorm[2] * tmpf;
	}
}

void TriMesh::para_to_point(FaceHandle fh,double alpha,double beta,double gamma, Point &p, Normal &n){
	Point destVertices;
	Normal destNormals;

	double tmpf;
	double coeff[10][3];
	double normCoeff[6][3];

	Point v1,v2,v3;
	Normal n1,n2,n3;

	FaceVertexIter fv_it=this->fv_iter(fh);

	v1=this->point(fv_it.handle());
	n1=this->normal(fv_it.handle());
	++fv_it;
	v2=this->point(fv_it.handle());
	n2=this->normal(fv_it.handle());
	++fv_it;
	v3=this->point(fv_it.handle());
	n3=this->normal(fv_it.handle());

	Point edge;
	Normal avgNorm;
	Point tangent1;
	Point tangent2;

	/***************************************************/
	/* Compute cubic Bezier control mesh for positions */
	/***************************************************/
	//    v2   0
	//         1 2
	//        3 4 5
	// v1  6 7 8 9  v3

	//3 corner vertices
	coeff[6][0] = v1[0];			coeff[6][1] = v1[1];			coeff[6][2] = v1[2];			
	coeff[0][0] = v2[0];			coeff[0][1] = v2[1];			coeff[0][2] = v2[2];			
	coeff[9][0] = v3[0];			coeff[9][1] = v3[1];			coeff[9][2] = v3[2];

	/***********/
	/* Edge #1 */
	/***********/
	edge = v2 - v1;	
	//E - (E.N)N
	tangent1 = edge;
	tangent1 -= n1 * DOTPROD3(tangent1,n1);
	coeff[3][0] = v1[0] + (tangent1[0] * 0.3333f);
	coeff[3][1] = v1[1] + (tangent1[1] * 0.3333f);
	coeff[3][2] = v1[2] + (tangent1[2] * 0.3333f);

	//E - (E.N)N
	tangent2 = edge;
	tangent2 -= n2 * DOTPROD3(tangent2, n2);
	coeff[1][0] = v2[0] - (tangent2[0] * 0.3333f);
	coeff[1][1] = v2[1] - (tangent2[1] * 0.3333f);
	coeff[1][2] = v2[2] - (tangent2[2] * 0.3333f);

	/***********/
	/* Edge #2 */
	/***********/
	edge = v3 - v2;
	//E - (E.N)N
	tangent1 = edge;
	tangent1 -= n2 * DOTPROD3(tangent1, n2);
	coeff[2][0] = v2[0] + (tangent1[0] * 0.3333f);
	coeff[2][1] = v2[1] + (tangent1[1] * 0.3333f);
	coeff[2][2] = v2[2] + (tangent1[2] * 0.3333f);
	//E - (E.N)N
	tangent2 = edge;
	tangent2 -= n3 * DOTPROD3(tangent2, n3);
	coeff[5][0] = v3[0] - (tangent2[0] * 0.3333f);
	coeff[5][1] = v3[1] - (tangent2[1] * 0.3333f);
	coeff[5][2] = v3[2] - (tangent2[2] * 0.3333f);

	/***********/
	/* Edge #3 */
	/***********/
	edge = v3 - v1;
	//E - (E.N)N
	tangent1 = edge;
	tangent1 -= n1 * DOTPROD3(tangent1, n1);

	coeff[7][0] = v1[0] + (tangent1[0] * 0.3333f);
	coeff[7][1] = v1[1] + (tangent1[1] * 0.3333f);
	coeff[7][2] = v1[2] + (tangent1[2] * 0.3333f);

	//E - (E.N)N
	tangent2 = edge;
	tangent2[0] -= n3[0] * DOTPROD3(tangent2, n3);

	coeff[8][0] = v3[0] - (tangent2[0] * 0.3333f);
	coeff[8][1] = v3[1] - (tangent2[1] * 0.3333f);
	coeff[8][2] = v3[2] - (tangent2[2] * 0.3333f);

	/****************/
	/* Middle Point */
	/****************/
	// Farin's Method
	//    v2  0
	//       1 2
	//      3 4 5
	// v1  6 7 8 9  v3
	coeff[4][0] = (coeff[1][0] + coeff[3][0] + coeff[2][0] + coeff[5][0] + coeff[7][0] + coeff[8][0]) / 4.0f;
	coeff[4][1] = (coeff[1][1] + coeff[3][1] + coeff[2][1] + coeff[5][1] + coeff[7][1] + coeff[8][1]) / 4.0f;
	coeff[4][2] = (coeff[1][2] + coeff[3][2] + coeff[2][2] + coeff[5][2] + coeff[7][2] + coeff[8][2]) / 4.0f;

	coeff[4][0] -= (coeff[0][0] + coeff[6][0] + coeff[9][0]) / 6.0f;
	coeff[4][1] -= (coeff[0][1] + coeff[6][1] + coeff[9][1]) / 6.0f;
	coeff[4][2] -= (coeff[0][2] + coeff[6][2] + coeff[9][2]) / 6.0f;

	/***********************************************************/
	/* Compute control mesh for quadratic normal interpolation */
	/***********************************************************/
	//    v2  0
	//       1 2
	// v1  3 4 5  v3

	//3 corner normals
	normCoeff[3][0] = n1[0];
	normCoeff[3][1] = n1[1];
	normCoeff[3][2] = n1[2];

	normCoeff[0][0] = n2[0];
	normCoeff[0][1] = n2[1];
	normCoeff[0][2] = n2[2];

	normCoeff[5][0] = n3[0];
	normCoeff[5][1] = n3[1];
	normCoeff[5][2] = n3[2];

	/*********************/
	/* Normal for edge 1 */
	/*********************/
	avgNorm = n1 + n2;
	//N - 2*(E.N)E
	edge = v2 - v1;
	tmpf = sqrt((edge[0]*edge[0]) + (edge[1]*edge[1]) + (edge[2]*edge[2]));
	if (tmpf > 0.0f)	{
		edge /= tmpf;

		tmpf = 2.0f * DOTPROD3(edge, avgNorm); //2*(E.N)

		avgNorm -= tmpf*edge;

		tmpf = 1.0f / sqrt((avgNorm[0]*avgNorm[0]) + (avgNorm[1]*avgNorm[1]) + (avgNorm[2]*avgNorm[2]));
		normCoeff[1][0] = avgNorm[0] * tmpf;
		normCoeff[1][1] = avgNorm[1] * tmpf;
		normCoeff[1][2] = avgNorm[2] * tmpf;
	}	else {//Don't reflect normal about plane
		tmpf = 1.0f / sqrt((avgNorm[0]*avgNorm[0]) + (avgNorm[1]*avgNorm[1]) + (avgNorm[2]*avgNorm[2]));
		normCoeff[1][0] = avgNorm[0] * tmpf;
		normCoeff[1][1] = avgNorm[1] * tmpf;
		normCoeff[1][2] = avgNorm[2] * tmpf;
	}

	/*********************/
	/* Normal for edge 2 */
	/*********************/
	avgNorm = n3 + n2;
	//N - 2*(E.N)E
	edge= v3 - v2;
	tmpf = sqrt((edge[0]*edge[0]) + (edge[1]*edge[1]) + (edge[2]*edge[2]));
	if (tmpf > 0.0f)	{
		edge /= tmpf;

		tmpf = 2.0f * DOTPROD3(edge, avgNorm); //2*(E.N)

		avgNorm -= tmpf*edge;

		tmpf = 1.0f / sqrt((avgNorm[0]*avgNorm[0]) + (avgNorm[1]*avgNorm[1]) + (avgNorm[2]*avgNorm[2]));
		normCoeff[2][0] = avgNorm[0] * tmpf;
		normCoeff[2][1] = avgNorm[1] * tmpf;
		normCoeff[2][2] = avgNorm[2] * tmpf;
	}	else {//Don't reflect normal about plane
		tmpf = 1.0f / sqrt((avgNorm[0]*avgNorm[0]) + (avgNorm[1]*avgNorm[1]) + (avgNorm[2]*avgNorm[2]));
		normCoeff[2][0] = avgNorm[0] * tmpf;
		normCoeff[2][1] = avgNorm[1] * tmpf;
		normCoeff[2][2] = avgNorm[2] * tmpf;
	}

	/*********************/
	/* Normal for edge 3 */
	/*********************/
	avgNorm= n1 + n3;
	//N - 2*(E.N)E
	edge = v1 - v3;
	tmpf = sqrt((edge[0]*edge[0]) + (edge[1]*edge[1]) + (edge[2]*edge[2]));
	if (tmpf > 0.0f)	{
		edge /= tmpf;

		tmpf = 2.0f * ((edge[0]*avgNorm[0]) + (edge[1]*avgNorm[1]) + (edge[2]*avgNorm[2])); //2*(E.N)

		avgNorm -= tmpf*edge;
		tmpf = 1.0f / sqrt((avgNorm[0]*avgNorm[0]) + (avgNorm[1]*avgNorm[1]) + (avgNorm[2]*avgNorm[2]));
		normCoeff[4][0] = avgNorm[0] * tmpf;
		normCoeff[4][1] = avgNorm[1] * tmpf;
		normCoeff[4][2] = avgNorm[2] * tmpf;
	}	else {//Don't reflect normal about plane
		tmpf = 1.0f / sqrt((avgNorm[0]*avgNorm[0]) + (avgNorm[1]*avgNorm[1]) + (avgNorm[2]*avgNorm[2]));
		normCoeff[4][0] = avgNorm[0] * tmpf;
		normCoeff[4][1] = avgNorm[1] * tmpf;
		normCoeff[4][2] = avgNorm[2] * tmpf;
	}

	/************************/  //   b2
	/* Compute new vertices */  //
	/************************/  // b1  b3

	EvaluateBezierPoint(coeff, normCoeff, alpha, beta, gamma, p, n);
}

void TriMesh::para_to_point_PN_triangle(FaceHandle fh, double u,double v,double w, Point &p, Normal &n){
	FaceVertexIter fv_it=fv_iter(fh);
	Point b300=point(fv_it.handle());
	Normal N1=normal(fv_it.handle());
	++fv_it;
	Point b030=point(fv_it.handle());
	Normal N2=normal(fv_it.handle());
	++fv_it;
	Point b003=point(fv_it.handle());
	Normal N3=normal(fv_it.handle());

	Point P1=b300; Point P2=b030; Point P3=b003;

	double w12=DOTPROD3(P2-P1, N1);
	double w21=DOTPROD3(P1-P2, N2);
	double w13=DOTPROD3(P3-P1, N1);
	double w31=DOTPROD3(P1-P3, N3);
	double w23=DOTPROD3(P3-P2, N2);
	double w32=DOTPROD3(P2-P3, N3);

	Point b210=(2.0*P1+P2-w12*N1)/3.0;
	Point b120=(2.0*P2+P1-w21*N2)/3.0;
	Point b021=(2.0*P2+P3-w23*N2)/3.0;
	Point b012=(2.0*P3+P2-w32*N3)/3.0;
	Point b102=(2.0*P3+P1-w31*N3)/3.0;
	Point b201=(2.0*P1+P3-w13*N1)/3.0;
	Point E=(b210+b120+b012+b021+b102+b201)/6.0;
	Point V=(P1+P2+P3)/3.0;
	Point b111=E+(E-V)/2.0;

	p = b300*w*w*w+b030*u*u*u+b003*v*v*v
		+3.0*b120*w*u*u+3.0*b210*w*w*u
		+3.0*b201*w*w*v+3.0*b102*w*v*v
		+3.0*b012*u*v*v+3.0*b021*u*u*v
		+6.0*b111*w*u*v;

	double v12=2.0*DOTPROD3(P2-P1, N1+N2)/DOTPROD3(P2-P1, P2-P1);
	double v23=2.0*DOTPROD3(P3-P2, N2+N3)/DOTPROD3(P3-P2, P3-P2);
	double v31=2.0*DOTPROD3(P1-P3, N3+N1)/DOTPROD3(P1-P3, P1-P3);

	Normal n110=N1+N2-v12*(P2-P1); GeometryTools::normalize(n110);
	Normal n011=N2+N3-v23*(P3-P2); GeometryTools::normalize(n011);
	Normal n101=N3+N1-v31*(P1-P3); GeometryTools::normalize(n101);

	n=N1*w*w+N2*u*u+N3*v*v+n110*w*u+n011*u*v+n101*w*v;
	GeometryTools::normalize(n);
}

TriMesh* TriMesh::read(const char* filename, OpenMesh::IO::Options* opt)
{
	TriMesh* mesh = new TriMesh;

	OpenMesh::IO::Options default_opt;
	if (!opt)
		opt = &default_opt;
	//opt+= OpenMesh::IO::Options::VertexNormal;

	if (!OpenMesh::IO::read_mesh(*mesh, filename, *opt)) {
		delete mesh;
		return NULL;
	}
	// If the file did not provide vertex normals, then calculate them
	if ( !opt->check(OpenMesh::IO::Options::VertexNormal) &&
		mesh->has_face_normals() && mesh->has_vertex_normals() ) {
			// let the mesh update the normals
			mesh->update_normals();
	}

	if (!mesh->has_vertex_colors())	{
		mesh->request_vertex_colors();
		mesh->set_vertex_color(255,255,153);
	}


	//const Color yellow(255,255,153);
	//const Color  blue(153, 153, 255);

	TriMesh::ConstFaceIter fit(mesh->faces_begin()), fEnd(mesh->faces_end());
	mesh->vertex_group0.clear();
	mesh->vertex_group1.clear();
	mesh->vertex_group2.clear();
	TriMesh::ConstFaceVertexIter fvit = mesh->cfv_iter(fit.handle());
	for(fvit = mesh->cfv_iter(fit.handle()); fit!=fEnd; ++fit)	{
		fvit = mesh->cfv_iter(fit.handle());
		mesh->vertex_group0.push_back(fvit.handle().idx());	
		mesh->vertex_group2.push_back(fvit.handle().idx());    ++fvit;
		mesh->vertex_group0.push_back(fvit.handle().idx());
		mesh->vertex_group2.push_back(fvit.handle().idx());     ++fvit;
		mesh->vertex_group0.push_back(fvit.handle().idx());
		mesh->vertex_group2.push_back(fvit.handle().idx());
	}


	mesh->average_edge_length=mesh->getAverageEdgeLength();

	mesh->init_boundary();


	return mesh;
}

bool TriMesh::save(const char* filename,  TriMesh* mesh, OpenMesh::IO::Options* opt)
{
	OpenMesh::IO::Options default_opt;
	if (!opt)
		opt = &default_opt;

	if (!OpenMesh::IO::write_mesh(*mesh, filename, *opt)) 
		return false;
	return true;
}

void TriMesh::getMeshStat(std::vector<double> &meshStat, TriMesh* meshCopy)
{
	// statistics 1: percentage of vertices having 6 neighbors. [0]
	meshStat.clear();
	TriMesh::VertexIter v_it;
	TriMesh::VertexVertexIter vv_it;
	int valence;
	int n_6n = 0; // number of vertices having 6 neighbors
	for(v_it=vertices_begin(); v_it!=vertices_end(); ++v_it)
	{
		valence = 0;
		for(vv_it=vv_iter(v_it); vv_it; ++vv_it)
			valence++;
		if(valence == 6)
			n_6n++;
	}
	meshStat.push_back(1.0*n_6n/n_vertices());
	// statistics 2: radius ratio, 2*r/R. [1] Tmin, [2] Tmean
	TriMesh::FaceIter f_it;
	std::vector<double> radiusRatio;
	/*for(f_it=faces_begin(); f_it!=faces_end(); ++f_it)
	{
		radiusRatio.push_back(2 * GeometryTools::radiusOfInscribedCircle(this, f_it.handle()) 
			/ GeometryTools::radiusOfCircumscribedCircle(this, f_it.handle()));
	}*/
	meshStat.push_back(*std::min_element(radiusRatio.begin(), radiusRatio.end()));
	meshStat.push_back(std::accumulate(radiusRatio.begin(), radiusRatio.end(), 0.0)/n_faces());
	radiusRatio.clear();
	// statistics 3: min angles in triangles, [3] minimum degree, [4] mean degree
	std::vector<double> angle;
	//for(f_it=faces_begin(); f_it!=faces_end(); ++f_it)
	//	angle.push_back(GeometryTools::minAngleInTriangle(this, f_it.handle()));
	meshStat.push_back(*std::min_element(angle.begin(), angle.end()) * 180.0 / 3.1415926);
	meshStat.push_back(std::accumulate(angle.begin(), angle.end(), 0.0)/n_faces() * 180.0 / 3.1415926);
	angle.clear();
	// statistics 4: relative RMS of edge lengths, [5] edge RRMS
	TriMesh::HalfedgeIter he_it, he_it_copy;
	double edgeRMS = 0;
	int num = 0;
	for(he_it=halfedges_begin(), he_it_copy=meshCopy->halfedges_begin(); he_it!=halfedges_end(); ++he_it, ++he_it_copy)
	{
		num++;
		TriMesh::Point p1_copy = meshCopy->point(meshCopy->from_vertex_handle(he_it_copy.handle()));
		TriMesh::Point p2_copy = meshCopy->point(meshCopy->to_vertex_handle(he_it_copy.handle()));
		TriMesh::Point p1 = point(from_vertex_handle(he_it.handle()));
		TriMesh::Point p2 = point(to_vertex_handle(he_it.handle()));
		edgeRMS += pow(abs((p1_copy-p2_copy).length() - (p1-p2).length())/(p1_copy-p2_copy).length(), 2);
	}
	edgeRMS = sqrt(edgeRMS / num);
	meshStat.push_back(edgeRMS);
	// statistics 5: relative RMS of volume magnitudes, [6]
	double volumeRMS = 0;
	double volume = 0;
	double volumeCopy = 0;
	TriMesh::FaceIter f_it_copy;
	TriMesh::FaceVertexIter fv_it, fv_it_copy;
	for(f_it=faces_begin(), f_it_copy=meshCopy->faces_begin(); f_it!=faces_end(); ++f_it, ++f_it_copy)
	{
		fv_it=fv_iter(f_it);
		TriMesh::Point p1 = point(fv_it.handle());
		++fv_it;
		TriMesh::Point p2 = point(fv_it.handle());
		++fv_it;
		TriMesh::Point p3 = point(fv_it.handle());
		volume += p3 | (p1%p2);
		fv_it_copy=meshCopy->fv_iter(f_it_copy);
		p1 = meshCopy->point(fv_it_copy.handle());
		++fv_it_copy;
		p2 = meshCopy->point(fv_it_copy.handle());
		++fv_it_copy;
		p3 = meshCopy->point(fv_it_copy.handle());
		volumeCopy += p3 | (p1%p2);
	}
	volumeRMS = abs(volume - volumeCopy)/volume;
	meshStat.push_back(volumeRMS);
}

TriMesh::VertexHandle TriMesh::approximatePointToVertex(Point center){
	TriMesh::VertexIter v_it,v_end(this->vertices_end()),center_it;
	double distance=99999999.9;
	Point current;
	for(v_it=this->vertices_begin();v_it!=v_end;++v_it){
		current=this->point(v_it);
		if(GeometryTools::distPoint2Point(current, center)<distance)
		{
			distance=GeometryTools::distPoint2Point(current,center);
			center_it=v_it;
		}
	}
	return center_it.handle();
}

bool TriMesh::selectInitialPoint(Point center){
	//start_vertex.push_back(approximatePointToVertex(center));
	//if(start_vertex.size()>0)
	//	start_vertex[0]=approximatePointToVertex(center);

	//if(is_valid_handle(start_vertex[start_vertex.size()-1]))
	//	return true;
	//else 
		return false;
}

void TriMesh::init_boundary(){
	boundary.clear();
	TriMesh::HalfedgeIter he_it=halfedges_begin();
	for ( ; he_it!=halfedges_end(); ++he_it)	{
		if (abs(point(to_vertex_handle(he_it))[2]-1.0)<eps_)
			continue;
		if (is_boundary(he_it.handle()))
			break;
	}
	if(he_it==halfedges_end())
		return;
	TriMesh::HalfedgeHandle cur_he=he_it.handle();
	boundary.push_back(to_vertex_handle(cur_he).idx());
	cur_he=next_halfedge_handle(cur_he);
	while (cur_he!=he_it.handle())	{
		boundary.push_back(to_vertex_handle(cur_he).idx());
		cur_he=next_halfedge_handle(cur_he);
	}
}