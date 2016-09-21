#include "GeometryTools.h"

double EPS=1e-14;

GeometryTools::GeometryTools(void)
{
}

GeometryTools::~GeometryTools(void)
{
}

void GeometryTools::normalize(Point& p){
	double len=sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
	if( len<EPS )
		return;
	p[0] = p[0]/len;
	p[1] = p[1]/len;
	p[2] = p[2]/len;
}

double GeometryTools::length2(Point p)
{
	return p[0]*p[0]+p[1]*p[1]+p[2]*p[2];
}

double GeometryTools::length(Point p)
{
	return sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
}

 double GeometryTools::dotProduct(Point p1, Point p2)
 {
	 return p1[0]*p2[0]+p1[1]*p2[1]+p1[2]*p2[2];
 }

 Point GeometryTools::crossProduct(Point p1, Point p2)
 {
	 return Point(p1[1]*p2[2] - p1[2]*p2[1], p1[2]*p2[0] - p1[0]*p2[2], p1[0]*p2[1] - p1[1]*p2[0]);
 }

bool GeometryTools::isInsideOfTriangle(Point p0, Point p1, Point p2, float p[])
{
	double v[2] = {p[0], p[1]};
	double v0[2] = {p0[0], p0[1]};
	double v1[2] = {p1[0]-p0[0], p1[1]-p0[1]};
	double v2[2] = {p2[0]-p0[0], p2[1]-p0[1]};
	double a = ((v[0]*v2[1]-v[1]*v2[0])-(v0[0]*v2[1]-v0[1]*v2[0]))/(v1[0]*v2[1]-v1[1]*v2[0]);
	double b = -((v[0]*v1[1]-v[1]*v1[0])-(v0[0]*v1[1]-v0[1]*v1[0]))/(v1[0]*v2[1]-v1[1]*v2[0]);
	return (a>0&&b>0&&(a+b)<1);
}

bool GeometryTools::isInsideOfTriangle(Point p0, Point p1, Point p2, Point p)
{
	double v[2] = {p[0], p[1]};
	double v0[2] = {p0[0], p0[1]};
	double v1[2] = {p1[0]-p0[0], p1[1]-p0[1]};
	double v2[2] = {p2[0]-p0[0], p2[1]-p0[1]};
	double a = ((v[0]*v2[1]-v[1]*v2[0])-(v0[0]*v2[1]-v0[1]*v2[0]))/(v1[0]*v2[1]-v1[1]*v2[0]);
	double b = -((v[0]*v1[1]-v[1]*v1[0])-(v0[0]*v1[1]-v0[1]*v1[0]))/(v1[0]*v2[1]-v1[1]*v2[0]);
	return (a>0&&b>0&&(a+b)<1);
}

double GeometryTools::distPoint2Point(Point p1, Point p2)
{
	return sqrt((p2[0]-p1[0])*(p2[0]-p1[0])+(p2[1]-p1[1])*(p2[1]-p1[1])+(p2[2]-p1[2])*(p2[2]-p1[2]));
}

double GeometryTools::squareDistPoint2Point(Point p1, Point p2)
{
	return (p2[0]-p1[0])*(p2[0]-p1[0])+(p2[1]-p1[1])*(p2[1]-p1[1])+(p2[2]-p1[2])*(p2[2]-p1[2]);
}

double GeometryTools::distPoint2Line(Point p1, Point p2, float p[])
{
	return abs((p2[0]-p1[0])*(p1[1]-p[1]) - (p1[0]-p[0])*(p2[1]-p1[1]))/sqrt((p2[0]-p1[0])*(p2[0]-p1[0])+(p2[1]-p1[1])*(p2[1]-p1[1]));
}
double GeometryTools::distPoint2Line(Point p1, Point p2, Point p){
	return abs((p2[0]-p1[0])*(p1[1]-p[1]) - (p1[0]-p[0])*(p2[1]-p1[1]))/sqrt((p2[0]-p1[0])*(p2[0]-p1[0])+(p2[1]-p1[1])*(p2[1]-p1[1]));
}

double GeometryTools::distPoint2Plane(double plane[4], Point p)
{
	return (plane[0]*p[0]+plane[1]*p[1]+plane[2]*p[2]+plane[3]) / sqrt(plane[0]*plane[0]+plane[1]*plane[1]+plane[2]*plane[2]);
}

void GeometryTools::spanPlane(Point p1, Point p2, Point p3, double result[4])
{
	Point plane = ((p2-p1)%(p3-p1)) / ((p2-p1).length() * (p3-p1).length());
	result[0] = plane[0];
	result[1] = plane[1];
	result[2] = plane[2];
	result[3] = -(result[0]*p1[0] + result[1]*p1[1] + result[2]*p1[2]);
}

void GeometryTools::projectPoint2Plane(double plane[4], Point p1, Point &p)
{
	double denom = plane[0]*plane[0] + plane[1]*plane[1] + plane[2]*plane[2];
	double mu = -(plane[3]+plane[0]*p1[0]+plane[1]*p1[1]+plane[2]*p1[2]) / denom;
	p[0] = p1[0] + mu*plane[0];
	p[1] = p1[1] + mu*plane[1];
	p[2] = p1[2] + mu*plane[2];
}

void GeometryTools::project(Point &p, Point normal){
	double dotprod=GeometryTools::dotProduct(p,normal);
	p[0]=p[0]-dotprod*normal[0];
	p[1]=p[1]-dotprod*normal[1];
	p[2]=p[2]-dotprod*normal[2];	
	GeometryTools::normalize(p);
}

int GeometryTools::isInsideOfPolygon(int nvert, float *vertx, float *verty, float testx, float testy)
{
	int i, j, c = 0;
	for (i = 0, j = nvert-1; i < nvert; j = i++) {
		if ( ((verty[i]>testy) != (verty[j]>testy)) &&
			(testx < (vertx[j]-vertx[i]) * (testy-verty[i]) / (verty[j]-verty[i]) + vertx[i]) )
			c = !c;
	}
	return c;
}

//double GeometryTools::radiusOfInscribedCircle(TriMesh *mesh, FaceHandle &fh)
//{
//	FaceVertexIter fv_it = mesh->fv_iter(fh);
//	Point p1 = mesh->point(fv_it);
//	++fv_it;
//	Point p2 = mesh->point(fv_it);
//	++fv_it;
//	Point p3 = mesh->point(fv_it);
//	double a = (p1-p2).length();
//	double b = (p2-p3).length();
//	double c = (p3-p1).length();
//	double s = (a+b+c)/2;
//	return sqrt((s-a)*(s-b)*(s-c)/s);
//}

//double GeometryTools:: radiusOfCircumscribedCircle(TriMesh *mesh, FaceHandle &fh)
//{
//	FaceVertexIter fv_it = mesh->fv_iter(fh);
//	Point p1 = mesh->point(fv_it);
//	++fv_it;
//	Point p2 = mesh->point(fv_it);
//	++fv_it;
//	Point p3 = mesh->point(fv_it);
//	double a = (p1-p2).length();
//	double b = (p2-p3).length();
//	double c = (p3-p1).length();
//	return a*b*c/sqrt((a+b+c)*(-a+b+c)*(a-b+c)*(a+b-c));
//}

//double GeometryTools::minAngleInTriangle(TriMesh *mesh, FaceHandle &fh)
//{
//	FaceVertexIter fv_it = mesh->fv_iter(fh);
//	Point p1 = mesh->point(fv_it);
//	++fv_it;
//	Point p2 = mesh->point(fv_it);
//	++fv_it;
//	Point p3 = mesh->point(fv_it);
//	double minAngle = 0;
//	double angle1 = acos((p1-p2)|(p3-p2) / ((p1-p2).length()*(p3-p2).length()));
//	double angle2 = acos((p2-p1)|(p3-p1) / ((p2-p1).length()*(p3-p1).length()));
//	double angle3 = acos((p1-p3)|(p2-p3) / ((p1-p3).length()*(p2-p3).length()));
//	minAngle = angle1<angle2 ? angle1:angle2;
//	minAngle = minAngle<angle3 ? minAngle:angle3;
//	return minAngle;
//}

bool GeometryTools::intersect_segment_to_segment(Point A, Point B, Point C, Point D, Point& p){
	//double demonminator= (A[0]-B[0])*(C[1]-D[1]) - (C[0]-D[0])*(A[1]-B[1]);
	if (abs( (A[0]-B[0])*(C[1]-D[1]) - (C[0]-D[0])*(A[1]-B[1]) )<EPS)
	    return false;
	p[0] = ( (A[0]*B[1]-B[0]*A[1])*(C[0]-D[0]) - (C[0]*D[1]-D[0]*C[1])*(A[0]-B[0]) )
		       / ( (A[0]-B[0])*(C[1]-D[1]) - (C[0]-D[0])*(A[1]-B[1]) );
	p[1] = ( (A[0]*B[1]-B[0]*A[1])*(C[1]-D[1]) - (C[0]*D[1]-D[0]*C[1])*(A[1]-B[1]) )
		      / ( (A[0]-B[0])*(C[1]-D[1]) - (C[0]-D[0])*(A[1]-B[1]) );
	p[2] = 0.0;

	if (distPoint2Point(p,A)+distPoint2Point(p,B) > distPoint2Point(A, B)+EPS)
	    return false;
	else if (distPoint2Point(p,C)+distPoint2Point(p,D) > distPoint2Point(C, D)+EPS)
		return false;
	else 
		return true;
}

bool GeometryTools::is_same_side(Point p1, Point p2, Point t1, Point t2){
	if ( (p2[0]-p1[0])*(t1[1]-p1[1]) > (p2[1]-p1[1])*(t1[0]-p1[0]) ){
		if ( (p2[0]-p1[0])*(t2[1]-p1[1]) - (p2[1]-p1[1])*(t2[0]-p1[0]) > -EPS)
			return true;
		else
			return false;
	} else{
		if (  (p2[1]-p1[1])*(t2[0]-p1[0]) - (p2[0]-p1[0])*(t2[1]-p1[1])  > -EPS)
			return true;
		else
			return false;
	}
}

bool GeometryTools::barycentrc_coord(Point curgeo, Point A, Point B, Point C, Point& coord){
	Point temp;
	if (!is_same_side(B, C, A, curgeo))
	{
		intersect_segment_to_segment(B, C, A, curgeo, temp);
		curgeo=temp;
	}
	if (!(A[0]==C[0]&&B[0]==C[0]))
	{
		if ((B[0]-C[0])*(A[1]-C[1]+A[2]-C[2])==(A[0]-C[0])*(B[1]-C[1]+B[2]-C[2]))
		{
			double ab=GeometryTools::distPoint2Point(A,B);
			double bc=GeometryTools::distPoint2Point(B,C);
			double ac=GeometryTools::distPoint2Point(C,A);

			double a=GeometryTools::distPoint2Point(B,C);
			double b=GeometryTools::distPoint2Point(B, curgeo);
			double c=GeometryTools::distPoint2Point(C,curgeo);

			double aa=GeometryTools::distPoint2Point(A,C);
			double bb=GeometryTools::distPoint2Point(A, curgeo);
			double cc=GeometryTools::distPoint2Point(C,curgeo);
			if( (ab+bc>ac)&&(ab+ac>bc)&&(ac+bc>ab)&&(a+b>c)&&(b+c>a)&&(a+c>b)&&(aa+bb>cc)&&(bb+cc>aa)&&(cc+aa>bb))
			{
				coord[0]=sqrt((a+b+c)/2.0*((a+b+c)/2.0-a)*((a+b+c)/2.0-b)*((a+b+c)/2.0-c))/sqrt((ab+bc+ac)/2.0*((ab+bc+ac)/2.0-ab)*((ab+bc+ac)/2.0-bc)*((ab+bc+ac)/2.0-ac));
				coord[1]=sqrt((aa+bb+cc)/2.0*((aa+bb+cc)/2.0-aa)*((aa+bb+cc)/2.0-bb)*((aa+bb+cc)/2.0-cc))/sqrt((ab+bc+ac)/2.0*((ab+bc+ac)/2.0-ab)*((ab+bc+ac)/2.0-bc)*((ab+bc+ac)/2.0-ac));
				coord[2]=1.0-coord[0]-coord[1];
				return true;
			}else{
				//not a triangle, there is two edge length sum equal to the third edge length
				//CString thes;
				//thes.Format("outer triangle area is zero when computing barycentric coord");
				//AfxMessageBox(thes, MB_OK|MB_ICONINFORMATION);
				return false;
			}
		}else{
			if (((B[0]-C[0])*(C[1]-curgeo[1]+C[2]-curgeo[2])-(C[0]-curgeo[0])*(B[1]-C[1]+B[2]-C[2]))==0)
				coord[0]=0.0;
			else
				coord[0]=((B[0]-C[0])*(C[1]-curgeo[1]+C[2]-curgeo[2])-(C[0]-curgeo[0])*(B[1]-C[1]+B[2]-C[2]))
				/((A[0]-C[0])*(B[1]-C[1]+B[2]-C[2])-(B[0]-C[0])*(A[1]-C[1]+A[2]-C[2]));
			if (((A[0]-C[0])*(C[1]-curgeo[1]+C[2]-curgeo[2])-(C[0]-curgeo[0])*(A[1]-C[1]+A[2]-C[2]))==0)
				coord[1]=0.0;
			else
				coord[1]=((A[0]-C[0])*(C[1]-curgeo[1]+C[2]-curgeo[2])-(C[0]-curgeo[0])*(A[1]-C[1]+A[2]-C[2]))
				/((B[0]-C[0])*(A[1]-C[1]+A[2]-C[2])-(A[0]-C[0])*(B[1]-C[1]+B[2]-C[2]));
			coord[2]=1.0-coord[0]-coord[1];
			return true;
		}
	}
	else{
		Point ba,bc,bd,cb,ca,cd;
		ba=A-B;
		bc=C-B;
		bd=curgeo-B;
		cb=B-C;
		ca=A-C;
		cd=curgeo-C;

		double test0=GeometryTools::squareDistPoint2Point(A,B)*GeometryTools::squareDistPoint2Point(C,B) - 
			GeometryTools::dotProduct(ba,bc)*GeometryTools::dotProduct(ba,bc);
		double test1=(GeometryTools::squareDistPoint2Point(curgeo,B)*GeometryTools::squareDistPoint2Point(C,B) - 
			GeometryTools::dotProduct(bd,bc)*GeometryTools::dotProduct(bd,bc));

		if ((GeometryTools::squareDistPoint2Point(curgeo,B)*GeometryTools::squareDistPoint2Point(C,B) - 
			GeometryTools::dotProduct(bd,bc)*GeometryTools::dotProduct(bd,bc))<EPS)
			coord[0]=0.0;
		else
			coord[0]=sqrt( GeometryTools::squareDistPoint2Point(curgeo,B)*GeometryTools::squareDistPoint2Point(C,B) - 
			GeometryTools::dotProduct(bd,bc)*GeometryTools::dotProduct(bd,bc) )
			/sqrt( GeometryTools::squareDistPoint2Point(A,B)*GeometryTools::squareDistPoint2Point(C,B) - 
			GeometryTools::dotProduct(ba,bc)*GeometryTools::dotProduct(ba,bc) );

		if ((GeometryTools::squareDistPoint2Point(curgeo,C)*GeometryTools::squareDistPoint2Point(C,A) - 
			GeometryTools::dotProduct(cd,ca)*GeometryTools::dotProduct(cd,ca))<EPS)
			coord[1]=0.0;
		else
			coord[1]=sqrt( GeometryTools::squareDistPoint2Point(curgeo,C)*GeometryTools::squareDistPoint2Point(C,A) - 
			GeometryTools::dotProduct(cd,ca)*GeometryTools::dotProduct(cd,ca))
			/sqrt( GeometryTools::squareDistPoint2Point(C,B)*GeometryTools::squareDistPoint2Point(C,A) - 
			GeometryTools::dotProduct(cb,ca)*GeometryTools::dotProduct(cb,ca) );
		coord[2]=1.0-coord[0]-coord[1];
		return true;
	}	
}