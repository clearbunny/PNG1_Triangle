#include "MeshViewer.h"
#include <glut.h>

typedef OpenMesh::Vec3d Point;
typedef OpenMesh::Vec3d Vec;
typedef OpenMesh::Vec3d Normal;

MeshViewer::MeshViewer(QWidget *parent) : QGLViewer(parent)
{
	wire_frame=false;
	draw_boundary=false;
	data = database::getInstance();
}

MeshViewer::~MeshViewer()
{
	//data->clear();
	database::Release();
}

void MeshViewer::BindTexture(){
	QImage glImg = QGLWidget::convertToGLFormat(*data->textureimage);  // flipped 32bit RGBA
	// Bind the texture...
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, glImg.width(), glImg.height(), 0,
		GL_RGBA, GL_UNSIGNED_BYTE, glImg.bits());
}

void MeshViewer::meshInit(const char* filename)
{
	delete[] data->themesh;
	data->themesh = TriMesh::read(filename);
	data->themesh->needBoundingBox();
	Point cen = data->themesh->getSceneCenter();
	float rad = (float)data->themesh->getSceneRadius();
	setSceneCenter(qglviewer::Vec((float)cen[0], (float)cen[1], (float)cen[2]));
	setSceneRadius(rad);
	showEntireScene();	
}

inline void MeshViewer::drawScene()
{
	if(data->themesh != NULL)	{
		drawMesh();

		if(draw_boundary)
			drawBoundary();
	}
}

inline void MeshViewer::draw()
{    
	this->setBackgroundColor(QColor(255,255,255));
	glDepthFunc(GL_LESS);
	glEnable(GL_DEPTH_TEST);
	glCullFace(GL_BACK);
	glEnable(GL_CULL_FACE);

	//glLightModeli(GL_LIGHT_MODEL_TWO_SIDE,1.0);

	//qglviewer::Vec cameraPos = camera()->position();
	//const GLfloat pos[4] = {cameraPos[0]*2.0, cameraPos[1]*2.0, cameraPos[2]*2.0, 1.0};
	//glLightfv(GL_LIGHT1, GL_POSITION, pos);

	//// Orientate light along view direction
	//glLightfv(GL_LIGHT1, GL_SPOT_DIRECTION, camera()->viewDirection());
	//glEnable(GL_LIGHT1);

	//// Light default parameters
	//const GLfloat light_ambient[4]  = {1.0, 1.0, 1.0, 1.0};
	//const GLfloat light_specular[4] = {1.0, 1.0, 1.0, 1.0};
	//const GLfloat light_diffuse[4]  = {1.0, 1.0, 1.0, 1.0};

	//glLightf( GL_LIGHT1, GL_SPOT_EXPONENT, 3.0);
	//glLightf( GL_LIGHT1, GL_SPOT_CUTOFF,   10.0);
	//glLightf( GL_LIGHT1, GL_CONSTANT_ATTENUATION,  0.1f);
	//glLightf( GL_LIGHT1, GL_LINEAR_ATTENUATION,    0.3f);
	//glLightf( GL_LIGHT1, GL_QUADRATIC_ATTENUATION, 0.3f);
	//glLightfv(GL_LIGHT1, GL_AMBIENT,  light_ambient);
	//glLightfv(GL_LIGHT1, GL_SPECULAR, light_specular);
	//glLightfv(GL_LIGHT1, GL_DIFFUSE,  light_diffuse);

	//glDisable(GL_DEPTH_TEST);

	drawScene();

	if(data->themesh!=NULL)
		drawMyText(VERTEX_FACE_NUMBER);
}

void MeshViewer::reRender(int textType)
{
	drawScene();
	//drawMyText(textType);
	updateGL();
}

inline void MeshViewer::mousePressEvent(QMouseEvent* e)
{
	bool handled = false;

	//DEFAULT CASE
	if (!handled)
	    QGLViewer::mousePressEvent(e); 
}

inline void MeshViewer::keyPressEvent(QKeyEvent *e)
{
	bool handled = false;
	if(e->key() == Qt::Key_W){
		if(wire_frame==true) 
			wire_frame=false; 
		else 
			wire_frame=true;
		reRender(WIRE_FRAME_SHADING);
		handled=true;
	}if(e->key() == Qt::Key_B){
		if(draw_boundary==true) 
			draw_boundary=false; 
		else 
			draw_boundary=true;
		reRender();
		handled=true;
	}
	

	//DEFAULT CASE
	if (!handled)
		QGLViewer::keyPressEvent(e);
}

inline void MeshViewer::mouseReleaseEvent(QMouseEvent* e)
{
	bool handled = false;

	//DEFAULT CASE
	if (!handled)
	    QGLViewer::mouseReleaseEvent(e);
}

inline void MeshViewer::mouseMoveEvent(QMouseEvent* e)
{
	bool handled = false;

	//DEFAULT CASE
    if (!handled)
	    QGLViewer::mouseMoveEvent(e);
}

inline void MeshViewer::mouseDoubleClickEvent(QMouseEvent *e)
{
	bool handled=false;	
    bool found;

	qglviewer::Vec vec_the = camera()->pointUnderPixel(e->pos(),found);
	Point the;
	the[0]=vec_the[0]; the[1]=vec_the[1]; the[2]=vec_the[2];
	if(!data->themesh->selectInitialPoint( the ) ){
		QMessageBox::information(NULL, "Select Initial Point Failed... ...", " Did Not Find the Center!");
		return;
	}
	
	reRender();
	handled=true;
	
	//DEFAULT CASE
	if(!handled)
		QGLViewer::mouseDoubleClickEvent(e);
}
inline void MeshViewer::drawMyText(int textType)
{
	switch(textType){
		case NO_MESH_ERROR: displayMessage("Error, no mesh loaded"); break;
		case WRONG_FILE_ERROR: displayMessage("Error, cannot load mesh from this file"); break;
		case VERTEX_FACE_NUMBER: displayMessage( QString("Vertex: ").append(QString::number(data->themesh->n_vertices())).append("        Face: ")
				.append(QString::number(data->themesh->n_faces())).append("        AEL: ").append(QString::number(data->themesh->average_edge_length)) ); break;
		}
}

inline void MeshViewer::drawMesh(){
	glEnableClientState(GL_VERTEX_ARRAY);
	glVertexPointer(3, GL_DOUBLE, 0, data->themesh->points());
	glEnableClientState(GL_NORMAL_ARRAY);
	glNormalPointer(GL_DOUBLE, 0, data->themesh->vertex_normals());
	glEnableClientState(GL_COLOR_ARRAY);
	glColorPointer(3, GL_UNSIGNED_BYTE, 0, data->themesh->vertex_colors());
	// draw solid mesh, with polygon offset
	if (data->themesh->vertex_group0.size()>0){
		glPolygonMode(GL_FRONT, GL_FILL);
		glEnable(GL_POLYGON_OFFSET_FILL);
		glPolygonOffset(2.5f,2.5f);
		glDrawElements(GL_TRIANGLES,data->themesh->vertex_group0.size(), GL_UNSIGNED_INT,&(data->themesh->vertex_group0[0]));
		glDisable(GL_POLYGON_OFFSET_FILL);
	}
	glDisableClientState(GL_COLOR_ARRAY);
	if(wire_frame){
		glPolygonMode(GL_FRONT, GL_LINE);
		glColor3d(0.0,0.0,0.0);
		glDrawElements(GL_TRIANGLES,data->themesh->vertex_group0.size(), GL_UNSIGNED_INT,&(data->themesh->vertex_group0[0]));
		glPolygonMode(GL_FRONT, GL_FILL);
	} 
	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_NORMAL_ARRAY);
	
}

inline void MeshViewer::drawBoundary(){
	glColor3d(0.8,0.8,1.0);
	TriMesh::Point the, next;
	for( int i=0; i<data->themesh->boundary.size(); i++){
		the=data->themesh->point(data->themesh->vertex_handle(data->themesh->boundary[i]));
		next=data->themesh->point(data->themesh->vertex_handle(data->themesh->boundary[(i+1)%data->themesh->boundary.size()]));
		drawCylinder(the, next);
	}
}




inline void MeshViewer::drawVertexNormal(){
	TriMesh::Point the,N;
	glColor3d(0.0, 0.0, 1.0);
	qglviewer::Vec from,to;
	
	TriMesh::ConstVertexIter v_it(data->themesh->vertices_begin()),v_end(data->themesh->vertices_end());
	for( ; v_it!=v_end; ++v_it){
		the=data->themesh->point(v_it.handle());
		N=data->themesh->normal(v_it.handle());
		from.x=the[0]; from.y=the[1], from.z=the[2];
		to.x= the[0]+N[0]*0.3;to.y= the[1]+N[1]*0.3;to.z= the[2]+N[2]*0.3;
		drawArrow(from, to, 0.005);
	}
	
}

inline void MeshViewer::drawSphere(TriMesh::Point the, double radius){
	glPushMatrix();
	glTranslated(the[0],the[1],the[2]);
	GLUquadricObj *qobj; 
	qobj = gluNewQuadric(); 
	gluQuadricDrawStyle( qobj, GLU_FILL ); 
	gluQuadricNormals( qobj, GLU_SMOOTH ); 
	gluSphere( qobj, radius, 30, 20 ); 
	glPopMatrix();
}

inline void MeshViewer::drawCylinder(TriMesh::Point p0, TriMesh::Point p1, double radius){
		GLUquadricObj * qobj;
		qobj=gluNewQuadric();
		gluQuadricCallback(qobj,GLU_ERROR,NULL);
		gluQuadricDrawStyle(qobj,GLU_FILL);
		gluQuadricOrientation(qobj,GLU_OUTSIDE);
		gluQuadricNormals(qobj,GLU_SMOOTH);
		glPushMatrix();
		TriMesh::Point z,w,n;
		w=p1-p0;
		GeometryTools::normalize(w);
		z[0]=0.0; 		z[1]=0.0;		z[2]=1.0;
		CROSSPROD3(n,z,w);
		double angle_z_to_w=acos(GeometryTools::dotProduct(z,w))*180.0/PI;
		glTranslated(p0[0],p0[1],p0[2]);
		glRotated(angle_z_to_w,n[0],n[1],n[2]);
		float mat_diffuse[]={1.0, 0.0, 0.0, 1.0f};
		glMaterialfv(GL_FRONT,GL_DIFFUSE,mat_diffuse);
		double height=GeometryTools::distPoint2Point(p0,p1);
		gluCylinder(qobj,radius,radius,height,60.0,80.0);

		glPopMatrix();
}

inline void MeshViewer::drawCone(TriMesh::Point p0, TriMesh::Point p1){
	GLUquadricObj * qobj;
	qobj=gluNewQuadric();
	gluQuadricCallback(qobj,GLU_ERROR,NULL);
	gluQuadricDrawStyle(qobj,GLU_FILL);
	gluQuadricOrientation(qobj,GLU_OUTSIDE);
	gluQuadricNormals(qobj,GLU_SMOOTH);
	glPushMatrix();
	TriMesh::Point z,w,n;
	w=p1-p0;
	GeometryTools::normalize(w);
	z[0]=0.0; 	z[1]=0.0; 	z[2]=1.0;
	n=GeometryTools::crossProduct(z,w);
	double angle_z_to_w=acos(GeometryTools::dotProduct(z,w))*180.0/PI;
	glTranslated(p0[0],p0[1],p0[2]);
	glRotated(angle_z_to_w,n[0],n[1],n[2]);
	float mat_diffuse[]={1.0, 0.0, 0.0, 1.0f};
	glMaterialfv(GL_FRONT,GL_DIFFUSE,mat_diffuse);
	double height=GeometryTools::distPoint2Point(p0,p1);
	gluCylinder(qobj,0.008,0.0,height,15.0,20.0);

	glPopMatrix();
}

inline void MeshViewer::drawFrame(TriMesh::VertexHandle vh, double length, double radius){
	Point the=data->themesh->point(vh);
	Normal n=data->themesh->normal(vh);
	OpenMesh::VPropHandleT<TriMesh::Point> principal_1;
	data->themesh->get_property_handle(principal_1, "principle_1");
	OpenMesh::VPropHandleT<TriMesh::Point> principal_2;
	data->themesh->get_property_handle(principal_2, "principle_2");
	Normal p1=data->themesh->property(principal_1,vh);
	Normal p2=data->themesh->property(principal_2,vh);

	drawFrame(the, n, p1,p2,length,radius);
	//p2=p2-Normal(0.03,0.0,0.0);
	//p2.normalize();
	//qglviewer::Vec from,to;
	//from.x=the[0]; from.y=the[1], from.z=the[2];

	//glColor3d(0.0,0.0,1.0);
	//n*=length*6.2; p1*=length*3.5; p2*=length*9;
	//to.x= the[0]-n[0];to.y= the[1]-n[1];to.z= the[2]-n[2];
	//drawArrow(from,to,radius);

	//to.x= the[0]+p1[0];to.y= the[1]+p1[1];to.z= the[2]+p1[2];
	//drawArrow(from,to,radius);
	//p1=6.2/3.5*Normal(0.0,0.0,-1.0)+Normal(1.0,0.0,0.0);
	//p1.normalize();
	//p1*=length*sqrt(6.2*6.2+3.5*3.5);
	//to.x= the[0]+p1[0];to.y= the[1]+p1[1];to.z= the[2]+p1[2];
	//drawArrow(from,to,radius);

	//to.x= the[0]-p2[0];to.y= the[1]-p2[1];to.z= the[2]-p2[2];
	//drawArrow(from,to,radius);

}

inline void MeshViewer::drawFrame(TriMesh::Point the, TriMesh::Normal N, TriMesh::Normal P1, TriMesh::Normal P2, double length, double radius){
	glColor3d(1.0,1.0,0.0);
	drawSphere(the,radius*3);
	qglviewer::Vec from,to;
	from.x=the[0]; from.y=the[1], from.z=the[2];

	N*=length*1.1; P1*=length; P2*=length;
	to.x= the[0]+N[0];to.y= the[1]+N[1];to.z= the[2]+N[2];
	glColor3d(0.0/255.0, 0.0/255.0,1.0);
	drawArrow(from,to,radius);

	to.x= the[0]+P1[0];to.y= the[1]+P1[1];to.z= the[2]+P1[2];
	glColor3d(1.0,102.0/255.0,102.0/255.0);
	drawArrow(from,to,radius);

	to.x= the[0]+P2[0];to.y= the[1]+P2[1];to.z= the[2]+P2[2];
	glColor3d(204.0/255.0,1.0,204.0/255.0);
	drawArrow(from,to,radius);
}
