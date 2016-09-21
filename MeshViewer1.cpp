#include "MeshViewer1.h"
#include <glut.h>

typedef OpenMesh::Vec3d Point;
typedef OpenMesh::Vec3d Vec;
typedef OpenMesh::Vec3d Normal;

MeshViewer1::MeshViewer1(QWidget *parent) : QGLViewer(parent)
{
	wire_frame=false;
	data = database::getInstance();
	this->setBackgroundColor(QColor(255,255,255));
}

MeshViewer1::~MeshViewer1()
{
	//data->clear();
	database::Release();
}

void MeshViewer1::BindTexture(){
	QImage glImg = QGLWidget::convertToGLFormat(*data->textureimage);  // flipped 32bit RGBA
	// Bind the texture...
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, glImg.width(), glImg.height(), 0,
		GL_RGBA, GL_UNSIGNED_BYTE, glImg.bits());
}

inline void MeshViewer1::drawScene()
{
	if(data->cubic_bezier != NULL)
	{
		drawCubicBezierMesh();

	}
}

inline void MeshViewer1::draw()
{    
	this->setBackgroundColor(QColor(255,255,255));
	glDepthFunc(GL_LESS);
	glEnable(GL_DEPTH_TEST);
	glCullFace(GL_BACK);
	glEnable(GL_CULL_FACE);
	drawScene();
	if(data->cubic_bezier!=NULL)
		drawMyText(VERTEX_FACE_NUMBER);
}

void MeshViewer1::reRender(int textType)
{
	drawScene();
	//drawMyText(textType);
	updateGL();
}

inline void MeshViewer1::mousePressEvent(QMouseEvent* e)
{
	bool handled = false;

	//DEFAULT CASE
	if (!handled)
	    QGLViewer::mousePressEvent(e); 
}

inline void MeshViewer1::keyPressEvent(QKeyEvent *e)
{
	bool handled = false;
	if(e->key() == Qt::Key_W){
		if(wire_frame==true) 
			wire_frame=false; 
		else 
			wire_frame=true;
		reRender(WIRE_FRAME_SHADING);
		handled=true;
	}else if(e->key() == Qt::Key_1){
		//data->PN_TRIANGLE(1);
		data->exact_pn_triangle(1);
		meshInit();
		handled=true;
	}else if(e->key() == Qt::Key_2){
		//data->PN_TRIANGLE(2);
		data->exact_pn_triangle(2);
		meshInit();
		handled=true;
	}else if(e->key() == Qt::Key_3){
		//data->PN_TRIANGLE(3);
		data->exact_pn_triangle(3);
		meshInit();
		handled=true;
	}else if(e->key() == Qt::Key_4){
		//data->PN_TRIANGLE(4);
		data->exact_pn_triangle(4);
		meshInit();
		handled=true;
	}else if(e->key() == Qt::Key_5){
		//data->PN_TRIANGLE(5);
		data->exact_pn_triangle(5);
		meshInit();
		handled=true;
	}else if(e->key() == Qt::Key_6){
		//data->PN_TRIANGLE(6);
		data->exact_pn_triangle(6);
		meshInit();
	}else if(e->key() == Qt::Key_7){
		//data->PN_TRIANGLE(7);
		data->exact_pn_triangle(7);
		meshInit();
		handled=true;
	}else if(e->key() == Qt::Key_8){
		//data->PN_TRIANGLE(8);
		//data->exact_pn_triangle(50);
		data->exact_torus(50);
		meshInit();
		handled=true;
	}else if(e->key() == Qt::Key_S){
		OpenMesh::IO::Options wopt;
		wopt.check(OpenMesh::IO::Options::VertexNormal);
		QString fileName = QFileDialog::getSaveFileName(this, tr("Save File"),"../models/untitled.off","All Files (*.*)");
		QString null_string="";
		if(fileName != null_string){	
			QByteArray bytearray = fileName.toAscii();	
			const char* filename = bytearray.data();
			if(!meshAssert(NO_MESH_ERROR))
				TriMesh::save(filename, data->cubic_bezier,&wopt);
		}
		handled=true;
	}

	//DEFAULT CASE
	if (!handled)
		QGLViewer::keyPressEvent(e);
}

inline void MeshViewer1::mouseReleaseEvent(QMouseEvent* e)
{
	bool handled = false;

	//DEFAULT CASE
	if (!handled)
	    QGLViewer::mouseReleaseEvent(e);
}

inline void MeshViewer1::mouseMoveEvent(QMouseEvent* e)
{
	bool handled = false;

	//DEFAULT CASE
    if (!handled)
	    QGLViewer::mouseMoveEvent(e);
}

inline void MeshViewer1::mouseDoubleClickEvent(QMouseEvent *e)
{
	bool handled=false;	
    bool found;

	qglviewer::Vec vec_the = camera()->pointUnderPixel(e->pos(),found);
	Point the;
	the[0]=vec_the[0]; the[1]=vec_the[1]; the[2]=vec_the[2];
		
	reRender();
	handled=true;
	
	//DEFAULT CASE
	if(!handled)
		QGLViewer::mouseDoubleClickEvent(e);
}
inline void MeshViewer1::drawMyText(int textType)
{
	switch(textType){
		case NO_MESH_ERROR: displayMessage("Error, no mesh loaded"); break;
		case WRONG_FILE_ERROR: displayMessage("Error, cannot load mesh from this file"); break;
	    case VERTEX_FACE_NUMBER: displayMessage( QString("Vertex: ").append(QString::number(data->cubic_bezier->n_vertices())).append("        Face: ")
									 .append(QString::number(data->cubic_bezier->n_faces())).append("        AEL: ").append(QString::number(data->cubic_bezier->average_edge_length)) ); break;
		}
}

void MeshViewer1::meshInit()
{
	Point cen = data->cubic_bezier->getSceneCenter();
	float rad = (float)data->cubic_bezier->getSceneRadius();
	setSceneCenter(qglviewer::Vec((float)cen[0], (float)cen[1], (float)cen[2]));
	setSceneRadius(rad);
	showEntireScene();

	if (!data->cubic_bezier->has_vertex_colors())	{
		data->cubic_bezier->request_vertex_colors();
		data->cubic_bezier->set_vertex_color(255, 255, 153);
	}

	data->cubic_bezier->average_edge_length=data->cubic_bezier->getAverageEdgeLength();

	TriMesh::ConstFaceIter fit(data->cubic_bezier->faces_begin()), fEnd(data->cubic_bezier->faces_end());
	data->cubic_bezier->vertex_group0.clear();
	data->cubic_bezier->vertex_group2.clear();
	for(; fit!=fEnd; ++fit)
	{
		TriMesh::ConstFaceVertexIter fvit = data->cubic_bezier->cfv_iter(fit.handle());
		data->cubic_bezier->vertex_group0.push_back(fvit.handle().idx());	
		data->cubic_bezier->vertex_group2.push_back(fvit.handle().idx());    ++fvit;
		data->cubic_bezier->vertex_group0.push_back(fvit.handle().idx());
		data->cubic_bezier->vertex_group2.push_back(fvit.handle().idx());     ++fvit;
		data->cubic_bezier->vertex_group0.push_back(fvit.handle().idx());
		data->cubic_bezier->vertex_group2.push_back(fvit.handle().idx());
	}

}


void MeshViewer1::drawCubicBezierMesh(){
	glEnableClientState(GL_VERTEX_ARRAY);
	glVertexPointer(3, GL_DOUBLE, 0, data->cubic_bezier->points());
	glEnableClientState(GL_NORMAL_ARRAY);
	glNormalPointer(GL_DOUBLE, 0, data->cubic_bezier->vertex_normals());
	glEnableClientState(GL_COLOR_ARRAY);
	glColorPointer(3, GL_UNSIGNED_BYTE, 0, data->cubic_bezier->vertex_colors());
	// draw solid mesh, with polygon offset
	if (data->cubic_bezier->vertex_group0.size()>0){
		glPolygonMode(GL_FRONT, GL_FILL);
		glEnable(GL_POLYGON_OFFSET_FILL);
		glPolygonOffset(2.5f,2.5f);
		glDrawElements(GL_TRIANGLES,data->cubic_bezier->vertex_group0.size(), GL_UNSIGNED_INT,&(data->cubic_bezier->vertex_group0[0]));
		glDisable(GL_POLYGON_OFFSET_FILL);
	}
	glDisableClientState(GL_COLOR_ARRAY);
	if(wire_frame){
		glPolygonMode(GL_FRONT, GL_LINE);
		glColor3d(0.0,0.0,0.0);
		glDrawElements(GL_TRIANGLES,data->cubic_bezier->vertex_group0.size(), GL_UNSIGNED_INT,&(data->cubic_bezier->vertex_group0[0]));
	} 
	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_NORMAL_ARRAY);
}

void MeshViewer1::drawSphere(TriMesh::Point the){
	glPushMatrix();
	glTranslated(the[0],the[1],the[2]);
	GLUquadricObj *qobj; 
	qobj = gluNewQuadric(); 
	gluQuadricDrawStyle( qobj, GLU_FILL ); 
	gluQuadricNormals( qobj, GLU_SMOOTH ); 
	gluSphere( qobj, 0.01, 15, 10 ); 
	glPopMatrix();
}

void MeshViewer1::drawCylinder(TriMesh::Point p0, TriMesh::Point p1){
		GLUquadricObj * qobj;
		qobj=gluNewQuadric();
		gluQuadricCallback(qobj,GLU_ERROR,NULL);
		gluQuadricDrawStyle(qobj,GLU_FILL);
		gluQuadricOrientation(qobj,GLU_OUTSIDE);
		gluQuadricNormals(qobj,GLU_SMOOTH);
		glPushMatrix();
		TriMesh::Point z,w,n;
		w=p1-p0;
		w.normalized();
		z[0]=0.0;
		z[1]=0.0;
		z[2]=1.0;
		n=GeometryTools::crossProduct(z,w);
		double angle_z_to_w=acos(GeometryTools::dotProduct(z,w))*180.0/PI;
		glTranslated(p0[0],p0[1],p0[2]);
		glRotated(angle_z_to_w,n[0],n[1],n[2]);
		float mat_diffuse[]={1.0, 0.0, 0.0, 1.0f};
		glMaterialfv(GL_FRONT,GL_DIFFUSE,mat_diffuse);
		double height=GeometryTools::distPoint2Point(p0,p1);
		gluCylinder(qobj,0.003,0.003,height,15.0,20.0);

		glPopMatrix();
}

void MeshViewer1::drawFrame(int vertex_idx){
	glPushMatrix();
	Point the=data->cubic_bezier->point(data->cubic_bezier->vertex_handle(vertex_idx));
	glTranslated(the[0],the[1],the[2]);
	Point z=data->cubic_bezier->normal(data->cubic_bezier->vertex_handle(vertex_idx));
	Point zz(0.0,0.0,1.0);
	Point axis=GeometryTools::crossProduct(zz,z);
	glRotated(acos(GeometryTools::dotProduct(zz,z))*180.0/PI, axis[0], axis[1], axis[2]);
	drawAxis(0.3);
	glPopMatrix();
}
