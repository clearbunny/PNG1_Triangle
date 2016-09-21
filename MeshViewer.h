/***********************************************************

Filename: MeshViewer.h

Description:

Display everything to screen, using QT4. For detail information, please
also check QT designer, QT's UI file, signal and slots usage. Currently
this class also maintains all data like mesh and region-of-interest. To 
change the user interface, please use QT designer to modify "viewerInterface.Qt4.ui".

Author: Peng Cheng

***********************************************************/

#ifndef   MESHVIEWER_H
#define  MESHVIEWER_H

#include <math.h>
#include "database.h"
#include "GlobalConstant.h"
#include "GeometryTools.h"
#include <qmessagebox.h>
#include <QGLViewer/camera.h>
#include <QtGui>
# include <QMenu>
# include <QKeyEvent>
#include <QFileDialog>

class MeshViewer : public QGLViewer
{
	Q_OBJECT

public:
	MeshViewer(QWidget *parent);
	~MeshViewer();	
	void BindTexture();

public slots: // slot functions for UI
	void reRender(int textType = 0); // keep current status and re-render, like glflush
	void setOpenFile(){
		QString fileName = QFileDialog::getOpenFileName(this, tr("Open File"),".", "Mesh Files (*.obj *.ply *.m *.off);;All files (*.*)");	
		QString null_string="";
		if(fileName!=null_string){
			QByteArray bytearray = fileName.toAscii();
			const char* filename = bytearray.data();
			meshInit(filename);
		}
		else	
			reRender(FILE_NAME_NULL);
	};
	void setSaveFile(){	
		
		data->themesh->save("Saved_Mesh.obj");
	};

public:
	bool wire_frame;
	bool draw_boundary;
	database *data;

signals:

protected :
	virtual void draw();
	virtual void drawMesh();
	virtual void drawBoundary();
	virtual void drawVertexNormal();
	virtual void drawSphere(TriMesh::Point the,double radius=0.01);
	virtual void drawCylinder(TriMesh::Point p0,  TriMesh::Point p1,double radius=0.001);
	virtual void drawCone(TriMesh::Point p0,  TriMesh::Point p1);
	virtual void drawFrame(TriMesh::VertexHandle vh, double length=0.2, double radius=0.004);
	virtual void drawFrame(TriMesh::Point the, TriMesh::Normal N, TriMesh::Normal P1, TriMesh::Normal P2, double length, double radius);
	virtual void keyPressEvent(QKeyEvent *e); // override key event
	virtual void mousePressEvent(QMouseEvent* e); // override mouse press event
	virtual void mouseReleaseEvent(QMouseEvent* e); // override mouse release event
	virtual void mouseMoveEvent(QMouseEvent* e); // override mouse move event
	virtual void mouseDoubleClickEvent(QMouseEvent *e);

private:
	void drawScene(); // method drawing all stuffs
	void drawMyText(int textType); // draw output message on left-bottom screen
	//void drawAnchorAndControl(); // draw anchor points as red spheres, control points as green spheres
	void meshInit(const char* filename); // initialize the data, call it when reloading
	bool meshAssert(int type)	{
		switch(type){
		case NO_MESH_ERROR: if(data->themesh == NULL){reRender(NO_MESH_ERROR); return true;} break;
		case WRONG_FILE_ERROR: if(data->themesh == NULL){reRender(WRONG_FILE_ERROR); return true;} break;
		}
		return false;
	};

private:
	
};

#endif // _MESHVIEWER_H___
