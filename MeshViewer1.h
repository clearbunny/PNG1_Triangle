

#ifndef   MESHVIEWER1_H
#define  MESHVIEWER1_H

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

class MeshViewer1 : public QGLViewer
{
	Q_OBJECT

public:
	MeshViewer1(QWidget *parent);
	~MeshViewer1();	
	void BindTexture();

public slots: // slot functions for UI
	void reRender(int textType = 0); // keep current status and re-render, like glflush
	
	void setSaveFile(){	
		QString fileName = QFileDialog::getSaveFileName(this, tr("Save File"),"../models/untitled.off","All Files (*.*)");
		QString null_string="";
		if(fileName != null_string){	
			QByteArray bytearray = fileName.toAscii();	
			const char* filename = bytearray.data();
			if(!meshAssert(NO_MESH_ERROR))
				TriMesh::save(filename, data->cubic_bezier);
		}
	};

	void generate_PNG1()
	{
		data->exact_pn_triangle(2);
		meshInit();
	};

public:
	bool wire_frame;
	database *data;

signals:

protected :
	virtual void draw();
	virtual void drawCubicBezierMesh();
	virtual void drawSphere(TriMesh::Point the);
	virtual void drawCylinder(TriMesh::Point p0,  TriMesh::Point p1);
	virtual void drawFrame(int vertex_idx);
	virtual void keyPressEvent(QKeyEvent *e); // override key event
	virtual void mousePressEvent(QMouseEvent* e); // override mouse press event
	virtual void mouseReleaseEvent(QMouseEvent* e); // override mouse release event
	virtual void mouseMoveEvent(QMouseEvent* e); // override mouse move event
	virtual void mouseDoubleClickEvent(QMouseEvent *e);

private:
	void drawScene(); // method drawing all stuffs
	void drawMyText(int textType); // draw output message on left-bottom screen
	//void drawAnchorAndControl(); // draw anchor points as red spheres, control points as green spheres
	void meshInit(); // initialize the data, call it when reloading
	bool meshAssert(int type)
	{
		switch(type){
		    case NO_MESH_ERROR: if(data->themesh == NULL){reRender(NO_MESH_ERROR); return true;} break;
		    case WRONG_FILE_ERROR: if(data->themesh == NULL){reRender(WRONG_FILE_ERROR); return true;} break;
		}
		return false;
	};

private:
	
};

#endif // _MESHVIEWER1_H___
