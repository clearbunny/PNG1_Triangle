#ifndef PNG1_TRIANGLE_H
#define PNG1_TRIANGLE_H

#include <QtGui/QMainWindow>
#include "ui_png1_triangle.h"

class PNG1_Triangle : public QMainWindow
{
	Q_OBJECT

public:
	PNG1_Triangle(QWidget *parent = 0, Qt::WFlags flags = 0);
	~PNG1_Triangle();

private:
	Ui::PNG1_TriangleClass ui;
};

#endif // PNG1_TRIANGLE_H
