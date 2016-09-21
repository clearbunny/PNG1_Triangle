#include "png1_triangle.h"
#include <QtGui/QApplication>

int main(int argc, char *argv[])
{
	QApplication a(argc, argv);
	PNG1_Triangle w;
	w.show();
	return a.exec();
}
