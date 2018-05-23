#include <QApplication>
#include "widget.h"

int main(int argc, char **argv)
{
  QApplication app(argc, argv);
  QGLFormat fmt = QGLFormat::defaultFormat();
  fmt.setSampleBuffers(true);
  fmt.setSamples(16); 
  QGLFormat::setDefaultFormat(fmt); 
  
  CGLWidget *widget = new CGLWidget;
  widget->loadMeshFromNetCDF(argv[1]);
  widget->buildKDTree();
  widget->show();
  
  return app.exec();
}
