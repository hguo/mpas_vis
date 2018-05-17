#ifndef _WIDGET_H
#define _WIDGET_H

#include "def.h"
#include <QGLWidget>
#include <QMouseEvent>
#include "trackball.h"

class CGLWidget : public QGLWidget {
  Q_OBJECT

public:
  CGLWidget(const QGLFormat& fmt = QGLFormat::defaultFormat());
  ~CGLWidget();

  void loadMeshFromNetCDF(const std::string& filename);

protected:
  void initializeGL();
  void resizeGL(int w, int h);
  void paintGL();
  
  void mousePressEvent(QMouseEvent*); 
  void mouseMoveEvent(QMouseEvent*);
  // void keyPressEvent(QKeyEvent*); 
  void wheelEvent(QWheelEvent*); 

private:
  CGLTrackball trackball;
  QMatrix4x4 projmatrix, mvmatrix; 
  
  const float fovy, znear, zfar; 
  const QVector3D eye, center, up;

protected:
  size_t nCells, nEdges, nVertices;
  double *latVertex, *lonVertex, *xVertex, *yVertex, *zVertex;
  const double radius = 6371220.;
};

#endif
