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
  const double radius = 6371220.;
  
  size_t nCells, nEdges, nVertices, nVertLevels;
  std::vector<double> latVertex, lonVertex, xVertex, yVertex, zVertex;
  std::vector<double> xCell, yCell, zCell;
  std::vector<int> indexToVertexID, verticesOnEdge;
  std::vector<double> velocityX, velocityY, velocityZ;
};

#endif
