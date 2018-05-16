#ifndef _WIDGET_H
#define _WIDGET_H

#include <QGLWidget>

class CGLWidget : public QGLWidget {
  Q_OBJECT

public:
  CGLWidget(const QGLFormat& fmt = QGLFormat::defaultFormat());
  ~CGLWidget();

protected:
  void initializeGL();
  void resizeGL(int w, int h);
  void paintGL();
};

#endif
