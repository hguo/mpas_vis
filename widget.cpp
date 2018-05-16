#include "widget.h"

CGLWidget::CGLWidget(const QGLFormat& fmt)
  : QGLWidget(fmt)
{
}

CGLWidget::~CGLWidget()
{
}

void CGLWidget::initializeGL()
{

}

void CGLWidget::resizeGL(int w, int h)
{

}

void CGLWidget::paintGL()
{
  glClearColor(1, 1, 1, 0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
}
