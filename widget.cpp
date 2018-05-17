#include "widget.h"
#include <OpenGL/glu.h>
#include <cfloat>
#include <netcdf.h>

CGLWidget::CGLWidget(const QGLFormat& fmt) : 
  QGLWidget(fmt),
  fovy(30.f), znear(0.1f), zfar(10.f), 
  eye(0, 0, 2.5), center(0, 0, 0), up(0, 1, 0)
{
}

CGLWidget::~CGLWidget()
{
}

void CGLWidget::loadMeshFromNetCDF(const std::string& filename)
{
  int ncid;
  int dimid_cells, dimid_edges, dimid_vertices;
  int varid_latVertex, varid_lonVertex, varid_xVertex, varid_yVertex, varid_zVertex;

  NC_SAFE_CALL( nc_open(filename.c_str(), NC_NOWRITE, &ncid) );

  NC_SAFE_CALL( nc_inq_dimid(ncid, "nCells", &dimid_cells) );
  NC_SAFE_CALL( nc_inq_dimid(ncid, "nEdges", &dimid_edges) );
  NC_SAFE_CALL( nc_inq_dimid(ncid, "nVertices", &dimid_vertices) );

  NC_SAFE_CALL( nc_inq_dimlen(ncid, dimid_cells, &nCells) );
  NC_SAFE_CALL( nc_inq_dimlen(ncid, dimid_edges, &nEdges) );
  NC_SAFE_CALL( nc_inq_dimlen(ncid, dimid_vertices, &nVertices) );

  NC_SAFE_CALL( nc_inq_varid(ncid, "latVertex", &varid_latVertex) );
  NC_SAFE_CALL( nc_inq_varid(ncid, "lonVertex", &varid_lonVertex) );
  NC_SAFE_CALL( nc_inq_varid(ncid, "xVertex", &varid_xVertex) );
  NC_SAFE_CALL( nc_inq_varid(ncid, "yVertex", &varid_yVertex) );
  NC_SAFE_CALL( nc_inq_varid(ncid, "zVertex", &varid_zVertex) );
 
  latVertex = (double*)malloc(sizeof(double)*nVertices);
  lonVertex = (double*)malloc(sizeof(double)*nVertices);
  xVertex = (double*)malloc(sizeof(double)*nVertices);
  yVertex = (double*)malloc(sizeof(double)*nVertices);
  zVertex = (double*)malloc(sizeof(double)*nVertices);

  const size_t start_vertices[1] = {0}, size_vertices[1] = {nVertices};
  NC_SAFE_CALL( nc_get_vara_double(ncid, varid_latVertex, start_vertices, size_vertices, latVertex) );
  NC_SAFE_CALL( nc_get_vara_double(ncid, varid_lonVertex, start_vertices, size_vertices, lonVertex) );
  NC_SAFE_CALL( nc_get_vara_double(ncid, varid_xVertex, start_vertices, size_vertices, xVertex) );
  NC_SAFE_CALL( nc_get_vara_double(ncid, varid_yVertex, start_vertices, size_vertices, yVertex) );
  NC_SAFE_CALL( nc_get_vara_double(ncid, varid_zVertex, start_vertices, size_vertices, zVertex) );

  NC_SAFE_CALL( nc_close(ncid) );

#if 0
  double xmin = DBL_MAX, xmax = -DBL_MAX, 
         ymin = DBL_MAX, ymax = -DBL_MAX,
         zmin = DBL_MAX, zmax = -DBL_MAX;

  fprintf(stderr, "%zu, %zu, %zu\n", nCells, nEdges, nVertices);
  for (int i=0; i<nVertices; i++) {
    fprintf(stderr, "lon=%f, lat=%f, x=%f, y=%f, z=%f\n", lonVertex[i], latVertex[i], xVertex[i], yVertex[i], zVertex[i]);

    xmin = std::min(xmin, xVertex[i]);
    xmax = std::max(xmax, xVertex[i]);
    ymin = std::min(ymin, yVertex[i]);
    ymax = std::max(ymax, yVertex[i]);
    zmin = std::min(zmin, zVertex[i]);
    zmax = std::max(zmax, zVertex[i]);
  }

  fprintf(stderr, "%f, %f, %f, %f, %f, %f\n", xmin, xmax, ymin, ymax, zmin, zmax);
#endif
}

void CGLWidget::initializeGL()
{
  trackball.init();
}

void CGLWidget::resizeGL(int w, int h)
{
  trackball.reshape(w, h);
  glViewport(0, 0, w, h);
  
  CHECK_GLERROR();
}

void CGLWidget::paintGL()
{
  glClearColor(1, 1, 1, 0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  
  projmatrix.setToIdentity(); 
  projmatrix.perspective(fovy, (float)width()/height(), znear, zfar); 
  mvmatrix.setToIdentity();
  mvmatrix.lookAt(eye, center, up);
  mvmatrix.rotate(trackball.getRotation());
  mvmatrix.scale(trackball.getScale());

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glLoadMatrixd(projmatrix.data()); 
  glMatrixMode(GL_MODELVIEW); 
  glLoadIdentity(); 
  glLoadMatrixd(mvmatrix.data()); 

  glColor3f(0, 0, 0);
  glBegin(GL_POINTS);
  for (int i=0; i<nVertices; i++) {
    glVertex3f(xVertex[i]/radius, yVertex[i]/radius, zVertex[i]/radius);
  }
  glEnd();

  CHECK_GLERROR();
}

void CGLWidget::mousePressEvent(QMouseEvent* e)
{
  trackball.mouse_rotate(e->x(), e->y()); 
}

void CGLWidget::mouseMoveEvent(QMouseEvent* e)
{
  trackball.motion_rotate(e->x(), e->y()); 
  updateGL(); 
}

void CGLWidget::wheelEvent(QWheelEvent* e)
{
  trackball.wheel(e->delta());
  updateGL(); 
}
