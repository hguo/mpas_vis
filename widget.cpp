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
  int dimid_cells, dimid_edges, dimid_vertices, dimid_vertLevels;
  int varid_latVertex, varid_lonVertex, varid_xVertex, varid_yVertex, varid_zVertex,
      varid_latCell, varid_lonCell, varid_xCell, varid_yCell, varid_zCell, 
      varid_verticesOnEdge, varid_indexToVertexID,
      varid_velocityX, varid_velocityY, varid_velocityZ;

  NC_SAFE_CALL( nc_open(filename.c_str(), NC_NOWRITE, &ncid) );

  NC_SAFE_CALL( nc_inq_dimid(ncid, "nCells", &dimid_cells) );
  NC_SAFE_CALL( nc_inq_dimid(ncid, "nEdges", &dimid_edges) );
  NC_SAFE_CALL( nc_inq_dimid(ncid, "nVertices", &dimid_vertices) );
  NC_SAFE_CALL( nc_inq_dimid(ncid, "nVertLevels", &dimid_vertLevels) );

  NC_SAFE_CALL( nc_inq_dimlen(ncid, dimid_cells, &nCells) );
  NC_SAFE_CALL( nc_inq_dimlen(ncid, dimid_edges, &nEdges) );
  NC_SAFE_CALL( nc_inq_dimlen(ncid, dimid_vertices, &nVertices) );
  NC_SAFE_CALL( nc_inq_dimlen(ncid, dimid_vertLevels, &nVertLevels) );

  NC_SAFE_CALL( nc_inq_varid(ncid, "indexToVertexID", &varid_indexToVertexID) );
  NC_SAFE_CALL( nc_inq_varid(ncid, "xCell", &varid_xCell) );
  NC_SAFE_CALL( nc_inq_varid(ncid, "yCell", &varid_yCell) );
  NC_SAFE_CALL( nc_inq_varid(ncid, "zCell", &varid_zCell) );
  NC_SAFE_CALL( nc_inq_varid(ncid, "latVertex", &varid_latVertex) );
  NC_SAFE_CALL( nc_inq_varid(ncid, "lonVertex", &varid_lonVertex) );
  NC_SAFE_CALL( nc_inq_varid(ncid, "xVertex", &varid_xVertex) );
  NC_SAFE_CALL( nc_inq_varid(ncid, "yVertex", &varid_yVertex) );
  NC_SAFE_CALL( nc_inq_varid(ncid, "zVertex", &varid_zVertex) );
  NC_SAFE_CALL( nc_inq_varid(ncid, "verticesOnEdge", &varid_verticesOnEdge) );
  NC_SAFE_CALL( nc_inq_varid(ncid, "velocityX", &varid_velocityX) );
  NC_SAFE_CALL( nc_inq_varid(ncid, "velocityY", &varid_velocityY) );
  NC_SAFE_CALL( nc_inq_varid(ncid, "velocityZ", &varid_velocityZ) );

  const size_t start_cells[1] = {0}, size_cells[1] = {nCells};
  xCell.resize(nCells);
  yCell.resize(nCells);
  zCell.resize(nCells);

  NC_SAFE_CALL( nc_get_vara_double(ncid, varid_xCell, start_cells, size_cells, &xCell[0]) );
  NC_SAFE_CALL( nc_get_vara_double(ncid, varid_yCell, start_cells, size_cells, &yCell[0]) );
  NC_SAFE_CALL( nc_get_vara_double(ncid, varid_zCell, start_cells, size_cells, &zCell[0]) );

  const size_t start_vertices[1] = {0}, size_vertices[1] = {nVertices};
  latVertex.resize(nVertices);
  lonVertex.resize(nVertices);
  xVertex.resize(nVertices);
  yVertex.resize(nVertices);
  zVertex.resize(nVertices);
  indexToVertexID.resize(nVertices);

  NC_SAFE_CALL( nc_get_vara_double(ncid, varid_latVertex, start_vertices, size_vertices, &latVertex[0]) );
  NC_SAFE_CALL( nc_get_vara_double(ncid, varid_lonVertex, start_vertices, size_vertices, &lonVertex[0]) );
  NC_SAFE_CALL( nc_get_vara_double(ncid, varid_xVertex, start_vertices, size_vertices, &xVertex[0]) );
  NC_SAFE_CALL( nc_get_vara_double(ncid, varid_yVertex, start_vertices, size_vertices, &yVertex[0]) );
  NC_SAFE_CALL( nc_get_vara_double(ncid, varid_zVertex, start_vertices, size_vertices, &zVertex[0]) );
  NC_SAFE_CALL( nc_get_vara_int(ncid, varid_indexToVertexID, start_vertices, size_vertices, &indexToVertexID[0]) );

  // for (int i=0; i<nVertices; i++) 
  //   fprintf(stderr, "%d, %d\n", i, indexToVertexID[i]);

  const size_t start_edges2[2] = {0, 0}, size_edges2[2] = {nEdges, 2};
  verticesOnEdge.resize(nEdges*2);

  NC_SAFE_CALL( nc_get_vara_int(ncid, varid_verticesOnEdge, start_edges2, size_edges2, &verticesOnEdge[0]) );

  // for (int i=0; i<nEdges; i++) 
  //   fprintf(stderr, "%d, %d\n", verticesOnEdge[i*2], verticesOnEdge[i*2+1]);


  const size_t start_time_cell_level[3] = {0, 0, 0}, size_time_cell_level[3] = {1, nCells, nVertLevels};
  velocityX.resize(nCells*nVertLevels);
  velocityY.resize(nCells*nVertLevels);
  velocityZ.resize(nCells*nVertLevels);

  NC_SAFE_CALL( nc_get_vara_double(ncid, varid_velocityX, start_time_cell_level, size_time_cell_level, &velocityX[0]) );
  NC_SAFE_CALL( nc_get_vara_double(ncid, varid_velocityY, start_time_cell_level, size_time_cell_level, &velocityY[0]) );
  NC_SAFE_CALL( nc_get_vara_double(ncid, varid_velocityZ, start_time_cell_level, size_time_cell_level, &velocityZ[0]) );

  NC_SAFE_CALL( nc_close(ncid) );
  
  fprintf(stderr, "%zu, %zu, %zu, %zu\n", nCells, nEdges, nVertices, nVertLevels);

#if 0
  double xmin = DBL_MAX, xmax = -DBL_MAX, 
         ymin = DBL_MAX, ymax = -DBL_MAX,
         zmin = DBL_MAX, zmax = -DBL_MAX;

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
  glPointSize(2.0);
  
#if 1 // vertices
  glBegin(GL_POINTS);
  for (int i=0; i<nVertices; i++) 
    glVertex3f(xVertex[i]/radius, yVertex[i]/radius, zVertex[i]/radius);
  glEnd();
#endif

#if 1 // cell centers
  for (int i=0; i<nCells; i++) {
    const double x = xCell[i] / radius, y = yCell[i] / radius, z = zCell[i] / radius;
    const double vx = velocityX[i], vy = velocityY[i], vz = velocityZ[i];
    const double alpha = 0.01;

    glBegin(GL_LINES);
    glVertex3f(x, y, z);
    glVertex3f(x + alpha*vx, y + alpha*vy, z + alpha*vz);
    glEnd();

    glBegin(GL_POINTS);
    glVertex3f(x, y, z); 
    glEnd();
  }
#endif

#if 1 // edges
  glBegin(GL_LINES);
  for (int i=0; i<nEdges; i++) {
    const int i0 = verticesOnEdge[i*2]-1, i1 = verticesOnEdge[i*2+1]-1;
    glVertex3f(xVertex[i0]/radius, yVertex[i0]/radius, zVertex[i0]/radius);
    glVertex3f(xVertex[i1]/radius, yVertex[i1]/radius, zVertex[i1]/radius);
  }
  glEnd();
#endif

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
