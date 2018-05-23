#include "widget.h"
#include <OpenGL/glu.h>
#include <cfloat>
#include <netcdf.h>

template <typename T>
static inline void cross_product(const T A[3], const T B[3], T C[3])
{
  C[0] = A[1]*B[2] - A[2]*B[1]; 
  C[1] = A[2]*B[0] - A[0]*B[2]; 
  C[2] = A[0]*B[1] - A[1]*B[0];
}

template <typename T>
static inline T dot_product(const T A[3], const T B[3])
{
  return A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
}

template <typename T>
static inline T barycentric_point2triangle(const T P0[3], const T P1[3], const T P2[3], const T P[3], T lambda[3])
{
  T u[3] = {P1[0] - P0[0], P1[1] - P0[1], P1[2] - P0[2]}, 
    v[3] = {P2[0] - P0[0], P2[1] - P0[1], P2[2] - P0[2]}, 
    w[3] = {P[0] - P0[0], P[1] - P0[1], P[2] - P0[2]};
  T n[3], uw[3], wv[3];

  cross_product(u, v, n); // n = u x v
  cross_product(u, w, uw); // uw = u x w
  cross_product(w, v, wv); // wv = w x v

  T n2 = dot_product(n, n); // n^2

  lambda[2] = dot_product(uw, n) / n2;
  lambda[1] = dot_product(wv, n) / n2;
  lambda[0] = 1 - lambda[2] - lambda[1];
}


CGLWidget::CGLWidget(const QGLFormat& fmt) : 
  QGLWidget(fmt),
  fovy(30.f), znear(0.1f), zfar(10.f), 
  eye(0, 0, 2.5), center(0, 0, 0), up(0, 1, 0)
{
}

CGLWidget::~CGLWidget()
{
}

void CGLWidget::buildKDTree()
{
  std::vector<int> cells(nCells);
  for (int i=0; i<nCells; i++) cells[i] = i;

  kdnode *root = new kdnode; 

  buildKDTreeRecursively(root, cells, 0);
}

void CGLWidget::buildKDTreeRecursively(kdnode *node, std::vector<int>& cells, int dir)
{
  const size_t size = cells.size();
  const size_t half_size = size / 2;

  if (size > 1) {
    std::stable_sort(cells.begin(), cells.end(), [this, dir](int c0, int c1) {return xyzCell[c0*3+dir] < xyzCell[c1*3+dir];});

    node->lo_val = xyzCell[cells.front()*3+dir];
    node->hi_val = xyzCell[cells.back()*3+dir];

    node->l = new kdnode;
    node->r = new kdnode;

    std::vector<int> lo(cells.begin(), cells.begin() + half_size);
    std::vector<int> hi(cells.begin() + half_size, cells.end());

    buildKDTreeRecursively(node->l, lo, (dir+1)%3);
    buildKDTreeRecursively(node->r, hi, (dir+1)%3);
  } else {
    node->cellId = cells[0];
    // TODO
  }
}

void CGLWidget::loadMeshFromNetCDF(const std::string& filename)
{
  int ncid;
  int dimid_cells, dimid_edges, dimid_vertices, dimid_vertLevels;
  int varid_latVertex, varid_lonVertex, varid_xVertex, varid_yVertex, varid_zVertex,
      varid_latCell, varid_lonCell, varid_xCell, varid_yCell, varid_zCell, 
      varid_verticesOnEdge, varid_cellsOnVertex, 
      varid_indexToVertexID, varid_indexToCellID,
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
  NC_SAFE_CALL( nc_inq_varid(ncid, "indexToCellID", &varid_indexToCellID) );
  NC_SAFE_CALL( nc_inq_varid(ncid, "latCell", &varid_latCell) );
  NC_SAFE_CALL( nc_inq_varid(ncid, "lonCell", &varid_lonCell) );
  NC_SAFE_CALL( nc_inq_varid(ncid, "xCell", &varid_xCell) );
  NC_SAFE_CALL( nc_inq_varid(ncid, "yCell", &varid_yCell) );
  NC_SAFE_CALL( nc_inq_varid(ncid, "zCell", &varid_zCell) );
  NC_SAFE_CALL( nc_inq_varid(ncid, "latVertex", &varid_latVertex) );
  NC_SAFE_CALL( nc_inq_varid(ncid, "lonVertex", &varid_lonVertex) );
  NC_SAFE_CALL( nc_inq_varid(ncid, "xVertex", &varid_xVertex) );
  NC_SAFE_CALL( nc_inq_varid(ncid, "yVertex", &varid_yVertex) );
  NC_SAFE_CALL( nc_inq_varid(ncid, "zVertex", &varid_zVertex) );
  NC_SAFE_CALL( nc_inq_varid(ncid, "verticesOnEdge", &varid_verticesOnEdge) );
  NC_SAFE_CALL( nc_inq_varid(ncid, "cellsOnVertex", &varid_cellsOnVertex) );
  NC_SAFE_CALL( nc_inq_varid(ncid, "velocityX", &varid_velocityX) );
  NC_SAFE_CALL( nc_inq_varid(ncid, "velocityY", &varid_velocityY) );
  NC_SAFE_CALL( nc_inq_varid(ncid, "velocityZ", &varid_velocityZ) );

  const size_t start_cells[1] = {0}, size_cells[1] = {nCells};

  indexToCellID.resize(nCells);
  NC_SAFE_CALL( nc_get_vara_int(ncid, varid_indexToCellID, start_cells, size_cells, &indexToCellID[0]) );
  for (int i=0; i<nCells; i++) {
    cellIndex[indexToCellID[i]] = i;
    // fprintf(stderr, "%d, %d\n", i, indexToCellID[i]);
  }

  std::vector<double> coord_cells;
  coord_cells.resize(nCells);
  xyzCell.resize(nCells*3);
  NC_SAFE_CALL( nc_get_vara_double(ncid, varid_xCell, start_cells, size_cells, &coord_cells[0]) );
  for (int i=0; i<nCells; i++) xyzCell[i*3] = coord_cells[i];
  NC_SAFE_CALL( nc_get_vara_double(ncid, varid_yCell, start_cells, size_cells, &coord_cells[0]) );
  for (int i=0; i<nCells; i++) xyzCell[i*3+1] = coord_cells[i];
  NC_SAFE_CALL( nc_get_vara_double(ncid, varid_zCell, start_cells, size_cells, &coord_cells[0]) );
  for (int i=0; i<nCells; i++) xyzCell[i*3+2] = coord_cells[i];

  const size_t start_vertices[1] = {0}, size_vertices[1] = {nVertices};
  latVertex.resize(nVertices);
  lonVertex.resize(nVertices);
  xVertex.resize(nVertices);
  yVertex.resize(nVertices);
  zVertex.resize(nVertices);
  indexToVertexID.resize(nVertices);

  NC_SAFE_CALL( nc_get_vara_int(ncid, varid_indexToVertexID, start_vertices, size_vertices, &indexToVertexID[0]) );
  NC_SAFE_CALL( nc_get_vara_double(ncid, varid_latVertex, start_vertices, size_vertices, &latVertex[0]) );
  NC_SAFE_CALL( nc_get_vara_double(ncid, varid_lonVertex, start_vertices, size_vertices, &lonVertex[0]) );
  NC_SAFE_CALL( nc_get_vara_double(ncid, varid_xVertex, start_vertices, size_vertices, &xVertex[0]) );
  NC_SAFE_CALL( nc_get_vara_double(ncid, varid_yVertex, start_vertices, size_vertices, &yVertex[0]) );
  NC_SAFE_CALL( nc_get_vara_double(ncid, varid_zVertex, start_vertices, size_vertices, &zVertex[0]) );

  for (int i=0; i<nVertices; i++) {
    vertexIndex[indexToVertexID[i]] = i;
    // fprintf(stderr, "%d, %d\n", i, indexToVertexID[i]);
  }

  const size_t start_edges2[2] = {0, 0}, size_edges2[2] = {nEdges, 2};
  verticesOnEdge.resize(nEdges*2);

  NC_SAFE_CALL( nc_get_vara_int(ncid, varid_verticesOnEdge, start_edges2, size_edges2, &verticesOnEdge[0]) );

  // for (int i=0; i<nEdges; i++) 
  //   fprintf(stderr, "%d, %d\n", verticesOnEdge[i*2], verticesOnEdge[i*2+1]);

  const size_t start_vertex_cell[2] = {0, 0}, size_vertex_cell[2] = {nVertices, 3};
  cellsOnVertex.resize(nVertices*3);

  NC_SAFE_CALL( nc_get_vara_int(ncid, varid_cellsOnVertex, start_vertex_cell, size_vertex_cell, &cellsOnVertex[0]) );
 

  const size_t start_time_cell_level[3] = {0, 0, 0}, size_time_cell_level[3] = {1, nCells, nVertLevels};
  velocityX.resize(nCells*nVertLevels);
  velocityY.resize(nCells*nVertLevels);
  velocityZ.resize(nCells*nVertLevels);

  NC_SAFE_CALL( nc_get_vara_double(ncid, varid_velocityX, start_time_cell_level, size_time_cell_level, &velocityX[0]) );
  NC_SAFE_CALL( nc_get_vara_double(ncid, varid_velocityY, start_time_cell_level, size_time_cell_level, &velocityY[0]) );
  NC_SAFE_CALL( nc_get_vara_double(ncid, varid_velocityZ, start_time_cell_level, size_time_cell_level, &velocityZ[0]) );

  NC_SAFE_CALL( nc_close(ncid) );
  
  fprintf(stderr, "%zu, %zu, %zu, %zu\n", nCells, nEdges, nVertices, nVertLevels);


  // derive velocity on verticies
  velocityXv.resize(nVertices*nVertLevels);
  velocityYv.resize(nVertices*nVertLevels);
  velocityZv.resize(nVertices*nVertLevels);
  
  for (int i=0; i<nVertices; i++) {
    if (cellsOnVertex[i*3] == 0 || cellsOnVertex[i*3+1] == 0 || cellsOnVertex[i*3+2] == 0) continue; // on boundary
    int c0 = cellIndex[cellsOnVertex[i*3]], c1 = cellIndex[cellsOnVertex[i*3+1]], c2 = cellIndex[cellsOnVertex[i*3+2]];

    double X[3][3] = {
      {xyzCell[c0*3], xyzCell[c0*3+1], xyzCell[c0*3+2]}, 
      {xyzCell[c1*3], xyzCell[c1*3+1], xyzCell[c1*3+2]}, 
      {xyzCell[c2*3], xyzCell[c2*3+1], xyzCell[c2*3+2]}, 
    };

    double P[3] = {xVertex[i], yVertex[i], zVertex[i]};

    double lambda[3];
    barycentric_point2triangle(X[0], X[1], X[2], P, lambda);

    velocityXv[i] = lambda[0] * velocityX[c0] + lambda[1] * velocityX[c1] + lambda[2] * velocityX[c2];
    velocityYv[i] = lambda[0] * velocityY[c0] + lambda[1] * velocityY[c1] + lambda[2] * velocityY[c2];
    velocityZv[i] = lambda[0] * velocityZ[c0] + lambda[1] * velocityZ[c1] + lambda[2] * velocityZ[c2];

    // fprintf(stderr, "vertex %d: %f, %f, %f\n", i, lambda[0], lambda[1], lambda[2]);
  }
  
  // for (int i=0; i<nVertices; i++) {
  //   fprintf(stderr, "%d, %d, %d\n", cellsOnVertex[i*3], cellsOnVertex[i*3+1], cellsOnVertex[i*3+2]);
  // }


#if 1
  double xmin = DBL_MAX, xmax = -DBL_MAX, 
         ymin = DBL_MAX, ymax = -DBL_MAX,
         zmin = DBL_MAX, zmax = -DBL_MAX, 
         latmin = DBL_MAX, latmax = -DBL_MAX, 
         lonmin = DBL_MAX, lonmax = -DBL_MAX;

  for (int i=0; i<nVertices; i++) {
    // fprintf(stderr, "lon=%f, lat=%f, x=%f, y=%f, z=%f\n", lonVertex[i], latVertex[i], xVertex[i], yVertex[i], zVertex[i]);

    xmin = std::min(xmin, xVertex[i]);
    xmax = std::max(xmax, xVertex[i]);
    ymin = std::min(ymin, yVertex[i]);
    ymax = std::max(ymax, yVertex[i]);
    zmin = std::min(zmin, zVertex[i]);
    zmax = std::max(zmax, zVertex[i]);
    latmin = std::min(latmin, latVertex[i]); 
    latmax = std::max(latmax, latVertex[i]);
    lonmin = std::min(lonmin, lonVertex[i]); 
    lonmax = std::max(lonmax, lonVertex[i]);
  }

  fprintf(stderr, "%f, %f, %f, %f, %f, %f, latmin=%f, latmax=%f, lonmin=%f, lonmax=%f\n", 
      xmin, xmax, ymin, ymax, zmin, zmax, latmin, latmax, lonmin, lonmax);
#endif
}

void CGLWidget::initializeGL()
{
  glewInit();

  trackball.init();

  sphere = gluNewQuadric();
  gluQuadricDrawStyle(sphere, GLU_FILL);
  gluQuadricTexture(sphere, TRUE);
  gluQuadricNormals(sphere, GLU_SMOOTH);
    
  glEnable(GL_POLYGON_OFFSET_FILL); 
  glPolygonOffset(1, 1);
    
  glEnable(GL_DEPTH_TEST);

  CHECK_GLERROR();
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

  if (sphere != NULL) {
    if (texEarth == 0) { 
      // QImage earthImage("./earth4096.jpg");
      // QImage earthImage("./water_16k.png");
      QImage earthImage("./water_8k.png");
      texEarth = bindTexture(earthImage); 
    }
    // fprintf(stderr, "%d\n", texEarth);

    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, texEarth);
 
    glCullFace(GL_BACK);
    glEnable(GL_CULL_FACE);

    glColor4f(1, 1, 1, 1);
    glPushMatrix();
    glRotatef(-90, 0, 0, 1);
    gluSphere(sphere, 1.0, 72, 60);
    glPopMatrix();

    glDisable(GL_CULL_FACE);
    glDisable(GL_TEXTURE_2D);
  }
  
  glColor3f(0, 0, 0);
  glPointSize(2.0);
  
#if 1 // vertex velocities
  for (int i=0; i<nVertices; i++) {
    const double x = xVertex[i] / radius, y = yVertex[i] / radius, z = zVertex[i] / radius;
    const double vx = velocityXv[i], vy = velocityYv[i], vz = velocityZv[i];
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

#if 1 // cell center velocities
  for (int i=0; i<nCells; i++) {
    const double x = xyzCell[i*3] / radius, y = xyzCell[i*3+1] / radius, z = xyzCell[i*3+2] / radius;
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
    const int i0 = vertexIndex[verticesOnEdge[i*2]], i1 = vertexIndex[verticesOnEdge[i*2+1]];
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
