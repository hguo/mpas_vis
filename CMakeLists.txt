cmake_minimum_required (VERSION 3.1.3)
project (mpas_vis)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
  
find_package (OpenGL REQUIRED)
find_package (Qt4 REQUIRED COMPONENTS QtCore QtGui QtOpenGL)
include (${QT_USE_FILE})

add_executable (main main.cpp widget.cpp)
set_target_properties (main PROPERTIES AUTOMOC ON)
target_link_libraries (main ${QT_LIBRARIES} ${OPENGL_LIBRARIES})
