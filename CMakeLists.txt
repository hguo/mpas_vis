cmake_minimum_required (VERSION 3.1.3)
project (mpas_vis)

set (CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_STANDARD_REQUIRED ON)

find_package (NetCDF REQUIRED)
include_directories (${NETCDF_INCLUDE_DIR})

find_package (GLEW REQUIRED)
include_directories (${GLEW_INCLUDE_DIR})

find_package (OpenGL REQUIRED)
find_package (Qt4 REQUIRED COMPONENTS QtCore QtGui QtOpenGL)
include (${QT_USE_FILE})

add_executable (main main.cpp widget.cpp trackball.cpp)
set_target_properties (main PROPERTIES AUTOMOC ON)
target_link_libraries (main ${QT_LIBRARIES} ${OPENGL_LIBRARIES} ${GLEW_LIBRARY} ${NETCDF_LIBRARIES})
