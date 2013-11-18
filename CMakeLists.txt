cmake_minimum_required (VERSION 2.6)
project (segmentation)

cmake_policy(SET CMP0015 NEW)
set(TARGET_NAME ${PROJECT_NAME})
set(DEPS_PKGCONFIG "eigen3 gsl opencv")
find_package(PkgConfig)


set (SOURCE
    src/regionGrowing.cpp
    src/kdtree.cpp
    src/plane.cpp
    src/hull.cpp
    src/spatialObject.cpp
)

FILE(GLOB HEADERS "./include/segmentation/*.h" "./include/segmentation/*.hpp")

include_directories(${CMAKE_SOURCE_DIR}/include)

add_definitions(
	-DTRANSPARENCY=64
	-std=gnu++0x	
	-frounding-math
)

pkg_check_modules(Eigen REQUIRED eigen3)
include_directories(${Eigen_INCLUDE_DIRS})

#need CGAL plus 3rd_PARTY_LIBRARIES!!! 
find_package(CGAL REQUIRED)
include_directories(${CGAL_INCLUDE_DIRS})
add_definitions ( ${CGAL_3RD_PARTY_DEFINITIONS}  ${CGAL_DEFINITIONS}  )
link_directories ( ${CGAL_LIBRARIES_DIR} ${CGAL_3RD_PARTY_LIBRARIES_DIRS} )

find_package( OpenCV REQUIRED )

pkg_check_modules(GSL gsl REQUIRED)
link_directories(GSL_LIBRARY_DIRS)

ADD_LIBRARY(${PROJECT_NAME} SHARED
  ${SOURCE}
)

target_link_libraries (segmentation ${CGAL_3RD_PARTY_LIBRARIES} ${CGAL_LIBRARY} ${GSL_LIBRARIES} ${OpenCV_LIBS})

CONFIGURE_FILE(./src/segmentation.pc.in segmentation.pc @ONLY)

INSTALL(FILES ${HEADERS} DESTINATION include/${PROJECT_NAME})
INSTALL(TARGETS segmentation DESTINATION lib)
INSTALL(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/segmentation.pc
  DESTINATION lib/pkgconfig)

