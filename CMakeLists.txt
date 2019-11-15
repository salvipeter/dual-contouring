cmake_minimum_required(VERSION 3.1.0)

project(dualcontour)

if(CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_CXX_FLAGS "-Wall -pedantic -std=c++17 -fPIC")
  set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -DDEBUG")
  set(CMAKE_CXX_FLAGS_RELEASE "-O3")
endif(CMAKE_COMPILER_IS_GNUCXX)

set(Transfinite_DIR "../transfinite/src") # https://github.com/salvipeter/transfinite
include_directories("${Transfinite_DIR}/geom")
link_directories("${Transfinite_DIR}/geom/Release" "${Transfinite_DIR}/../release/geom")

add_library(dualcontour STATIC dc.cc quadmesh.cc)

target_link_libraries(dualcontour "geom")
