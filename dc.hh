#pragma once

#include <functional>

#include <geometry.hh>

#include "quadmesh.hh"

// Constructs a quad mesh approximating the isosurface of a given scalar function

QuadMesh
dualContouring(std::function<double(const Geometry::Point3D &)> f,
               double isolevel, 
               const std::array<Geometry::Point3D, 2> &bounding_box, 
               const std::array<size_t, 3> &resolution);
