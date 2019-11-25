#pragma once

#include <functional>

#include "quadmesh.hh"

// Constructs a quad mesh approximating the isosurface of a given scalar function

namespace DualContouring {

QuadMesh isosurface(std::function<double(const Point3D &)> f, double isolevel,
                    const std::array<Point3D, 2> &bounding_box,
                    const std::array<size_t, 3> &resolution);

QuadMesh isosurface(const std::vector<double> &distances, double isolevel,
                    const std::array<Point3D, 2> &bounding_box,
                    const std::array<size_t, 3> &resolution);

}
