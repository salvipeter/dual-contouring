#pragma once

#include <vector>

#include "vector3d.hh"

namespace DualContouring {

struct QuadMesh {
  using Quad = std::array<size_t, 4>;
  std::vector<Point3D> points;
  std::vector<Quad> quads;

  void addPoint(const Point3D &p);
  void addQuad(size_t a, size_t b, size_t c, size_t d);
  void writeOBJ(std::string filename) const;
};

}
