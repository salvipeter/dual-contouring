#pragma once

#include "geometry.hh"

struct QuadMesh {
  using Quad = std::array<size_t, 4>;
  Geometry::PointVector points;
  std::vector<Quad> quads;

  void addPoint(const Geometry::Point3D &p);
  void addQuad(size_t a, size_t b, size_t c, size_t d);
  void writeOBJ(std::string filename) const;
};
