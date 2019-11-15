#include "quadmesh.hh"

#include <fstream>

namespace DualContouring {

void
QuadMesh::addPoint(const Point3D &p) {
  points.push_back(p);
}

void
QuadMesh::addQuad(size_t a, size_t b, size_t c, size_t d) {
  quads.push_back({a, b, c, d});
}

void
QuadMesh::writeOBJ(std::string filename) const {
  std::ofstream f(filename);
  f.exceptions(std::ios::failbit | std::ios::badbit);
  for (const auto &p : points)
    f << "v " << p[0] << ' ' << p[1] << ' ' << p[2] << std::endl;
  for (const auto &t : quads)
    f << "f " << t[0] << ' ' << t[1] << ' ' << t[2] << ' ' << t[3] << std::endl;
}

}
