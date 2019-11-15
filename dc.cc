#include "dc.hh"

using namespace Geometry;

namespace {

  bool computeCell(const Point3D &origin, const std::array<Vector3D, 3> &delta,
                   const std::array<double, 8> &vertices, Point3D &surface_point) {
    bool negative = false, positive = false;
    for (size_t i = 0; i < 8; ++i)
      if (vertices[i] < 0)
        negative = true;
      else if (vertices[i] > 0)
        positive = true;
    if (!negative || !positive)
      return false;

    // Placeholder: just take the center point
    surface_point = origin + delta[0] * 0.5 + delta[1] * 0.5 + delta[2] * 0.5;
    return true;
  }

}

QuadMesh
dualContouring(std::function<double(const Point3D &)> f,
               double isolevel, 
               const std::array<Point3D, 2> &bounding_box, 
               const std::array<size_t, 3> &resolution) {
  QuadMesh mesh;
  std::vector<size_t> cells;
  cells.reserve(resolution[0] * resolution[1] * resolution[2]);

  auto axis = bounding_box[1] - bounding_box[0];
  std::array<Vector3D, 3> delta = { {
    { axis[0] / resolution[0], 0, 0 },
    { 0, axis[1] / resolution[1], 0 },
    { 0, 0, axis[2] / resolution[2] }
  } };

  // Compute distances at the cell corners
  std::vector<double> values;
  values.reserve((resolution[0] + 1) * (resolution[1] + 1) * (resolution[2] + 1));
  for (size_t i = 0; i <= resolution[0]; ++i) {
    double u = (double)i / resolution[0];
    for (size_t j = 0; j <= resolution[1]; ++j) {
      double v = (double)j / resolution[1];
      for (size_t k = 0; k <= resolution[2]; ++k) {
        double w = (double)k / resolution[2];
        values.push_back(f({u, v, w}));
      }
    }
  }

  // Compute representative points for each cell
  size_t ni = (resolution[1] + 1) * (resolution[2] + 1), nj = (resolution[2] + 1);
  size_t point_index = 1;
  std::array<double, 8> vertices;
  for (size_t i = 0; i < resolution[0]; ++i) {
    size_t index1 = i * ni;
    for (size_t j = 0; j < resolution[1]; ++j) {
      size_t index2 = index1 + j * nj;
      for (size_t k = 0; k < resolution[2]; ++k) {
        size_t index = index2 + k;

        // Gather the distance values for cell (i,j,k)
        for (size_t di = 0, vi = 0; di <= 1; ++di)
          for (size_t dj = 0; dj <= 1; ++dj)
            for (size_t dk = 0; dk <= 1; ++dk, ++vi)
              vertices[vi] = values[index+di*ni+dj*nj+dk];

        Point3D origin = bounding_box[0] + delta[0] * i + delta[1] * j + delta[2] * k;
        Point3D surface_point;
        bool found = computeCell(origin, delta, vertices, surface_point);
        if (found) {
          mesh.addPoint(surface_point);
          cells.push_back(point_index++);
        } else
          cells.push_back(0);
      }
    }
  }

  // Add the quads
  constexpr std::array<size_t, 3> c1 = { 1, 2, 0 }, c2 = { 2, 0, 1 };
  std::array<size_t, 3> ns = { resolution[1] * resolution[2], resolution[2], 1 };
  for (size_t c = 0; c < 3; ++c) {
    size_t ni = ns[c], nj = ns[c1[c]], nk = ns[c2[c]];
    for (size_t i = 0; i < resolution[c]; ++i) {
      for (size_t j = 1; j < resolution[c1[c]]; ++j) {
        for (size_t k = 1; k < resolution[c2[c]]; ++k) {
          size_t index = i * ni + j * nj + k * nk;
          size_t a = cells[index], b = cells[index-nj], c = cells[index-nj-nk], d = cells[index-nk];
          if (a * b * c * d != 0)
            mesh.addQuad(a, b, c, d);
        }
      }
    }
  }

  return mesh;
}
