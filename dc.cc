#include "dc.hh"

namespace DualContouring {

const double SMALL_NUMBER = 1e-15;

bool computeCell(const Point3D &origin, const std::array<Vector3D, 3> &dirs,
                 const std::array<double, 8> &vertices, Point3D &surface_point) {
  bool negative = false, positive = false;
  for (size_t i = 0; i < 8; ++i)
    if (vertices[i] < 0)
      negative = true;
    else
      positive = true;
  if (!negative || !positive)
    return false;

  // Find all edge intersections, and take their mean
  static constexpr std::array<size_t, 24> edges =
    { 0, 4, 1, 5, 2, 6, 3, 7, 0, 2, 1, 3, 4, 6, 5, 7, 0, 1, 2, 3, 4, 5, 6, 7 };

  Point3D mean(0, 0, 0);
  size_t count = 0;
  for (size_t i = 0; i < 12; ++i) {
    size_t v1 = edges[2*i], v2 = edges[2*i+1];
    double a = vertices[v1], b = vertices[v2];
    if (a * b > 0)
      continue;
    double denom = std::abs(b - a);
    double x = denom < SMALL_NUMBER ? 0.5 : std::abs(a) / denom;
    mean += dirs[0] * (v1 / 4) + dirs[1] * ((v1 % 4) / 2) + dirs[2] * (v1 % 2) + dirs[i/4] * x;
    ++count;
  }
  surface_point = origin + mean / count;
  return true;
}

std::vector<size_t> addPoints(QuadMesh &mesh, const std::vector<double> &values,
                              const Point3D &corner, const Vector3D &delta,
                              const std::array<size_t, 3> &resolution) {
  std::vector<size_t> cells;
  cells.reserve(resolution[0] * resolution[1] * resolution[2]);

  std::array<Vector3D, 3> dirs =
    { { { delta[0], 0, 0 },
        { 0, delta[1], 0 },
        { 0, 0, delta[2] } } };

  size_t mi = (resolution[1] + 1) * (resolution[2] + 1), mj = (resolution[2] + 1);
  size_t point_index = 1;
  std::array<double, 8> vertices;
  for (size_t i = 0; i < resolution[0]; ++i) {
    size_t index1 = i * mi;
    for (size_t j = 0; j < resolution[1]; ++j) {
      size_t index2 = index1 + j * mj;
      for (size_t k = 0; k < resolution[2]; ++k) {
        size_t index = index2 + k;

        // Gather the distance values for cell (i,j,k)
        for (size_t di = 0, vi = 0; di <= 1; ++di)
          for (size_t dj = 0; dj <= 1; ++dj)
            for (size_t dk = 0; dk <= 1; ++dk, ++vi)
              vertices[vi] = values[index+di*mi+dj*mj+dk];

        Point3D origin = corner + Vector3D(delta[0] * i, delta[1] * j, delta[2] * k);
        Point3D surface_point;
        if (computeCell(origin, dirs, vertices, surface_point)) {
          mesh.addPoint(surface_point);
          cells.push_back(point_index++);
        } else
          cells.push_back(0);
      }
    }
  }

  return cells;
}

void addQuads(QuadMesh &mesh, const std::vector<double> &values, const std::vector<size_t> &cells,
              const std::array<size_t, 3> &resolution) {
  static constexpr std::array<size_t, 3> c1 = { 1, 2, 0 }, c2 = { 2, 0, 1 };
  std::array<size_t, 3> ns = { resolution[1] * resolution[2], resolution[2], 1 };
  std::array<size_t, 3> ms = { (resolution[1] + 1) * (resolution[2] + 1), resolution[2] + 1, 1 };
  for (size_t c = 0; c < 3; ++c) {
    size_t ni = ns[c], nj = ns[c1[c]], nk = ns[c2[c]];
    size_t mi = ms[c], mj = ms[c1[c]], mk = ms[c2[c]];
    for (size_t i = 0; i < resolution[c]; ++i) {
      for (size_t j = 1; j < resolution[c1[c]]; ++j) {
        for (size_t k = 1; k < resolution[c2[c]]; ++k) {
          size_t index = i * ni + j * nj + k * nk;
          size_t a = cells[index], b = cells[index-nj], c = cells[index-nj-nk], d = cells[index-nk];
          if (a * b * c * d == 0)
            continue;
          size_t vi = i * mi + j * mj + k * mk;
          double v1 = values[vi], v2 = values[vi+mi];
          if (v1 * v2 > 0)
            continue;
          if (v1 < 0)
            mesh.addQuad(a, b, c, d);
          else
            mesh.addQuad(d, c, b, a);
        }
      }
    }
  }
}

QuadMesh isosurface(std::function<double(const Point3D &)> f, double isolevel,
                    const std::array<Point3D, 2> &bounding_box,
                    const std::array<size_t, 3> &resolution) {
  auto axis = bounding_box[1] - bounding_box[0];
  Vector3D delta(axis[0] / resolution[0], axis[1] / resolution[1], axis[2] / resolution[2]);

  // Compute distances at the cell corners
  std::vector<double> values;
  values.reserve((resolution[0] + 1) * (resolution[1] + 1) * (resolution[2] + 1));
  for (size_t i = 0; i <= resolution[0]; ++i) {
    double x = delta[0] * i;
    for (size_t j = 0; j <= resolution[1]; ++j) {
      double y = delta[1] * j;
      for (size_t k = 0; k <= resolution[2]; ++k) {
        double z = delta[2] * k;
        double v = f(bounding_box[0] + Vector3D(x, y, z)) - isolevel;
        values.push_back(std::abs(v) < SMALL_NUMBER ? (v < 0 ? -SMALL_NUMBER : SMALL_NUMBER) : v);
      }
    }
  }

  QuadMesh mesh;
  auto cells = addPoints(mesh, values, bounding_box[0], delta, resolution);
  addQuads(mesh, values, cells, resolution);

  return mesh;
}

QuadMesh isosurface(const std::vector<double> &distances, double isolevel,
                    const std::array<Point3D, 2> &bounding_box,
                    const std::array<size_t, 3> &resolution) {
  auto axis = bounding_box[1] - bounding_box[0];
  Vector3D delta(axis[0] / resolution[0], axis[1] / resolution[1], axis[2] / resolution[2]);

  std::vector<double> values = distances;
  for (double &v : values) {
    v -= isolevel;
    if (std::abs(v) < SMALL_NUMBER)
      v = v < 0 ? -SMALL_NUMBER : SMALL_NUMBER;
  }

  QuadMesh mesh;
  auto cells = addPoints(mesh, values, bounding_box[0], delta, resolution);
  addQuads(mesh, values, cells, resolution);

  return mesh;
}

}
