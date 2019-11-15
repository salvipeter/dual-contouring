#pragma once

#include <array>

namespace DualContouring {

  struct Vector3D {
    std::array<double, 3> data;

    Vector3D(double x = 0, double y = 0, double z = 0) : data({x, y, z}) { }
    Vector3D(const double *raw) : data({raw[0], raw[1], raw[2]}) { }
    double operator[](size_t i) const { return data[i]; }
    double &operator[](size_t i) { return data[i]; }
    Vector3D &operator+=(const Vector3D &v) {
      data[0] += v.data[0]; data[1] += v.data[1]; data[2] += v.data[2];
      return *this;
    }
    Vector3D operator+(const Vector3D &v) const {
      return { data[0] + v.data[0], data[1] + v.data[1], data[2] + v.data[2] };
    }
    Vector3D operator-(const Vector3D &v) const {
      return { data[0] - v.data[0], data[1] - v.data[1], data[2] - v.data[2] };
    }
    Vector3D operator*(double x) const {
      return { data[0] * x, data[1] * x, data[2] * x };
    }
    Vector3D operator/(double x) const { return operator*(1 / x); }
  };

  using Point3D = Vector3D;

}
