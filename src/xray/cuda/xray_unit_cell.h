#ifndef AMBER_XRAY_UNIT_CELL_H
#define AMBER_XRAY_UNIT_CELL_H

#include <array>
#include <cmath>

namespace {

[[nodiscard]] inline double to_radians(double degrees) {
  return degrees * M_PI / 180.0;
}

[[nodiscard]] inline double to_degrees(double radians) {
  return radians / M_PI * 180.0;
}
} // namespace

namespace xray {

struct Sym33 {
  double xx;
  double yy;
  double zz;
  double xy;
  double xz;
  double yz;
};

class UnitCell {
public:
  UnitCell() = default;
  UnitCell(const UnitCell &) = default;
  UnitCell &operator=(const UnitCell &) = default;

  UnitCell(double a, double b, double c, double alpha_deg, double beta_deg,
           double gamma_deg);

  [[nodiscard]] inline double to_squared_orth_norm(double fx, double fy,
                                                   double fz) const {
    const auto &m = m_metric_tensor;
    return m.xx * fx * fx + m.yy * fy * fy + m.zz * fz * fz +
           2 * (m.xy * fx * fy + m.xz * fx * fz + m.yz * fy * fz);
  }

  [[nodiscard]] inline std::array<double, 3> to_orth(double fx, double fy,
                                                     double fz) const {
    auto &m = m_frac_to_orth;
    return {
        m[0][0] * fx + m[0][1] * fy + m[0][2] * fz,
        m[1][0] * fx + m[1][1] * fy + m[1][2] * fz,
        m[2][0] * fx + m[2][1] * fy + m[2][2] * fz,
    };
  }
  [[nodiscard]] inline std::array<double, 3> to_frac(double x, double y,
                                                     double z) const {
    auto &m = m_orth_to_frac;
    return {
        m[0][0] * x + m[0][1] * y + m[0][2] * z,
        m[1][0] * x + m[1][1] * y + m[1][2] * z,
        m[2][0] * x + m[2][1] * y + m[2][2] * z,
    };
  }

  [[nodiscard]] inline const Sym33 &metric_tensor() const {
    return m_metric_tensor;
  }

  [[nodiscard]] inline double get_volume() const {
    return m_volume;
  }

private:
  double m_a = 0;
  double m_b = 0;
  double m_c = 0;
  double m_alpha_deg = 0;
  double m_beta_deg = 0;
  double m_gamma_deg = 0;

  std::array<std::array<double, 3>, 3> m_frac_to_orth = {};
  std::array<std::array<double, 3>, 3> m_orth_to_frac = {};
  Sym33 m_metric_tensor = Sym33{0, 0, 0, 0, 0, 0};
  double m_volume = 0;
};

} // namespace xray

#endif // AMBER_XRAY_UNIT_CELL_H
