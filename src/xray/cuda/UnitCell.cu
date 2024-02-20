#include "UnitCell.h"

xray::UnitCell::UnitCell(double a, double b, double c, double alpha_deg,
                         double beta_deg, double gamma_deg)
    : m_a(a), m_b(b), m_c(c), m_alpha_deg(alpha_deg), m_beta_deg(beta_deg),
      m_gamma_deg(gamma_deg) {
  bool is_cubic = (alpha_deg == 90 && beta_deg == 90 && gamma_deg == 90);

  double cos_alpha = 0;
  double cos_beta = 0;
  double cos_gamma = 0;

  if (is_cubic) {
    m_frac_to_orth = {{{a, 0., 0.}, {0., b, 0.}, {0., 0., c}}};

    m_orth_to_frac = {{{1 / a, 0., 0.}, {0., 1 / b, 0.}, {0., 0., 1 / c}}};
  } else {
    cos_alpha = std::cos(to_radians(alpha_deg));
    cos_beta = std::cos(to_radians(beta_deg));
    cos_gamma = std::cos(to_radians(gamma_deg));
    double sin_alpha = std::sin(to_radians(alpha_deg));
    double sin_beta = std::sin(to_radians(beta_deg));
    double sin_gamma = std::sin(to_radians(gamma_deg));
    double cos_alphar_sin_beta = (cos_beta * cos_gamma - cos_alpha) / sin_gamma;
    double cos_alphar = cos_alphar_sin_beta / sin_beta;
    //      double cos_betar =
    //          (cos_alpha * cos_gamma - cos_beta) / (sin_alpha * sin_gamma);
    //      double cos_gammar =
    //          (cos_alpha * cos_beta - cos_gamma) / (sin_alpha * sin_beta);
    double sin_alphar = std::sqrt(1.0 - cos_alphar * cos_alphar);
    m_frac_to_orth = {{
        {a, b * cos_gamma, c * cos_beta},
        {0., b * sin_gamma, -c * cos_alphar_sin_beta},
        {0., 0., c * sin_beta * sin_alphar},

    }};
    double o12 = -cos_gamma / (sin_gamma * a);
    double o13 = -(cos_gamma * cos_alphar_sin_beta + cos_beta * sin_gamma) /
                 (sin_alphar * sin_beta * sin_gamma * a);
    double o23 = cos_alphar / (sin_alphar * sin_gamma * b);
    m_orth_to_frac = {{{1 / a, o12, o13},
                       {0., 1 / m_frac_to_orth[1][1], o23},
                       {0., 0., 1 / m_frac_to_orth[2][2]}}};
  }
  m_metric_tensor = {a * a,
                     b * b,
                     c * c,
                     a * b * cos_gamma,
                     a * c * cos_beta,
                     b * c * cos_alpha};

  m_volume = a * b * c * std::sqrt(1 - cos_alpha * cos_alpha
                                   - cos_beta * cos_beta - cos_gamma * cos_gamma
                                   + 2 * cos_alpha * cos_beta * cos_gamma);
}
