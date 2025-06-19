#include "force2D.h"

LJ_hWalls::LJ_hWalls(double Ly, bool flag_attractive_lower, bool flag_attractive_upper,
                                    double eps, double sigma)
  : Ly_(Ly), eps_6_(eps * 6), sigma_(sigma), half_sigma_(0.5 * sigma) {
  if (flag_attractive_lower) {
    r_cut_[0] = 2.5 * sigma;
  } else {
    r_cut_[0] = std::pow(2., 1. / 6) * sigma;
  }

  if (flag_attractive_upper) {
    r_cut_[1] = 2.5 * sigma;
  } else {
    r_cut_[1] = std::pow(2., 1. / 6) * sigma;
  }
  std::cout << "sigma=" << sigma << "\tr_cut=" << r_cut_[0] << ", " << r_cut_[1] << std::endl;
}


double LJ_hWalls::get_LJ_force(double r) const {
  double inverse_r = 1. / r;
  double pow6 = std::pow(sigma_ * inverse_r, 6);
  double pow12 = pow6 * pow6;
  double f = eps_6_ * (2 * pow12 - pow6) * inverse_r;
  return f;
}
