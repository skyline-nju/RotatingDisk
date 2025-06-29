#pragma once
#include "vect.h"
#include <cmath>
#include "config.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


// F(r) = k * (r - sigma)
class SpringForce {
public:
  SpringForce(double r_cut, double k, double sigma = 1.)
    : r_cut_(r_cut), r_cut_square_(r_cut* r_cut), k_(k), sigma_(sigma) {}

  template <typename TPar>
  void interact(TPar& p1, TPar& p2, const Vec_2<double>& r12_vec, double r12_square) const;

  template <typename TPar>
  void operator ()(TPar& p1, TPar& p2) const;

  template <typename TPar, typename BdyCondi>
  void operator ()(TPar& p1, TPar& p2, const BdyCondi& bc) const;

private:
  double r_cut_;
  double r_cut_square_;
  double k_;
  double sigma_;
};


template<typename TPar>
void SpringForce::interact(TPar& p1, TPar& p2, const Vec_2<double>& r12_vec, double r12_square) const {
  double r12 = std::sqrt(r12_square);
  Vec_2<double> r12_hat = r12_vec / r12;
  Vec_2<double> f12 = k_ * (r12 - sigma_) * r12_hat;
  p1.force += f12;
  p2.force -= f12;
}

template<typename TPar>
void SpringForce::operator()(TPar& p1, TPar& p2) const {
  Vec_2<double> r12 = p2.pos - p1.pos;
  double r12_square = r12.square();
  if (r12_square < r_cut_square_) {
    interact(p1, p2, r12, r12_square);
  }
}


template<typename TPar, typename BdyCondi>
void SpringForce::operator()(TPar& p1, TPar& p2, const BdyCondi& bc) const {
  Vec_2<double> r12 = p2.pos - p1.pos;
  bc.untangle(r12);
  double r12_square = r12.square();
  if (r12_square < r_cut_square_) {
    interact(p1, p2, r12, r12_square);
  }
}



// U(r) = eps * [(sigma/r)^12 - (sigma/r)^6]
// F(r) = 6eps * [2 * (sigma/r)^12 - (sigma/r)^6] / r
class LJ_Force {
public:
  LJ_Force(double r_cut, double eps = 1., double sigma = 1.)
    : r_cut_(r_cut), r_cut_square_(r_cut* r_cut), eps_6_(eps * 6), sigma_(sigma) {}

  template <typename TPar>
  void interact(TPar& p1, TPar& p2, const Vec_2<double>& r12_vec, double r12_square) const;

  template <typename TPar>
  void operator ()(TPar& p1, TPar& p2) const;

  template <typename TPar, typename BdyCondi>
  void operator ()(TPar& p1, TPar& p2, const BdyCondi& bc) const;

private:
  double r_cut_;
  double r_cut_square_;
  double eps_6_;
  double sigma_;
};


template<typename TPar>
void LJ_Force::interact(TPar& p1, TPar& p2, const Vec_2<double>& r12_vec, double r12_square) const {
  double r12 = std::sqrt(r12_square);
  double inverse_r = 1. / r12;
  double pow6 = std::pow(sigma_ * inverse_r, 6);
  double pow12 = pow6 * pow6;
  Vec_2<double> f12 = eps_6_ * (2 * pow12 - pow6) * inverse_r * inverse_r * r12_vec;
  p1.force -= f12;
  p2.force += f12;
}


template<typename TPar>
void LJ_Force::operator()(TPar& p1, TPar& p2) const {
  Vec_2<double> r12 = p2.pos - p1.pos;
  double r12_square = r12.square();
  if (r12_square < r_cut_square_) {
    interact(p1, p2, r12, r12_square);
  }
}


template<typename TPar, typename BdyCondi>
void LJ_Force::operator()(TPar& p1, TPar& p2, const BdyCondi& bc) const {
  Vec_2<double> r12 = p2.pos - p1.pos;
  bc.untangle(r12);
  double r12_square = r12.square();
  if (r12_square < r_cut_square_) {
    interact(p1, p2, r12, r12_square);
  }
}

class LJ_hWalls {
public:
  LJ_hWalls(double Ly, bool flag_attractive_lower, bool flag_attractive_upper,
            double eps, double sigma=1.0);

  template <typename TPar>
  void interact(TPar &p) const;

  double get_LJ_force(double dx) const;

private:
  double Ly_;
  double r_cut_[2];
  double eps_6_;
  double sigma_;
  double half_sigma_;
};

template<typename TPar>
void LJ_hWalls::interact(TPar& p) const {
  double dy0 = p.pos.y + half_sigma_;
#ifdef DEBUG_BDY
    if (dy0 < 0.) {
      std::cout << "Error, y=" << p.pos.y << std::endl;
    }
#endif
  if (dy0 < r_cut_[0]) {
    double fy = get_LJ_force(dy0);
    p.force.y += fy;
  } else {
    double dy1 = Ly_ - p.pos.y + half_sigma_;
#ifdef DEBUG_BDY
    if (dy1 < 0.) {
      std::cout << "Error, y=" << p.pos.y << std::endl;
    }
#endif
    if (dy1 < r_cut_[1]) {
      p.force.y -= get_LJ_force(dy1);
    }
  }
}
