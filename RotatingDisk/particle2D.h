#pragma once
#include "vect.h"
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

class BD_2 {
public:
  BD_2() : pos(), force() {}
  BD_2(const Vec_2<double>& pos0) : pos(pos0), force() {}

  template <typename TRan>
  BD_2(TRan& myran, const Vec_2<double>& l, const Vec_2<double>& o);

  double get_theta() const { return 0.; }

  Vec_2<double> pos;
  Vec_2<double> force;
};


template<typename TRan>
BD_2::BD_2(TRan& myran, const Vec_2<double>& l, const Vec_2<double>& o) {
  pos.x = myran.doub() * l.x + o.x;
  pos.y = myran.doub() * l.y + o.y;
  force.x = 0;
  force.y = 0;
}


