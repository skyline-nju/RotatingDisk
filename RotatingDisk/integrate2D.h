#pragma once
#include <cmath>
#include "particle2D.h"
#include "node.h"


class BrownianDynamicsEM{
public:
  BrownianDynamicsEM(double h, double Dt)
    : h_(h), Dt_(24 * Dt * h) {}

  template <class TPar, class TDomain, class TRan>
  void update(TPar& p, const TDomain& dm, TRan& myran) const;

protected:
  double h_;
  double Dt_;
};

template<class TPar, class TDomain, class TRan>
void BrownianDynamicsEM::update(TPar& p, const TDomain& dm, TRan& myran) const {
  p.pos += p.force * h_;
  p.pos.x += (myran.doub() - 0.5) * Dt_;
  p.pos.y += (myran.doub() - 0.5) * Dt_;

  dm.tangle(p.pos);

  //if (p.pos.x >= dm.gl_l().x || p.pos.x < 0) {
  //  std::cout << p.force << std::endl;
  //}
  p.force.x = 0.;
  p.force.y = 0.;
}
