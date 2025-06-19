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

  template<class TPar, class TDomain, class TRan, class TCellList>
  void update_par_cellList(TPar& p, const TDomain& dm, TRan& myran, TCellList& cl) const;


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


template<class TPar, class TDomain, class TRan, class TCellList>
void BrownianDynamicsEM::update_par_cellList(TPar& p, const TDomain& dm,
                                             TRan& myran, TCellList& cl) const {
  auto ic_old = cl.get_ic(p);
  update(p, dm, myran);
  auto ic_new = cl.get_ic(p);
  if (ic_old != ic_new) {
    cl.update(p, ic_old, ic_new);
  }
}
