#pragma once

#include "domain2D.h"
#include "cellList2D.h"
#include "force2D.h"
#include "integrate2D.h"
#include "exporter2D.h"
#include "comn.h"
#include "particle2D.h"

template <typename TRan, typename TPar>
void ini_rand(std::vector<TPar>& p_arr, TRan& myran, int n_par, const Vec_2<double>& gl_l) {
  p_arr.reserve(n_par);
  for (size_t i = 0; i < n_par; i++) {
    p_arr.emplace_back(myran, gl_l, Vec_2<double>());
  }
}


template <typename TRan, typename TPar>
void ini_rand_wo_overlap(std::vector<TPar>& p_arr, TRan& myran,
                         int n_par, const Vec_2<double>& gl_l,
                         double sigma = 1.) {
  typedef BiNode<BD_2> node_t;
  ini_rand(p_arr, myran, n_par, gl_l);

  double Dt = 0.;
  double k = 10;
  double h = 0.05;
  int n_step = 10000;

  double r_cut = sigma;
  Grid_2 grid(gl_l, r_cut);
  PeriodicDomain_2 pdm(gl_l);
  CellListNode_2<node_t> cl(pdm, grid);

  // ini integrator
  BrownianDynamicsEM integrator(h, Dt);

  // cal force
  SpringForce kernal(r_cut, k, sigma);
  auto f1 = [&kernal](node_t* p1, node_t* p2) {
    kernal(*p1, *p2);
  };
  auto f2 = [&kernal, &pdm](node_t* p1, node_t* p2) {
    kernal(*p1, *p2, pdm);
  };

  cl.create(p_arr);
  for (int t = 1; t <= n_step; t++) {
    cl.for_each_pair(f1, f2);
    for (int i = 0; i < n_par; i++) {
      integrator.update(p_arr[i], pdm, myran);
    }
    cl.recreate(p_arr);
  }
}

template <typename TPar, typename TSnap>
void ini_from_snap(std::vector<TPar>& p_arr, TSnap& snap) {
  snap.read_last_frame(p_arr);
}


template <typename TPar, typename TRan, typename TSnap>
void ini(const std::string& ini_mode, std::vector<TPar>& p_arr, TRan& myran,
  int n_par, const Vec_2<double>& gl_l, TSnap& snap,
  double sigma = 1.) {
  if (ini_mode == "rand") {
    ini_rand_wo_overlap(p_arr, myran, n_par, gl_l, sigma);
  } else if (ini_mode == "resume") {
    ini_from_snap(p_arr, snap);
  } else {
    std::cout << "Error, ini_mode must be one of rand or resume" << std::endl;
    exit(1);
  }
}

