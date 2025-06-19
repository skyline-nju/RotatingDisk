#include <iostream>
#include "LJ_fluids.h"

int main(int argc, char* argv[]) {
  double Lx = 100;
  double Ly = 100;

  double rho0 = 0.25;
  double Dt = 1;
  double eps = 0.5;
  double h = 0.0025;
  int n_step = 100000;
  //int snap_interval = int(round(200 / h * 0.1));
  double snap_log_sep = 1000;

  int seed = 123;
  std::string ini_mode = "resume";  // should be "rand" or "resume"

#ifdef _MSC_VER
  char folder[] = "data\\";
#else
  char folder[] = "/home/sli/yduan/ChiralBPs/LJ/";
#endif

  run_LJ(Lx, Ly, rho0, eps, h, n_step, snap_log_sep, seed, ini_mode, folder);

}
