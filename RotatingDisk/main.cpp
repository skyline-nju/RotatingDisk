#include <iostream>
#include "LJ_fluids.h"

int main(int argc, char* argv[]) {
  double Lx = 100;
  double Ly = 100;

  double rho0 = 0.25;
  double Dt = 1;
  double eps = 0.5;
  double eps_W = 1;
  double h = 0.0025;
  int n_step = 1000000;
  //int snap_interval = int(round(200 / h * 0.1));
  double snap_log_sep = 10000;

  int seed = 1234;
  std::string ini_mode = "resume";  // should be "rand" or "resume"

#ifdef _MSC_VER
  char folder[] = "data\\";
#else
  char folder[] = "/home/sli/yduan/ChiralBPs/LJwall/";
#endif

  run_LJ_hWall(Lx, Ly, rho0, eps, eps_W, 
               h, n_step, snap_log_sep, seed,
               ini_mode, folder);
}
