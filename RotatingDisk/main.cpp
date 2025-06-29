#include <iostream>
#include "LJ_fluids.h"

int main(int argc, char* argv[]) {
  double Lx = 150;
  double Ly = 300;

#ifdef _MSC_VER
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
  char folder[] = "data\\";
#else
  double rho0 = atof(argv[1]);
  double Dt = 1;
  double eps = atof(argv[2]);
  double eps_W = atof(argv[3]);
  double h = atof(argv[4]);
  int n_step = atoi(argv[5]);
  double snap_log_sep = atof(argv[6]);
  int seed = atoi(argv[7]);
  std::string ini_mode = argv[8];  // should be "rand" or "resume"
  char folder[] = "/home/sli/yduan/ChiralBPs/LJ_rep_wall/L150_300_debug/";
#endif

  run_LJ_hWall(Lx, Ly, rho0, eps, eps_W, 
               h, n_step, snap_log_sep, seed,
               ini_mode, folder);
}
