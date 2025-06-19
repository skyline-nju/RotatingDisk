#include <iostream>
#include "ini.h"

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
  std::string ini_mode = "rand";  // should be "rand" or "resume"

  int snap_interval = 1;
  if (snap_log_sep >= 1) {
    snap_interval = int(snap_log_sep);
    snap_log_sep = -1;
  }
  int n_par = int(rho0 * Lx * Ly);

  typedef BiNode<BD_2> node_t;
  Ranq2 myran(seed);
  Vec_2<double> gl_l(Lx, Ly);
  double r_cut = 2.5;
  Grid_2 grid(gl_l, r_cut);
  PeriodicDomain_2 pdm(gl_l);
  CellListNode_2<node_t> cl(pdm, grid);
  std::vector<node_t> p_arr;

  // ini integrator
  BrownianDynamicsEM integrator(h, Dt);

  // cal force
  LJ_Force kernal(r_cut, eps);
  //SpringForce kernal(r_cut, 1);
  auto f1 = [&kernal](node_t* p1, node_t* p2) {
    kernal(*p1, *p2);
  };
  auto f2 = [&kernal, &pdm](node_t* p1, node_t* p2) {
    kernal(*p1, *p2, pdm);
  };

  // set output
  char basename[255];
  char log_file[255];
  char op_file[255];
  char gsd_file[255];
#ifdef _MSC_VER
  char folder[] = "data\\";
#else
  char folder[] = "/home/sli/yduan/ChiralBPs/LJ/";
#endif

  char log_folder[255];
  snprintf(log_folder, 255, "%slog%s", folder, delimiter.c_str());
  mkdir(log_folder);

  snprintf(basename, 255, "L%g_%g_r%g_e%g_h%g_S%d",
    Lx, Ly, rho0, eps, h, seed);
  snprintf(gsd_file, 255, "%s%s.gsd", folder, basename);

  int start = 0;

  exporter::Snap_GSD_2 gsd(gsd_file, n_step, snap_interval, start, h, snap_log_sep, gl_l, ini_mode);

  int log_interval = 10000;
  int op_interval = 100;
  snprintf(log_file, 255, "%s%s_t%d.dat", log_folder, basename, start);
  exporter::LogExporter log(log_file, start, n_step, log_interval, n_par);

  // ini particles
  ini(ini_mode, p_arr, myran, n_par, gl_l, gsd, 1.);

  cl.create(p_arr);

  for (int t = 1; t <= n_step; t++) {
    //std::cout << "t = " << t << std::endl;
    cl.for_each_pair(f1, f2);
    for (int i = 0; i < n_par; i++) {
      integrator.update(p_arr[i], pdm, myran);
      //integrator.update_par_cellList(p_arr[i], pdm, myran, cl);
    }
    cl.recreate(p_arr);
    gsd.dump(t, p_arr);
    log.record(t);
  }



}
