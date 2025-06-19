#include "LJ_fluids.h"

void run_LJ(double Lx, double Ly, double rho0, double eps,
            double h, int n_step, double snap_log_sep,
            int seed, const std::string& ini_mode, char folder[]) {
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
  double Dt = 1.;
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


  char log_folder[255];
  snprintf(log_folder, 255, "%slog%s", folder, delimiter.c_str());
  mkdir(log_folder);

  snprintf(basename, 255, "L%g_%g_r%.4f_e%.4f_h%g_S%d",
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
      //integrator.update(p_arr[i], pdm, myran);
      integrator.update_par_cellList(p_arr[i], pdm, myran, cl);
    }
    //cl.recreate(p_arr);
    gsd.dump(t, p_arr);
    log.record(t);
  }
}


void run_LJ_hWall(double Lx, double Ly, double rho0, double eps, double eps_W,
                  double h, int n_step, double snap_log_sep, int seed,
                  const std::string& ini_mode, char folder[]) {
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
  Domain_w_hWalls<LJ_hWalls> pdm(gl_l, true, false, eps_W, 1.);
  CellListNode_2<node_t> cl(pdm, grid);
  std::vector<node_t> p_arr;

  // ini integrator
  double Dt = 1.;
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


  char log_folder[255];
  snprintf(log_folder, 255, "%slog%s", folder, delimiter.c_str());
  mkdir(log_folder);

  snprintf(basename, 255, "L%g_%g_r%.4f_e%.4f_%.4f_h%g_S%d",
    Lx, Ly, rho0, eps, eps_W, h, seed);
  snprintf(gsd_file, 255, "%s%s.gsd", folder, basename);

  int start = 0;

  exporter::Snap_GSD_2 gsd(gsd_file, n_step, snap_interval, start, h, snap_log_sep, gl_l, ini_mode);

  int log_interval = 10000;
  int op_interval = 100;
  snprintf(log_file, 255, "%s%s_t%d.dat", log_folder, basename, start);
  exporter::LogExporter log(log_file, start, n_step, log_interval, n_par);

  // ini particles
  ini(ini_mode, p_arr, myran, n_par, gl_l, gsd, 1., true);

  cl.create(p_arr);

  for (int t = 1; t <= n_step; t++) {
    //std::cout << "t = " << t << std::endl;
    cl.for_each_pair(f1, f2);
    for (int i = 0; i < n_par; i++) {
      //integrator.update(p_arr[i], pdm, myran);
      integrator.update_par_cellList(p_arr[i], pdm, myran, cl);
    }
    //cl.recreate(p_arr);
    gsd.dump(t, p_arr);
    log.record(t);
  }

}