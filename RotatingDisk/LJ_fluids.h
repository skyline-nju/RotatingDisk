#pragma once

#include "ini.h"


void run_LJ(double Lx, double Ly, double rho0, double eps,
            double h, int n_step, double snap_log_sep,
            int seed, const std::string& ini_mode, char folder[]);



void run_LJ_hWall(double Lx, double Ly, double rho0, double eps, double eps_W,
                 double h, int n_step, double snap_log_sep,
                 int seed, const std::string& ini_mode, char folder[]);