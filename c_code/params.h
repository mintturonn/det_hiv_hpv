#ifndef PARAMS_H
#define PARAMS_H

#include "matrix_8d.h"

class params;

class params {
  public:
	  double* bp0;
	  double* dp0;
	  double* hp0;
	  double* pp0;
    int threads;

    #define DIM_h1 1
    #define DIM_h2 2
    #define DIM_h3 3
    #define DIM_i  4
    #define DIM_s  5
    #define DIM_r  6
    #define DIM_a  7
    #define DIM_v1 8

    int vaccscen_setup;
    double vaccscen_agef1 = 0.0;
    double vaccscen_agef2hp = 0.0;
    double vaccscen_agef3hp = 0.0;
    double vaccscen_agef4hp = 0.0;

    int screenscen_setup;
    int screenscen_yr1 = 0;
    int screenscen_yr2 = 0;
    int screenscen_yr3 = 0;
    DArray screenscen_scen_15to24;
    DArray screenscen_scen_25to34;
    DArray screenscen_scen_35to49;
    DArray screenscen_scen_15to24hn;
    DArray screenscen_scen_25to34hn;
    DArray screenscen_scen_35to49hn;
    DArray screenscen_scen_15to24hp;
    DArray screenscen_scen_25to34hp;
    DArray screenscen_scen_35to49hp;
    double txsuxx_cinhn;
    double txsuxx_cinhp;
    double txsuxx_cchn;
    double txsuxx_cchp;
    DArray screenscen_cintx;
    DArray screenscen_cctx;

    params(char* params_file, char* scenario_file);
};

#endif
