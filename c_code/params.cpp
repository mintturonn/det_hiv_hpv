#include <stdio.h>
#include <stdlib.h>
#include "params.h"

double* read_array(FILE* f) {
  int length;
  fread(&length, 1, 4, f);
  double* arr = new double[length + 1];

  // It will be easier if the arrays line up with matlab,
  // which starts at index 1, so...

  arr[0] = -999;

  for (int i = 1; i <= length; i++) {
    fread(&arr[i], 1, 8, f);
  }
  return arr;
}

params::params(char* params_file, char* scenario_file) {
  threads = 1;

  FILE* f = fopen(params_file, "rb");
  if (!f) {
    printf("Parameter File not found, %s\n", params_file);
    exit(-1);
  }

  bp0 = read_array(f);
  dp0 = read_array(f);
  hp0 = read_array(f);
  pp0 = read_array(f);

  fclose(f);


  vaccscen_setup = 0;
  vaccscen_agef1 = 0;
  vaccscen_agef2hp = 0;
  vaccscen_agef3hp = 0;
  vaccscen_agef4hp = 0;


  screenscen_setup = 0;
  screenscen_scen_15to24 = CREATE_DARRAY(4, screenscen_scen_15to24, 0.0);
  screenscen_scen_25to34 = CREATE_DARRAY(4, screenscen_scen_25to34, 0.0);
  screenscen_scen_35to49 = CREATE_DARRAY(4, screenscen_scen_35to49, 0.0);
  screenscen_scen_15to24hn = CREATE_DARRAY(4, screenscen_scen_15to24hn, 0.0);
  screenscen_scen_25to34hn = CREATE_DARRAY(4, screenscen_scen_25to34hn, 0.0);
  screenscen_scen_35to49hn = CREATE_DARRAY(4, screenscen_scen_35to49hn, 0.0);
  screenscen_scen_15to24hp = CREATE_DARRAY(4, screenscen_scen_15to24hp, 0.0);
  screenscen_scen_25to34hp = CREATE_DARRAY(4, screenscen_scen_25to34hp, 0.0);
  screenscen_scen_35to49hp = CREATE_DARRAY(4, screenscen_scen_35to49hp, 0.0);

  screenscen_yr1 = 0;
  screenscen_yr1 = 1;
  screenscen_yr1 = 2;

  txsuxx_cinhn = 0.0;
  txsuxx_cinhp = 0.0;
  txsuxx_cchn = 0.0;
  txsuxx_cchp = 0.0;
  screenscen_cintx = CREATE_DARRAY(2, screenscen_cintx, 0.0)
  screenscen_cctx = CREATE_DARRAY(2, screenscen_cctx, 0.0)

  f = fopen(scenario_file, "rb");
  if (!f) {
    printf("Scenario File not found, %s\n", scenario_file);
    exit(-1);
  }

  fread(&vaccscen_setup, 4, 1, f);
  if (vaccscen_setup >= 1) {
    fread(&vaccscen_agef1, 8, 1, f);
    fread(&vaccscen_agef2hp, 8, 1, f);
    fread(&vaccscen_agef3hp, 8, 1, f);
    fread(&vaccscen_agef4hp, 8, 1, f);
  }

  fread(&screenscen_setup, 4, 1, f);
  if (screenscen_setup > 0) {
    fread(&screenscen_yr1, 4, 1, f);
    fread(&screenscen_yr2, 4, 1, f);
    fread(&screenscen_yr3, 4, 1, f);
    fread(&screenscen_cintx[1], 8, 2, f);
    fread(&screenscen_cctx[1], 8, 2, f);
  }
  if (screenscen_setup == 1) {
    for (int si = 1; si <= 3; si++) {
      fread(&screenscen_scen_15to24[si], 8, 1, f);
      fread(&screenscen_scen_25to34[si], 8, 1, f);
      fread(&screenscen_scen_35to49[si], 8, 1, f);
    }
  } else if (screenscen_setup == 2) {
    for (int si = 1; si <= 3; si++) {
      fread(&screenscen_scen_15to24hn[si], 8, 1, f);
      fread(&screenscen_scen_25to34hn[si], 8, 1, f);
      fread(&screenscen_scen_35to49hn[si], 8, 1, f);
      fread(&screenscen_scen_15to24hp[si], 8, 1, f);
      fread(&screenscen_scen_25to34hp[si], 8, 1, f);
      fread(&screenscen_scen_35to49hp[si], 8, 1, f);
    }
  }

  fclose(f);
}
