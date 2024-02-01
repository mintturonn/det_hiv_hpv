#include "hpv_hiv_model_vacc.h"
#include <stdio.h>

/*
% 5 November 2019 Minttu R�nn
% South Africa model of HIV and HPV transmission and development of
% cervical cancer

% 3 HPV types - SIR with a waning immunity(3 comps), and cervical cancer development(2 comps)
% HIV infection - 4 compartment SI model
% Two sex model
% Three risk groups
% 13 age groups
% vaccination status(9 - valent)

% time - varying processes
% Population growth
% Perinatally infected intorduced as HIV + to 9 - 14 yos compartment
% Increase in ART use over time
% Increase in voluntary male circumcision(VMC) over time
% Increase in condom use over time
*/

void hiv_hpv_model_vacc(params* p, char* pop_file, char* output_file) {
  state* S = new state();

  /*
    h1 = dp0(1); % 5; % hpv1(16 / 18)
    h2 = dp0(2); % 5; % hvp2(other vaccine types)
    h3 = dp0(3); % 5; % hvp3 non - vaccine types
    i = dp0(4); % 4; % hiv
    s = dp0(5); % 2; % sex% F = 1, M = 2
    r = dp0(6); % 3; % risk group(3)
    a = dp0(7); % 5; % age; 1 = 9 - 14, 2 = 15 - 25, 3 = 26 - 35, 4 = 36 - 49, 5 = 50 - 70
    v1 = dp0(8); % 2; % vaccine; 1 = not 2 = vaccinated // when calibrating this is turned off
  */
  #define _h1 (int) round(p->dp0[DIM_h1])
  #define _h2 (int) round(p->dp0[DIM_h2])
  #define _h3 (int) round(p->dp0[DIM_h3])
  #define _i  (int) round(p->dp0[DIM_i])
  #define _s  (int) round(p->dp0[DIM_s])
  #define _r  (int) round(p->dp0[DIM_r])
  #define _a  (int) round(p->dp0[DIM_a])
  #define _v1 (int) round(p->dp0[DIM_v1])

/*
  [hivS_id, hiv_id, hiv2_id, acute_id, art_id, noart_id, vacc_id, hpv1_id, hpv2_id, ...
  hpv3_id, hpv1S_id, hpv2S_id, hpv3S_id, hpv12_id, hpv1R_id, hpv2R_id, hpv3R_id, ...
  hpv1I_id, hpv2I_id, hpv3I_id, hpv1NC_id, hpv2NC_id, hpv3NC_id, ...
  hpv1CC_id, hpv2CC_id, hpv3CC_id, hpv1CIN_id, hpv2CIN_id, hpv3CIN_id, CIN_id, CC_id,...
  sexf_id, sexm_id, risk_id, age_id, age_id2] = deal(zeros(h1, h2, h3, i, s, r, a, v1));
*/

  int* dims = new int[9] {0, _h1,_h2,_h3,_i, _s, _r, _a, _v1};

  S->hivS_id = cNEW_8D(S->hivS_id);
  S->hiv_id = cNEW_8D(S->hiv_id)
  S->hiv2_id = cNEW_8D(S->hiv2_id)
  cMAT_8D acute_id =   cNEW_8D(acute_id)
  S->art_id = cNEW_8D(S->art_id)
  S->noart_id = cNEW_8D(S->noart_id)
  S->vacc_id = cNEW_8D(S->vacc_id)
  S->hpv1_id = cNEW_8D(S->hpv1_id)
  S->hpv2_id = cNEW_8D(S->hpv2_id)
  S->hpv3_id = cNEW_8D(S->hpv3_id)
  S->hpv1S_id = cNEW_8D(S->hpv1S_id)
  S->hpv2S_id = cNEW_8D(S->hpv2S_id)
  S->hpv3S_id = cNEW_8D(S->hpv3S_id)
  S->hpv12_id = cNEW_8D(S->hpv12_id)
  S->hpv1R_id = cNEW_8D(S->hpv1R_id)
  S->hpv2R_id = cNEW_8D(S->hpv2R_id)
  S->hpv3R_id = cNEW_8D(S->hpv3R_id)
  S->hpv1I_id = cNEW_8D(S->hpv1I_id)
  S->hpv2I_id = cNEW_8D(S->hpv2I_id)
  S->hpv3I_id = cNEW_8D(S->hpv3I_id)
  cMAT_8D hpv1NC_id = cNEW_8D(hpv1NC_id)
  cMAT_8D hpv2NC_id = cNEW_8D(hpv2NC_id)
  cMAT_8D hpv3NC_id = cNEW_8D(hpv3NC_id)
  S->hpv1CC_id = cNEW_8D(S->hpv1CC_id)
  S->hpv2CC_id = cNEW_8D(S->hpv2CC_id)
  S->hpv3CC_id = cNEW_8D(S->hpv3CC_id)
  S->hpv1CIN_id = cNEW_8D(S->hpv1CIN_id)
  S->hpv2CIN_id = cNEW_8D(S->hpv2CIN_id)
  S->hpv3CIN_id = cNEW_8D(S->hpv3CIN_id)
  S->CIN_id = cNEW_8D(S->CIN_id)
  S->CC_id = cNEW_8D(S->CC_id)      // #new_age_structre
  S->sexf_id = cNEW_8D(S->sexf_id)
  S->sexm_id = cNEW_8D(S->sexm_id)
  S->risk_id = cNEW_8D(S->risk_id)
  S->age_id = cNEW_8D(S->age_id)
  S->age_id2 = cNEW_8D(S->age_id2)

  //    hivS_id (:, : , : , 1,     1, : , : , : ) = 1;% Susceptible for HIV - F
  //    hiv_id  (:, : , : , 2 : 4, 1, : , : , : ) = 1;% infected with HIV - F
  //    acute_id(:, : , : , 2,     1, : , : , : ) = 1;% infected with HIV - F
  //    hiv2_id (:, : , : , 3,     1, : , : , : ) = 1;% HIV infected, not on ART - F
  //    art_id  (:, : , : , 4,     1, : , : , : ) = 1;% on ART - F
  //    noart_id(:, : , : , 2 : 3, 1, : , : , : ) = 1; % not on ART - F


  LOOP_H1 LOOP_H2 LOOP_H3 LOOP_R LOOP_A LOOP_V1 {
    SET(S->hivS_id,   _dh1, _dh2, _dh3, 1, 1, _dr, _da, _dv1, 1)
    SET(S->hiv_id,    _dh1, _dh2, _dh3, 2, 1, _dr, _da, _dv1, 1)
    SET(S->hiv_id,    _dh1, _dh2, _dh3, 3, 1, _dr, _da, _dv1, 1)
    SET(S->hiv_id,    _dh1, _dh2, _dh3, 4, 1, _dr, _da, _dv1, 1)
    SET(   acute_id,  _dh1, _dh2, _dh3, 2, 1, _dr, _da, _dv1, 1)
    SET(S->hiv2_id,   _dh1, _dh2, _dh3, 3, 1, _dr, _da, _dv1, 1)
    SET(S->art_id,    _dh1, _dh2, _dh3, 4, 1, _dr, _da, _dv1, 1)
    SET(S->noart_id,  _dh1, _dh2, _dh3, 2, 1, _dr, _da, _dv1, 1)
    SET(S->noart_id,  _dh1, _dh2, _dh3, 3, 1, _dr, _da, _dv1, 1)
  }

  //    hivS_id (1 : 3, 1 : 3, 1 : 3, 1,     2, : , : , : ) = 1;% Susceptible for HIV - M
  //     hiv_id (1 : 3, 1 : 3, 1 : 3, 2 : 4, 2, : , : , : ) = 1;% infected with HIV - M
  //    acute_id(1 : 3, 1 : 3, 1 : 3, 2,     2, : , : , : ) = 1;% infected with HIV - M
  //    hiv2_id (1 : 3, 1 : 3, 1 : 3, 3,     2, : , : , : ) = 1;% HIV infected, not on ART - M
  //    art_id  (1 : 3, 1 : 3, 1 : 3, 4,     2, : , : , : ) = 1;% on ART - M
  //    noart_id(1 : 3, 1 : 3, 1 : 3, 2 : 3, 2, : , : , : ) = 1; % not on ART - M

  LOOP_R LOOP_A LOOP_V1
    for (int _dh1=1; _dh1<=3; _dh1++)
      for (int _dh2 = 1; _dh2 <= 3; _dh2++)
        for (int _dh3 = 1; _dh3 <= 3; _dh3++) {
          SET(S->hivS_id,  _dh1, _dh2, _dh3, 1, 2, _dr, _da, _dv1, 1)
          SET(S->hiv_id,   _dh1, _dh2, _dh3, 2, 2, _dr, _da, _dv1, 1)
          SET(S->hiv_id,   _dh1, _dh2, _dh3, 3, 2, _dr, _da, _dv1, 1)
          SET(S->hiv_id,   _dh1, _dh2, _dh3, 4, 2, _dr, _da, _dv1, 1)
          SET(   acute_id, _dh1, _dh2, _dh3, 2, 2, _dr, _da, _dv1, 1)
          SET(S->hiv2_id,  _dh1, _dh2, _dh3, 3, 2, _dr, _da, _dv1, 1)
          SET(S->art_id,   _dh1, _dh2, _dh3, 4, 2, _dr, _da, _dv1, 1)
          SET(S->noart_id, _dh1, _dh2, _dh3, 2, 2, _dr, _da, _dv1, 1)
          SET(S->noart_id, _dh1, _dh2, _dh3, 3, 2, _dr, _da, _dv1, 1)

        }

  //    hpv1S_id(1, :, : , : , 1, : , : , : ) = 1;% susceptible to 16 / 18vtHPV
  //    hpv2S_id(:, 1, : , : , 1, : , : , : ) = 1;% susceptible to oth vtHPV
  //    hpv3S_id(:, : , 1, : , 1, : , : , : ) = 1;% susceptible to nvtHPV
  //    hpv1I_id(2, :, : , : , 1, : , : , : ) = 1;% can develop into CIN2 / 3 16 / 18vtHPV
  //    hpv2I_id(:, 2, : , : , 1, : , : , : ) = 1;% can develop into CIN2 / 3 oth vtHPV
  //    hpv3I_id(:, : , 2, : , 1, : , : , : ) = 1;% can develop into CIN2 / 3 nvtHPV
  //  hpv1CIN_id(4, :, : , : , 1, : , : , : ) = 1;% CIN2 + 16 / 18HPV
  //  hpv2CIN_id(:, 4, : , : , 1, : , : , : ) = 1;% CIN2 + oth vtHPV
  //  hpv3CIN_id(:, : , 4, : , 1, : , : , : ) = 1;% CIN2 + nvtHPV
  //      CIN_id(4, :, : , : , 1, : , : , :) = 1;% CIN2 + 16 / 18HPV
  //      CIN_id(:, 4, : , : , 1, : , : , :) = 1;% CIN2 + oth vtHPV
  //      CIN_id(:, : , 4, : , 1, : , : , :) = 1;% CIN2 + nvtHPV
  //   hpv1NC_id(4, :, : , : , 1, : , : , : ) = 1;% Screening would also detect CC, but only able to treat CIN2 +
  //   hpv2NC_id(:, 4, : , : , 1, : , : , : ) = 1;
  //   hpv3NC_id(:, : , 4, : , 1, : , : , : ) = 1;
  //   hpv1CC_id(5, :, : , : , 1, : , : , : ) = 1;% CC mortalityand treatment, same reagrdless of CIN2 type
  //   hpv2CC_id(:, 5, : , : , 1, : , : , : ) = 1;
  //   hpv3CC_id(:, : , 5, : , 1, : , : , : ) = 1;

  //    #new_age_structre
  //    CC_id(5, :, : , : , 1, : , : , : ) = 1;% CIN2 + 16 / 18HPV
  //    CC_id(:, 5, : , : , 1, : , : , : ) = 1; % CIN2 + oth vtHPV
  //    CC_id(:, : , 5, : , 1, : , : , : ) = 1; % CIN2 + nvtHPV




  LOOP_I LOOP_R LOOP_A LOOP_V1 {
    LOOP_H3 {
      LOOP_H2 {
        SET(S->hpv1S_id,   1, _dh2, _dh3, _di, 1, _dr, _da, _dv1, 1)
        SET(S->hpv1I_id,   2, _dh2, _dh3, _di, 1, _dr, _da, _dv1, 1)
        SET(S->hpv1CIN_id, 4, _dh2, _dh3, _di, 1, _dr, _da, _dv1, 1)
        SET(S->CIN_id,     4, _dh2, _dh3, _di, 1, _dr, _da, _dv1, 1)
        SET(S->CC_id,      5, _dh2, _dh3, _di, 1, _dr, _da, _dv1, 1) // #new_age_structure
        SET(   hpv1NC_id,  4, _dh2, _dh3, _di, 1, _dr, _da, _dv1, 1)
        SET(S->hpv1CC_id,  5, _dh2, _dh3, _di, 1, _dr, _da, _dv1, 1)
      }
      LOOP_H1 {
        SET(S->hpv2S_id,   _dh1, 1, _dh3, _di, 1, _dr, _da, _dv1, 1)
        SET(S->hpv2I_id,   _dh1, 2, _dh3, _di, 1, _dr, _da, _dv1, 1)
        SET(S->hpv2CIN_id, _dh1, 4, _dh3, _di, 1, _dr, _da, _dv1, 1)
        SET(S->CIN_id,     _dh1, 4, _dh3, _di, 1, _dr, _da, _dv1, 1)
        SET(S-> CC_id,     _dh1, 5, _dh3, _di, 1, _dr, _da, _dv1, 1) // #new_age_structure
        SET(   hpv2NC_id,  _dh1, 4, _dh3, _di, 1, _dr, _da, _dv1, 1)
        SET(S->hpv2CC_id,  _dh1, 5, _dh3, _di, 1, _dr, _da, _dv1, 1)
      }
    }

    LOOP_H1 LOOP_H2 {
      SET(S->hpv3S_id,   _dh1, _dh2, 1, _di, 1, _dr, _da, _dv1, 1)
      SET(S->hpv3I_id,   _dh1, _dh2, 2, _di, 1, _dr, _da, _dv1, 1)
      SET(S->hpv3CIN_id, _dh1, _dh2, 4, _di, 1, _dr, _da, _dv1, 1)
      SET(S->CIN_id,     _dh1, _dh2, 4, _di, 1, _dr, _da, _dv1, 1)
      SET(S->CC_id,      _dh1, _dh2, 5, _di, 1, _dr, _da, _dv1, 1) // #new_age_structure
      SET(   hpv3NC_id,  _dh1, _dh2, 4, _di, 1, _dr, _da, _dv1, 1)
      SET(S->hpv3CC_id,  _dh1, _dh2, 5, _di, 1, _dr, _da, _dv1, 1)
    }

    //    hpv1S_id(1,     1 : 3, 1 : 3, : , 2, : , : , : ) = 1;% susceptible to 16 / 18vtHPV
    //    hpv2S_id(1 : 3, 1,     1 : 3, : , 2, : , : , : ) = 1;% susceptible to oth vtHPV
    //    hpv3S_id(1 : 3, 1 : 3, 1,     : , 2, : , : , : ) = 1;% susceptible to nvtHPV

    for (int _dh3 = 1; _dh3 <= 3; _dh3++) {
      for (int _dh2 = 1; _dh2 <= 3; _dh2++) SET(S->hpv1S_id, 1, _dh2, _dh3, _di, 2, _dr, _da, _dv1, 1)
      for (int _dh1 = 1; _dh1 <= 3; _dh1++) SET(S->hpv2S_id, _dh1, 1, _dh3, _di, 2, _dr, _da, _dv1, 1)
    }
    for (int _dh1 = 1; _dh1 <= 3; _dh1++)
      for (int _dh2 = 1; _dh2 <= 3; _dh2++) SET(S->hpv3S_id, _dh1, _dh2, 1, _di, 2, _dr, _da, _dv1, 1)
  }

  //    hpv1_id(2, :, : , : , : , : , : , : ) = 1;% infected  with 16 / 18HPV no CIN2 + or CC
  //    hpv2_id(:, 2, : , : , : , : , : , : ) = 1;% infected  with oth vtHPV  no CIN2 + or CC
  //    hpv3_id(:, : , 2, : , : , : , : , : ) = 1;% infected  with nvtHPV  no CIN2 + or CC
  //    hpv1R_id(3, :, : , : , : , : , : , : ) = 1;% Temporal immunity with 16 / 18HPV
  //    hpv2R_id(:, 3, : , : , : , : , : , : ) = 1;% Temporal immunity oth oth vtHPV
  //    hpv3R_id(:, : , 3, : , : , : , : , : ) = 1;% Temporal immunity oth nvtHPV


  LOOP_I LOOP_S LOOP_R LOOP_A LOOP_V1 {
    LOOP_H3 {
      LOOP_H2 {
        SET(S->hpv1_id,  2, _dh2, _dh3, _di, _ds, _dr, _da, _dv1, 1)
        SET(S->hpv1R_id, 3, _dh2, _dh3, _di, _ds, _dr, _da, _dv1, 1.0)
      }
      LOOP_H1 {
        SET(S->hpv2_id,  _dh1, 2, _dh3, _di, _ds, _dr, _da, _dv1, 1)
        SET(S->hpv2R_id, _dh1, 3, _dh3, _di, _ds, _dr, _da, _dv1, 1)
      }
    }
    LOOP_H1 LOOP_H2 {
      SET(S->hpv3_id,  _dh1, _dh2, 2, _di, _ds, _dr, _da, _dv1, 1)
      SET(S->hpv3R_id, _dh1, _dh2, 3, _di, _ds, _dr, _da, _dv1, 1)
    }
  }

   //    hpv12_id(([2, 4, 5]), :, : , : , 1, : , : , : ) = 1;% infected with vtHPV
   //    hpv12_id(:, ([2, 4, 5]), : , : , 1, : , : , : ) = 1;
   //    hpv12_id(2, :, : , : , 2, : , : , : ) = 1;% infected with vtHPV
   //    hpv12_id(:, 2, : , : , 2, : , : , : ) = 1;

  LOOP_H3 LOOP_I LOOP_R LOOP_A LOOP_V1 {
    LOOP_H2 {
      SET(S->hpv12_id, 2, _dh2, _dh3, _di, 1, _dr, _da, _dv1, 1)
      SET(S->hpv12_id, 4, _dh2, _dh3, _di, 1, _dr, _da, _dv1, 1)
      SET(S->hpv12_id, 5, _dh2, _dh3, _di, 1, _dr, _da, _dv1, 1)
      SET(S->hpv12_id, 2, _dh2, _dh3, _di, 2, _dr, _da, _dv1, 1)
    }
    LOOP_H1 {
      SET(S->hpv12_id, _dh1, 2, _dh3, _di, 1, _dr, _da, _dv1, 1)
      SET(S->hpv12_id, _dh1, 4, _dh3, _di, 1, _dr, _da, _dv1, 1)
      SET(S->hpv12_id, _dh1, 5, _dh3, _di, 1, _dr, _da, _dv1, 1)
      SET(S->hpv12_id, _dh1, 2, _dh3, _di, 2, _dr, _da, _dv1, 1)
    }
  }


  //    % here careful, need to have any infected, so on top of row + column

  //    hpv123_id = hpv12_id;

  S->hpv123_id = cNEW_COPY_8D(S->hpv123_id, S->hpv12_id)

  //    hpv123_id(:, : , ([2, 4, 5]), : , 1, : , : , : ) = 1;% infected with vtHPV or nvtHPV
  //    hpv123_id(:, : , 2,           : , 2, : , : , : ) = 1;

  LOOP_H1 LOOP_H2 LOOP_I LOOP_R LOOP_A LOOP_V1 {
    SET(S->hpv123_id, _dh1, _dh2, 2, _di, 1, _dr, _da, _dv1, 1)
    SET(S->hpv123_id, _dh1, _dh2, 4, _di, 1, _dr, _da, _dv1, 1)
    SET(S->hpv123_id, _dh1, _dh2, 5, _di, 1, _dr, _da, _dv1, 1)
    SET(S->hpv123_id, _dh1, _dh2, 2, _di, 2, _dr, _da, _dv1, 1)
  }

  //    sexf_id(:, : , : , : , 1, : , : , : ) = 1;% women
  //    sexm_id(:, : , : , : , 2, : , : , : ) = 1;% men
  //    vacc_id(:, : , : , : , : , : , : , 2) = 1;% vaccinated
  //    risk_id(:, : , : , : , : , 1, : , : ) = 1;% LR
  //    risk_id(:, : , : , : , : , 2, : , : ) = 2;% MR
  //    risk_id(:, : , : , : , : , 3, : , : ) = 3;% HR

  LOOP_H1 LOOP_H2 LOOP_H3 LOOP_I LOOP_A LOOP_V1 {
    LOOP_R {
      SET(S->sexf_id, _dh1, _dh2, _dh3, _di, 1, _dr, _da, _dv1, 1)
      SET(S->sexm_id, _dh1, _dh2, _dh3, _di, 2, _dr, _da, _dv1, 1)
      LOOP_S SET(S->vacc_id, _dh1, _dh2, _dh3, _di, _ds, _dr, _da, 2, 1)
    }
    LOOP_S {
      SET(S->risk_id, _dh1, _dh2, _dh3, _di, _ds, 1, _da, _dv1, 1)
      SET(S->risk_id, _dh1, _dh2, _dh3, _di, _ds, 2, _da, _dv1, 2)
      SET(S->risk_id, _dh1, _dh2, _dh3, _di, _ds, 3, _da, _dv1, 3)
    }
  }
    // old_age_structure:
    //  age_id(:, : , : , : , : , : , 1, : ) = 1;
    //  age_id(:, : , : , : , : , : , 2, : ) = 2;
    //  age_id(:, : , : , : , : , : , 3, : ) = 3;
    //  age_id(:, : , : , : , : , : , 4, : ) = 4;
    //  age_id(:, : , : , : , : , : , 5, : ) = 5;

    // LOOP_H1 LOOP_H2 LOOP_H3 LOOP_I LOOP_S LOOP_R LOOP_V1{
    // for (int _da = 1; _da <= 5; _da++) SET(S->age_id, _dh1, _dh2, _dh3, _di, _ds, _dr, _da, _dv1, (double)_da);
    //}


    // new_age_structure:
    // age_id(:, :, :, :, :, :, 1,      : ) = 1;
    // age_id(:, :, :, :, :, :, 2 : 3,  : ) = 2;
    // age_id(:, :, :, :, :, :, 4 : 5,  : ) = 3;
    // age_id(:, :, :, :, :, :, 6 : 8,  : ) = 4;
    // age_id(:, :, :, :, :, :, 9 : 13, : ) = 5;

    LOOP_H1 LOOP_H2 LOOP_H3 LOOP_I LOOP_S LOOP_R LOOP_V1 {
      for (int _da = 1; _da <= 1;  _da++) GET_ALL(S->age_id) = 1;
      for (int _da = 2; _da <= 3;  _da++) GET_ALL(S->age_id) = 2;
      for (int _da = 4; _da <= 5;  _da++) GET_ALL(S->age_id) = 3;
      for (int _da = 6; _da <= 8;  _da++) GET_ALL(S->age_id) = 4;
      for (int _da = 9; _da <= 13; _da++) GET_ALL(S->age_id) = 5;
    }

    // this is not the best formulation but I need to understand the
    // SET better to make it work
    LOOP_H1 LOOP_H2 LOOP_H3 LOOP_I LOOP_S LOOP_R LOOP_V1 {
      for (int _da = 1; _da <= 1;  _da++) GET_ALL(S->age_id2) = 1;
      for (int _da = 2; _da <= 2;  _da++) GET_ALL(S->age_id2) = 2;
      for (int _da = 3; _da <= 3;  _da++) GET_ALL(S->age_id2) = 3;
      for (int _da = 4; _da <= 4;  _da++) GET_ALL(S->age_id2) = 4;
      for (int _da = 5; _da <= 5;  _da++) GET_ALL(S->age_id2) = 5;
      for (int _da = 6; _da <= 6;  _da++) GET_ALL(S->age_id2) = 6;
      for (int _da = 7; _da <= 7;  _da++) GET_ALL(S->age_id2) = 7;
      for (int _da = 8; _da <= 8;  _da++) GET_ALL(S->age_id2) = 8;
      for (int _da = 9; _da <= 9;  _da++) GET_ALL(S->age_id2) = 9;
      for (int _da = 10;_da <= 10; _da++) GET_ALL(S->age_id2) = 10;
      for (int _da = 11;_da <= 11; _da++) GET_ALL(S->age_id2) = 11;
      for (int _da = 12;_da <= 12; _da++) GET_ALL(S->age_id2) = 12;
      for (int _da = 13;_da <= 13; _da++) GET_ALL(S->age_id2) = 13;
    }

   // This did not work, I will investigate further
   // LOOP_H1 LOOP_H2 LOOP_H3 LOOP_I LOOP_S LOOP_R LOOP_V1 {
   // SET(S->age_id2, _dh1, _dh2, _dh3, _di, _ds, _dr, 1, _dv1, 1)
   // SET(S->age_id2, _dh1, _dh2, _dh3, _di, _ds, _dr, 2, _dv1, 2)
   // SET(S->age_id2, _dh1, _dh2, _dh3, _di, _ds, _dr, 3, _dv1, 3)
   // SET(S->age_id2, _dh1, _dh2, _dh3, _di, _ds, _dr, 4, _dv1, 4)
   // SET(S->age_id2, _dh1, _dh2, _dh3, _di, _ds, _dr, 5, _dv1, 5)
   // SET(S->age_id2, _dh1, _dh2, _dh3, _di, _ds, _dr, 6, _dv1, 6)
   // SET(S->age_id2, _dh1, _dh2, _dh3, _di, _ds, _dr, 7, _dv1, 7)
   // SET(S->age_id2, _dh1, _dh2, _dh3, _di, _ds, _dr, 8, _dv1, 8)
   // SET(S->age_id2, _dh1, _dh2, _dh3, _di, _ds, _dr, 9, _dv1, 9)
   // SET(S->age_id2, _dh1, _dh2, _dh3, _di, _ds, _dr, 10, dv1, 10)
   // SET(S->age_id2, _dh1, _dh2, _dh3, _di, _ds, _dr, 11, dv1, 11)
   // SET(S->age_id2, _dh1, _dh2, _dh3, _di, _ds, _dr, 12, dv1, 12)
   // SET(S->age_id2, _dh1, _dh2, _dh3, _di, _ds, _dr, 13, dv1, 13)
   // }


  //#define TEST_OUTPUT 1

  #ifdef TEST_OUTPUT
  FILE* fout = fopen("debug.bin", "wb");

  // Output for debugging 1

  LOOP_ALL fwrite(&GET_ALL(S->hivS_id), 1, 1, fout);
  LOOP_ALL fwrite(&GET_ALL(S->hiv_id), 1, 1, fout);
  LOOP_ALL fwrite(&GET_ALL(   acute_id), 1, 1, fout);
  LOOP_ALL fwrite(&GET_ALL(S->hiv2_id), 1, 1, fout);
  LOOP_ALL fwrite(&GET_ALL(S->art_id), 1, 1, fout);
  LOOP_ALL fwrite(&GET_ALL(S->noart_id), 1, 1, fout);
  LOOP_ALL fwrite(&GET_ALL(S->hpv1S_id), 1, 1, fout);
  LOOP_ALL fwrite(&GET_ALL(S->hpv2S_id), 1, 1, fout);
  LOOP_ALL fwrite(&GET_ALL(S->hpv3S_id), 1, 1, fout);
  LOOP_ALL fwrite(&GET_ALL(S->hpv1I_id), 1, 1, fout);
  LOOP_ALL fwrite(&GET_ALL(S->hpv2I_id), 1, 1, fout);
  LOOP_ALL fwrite(&GET_ALL(S->hpv3I_id), 1, 1, fout);
  LOOP_ALL fwrite(&GET_ALL(S->hpv1_id), 1, 1, fout);
  LOOP_ALL fwrite(&GET_ALL(S->hpv2_id), 1, 1, fout);
  LOOP_ALL fwrite(&GET_ALL(S->hpv3_id), 1, 1, fout);
  LOOP_ALL fwrite(&GET_ALL(S->hpv1R_id), 1, 1, fout);
  LOOP_ALL fwrite(&GET_ALL(S->hpv2R_id), 1, 1, fout);
  LOOP_ALL fwrite(&GET_ALL(S->hpv3R_id), 1, 1, fout);
  LOOP_ALL fwrite(&GET_ALL(S->hpv1CIN_id), 1, 1, fout);
  LOOP_ALL fwrite(&GET_ALL(S->hpv2CIN_id), 1, 1, fout);
  LOOP_ALL fwrite(&GET_ALL(S->hpv3CIN_id), 1, 1, fout);
  LOOP_ALL fwrite(&GET_ALL(S->CIN_id), 1, 1, fout);
  LOOP_ALL fwrite(&GET_ALL(S->hpv1CC_id), 1, 1, fout);
  LOOP_ALL fwrite(&GET_ALL(S->hpv2CC_id), 1, 1, fout);
  LOOP_ALL fwrite(&GET_ALL(S->hpv3CC_id), 1, 1, fout);
  LOOP_ALL fwrite(&GET_ALL(   hpv1NC_id), 1, 1, fout);
  LOOP_ALL fwrite(&GET_ALL(   hpv2NC_id), 1, 1, fout);
  LOOP_ALL fwrite(&GET_ALL(   hpv3NC_id), 1, 1, fout);
  LOOP_ALL fwrite(&GET_ALL(S->hpv12_id), 1, 1, fout);
  LOOP_ALL fwrite(&GET_ALL(S->hpv123_id), 1, 1, fout);
  LOOP_ALL fwrite(&GET_ALL(S->vacc_id), 1, 1, fout);
  LOOP_ALL fwrite(&GET_ALL(S->sexf_id), 1, 1, fout);
  LOOP_ALL fwrite(&GET_ALL(S->sexm_id), 1, 1, fout);
  LOOP_ALL fwrite(&GET_ALL(S->risk_id), 1, 1, fout);
  LOOP_ALL fwrite(&GET_ALL(S->age_id), 1, 1, fout);
  #endif


  /*
  % %Parameters
    %%%%%%%%%%%%%%%
    [betahpv1, betahpv2, betahpv3, betahivacute, betahiv2, betahivart, betahivacute_ai, betahiv2_ai, betahivart_ai, ...
    delta1, delta2, delta3, vxxwane, nacts, nacts_comm, p_anal, pcomm_anal, ...
    mu, nu, omega1, omega2, omega3, eta, tau, omikron, theta1, theta2, theta3, alpha, ...
    sigma1, sigma2, sigma3, psi1, psi2, psi3, pi1, pi2, pi3, ups1, ups2, ups3, hups1, hups2, hups3] = deal(zeros(h1, h2, h3, i, s, r, a, v1));
  */

  S->betahpv1 = dNEW_8D(S->betahpv1)
  S->betahpv2 = dNEW_8D( S->betahpv2 )
  S->betahpv3 = dNEW_8D( S->betahpv3 )
  S->betahivacute = dNEW_8D(S->betahivacute);
  S->betahiv2 = dNEW_8D(S->betahiv2)
  S->betahivart = dNEW_8D(S->betahivart)
  S->betahivacute_ai = dNEW_8D(S->betahivacute_ai)
  S->betahiv2_ai = dNEW_8D(S->betahiv2_ai)
  S->betahivart_ai = dNEW_8D(S->betahivart_ai)
  S->delta1 = dNEW_8D(S->delta1)
  S->delta2 = dNEW_8D(S->delta2)
  S->delta3 = dNEW_8D(S->delta3)
  S->vxxwane = dNEW_8D(S->vxxwane)
  S->nacts = dNEW_8D(S->nacts)
  S->nacts_comm = dNEW_8D(S->nacts_comm)
  S->p_anal = dNEW_8D(S->p_anal)
  S->pcomm_anal = dNEW_8D(S->pcomm_anal)
  S->mu = dNEW_8D(S->mu)
  S->nu = dNEW_8D(S->nu)
  S->omega1 = dNEW_8D(S->omega1)
  S->omega2 = dNEW_8D(S->omega2)
  S->omega3 = dNEW_8D(S->omega3)
  S->eta = dNEW_8D(S->eta)
  S->tau = dNEW_8D(S->tau)
  S->omikron = dNEW_8D(S->omikron)
  S->theta1 = dNEW_8D(S->theta1)
  S->theta2 = dNEW_8D(S->theta2)
  S->theta3 = dNEW_8D(S->theta3)
  S->alpha = dNEW_8D(S->alpha)
  S->sigma1 = dNEW_8D(S->sigma1)
  S->sigma2 = dNEW_8D(S->sigma2)
  S->sigma3 = dNEW_8D(S->sigma3)
  S->psi1 = dNEW_8D(S->psi1)
  S->psi2 = dNEW_8D(S->psi2)
  S->psi3 = dNEW_8D(S->psi3)
  S->pi1 = dNEW_8D(S->pi1)
  S->pi2 = dNEW_8D(S->pi2)
  S->pi3 = dNEW_8D(S->pi3)
  S->ups1 = dNEW_8D(S->ups1)
  S->ups2 = dNEW_8D(S->ups2)
  S->ups3 = dNEW_8D(S->ups3)
  S->hups1 = dNEW_8D(S->hups1)
  S->hups2 = dNEW_8D(S->hups2)
  S->hups3 = dNEW_8D(S->hups3)

  /*
    %%%%%%%%%%%%%%%

    %% TRANSMISSION PROBABILITIES

    % HIV related params -- need to connect beta with the infected person as
    % infectiousness is different by HIV state BUT susceptibility is also
    % different by HPV status -- more stratification within the diffeq equation

    % mcirc added on top of this in the niterations
  */

  //    % Acute HIV - VI

  //    betahivacute(hivS_id == 1 & sexf_id == 1) = hp0(2) * hp0(1) * hp0(5);% from acute men
  //    betahivacute(hivS_id == 1 & sexf_id == 1 & hpv123_id == 1) = hp0(2) * hp0(1) * hp0(5) * hp0(7);

  //    % Acute HIV - AI // HPV infection assumed to not play a role in increased
  //    % transmission

  //    betahivacute_ai(hivS_id == 1 & sexf_id == 1) = hp0(2) * hp0(1) * hp0(4) * hp0(5);% from acute men, woman RAI

  //    % HIV after acute, not on ART

  //    betahiv2(hivS_id == 1 & sexf_id == 1) = hp0(2) * hp0(1);%
  //    betahiv2(hivS_id == 1 & sexf_id == 1 & hpv123_id == 1) = hp0(2) * hp0(1) * hp0(7);%

  //    betahiv2_ai(hivS_id == 1 & sexf_id == 1) = hp0(2) * hp0(1);% woman RAI, man not on ART

  //    % HIV on ART

  //    betahivart(hivS_id == 1 & sexf_id == 1) = hp0(2) * hp0(1) * (1 - hp0(6));% infectivity from men on ART
  //    betahivart(hivS_id == 1 & sexf_id == 1 & hpv123_id == 1) = hp0(2) * hp0(1) * (1 - hp0(6)) * hp0(7);% i

  //    betahivart_ai(hivS_id == 1 & sexf_id == 1) = hp0(1) * hp0(4) * (1 - hp0(6));% infectivity from men on ART, RAI

  //    %**********%

  // These are factored out of the loop to make the setting faster...

  double _p1_p2 = p->pp0[1] * p->pp0[2];
  double _p1_p2_p4 = _p1_p2 * p->pp0[4];
  double _p1_p3 = p->pp0[1] * p->pp0[3];
  double _p1_p3_p4 = _p1_p3 * p->pp0[4];
  double _p1_p4 = p->pp0[1] * p->pp0[4];

  double _h1_h2 = p->hp0[1] * p->hp0[2];
  double _h1_h2_h4 = _h1_h2 * p->hp0[4];
  double _h1_h2_h4_h5 = _h1_h2_h4 * p->hp0[5];
  double _h1_h2_h5 = _h1_h2 * p->hp0[5];
  double _h1_h2_h5_h7 = _h1_h2_h5 * p->hp0[7];
  double _h1_h2_m6 = _h1_h2 * (1.0 - p->hp0[6]);
  double _h1_h2_m6_h7 = _h1_h2 * (1.0 - p->hp0[6]) * p->hp0[7];
  double _h1_h2_h7 = _h1_h2 * p->hp0[7];
  double _h1_h4 = p->hp0[1] * p->hp0[4];
  double _h1_h4_m6 = _h1_h4 * (1 - p->hp0[6]);

  LOOP_ALL {
    if IS(S->hivS_id, 1) {
      if IS(S->sexf_id, 1) {
        GET_ALL(S->betahivacute) = _h1_h2_h5;
        GET_ALL(S->betahivacute_ai) = _h1_h2_h4_h5;
        GET_ALL(S->betahiv2) = _h1_h2;
        GET_ALL(S->betahiv2_ai) = _h1_h2;
        GET_ALL(S->betahivart) = _h1_h2_m6;
        GET_ALL(S->betahivart_ai) = _h1_h4_m6;

        if IS(S->hpv123_id, 1) {
          GET_ALL(S->betahivacute) = _h1_h2_h5_h7;
          GET_ALL(S->betahiv2) = _h1_h2_h7;
          GET_ALL(S->betahivart) = _h1_h2_m6_h7;
        }
      }
    }

    bool _hiv_id = IS(S->hiv_id, 1);

    //    % HPV beta
    //    betahpv1(hpv1S_id == 1) = pp0(1);% vtHPV acquisition_act prob
    //    betahpv1(hpv1S_id == 1 & hiv_id == 1) = pp0(1) * pp0(4);% vvtHPV acquisition_act prob

    if IS(S->hpv1S_id, 1) {
      GET_ALL(S->betahpv1) = p->pp0[1];
      if (_hiv_id) GET_ALL(S->betahpv1) = _p1_p4;
    }

    //    betahpv2(hpv2S_id == 1) = pp0(1) * pp0(2);% vtHPV acquisition_act prob
    //    betahpv2(hpv2S_id == 1 & hiv_id == 1) = pp0(1) * pp0(2) * pp0(4);% vvtHPV acquisition_act prob

    if IS(S->hpv2S_id, 1) {
      GET_ALL(S->betahpv2) = _p1_p2;
      if (_hiv_id) GET_ALL(S->betahpv2) = _p1_p2_p4;
    }

    //    betahpv3(hpv3S_id == 1) = pp0(1) * pp0(3);% vtHPV acquisition_act prob
    //    betahpv3(hpv3S_id == 1 & hiv_id == 1) = pp0(1) * pp0(3) * pp0(4);% vvtHPV acquisition_act prob

    if IS(S->hpv3S_id, 1) {
      GET_ALL(S->betahpv3) = _p1_p3;
      if (_hiv_id) GET_ALL(S->betahpv3) = _p1_p3_p4;
    }

  }

  //    coneffhiv = hp0(9);% efficacy of condom use against hiv

  // Used once as a constant - see #define in diffeq_vacc.cpp

  //    %% HPV parameters
  //    % clearance of HPV
  //    % faster clearance for hpv2and hpv3 than for hpv1

  //    % proportion who develop immunity
  //    em = pp0(7);
  //  % proportion who clear infection once treated
  //    zeta = pp0(8);
  //  % proportion who clear infection if natural regression from CIN2 / CIN3
  //    qu = pp0(9);

  #define qu p->pp0[9]

  double _p5_p10 = p->pp0[5]  * p->pp0[10];
  double _p5_p10_p11 = p->pp0[10] * p->pp0[11] * p->pp0[5];
  double _p5_p10_p12 = p->pp0[10] * p->pp0[12] * p->pp0[5];
  double _p5_p10_p12_p32 = _p5_p10_p12 * p->pp0[32];
  double _p5_p10_p11_p32 = p->pp0[10] * p->pp0[11] * p->pp0[5] * p->pp0[32];
  double _p5_p10_p32 = _p5_p10 * p->pp0[32];
  double _p10_p11 = p->pp0[10] * p->pp0[11];
  double _p10_p12 = p->pp0[10] * p->pp0[12];
  double _p10_p12_p32 = _p10_p12 * p->pp0[32];

  double _p10_p11_p32 = _p10_p11 * p->pp0[32];
  double _p10_p32 = p->pp0[10] * p->pp0[32];

  // Defines for em, zeta, qu, moved to where they are needed...

  //  % Clearance of HPV -- HIV affects, ART may also affect ADD A PARAMETER HERE
  //    sigma1(hpv1_id == 1 & sexf_id == 1) = pp0(10);
  //    sigma1(hpv1_id == 1 & hiv_id == 1 & sexf_id == 1) = pp0(10) * pp0(5);
  //    sigma1(hpv1_id == 1 & sexm_id == 1) = pp0(10) * pp0(32);
  //    sigma1(hpv1_id == 1 & hiv_id == 1 & sexm_id == 1) = pp0(10) * pp0(5) * pp0(32);

  //  sigma2(hpv2_id == 1 & sexf_id == 1) = pp0(10) * pp0(11);
  //  sigma2(hpv2_id == 1 & hiv_id == 1 & sexf_id == 1) = pp0(10) * pp0(11) * pp0(5);
  //  sigma2(hpv2_id == 1 & sexm_id == 1) = pp0(10) * pp0(11) * pp0(32);
  //  sigma2(hpv2_id == 1 & hiv_id == 1 & sexm_id == 1) = pp0(10) * pp0(11) * pp0(5) * pp0(32);

  //  sigma3(hpv3_id == 1 & sexf_id == 1) = pp0(10) * pp0(12);
  //  sigma3(hpv3_id == 1 & hiv_id == 1 & sexf_id == 1) = pp0(10) * pp0(12) * pp0(5);

  //  sigma3(hpv3_id == 1 & sexm_id == 1) = pp0(10) * pp0(12) * pp0(32);
  //  sigma3(hpv3_id == 1 & hiv_id == 1 & sexm_id == 1) = pp0(10) * pp0(12) * pp0(5) * pp0(32);


  LOOP_ALL {
    if IS(S->sexf_id, 1) {
      if IS(S->hpv1_id, 1) {
        GET_ALL(S->sigma1) = p->pp0[10];
        if IS(S->hiv_id, 1) GET_ALL(S->sigma1) = _p5_p10;
      }
      if IS(S->hpv2_id, 1) {
        GET_ALL(S->sigma2) = _p10_p11;
        if IS(S->hiv_id, 1) GET_ALL(S->sigma2) = _p5_p10_p11;
      }
      if IS(S->hpv3_id, 1) {
        GET_ALL(S->sigma3) = _p10_p12;
        if IS(S->hiv_id, 1) GET_ALL(S->sigma3) = _p5_p10_p12;
       }
    }
    if IS(S->sexm_id, 1) {
      if IS(S->hpv1_id, 1) {
        GET_ALL(S->sigma1) = _p10_p32;
        if IS(S->hiv_id, 1) GET_ALL(S->sigma1) = _p5_p10_p32;
      }
      if IS(S->hpv2_id, 1) {
        GET_ALL(S->sigma2) = _p10_p11_p32;
        if IS(S->hiv_id, 1) GET_ALL(S->sigma2) = _p5_p10_p11_p32;
      }
      if IS(S->hpv3_id, 1) {
        GET_ALL(S->sigma3) = _p10_p12_p32;
        if IS(S->hiv_id, 1) GET_ALL(S->sigma3) = _p5_p10_p12_p32;
      }
    }

    //  % Waning of natural immunity,
    //    % same across age varies by HPV typesand by HIV sttus
    //    delta1(hpv1R_id == 1) = pp0(13); %
    //    delta1(hpv1R_id == 1 & hiv_id == 1) = pp0(13) * pp0(16); %

    double _p13_p16 = p->pp0[13] * p->pp0[16];
    double _p13_p14 = p->pp0[13] * p->pp0[14];
    double _p13_p14_p16 = _p13_p14 * p->pp0[16];
    double _p13_p15 = p->pp0[13] * p->pp0[15];
    double _p13_p15_p16 = _p13_p15 * p->pp0[16];
    double _p17_p18 = p->pp0[17] * p->pp0[18];
    double _p17_p18_p48 = _p17_p18 * p->pp0[48];
    double _p17_p19 = p->pp0[17] * p->pp0[19];
    double _p17_p19_p48 = _p17_p19 * p->pp0[48];
    double _p17_p48 = p->pp0[17] * p->pp0[48];


    if IS(S->hpv1R_id, 1) {
      GET_ALL(S->delta1) = p->pp0[13];
      if IS(S->hiv_id, 1) GET_ALL(S->delta1) = _p13_p16;
    }

    //    delta2(hpv2R_id == 1) = pp0(13) * pp0(14); %
    //    delta2(hpv2R_id == 1 & hiv_id == 1) = pp0(13) * pp0(14) * pp0(16); %

    if IS(S->hpv2R_id, 1) {
      GET_ALL(S->delta2) = _p13_p14;
      if IS(S->hiv_id, 1) GET_ALL(S->delta2) = _p13_p14_p16;
    }

    //    delta3(hpv3R_id == 1) = pp0(13) * pp0(15); %
    //    delta3(hpv3R_id == 1 & hiv_id == 1) = pp0(13) * pp0(15) * pp0(16); %

    if IS(S->hpv3R_id, 1) {
      GET_ALL(S->delta3) = _p13_p15;
      if IS(S->hiv_id, 1) GET_ALL(S->delta3) = _p13_p15_p16;
    }

    //    % Add impact of HIV -- ART may affect ADD A PARAMETER
    //    psi1(hpv1I_id == 1) = pp0(17);% progression to CIN2 / CIN3
    //    psi1(hpv1I_id == 1 & hiv_id == 1) = pp0(17) * pp0(48);% progression to CIN2 / CIN3

    if IS(S->hpv1I_id, 1) {
      GET_ALL(S->psi1) = p->pp0[17];
      if IS(S->hiv_id, 1) GET_ALL(S->psi1) = _p17_p48;
    }

    //    psi2(hpv2I_id == 1) = pp0(17) * pp0(18);% progression to CIN2 / CIN3
    //    psi2(hpv2I_id == 1 & hiv_id == 1) = pp0(17) * pp0(18) * pp0(48);% progression to CIN2 / CIN3

    if IS(S->hpv2I_id, 1) {
      GET_ALL(S->psi2) = _p17_p18;
      if IS(S->hiv_id, 1) GET_ALL(S->psi2) = _p17_p18_p48;
    }

    //    psi3(hpv3I_id == 1) = pp0(17) * pp0(19);% progression to CIN2 / CIN3
    //    psi3(hpv3I_id == 1 & hiv_id == 1) = pp0(17) * pp0(19) * pp0(48);% progression to CIN2 / CIN3

    if IS(S->hpv3I_id, 1) {
      GET_ALL(S->psi3) = _p17_p19;
      if IS(S->hiv_id, 1) GET_ALL(S->psi3) = _p17_p19_p48;
    }

    //    % **********%

    double _p20_p23 = p->pp0[20] * p->pp0[23];
    double _p20_o_p24 = p->pp0[20] / p->pp0[24];
    double _p20_p23_o_p24 = _p20_p23 / p->pp0[24];
    double _p20_p21 = p->pp0[20] * p->pp0[21];
    double _p20_p21_p23 = _p20_p21 * p->pp0[23];
    double _p20_p21_o_p24 = p->pp0[20] * p->pp0[21] / p->pp0[24];
    double _p20_p21_p23_o_p24 = p->pp0[20] * p->pp0[21] * p->pp0[23] / p->pp0[24];
    double _p20_p22 = p->pp0[20] * p->pp0[22];
    double _p20_p22_p23 = _p20_p22 * p->pp0[23];
    double _p20_p22_p23_o_p24 = _p20_p22_p23 / p->pp0[24];
    double _p20_p22_o_p24 = _p20_p22 / p->pp0[24];
    double _p26 = p->pp0[26];
    double _p27 =  p->pp0[27]; // 60+
    double _p26_p28 = p->pp0[26] * p->pp0[28]; // 50+
    double _p25 =  p->pp0[25]; // 70+
    double _p29_p30 = p->pp0[29] * p->pp0[30];
    double _p29_p31 = p->pp0[29] * p->pp0[31];
    double _h18_h19 = p->hp0[18] * p->hp0[19];

    //    omega1(hpv1CIN_id == 1 & age_id < 4) = pp0(20);% regression from CIN2 / CIN3
    //    omega1(hpv1CIN_id == 1 & hiv_id == 1 & age_id < 4) = pp0(20) * pp0(23);% regression from CIN2 / CIN3
    //    omega1(hpv1CIN_id == 1 & age_id > 3) = pp0(20) / pp0(24);% regression from CIN2 / CIN3
    //    omega1(hpv1CIN_id == 1 & hiv_id == 1 & age_id > 3) = pp0(20) * pp0(23) / pp0(24);% regression from CIN2 / CIN3
    //    pi1(hpv1CIN_id == 1 & age_id < 4) = pp0(26);% development of CC in younger people
    //    pi1(hpv1CIN_id == 1 & age_id == 4) = pp0(26) * pp0(27);% development of CC
    //    pi1(hpv1CIN_id == 1 & age_id == 5) = pp0(26) * pp0(28);% development of CC


    if IS(S->hpv1CIN_id, 1) {
      if LT(S->age_id, 3) {
         GET_ALL(S->pi1) = p->pp0[26]/10;
      }
      if LT(S->age_id, 4) {
        GET_ALL(S->omega1) = p->pp0[20];
        if IS(S->hiv_id, 1) GET_ALL(S->omega1) = _p20_p23;
        if IS(S->age_id, 3) GET_ALL(S->pi1) = p->pp0[26];
       // GET_ALL(S->pi1) = p->pp0[26];
      }
      if GT(S->age_id, 3) {
        GET_ALL(S->omega1) = _p20_o_p24;
        if IS(S->hiv_id, 1) GET_ALL(S->omega1) = _p20_p23_o_p24;
        if IS(S->age_id, 4) GET_ALL(S->pi1) = _p26;  // times 2 (11-11-2020) - removed (22-04-2021)
       // else if IS(S->age_id, 5) GET_ALL(S->pi1) = _p26_p28;
      }
      if GT(S->age_id, 4) {
        GET_ALL(S->omega1) = _p20_o_p24/20;
        if IS(S->hiv_id, 1) GET_ALL(S->omega1) = _p20_p23_o_p24/20;
      }
    if GT(S->age_id2, 11) {
            GET_ALL(S->omega1) = _p20_o_p24/50;
            if IS(S->hiv_id, 1) GET_ALL(S->omega1) = _p20_p23_o_p24/50;
          }
   }

    //    omega2(hpv2CIN_id == 1 & age_id < 4) = pp0(20) * pp0(21);% regression from CIN2 / CIN3
    //    omega2(hpv2CIN_id == 1 & hiv_id == 1 & age_id < 4) = pp0(20) * pp0(21) * pp0(23);% regression from CIN2 / CIN3

    //    omega2(hpv2CIN_id == 1 & age_id > 3) = pp0(20) * pp0(21) / pp0(24); % regression from CIN2 / CIN3
    //    omega2(hpv2CIN_id == 1 & hiv_id == 1 & age_id > 3) = pp0(20) * pp0(21) * pp0(23) / pp0(24); % regression from CIN2 / CIN3
    //    pi2(hpv2CIN_id == 1 & age_id < 4) = pp0(26);% pp0(27);
    //    pi2(hpv2CIN_id == 1 & age_id == 4) = pp0(26) * pp0(27);% pp0(27);
    //    pi2(hpv2CIN_id == 1 & age_id == 5) = pp0(26) * pp0(28); %


    if IS(S->hpv2CIN_id, 1) {
      if LT(S->age_id, 3) {
         GET_ALL(S->pi2) = p->pp0[26]/10;
      }
      if LT(S->age_id, 4) {
        GET_ALL(S->omega2) = _p20_p21;
        if IS(S->hiv_id, 1) GET_ALL(S->omega2) = _p20_p21_p23;
        if IS(S->age_id, 3) GET_ALL(S->pi2) = p->pp0[26];
      //  GET_ALL(S->pi2) = p->pp0[26];
      }
      if GT(S->age_id, 3) {
        GET_ALL(S->omega2) = _p20_p21_o_p24;
        if IS(S->hiv_id, 1) GET_ALL(S->omega2) = _p20_p21_p23_o_p24;
        if IS(S->age_id, 4) GET_ALL(S->pi2) = _p26;  // times 2 (11-11-2020)  - removed (22-04-2021)
      //  else if IS(S->age_id, 5) GET_ALL(S->pi2) = _p26_p28;
      }
      if GT(S->age_id, 4) {
            GET_ALL(S->omega2) = _p20_p21_o_p24/20;
            if IS(S->hiv_id, 1) GET_ALL(S->omega2) = _p20_p21_p23_o_p24/20;
          }
    if GT(S->age_id2, 11) {
             GET_ALL(S->omega2) = _p20_p21_o_p24/50;
             if IS(S->hiv_id, 1) GET_ALL(S->omega2) = _p20_p21_p23_o_p24/50;
           }

     }

    //    omega3(hpv3CIN_id == 1 & age_id < 4) = pp0(20) * pp0(22);% regression from CIN2 / CIN3
    //    omega3(hpv3CIN_id == 1 & hiv_id == 1 & age_id < 4) = pp0(20) * pp0(22) * pp0(23);% regression from CIN2 / CIN3
    //    omega3(hpv3CIN_id == 1 & age_id > 3) = pp0(20) * pp0(22) / pp0(24);% regression from CIN2 / CIN3
    //    omega3(hpv3CIN_id == 1 & hiv_id == 1 & age_id > 3) = pp0(20) * pp0(22) * pp0(23) / pp0(24);% regression from CIN2 / CIN3
    //    pi3(hpv3CIN_id == 1 & age_id < 4) = pp0(26);% pp0(28);
    //    pi3(hpv3CIN_id == 1 & age_id == 4) = pp0(26) * pp0(27);% pp0(28);
    //    pi3(hpv3CIN_id == 1 & age_id == 5) = pp0(26) * pp0(28); % p

    if IS(S->hpv3CIN_id, 1) {
        if LT(S->age_id, 3) {
           GET_ALL(S->pi3) = p->pp0[26]/10;
        }
        if LT(S->age_id, 4) {
        GET_ALL(S->omega3) = _p20_p22;
        if IS(S->hiv_id, 1) GET_ALL(S->omega3) = _p20_p22_p23;
        if IS(S->age_id, 3) GET_ALL(S->pi3) = p->pp0[26];
     //   GET_ALL(S->pi3) = p->pp0[26];
      }
      if GT(S->age_id, 3) {
        GET_ALL(S->omega3) = _p20_p22_o_p24;
        if IS(S->hiv_id, 1) GET_ALL(S->omega3) = _p20_p22_p23_o_p24;
        if IS(S->age_id, 4) GET_ALL(S->pi3) = _p26; // times 2 (11-11-2020)  - removed (22-04-2021)
       // else if IS(S->age_id, 5) GET_ALL(S->pi3) = _p26_p28;
      }
      if GT(S->age_id, 4) {
        GET_ALL(S->omega3) = _p20_p22_o_p24/20;
        if IS(S->hiv_id, 1) GET_ALL(S->omega3) = _p20_p22_p23_o_p24/20;
      }
    if GT(S->age_id2, 11) {
        GET_ALL(S->omega3) = _p20_p22_o_p24/50;
        if IS(S->hiv_id, 1) GET_ALL(S->omega3) = _p20_p22_p23_o_p24/50;
      }
    }

    //////////////////////////////////////////////////////////////////

       if IS(S->hpv1CIN_id, 1) {
            if IS(S->age_id2,  9)  GET_ALL(S->pi1) = _p26_p28;
            else if IS(S->age_id2, 10)  GET_ALL(S->pi1) = _p26_p28;
            else if IS(S->age_id2, 11)  GET_ALL(S->pi1) = _p26_p28*_p27;
            else if IS(S->age_id2, 12) GET_ALL(S->pi1)  = _p26_p28*_p27;
           	else if IS(S->age_id2, 13) GET_ALL(S->pi1)  = _p26_p28*_p25;
           }

       if IS(S->hpv2CIN_id, 1) {
             if IS(S->age_id2, 9) GET_ALL(S->pi2) = _p26_p28;
             else if IS(S->age_id2, 10) GET_ALL(S->pi2) = _p26_p28;
             else if IS(S->age_id2, 11) GET_ALL(S->pi2) = _p26_p28*_p27;
             else if IS(S->age_id2, 12) GET_ALL(S->pi2) = _p26_p28*_p27;
             else if IS(S->age_id2, 13) GET_ALL(S->pi2) = _p26_p28*_p25;
           }

       if IS(S->hpv3CIN_id, 1) {
             if IS(S->age_id2, 9) GET_ALL(S->pi3) = _p26_p28;
             else if IS(S->age_id2, 10) GET_ALL(S->pi3) = _p26_p28;
             else if IS(S->age_id2, 11) GET_ALL(S->pi3) = _p26_p28*_p27;
             else if IS(S->age_id2, 12) GET_ALL(S->pi3) = _p26_p28*_p27;
             else if IS(S->age_id2, 13) GET_ALL(S->pi3) = _p26_p28*_p25;
           }

       //////////////////////////////////////////////////////////////////

    //  ups1(hpv1CC_id == 1 & age_id < 4) = pp0(29);% CC death
    //  ups1(hpv1CC_id == 1 & age_id == 4) = pp0(29) * pp0(30);% CC death
    //  ups1(hpv1CC_id == 1 & age_id == 5) = pp0(29) * pp0(31);% CC death


    if IS(S->hpv1CC_id, 1) {
      if LT(S->age_id, 4) GET_ALL(S->ups1) = p->pp0[29];
      else if IS(S->age_id, 4) GET_ALL(S->ups1) = _p29_p30;
      else if IS(S->age_id, 5) GET_ALL(S->ups1) = _p29_p31;
    }

    //  ups2(hpv2CC_id == 1 & age_id < 4) = pp0(29);
    //  ups2(hpv2CC_id == 1 & age_id == 4) = pp0(29) * pp0(30);
    //  ups2(hpv2CC_id == 1 & age_id == 5) = pp0(29) * pp0(31);

    if IS(S->hpv2CC_id, 1) {
      if LT(S->age_id, 4) GET_ALL(S->ups2) = p->pp0[29];
      else if IS(S->age_id, 4) GET_ALL(S->ups2) = _p29_p30;
      else if IS(S->age_id, 5) GET_ALL(S->ups2) = _p29_p31;

    }

    //  ups3(hpv3CC_id == 1 & age_id < 4) = pp0(29);
    //  ups3(hpv3CC_id == 1 & age_id == 4) = pp0(29) * pp0(30);
    //  ups3(hpv3CC_id == 1 & age_id == 5) = pp0(29) * pp0(31);

    if IS(S->hpv3CC_id, 1) {
      if LT(S->age_id, 4) GET_ALL(S->ups3) = p->pp0[29];
      else if IS(S->age_id, 4) GET_ALL(S->ups3) = _p29_p30;
      else if IS(S->age_id, 5) GET_ALL(S->ups3) = _p29_p31;

    }

    // Some extras to import into this loop:

    //    nu(noart_id == 1) = hp0(18);%% AIDS deaths in acuteand non - acute stage
    //    nu(hiv2_id == 1) = hp0(18);%
    //    nu(art_id == 1) = hp0(19) * hp0(18); % AIDS deaths when on ART
    //    eta(acute_id == 1) = hp0(10);% acute phase

    if IS(S->noart_id, 1) GET_ALL(S->nu) = p->hp0[18];
    if IS(S->hiv2_id, 1) GET_ALL(S->nu) = p->hp0[18];
    if IS(S->art_id, 1) GET_ALL(S->nu) = _h18_h19;
    if IS(acute_id, 1) GET_ALL(S->eta) = p->hp0[10];

  }

  //    % all three 'ups' will be applied in the differential equations, need to make
  //    % sure CC deaths are not overestimated where there is > 1 cancer, which would not happen in reality
  //    % however, the proportion of people in these
  //    % comps should be very small so likely not very influential

  //  ups1(5, 5, 5, :, 1, : , : , : ) = ups1(5, 5, 5, :, 1, : , : , : ) / 3;
  //  ups2(5, 5, 5, :, 1, : , : , : ) = ups2(5, 5, 5, :, 1, : , : , : ) / 3;
  //  ups3(5, 5, 5, :, 1, : , : , : ) = ups3(5, 5, 5, :, 1, : , : , : ) / 3;


  //  ups1(5, :, 5, : , 1, : , : , : ) = ups1(5, :, 5, : , 1, : , : , : ) / 2;
  //  ups3(5, :, 5, : , 1, : , : , : ) = ups3(5, :, 5, : , 1, : , : , : ) / 2;

  //  ups1(5, 5, :, : , 1, : , : , : ) = ups1(5, 5, :, : , 1, : , : , : ) / 2;
  //  ups2(5, 5, :, : , 1, : , : , : ) = ups2(5, 5, :, : , 1, : , : , : ) / 2;

  //  ups2(:, 5, 5, : , 1, : , : , : ) = ups2(:, 5, 5, : , 1, : , : , : ) / 2;
  //  ups3(:, 5, 5, : , 1, : , : , : ) = ups3(:, 5, 5, : , 1, : , : , : ) / 2;

  LOOP_I LOOP_R LOOP_A LOOP_V1 {
    GET(S->ups1, 5, 5, 5, _di, 1, _dr, _da, _dv1) /= 3.0;
    GET(S->ups2, 5, 5, 5, _di, 1, _dr, _da, _dv1) /= 3.0;
    GET(S->ups3, 5, 5, 5, _di, 1, _dr, _da, _dv1) /= 3.0;

    LOOP_H2 {
      GET(S->ups1, 5, _dh2, 5, _di, 1, _dr, _da, _dv1) /= 2.0;
      GET(S->ups3, 5, _dh2, 5, _di, 1, _dr, _da, _dv1) /= 2.0;
    }

    LOOP_H3 {
      GET(S->ups1, 5, 5, _dh3, _di, 1, _dr, _da, _dv1) /= 2.0;
      GET(S->ups2, 5, 5, _dh3, _di, 1, _dr, _da, _dv1) /= 2.0;
    }
    LOOP_H1 {
      GET(S->ups2, _dh1, 5, 5, _di, 1, _dr, _da, _dv1) /= 2.0;
      GET(S->ups3, _dh1, 5, 5, _di, 1, _dr, _da, _dv1) /= 2.0;
    }
  }

  //      % %VACCINATION parameters

  //    phi = pp0(43);% vaccine efficacy for hpv1
  //    phi_hiv = pp0(44);% vaccine efficacy for hpv1
  //    phi2 = pp0(45);% Vaccine cross - protection for nvtHPV

  // Moving the above 3 into diffeq_vacc where they are used.

  // vaccination waning - DEFAULT 0
       LOOP_H1 LOOP_H2 LOOP_H3 LOOP_I LOOP_S LOOP_R LOOP_A
         GET(S->vxxwane, _dh1, _dh2, _dh3, _di, _ds, _dr, _da, 2) = p->pp0[33];


  //    %% Behavior
  S->pcr = dNEW_3D_SRA(S->pcr);

  //  % 9 - 14 years old women
  //  pcr(1, 1, 1) = bp0(7);
  //  pcr(1, 2, 1) = bp0(8);
  //  pcr(1, 3, 1) = bp0(9);

  GET_3D_SRA(S->pcr, 1, 1, 1) = p->bp0[7];
  GET_3D_SRA(S->pcr, 1, 2, 1) = p->bp0[8];
  GET_3D_SRA(S->pcr, 1, 3, 1) = p->bp0[9];

  //  % 9 - 14 years old men
  //    pcr(2, 1, 1) = bp0(10);
  //  pcr(2, 2, 1) = bp0(11);
  //  pcr(2, 3, 1) = bp0(12);

  GET_3D_SRA(S->pcr, 2, 1, 1) = p->bp0[10];
  GET_3D_SRA(S->pcr, 2, 2, 1) = p->bp0[11];
  GET_3D_SRA(S->pcr, 2, 3, 1) = p->bp0[12];

//  % 15 - 24 years old women
//  new_age_structure changes 2 -> 2:3

  //  pcr(1, 1, 2:3) = bp0(13) * bp0(19);
  //  pcr(1, 2, 2:3) = bp0(14) * bp0(20);
  //  pcr(1, 3, 2:3) = bp0(15);

  //  % 15 - 24 years old men
  //  pcr(2, 1, 2:3) = bp0(16) * bp0(22);
  //  pcr(2, 2, 2:3) = bp0(17) * bp0(23);
  //  pcr(2, 3, 2:3) = bp0(18) * bp0(24);

  for (int _da = 2; _da<=3; _da++) {
    GET_3D_SRA(S->pcr, 1, 1, _da) = p->bp0[13] * p->bp0[19];
    GET_3D_SRA(S->pcr, 1, 2, _da) = p->bp0[14] * p->bp0[20];
    GET_3D_SRA(S->pcr, 1, 3, _da) = p->bp0[15];
    GET_3D_SRA(S->pcr, 2, 1, _da) = p->bp0[16] * p->bp0[22];
    GET_3D_SRA(S->pcr, 2, 2, _da) = p->bp0[17] * p->bp0[23];
    GET_3D_SRA(S->pcr, 2, 3, _da) = p->bp0[18] * p->bp0[24];
  }

  //  % 25 - 34 years old women
  //  new_age_structure changes 3 -> 4:5
  //  pcr(1, 1, 4:5) = bp0(19);
  //  pcr(1, 2, 4:5) = bp0(20);
  //  pcr(1, 3, 4:5) = bp0(21);

  //  % 25 - 34 years old men
  //  pcr(2, 1, 4:5) = bp0(22);
  //  pcr(2, 2, 4:5) = bp0(23);
  //  pcr(2, 3, 4:5) = bp0(24);

  for (int _da = 4; _da <= 5; _da++) {
    GET_3D_SRA(S->pcr, 1, 1, _da) = p->bp0[19];
    GET_3D_SRA(S->pcr, 1, 2, _da) = p->bp0[20];
    GET_3D_SRA(S->pcr, 1, 3, _da) = p->bp0[21];
    GET_3D_SRA(S->pcr, 2, 1, _da) = p->bp0[22];
    GET_3D_SRA(S->pcr, 2, 2, _da) = p->bp0[23];
    GET_3D_SRA(S->pcr, 2, 3, _da) = p->bp0[24];
  }

  //  % 35 - 49 years old women
  //  new_age_structure changes 4 -> 6:8
  //  pcr(1, 1, 6:8) = bp0(25) * bp0(19);
  //  pcr(1, 2, 6:8) = bp0(26) * bp0(20);
  //  pcr(1, 3, 6:8) = bp0(27);

  //  % 35 - 49 years old men
  //  pcr(2, 1, 6:8) = bp0(28) * bp0(22);
  //  pcr(2, 2, 6:8) = bp0(29) * bp0(23);
  //  pcr(2, 3, 6:8) = bp0(30) * bp0(24);

  for (int _da = 6; _da <= 8; _da++) {
    GET_3D_SRA(S->pcr, 1, 1, _da) = p->bp0[25] * p->bp0[19];
    GET_3D_SRA(S->pcr, 1, 2, _da) = p->bp0[26] * p->bp0[20];
    GET_3D_SRA(S->pcr, 1, 3, _da) = p->bp0[27];
    GET_3D_SRA(S->pcr, 2, 1, _da) = p->bp0[28] * p->bp0[22];
    GET_3D_SRA(S->pcr, 2, 2, _da) = p->bp0[29] * p->bp0[23];
    GET_3D_SRA(S->pcr, 2, 3, _da) = p->bp0[30] * p->bp0[24];
  }

//  % 50 - 74 years old women -> now becomes 50-59 and 60-74.
// age group 5 becomes 9:10, and 11:13 with new params...

//  % 50 - 59 years old women
//  pcr(1, 1, 9:10) = bp0(31) * bp0(19);
//  pcr(1, 2, 9:10) = bp0(32) * bp0(20);
//  pcr(1, 3, 9:10) = bp0(33);

//  % 50 - 59 years old men
//  pcr(2, 1, 9:10) = bp0(34) * bp0(22);
//  pcr(2, 2, 9:10) = bp0(35) * bp0(23);
//  pcr(2, 3, 9:10) = bp0(36) * bp0(24);

  for (int _da = 9; _da <= 10; _da++) {
    GET_3D_SRA(S->pcr, 1, 1, _da) = p->bp0[31] * p->bp0[19];
    GET_3D_SRA(S->pcr, 1, 2, _da) = p->bp0[32] * p->bp0[20];
    GET_3D_SRA(S->pcr, 1, 3, _da) = p->bp0[33];
    GET_3D_SRA(S->pcr, 2, 1, _da) = p->bp0[34] * p->bp0[22];
    GET_3D_SRA(S->pcr, 2, 2, _da) = p->bp0[35] * p->bp0[23];
    GET_3D_SRA(S->pcr, 2, 3, _da) = p->bp0[36] * p->bp0[24];
  }

//    % 60 - 74 years old women // HERE NOW ASSUMED THAT THE OLDEST AGES ALL HAVE
//    % THE LR PCR multiplier
//    pcr(1, 1, 11:13) = bp0(31) * bp0(19);
//    pcr(1, 2, 11:13) = bp0(31) * bp0(19);
//    pcr(1, 3, 11:13) = bp0(31) * bp0(19);

//  % 60 - 74 years old men // HERE NOW ASSUMED THAT THE OLDEST AGES ALL HAVE
//    % THE LR PCR
//    pcr(2, 1, 11:13) = bp0(34) * bp0(22);
//    pcr(2, 2, 11:13) = bp0(34) * bp0(22);
//    pcr(2, 3, 11:13) = bp0(34) * bp0(22);

  for (int _da = 11; _da <= 13; _da++) {
    GET_3D_SRA(S->pcr, 1, 1, _da) = p->bp0[31] * p->bp0[19];
    GET_3D_SRA(S->pcr, 1, 2, _da) = p->bp0[31] * p->bp0[19];
    GET_3D_SRA(S->pcr, 1, 3, _da) = p->bp0[31] * p->bp0[19];
    GET_3D_SRA(S->pcr, 2, 1, _da) = p->bp0[34] * p->bp0[22];
    GET_3D_SRA(S->pcr, 2, 2, _da) = p->bp0[34] * p->bp0[22];
    GET_3D_SRA(S->pcr, 2, 3, _da) = p->bp0[34] * p->bp0[22];
  }

  //  % commercial partnerships

  //  pcr_comm = zeros(s, r, a);

  S->pcr_comm = dNEW_3D_SRA(S->pcr_comm);

  // old age structure
  //  pcr_comm(1, 3, 2) = bp0(37);
  //  pcr_comm(1, 3, 3) = bp0(38);
  //  pcr_comm(1, 3, 4) = bp0(39);
  //  pcr_comm(1, 3, 5) = bp0(40);

  // new age structure
  //  pcr_comm(1, 3, 2:3) = bp0(37);
  //  pcr_comm(1, 3, 4:5) = bp0(38);
  //  pcr_comm(1, 3, 6:7) = bp0(39);
  //  pcr_comm(1, 3, 9:10) = bp0(40);
  //  pcr_comm(1, 3, 11:13) = 0;% assume FSW > 65 retire

  for (int _da = 2;  _da <= 3;  _da++) GET_3D_SRA(S->pcr_comm, 1, 3, _da) = p->bp0[37];
  for (int _da = 4;  _da <= 5;  _da++) GET_3D_SRA(S->pcr_comm, 1, 3, _da) = p->bp0[38];
  for (int _da = 6;  _da <= 7;  _da++) GET_3D_SRA(S->pcr_comm, 1, 3, _da) = p->bp0[39];
  for (int _da = 9;  _da <= 10; _da++) GET_3D_SRA(S->pcr_comm, 1, 3, _da) = p->bp0[40];
  for (int _da = 11; _da <= 11; _da++) GET_3D_SRA(S->pcr_comm, 1, 3, _da) = 0;

  // old age structure
  //  pcr_comm(2, 3, 2) = bp0(41);
  //  pcr_comm(2, 3, 3) = bp0(42);
  //  pcr_comm(2, 3, 4) = bp0(43);
  //  pcr_comm(2, 3, 5) = bp0(44);

  // new age structure
  //  pcr_comm(2, 3, 2:3) = bp0(41);
  //  pcr_comm(2, 3, 4:5) = bp0(42);
  //  pcr_comm(2, 3, 6:8) = bp0(43);
  //  pcr_comm(2, 3, 9:13) = bp0(44);

  for (int _da = 2; _da <= 3;  _da++) GET_3D_SRA(S->pcr_comm, 2, 3, _da) = p->bp0[41];
  for (int _da = 4; _da <= 5;  _da++) GET_3D_SRA(S->pcr_comm, 2, 3, _da) = p->bp0[42];
  for (int _da = 6; _da <= 8;  _da++) GET_3D_SRA(S->pcr_comm, 2, 3, _da) = p->bp0[43];
  for (int _da = 9; _da <= 13; _da++) GET_3D_SRA(S->pcr_comm, 2, 3, _da) = p->bp0[44];

  //  % number of acts
  //  % 9 - 14 years old women
  //  nacts(risk_id == 1 & sexf_id == 1 & age_id == 1) = bp0(45) / 2;
  //  nacts(risk_id == 2 & sexf_id == 1 & age_id == 1) = bp0(46) / 2;
  //  nacts(risk_id == 3 & sexf_id == 1 & age_id == 1) = bp0(47) / 2;

  //  % 9 - 14 years old men
  //  nacts(risk_id == 1 & sexm_id == 1 & age_id == 1) = bp0(45) / 2;
  //  nacts(risk_id == 2 & sexm_id == 1 & age_id == 1) = bp0(46) / 2;
  //  nacts(risk_id == 3 & sexm_id == 1 & age_id == 1) = bp0(47) / 2;

  //  % 15 - 24 years old women
  //  nacts(risk_id == 1 & sexf_id == 1 & age_id == 2) = bp0(45). / pcr(1, 1, 2);
  //  nacts(risk_id == 2 & sexf_id == 1 & age_id == 2) = bp0(46). / pcr(1, 2, 2);
  //  nacts(risk_id == 3 & sexf_id == 1 & age_id == 2) = bp0(47). / pcr(1, 3, 2);

  //   % 15 - 24 years old men
  //  nacts(risk_id == 1 & sexm_id == 1 & age_id == 2) = bp0(45). / pcr(2, 1, 2);
  //  nacts(risk_id == 2 & sexm_id == 1 & age_id == 2) = bp0(46). / pcr(2, 2, 2);
  //  nacts(risk_id == 3 & sexm_id == 1 & age_id == 2) = bp0(47). / pcr(2, 3, 2);

  //  % 25 - 34 years old women
  //  nacts(risk_id == 1 & sexf_id == 1 & age_id == 3) = bp0(45). / pcr(1, 1, 3);
  //  nacts(risk_id == 2 & sexf_id == 1 & age_id == 3) = bp0(46). / pcr(1, 2, 3);
  //  nacts(risk_id == 3 & sexf_id == 1 & age_id == 3) = bp0(47). / pcr(1, 3, 3);

  //  % 25 - 34 years old men
  //  nacts(risk_id == 1 & sexm_id == 1 & age_id == 3) = bp0(45). / pcr(2, 1, 3);
  //  nacts(risk_id == 2 & sexm_id == 1 & age_id == 3) = bp0(46). / pcr(2, 2, 3);
  //  nacts(risk_id == 3 & sexm_id == 1 & age_id == 3) = bp0(47). / pcr(2, 3, 3);

  //  % 35 - 49 years old women
  //  nacts(risk_id == 1 & sexf_id == 1 & age_id == 4) = bp0(45). / pcr(1, 1, 4);
  //  nacts(risk_id == 2 & sexf_id == 1 & age_id == 4) = bp0(46). / pcr(1, 2, 4);
  //  nacts(risk_id == 3 & sexf_id == 1 & age_id == 4) = bp0(47). / pcr(1, 3, 4);

  //  % 35 - 49 years old men
  //  nacts(risk_id == 1 & sexm_id == 1 & age_id == 4) = bp0(45). / pcr(2, 1, 4);
  //  nacts(risk_id == 2 & sexm_id == 1 & age_id == 4) = bp0(46). / pcr(2, 2, 4);
  //  nacts(risk_id == 3 & sexm_id == 1 & age_id == 4) = bp0(47). / pcr(2, 3, 4);

  //  % 50 - 74 years old women
  //  nacts(risk_id == 1 & sexf_id == 1 & age_id == 5) = bp0(45). / pcr(1, 1, 5);
  //  nacts(risk_id == 2 & sexf_id == 1 & age_id == 5) = bp0(46). / pcr(1, 2, 5);
  //  nacts(risk_id == 3 & sexf_id == 1 & age_id == 5) = bp0(47). / pcr(1, 3, 5);

  //  % 50 - 74 years old men
  //  nacts(risk_id == 1 & sexm_id == 1 & age_id == 5) = bp0(45). / pcr(2, 1, 5);
  //  nacts(risk_id == 2 & sexm_id == 1 & age_id == 5) = bp0(46). / pcr(2, 2, 5);
  //  nacts(risk_id == 3 & sexm_id == 1 & age_id == 5) = bp0(47). / pcr(2, 3, 5);

  double _b49_b50 = p->bp0[49] * p->bp0[50];

  LOOP_ALL {
    if IS(S->age_id, 1) {
      if (IS(S->sexf_id, 1) || IS(S->sexm_id, 1)) {
        if IS(S->risk_id, 1) GET_ALL(S->nacts) = p->bp0[45] / 2.0;
        else if IS(S->risk_id, 2) GET_ALL(S->nacts) = p->bp0[46] / 2.0;
        else if IS(S->risk_id, 3) GET_ALL(S->nacts) = p->bp0[47] / 2.0;
      }

    }
    else {
      if IS(S->sexf_id, 1) GET_ALL(S->nacts) =      p->bp0[45] / GET_3D_SRA(S->pcr, 1, (int) GET_ALL(S->risk_id), (int) GET_ALL(S->age_id));
      else if IS(S->sexm_id, 1) GET_ALL(S->nacts) = p->bp0[45] / GET_3D_SRA(S->pcr, 2, (int) GET_ALL(S->risk_id), (int) GET_ALL(S->age_id));
    }

    //  nacts(risk_id == 2) = bp0(46);
    //  nacts(risk_id == 3) = bp0(47);
    //  nacts_comm(risk_id == 3) = bp0(48);

    // The above overwrites what we've just done to nacts...?

    //  p_anal(risk_id == 1) = bp0(49);
    //  p_anal(risk_id == 2) = bp0(50) * bp0(49);% RR to LR
    //  p_anal(risk_id == 3) = bp0(51);
    //  pcomm_anal(risk_id == 3) = bp0(52);

    if IS(S->risk_id, 1) {
      GET_ALL(S->p_anal) = p->bp0[49];
    }

    if IS(S->risk_id, 2) {
      GET_ALL(S->nacts) = p->bp0[46];
      GET_ALL(S->p_anal) = _b49_b50;
    }
    if IS(S->risk_id, 3) {
      GET_ALL(S->nacts) = p->bp0[47];
      GET_ALL(S->nacts_comm) = p->bp0[48];
      GET_ALL(S->p_anal) = p->bp0[51];
      GET_ALL(S->pcomm_anal) = p->bp0[52];
    }

    //  %% Demography

    // old age structure:
    //    alpha(age_id == 1) = dp0(9);%% ageing
    //    alpha(age_id == 2) = dp0(10);%% ageing
    //    alpha(age_id == 3) = dp0(11);%% ageing
    //    alpha(age_id == 4) = dp0(12);%% ageing
    //    alpha(age_id == 5) = dp0(13);%% ageing

    // new age structure:
    //    alpha(age_id == 1) = dp0(9);%% ageing
    //    alpha(age_id == 2) = dp0(10);%% ageing
    //    alpha(age_id == 3) = dp0(10);%% ageing
    //    alpha(age_id == 4) = dp0(10);%% ageing
    //    alpha(age_id == 5) = dp0(10);%% ageing

    if IS(S->age_id, 1) GET_ALL(S->alpha) = p->dp0[9];
    else if ((GT(S->age_id, 1)) && (LT(S->age_id, 6))) GET_ALL(S->alpha) = p->dp0[10];

    //  mu(age_id == 1 & sexf_id == 1) = dp0(17);% bakcground deaths, women
    //  mu(age_id == 2 & sexf_id == 1) = dp0(18);
    //  mu(age_id == 3 & sexf_id == 1) = dp0(19);
    //  mu(age_id == 4 & sexf_id == 1) = dp0(20);
    //  mu(age_id == 5 & sexf_id == 1) = dp0(21);

    if IS(S->sexf_id, 1) GET_ALL(S->mu) = p->dp0[16 + (int)round(GET_ALL(S->age_id2))];

//  mu(age_id == 1 & sexm_id == 1) = dp0(17) * dp0(26);% bakcground deaths, men
//  mu(age_id == 2 & sexm_id == 1) = dp0(18) * dp0(26);
//  mu(age_id == 3 & sexm_id == 1) = dp0(19) * dp0(26);
//  mu(age_id == 4 & sexm_id == 1) = dp0(20) * dp0(26);
//  mu(age_id == 5 & sexm_id == 1) = dp0(21) * dp0(26);

    if IS(S->sexm_id, 1) GET_ALL(S->mu) = p->dp0[16 + (int)round(GET_ALL(S->age_id2))] * p->dp0[30];
  }
//  P = zeros(h1, h2, h3, i, s, r, a, v1);
//  % h1, h2, h3, i, s, r, a, v1

  dMAT_8D P = dNEW_8D(P);

//    % I have implemented initial population in two ways, either by
//    % starting from close to HPV equilibrium(implemented now, below), and keeping
//    % population in risk groups constant ...

//    prsk3_f = bp0(6);% atm constant
//    prsk2_f = bp0(5);% atm constant
//    prsk3_m = bp0(3); % ''
//    prsk2_m = bp0(2); % ''

  #define prsk3_f p->bp0[6]
  #define prsk2_f p->bp0[5]
  #define prsk3_m p->bp0[3]
  #define prsk2_m p->bp0[2]

// basepop = load('popstart.mat');
// P(:, : , : , : , : , : , : , 1) = basepop.popn;

  DArray basepop_1950 = CREATE_DARRAY(_SIZE_A, basepop_1950, 0.0);

  read_popfile(pop_file, basepop_1950);

  // new_age_structure: pop_file will now contain 13 numbers (basepop.pop1950)

  // New work needed to be done on P now, as it varies with parameters.



  //% ... or by varying the distribution of risk groups, and seeding HPV at a
  //% low level(will be further from equilibrium)

  //% P(1, 1, 1, 1, 1, 1, 1, 1) = 911999;
  //% P(1, 1, 1, 1, 1, 1, 2, 1) = 1230001;
  //% P(1, 1, 1, 1, 1, 1, 3, 1) = 980000;
  //% P(1, 1, 1, 1, 1, 1, 4, 1) = 1020000;
  //% P(1, 1, 1, 1, 1, 1, 5, 1) = 838965;
  //%
  //% P(1, 1, 1, 1, 2, 1, 1, 1) = 913800;
  //% P(1, 1, 1, 1, 2, 1, 2, 1) = 1306135;
  //% P(1, 1, 1, 1, 2, 1, 3, 1) = 1034000;
  //% P(1, 1, 1, 1, 2, 1, 4, 1) = 1085999;
  //% P(1, 1, 1, 1, 2, 1, 5, 1) = 751906;

  //% P(1, 1, 1, 1, 1, 3, :, 1) = prsk3_f.*P(1, 1, 1, 1, 1, 1, :, 1);
  //% P(1, 1, 1, 1, 1, 2, :, 1) = prsk2_f.*P(1, 1, 1, 1, 1, 1, :, 1);
  //% P(1, 1, 1, 1, 1, 1, :, 1) = P(1, 1, 1, 1, 1, 1, :, 1) - (prsk3_f + prsk2_f).*P(1, 1, 1, 1, 1, 1, :, 1);
  //%
  //% P(1, 1, 1, 1, 2, 3, :, 1) = prsk3_m.*P(1, 1, 1, 1, 2, 1, :, 1);
  //% P(1, 1, 1, 1, 2, 2, :, 1) = prsk2_m.*P(1, 1, 1, 1, 2, 1, :, 1);
  //% P(1, 1, 1, 1, 2, 1, :, 1) = P(1, 1, 1, 1, 2, 1, :, 1) - (prsk3_m + prsk2_m).*P(1, 1, 1, 1, 2, 1, :, 1);

  //% Check initial conditions are correct
  //  % sum(sum(sum(sum(sum(sum(sum(P(:, : , : , : , 1, : , : , : ))))))))
  //  % 9475286
  //  % sum(sum(sum(sum(sum(sum(sum(P(:, : , : , : , 1, 1, : , : )))))))) / sum(sum(sum(sum(sum(sum(sum(P(:, : , : , : , 1, : , : , : ))))))))
  //  % sum(sum(sum(sum(sum(sum(sum(P(:, : , : , : , 1, 2, : , : )))))))) / sum(sum(sum(sum(sum(sum(sum(P(:, : , : , : , 1, : , : , : ))))))))
  //  % sum(sum(sum(sum(sum(sum(sum(P(:, : , : , : , 1, 3, : , : )))))))) / sum(sum(sum(sum(sum(sum(sum(P(:, : , : , : , 1, : , : , : ))))))))

  //  % sum(sum(sum(sum(sum(sum(sum(P(:, : , : , : , 2, : , : , : ))))))))
  //  % 9471386
  //  % sum(sum(sum(sum(sum(sum(sum(P(:, : , : , : , 2, 1, : , : )))))))) / sum(sum(sum(sum(sum(sum(sum(P(:, : , : , : , 2, : , : , : ))))))))
  //  % sum(sum(sum(sum(sum(sum(sum(P(:, : , : , : , 2, 2, : , : )))))))) / sum(sum(sum(sum(sum(sum(sum(P(:, : , : , : , 2, : , : , : ))))))))
  //  % sum(sum(sum(sum(sum(sum(sum(P(:, : , : , : , 2, 3, : , : )))))))) / sum(sum(sum(sum(sum(sum(sum(P(:, : , : , : , 2, : , : , : ))))))))
  //  %
  //  %%seed HPV

  LOOP_A {
    // P(1, 1, 1, 1, 1, 1, :, 1) = (1 - prsk2_f - prsk3_f) * basepop.pop1950 / s;
    // P(1, 1, 1, 1, 1, 2, :, 1) = prsk2_f * basepop.pop1950 / s;
    // P(1, 1, 1, 1, 1, 3, :, 1) = prsk3_f * basepop.pop1950 / s;

    GET(P, 1, 1, 1, 1, 1, 1, _da, 1) = (1 - prsk2_f - prsk3_f) * basepop_1950[_da] / _SIZE_S;
    GET(P, 1, 1, 1, 1, 1, 2, _da, 1) = prsk2_f * basepop_1950[_da] / _SIZE_S;
    GET(P, 1, 1, 1, 1, 1, 3, _da, 1) = prsk3_f * basepop_1950[_da] / _SIZE_S;

    // P(1, 1, 1, 1, 2, 1, :, 1) = (1 - prsk2_m - prsk3_m) * basepop.pop1950 / s;
    // P(1, 1, 1, 1, 2, 2, :, 1) = prsk2_m * basepop.pop1950 / s;
    // P(1, 1, 1, 1, 2, 3, :, 1) = prsk3_m * basepop.pop1950 / s;


    GET(P, 1, 1, 1, 1, 2, 1, _da, 1) = (1 - prsk2_m - prsk3_m) * basepop_1950[_da] / _SIZE_S;
    GET(P, 1, 1, 1, 1, 2, 2, _da, 1) = prsk2_m * basepop_1950[_da] / _SIZE_S;
    GET(P, 1, 1, 1, 1, 2, 3, _da, 1) = prsk3_m * basepop_1950[_da] / _SIZE_S;

    LOOP_R {
      LOOP_S {

        // P(2, 1, 1, 1, :, : , : , 1) = 0.08.*P(1, 1, 1, 1, :, : , : ,1);
        // P(1, 2, 1, 1, :, : , : , 1) = 0.08.*P(1, 1, 1, 1, :, : , : ,1);
        // P(1, 1, 2, 1, :, : , : , 1) = 0.08.*P(1, 1, 1, 1, :, : , : ,1);
        // P(1, 1, 1, 1, :, : , : , 1) = P(1, 1, 1, 1, :, : , : , 1) - 3 * 0.08 * P(1, 1, 1, 1, :, : , : , 1);

        GET(P, 2, 1, 1, 1, _ds, _dr, _da, 1) = 0.08 * GET(P, 1, 1 ,1, 1, _ds, _dr, _da, 1);
        GET(P, 1, 2, 1, 1, _ds, _dr, _da, 1) = 0.08 * GET(P, 1, 1, 1, 1, _ds, _dr, _da, 1);
        GET(P, 1, 1, 2, 1, _ds, _dr, _da, 1) = 0.08 * GET(P, 1, 1, 1, 1, _ds, _dr, _da, 1);
        GET(P, 1, 1, 1, 1, _ds, _dr, _da, 1) = GET(P, 1, 1, 1, 1, _ds, _dr, _da, 1) - 3 * 0.08 * GET(P, 1, 1, 1, 1, _ds, _dr, _da, 1);

        // P(3, 1, 1, 1, :, : , : , 1) = 0.03.*P(1, 1, 1, 1, :, : , : , 1);
        // P(1, 3, 1, 1, :, : , : , 1) = 0.03.*P(1, 1, 1, 1, :, : , : , 1);
        // P(1, 1, 3, 1, :, : , : , 1) = 0.03.*P(1, 1, 1, 1, :, : , : , 1);
        // P(1, 1, 1, 1, :, : , : , 1) = P(1, 1, 1, 1, :, : , : , 1) - 3 * 0.03 * P(1, 1, 1, 1, :, : , : , 1);

        GET(P, 3, 1, 1, 1, _ds, _dr, _da, 1) = 0.03 * GET(P, 1, 1, 1, 1, _ds, _dr, _da, 1);
        GET(P, 1, 3, 1, 1, _ds, _dr, _da, 1) = 0.03 * GET(P, 1, 1, 1, 1, _ds, _dr, _da, 1);
        GET(P, 1, 1, 3, 1, _ds, _dr, _da, 1) = 0.03 * GET(P, 1, 1, 1, 1, _ds, _dr, _da, 1);
        GET(P, 1, 1, 1, 1, _ds, _dr, _da, 1) = GET(P, 1, 1, 1, 1, _ds, _dr, _da, 1) - 3 * 0.03 * GET(P, 1, 1, 1, 1, _ds, _dr, _da, 1);

      } // end of S loop

      // P(4, 1, 1, 1, 1, :, : , 1) = 0.01.*P(1, 1, 1, 1, 1, :, : , 1);
      // P(1, 4, 1, 1, 1, :, : , 1) = 0.01.*P(1, 1, 1, 1, 1, :, : , 1);
      // P(1, 1, 4, 1, 1, :, : , 1) = 0.01.*P(1, 1, 1, 1, 1, :, : , 1);
      // P(1, 1, 1, 1, 1, :, : , 1) = P(1, 1, 1, 1, 1, :, : , 1) - 3 * 0.01 * P(1, 1, 1, 1, 1, :, : , 1);

      GET(P, 4, 1, 1, 1, 1, _dr, _da, 1) = 0.01 * GET(P, 1, 1, 1, 1, 1, _dr, _da, 1);
      GET(P, 1, 4, 1, 1, 1, _dr, _da, 1) = 0.01 * GET(P, 1, 1, 1, 1, 1, _dr, _da, 1);
      GET(P, 1, 1, 4, 1, 1, _dr, _da, 1) = 0.01 * GET(P, 1, 1, 1, 1, 1, _dr, _da, 1);
      GET(P, 1, 1, 1, 1, 1, _dr, _da, 1) = GET(P, 1, 1, 1, 1, 1, _dr, _da, 1) - 3 * 0.01 * GET(P, 1, 1, 1, 1, 1, _dr, _da, 1);

      // P(5, 1, 1, 1, 1, :, : , 1) = 0.001.*P(1, 1, 1, 1, 1, :, : , 1);
      // P(1, 5, 1, 1, 1, :, : , 1) = 0.001.*P(1, 1, 1, 1, 1, :, : , 1);
      // P(1, 1, 5, 1, 1, :, : , 1) = 0.001.*P(1, 1, 1, 1, 1, :, : , 1);
      // P(1, 1, 1, 1, 1, :, : , 1) = P(1, 1, 1, 1, 1, :, : , 1) - 3 * 0.001 * P(1, 1, 1, 1, 1, :, : , 1);

      GET(P, 5, 1, 1, 1, 1, _dr, _da, 1) = 0.001 * GET(P, 1, 1, 1, 1, 1, _dr, _da, 1);
      GET(P, 1, 5, 1, 1, 1, _dr, _da, 1) = 0.001 * GET(P, 1, 1, 1, 1, 1, _dr, _da, 1);
      GET(P, 1, 1, 5, 1, 1, _dr, _da, 1) = 0.001 * GET(P, 1, 1, 1, 1, 1, _dr, _da, 1);
      GET(P, 1, 1, 1, 1, 1, _dr, _da, 1) = GET(P, 1, 1, 1, 1, 1, _dr, _da, 1) - 3 * 0.001 * GET(P, 1, 1, 1, 1, 1, _dr, _da, 1);
    } // End R
  }

//  %% Call function where Euler solver

//    csize = h1 * h2 * h3 * i * s * r * a * v1;

//  int csize = _h1 * _h2 * _h3 * _i * _s * _r * _a * _v1;

//  pop0 = P;
//  %% matrix P into a new initial population column vector called pop0 - used in ode45 command

  dMAT_8D pop0 = dNEW_COPY_8D(pop0, P);


//  T0 = dp0(31); % 1950;
//  T1 = dp0(32); % 150; %
//  tstep = dp0(33);% timestep
//  tspan = (round(tstep / 2)) : tstep : tstep * T1;

  // - Move this to diffeq_vacc, as only need it there
  // - #define T0 p->dp0[31]

  #define T1 (int) round(p->dp0[32])
  #define tstep (int) round(p->dp0[33])

  // tpsan is never used.


//    [popbyrisk_f, popbyrisk_m, hiv_9to14, hiv_15to49, hiv_15to74, hivart_9to14, hivart_15to49, hivart_15to74, ...
//    art_15plus, hiv_fage, hiv_mage, art_age, pop1980, pop2018, pop2120, hiv_fsw, hiv_mcli, ...
//    psize, psize_age_f, psize_age_f_neg, psize_age_f_pos, psize_age_f_art, ...
//    psize_age_m, psize_age_m_neg, psize_age_m_pos, psize_age_m_art, vacc_age_f, vacc_age_m, vacc_age_f_pos, vacc_age_m_pos, ...
//    hivd_f, hivd_m, ccd_f, ccd_f_neg, ccd_f_pos, ccd_f_art, othd_f, othd_m, hpv1cum_pop, hpv2cum_pop, hpv3cum_pop, ...
//    hpv1hncum_pop, hpv2hncum_pop, hpv3hncum_pop, ...
//    hivsus_pop, hivinc_pop_cum, ...
//    ccinc_f, ccinc_f_neg, ccinc_f_pos, ccinc_f_art, ccinc_nvtf, ccinc_nvtf_neg, ccinc_nvtf_pos, ccinc_nvtf_art, cinin_f, cinin_f_neg, cinin_f_pos, cinin_f_art, ...
//    hpv1618_f, hpv1618_m, nvthpv_f, nvthpv_m, hpv9vt_f, hpv9vt_m, hpv_f, hpv_m, cinpr_f, ccpr_f, screencov_yr, ...
//    errmessage, testtime, testtime1] = diffeq_vacc(T1, tstep, pop0, dp0, bp0, hp0, pp0, vaccscen, screenscen);

  /* TEST 2  */
  #ifdef TEST_OUTPUT
  LOOP_ALL fwrite(&GET_ALL(S->betahpv1), 8, 1, fout);
  LOOP_ALL fwrite(&GET_ALL(S->betahpv2), 8, 1, fout);
  LOOP_ALL fwrite(&GET_ALL(S->betahpv3), 8, 1, fout);
  LOOP_ALL fwrite(&GET_ALL(S->betahivacute), 8, 1, fout);
  LOOP_ALL fwrite(&GET_ALL(S->betahiv2), 8, 1, fout);
  LOOP_ALL fwrite(&GET_ALL(S->betahivart), 8, 1, fout);
  LOOP_ALL fwrite(&GET_ALL(S->betahivacute_ai), 8, 1, fout);
  LOOP_ALL fwrite(&GET_ALL(S->betahiv2_ai), 8, 1, fout);
  LOOP_ALL fwrite(&GET_ALL(S->betahivart_ai), 8, 1, fout);
  LOOP_ALL fwrite(&GET_ALL(S->delta1), 8, 1, fout);
  LOOP_ALL fwrite(&GET_ALL(S->delta2), 8, 1, fout);
  LOOP_ALL fwrite(&GET_ALL(S->delta3), 8, 1, fout);
  LOOP_ALL fwrite(&GET_ALL(S->vxxwane), 8, 1, fout);
  LOOP_ALL fwrite(&GET_ALL(S->nacts), 8, 1, fout);
  LOOP_ALL fwrite(&GET_ALL(S->nacts_comm), 8, 1, fout);
  LOOP_ALL fwrite(&GET_ALL(S->p_anal), 8, 1, fout);
  LOOP_ALL fwrite(&GET_ALL(S->pcomm_anal), 8, 1, fout);
  LOOP_ALL fwrite(&GET_ALL(S->mu), 8, 1, fout);
  LOOP_ALL fwrite(&GET_ALL(S->nu), 8, 1, fout);
  LOOP_ALL fwrite(&GET_ALL(S->omega1), 8, 1, fout);
  LOOP_ALL fwrite(&GET_ALL(S->omega2), 8, 1, fout);
  LOOP_ALL fwrite(&GET_ALL(S->omega3), 8, 1, fout);
  LOOP_ALL fwrite(&GET_ALL(S->eta), 8, 1, fout);
  LOOP_ALL fwrite(&GET_ALL(S->tau), 8, 1, fout);
  LOOP_ALL fwrite(&GET_ALL(S->omikron), 8, 1, fout);
  LOOP_ALL fwrite(&GET_ALL(S->theta1), 8, 1, fout);
  LOOP_ALL fwrite(&GET_ALL(S->theta2), 8, 1, fout);
  LOOP_ALL fwrite(&GET_ALL(S->theta3), 8, 1, fout);
  LOOP_ALL fwrite(&GET_ALL(S->alpha), 8, 1, fout);
  LOOP_ALL fwrite(&GET_ALL(S->sigma1), 8, 1, fout);
  LOOP_ALL fwrite(&GET_ALL(S->sigma2), 8, 1, fout);
  LOOP_ALL fwrite(&GET_ALL(S->sigma3), 8, 1, fout);
  LOOP_ALL fwrite(&GET_ALL(S->psi1), 8, 1, fout);
  LOOP_ALL fwrite(&GET_ALL(S->psi2), 8, 1, fout);
  LOOP_ALL fwrite(&GET_ALL(S->psi3), 8, 1, fout);
  LOOP_ALL fwrite(&GET_ALL(S->pi1), 8, 1, fout);
  LOOP_ALL fwrite(&GET_ALL(S->pi2), 8, 1, fout);
  LOOP_ALL fwrite(&GET_ALL(S->pi3), 8, 1, fout);
  LOOP_ALL fwrite(&GET_ALL(S->ups1), 8, 1, fout);
  LOOP_ALL fwrite(&GET_ALL(S->ups2), 8, 1, fout);
  LOOP_ALL fwrite(&GET_ALL(S->ups3), 8, 1, fout);
  LOOP_ALL fwrite(&GET_ALL(S->hups1), 8, 1, fout);
  LOOP_ALL fwrite(&GET_ALL(S->hups2), 8, 1, fout);
  LOOP_ALL fwrite(&GET_ALL(S->hups3), 8, 1, fout);
  LOOP_ALL fwrite(&GET_ALL(P), 8, 1, fout);
  LOOP_ALL fwrite(&GET_ALL(pop0), 8, 1, fout);

  fclose(fout);
  #endif

  results* R = diffeq_vacc(T1, tstep, pop0, p, S, dims);

//    % Defining outputs

//    for y = 2:dp0(23)+1

  int* dims_y_a =  new int[3] {0, (int)round(p->dp0[32]) + 1, _SIZE_A};
  R->deatht = CREATE_DMATRIX(dims_y_a, R->deatht, 0.0)
  R->hivdftot = CREATE_DMATRIX(dims_y_a, R->hivdftot, 0.0)
  R->hivdmtot = CREATE_DMATRIX(dims_y_a, R->hivdmtot, 0.0)
  R->ccdftot = CREATE_DMATRIX(dims_y_a, R->ccdftot, 0.0)
  R->ccdftot_neg = CREATE_DMATRIX(dims_y_a, R->ccdftot_neg, 0.0)
  R->ccdftot_pos = CREATE_DMATRIX(dims_y_a, R->ccdftot_pos, 0.0)
  R->ccdftot_art = CREATE_DMATRIX(dims_y_a, R->ccdftot_art, 0.0)
  R->ccinc = CREATE_DMATRIX(dims_y_a, R->ccinc, 0.0)
  R->ccinc_neg = CREATE_DMATRIX(dims_y_a, R->ccinc_neg, 0.0)
  R->ccinc_pos = CREATE_DMATRIX(dims_y_a, R->ccinc_pos, 0.0)
  R->ccinc_art = CREATE_DMATRIX(dims_y_a, R->ccinc_art, 0.0)
  R->ccinc_nvt = CREATE_DMATRIX(dims_y_a, R->ccinc_nvt, 0.0)
  R->ccinc_nvt_neg = CREATE_DMATRIX(dims_y_a, R->ccinc_nvt_neg, 0.0)
  R->ccinc_nvt_pos = CREATE_DMATRIX(dims_y_a, R->ccinc_nvt_pos, 0.0)
  R->ccinc_nvt_art = CREATE_DMATRIX(dims_y_a, R->ccinc_nvt_art, 0.0)
  R->cin2in = CREATE_DMATRIX(dims_y_a, R->cin2in, 0.0)
  R->cin2in_neg = CREATE_DMATRIX(dims_y_a, R->cin2in_neg, 0.0)
  R->cin2in_pos = CREATE_DMATRIX(dims_y_a, R->cin2in_pos, 0.0)
  R->cin2in_art = CREATE_DMATRIX(dims_y_a, R->cin2in_art, 0.0)
  int* dims_y_s_a_v = new int[5] {0, (int)round(p->dp0[32] + 1), dims[DIM_s], dims[DIM_a], dims[DIM_v1]};
  R->hivinctot = CREATE_DMATRIX4D(dims_y_s_a_v, R->hivinctot, 0.0)
  int* dims_y32 = new int[4] {0, (int)round(p->dp0[32] + 1), 3, 2};

 // R->screencovf = CREATE_DMATRIX3D(dims_y32, R->screencovf, 0.0)
  R->screencovf_neg = CREATE_DMATRIX3D(dims_y32, R->screencovf_neg, 0.0)
  R->screencovf_pos = CREATE_DMATRIX3D(dims_y32, R->screencovf_pos, 0.0)
  R->screencovf_art = CREATE_DMATRIX3D(dims_y32, R->screencovf_art, 0.0)

  R->hpv1hncum = CREATE_DMATRIX(R->dim_y_sav, R->hpv1hncum, 0.0)
  R->hpv2hncum = CREATE_DMATRIX(R->dim_y_sav, R->hpv2hncum, 0.0)
  R->hpv3hncum = CREATE_DMATRIX(R->dim_y_sav, R->hpv3hncum, 0.0)
  R->hpv1cum = CREATE_DMATRIX(R->dim_y_sav, R->hpv1cum, 0.0)
  R->hpv2cum = CREATE_DMATRIX(R->dim_y_sav, R->hpv2cum, 0.0)
  R->hpv3cum = CREATE_DMATRIX(R->dim_y_sav, R->hpv3cum, 0.0)

  for (int y = 2; y <= (p->dp0[32] + 1); y++) {
    // Inserted in new age structure:
  // hpv1hncum(y, :) = hpv1hncum_pop(y, :) - hpv1hncum_pop(y - 1, :);
  // hpv2hncum(y, :) = hpv2hncum_pop(y, :) - hpv2hncum_pop(y - 1, :);
  // hpv3hncum(y, :) = hpv3hncum_pop(y, :) - hpv3hncum_pop(y - 1, :);
  // hpv1cum(y, :) = hpv1cum_pop(y, :) - hpv1cum_pop(y - 1, :);
  // hpv2cum(y, :) = hpv2cum_pop(y, :) - hpv2cum_pop(y - 1, :);
  // hpv3cum(y, :) = hpv3cum_pop(y, :) - hpv3cum_pop(y - 1, :);
    for (int _sav=1; _sav <= _SIZE_S * _SIZE_A * _SIZE_V1; _sav++) {
      R->hpv1hncum[y][_sav] = R->hpv1hncum_pop[y][_sav] - R->hpv1hncum_pop[y - 1][_sav];
      R->hpv2hncum[y][_sav] = R->hpv2hncum_pop[y][_sav] - R->hpv2hncum_pop[y - 1][_sav];
      R->hpv3hncum[y][_sav] = R->hpv3hncum_pop[y][_sav] - R->hpv3hncum_pop[y - 1][_sav];

      R->hpv1cum[y][_sav] = R->hpv1cum_pop[y][_sav] - R->hpv1cum_pop[y - 1][_sav];
      R->hpv2cum[y][_sav] = R->hpv2cum_pop[y][_sav] - R->hpv2cum_pop[y - 1][_sav];
      R->hpv3cum[y][_sav] = R->hpv3cum_pop[y][_sav] - R->hpv3cum_pop[y - 1][_sav];
    }

    LOOP_A {

//    deatht(y, :) = hivd_f(y, :)     + hivd_m(y, :)     + ccd_f(y, :)     + othd_f(y, :)     + othd_m(y, :)
//                 - hivd_f(y - 1, :) - hivd_m(y - 1, :) - ccd_f(y - 1, :) - othd_f(y - 1, :) - othd_m(y - 1, :);

      R->deatht[y][_da] = R->hivd_f[y][_da]     + R->hivd_m[y][_da]     + R->ccd_f[y][_da]     + R->othd_f[y][_da]     + R->othd_m[y][_da]
                      - R->hivd_f[y - 1][_da] - R->hivd_m[y - 1][_da] - R->ccd_f[y - 1][_da] - R->othd_f[y - 1][_da] - R->othd_m[y - 1][_da];

//    hivdftot(y, :) = hivd_f(y, :) - hivd_f(y - 1, :);
//    hivdmtot(y, :) = hivd_m(y, :) - hivd_m(y - 1, :);

      R->hivdftot[y][_da] = R->hivd_f[y][_da] - R->hivd_f[y - 1][_da];
      R->hivdmtot[y][_da] = R->hivd_m[y][_da] - R->hivd_m[y - 1][_da];

//    ccdftot(y, :) = ccd_f(y, :) - ccd_f(y - 1, :);
//    ccdftot_neg(y, :) = ccd_f_neg(y, :) - ccd_f_neg(y - 1, :);
//    ccdftot_pos(y, :) = ccd_f_pos(y, :) - ccd_f_pos(y - 1, :);
//    ccdftot_art(y, :) = ccd_f_art(y, :) - ccd_f_art(y - 1, :);

      R->ccdftot[y][_da] = R->ccd_f[y][_da] - R->ccd_f[y - 1][_da];
      R->ccdftot_neg[y][_da] = R->ccd_f_neg[y][_da] - R->ccd_f_neg[y - 1][_da];
      R->ccdftot_pos[y][_da] = R->ccd_f_pos[y][_da] - R->ccd_f_pos[y - 1][_da];
      R->ccdftot_art[y][_da] = R->ccd_f_art[y][_da] - R->ccd_f_art[y - 1][_da];

      //    ccinc(y, :) = ccinc_f(y, :) - ccinc_f(y - 1, :);
      //    ccinc_neg(y, :) = ccinc_f_neg(y, :) - ccinc_f_neg(y - 1, :);
      //    ccinc_pos(y, :) = ccinc_f_pos(y, :) - ccinc_f_pos(y - 1, :);
      //    ccinc_art(y, :) = ccinc_f_art(y, :) - ccinc_f_art(y - 1, :);

      R->ccinc[y][_da] = R->ccinc_f[y][_da] - R->ccinc_f[y - 1][_da];
      R->ccinc_neg[y][_da] = R->ccinc_f_neg[y][_da] - R->ccinc_f_neg[y - 1][_da];
      R->ccinc_pos[y][_da] = R->ccinc_f_pos[y][_da] - R->ccinc_f_pos[y - 1][_da];
      R->ccinc_art[y][_da] = R->ccinc_f_art[y][_da] - R->ccinc_f_art[y - 1][_da];

      //    ccinc_nvt(y, :) = ccinc_nvtf(y, :) - ccinc_nvtf(y - 1, :);
      //    ccinc_nvt_neg(y, :) = ccinc_nvtf_neg(y, :) - ccinc_nvtf_neg(y - 1, :);
      //    ccinc_nvt_pos(y, :) = ccinc_nvtf_pos(y, :) - ccinc_nvtf_pos(y - 1, :);
      //    ccinc_nvt_art(y, :) = ccinc_nvtf_art(y, :) - ccinc_nvtf_art(y - 1, :);

      R->ccinc_nvt[y][_da] =     R->ccinc_nvtf[y][_da] -     R->ccinc_nvtf[y - 1][_da];
      R->ccinc_nvt_neg[y][_da] = R->ccinc_nvtf_neg[y][_da] - R->ccinc_nvtf_neg[y - 1][_da];
      R->ccinc_nvt_pos[y][_da] = R->ccinc_nvtf_pos[y][_da] - R->ccinc_nvtf_pos[y - 1][_da];
      R->ccinc_nvt_art[y][_da] = R->ccinc_nvtf_art[y][_da] - R->ccinc_nvtf_art[y - 1][_da];

//    cin2in(y, :) = cinin_f(y, :) - cinin_f(y - 1, :);
//    cin2in_neg(y, :) = cinin_f_neg(y, :) - cinin_f_neg(y - 1, :);
//    cin2in_pos(y, :) = cinin_f_pos(y, :) - cinin_f_pos(y - 1, :);
//    cin2in_art(y, :) = cinin_f_art(y, :) - cinin_f_art(y - 1, :);

      R->cin2in[y][_da] = R->cinin_f[y][_da] - R->cinin_f[y - 1][_da];
      R->cin2in_neg[y][_da] = R->cinin_f_neg[y][_da] - R->cinin_f_neg[y - 1][_da];
      R->cin2in_pos[y][_da] = R->cinin_f_pos[y][_da] - R->cinin_f_pos[y - 1][_da];
      R->cin2in_art[y][_da] = R->cinin_f_art[y][_da] - R->cinin_f_art[y - 1][_da];
    }

//  hivinctot(y, :, : , : ) = hivinc_pop_cum(y, :, : , : ) - hivinc_pop_cum(y - 1, :, : , : );

    LOOP_S LOOP_A LOOP_V1
         R->hivinctot[y][_ds][_da][_dv1] = R->hivinc_pop_cum[y][_ds][_da][_dv1] - R->hivinc_pop_cum[y - 1][_ds][_da][_dv1];

//  screencovf(y, :, : ) = screencov_yr(y, :, : ) - screencov_yr(y - 1, :, : );

    for (int j=1; j<=3; j++)
      for (int k=1; k<=2; k++)
        R->screencovf_neg[y][j][k] = R->screencov_neg_yr[y][j][k] - R->screencov_neg_yr[y - 1][j][k];
    for (int j=1; j<=3; j++)
      for (int k=1; k<=2; k++)
    	R->screencovf_pos[y][j][k] = R->screencov_pos_yr[y][j][k] - R->screencov_pos_yr[y - 1][j][k];
    for (int j=1; j<=3; j++)
      for (int k=1; k<=2; k++)
    	R->screencovf_art[y][j][k] = R->screencov_art_yr[y][j][k] - R->screencov_art_yr[y - 1][j][k];


 // end

  }



  // So, these are the final outputs - what do we want to do with them. Should deal with that here I think, rather
  // than building another object and returning it all to main.cpp...
  // Let's write them to a bin file, and then write similar code in matlab to retrieve them.

  FILE* f = fopen(output_file, "wb");

//  output.popbyrisk_f = popbyrisk_f;
//  output.popbyrisk_m = popbyrisk_m;

  writeMatrix(f, R->popbyrisk_f, R->dim_g_r);
  writeMatrix(f, R->popbyrisk_m, R->dim_g_r);

//  output.hivart_9to14 = hivart_9to14;
//  output.hivart_15to49 = hivart_15to49;
//  output.hivart_15to74 = hivart_15to74;

  writeMatrix(f, R->hivart_9to14, R->dim_gg_s);
  writeMatrix(f, R->hivart_15to49, R->dim_gg_s);
  writeMatrix(f, R->hivart_15to74, R->dim_gg_s);


  //output.hiv_9to14 = hiv_9to14;
  //output.hiv_15to49 = hiv_15to49;
  //output.hiv_15to74 = hiv_15to74;
  //output.art_15plus = art_15plus;

  writeMatrix(f, R->hiv_9to14, R->dim_gg_s);
  writeMatrix(f, R->hiv_15to49, R->dim_gg_s);
  writeMatrix(f, R->hiv_15to74, R->dim_gg_s);
  writeMatrix(f, R->art_15plus, R->dim_gg_s);

  //output.hiv_fage = hiv_fage;
  //output.hiv_mage = hiv_mage;
  //output.art_age = art_age;

  writeMatrix(f, R->hiv_fage, R->dim_gg_a5);
  writeMatrix(f, R->hiv_mage, R->dim_gg_a5);
  writeMatrix(f, R->art_age, R->dim_gg_a5);

  //output.pop1980 = pop1980;
  //output.pop2018 = pop2018;
  //output.pop2120 = pop2120;

  write8D(f, R->pop2000);
  write8D(f, R->pop2018);
  write8D(f, R->pop2120);

  //output.psize = psize;

  writeMatrix(f, R->psize, R->dim_g_s);

//  output.psize_age_f = psize_age_f;
//  output.psize_age_f_neg = psize_age_f_neg;
//  output.psize_age_f_pos = psize_age_f_pos;
//  output.psize_age_f_art = psize_age_f_art;
//  output.psize_age_m = psize_age_m;
//  output.psize_age_m_neg = psize_age_m_neg;
//  output.psize_age_m_pos = psize_age_m_pos;
//  output.psize_age_m_art = psize_age_m_art;

  writeMatrix(f, R->psize_age_f, R->dim_g_a);
  writeMatrix(f, R->psize_age_f_neg, R->dim_g_a);
  writeMatrix(f, R->psize_age_f_pos, R->dim_g_a);
  writeMatrix(f, R->psize_age_f_art, R->dim_g_a);
  writeMatrix(f, R->psize_age_m, R->dim_g_a);
  writeMatrix(f, R->psize_age_m_neg, R->dim_g_a);
  writeMatrix(f, R->psize_age_m_pos, R->dim_g_a);
  writeMatrix(f, R->psize_age_m_art, R->dim_g_a);

  //output.vacc_age_f = vacc_age_f;
  //output.vacc_age_m = vacc_age_m;
  //output.vacc_age_f_pos = vacc_age_f_pos;
  //output.vacc_age_m_pos = vacc_age_m_pos;

  writeMatrix(f, R->vacc_age_f, R->dim_gg_a);
  writeMatrix(f, R->vacc_age_m, R->dim_gg_a);
  writeMatrix(f, R->vacc_age_f_pos, R->dim_gg_a);
  writeMatrix(f, R->vacc_age_f_art, R->dim_gg_a);
  writeMatrix(f, R->vacc_age_m_pos, R->dim_gg_a);

  // output.deatht = deatht;
  // output.hivdftot = hivdftot;
  // output.hivdmtot = hivdmtot;
  //output.ccdftot = ccdftot;
  //output.ccdftot_neg = ccdftot_neg;
  //output.ccdftot_pos = ccdftot_pos;
  //output.ccdftot_art = ccdftot_art;

  writeMatrix(f, R->deatht, dims_y_a);
  writeMatrix(f, R->hivdftot, dims_y_a);
  writeMatrix(f, R->hivdmtot, dims_y_a);

  writeMatrix(f, R->ccdftot, dims_y_a);
  writeMatrix(f, R->ccdftot_neg, dims_y_a);
  writeMatrix(f, R->ccdftot_pos, dims_y_a);
  writeMatrix(f, R->ccdftot_art, dims_y_a);

  //output.hpv1cum_pop = hpv1cum_pop;
  //output.hpv2cum_pop = hpv2cum_pop;
  //output.hpv3cum_pop = hpv3cum_pop;
  //output.hpv1hncum_pop = hpv1hncum_pop;
  //output.hpv2hncum_pop = hpv2hncum_pop;
  //output.hpv3hncum_pop = hpv3hncum_pop;

  //output.hivsus_pop = hivsus_pop;
  //output.hpv1hnsus_pop = hpv1hnsus_pop;
  //output.hpv2hnsus_pop = hpv2hnsus_pop;
  //output.hpv3hnsus_pop = hpv3hnsus_pop;

  writeMatrix(f, R->hpv1cum_pop, R->dim_y_sav);
  writeMatrix(f, R->hpv2cum_pop, R->dim_y_sav);
  writeMatrix(f, R->hpv3cum_pop, R->dim_y_sav);
  writeMatrix(f, R->hpv1hncum_pop, R->dim_y_sav);
  writeMatrix(f, R->hpv2hncum_pop, R->dim_y_sav);
  writeMatrix(f, R->hpv3hncum_pop, R->dim_y_sav);

  //output.hivsus_pop = hivsus_pop;
  writeMatrix5D(f, R->hivsus_pop, R->dim_gg_1sav);

  writeMatrix(f, R->hpv1hncum, R->dim_y_sav);
  writeMatrix(f, R->hpv2hncum, R->dim_y_sav);
  writeMatrix(f, R->hpv3hncum, R->dim_y_sav);
  writeMatrix(f, R->hpv1cum, R->dim_y_sav);
  writeMatrix(f, R->hpv2cum, R->dim_y_sav);
  writeMatrix(f, R->hpv3cum, R->dim_y_sav);

  //output.ccinc = ccinc;
  //output.ccinc_neg = ccinc_neg;
  //output.ccinc_pos = ccinc_pos;
  //output.ccinc_art = ccinc_art;
  //output.ccinc_nvt = ccinc_nvt;
  //output.ccinc_nvt_neg = ccinc_nvt_neg;
  //output.ccinc_nvt_pos = ccinc_nvt_pos;
  //output.ccinc_nvt_art = ccinc_nvt_art;
  //output.cin = cin2in;
  //output.cin_neg = cin2in_neg;
  //output.cin_pos = cin2in_pos;
  //output.cin_art = cin2in_art;

  writeMatrix(f, R->ccinc, dims_y_a);
  writeMatrix(f, R->ccinc_neg, dims_y_a);
  writeMatrix(f, R->ccinc_pos, dims_y_a);
  writeMatrix(f, R->ccinc_art, dims_y_a);
  writeMatrix(f, R->ccinc_nvt, dims_y_a);
  writeMatrix(f, R->ccinc_nvt_neg, dims_y_a);
  writeMatrix(f, R->ccinc_nvt_pos, dims_y_a);
  writeMatrix(f, R->ccinc_nvt_art, dims_y_a);
  writeMatrix(f, R->cin2in, dims_y_a);
  writeMatrix(f, R->cin2in_neg, dims_y_a);
  writeMatrix(f, R->cin2in_pos, dims_y_a);
  writeMatrix(f, R->cin2in_art, dims_y_a);

  //output.hivinctot = hivinctot;

  writeMatrix4D(f, R->hivinctot, dims_y_s_a_v);

  //output.hpv1618_f = hpv1618_f;
  //output.hpv1618_m = hpv1618_m;
  //output.nvthpv_f = nvthpv_f;
  //output.nvthpv_m = nvthpv_m;

  writeMatrix3D(f, R->hpv1618_f, R->dim_g_i_a);
  writeMatrix3D(f, R->hpv1618_m, R->dim_g_i_a);
  writeMatrix3D(f, R->nvthpv_f, R->dim_g_i_a);
  writeMatrix3D(f, R->nvthpv_m, R->dim_g_i_a);

  //output.hpv9vt_f = hpv9vt_f;
  //output.hpv9vt_m = hpv9vt_m;

  writeMatrix3D(f, R->hpv9vt_f, R->dim_g_i_a);
  writeMatrix3D(f, R->hpv9vt_m, R->dim_g_i_a);

  //output.hpv_f = hpv_f;
  //output.hpv_m = hpv_m;
  //output.cinpr_f = cinpr_f;
  writeMatrix3D(f, R->hpv_f, R->dim_g_i_a);
  writeMatrix3D(f, R->hpv_m, R->dim_g_i_a);
  writeMatrix3D(f, R->cinpr_f, R->dim_g_i_a);
  writeMatrix3D(f, R->ccpr_f, R->dim_g_i_a);

  // output.screencovf = screencovf;

 // writeMatrix3D(f, R->screencovf, dims_y32);
  writeMatrix3D(f, R->screencovf_neg, dims_y32);
  writeMatrix3D(f, R->screencovf_pos, dims_y32);
  writeMatrix3D(f, R->screencovf_art, dims_y32);

  //output.hiv_fsw = hiv_fsw;
  //output.hiv_mcli = hiv_mcli;
  int _hyrid_length = 1 + ((T1 + 1) - (1985 - (int)(p->dp0[31]) + 1));

  writeArray(f, R->hiv_fsw, _hyrid_length);
  writeArray(f, R->hiv_mcli, _hyrid_length);

  //output.errmessage = errmessage;
  //output.testtime = testtime;
  //output.testtime1yr = testtime1;

  writeIntArray(f, R->errmessage, 3);
  writeArray(f, R->testtime, _hyrid_length);
  writeArray(f,R->testtime1yr, T1 + 1);

  fclose(f);

}
