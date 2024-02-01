#ifndef STATE_H
#define STATE_H

#include "matrix_8d.h"
#include "lmatrix_8d.h"

class state;

class state {
public:
  cMAT_8D age_id;
  cMAT_8D age_id2;
  dMAT_8D alpha;
  cMAT_8D art_id;
  dMAT_8D betahivacute;
  dMAT_8D betahivart;
  dMAT_8D betahiv2;
  dMAT_8D betahivacute_ai;
  dMAT_8D betahiv2_ai;
  dMAT_8D betahivart_ai;
  dMAT_8D betahpv1;
  dMAT_8D betahpv2;
  dMAT_8D betahpv3;
  cMAT_8D CC_id;   // #new_age_structre
  cMAT_8D CIN_id;
  dMAT_8D delta1;
  dMAT_8D delta2;
  dMAT_8D delta3;
  dMAT_8D eta;
  cMAT_8D hiv_id;
  cMAT_8D hiv2_id;
  cMAT_8D hivS_id;
  cMAT_8D hpv123_id;
  cMAT_8D hpv12_id;
  cMAT_8D hpv1_id;
  cMAT_8D hpv2_id;
  cMAT_8D hpv3_id;
  cMAT_8D hpv1CIN_id;
  cMAT_8D hpv2CIN_id;
  cMAT_8D hpv3CIN_id;
  cMAT_8D hpv1CC_id;
  cMAT_8D hpv2CC_id;
  cMAT_8D hpv3CC_id;
  cMAT_8D hpv1I_id;
  cMAT_8D hpv2I_id;
  cMAT_8D hpv3I_id;
  cMAT_8D hpv1S_id;
  cMAT_8D hpv2S_id;
  cMAT_8D hpv3S_id;
  cMAT_8D hpv1R_id;
  cMAT_8D hpv2R_id;
  cMAT_8D hpv3R_id;
  dMAT_8D hups1;
  dMAT_8D hups2;
  dMAT_8D hups3;
  dMAT_8D nacts;
  dMAT_8D nacts_comm;
  cMAT_8D noart_id;
  dMAT_8D mu;
  dMAT_8D nu;
  dMAT_8D omega1;
  dMAT_8D omega2;
  dMAT_8D omega3;
  dMAT_8D omikron;
  dMAT_8D p_anal;
  dMAT_8D pcomm_anal;
  dMAT_3D_SRA pcr_comm;
  dMAT_3D_SRA pcr;
  dMAT_8D pi1;
  dMAT_8D pi2;
  dMAT_8D pi3;
  dMAT_8D psi1;
  dMAT_8D psi2;
  dMAT_8D psi3;
  cMAT_8D risk_id;
  cMAT_8D sexf_id;
  cMAT_8D sexm_id;
  dMAT_8D sigma1;
  dMAT_8D sigma2;
  dMAT_8D sigma3;
  dMAT_8D tau;
  dMAT_8D theta1;
  dMAT_8D theta2;
  dMAT_8D theta3;
  dMAT_8D ups1;
  dMAT_8D ups2;
  dMAT_8D ups3;
  cMAT_8D vacc_id;
  dMAT_8D vxxwane;

  ~state();
  state();
};

#endif
