#include "state.h"
#include <cstddef>

state::state() {
  age_id = (cMAT_8D) NULL;
  age_id2 = (cMAT_8D) NULL;
  alpha = (dMAT_8D)NULL;
  art_id = (cMAT_8D) NULL;
  betahivacute = (dMAT_8D) NULL;
  betahivart = (dMAT_8D) NULL;
  betahiv2 = (dMAT_8D)NULL;
  betahivacute_ai = (dMAT_8D)NULL;
  betahiv2_ai = (dMAT_8D)NULL;
  betahivart_ai = (dMAT_8D)NULL;
  betahpv1 = (dMAT_8D) NULL;
  betahpv2 = (dMAT_8D) NULL;
  betahpv3 = (dMAT_8D) NULL;
  CC_id = (cMAT_8D) NULL;   // #new_age_structre
  CIN_id = (cMAT_8D) NULL;
  delta1 = (dMAT_8D) NULL;
  delta2 = (dMAT_8D) NULL;
  delta3 = (dMAT_8D) NULL;
  eta = (dMAT_8D)NULL;
  hiv_id = (cMAT_8D) NULL;
  hiv2_id = (cMAT_8D) NULL;
  hivS_id = (cMAT_8D) NULL;
  hpv123_id = (cMAT_8D) NULL;
  hpv1_id = (cMAT_8D) NULL;
  hpv2_id = (cMAT_8D) NULL;
  hpv3_id = (cMAT_8D) NULL;
  hpv12_id = (cMAT_8D) NULL;
  hpv1I_id = (cMAT_8D) NULL;
  hpv2I_id = (cMAT_8D) NULL;
  hpv3I_id = (cMAT_8D) NULL;
  hpv1R_id = (cMAT_8D) NULL;
  hpv2R_id = (cMAT_8D) NULL;
  hpv3R_id = (cMAT_8D) NULL;
  hpv1S_id = (cMAT_8D) NULL;
  hpv2S_id = (cMAT_8D) NULL;
  hpv3S_id = (cMAT_8D) NULL;
  hpv1CIN_id = (cMAT_8D) NULL;
  hpv2CIN_id = (cMAT_8D) NULL;
  hpv3CIN_id = (cMAT_8D) NULL;
  hpv1CC_id= (cMAT_8D)NULL;
  hpv2CC_id = (cMAT_8D)NULL;
  hpv3CC_id = (cMAT_8D)NULL;
  hups1 = (dMAT_8D) NULL;
  hups2 = (dMAT_8D) NULL;
  hups3 = (dMAT_8D) NULL;
  nacts = (dMAT_8D) NULL;
  nacts_comm = (dMAT_8D) NULL;
  mu = (dMAT_8D) NULL;
  noart_id = (cMAT_8D) NULL;
  nu = (dMAT_8D)NULL;
  omega1 = (dMAT_8D) NULL;
  omega2 = (dMAT_8D) NULL;
  omega3 = (dMAT_8D) NULL;
  omikron = (dMAT_8D) NULL;
  p_anal = (dMAT_8D) NULL;
  pcomm_anal = (dMAT_8D) NULL;
  pcr_comm = (dMAT_3D_SRA) NULL;
  pcr = (dMAT_3D_SRA) NULL;
  pi1 = (dMAT_8D)NULL;
  pi2 = (dMAT_8D)NULL;
  pi3 = (dMAT_8D)NULL;
  psi1 = (dMAT_8D) NULL;
  psi2 = (dMAT_8D) NULL;
  psi3 = (dMAT_8D) NULL;
  risk_id = (cMAT_8D) NULL;
  sexf_id = (cMAT_8D) NULL;
  sexm_id = (cMAT_8D) NULL;
  sigma1 = (dMAT_8D) NULL;
  sigma2 = (dMAT_8D) NULL;
  sigma3 = (dMAT_8D) NULL;
  tau = (dMAT_8D) NULL;
  theta1 = (dMAT_8D)NULL;
  theta2 = (dMAT_8D)NULL;
  theta3 = (dMAT_8D)NULL;
  ups1 = (dMAT_8D)NULL;
  ups2 = (dMAT_8D)NULL;
  ups3 = (dMAT_8D)NULL;
  vacc_id = (cMAT_8D)NULL;
 // vxxwane = (dMAT_8D)NULL;
}