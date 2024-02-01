#ifndef RESULTS_H
#define RESULTS_H

#include "matrix_8d.h"
#include "lmatrix_8d.h"

#include <cstddef>
class result;

class results {
  public:

// Results all from the function call to diffeq_vacc...

//     [popbyrisk_f, popbyrisk_m, hiv_9to14, hiv_15to49, hiv_15to74, hivart_9to14, hivart_15to49, hivart_15to74, ...

    DMatrix popbyrisk_f;
    DMatrix popbyrisk_m;
    DMatrix hiv_9to14;
    DMatrix hiv_15to49;
    DMatrix hiv_15to74;
    DMatrix hivart_9to14;
    DMatrix hivart_15to49;
    DMatrix hivart_15to74;

//    art_15plus, hiv_fage, hiv_mage, art_age, pop1980, pop2018, pop2120, hiv_fsw, hiv_mcli, ...


    DMatrix art_15plus;
    DMatrix hiv_fage;
    DMatrix hiv_mage;
    DMatrix art_age;
    dMAT_8D pop2000;
    dMAT_8D pop2018;
    dMAT_8D pop2120;
    DArray hiv_fsw;
    DArray hiv_mcli;

//    psize, psize_age_f, psize_age_f_neg, psize_age_f_pos, psize_age_f_art, ...

    DMatrix psize;
    DMatrix psize_age_f;
    DMatrix psize_age_f_neg;
    DMatrix psize_age_f_pos;
    DMatrix psize_age_f_art;

//    psize_age_m, psize_age_m_neg, psize_age_m_pos, psize_age_m_art, vacc_age_f, vacc_age_m, vacc_age_f_pos, vacc_age_m_pos, ...

    DMatrix psize_age_m;
    DMatrix psize_age_m_neg;
    DMatrix psize_age_m_pos;
    DMatrix psize_age_m_art;
    DMatrix vacc_age_f;
    DMatrix vacc_age_m;
    DMatrix vacc_age_f_pos;
    DMatrix vacc_age_f_art;
    DMatrix vacc_age_m_pos;

//       hivd_f, hivd_m, ccd_f, ccd_f_neg, ccd_f_pos, ccd_f_art, othd_f, othd_m, ...

    DMatrix hivd_f;
    DMatrix hivd_m;
    DMatrix ccd_f;
    DMatrix ccd_f_neg;
    DMatrix ccd_f_pos;
    DMatrix ccd_f_art;
    DMatrix othd_f;
    DMatrix othd_m;

//    hpv1hncum_pop, hpv2hncum_pop, hpv3hncum_pop, hpv1cum_pop, hpv2cum_pop, hpv3cum_pop, ...

    DMatrix hpv1hncum_pop;
    DMatrix hpv2hncum_pop;
    DMatrix hpv3hncum_pop;
    DMatrix hpv1cum_pop;
    DMatrix hpv2cum_pop;
    DMatrix hpv3cum_pop;

//    hivsus_pop, hivinc_pop_cum, ...

    DMatrix5D hivsus_pop;
    DMatrix4D hivinc_pop_cum;

    //ccinc_f, ccinc_f_neg, ccinc_f_pos, ccinc_f_art, ccinc_nvtf, ccinc_nvtf_neg, ccinc_nvtf_pos, ccinc_nvtf_art, cinin_f, cinin_f_neg, cinin_f_pos, cinin_f_art, ...


    DMatrix ccinc_f;
    DMatrix ccinc_f_neg;
    DMatrix ccinc_f_pos;
    DMatrix ccinc_f_art;
    DMatrix ccinc_nvtf;
    DMatrix ccinc_nvtf_neg;
    DMatrix ccinc_nvtf_pos;
    DMatrix ccinc_nvtf_art;
    DMatrix cinin_f;
    DMatrix cinin_f_neg;
    DMatrix cinin_f_pos;
    DMatrix cinin_f_art;

//    hpv1618_f, hpv1618_m, nvthpv_f, nvthpv_m, hpv9vt_f, hpv9vt_m, hpv_f, hpv_m, cinpr_f, ccpr_f, screencov_yr, ...

    DMatrix3D hpv1618_f;
    DMatrix3D hpv1618_m;
    DMatrix3D nvthpv_f;
    DMatrix3D nvthpv_m;
    DMatrix3D hpv9vt_f;
    DMatrix3D hpv9vt_m;
    DMatrix3D hpv_f;
    DMatrix3D hpv_m;
    DMatrix3D cinpr_f;
    DMatrix3D ccpr_f;
  //  DMatrix3D screencov_yr;
    DMatrix3D screencov_neg_yr;
    DMatrix3D screencov_pos_yr;
    DMatrix3D screencov_art_yr;

//    errmessage, testtime, testtime1] = diffeq_vacc(T1, tstep, pop0, dp0, bp0, hp0, pp0, vaccscen, screenscen);
//  Noting that testtime1 is renamed to testtime1yr in the actiona diffeq_vacc function header.

    int* errmessage;
    DArray testtime;
    DArray testtime1yr;

// These, additionally, are global in the original code, so should be encapsulated somewhere - here seems
// as good a place as any.

    DMatrix deatht;
    DMatrix hivdftot;
    DMatrix hivdmtot;
    DMatrix ccdftot;
    DMatrix ccdftot_neg;
    DMatrix ccdftot_pos;
    DMatrix ccdftot_art;
    DMatrix ccinc;
    DMatrix ccinc_neg;
    DMatrix ccinc_pos;
    DMatrix ccinc_art;
    DMatrix ccinc_nvt;
    DMatrix ccinc_nvt_neg;
    DMatrix ccinc_nvt_pos;
    DMatrix ccinc_nvt_art;
    //DMatrix cin;
    //DMatrix cin_neg;
    //DMatrix cin_pos;
    //DMatrix cin_art;
    DMatrix cin2in;
    DMatrix cin2in_neg;
    DMatrix cin2in_pos;
    DMatrix cin2in_art;
    DMatrix4D hivinctot;
    //DMatrix3D screencovf;
    DMatrix3D screencovf_neg;
    DMatrix3D screencovf_pos;
    DMatrix3D screencovf_art;

    DMatrix hpv1hncum;
    DMatrix hpv2hncum;
    DMatrix hpv3hncum;
    DMatrix hpv1cum;
    DMatrix hpv2cum;
    DMatrix hpv3cum;

    int* dim_g_a;
    int* dim_g_i_a;
    int* dim_g_r;
    int* dim_g_s;
    int* dim_g_r_a5;
    int* dim_gg_s;
    int* dim_gg_1sav;
    int* dim_gg_sav;
    int* dim_gg_a;
    int* dim_gg_a5;
    int* dim_y_sav;

    results();
};

#endif
