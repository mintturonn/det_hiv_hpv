#include "diffeq_vacc.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

results* diffeq_vacc(int tend, int tstep, dMAT_8D pop0, params* p, state* S, int* dims) {

  results* R = new results();

  //% Population structure
  //    popn = pop0;% simulation, ode45; pop_col = pop0 as defined by[t, p]
  //    % size(pop0)
  //    niterations = tend * tstep;
  //  dt = 1 / tstep;

  dMAT_8D popn = dNEW_COPY_8D(popn, pop0);
  dMAT_8D spare = dNEW_8D(spare);
  dMAT_8D popn_sum_h1 = dNEW_8D(popn_sum_h1);
  dMAT_8D popn_sum_h1_h2 = dNEW_8D(popn_sum_h1_h2);
  dMAT_8D popn_sum_h1_h2_h3 = dNEW_8D(popn_sum_h1_h2_h3);
  dMAT_8D popn_sum_h1_h2_h3_i = dNEW_8D(popn_sum_h1_h2_h3_i);
  dMAT_8D tmpsum_h1_h2 = dNEW_8D(tmpsum_h1_h2);
  dMAT_8D tmpsum_h1_h2_h3 = dNEW_8D(tmpsum_h1_h2_h3);
  dMAT_8D tmpsum_h1_h2_h3_i = dNEW_8D(tmpsum_h1_h2_h3_i);

  DArray spare_v1 = CREATE_DARRAY(_SIZE_V1, spare_v1, 0);
  dMAT_2D_AV spare_av = dNEW_2D_AV(spare_av);
  dMAT_3D_IAV spare_iav = dNEW_3D_IAV(spare_iav);
  dMAT_4D_H3IAV spare_h3iav = dNEW_4D_H3IAV(spare_h3iav);
  dMAT_4D_H2H3IR spare_h2h3ir = dNEW_4D_H2H3IR(spare_h2h3ir);
  dMAT_4D_H1H3IR spare_h1h3ir = dNEW_4D_H1H3IR(spare_h1h3ir);
  dMAT_4D_H1H2IR spare_h1h2ir = dNEW_4D_H1H2IR(spare_h1h2ir);


  R->pop2000 = dNEW_8D(R->pop2000);
  R->pop2018 = dNEW_8D(R->pop2018);
  R->pop2120 = dNEW_8D(R->pop2120);

  int niterations = (int)(tend * tstep);

  double dt = 1.0 / tstep;

  //  [foihiv, foihpv1, foihpv2, foihpv3] = deal(zeros(h1, h2, h3, i, s, r, a, v1));

  dMAT_8D foihiv = dNEW_8D(foihiv);
  dMAT_8D foihpv1 = dNEW_8D(foihpv1);
  dMAT_8D foihpv2 = dNEW_8D(foihpv2);
  dMAT_8D foihpv3 = dNEW_8D(foihpv3);

  //  [tempfoih1, tempfoih2, tempfoih3, tempfoih1_comm, tempfoih2_comm, tempfoih3_comm, ...
  //    tempmoih1, tempmoih2, tempmoih3, tempmoih1_comm, tempmoih2_comm, tempmoih3_comm, ...
  //    tempfoihpv1, tempfoihpv2, tempfoihpv3, tempmoihpv1, tempmoihpv2, tempmoihpv3] = deal(zeros(r * a, r * a));

  dMAT_2D_RA6RA6 tempfoih1 = dNEW_2D_RA6RA6(tempfoih1)
  dMAT_2D_RA6RA6 tempfoih2 = dNEW_2D_RA6RA6(tempfoih2)
  dMAT_2D_RA6RA6 tempfoih3 = dNEW_2D_RA6RA6(tempfoih3)

  // These are declared up here in the Matlab code above, but when used, have an immediately different shape...

  //DMatrix tempfoih1_comm = CREATE_DMATRIX(dims_ra, tempfoih1_comm, 0.0)
  //DMatrix tempfoih2_comm = CREATE_DMATRIX(dims_ra, tempfoih2_comm, 0.0)
  //DMatrix tempfoih3_comm = CREATE_DMATRIX(dims_ra, tempfoih3_comm, 0.0)
  //DMatrix tempmoih1_comm = CREATE_DMATRIX(dims_ra, tempmoih1_comm, 0.0)
  //DMatrix tempmoih2_comm = CREATE_DMATRIX(dims_ra, tempmoih2_comm, 0.0)
  //DMatrix tempmoih3_comm = CREATE_DMATRIX(dims_ra, tempmoih3_comm, 0.0)

  dMAT_2D_RA6RA6 tempmoih1 = dNEW_2D_RA6RA6(tempmoih1);
  dMAT_2D_RA6RA6 tempmoih2 = dNEW_2D_RA6RA6(tempmoih2);
  dMAT_2D_RA6RA6 tempmoih3 = dNEW_2D_RA6RA6(tempmoih3);
  dMAT_2D_RA6RA6 tempfoihpv1 = dNEW_2D_RA6RA6(tempfoihpv1);
  dMAT_2D_RA6RA6 tempfoihpv2 = dNEW_2D_RA6RA6(tempfoihpv2);
  dMAT_2D_RA6RA6 tempfoihpv3 = dNEW_2D_RA6RA6(tempfoihpv3);
  dMAT_2D_RA6RA6 tempmoihpv1 = dNEW_2D_RA6RA6(tempmoihpv1);
  dMAT_2D_RA6RA6 tempmoihpv2 = dNEW_2D_RA6RA6(tempmoihpv2);
  dMAT_2D_RA6RA6 tempmoihpv3 = dNEW_2D_RA6RA6(tempmoihpv3);

  dMAT_2D_RARA6 tempfoih1_full = dNEW_2D_RARA6(tempfoih1_full);
  dMAT_2D_RARA6 tempfoih2_full = dNEW_2D_RARA6(tempfoih2_full);
  dMAT_2D_RARA6 tempfoih3_full = dNEW_2D_RARA6(tempfoih3_full);
  dMAT_2D_RARA6 tempmoih1_full = dNEW_2D_RARA6(tempmoih1_full);
  dMAT_2D_RARA6 tempmoih2_full = dNEW_2D_RARA6(tempmoih2_full);
  dMAT_2D_RARA6 tempmoih3_full = dNEW_2D_RARA6(tempmoih3_full);
  dMAT_2D_RARA6 tempfoihpv1_full = dNEW_2D_RARA6(tempfoihpv1_full);
  dMAT_2D_RARA6 tempfoihpv2_full = dNEW_2D_RARA6(tempfoihpv2_full);
  dMAT_2D_RARA6 tempfoihpv3_full = dNEW_2D_RARA6(tempfoihpv3_full);
  dMAT_2D_RARA6 tempmoihpv1_full = dNEW_2D_RARA6(tempmoihpv1_full);
  dMAT_2D_RARA6 tempmoihpv2_full = dNEW_2D_RARA6(tempmoihpv2_full);
  dMAT_2D_RARA6 tempmoihpv3_full = dNEW_2D_RARA6(tempmoihpv3_full);


  //  errmessage = zeros(1, 3);
  //  [popbyhiv_agef, popbyhiv_agem, hiv_agefsw, hiv_ageclient] = deal(zeros(niterations, i * a));

  R->errmessage = new int[4]{ 0, 0, 0, 0 };

  //[prevsall] = deal(zeros(niterations, s * r * a));

  /***** How is deal(zeros(x,y)) different from zeros(x,y) ? *****/

  //int* dims_n_sra = new int[3]{ 0, niterations, dims[DIM_s] * dims[DIM_r] * dims[DIM_a] };
  //DMatrix3D prevsall = CREATE_DMATRIX(dims_n_sra, prevsall, 0.0)  -  never used.

  //delete[] dims_n_sra;

  //  [hivd_ft, hivd_mt, othd_ft, othd_mt, ...
  //    ccinc_ft, ccinc_ft_art, ccinc_ft_pos, ccinc_ft_neg, ...
  //    ccinc_nvtft, ccinc_nvtft_art, ccinc_nvtft_pos, ccinc_nvtft_neg, ...
  //    cin_ft, cin_ft_art, cin_ft_pos, cin_ft_neg, ...
  //    ccd_ft, ccd_ft_art, ccd_ft_pos, ccd_ft_neg] = deal(zeros(niterations, a));

  dMAT_2D_NA hivd_ft = dNEW_2D_NA(hivd_ft);
  dMAT_2D_NA hivd_mt = dNEW_2D_NA(hivd_mt);
  dMAT_2D_NA othd_ft = dNEW_2D_NA(othd_ft);
  dMAT_2D_NA othd_mt = dNEW_2D_NA(othd_mt);
  dMAT_2D_NA ccinc_ft = dNEW_2D_NA(ccinc_ft);
  dMAT_2D_NA ccinc_ft_art = dNEW_2D_NA(ccinc_ft_art);
  dMAT_2D_NA ccinc_ft_pos = dNEW_2D_NA(ccinc_ft_pos);
  dMAT_2D_NA ccinc_ft_neg = dNEW_2D_NA(ccinc_ft_neg);
  dMAT_2D_NA ccinc_nvtft = dNEW_2D_NA(ccinc_nvtft);
  dMAT_2D_NA ccinc_nvtft_art = dNEW_2D_NA(ccinc_nvtft_art);
  dMAT_2D_NA ccinc_nvtft_pos = dNEW_2D_NA(ccinc_nvtft_pos);
  dMAT_2D_NA ccinc_nvtft_neg = dNEW_2D_NA(ccinc_nvtft_neg);
  dMAT_2D_NA cin_ft = dNEW_2D_NA(cin_ft);
  dMAT_2D_NA cin_ft_art = dNEW_2D_NA(cin_ft_art);
  dMAT_2D_NA cin_ft_pos = dNEW_2D_NA(cin_ft_pos);
  dMAT_2D_NA cin_ft_neg = dNEW_2D_NA(cin_ft_neg);
  dMAT_2D_NA ccd_ft = dNEW_2D_NA(ccd_ft);
  dMAT_2D_NA ccd_ft_art = dNEW_2D_NA(ccd_ft_art);
  dMAT_2D_NA ccd_ft_pos = dNEW_2D_NA(ccd_ft_pos);
  dMAT_2D_NA ccd_ft_neg = dNEW_2D_NA(ccd_ft_neg);

  //  [hivinc_pop_cumt] = deal(zeros(niterations, s, a, v1));
  //  [screencov] = deal(zeros(niterations, 3, 2));

  dMAT_4D_NSAV1 hivinc_pop_cumt = dNEW_4D_NSAV1(hivinc_pop_cumt);
  dMAT_4D_NSAV1 hpv1cum_popt = dNEW_4D_NSAV1(hpv1cum_popt);
  dMAT_4D_NSAV1 hpv2cum_popt = dNEW_4D_NSAV1(hpv2cum_popt);
  dMAT_4D_NSAV1 hpv3cum_popt = dNEW_4D_NSAV1(hpv3cum_popt);
  dMAT_4D_NSAV1 hpv1hncum_popt = dNEW_4D_NSAV1(hpv1hncum_popt);
  dMAT_4D_NSAV1 hpv2hncum_popt = dNEW_4D_NSAV1(hpv2hncum_popt);
  dMAT_4D_NSAV1 hpv3hncum_popt = dNEW_4D_NSAV1(hpv3hncum_popt);

 // dMAT_3D_N32 screencov = dNEW_3D_N32(screencov);
  dMAT_3D_N32 screencov_art = dNEW_3D_N32(screencov_art);
  dMAT_3D_N32 screencov_pos = dNEW_3D_N32(screencov_pos);
  dMAT_3D_N32 screencov_neg = dNEW_3D_N32(screencov_neg);

  //  [agein1, agein2, agein3, conage] = deal(zeros(h1, h2, h3, i, s, r, a, v1));
  //  [ , ageinm] = deal(zeros(3, 1));

    //agein gets defined lower down in 7 dimensions.

  dMAT_8D agein1 = dNEW_8D(agein1);
  dMAT_8D agein2 = dNEW_8D(agein2);
  dMAT_8D agein3 = dNEW_8D(agein3);
  dMAT_8D agein4 = dNEW_8D(agein4);
  dMAT_8D agein5 = dNEW_8D(agein5);
  dMAT_8D agein6 = dNEW_8D(agein6);
  dMAT_8D conage = dNEW_8D(conage);
  dMAT_8D betahpv1pp_comm = dNEW_8D(betahpv1pp_comm);
  dMAT_8D betahpv2pp_comm = dNEW_8D(betahpv2pp_comm);
  dMAT_8D betahpv3pp_comm = dNEW_8D(betahpv3pp_comm);

  //  %% Time - varying parameters
  //    % Pop growth : how many new people enter the model per year
  //    gr_yr = T0; % year when growth in popsize starts
  //    gr_start = 0; %
  //    gr_min = dp0(15);
  //    gr_max = dp0(14);
  //    gr_grate = dp0(16);
  //    gr_plat = (2040 - T0) * tstep;% year after which forced to be stable
  //    gr_tv = tvarying_plat(niterations, gr_start, gr_min, gr_max, gr_grate, gr_plat);

  // gr_yr is never used

#define T0 p->dp0[31]

#define gr_start 0
#define gr_min p->dp0[15]
#define gr_max p->dp0[14]
#define gr_grate p->dp0[16]
  const int gr_plat = (const int)round((2040 - T0) * tstep);
  double* gr_tv = tvarying_plat(niterations, gr_start, gr_min, gr_max, gr_grate, gr_plat);

  //  % population size to remain constant
  //    % gr_tv = ones(niterations, 1);

  //  % exponential population growth
  //    % gr_tv = ones(niterations, 1);
  //  % gr_tv = gr_tv + dp0(15) - 1;

  //  % Perinatal HIV infection
  //  a1 = round((hp0(30) - T0) * tstep, 0);
  //  a2 = round((hp0(31) - T0) * tstep, 0);
  //  a3 = round((hp0(32) - T0) * tstep, 0);
  //  a4 = round((hp0(33) - T0) * tstep, 0);
  //  x0 = hp0(25);
  //  x1 = hp0(26);
  //  x2 = hp0(27);
  //  x3 = x2 * hp0(28);
  //  x4 = x3 * hp0(29);

  int a1 = (int)round((p->hp0[30] - T0) * tstep);
  int a2 = (int)round((p->hp0[31] - T0) * tstep);
  int a3 = (int)round((p->hp0[32] - T0) * tstep);
  int a4 = (int)round((p->hp0[33] - T0) * tstep);
#define x0 p->hp0[25]
#define x1 p->hp0[26]
#define x2 p->hp0[27]
  double x3 = x2 * p->hp0[28];
  double x4 = x3 * p->hp0[29];

  //  [perin, perin2] = hivperinatal(niterations, a1, a2, a3, a4, x0, x1, x2, x3, x4);

  // This is a single-use function, and easier I think to just inline it here.

  //    function[perin, perin2] = hivperinatal(niterations, a1, a2, a3, a4, x0, x1, x2, x3, x4)

  //    perin = zeros(niterations, 1);
  //    perin(1:a1) = x0;
  //    perin(a1:a2) = x1;

  //    perin2 = zeros(niterations, 1);
  //    perin2(a2 + 1:a3) = x2;
  //    perin2(a3 + 1:a4) = x3;
  //    perin2(a4 + 1:end) = x4;

  //    end

  DArray perin = CREATE_DARRAY_START(niterations, perin) CREATE_DARRAY_END
  for (int i = 1; i < a1; i++) perin[i] = x0;
  for (int i = a1; i <= a2; i++) perin[i] = x1;
  for (int i = a2 + 1; i <= niterations; i++) perin[i] = 0;

  DArray perin2 = CREATE_DARRAY_START(niterations, perin2) CREATE_DARRAY_END
  for (int i = 1; i <= a2; i++) perin2[i] = 0;
  for (int i = a2 + 1; i <= a3; i++) perin2[i] = x2;
  for (int i = a3 + 1; i <= a4; i++) perin2[i] = x3;
  for (int i = a4 + 1; i <= niterations; i++) perin2[i] = x4;

  //   % if you want to turn off ART completely, you need to do this too
  //   % perin2(:) = 0;

  //   % if perinatal HIV needs to be turned off :
  //   % [perin, perin2] = deal(zeros(niterations, 1));

  //   % VMC: the proportion of MEN who are VMC
  //     mcirc_yr = round(hp0(20), 0); % year when growth in VMC starts
  //     mcirc_yr2 = round(hp0(21), 0);
  //     mcirc_start = (mcirc_yr - T0) * tstep; % year when growth in VMC starts
  //     mcirc_start2 = (mcirc_yr2 - T0) * tstep; %
  //     mcirc_min = hp0(22);% implement in betahiv like this 0.08 * 0.6 + (1 - 0.08)
  //     mcirc_max = hp0(23);
  //     mcirc_tv = tvarying_artcon(niterations, mcirc_start, mcirc_start2, mcirc_min, mcirc_max);

  int mcirc_yr = (int)round(p->hp0[20]);
  int mcirc_yr2 = (int)round(p->hp0[21]);
  int mcirc_start = (int)round((mcirc_yr - T0) * tstep);
  int mcirc_start2 = (int)round((mcirc_yr2 - T0) * tstep);
#define mcirc_min p->hp0[22]
#define mcirc_max p->hp0[23]
  DArray mcirc_tv = tvarying_artcon(niterations, mcirc_start, mcirc_start2, mcirc_min, mcirc_max);


  //    % VMC as an average(efficacy of VMC 40 %), applied to betahiv as a
  //    % multiplier
  //    mcirc_effic = hp0(8);
  //    mcirc = mcirc_tv .* (1 - mcirc_effic) + (1 - mcirc_tv);

#define mcirc_effic p->hp0[8]
  DArray mcirc = CREATE_DARRAY_START(niterations, mcirc_tv) CREATE_DARRAY_END
  for (int i = 1; i <= niterations; i++)
    mcirc[i] = mcirc_tv[i] * (1 - mcirc_effic) + (1 - mcirc_tv[i]);

  //  % Condom use, proportion of acts protected - non - commercial
  //  conuse_yr = round(bp0(53), 0);
  //  conuse_yr2 = round(bp0(54), 0);
  //  conuse_start = (conuse_yr - T0) * tstep; %
  //  conuse_start2 = (conuse_yr2 - T0) * tstep; %
  //  conuse_min = bp0(55);%
  //  conuse_max = bp0(56);
  //  conuse_tv = tvarying_artcon(niterations, conuse_start, conuse_start2, conuse_min, conuse_max);

#define conuse_yr ((int) round(p->bp0[53]))
#define conuse_yr2 ((int) round(p->bp0[54]))
  int conuse_start = (int)round((conuse_yr - T0) * tstep);
  int conuse_start2 = (int)round((conuse_yr2 - T0) * tstep);
#define conuse_min p->bp0[55]
#define conuse_max p->bp0[56]
  DArray conuse_tv = tvarying_artcon(niterations, conuse_start, conuse_start2, conuse_min, conuse_max);

  //  % Condom use, proportion of acts protected - commercial
  //    conuse_commmin = bp0(57);
  //    conuse_commmax = bp0(58);
  //    conuse_commtv = tvarying_artcon(niterations, conuse_start, conuse_start2, conuse_commmin, conuse_commmax);

#define conuse_commmin p->bp0[57]
#define conuse_commmax p->bp0[58]
  DArray conuse_commtv = tvarying_artcon(niterations, conuse_start, conuse_start2, conuse_commmin, conuse_commmax);


  //  conage(age_id == 1) = bp0(73) * bp0(74);
  //  conage(age_id == 2) = bp0(74);
  //  conage(age_id == 3) = bp0(75) * bp0(74);
  //  conage(age_id == 4) = bp0(76) * bp0(74);
  //  conage(age_id == 5) = bp0(77) * bp0(74);

  LOOP_ALL {
    int _age = (int)round(GET_ALL(S->age_id));
    if (_age == 2) GET_ALL(conage) = p->bp0[74];
    else GET_ALL(conage) = p->bp0[72 + _age] * p->bp0[74];
  }


    //  % ART
    //    tau_yr = round(hp0(11), 0); % year when growth ART initiaion starts
    //    tau_yr2 = round(hp0(12), 0); % year when growth ART initiaion starts
    //    tau_start = (tau_yr - T0) * tstep; % niteration in ART initiation starts
    //    tau_start2 = (tau_yr2 - T0) * tstep; %
    //    tau_max = hp0(13);
    //    tau_tv = tvarying_artcon(niterations, tau_start, tau_start2, 0, tau_max);

#define tau_yr ((int) round(p->hp0[11]))
#define tau_yr2 ((int) round(p->hp0[12]))
  int tau_start = (int)round((tau_yr - T0) * tstep);
  int tau_start2 = (int)round((tau_yr2 - T0) * tstep);
#define tau_max p->hp0[13]
  DArray tau_tv = tvarying_artcon(niterations, tau_start, tau_start2, 0.0, tau_max);

  //  % Relative rates based on sexand risk
  //    tau_coff = zeros(niterations, 1);
  //  tau_coff_yr = (pp0(46) - T0) * tstep;
  //  tau_coff(1:tau_coff_yr) = pp0(6) * pp0(47);
  //  tau_coff(tau_coff_yr:end) = pp0(6);

  DArray tau_coff = CREATE_DARRAY(niterations, tau_coff, p->pp0[6] * p->pp0[47]);
  int tau_coff_yr = (int)round((p->pp0[46] - T0) * tstep);
  for (int i = tau_coff_yr; i <= niterations; i++) tau_coff[i] = p->pp0[6];

  //  % Cervical cancer screening
  //  tyr1 = (pp0(68) - T0) * tstep;
  //  tyr2 = (pp0(69) - T0) * tstep;
  //  tyr3 = (pp0(70) - T0) * tstep;
  //  tyr4 = (pp0(71) - T0) * tstep;
  //  tyr5 = (pp0(72) - T0) * tstep;

  int tyr1 = (int)round((p->pp0[68] - T0) * tstep);
  int tyr2 = (int)round((p->pp0[69] - T0) * tstep);
  int tyr3 = (int)round((p->pp0[70] - T0) * tstep);
  int tyr4 = (int)round((p->pp0[71] - T0) * tstep);
  int tyr5 = (int)round((p->pp0[72] - T0) * tstep);

  //  %% the following pertains to CC screening scenarios in the interventions
  //    if screenscen.setup == 0

  //      theta15to24hn = tvarying_theta(niterations, tyr1, tyr2, tyr3, tyr4, tyr5, ...
  //        pp0(49), pp0(50), pp0(51), pp0(52), pp0(53));

  //      theta15to24hp = theta15to24hn;

  //      theta25to34hn = tvarying_theta(niterations, tyr1, tyr2, tyr3, tyr4, tyr5, ...
  //        pp0(54), pp0(55), pp0(56), pp0(57), pp0(58));

  //      theta25to34hp = theta25to34hn;

  //      theta35to49hn = tvarying_theta(niterations, tyr1, tyr2, tyr3, tyr4, tyr5, ...
  //        pp0(59), pp0(60), pp0(61), pp0(62), pp0(63));

  //      theta35to49hp = theta35to49hn;
  //  end

  DArray theta25to29hn = (DArray) NULL;
  DArray theta25to29hp = (DArray) CREATE_DARRAY(niterations, theta25to29hp, 0.0)
 // DArray theta30to34hn = (DArray) CREATE_DARRAY(niterations, theta30to34hn, 0.0)
  DArray theta30to34hp = (DArray) CREATE_DARRAY(niterations, theta30to34hp, 0.0)
  DArray theta35to39hn = (DArray) NULL;
  DArray theta35to39hp = (DArray) CREATE_DARRAY(niterations, theta35to39hp, 0.0)
 // DArray theta40to44hn = (DArray) CREATE_DARRAY(niterations, theta40to44hn, 0.0)
  DArray theta40to44hp = (DArray) CREATE_DARRAY(niterations, theta40to44hp, 0.0)
  DArray theta45to49hn = (DArray) NULL;
  DArray theta45to49hp = (DArray) CREATE_DARRAY(niterations, theta45to49hp, 0.0)
// this is baseline
  if ((int)round(p->screenscen_setup) == 0) {
    theta25to29hn = tvarying_theta(niterations, tyr1, tyr2, tyr3, tyr4, tyr5,
    p->pp0[49], p->pp0[50], p->pp0[51], p->pp0[52], p->pp0[53]);
    theta35to39hn = tvarying_theta(niterations, tyr1, tyr2, tyr3, tyr4, tyr5,
    p->pp0[54], p->pp0[55], p->pp0[56], p->pp0[57], p->pp0[58]);
    theta45to49hn = tvarying_theta(niterations, tyr1, tyr2, tyr3, tyr4, tyr5,
    p->pp0[59], p->pp0[60], p->pp0[61], p->pp0[62], p->pp0[63]);
    COPY_DARRAY(theta25to29hn, theta25to29hp, niterations)
    COPY_DARRAY(theta35to39hn, theta35to39hp, niterations)
    COPY_DARRAY(theta45to49hn, theta45to49hp, niterations)
  }

  //  if screenscen.setup == 1

  if ((int)round(p->screenscen_setup) == 1) {

    //    tyr6 = (screenscen.yr1 - T0) * tstep;
    //    tyr7 = (screenscen.yr2 - T0) * tstep;
    //    tyr8 = (screenscen.yr3 - T0) * tstep;

    int tyr6 = (int)round((p->screenscen_yr1 - T0) * tstep);
    int tyr7 = (int)round((p->screenscen_yr2 - T0) * tstep);
    int tyr8 = (int)round((p->screenscen_yr3 - T0) * tstep);

    //  theta15to24hn = tvarying_theta_add(niterations, tyr1, tyr2, tyr3, tyr4, tyr5, tyr6, tyr7, tyr8, ...
    //    pp0(49), pp0(50), pp0(51), pp0(52), pp0(53), screenscen.scen_15to24(1), screenscen.scen_15to24(2), screenscen.scen_15to24(3));

    theta25to29hn = tvarying_theta_add(niterations, tyr1, tyr2, tyr3, tyr4, tyr5, tyr6, tyr7, tyr8,
      p->pp0[49], p->pp0[50], p->pp0[51], p->pp0[52], p->pp0[53], p->screenscen_scen_15to24[1], p->screenscen_scen_15to24[2], p->screenscen_scen_15to24[3]);

    //  theta15to24hp = theta15to24hn;

    COPY_DARRAY(theta25to29hn, theta25to29hp, niterations)


      //  theta25to34hn = tvarying_theta_add(niterations, tyr1, tyr2, tyr3, tyr4, tyr5, tyr6, tyr7, tyr8, ...
      //    pp0(54), pp0(55), pp0(56), pp0(57), pp0(58), screenscen.scen_25to34(1), screenscen.scen_25to34(2), screenscen.scen_25to34(3));

    theta35to39hn = tvarying_theta_add(niterations, tyr1, tyr2, tyr3, tyr4, tyr5, tyr6, tyr7, tyr8,
      p->pp0[54], p->pp0[55], p->pp0[56], p->pp0[57], p->pp0[58], p->screenscen_scen_25to34[1], p->screenscen_scen_25to34[2], p->screenscen_scen_25to34[3]);

    //  theta25to34hp = theta25to34hn;

    COPY_DARRAY(theta35to39hn, theta35to39hp, niterations)

      //  theta35to49hn = tvarying_theta_add(niterations, tyr1, tyr2, tyr3, tyr4, tyr5, tyr6, tyr7, tyr8, ...
      //    pp0(59), pp0(60), pp0(61), pp0(62), pp0(62), screenscen.scen_35to49(1), screenscen.scen_35to49(2), screenscen.scen_35to49(3));

    theta45to49hn = tvarying_theta_add(niterations, tyr1, tyr2, tyr3, tyr4, tyr5, tyr6, tyr7, tyr8,
      p->pp0[59], p->pp0[60], p->pp0[61], p->pp0[62], p->pp0[63], p->screenscen_scen_35to49[1], p->screenscen_scen_35to49[2], p->screenscen_scen_35to49[3]);

    //  theta35to49hp = theta35to49hn;
    COPY_DARRAY(theta45to49hn, theta45to49hp, niterations)
  }
  //  end

  //  if screenscen.setup == 2

  if ((int)round(p->screenscen_setup) == 2) {

    //  tyr6 = (screenscen.yr1 - T0) * tstep;
    //  tyr7 = (screenscen.yr2 - T0) * tstep;
    //  tyr8 = (screenscen.yr3 - T0) * tstep;

    int tyr6 = (int)round((p->screenscen_yr1 - T0) * tstep);
    int tyr7 = (int)round((p->screenscen_yr2 - T0) * tstep);
    int tyr8 = (int)round((p->screenscen_yr3 - T0) * tstep);

    //  theta15to24hn = tvarying_theta_add(niterations, tyr1, tyr2, tyr3, tyr4, tyr5, tyr6, tyr7, tyr8, ...
    //    pp0(49), pp0(50), pp0(51), pp0(52), pp0(53), screenscen.scen_15to24hn(1), screenscen.scen_15to24hn(2), screenscen.scen_15to24hn(3));

    theta25to29hn = tvarying_theta_add(niterations, tyr1, tyr2, tyr3, tyr4, tyr5, tyr6, tyr7, tyr8,
      p->pp0[49], p->pp0[50], p->pp0[51], p->pp0[52], p->pp0[53], p->screenscen_scen_15to24hn[1], p->screenscen_scen_15to24hn[2], p->screenscen_scen_15to24hn[3]);

    //  theta15to24hp = tvarying_theta_add(niterations, tyr1, tyr2, tyr3, tyr4, tyr5, tyr6, tyr7, tyr8, ...
    //    pp0(49), pp0(50), pp0(51), pp0(52), pp0(53), screenscen.scen_15to24hp(1), screenscen.scen_15to24hp(2), screenscen.scen_15to24hp(3));

    theta25to29hp = tvarying_theta_add(niterations, tyr1, tyr2, tyr3, tyr4, tyr5, tyr6, tyr7, tyr8,
      p->pp0[49], p->pp0[50], p->pp0[51], p->pp0[52], p->pp0[53], p->screenscen_scen_15to24hp[1], p->screenscen_scen_15to24hp[2], p->screenscen_scen_15to24hp[3]);

    //  theta25to34hn = tvarying_theta_add(niterations, tyr1, tyr2, tyr3, tyr4, tyr5, tyr6, tyr7, tyr8, ...
    //    pp0(54), pp0(55), pp0(56), pp0(57), pp0(58), screenscen.scen_25to34hn(1), screenscen.scen_25to34hn(2), screenscen.scen_25to34hn(3));

    theta35to39hn = tvarying_theta_add(niterations, tyr1, tyr2, tyr3, tyr4, tyr5, tyr6, tyr7, tyr8,
      p->pp0[54], p->pp0[55], p->pp0[56], p->pp0[57], p->pp0[58], p->screenscen_scen_25to34hn[1], p->screenscen_scen_25to34hn[2], p->screenscen_scen_25to34hn[3]);

    //  theta25to34hp = tvarying_theta_add(niterations, tyr1, tyr2, tyr3, tyr4, tyr5, tyr6, tyr7, tyr8, ...
    //    pp0(54), pp0(55), pp0(56), pp0(57), pp0(58), screenscen.scen_25to34hp(1), screenscen.scen_25to34hp(2), screenscen.scen_25to34hp(3));

    theta35to39hp = tvarying_theta_add(niterations, tyr1, tyr2, tyr3, tyr4, tyr5, tyr6, tyr7, tyr8,
      p->pp0[54], p->pp0[55], p->pp0[56], p->pp0[57], p->pp0[58], p->screenscen_scen_25to34hp[1], p->screenscen_scen_25to34hp[2], p->screenscen_scen_25to34hp[3]);

    COPY_DARRAY(theta35to39hp, theta30to34hp, niterations)
    //  theta35to49hn = tvarying_theta_add(niterations, tyr1, tyr2, tyr3, tyr4, tyr5, tyr6, tyr7, tyr8, ...
    //    pp0(59), pp0(60), pp0(61), pp0(62), pp0(62), screenscen.scen_35to49hn(1), screenscen.scen_35to49hn(2), screenscen.scen_35to49hn(3));

    theta45to49hn = tvarying_theta_add(niterations, tyr1, tyr2, tyr3, tyr4, tyr5, tyr6, tyr7, tyr8,
      p->pp0[59], p->pp0[60], p->pp0[61], p->pp0[62], p->pp0[63], p->screenscen_scen_35to49hn[1], p->screenscen_scen_35to49hn[2], p->screenscen_scen_35to49hn[3]);

    //  theta35to49hp = tvarying_theta_add(niterations, tyr1, tyr2, tyr3, tyr4, tyr5, tyr6, tyr7, tyr8, ...
    //    pp0(59), pp0(60), pp0(61), pp0(62), pp0(62), screenscen.scen_35to49hp(1), screenscen_scen_35to49hp(2), screenscen.scen_35to49hp(3));

    theta45to49hp = tvarying_theta_add(niterations, tyr1, tyr2, tyr3, tyr4, tyr5, tyr6, tyr7, tyr8,
      p->pp0[59], p->pp0[60], p->pp0[61], p->pp0[62], p->pp0[63], p->screenscen_scen_35to49hp[1], p->screenscen_scen_35to49hp[2], p->screenscen_scen_35to49hp[3]);

    COPY_DARRAY(theta45to49hp, theta40to44hp, niterations)
  }
  //  end

  /**************** END ***********/

  //  %% Mixing by age

  //  f_am = [bp0(61) bp0(62) bp0(63) bp0(64) bp0(65)];
  //  f_amx = bp0(66);
  //  m_am = [bp0(67) bp0(68) bp0(69) bp0(70) bp0(71)];
  //  m_amx = bp0(72);

  DArray f_am = new double[6]{ 0, p->bp0[61], p->bp0[62], p->bp0[63], p->bp0[64], p->bp0[65] };
#define f_amx p->bp0[66]
  DArray m_am = new double[6]{ 0, p->bp0[67], p->bp0[68], p->bp0[69], p->bp0[70], p->bp0[71] };
#define m_amx p->bp0[72]

  //  eaf = [f_am(1) - 3 * f_amx   1 - f_am(1)           f_amx             f_amx          f_amx              0; ...
  //         0                     f_am(2) - 2 * f_amx   1 - f_am(2)       f_amx          f_amx/2            f_amx/2; ...
  //         0                     0                     f_am(3) - f_amx   1 - f_am(3)    f_amx/2            f_amx/2; ...
  //         0                     0                     0                 f_am(4)        (1 - f_am(4))/2    (1-f_am(4))/2; ...
  //         0                     0                     0                 0              f_am(5)            1-f_am(5)];
  //         0                     0                     0                 0              0                  0 ] ;


  DMatrix eaf = new DArray[7];
  eaf[1] = new double[7]{ 0, f_am[1] - 3 * f_amx, 1 - f_am[1],         f_amx,            f_amx,       f_amx,              0};
  eaf[2] = new double[7]{ 0, 0,                   f_am[2] - 2 * f_amx, 1 - f_am[2],      f_amx,       f_amx / 2.0,        f_amx / 2.0};
  eaf[3] = new double[7]{ 0, 0,                   0,                   f_am[3] - f_amx,  1 - f_am[3], f_amx / 2.0,        f_amx  / 2.0};
  eaf[4] = new double[7]{ 0, 0,                   0,                   0,                f_am[4],     (1 - f_am[4])/2.0,  (1 - f_am[4]) / 2.0};
  eaf[5] = new double[7]{ 0, 0,                   0,                   0,                0,           f_am[5],            (1 - f_am[5])};
  eaf[6] = new double[7]{ 0, 0,                   0,                   0,                0,           0,                  1};

// eam = [m_am(1)       0             0                0                    0                    0; ...
//        1 - m_am(2)   m_am(2)       0                0                    0                    0; ...
//        m_amx         1 - m_am(3)   m_am(3) - m_amx  0                    0                    0; ...
//        m_amx         m_amx         1 - m_am(4)      m_am(4) - 2 * m_amx  0                    0; ...
//        m_amx         m_amx         m_amx            1 - m_am(5)          m_am(5) - 3 * m_amx  0; ...
//        0             m_amx         m_amx            m_amx                1 - m_am(5)          m_am(5) - 3 * m_amx];

  DMatrix eam = new DArray[7];
  eam[1] = new double[7]{ 0,    m_am[1],     0,           0,               0,                   0,                         0};
  eam[2] = new double[7]{ 0,    1 - m_am[2], m_am[2],     0,               0,                   0,                         0};
  eam[3] = new double[7]{ 0,    m_amx,       1 - m_am[3], m_am[3] - m_amx, 0,                   0,                         0};
  eam[4] = new double[7]{ 0,    m_amx,       m_amx,       1 - m_am[4],     m_am[4] - 2 * m_amx, 0,                         0};
  eam[5] = new double[7]{ 0,    m_amx,       m_amx,       m_amx,           1 - m_am[5],         m_am[5] - 3 * m_amx,       0};
  eam[6] = new double[7]{ 0,    0,           m_amx,       m_amx,           m_amx,               1 - m_am[5],               m_am[5] - 3 * m_amx};
  //  % identity matrix for risk
  //    ed = [1 0 0; 0 1 0; 0 0.5 0.5];
  //    e = bp0(78);

  DMatrix ed = new DArray[4];
  ed[1] = new double[4]{ 0, 1, 0,   0 };
  ed[2] = new double[4]{ 0, 0, 1,   0 };
  ed[3] = new double[4]{ 0, 0, 0.5, 0.5 };

  double e = p->bp0[78];

  /**** Below is quite manky - double-check this later ***/

//  %% time indexing
//    yrid = tstep:tstep:(tend * tstep);
//    yrid = [1 yrid];
//    hyrid = yrid((1985 - dp0(22) + 1) : end);% when HIV introduced
//    h5yrid = yrid((1985 - dp0(22) + 1) : 5 : end);

  DArray yrid = CREATE_DARRAY_START(tend + 1, yrid) CREATE_DARRAY_END
  yrid[1] = 1;
  yrid[2] = tstep;
  for (int i = 3; i <= tend + 1; i++) yrid[i] = yrid[i - 1] + tstep;

  //   int _hyrid_length = 1 + ((tend + 1) - (1985 - (int)(p->dp0[31]) + 1));
  int _hyrid_length = 1 + ((tend + 1) - ((int)(p->hp0[14]) - (int)(p->dp0[31]) + 1));
  DArray hyrid = CREATE_DARRAY_START(_hyrid_length, hyrid) CREATE_DARRAY_END
  // for (int i = 0; i < _hyrid_length; i++) hyrid[i + 1] = yrid[(1985 - (int)(p->dp0[31]) + 1) + i];
  for (int i = 0; i < _hyrid_length; i++) hyrid[i + 1] = yrid[((int)(p->hp0[14]) - (int)(p->dp0[31]) + 1) + i];

  int _h5yrid_length = (int)round(1 + ((tend * tstep) - hyrid[1]) / (tstep * 5));
  DArray h5yrid = CREATE_DARRAY_START(_h5yrid_length, h5yrid) CREATE_DARRAY_END
  h5yrid[1] = hyrid[1];
  for (int i = 2; i <= _h5yrid_length; i++) h5yrid[i] = h5yrid[i - 1] + (tstep * 5);

  //  % initiating counters for outputs
  //  g = 1;
  //  gg = 1;
  //  ggg = 1;

  int g = 1;
  int gg = 1;
  int ggg = 1;

  //  %% vaccination initial
  //     vxx = deal(zeros(h1, h2, h3, i, s, r, a, v1));

  dMAT_8D vxx = dNEW_8D(vxx);
//  dMAT_8D vxxwane = dNEW_8D(vxxwane);

  /*************************************************************/
  /* MEMORY MANAGEMENT HERE                                    */
  /*
  /* Here we want to put all the memory creation stuff that    */
  /* happens in the main loop. Don't want to create it/destroy */
  /* it in every iteration                                     */
  /*************************************************************/

  dMAT_3D_SRA6 popbysexriskage = dNEW_3D_SRA6(popbysexriskage);
  dMAT_3D_SRA6 pships_com = dNEW_3D_SRA6(pships_com);
  dMAT_3D_SRA6 pcr_comm_bal = dNEW_3D_SRA6(pcr_comm_bal);
  dMAT_3D_SRA6 pships_com_bal = dNEW_3D_SRA6(pships_com_bal);
  dMAT_3D_SRA6 pshipsprop_age_com = dNEW_3D_SRA6(pshipsprop_age_com);
  dMAT_3D_SRA6 pships_main = dNEW_3D_SRA6(pships_main);
  dMAT_3D_SRA6 pships_ageriskprop = dNEW_3D_SRA6(pships_ageriskprop);

  dMAT_4D_ISRA6 popbyhiv = dNEW_4D_ISRA6(popbyhiv);
  dMAT_4D_H3SRA6 popbyhpv3 = dNEW_4D_H3SRA6(popbyhpv3);
  dMAT_4D_H2SRA6 popbyhpv2 = dNEW_4D_H2SRA6(popbyhpv2);
  dMAT_4D_H1SRA6 popbyhpv1 = dNEW_4D_H1SRA6(popbyhpv1);

//  DMatrix4D popbyhiv = CREATE_DMATRIX4D(dims_isra, popbyhiv, 0.0);
//  DMatrix4D popbyhpv3 = CREATE_DMATRIX4D(dims_h3sra, popbyhpv3, 0.0);
//  DMatrix4D popbyhpv2 = CREATE_DMATRIX4D(dims_h2sra, popbyhpv2, 0.0);
//  DMatrix4D popbyhpv1 = CREATE_DMATRIX4D(dims_h1sra, popbyhpv1, 0.0);

  int* dims_sa = new int[3]{ 0, _SIZE_S, _SIZE_A6};
  DMatrix pships_age_sum = CREATE_DMATRIX(dims_sa, pships_age_sum, 0.0)
  delete[] dims_sa;

  DArray pships_com_bysex = CREATE_DARRAY(_SIZE_S, pships_com_bysex, 0.0)
  DArray pships_com_bysex_bal = CREATE_DARRAY(_SIZE_S, pships_com_bysex_bal, 0.0)
  DArray pshipsprop_age_f = CREATE_DARRAY(_SIZE_R * _SIZE_A6, pshipsprop_age_f, 0.0);
  DArray pshipsprop_age_m = CREATE_DARRAY(_SIZE_R * _SIZE_A6, pshipsprop_age_f, 0.0);

  dMAT_2D_RA6RA6 rhof_at = dNEW_2D_RA6RA6(rhof_at);
  dMAT_2D_RA6RA6 rhom_at = dNEW_2D_RA6RA6(rhom_at);
  dMAT_2D_RA6RA6 rhof = dNEW_2D_RA6RA6(rhof);
  dMAT_2D_RA6RA6 rhom = dNEW_2D_RA6RA6(rhom);
  dMAT_2D_RA6RA6 cpn_f = dNEW_2D_RA6RA6(cpn_f);
  dMAT_2D_RA6RA6 cpn_m = dNEW_2D_RA6RA6(cpn_m);
  dMAT_2D_RA6RA6 bal = dNEW_2D_RA6RA6(bal);
  dMAT_2D_RA6RA6 bal_f = dNEW_2D_RA6RA6(bal_f);
  dMAT_2D_RA6RA6 c_f = dNEW_2D_RA6RA6(c_f);
  dMAT_2D_RA6RA6 c_m = dNEW_2D_RA6RA6(c_m);

  DArray tempfoih1_comm = CREATE_DARRAY(_SIZE_R * _SIZE_A6, tempfoih1_comm, 0.0)
  DArray tempfoih2_comm = CREATE_DARRAY(_SIZE_R * _SIZE_A6, tempfoih2_comm, 0.0)
  DArray tempfoih3_comm = CREATE_DARRAY(_SIZE_R * _SIZE_A6, tempfoih3_comm, 0.0)
  DArray tempmoih1_comm = CREATE_DARRAY(_SIZE_R * _SIZE_A6, tempmoih1_comm, 0.0)
  DArray tempmoih2_comm = CREATE_DARRAY(_SIZE_R * _SIZE_A6, tempmoih2_comm, 0.0)
  DArray tempmoih3_comm = CREATE_DARRAY(_SIZE_R * _SIZE_A6, tempmoih3_comm, 0.0)
  DArray tempfoihpv1_comm = CREATE_DARRAY(_SIZE_R * _SIZE_A6, tempfoihpv1_comm, 0.0)
  DArray tempfoihpv2_comm = CREATE_DARRAY(_SIZE_R * _SIZE_A6, tempfoihpv2_comm, 0.0)
  DArray tempfoihpv3_comm = CREATE_DARRAY(_SIZE_R * _SIZE_A6, tempfoihpv3_comm, 0.0)
  DArray tempmoihpv1_comm = CREATE_DARRAY(_SIZE_R * _SIZE_A6, tempmoihpv1_comm, 0.0)
  DArray tempmoihpv2_comm = CREATE_DARRAY(_SIZE_R * _SIZE_A6, tempmoihpv2_comm, 0.0)
  DArray tempmoihpv3_comm = CREATE_DARRAY(_SIZE_R * _SIZE_A6, tempmoihpv3_comm, 0.0)

  DArray tempfoih1_comm_full = CREATE_DARRAY(_SIZE_R * _SIZE_A, tempfoih1_comm_full, 0.0)
  DArray tempfoih2_comm_full = CREATE_DARRAY(_SIZE_R * _SIZE_A, tempfoih2_comm_full, 0.0)
  DArray tempfoih3_comm_full = CREATE_DARRAY(_SIZE_R * _SIZE_A, tempfoih3_comm_full, 0.0)
  DArray tempmoih1_comm_full = CREATE_DARRAY(_SIZE_R * _SIZE_A, tempmoih1_comm_full, 0.0)
  DArray tempmoih2_comm_full = CREATE_DARRAY(_SIZE_R * _SIZE_A, tempmoih2_comm_full, 0.0)
  DArray tempmoih3_comm_full = CREATE_DARRAY(_SIZE_R * _SIZE_A, tempmoih3_comm_full, 0.0)
  DArray tempfoihpv1_comm_full = CREATE_DARRAY(_SIZE_R * _SIZE_A, tempfoihpv1_comm_full, 0.0)
  DArray tempfoihpv2_comm_full = CREATE_DARRAY(_SIZE_R * _SIZE_A, tempfoihpv2_comm_full, 0.0)
  DArray tempfoihpv3_comm_full = CREATE_DARRAY(_SIZE_R * _SIZE_A, tempfoihpv3_comm_full, 0.0)
  DArray tempmoihpv1_comm_full = CREATE_DARRAY(_SIZE_R * _SIZE_A, tempmoihpv1_comm_full, 0.0)
  DArray tempmoihpv2_comm_full = CREATE_DARRAY(_SIZE_R * _SIZE_A, tempmoihpv2_comm_full, 0.0)
  DArray tempmoihpv3_comm_full = CREATE_DARRAY(_SIZE_R * _SIZE_A, tempmoihpv3_comm_full, 0.0)


    // 3 is the size of ed - is this DIM_r?

  int* dims_3a = new int[3]{ 0, 3 * _SIZE_A6, 3 * _SIZE_A6 };
  DMatrix ed_at = CREATE_DMATRIX(dims_3a, ed_at, 0.0)

  int dim_ra = _SIZE_R * _SIZE_A6;
  DArray womsize = CREATE_DARRAY(dim_ra, womsize, 0.0)
  DArray mensize = CREATE_DARRAY(dim_ra, mensize, 0.0)
  DArray pcr_wom = CREATE_DARRAY(dim_ra, pcr_wom, 0.0)
  DArray pcr_men = CREATE_DARRAY(dim_ra, pcr_men, 0.0)

  dMAT_8D betahivacute_esc = dNEW_8D(betahivacute_esc)
  dMAT_8D betahiv2_esc = dNEW_8D(betahiv2_esc)
  dMAT_8D betahivart_esc = dNEW_8D(betahivart_esc)
  dMAT_8D betahivacute_esc_comm = dNEW_8D(betahivacute_esc_comm)
  dMAT_8D betahiv2_esc_comm = dNEW_8D(betahiv2_esc_comm)
  dMAT_8D betahivart_esc_comm = dNEW_8D(betahivart_esc_comm)
  dMAT_8D betahivpp_acute = dNEW_8D(betahivpp_acute)
  dMAT_8D betahivpp_hiv2 = dNEW_8D(betahivpp_hiv2)
  dMAT_8D betahivpp_art = dNEW_8D(betahivpp_art)
  dMAT_8D betahivpp_acute_comm = dNEW_8D(betahivpp_acute_comm)
  dMAT_8D betahivpp_hiv2_comm = dNEW_8D(betahivpp_hiv2_comm)
  dMAT_8D betahivpp_art_comm = dNEW_8D(betahivpp_art_comm)

  dMAT_2D_S_RA6 prev_acute = dNEW_2D_S_RA6(prev_acute);
  dMAT_2D_S_RA6 prev_noart = dNEW_2D_S_RA6(prev_noart);
  dMAT_2D_S_RA6 prev_art = dNEW_2D_S_RA6(prev_art);
  dMAT_2D_S_RA6 prev_hpv1 = dNEW_2D_S_RA6(prev_hpv1);
  dMAT_2D_S_RA6 prev_hpv2 = dNEW_2D_S_RA6(prev_hpv2);
  dMAT_2D_S_RA6 prev_hpv3 = dNEW_2D_S_RA6(prev_hpv3);

  dMAT_8D pdt1 = dNEW_8D(pdt1);
  dMAT_8D pdt2 = dNEW_8D(pdt2);
  dMAT_8D pdt3 = dNEW_8D(pdt3);
  dMAT_8D pdt4 = dNEW_8D(pdt4);
  dMAT_8D pdt5 = dNEW_8D(pdt5);
  dMAT_8D pdt6 = dNEW_8D(pdt6);
  dMAT_8D pdt7 = dNEW_8D(pdt7);

  R->dim_y_sav = new int[3]{ 0, (int)round(p->dp0[32] + 1), _SIZE_S * _SIZE_A * _SIZE_V1 };
  R->hpv1cum_pop = CREATE_DMATRIX(R->dim_y_sav, R->hpv1cum_pop, 0.0);
  R->hpv2cum_pop = CREATE_DMATRIX(R->dim_y_sav, R->hpv2cum_pop, 0.0);
  R->hpv3cum_pop = CREATE_DMATRIX(R->dim_y_sav, R->hpv3cum_pop, 0.0);
  R->hpv1hncum_pop = CREATE_DMATRIX(R->dim_y_sav, R->hpv1hncum_pop, 0.0);
  R->hpv2hncum_pop = CREATE_DMATRIX(R->dim_y_sav, R->hpv2hncum_pop, 0.0);
  R->hpv3hncum_pop = CREATE_DMATRIX(R->dim_y_sav, R->hpv3hncum_pop, 0.0);


  dMAT_8D ccinc_pop = dNEW_8D(ccinc_pop);
  dMAT_8D cin_pop = dNEW_8D(cin_pop);
  dMAT_8D hivdeaths_pop = dNEW_8D(hivdeaths_pop);
  dMAT_8D ccdeaths_pop = dNEW_8D(ccdeaths_pop);
  dMAT_8D othdeaths_pop = dNEW_8D(othdeaths_pop);

  R->dim_gg_sav = new int[5]{ 0, _hyrid_length, dims[DIM_s], dims[DIM_a], dims[DIM_v1] };
  R->dim_gg_1sav = new int[6]{ 0, _hyrid_length, 1, dims[DIM_s], dims[DIM_a], dims[DIM_v1] };
  R->hivsus_pop = CREATE_DMATRIX5D(R->dim_gg_1sav, R->hivsus_pop, 0.0); // or is it 5d?
  R->dim_gg_s = new int[3]{ 0, _hyrid_length, dims[DIM_s] };
  R->hiv_15to49 = CREATE_DMATRIX(R->dim_gg_s, R->hiv_15to49, 0.0);
  R->hiv_15to74 = CREATE_DMATRIX(R->dim_gg_s, R->hiv_15to74, 0.0);
  R->hiv_9to14 = CREATE_DMATRIX(R->dim_gg_s, R->hiv_9to14, 0.0);
  R->hivart_15to49 = CREATE_DMATRIX(R->dim_gg_s, R->hivart_15to49, 0.0);
  R->hivart_15to74 = CREATE_DMATRIX(R->dim_gg_s, R->hivart_15to74, 0.0);
  R->hivart_9to14 = CREATE_DMATRIX(R->dim_gg_s, R->hivart_9to14, 0.0);
  R->art_15plus = CREATE_DMATRIX(R->dim_gg_s, R->art_15plus, 0.0);
  DArray sum1_s = CREATE_DARRAY(dims[DIM_s], sum1_s, 0);
  DArray sum2_s = CREATE_DARRAY(dims[DIM_s], sum2_s, 0);
  DArray sum1_a = CREATE_DARRAY(dims[DIM_a], sum1_a, 0);
  DArray sum2_a = CREATE_DARRAY(dims[DIM_a], sum2_a, 0);

  R->dim_gg_a5 = new int[3] {0, _hyrid_length, 5};
  R->art_age = CREATE_DMATRIX(R->dim_gg_a5, R->art_age, 0.0);
  R->hiv_fage = CREATE_DMATRIX(R->dim_gg_a5, R->hiv_fage, 0.0);
  R->hiv_mage = CREATE_DMATRIX(R->dim_gg_a5, R->hiv_mage, 0.0);

  R->dim_gg_a = new int[3]{ 0, _hyrid_length, dims[DIM_a] };
  R->vacc_age_f = CREATE_DMATRIX(R->dim_gg_a, R->vacc_age_f, 0.0);
  R->vacc_age_m = CREATE_DMATRIX(R->dim_gg_a, R->vacc_age_m, 0.0);
  R->vacc_age_f_pos = CREATE_DMATRIX(R->dim_gg_a, R->vacc_age_f_pos, 0.0);
  R->vacc_age_f_art = CREATE_DMATRIX(R->dim_gg_a, R->vacc_age_f_art, 0.0);
  R->vacc_age_m_pos = CREATE_DMATRIX(R->dim_gg_a, R->vacc_age_m_pos, 0.0);
  R->hiv_fsw = CREATE_DARRAY(_hyrid_length, R->hiv_fsw, 0.0)
  R->hiv_mcli = CREATE_DARRAY(_hyrid_length, R->hiv_mcli, 0.0)
  R->testtime = CREATE_DARRAY(_hyrid_length, R->testtime, 0.0)

  R->dim_g_s = new int[3]{ 0, tend + 1, dims[DIM_s] };
  R->dim_g_a = new int[3]{ 0, tend + 1, dims[DIM_a] };
  R->psize = CREATE_DMATRIX(R->dim_g_s, R->psize, 0.0);
  R->psize_age_m = CREATE_DMATRIX(R->dim_g_a, R->psize_age_m, 0.0);
  R->psize_age_m_neg = CREATE_DMATRIX(R->dim_g_a, R->psize_age_m_neg, 0.0);
  R->psize_age_m_pos = CREATE_DMATRIX(R->dim_g_a, R->psize_age_m_pos, 0.0);
  R->psize_age_m_art = CREATE_DMATRIX(R->dim_g_a, R->psize_age_m_art, 0.0);
  R->psize_age_f = CREATE_DMATRIX(R->dim_g_a, R->psize_age_f, 0.0);
  R->psize_age_f_neg = CREATE_DMATRIX(R->dim_g_a, R->psize_age_f_neg, 0.0);
  R->psize_age_f_pos = CREATE_DMATRIX(R->dim_g_a, R->psize_age_f_pos, 0.0);
  R->psize_age_f_art = CREATE_DMATRIX(R->dim_g_a, R->psize_age_f_art, 0.0);

  R->dim_g_i_a = new int[4]{ 0, tend + 1, _SIZE_I, _SIZE_A };
  R->hpv1618_f = CREATE_DMATRIX3D(R->dim_g_i_a, R->hpv1618_f, 0.0)
  R->hpv1618_m = CREATE_DMATRIX3D(R->dim_g_i_a, R->hpv1618_m, 0.0)
  R->nvthpv_f = CREATE_DMATRIX3D(R->dim_g_i_a, R->nvthpv_f, 0.0)
  R->nvthpv_m = CREATE_DMATRIX3D(R->dim_g_i_a, R->nvthpv_m, 0.0)

  R->hpv9vt_f = CREATE_DMATRIX3D(R->dim_g_i_a, R->hpv9vt_f, 0.0)
  R->hpv9vt_m = CREATE_DMATRIX3D(R->dim_g_i_a, R->hpv9vt_m, 0.0)
  R->hpv_f = CREATE_DMATRIX3D(R->dim_g_i_a, R->hpv_f, 0.0)
  R->hpv_m = CREATE_DMATRIX3D(R->dim_g_i_a, R->hpv_m, 0.0)
  R->cinpr_f = CREATE_DMATRIX3D(R->dim_g_i_a, R->cinpr_f, 0.0)
  R->ccpr_f = CREATE_DMATRIX3D(R->dim_g_i_a, R->ccpr_f, 0.0)

  R->dim_g_r_a5 = new int[4]{ 0, tend + 1, _SIZE_R, 5};

//  R->hpv9vt_f = CREATE_DMATRIX3D(R->dim_g_r_a5, R->hpv9vt_f, 0.0)
//  R->hpv9vt_m = CREATE_DMATRIX3D(R->dim_g_r_a5, R->hpv9vt_m, 0.0)
//  R->hpv_f = CREATE_DMATRIX3D(R->dim_g_r_a5, R->hpv_f, 0.0)
//  R->hpv_m = CREATE_DMATRIX3D(R->dim_g_r_a5, R->hpv_m, 0.0)
//  R->cinpr_f = CREATE_DMATRIX3D(R->dim_g_r_a5, R->cinpr_f, 0.0)
 // R->ccpr_f = CREATE_DMATRIX3D(R->dim_g_r_a5, R->ccpr_f, 0.0)

  R->testtime1yr = CREATE_DARRAY(tend + 1, R->testtime1yr, 0.0)

  R->dim_g_r = new int[3]{ 0, _h5yrid_length, dims[DIM_r] };
  R->popbyrisk_f = CREATE_DMATRIX(R->dim_g_r, R->popbyrisk_f, 0.0)
  R->popbyrisk_m = CREATE_DMATRIX(R->dim_g_r, R->popbyrisk_m, 0.0)

  DArray ageinf = CREATE_DARRAY(3, ageinf, 0.0)
  DArray ageinm = CREATE_DARRAY(3, ageinm, 0.0)
  double _h15_16 = p->hp0[15] * p->hp0[16];


  //  for n = 1:niterations
  // int timer0 = -clock();
  for (int n = 1; n <= niterations; n++) {

    //    % %Introduction of HIV
    //    %
    //    if (n == tstep * (1985 - dp0(22))) %
    //      %(h1, h2, h3, i, s, r, a, v1)
    //      % keyboard
    //      popn(1, 1, 1, 2, :, 1, 2, 1) = hp0(24) * 0.5;%
    //      popn(2, 2, 2, 2, :, 1, 2, 1) = hp0(24) * 0.5;%

    //    end

        /************ Above, two lines are exactly repeated *******************/

    if (n == (int)round(tstep * ((int)(p->hp0[14]) - p->dp0[31]))) {
      for (int _ds = 1; _ds <= dims[DIM_s]; _ds++) {
        GET(popn, 1, 1, 1, 2, _ds, 1, 2, 1) = p->hp0[24] * 0.5;
        GET(popn, 2, 2, 2, 2, _ds, 1, 2, 1) = p->hp0[24] * 0.5;
      }
      // popn might be updated here...

      sum_h1(popn, popn_sum_h1);
      sum_h2(popn_sum_h1, popn_sum_h1_h2);
      sum_h3(popn_sum_h1_h2, popn_sum_h1_h2_h3);
      sum_i(popn_sum_h1_h2_h3, popn_sum_h1_h2_h3_i);
    }

    // We often do popn(:, :, :) or even popn(:, :, :, :) - so calculate these once per timestep.
    // Initialise here - but we actually need it at the end of the timestep.

    if (n == 1) {
      sum_h1(popn, popn_sum_h1);
      sum_h2(popn_sum_h1, popn_sum_h1_h2);
      sum_h3(popn_sum_h1_h2, popn_sum_h1_h2_h3);
      sum_i(popn_sum_h1_h2_h3, popn_sum_h1_h2_h3_i);
    }


    //    %% VACCINATION(intervention scenarios)
    //      if vaccscen.setup > 0
    //        if (n > tstep* (2020 - dp0(22))& n < tstep * (2023 - dp0(22)) + 1) % start vaccination at 2019

    if ((int)round(p->vaccscen_setup) > 0) {
      if (n > tstep* (2020 - p->dp0[31]) && n < tstep * (2022 - p->dp0[31]) + 1) {

        //          % if n > tstep* (2020 - dp0(22)) %
        //          %first year burn in period to assure rapid vaccination scale up
        //          % (h1, h2, h3, i, s, r, a, v1)
        //          % desired vacc coverage pp0(33) - (1 / 1) * LN(1 - K31)

        //          vxx(:, : , : , : , 1, : , 1, 1) = vaccscen.agef1; %; -log(1 - pp0(33));

        LOOP_H1 LOOP_H2 LOOP_H3 LOOP_I LOOP_R
          GET(vxx, _dh1, _dh2, _dh3, _di, 1, _dr, 1, 1) = p->vaccscen_agef1;

        //          if vaccscen.setup == 2
        //            vxx(:, : , : , 2 : 4, 1, : , 2, 1) = vaccscen.agef2hp / 2;%% N.B.divided by 2 here to reduce speed of scale up
        //            vxx(:, : , : , 2 : 4, 1, : , 3, 1) = vaccscen.agef3hp / 2;
        //            vxx(:, : , : , 2 : 4, 1, : , 4, 1) = vaccscen.agef4hp / 2;
        //          end

        if ((int)round(p->vaccscen_setup == 2)) {
          LOOP_H1 LOOP_H2 LOOP_H3 LOOP_R
            for (int _di = 2; _di <= 4; _di++) {
              GET(vxx, _dh1, _dh2, _dh3, _di, 1, _dr, 2, 1) = p->vaccscen_agef2hp ;
              GET(vxx, _dh1, _dh2, _dh3, _di, 1, _dr, 3, 1) = p->vaccscen_agef3hp ;
              GET(vxx, _dh1, _dh2, _dh3, _di, 1, _dr, 4, 1) = p->vaccscen_agef4hp ;
            }
        }
        //        end
      }

      //
	   if ((int)round(p->vaccscen_setup == 2)) {
		   if (p->vaccscen_agef2hp < 0.5 && n > tstep* (2023 - p->dp0[31]) ) {

			LOOP_H1 LOOP_H2 LOOP_H3 LOOP_R
			  for (int _di = 2; _di <= 4; _di++) {
				GET(vxx, _dh1, _dh2, _dh3, _di, 1, _dr, 2, 1) =  0;
				GET(vxx, _dh1, _dh2, _dh3, _di, 1, _dr, 3, 1) =  0;
				GET(vxx, _dh1, _dh2, _dh3, _di, 1, _dr, 4, 1) =  0;

				// if (n == tstep* (2023 - p->dp0[31]) +10)	{
				//		printf("test1");
				//	  }
		   }
		  }
		   if (p->vaccscen_agef2hp > 0.5 && n > tstep* (2026 - p->dp0[31]) && n < tstep * (2057 - p->dp0[31]) + 1) {
            LOOP_H1 LOOP_H2 LOOP_H3 LOOP_R
              for (int _di = 2; _di <= 4; _di++) {
                GET(vxx, _dh1, _dh2, _dh3, _di, 1, _dr, 2, 1) =  p->vaccscen_agef2hp /10;
                GET(vxx, _dh1, _dh2, _dh3, _di, 1, _dr, 3, 1) =  p->vaccscen_agef3hp /10;
                GET(vxx, _dh1, _dh2, _dh3, _di, 1, _dr, 4, 1) =  p->vaccscen_agef4hp /10;

				// if (n == tstep* (2026 - p->dp0[31]) +10)	{
				//		printf("test2");
				//	  }
		   }
         }

		   if (n > tstep* (2057 - p->dp0[31])) {
            LOOP_H1 LOOP_H2 LOOP_H3 LOOP_R
              for (int _di = 2; _di <= 4; _di++) {
                GET(vxx, _dh1, _dh2, _dh3, _di, 1, _dr, 2, 1) = 0;// p->vaccscen_agef2hp /10;
                GET(vxx, _dh1, _dh2, _dh3, _di, 1, _dr, 3, 1) = 0;// p->vaccscen_agef3hp /10;
                GET(vxx, _dh1, _dh2, _dh3, _di, 1, _dr, 4, 1) = 0;// p->vaccscen_agef4hp /10;

				// if (n == tstep* (2057 - p->dp0[31]) +10)	{
				//		printf("test3");
				//	  }
              }
          }
	   }
      //        if (n > tstep* (2023 - dp0(22)))
      if (n > tstep* (2022 - p->dp0[31])) {
        //          % after 3 first years burn - in start reducing vaccination impact
        //
        //          Nvfage1 = reshape(sum(sum(sum(sum(sum(popn(:, : , : , : , 1, : , 1, : ), 1), 2), 3), 4), 6), 1, 2);

        double Nvfage1_1 = 0.0;
        double Nvfage1_2 = 0.0;

        LOOP_R {
          Nvfage1_1 += GET(popn_sum_h1_h2_h3_i, 1, 1, 1, 1, 1, _dr, 1, 1);
          Nvfage1_2 += GET(popn_sum_h1_h2_h3_i, 1, 1, 1 ,1 ,1 ,_dr, 1, 2);
        }
        double Nvfage1sum = Nvfage1_1 + Nvfage1_2;

        //          vaccrate1 = (vaccscen.agef1 * sum(Nvfage1) - Nvfage1(1, 2)) / Nvfage1(1, 1);
        //          if vaccrate1 < 0
        //            vaccrate1 = 0;
        //            % keyboard
        //          end

        double vaccrate1 = (p->vaccscen_agef1 * Nvfage1sum - Nvfage1_2) / Nvfage1_1;
        if (vaccrate1 < 0) vaccrate1 = 0;

        //          vxx(:, : , : , : , 1, : , 1, 1) = vaccrate1;
        LOOP_H1 LOOP_H2 LOOP_H3 LOOP_I LOOP_R
          GET(vxx, _dh1, _dh2, _dh3, _di, 1, _dr, 1, 1) = vaccrate1;


        //          if vaccscen.setup == 2

//// NOT IN USE ////////////
      //  if ((int)round(p->vaccscen_setup) == 9999) {
          //            if vaccscen.agef2hp > 0
          //              Nvfage2hp = reshape(sum(sum(sum(sum(sum(popn(:, : , : , 2 : 4, 1, : , 2, : ), 1), 2), 3), 4), 6), 1, 2);
          //              vaccrate2 = (vaccscen.agef2hp * sum(Nvfage2hp) - Nvfage2hp(1, 2)) / Nvfage2hp(1, 1);
          //              if vaccrate2 < 0
          //                vaccrate2 = 0;
          //                % keyboard
          //              end
          //              vxx(:, : , : , 3 : 4, 1, : , 2, 1) = vaccrate2;
          //            end

     //     if ((int)round(p->vaccscen_agef2hp) > 9999) {
     //       LOOP_LINEAR GET_LIN_ALL(spare) = GET_LIN_ALL(popn_sum_h1_h2_h3);

     //       for (int _di = 3; _di <= 3; _di++) LOOP_R LOOP_V1
     //         GET(spare, 1, 1, 1, 2, 1, _dr, 2, _dv1) += GET(spare, 1, 1, 1, _di, 1, _dr, 2, _dv1);

     //       for (int _dr = 2; _dr <= dims[DIM_r]; _dr++) LOOP_V1
     //         GET(spare, 1, 1, 1, 2, 1, 1, 2, _dv1) += GET(spare, 1, 1, 1, 2, 1, _dr, 2, _dv1);

     //       double Nvfage2hp_11 = GET(spare, 1, 1, 1, 2, 1, 1, 2, 1);
     //       double Nvfage2hp_21 = GET(spare, 1, 1, 1, 2, 1, 1, 2, 2);
     //       double Nvfage2hp_sum1 = Nvfage2hp_11 + Nvfage2hp_21;

     //       double vaccrate2_1 = (p->vaccscen_agef2hp * Nvfage2hp_sum1 - Nvfage2hp_21) / Nvfage2hp_11;
     //       if (vaccrate2_1 < 0) vaccrate2_1 = 0;
     //       LOOP_H1 LOOP_H2 LOOP_H3 LOOP_R
     //         for (int _di = 3; _di <= 3; _di++)
     //          GET(vxx, _dh1, _dh2, _dh3, _di, 1, _dr, 2, 1) = vaccrate2_1;


            /////

     //       for (int _di = 4; _di <= 4; _di++) LOOP_R LOOP_V1
     //                  GET(spare, 1, 1, 1, 2, 1, _dr, 2, _dv1) += GET(spare, 1, 1, 1, _di, 1, _dr, 2, _dv1);

     //                for (int _dr = 2; _dr <= dims[DIM_r]; _dr++) LOOP_V1
     //                  GET(spare, 1, 1, 1, 2, 1, 1, 2, _dv1) += GET(spare, 1, 1, 1, 2, 1, _dr, 2, _dv1);

     //                double Nvfage2hp_1 = GET(spare, 1, 1, 1, 2, 1, 1, 2, 1);
     //                double Nvfage2hp_2 = GET(spare, 1, 1, 1, 2, 1, 1, 2, 2);
     //                double Nvfage2hp_sum = Nvfage2hp_1 + Nvfage2hp_2;

      //               double vaccrate2 = (p->vaccscen_agef2hp * Nvfage2hp_sum - Nvfage2hp_2) / Nvfage2hp_1;
      //               if (vaccrate2 < 0) vaccrate2 = 0;
      //               LOOP_H1 LOOP_H2 LOOP_H3 LOOP_R
      //                 for (int _di = 4; _di <= 4; _di++)
      //                   GET(vxx, _dh1, _dh2, _dh3, _di, 1, _dr, 2, 1) = vaccrate2;

      //    }

          //            if vaccscen.agef3hp > 0
          //              Nvfage3hp = reshape(sum(sum(sum(sum(sum(popn(:, : , : , 2 : 4, 1, : , 3, : ), 1), 2), 3), 4), 6), 1, 2);
          //              vaccrate3 = (vaccscen.agef3hp * sum(Nvfage3hp) - Nvfage3hp(1, 2)) / Nvfage3hp(1, 1);
          //              if vaccrate3 < 0
          //                vaccrate3 = 0;
          //                % keyboard
          //              end
          //              vxx(:, : , : , 2 : 4, 1, : , 3, 1) = vaccrate3;
          //            end

     //     if ((int)round(p->vaccscen_agef3hp) > 9999) {

     //       LOOP_LINEAR GET_LIN_ALL(spare) = GET_LIN_ALL(popn_sum_h1_h2_h3);

      //      for (int _di = 3; _di <= 3; _di++) LOOP_R LOOP_V1
      //        GET(spare, 1, 1, 1, 2, 1, _dr, 3, _dv1) += GET(spare, 1, 1, 1, _di, 1, _dr, 3, _dv1);

     //       for (int _dr = 2; _dr <= dims[DIM_r]; _dr++) LOOP_V1
     //         GET(spare, 1, 1, 1, 2, 1, 1, 3, _dv1) += GET(spare, 1, 1, 1, 2, 1, _dr, 3, _dv1);

     //       double Nvfage3hp_11 = GET(spare, 1, 1, 1, 2, 1, 1, 3, 1);
     //       double Nvfage3hp_21 = GET(spare, 1, 1, 1, 2, 1, 1, 3, 2);
     //       double Nvfage3hp_sum1 = Nvfage3hp_11 + Nvfage3hp_21;

    //        double vaccrate3_1 = (p->vaccscen_agef3hp * Nvfage3hp_sum1 - Nvfage3hp_21) / Nvfage3hp_11;
    //        if (vaccrate3_1 < 0) vaccrate3_1 = 0;
    //        LOOP_H1 LOOP_H2 LOOP_H3 LOOP_R
    //          for (int _di = 3; _di <= 3; _di++)
    //            GET(vxx, _dh1, _dh2, _dh3, _di, 1, _dr, 3, 1) = vaccrate3_1 ;

            ////

    //        LOOP_LINEAR GET_LIN_ALL(spare) = GET_LIN_ALL(popn_sum_h1_h2_h3);

     //                  for (int _di = 4; _di <= 4; _di++) LOOP_R LOOP_V1
    //                     GET(spare, 1, 1, 1, 2, 1, _dr, 3, _dv1) += GET(spare, 1, 1, 1, _di, 1, _dr, 3, _dv1);

     //                  for (int _dr = 2; _dr <= dims[DIM_r]; _dr++) LOOP_V1
     //                    GET(spare, 1, 1, 1, 2, 1, 1, 3, _dv1) += GET(spare, 1, 1, 1, 2, 1, _dr, 3, _dv1);

     //                  double Nvfage3hp_1 = GET(spare, 1, 1, 1, 2, 1, 1, 3, 1);
     //                  double Nvfage3hp_2 = GET(spare, 1, 1, 1, 2, 1, 1, 3, 2);
     //                  double Nvfage3hp_sum = Nvfage3hp_1 + Nvfage3hp_2;

      //                 double vaccrate3 = (p->vaccscen_agef3hp * Nvfage3hp_sum - Nvfage3hp_2) / Nvfage3hp_1;
      //                 if (vaccrate3 < 0) vaccrate3 = 0;
      //                 LOOP_H1 LOOP_H2 LOOP_H3 LOOP_R
      //                   for (int _di = 4; _di <= 4; _di++)
      //                     GET(vxx, _dh1, _dh2, _dh3, _di, 1, _dr, 3, 1) = vaccrate3 ;

      //    }

          //            if vaccscen.agef4hp > 0
          //              Nvfage4hp = reshape(sum(sum(sum(sum(sum(popn(:, : , : , 2 : 4, 1, : , 4, : ), 1), 2), 3), 4), 6), 1, 2);
          //              vaccrate4 = (vaccscen.agef4hp * sum(Nvfage4hp) - Nvfage4hp(1, 2)) / Nvfage4hp(1, 1);

          //              if vaccrate4 < 0
          //                vaccrate4 = 0;
          //                % keyboard
          //              end
          //              vxx(:, : , : , 2 : 4, 1, : , 4, 1) = vaccrate4;
          //            end

       //   if ((int)round(p->vaccscen_agef4hp > 9999)) {

       //     LOOP_LINEAR GET_LIN_ALL(spare) = GET_LIN_ALL(popn_sum_h1_h2_h3);

       //     for (int _di = 2; _di <= 4; _di++) LOOP_R LOOP_V1
       //       GET(spare, 1, 1, 1, 2, 1, _dr, 4, _dv1) += GET(spare, 1, 1, 1, _di, 1, _dr, 4, _dv1);

       //     for (int _dr = 2; _dr <= dims[DIM_r]; _dr++) LOOP_V1
       //       GET(spare, 1, 1, 1, 2, 1, 1, 4, _dv1) += GET(spare, 1, 1, 1, 2, 1, _dr, 4, _dv1);

      //      double Nvfage4hp_1 = GET(spare, 1, 1, 1, 2, 1, 1, 4, 1);
      //      double Nvfage4hp_2 = GET(spare, 1, 1, 1, 2, 1, 1, 4, 2);
      //      double Nvfage4hp_sum = Nvfage4hp_1 + Nvfage4hp_2;
      //       double vaccrate4 = (p->vaccscen_agef4hp * Nvfage4hp_sum - Nvfage4hp_2) / Nvfage4hp_1;
      //      if (vaccrate4 < 0) vaccrate4 = 0;
     //       LOOP_H1 LOOP_H2 LOOP_H3 LOOP_R
     //         for (int _di = 3; _di <= 4; _di++)
     //           GET(vxx, _dh1, _dh2, _dh3, _di, 1, _dr, 4, 1) = vaccrate4 ;
     //     }
///////////////////////////////////////////////////////////////////////////////////
          //          end ** END if vaccscen.setup == 2
          //        end ** END if (n > step...)
          //      end ** END if vaccscesn.setup > 0

       // }
      }
    }


    //      %%
    //      % poptot = sum(sum(sum(sum(sum(sum(sum(popn(:, : , : , : , : , : , : ))))))));

    //      new_age_structure

    // popbysexriskage(:, :, 1) = reshape(sum(sum(sum(sum(sum(sum(popn(:, :, :, :, 1:s, 1:r, 1,     : ), 1), 2), 3), 4), 7), 8), s, r, 1);
    // popbysexriskage(:, :, 2) = reshape(sum(sum(sum(sum(sum(sum(popn(:, :, :, :, 1:s, 1:r, 2:3,   : ), 1), 2), 3), 4), 7), 8), s, r, 1);
    // popbysexriskage(:, :, 3) = reshape(sum(sum(sum(sum(sum(sum(popn(:, :, :, :, 1:s, 1:r, 4:5,   : ), 1), 2), 3), 4), 7), 8), s, r, 1);
    // popbysexriskage(:, :, 4) = reshape(sum(sum(sum(sum(sum(sum(popn(:, :, :, :, 1:s, 1:r, 6:8,   : ), 1), 2), 3), 4), 7), 8), s, r, 1);
    // popbysexriskage(:, :, 5) = reshape(sum(sum(sum(sum(sum(sum(popn(:, :, :, :, 1:s, 1:r, 9:10,  : ), 1), 2), 3), 4), 7), 8), s, r, 1);
    // popbysexriskage(:, :, 6) = reshape(sum(sum(sum(sum(sum(sum(popn(:, :, :, :, 1:s, 1:r, 11:13, : ), 1), 2), 3), 4), 7), 8), s, r, 1);

    sum_1234_8_sra_a6(spare_v1, popn_sum_h1_h2_h3_i, 1, _SIZE_S, 1, _SIZE_R, 1, 1, popbysexriskage, 1);
    sum_1234_8_sra_a6(spare_v1, popn_sum_h1_h2_h3_i, 1, _SIZE_S, 1, _SIZE_R, 2, 3, popbysexriskage, 2);
    sum_1234_8_sra_a6(spare_v1, popn_sum_h1_h2_h3_i, 1, _SIZE_S, 1, _SIZE_R, 4, 5, popbysexriskage, 3);
    sum_1234_8_sra_a6(spare_v1, popn_sum_h1_h2_h3_i, 1, _SIZE_S, 1, _SIZE_R, 6, 8, popbysexriskage, 4);
    sum_1234_8_sra_a6(spare_v1, popn_sum_h1_h2_h3_i, 1, _SIZE_S, 1, _SIZE_R, 9, 10, popbysexriskage, 5);
    sum_1234_8_sra_a6(spare_v1, popn_sum_h1_h2_h3_i, 1, _SIZE_S, 1, _SIZE_R, 11, 13, popbysexriskage, 6);

    // popbyhiv(:, : , : , 1) = reshape(sum(sum(sum(sum(sum(popn(:, : , : , 1 : i, 1 : s, 1 : r, 1, : ), 1), 2), 3), 7), 8), i, s, r, 1);
    // popbyhiv(:, : , : , 2) = reshape(sum(sum(sum(sum(sum(popn(:, : , : , 1 : i, 1 : s, 1 : r, 2 : 3, : ), 1), 2), 3), 7), 8), i, s, r, 1);
    // popbyhiv(:, : , : , 3) = reshape(sum(sum(sum(sum(sum(popn(:, : , : , 1 : i, 1 : s, 1 : r, 4 : 5, : ), 1), 2), 3), 7), 8), i, s, r, 1);
    // popbyhiv(:, : , : , 4) = reshape(sum(sum(sum(sum(sum(popn(:, : , : , 1 : i, 1 : s, 1 : r, 6 : 8, : ), 1), 2), 3), 7), 8), i, s, r, 1);
    // popbyhiv(:, : , : , 5) = reshape(sum(sum(sum(sum(sum(popn(:, : , : , 1 : i, 1 : s, 1 : r, 9 : 10, : ), 1), 2), 3), 7), 8), i, s, r, 1);
    // popbyhiv(:, : , : , 6) = reshape(sum(sum(sum(sum(sum(popn(:, : , : , 1 : i, 1 : s, 1 : r, 11 : 13, : ), 1), 2), 3), 7), 8), i, s, r, 1);

    sum123_8_isr_a(spare_v1, popn_sum_h1_h2_h3, 1, _SIZE_I, 1, _SIZE_S, 1, _SIZE_R, 1, 1, popbyhiv, 1);
    sum123_8_isr_a(spare_v1, popn_sum_h1_h2_h3, 1, _SIZE_I, 1, _SIZE_S, 1, _SIZE_R, 2, 3, popbyhiv, 2);
    sum123_8_isr_a(spare_v1, popn_sum_h1_h2_h3, 1, _SIZE_I, 1, _SIZE_S, 1, _SIZE_R, 4, 5, popbyhiv, 3);
    sum123_8_isr_a(spare_v1, popn_sum_h1_h2_h3, 1, _SIZE_I, 1, _SIZE_S, 1, _SIZE_R, 6, 8, popbyhiv, 4);
    sum123_8_isr_a(spare_v1, popn_sum_h1_h2_h3, 1, _SIZE_I, 1, _SIZE_S, 1, _SIZE_R, 9, 10, popbyhiv, 5);
    sum123_8_isr_a(spare_v1, popn_sum_h1_h2_h3, 1, _SIZE_I, 1, _SIZE_S, 1, _SIZE_R, 11, 13, popbyhiv, 6);

    // popbyhpv3(:,:,:,1) = reshape(sum(sum(sum(sum(sum(popn(:, :, 1:h3, :, 1:s, 1:r, 1,   :),1),2),4),7),8), h3, s, r, 1);
    // popbyhpv3(:,:,:,2) = reshape(sum(sum(sum(sum(sum(popn(:, : , 1 : h3, : , 1 : s, 1 : r, 2 : 3, : ), 1), 2), 4), 7), 8), h3, s, r, 1);
    // popbyhpv3(:,:,:,3) = reshape(sum(sum(sum(sum(sum(popn(:, : , 1 : h3, : , 1 : s, 1 : r, 4 : 5, : ), 1), 2), 4), 7), 8), h3, s, r, 1);
    // popbyhpv3(:,:,:,4) = reshape(sum(sum(sum(sum(sum(popn(:, : , 1 : h3, : , 1 : s, 1 : r, 6 : 8, : ), 1), 2), 4), 7), 8), h3, s, r, 1);
    // popbyhpv3(:,:,:,5) = reshape(sum(sum(sum(sum(sum(popn(:, : , 1 : h3, : , 1 : s, 1 : r, 9 : 10, : ), 1), 2), 4), 7), 8), h3, s, r, 1);
    // popbyhpv3(:,:,:,6) = reshape(sum(sum(sum(sum(sum(popn(:, : , 1 : h3, : , 1 : s, 1 : r, 11 : 13, : ), 1), 2), 4), 7), 8), h3, s, r, 1);

    sum12_48_h3sr_a(spare_av, popn_sum_h1_h2, 1, _SIZE_H3, 1, _SIZE_S, 1, _SIZE_R, 1, 1, popbyhpv3, 1);
    sum12_48_h3sr_a(spare_av, popn_sum_h1_h2, 1, _SIZE_H3, 1, _SIZE_S, 1, _SIZE_R, 2, 3, popbyhpv3, 2);
    sum12_48_h3sr_a(spare_av, popn_sum_h1_h2, 1, _SIZE_H3, 1, _SIZE_S, 1, _SIZE_R, 4, 5, popbyhpv3, 3);
    sum12_48_h3sr_a(spare_av, popn_sum_h1_h2, 1, _SIZE_H3, 1, _SIZE_S, 1, _SIZE_R, 6, 8, popbyhpv3, 4);
    sum12_48_h3sr_a(spare_av, popn_sum_h1_h2, 1, _SIZE_H3, 1, _SIZE_S, 1, _SIZE_R, 9, 10, popbyhpv3, 5);
    sum12_48_h3sr_a(spare_av, popn_sum_h1_h2, 1, _SIZE_H3, 1, _SIZE_S, 1, _SIZE_R, 11, 13, popbyhpv3, 6);

    // popbyhpv2(:, : , : , 1) = reshape(sum(sum(sum(sum(sum(popn(:, 1 : h2, : , : , 1 : s, 1 : r, 1, : ), 1), 3), 4), 7), 8), h2, s, r, 1);
    // popbyhpv2(:, : , : , 2) = reshape(sum(sum(sum(sum(sum(popn(:, 1 : h2, : , : , 1 : s, 1 : r, 2 : 3, : ), 1), 3), 4), 7), 8), h2, s, r, 1);
    // popbyhpv2(:, : , : , 3) = reshape(sum(sum(sum(sum(sum(popn(:, 1 : h2, : , : , 1 : s, 1 : r, 4 : 5, : ), 1), 3), 4), 7), 8), h2, s, r, 1);
    // popbyhpv2(:, : , : , 4) = reshape(sum(sum(sum(sum(sum(popn(:, 1 : h2, : , : , 1 : s, 1 : r, 6 : 8, : ), 1), 3), 4), 7), 8), h2, s, r, 1);
    // popbyhpv2(:, : , : , 5) = reshape(sum(sum(sum(sum(sum(popn(:, 1 : h2, : , : , 1 : s, 1 : r, 9 : 10, : ), 1), 3), 4), 7), 8), h2, s, r, 1);
    // popbyhpv2(:, : , : , 6) = reshape(sum(sum(sum(sum(sum(popn(:, 1 : h2, : , : , 1 : s, 1 : r, 11 : 13, : ), 1), 3), 4), 7), 8), h2, s, r, 1);

    sum1_348_h2sr_a(spare_iav, popn_sum_h1, 1, _SIZE_H2, 1, _SIZE_S, 1, _SIZE_R, 1, 1, popbyhpv2, 1);
    sum1_348_h2sr_a(spare_iav, popn_sum_h1, 1, _SIZE_H2, 1, _SIZE_S, 1, _SIZE_R, 2, 3, popbyhpv2, 2);
    sum1_348_h2sr_a(spare_iav, popn_sum_h1, 1, _SIZE_H2, 1, _SIZE_S, 1, _SIZE_R, 4, 5, popbyhpv2, 3);
    sum1_348_h2sr_a(spare_iav, popn_sum_h1, 1, _SIZE_H2, 1, _SIZE_S, 1, _SIZE_R, 6, 8, popbyhpv2, 4);
    sum1_348_h2sr_a(spare_iav, popn_sum_h1, 1, _SIZE_H2, 1, _SIZE_S, 1, _SIZE_R, 9, 10, popbyhpv2, 5);
    sum1_348_h2sr_a(spare_iav, popn_sum_h1, 1, _SIZE_H2, 1, _SIZE_S, 1, _SIZE_R, 11, 13, popbyhpv2, 6);

    // popbyhpv1(:, : , : , 1) = reshape(sum(sum(sum(sum(sum(popn(1:h1, : , : , : , 1 : s, 1 : r, 1, : ), 2), 3), 4), 7), 8), h1, s, r, 1);
    // popbyhpv1(:, : , : , 2) = reshape(sum(sum(sum(sum(sum(popn(1:h1, : , : , : , 1 : s, 1 : r, 2 : 3, : ), 2), 3), 4), 7), 8), h1, s, r, 1);
    // popbyhpv1(:, : , : , 3) = reshape(sum(sum(sum(sum(sum(popn(1:h1, : , : , : , 1 : s, 1 : r, 4 : 5, : ), 2), 3), 4), 7), 8), h1, s, r, 1);
    // popbyhpv1(:, : , : , 4) = reshape(sum(sum(sum(sum(sum(popn(1:h1, : , : , : , 1 : s, 1 : r, 6 : 8, : ), 2), 3), 4), 7), 8), h1, s, r, 1);
    // popbyhpv1(:, : , : , 5) = reshape(sum(sum(sum(sum(sum(popn(1:h1, : , : , : , 1 : s, 1 : r, 9 : 10, : ), 2), 3), 4), 7), 8), h1, s, r, 1);
    // popbyhpv1(:, : , : , 6) = reshape(sum(sum(sum(sum(sum(popn(1:h1, : , : , : , 1 : s, 1 : r, 11 : 13, : ), 2), 3), 4), 7), 8), h1, s, r, 1);

    sum_2348_h1sr_a(spare_h3iav, popn, 1, _SIZE_H1, 1, _SIZE_S, 1, _SIZE_R, 1, 1, popbyhpv1, 1);
    sum_2348_h1sr_a(spare_h3iav, popn, 1, _SIZE_H1, 1, _SIZE_S, 1, _SIZE_R, 2, 3, popbyhpv1, 2);
    sum_2348_h1sr_a(spare_h3iav, popn, 1, _SIZE_H1, 1, _SIZE_S, 1, _SIZE_R, 4, 5, popbyhpv1, 3);
    sum_2348_h1sr_a(spare_h3iav, popn, 1, _SIZE_H1, 1, _SIZE_S, 1, _SIZE_R, 6, 8, popbyhpv1, 4);
    sum_2348_h1sr_a(spare_h3iav, popn, 1, _SIZE_H1, 1, _SIZE_S, 1, _SIZE_R, 9, 10, popbyhpv1, 5);
    sum_2348_h1sr_a(spare_h3iav, popn, 1, _SIZE_H1, 1, _SIZE_S, 1, _SIZE_R, 11, 13, popbyhpv1, 6);

    //      % if any
    //      % (any(any((popbyhpv3(:, 2, : , : )) < 0) == 1) == 1) == 1
    //      % disp(n)
    //      % keyboard
    //      % end

    //        %% INDICATOR THAT THINGS ARE GOING WRONG
    //          if isnan(popbysexriskage(1, 1, 2))
    //            % disp(n)
    //            % keyboard
    //            % disp('too large timestep')
    //            errmessage(1, 1) = -666;
    //            [popbyrisk_f, popbyrisk_m, , ...
    //            art_15plus, hiv_fage, hiv_mage, art_age, popbyhivhpvsexage, psize, psize_age, ...
    //            hivd_f, hivd_m, ccd_f, othd_f, othd_m, hpv1cum_pop, hpv2cum_pop, hpv3cum_pop, ...
    //            hpv1hncum_pop, hpv2hncum_pop, hpv3hncum_pop, hpv1sus_pop, hpv2sus_pop, hpv3sus_pop, ...
    //            hpv1hnsus_pop, hpv2hnsus_pop, hpv3hnsus_pop, hivsus_pop, ccinc_f, hivinc_pop_cum] = deal(0);

    //            break
    //          end

        // popbysexriskage where s=1, r=1, a=2, and all others are 1.
    if (isnan(GET_3D_SRA6(popbysexriskage, 1, 1, 2))) {
      R->errmessage[1] = -666;
      // Will come back to all those zeroes later.
      // printf("IS_NAN OCCURRED at time %d\n", n);
      break;
    }

    // this tests whether HIV among 15-49 yos is within plausible range, otherwise model run terminated
    //if (n==tstep * 45 +1 && hiv_15to49[11][1] < 0.05) {
    //	if (hiv_15to49[11][1] < 0.05 || hiv_15to49[11][1] > 0.13 || hiv_15to49[11][2] < 0.05 || hiv_15to49[11][2] > 0.13) {
    //  R->errmessage[2] = -333;
    //  break;
    // }
   // }

    // EARLY BREAK NOW DISABLED!
  //  if (n==tstep * 45 + 1 && (R->hiv_15to49[11][1] < 0.025 || R->hiv_15to49[11][1] > 0.14 || R->hiv_15to49[11][2] < 0.025 || R->hiv_15to49[11][2]  > 0.1)) {
     // R->errmessage[3] = -333;
  //    R->errmessage[2] = -111;

  //    sum_123_4678_ia(spare, popn_sum_h1_h2_h3, 3, 4, 2, 8, sum1_s);
  //    sum_123_4678_ia(spare, popn_sum_h1_h2_h3, 1, 4, 2, 8, sum2_s);

  //    R-> errmessage[1] = sum1_s[1] / sum2_s[1];
  //    R-> errmessage[2] = sum1_s[2] / sum2_s[2];

   //   printf(test);

   //   break;
  //  }


//////////////////////////////// early break in 2021 for calibration //////////////////////////////////////////////////////////
    // THIS BREAKS THE RUN IN 2021
   if (n==tstep * 72) {
      R->errmessage[3] = -333;
      R->errmessage[2] = -999;

    break;
  }


    //    %% Partnerships, mixing& balancing

    //      % [popbysexriskage] by sex(rows), risk(columns) and age(3d);


    //    %% commercial partners
    //    pships_com = popbysexriskage .* pcr_comm(:, :, [1 2 4 6 9 11]);
    //    pcr_comm_bal = pcr_comm(:, :, [1 2 4 6 9 11]);
    //    % total partnerships by sex, age, and risk
    //    pships_main = popbysexriskage .* pcr(:, :, [1 2 4 6 9 11]);

    LOOP_S LOOP_R {
      int index_sra = GET_3D_SRA_INDEX(_ds, _dr, 1) - 1;
      int index_sra6 = GET_3D_SRA6_INDEX(_ds, _dr, 1);

      pships_com[index_sra6]     = popbysexriskage[index_sra6] * S->pcr_comm[index_sra + 1];
      pships_main[index_sra6]    = popbysexriskage[index_sra6] * S->pcr[index_sra + 1];
      pcr_comm_bal[index_sra6++] = S->pcr_comm[index_sra + 1];

      pships_com[index_sra6]     = popbysexriskage[index_sra6] * S->pcr_comm[index_sra + 2];
      pships_main[index_sra6]    = popbysexriskage[index_sra6] * S->pcr[index_sra + 2];
      pcr_comm_bal[index_sra6++] = S->pcr_comm[index_sra + 2];

      pships_com[index_sra6]     = popbysexriskage[index_sra6] * S->pcr_comm[index_sra + 4];
      pships_main[index_sra6]    = popbysexriskage[index_sra6] * S->pcr[index_sra + 4];
      pcr_comm_bal[index_sra6++] = S->pcr_comm[index_sra + 4];

      pships_com[index_sra6]     = popbysexriskage[index_sra6] * S->pcr_comm[index_sra + 6];
      pships_main[index_sra6]    = popbysexriskage[index_sra6] * S->pcr[index_sra + 6];
      pcr_comm_bal[index_sra6++] = S->pcr_comm[index_sra + 6];

      pships_com[index_sra6]     = popbysexriskage[index_sra6] * S->pcr_comm[index_sra + 9];
      pships_main[index_sra6]    = popbysexriskage[index_sra6] * S->pcr[index_sra + 9];
      pcr_comm_bal[index_sra6++] = S->pcr_comm[index_sra + 9];

      pships_com[index_sra6]   = popbysexriskage[index_sra6] * S->pcr_comm[index_sra + 11];
      pships_main[index_sra6]  = popbysexriskage[index_sra6] * S->pcr[index_sra + 11];
      pcr_comm_bal[index_sra6] = S->pcr_comm[index_sra + 11];
    }

    //    pships_com_bysex = sum(sum(pships_com, 2), 3);

    LOOP_S {
      pships_com_bysex[_ds] = 0;
      LOOP_R LOOP_A6
        pships_com_bysex[_ds] += GET_3D_SRA6(pships_com, _ds, _dr, _da);
    }

    //    % no clients in the youngest age group(assumption), so only age groups
    //    % 2 + are considered c = cN / N
    //    pcr_client_men = pships_com_bysex(1)/sum(popbysexriskage(2,3,2:end)); % starts from 15+



    double pcr_client_men = pships_com_bysex[1] / (GET_3D_SRA6(popbysexriskage, 2, 3, 2) +
                                                   GET_3D_SRA6(popbysexriskage, 2, 3, 3) +
                                                   GET_3D_SRA6(popbysexriskage, 2, 3, 4) +
                                                   GET_3D_SRA6(popbysexriskage, 2, 3, 5) +
                                                   GET_3D_SRA6(popbysexriskage, 2, 3, 6));

    // pcr_comm_bal(2, 3, 2:length(eaf(1, :))) = pcr_client_men; % 6 hardcoded

    int index_232 = GET_3D_SRA6_INDEX(2, 3, 2) - 2;
    for (int _da = 2; _da <= _SIZE_A6; _da++) pcr_comm_bal[index_232 + _da] = pcr_client_men;

    //    pships_com_bal = popbysexriskage.*pcr_comm_bal;

    //      %% main partners
    //        % pships_risk = sum(pships, 3). / sum(sum(pships, 3), 2);


    LOOP_LINEAR_3D_SRA6
      GET_LIN_ALL(pships_com_bal) = GET_LIN_ALL(popbysexriskage) * GET_LIN_ALL(pcr_comm_bal);

    //    pships_com_bysex_bal = sum(sum(pships_com_bal, 2), 3);

    LOOP_S {
      pships_com_bysex_bal[_ds] = 0;
      LOOP_R LOOP_A6
        pships_com_bysex_bal[_ds] += GET_3D_SRA6(pships_com_bal, _ds, _dr, _da);
    }

    //    pshipsprop_age_com = pships_com_bal . / pships_com_bysex_bal;
    //    bal_check = pships_com_bysex_bal(1) / pships_com_bysex_bal(2);

    LOOP_S LOOP_R LOOP_A6
      GET_3D_SRA6(pshipsprop_age_com, _ds, _dr, _da) = GET_3D_SRA6(pships_com_bal, _ds, _dr, _da) / pships_com_bysex_bal[_ds];

    //      double bal_check = pships_com_bysex_bal[1] / pships_com_bysex_bal[2]
    //      Since the check below is commented, let's leave this out.

    //    % if ~bal_check
    //      % disp('commercial partner trouble');
    //    % keyboard;
    //    % end

    //      % total partnership by sexand age
    //        pships_agesum = sum(pships_main, 2);

    LOOP_S LOOP_A6 {
      pships_age_sum[_ds][_da] = GET_3D_SRA6(pships_main, _ds, 1, _da);
      for (int _dr = 2; _dr <= _SIZE_R; _dr++)
        pships_age_sum[_ds][_da] += GET_3D_SRA6(pships_main, _ds, _dr, _da);
    }

    //      % proportional partnerships in sexand age(proportional to risk)
    //        pships_ageriskprop = pships_main . / pships_agesum;

    LOOP_S LOOP_R LOOP_A6
      GET_3D_SRA6(pships_ageriskprop, _ds, _dr, _da) = GET_3D_SRA6(pships_main, _ds, _dr, _da) / pships_age_sum[_ds][_da];

    //      % proportional partnerships available from women and men
    //        pshipsprop_age_f = reshape(pships_ageriskprop(1, :, : ), 1, r * length(eaf(1, :)));
    //        pshipsprop_age_m = reshape(pships_ageriskprop(2, :, : ), 1, r * length(eaf(1, :)));
    // - length is _SIZE_R * _SIZE_A

    LOOP_A6 LOOP_R {
      int index = ((_da - 1) * _SIZE_R) + _dr;
      pshipsprop_age_f[index] = GET_3D_SRA6(pships_ageriskprop, 1, _dr, _da);
      pshipsprop_age_m[index] = GET_3D_SRA6(pships_ageriskprop, 2, _dr, _da);
    }

    //    % this replicates the mixing matrix to be the size of a* r
    //    rhof_at = repelem(eaf, r, r);
    //    rhom_at = repelem(eam, r, r);

    // Example for repelem:
    // Size of r is 3 - so rhof_at(1,1) to rhof_at(3,3) are all the same, and equal to eaf(1,1)
    //                     rhof_at(4,4) to rhof_at(6,6) are all the same, and equal to eaf(2,2)

    for (int _y = 1; _y <= 6; _y++) {
      int base_y = ((_y - 1) * _SIZE_R);
      for (int _x = 1; _x <= 6; _x++) {
        int base_x = ((_x - 1) * _SIZE_R);

        for (int _yy = 1; _yy <= _SIZE_R; _yy++)
          for (int _xx = 1; _xx <= _SIZE_R; _xx++) {
            int index = GET_2D_RA6RA6_INDEX(base_y + _yy, base_x + _xx);
            rhof_at[index] = eaf[_y][_x];
            rhom_at[index] = eam[_y][_x];
          }
      }
    }

    //    repmat repeats the whole matrix.
    //    So when ed is size (3,3),
    //    ed_at[1,2,3][1] = ed_at[4,5,6][1] = ed[1,2,3][1] on both axes.

    //    ed_at = repmat(ed, length(eaf(1,:)), length(eaf(1,:)));
    //
    //    length(eaf(1,:))) = 6 (A6)

    for (int _y = 1; _y <= _SIZE_A6 * 3; _y++)
      for (int _x = 1; _x <= _SIZE_A6 * 3; _x++)
        ed_at[_y][_x] = ed[1 + ((_y - 1) % 3)][1 + ((_x - 1) % 3)];

    //    rhof = rhof_at .* (e.*ed_at + (1 - e).*pshipsprop_age_m);
    //    rhom = rhom_at .* (e.*ed_at + (1 - e).*pshipsprop_age_f);

    //    % rhof = rhof_at.*pshipsprop_age_m;
    //    % rhom = rhom_at.*pshipsprop_age_f;

    for (int _y = 1; _y <= _SIZE_A6 * _SIZE_R; _y++)
      for (int _x = 1; _x <= _SIZE_A6 * _SIZE_R; _x++) {
        int index = GET_2D_RA6RA6_INDEX(_y, _x);
        rhof[index] = rhof_at[index] * (e * ed_at[_y][_x] + (1 - e) * pshipsprop_age_m[_x]);
        rhom[index] = rhom_at[index] * (e * ed_at[_y][_x] + (1 - e) * pshipsprop_age_f[_x]);
      }

    //    % balance men to women
    //      womsize = reshape(popbysexriskage(1,:,:), r*length(eaf(1,:)) ,1) ;
    //      mensize = reshape(popbysexriskage(2,:,:), r*length(eaf(1,:)) ,1) ;
    //      pcr_wom = reshape(pcr(1,:,[1 2 4 6 9 11]), r*length(eaf(1,:)), 1);
    //      pcr_men = reshape(pcr(2,:,[1 2 4 6 9 11]), r*length(eaf(1,:)), 1);

    LOOP_R {
      LOOP_A6 {
        womsize[((_da - 1) * _SIZE_R) + _dr] = GET_3D_SRA6(popbysexriskage, 1, _dr, _da);
        mensize[((_da - 1) * _SIZE_R) + _dr] = GET_3D_SRA6(popbysexriskage, 2, _dr, _da);
      }
      pcr_wom[((1 - 1) * _SIZE_R) + _dr] = GET_3D_SRA(S->pcr, 1, _dr, 1);
      pcr_men[((1 - 1) * _SIZE_R) + _dr] = GET_3D_SRA(S->pcr, 2, _dr, 1);
      pcr_wom[((2 - 1) * _SIZE_R) + _dr] = GET_3D_SRA(S->pcr, 1, _dr, 2);
      pcr_men[((2 - 1) * _SIZE_R) + _dr] = GET_3D_SRA(S->pcr, 2, _dr, 2);
      pcr_wom[((3 - 1) * _SIZE_R) + _dr] = GET_3D_SRA(S->pcr, 1, _dr, 4);
      pcr_men[((3 - 1) * _SIZE_R) + _dr] = GET_3D_SRA(S->pcr, 2, _dr, 4);
      pcr_wom[((4 - 1) * _SIZE_R) + _dr] = GET_3D_SRA(S->pcr, 1, _dr, 6);
      pcr_men[((4 - 1) * _SIZE_R) + _dr] = GET_3D_SRA(S->pcr, 2, _dr, 6);
      pcr_wom[((5 - 1) * _SIZE_R) + _dr] = GET_3D_SRA(S->pcr, 1, _dr, 9);
      pcr_men[((5 - 1) * _SIZE_R) + _dr] = GET_3D_SRA(S->pcr, 2, _dr, 9);
      pcr_wom[((6 - 1) * _SIZE_R) + _dr] = GET_3D_SRA(S->pcr, 1, _dr, 11);
      pcr_men[((6 - 1) * _SIZE_R) + _dr] = GET_3D_SRA(S->pcr, 2, _dr, 11);
    }


    //    % final balance to follow the rule(cpN = cpN)
    //      cpn_f = (pcr_wom.*rhof.*womsize);
    //      cpn_m = (pcr_men.*rhom.*mensize);

    for (int _da = 1; _da <= _SIZE_R * _SIZE_A6; _da++)
      for (int _dr = 1; _dr <= _SIZE_R * _SIZE_A6; _dr++) {
        int index = GET_2D_RA6RA6_INDEX(_da,_dr);
        cpn_f[index] = rhof[index] * pcr_wom[_da] * womsize[_da];
        cpn_m[index] = rhom[index] * pcr_men[_da] * mensize[_da];
      }

    //    for q = 1:length(cpn_f(:, 1))
    //      for w = 1 : length(cpn_f(1, :))
    //        bal(w, q) = cpn_f(q, w) / cpn_m(w, q);
    //      end
    //    end
    //    bal(isnan(bal)) = 0;

    for (int q = 1; q <= _SIZE_R * _SIZE_A6; q++)
      for (int w = 1; w <= _SIZE_R * _SIZE_A6; w++) {
        int index_wq = GET_2D_RA6RA6_INDEX(w, q);
        int index_qw = GET_2D_RA6RA6_INDEX(q, w);
        bal[index_wq] = cpn_f[index_qw] / cpn_m[index_wq];
        if (isnan(bal[index_wq])) bal[index_wq] = 0;
      }


    for (int q = 1; q <= _SIZE_R * _SIZE_A6; q++)
      for (int w = 1; w <= _SIZE_R * _SIZE_A6; w++)
        bal_f[GET_2D_RA6RA6_INDEX(w,q)] = bal[GET_2D_RA6RA6_INDEX(q,w)];

    //    c_m = pcr_men.*sqrt(bal);
    //    c_m(isinf(c_m)) = 0;
    //    c_f = pcr_wom.*sqrt(bal_f. ^ -(1 - 0.5);
    //    c_f(isinf(c_f)) = 0;

    for (int _d1 = 1; _d1 <= _SIZE_R * _SIZE_A6; _d1++)
      for (int _d2 = 1; _d2 <= _SIZE_R * _SIZE_A6; _d2++) {
        int index = GET_2D_RA6RA6_INDEX(_d1, _d2);
        c_m[index] = pcr_men[_d1] * sqrt(bal[index]);
        if (isinf(c_m[index])) c_m[index] = 0;
        c_f[index] = pcr_wom[_d1] * (1.0 / sqrt(bal_f[index]));
        if (isinf(c_f[index])) c_f[index] = 0;
      }


    //    % below is to check
    //    % pn_f = (rhof.*womsize);
    //    % cpn_m = (pcr_men.*rhom.*mensize);
    //    %
    //      % for q = 1:length(cpn_m(:, 1))
    //      % for w = 1 : length(cpn_m(1, :))
    //      %
    //      %c_f(w, q) = cpn_m(q, w) / pn_f(w, q);
    //    %
    //      % end
    //      % end
    //      % c_f(isnan(c_f)) = 0;
    //    % c_f(isinf(c_f)) = 0;
    //    %
    //      % c_m = pcr_men;
    //    % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //      % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //      % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //      % check
    //      % newpshipsf = womsize.*c_f.*rhof;
    //    % newpshipsm = mensize.*c_m.*rhom;
    //    %%
    //      %%
    //      %%
    //      % for q = 1:length(c_f(:, 1))
    //      % for w = 1 : length(c_f(1, :))
    //      %
    //      %test(w, q) = newpshipsf(q, w) . / newpshipsm(w, q);
    //    %
    //      % end
    //      % end
    //      %
    //      %
    //      %
    //      % this produces a matrix of ones when balancing works :
    //    % newpshipsf . / newpshipsm
    //      %
    //      %%used for checking
    //      % if any(any(newpshipsf . / newpshipsm)~= 1) == 1
    //      % disp("problem with balancing")
    //      % keyboard
    //      %
    //      %end

    //      % %Time - varying-- > scale up of ART by reducing stopping ART in a step - wise fashion

    //      if n < tstep * (2017.5 - dp0(22)) + 1

    // if (n < tstep * (2017.5 - p->dp0[31]) + 1) {
    /////////////////////////////////////////// HIV sens
    if (n < tstep * (2500 - p->dp0[31]) + 1) {

      //        omikron(art_id == 1 & sexf_id == 1 & age_id == 1) = hp0(15) * hp0(34); % stopping ART in women
      //        omikron(art_id == 1 & sexf_id == 1 & age_id == 2) = hp0(15) * hp0(35); % stopping ART in women
      //        omikron(art_id == 1 & sexf_id == 1 & age_id == 3) = hp0(15) * hp0(36); % stopping ART in women
      //        omikron(art_id == 1 & sexf_id == 1 & age_id == 4) = hp0(15) * hp0(37); % stopping ART in women
      //        omikron(art_id == 1 & sexf_id == 1 & age_id == 5) = hp0(15) * hp0(38); % stopping ART in women

      LOOP_LINEAR{
        if IS_LIN(S->art_id, 1) {
          if IS_LIN(S->sexf_id, 1) {

            //        omikron(art_id == 1 & sexm_id == 1 & age_id == 1) = hp0(15) * hp0(16) * hp0(34); % stopping ART in men(faster)
            //        omikron(art_id == 1 & sexm_id == 1 & age_id == 2) = hp0(15) * hp0(16) * hp0(35); % stopping ART in men(faster)
            //        omikron(art_id == 1 & sexm_id == 1 & age_id == 3) = hp0(15) * hp0(16) * hp0(36); % stopping ART in men(faster)
            //        omikron(art_id == 1 & sexm_id == 1 & age_id == 4) = hp0(15) * hp0(16) * hp0(37); % stopping ART in men(faster)
            //        omikron(art_id == 1 & sexm_id == 1 & age_id == 5) = hp0(15) * hp0(16) * hp0(38); % stopping ART in men(faster)

            if IS_LIN(S->age_id, 1) GET_LIN_ALL(S->omikron) = p->hp0[15] * p->hp0[34];
            else if IS_LIN(S->age_id, 2) GET_LIN_ALL(S->omikron) = p->hp0[15] * p->hp0[35];
            else if IS_LIN(S->age_id, 3) GET_LIN_ALL(S->omikron) = p->hp0[15] * p->hp0[36];
            else if IS_LIN(S->age_id, 4) GET_LIN_ALL(S->omikron) = p->hp0[15] * p->hp0[37];
            else if IS_LIN(S->age_id, 5) GET_LIN_ALL(S->omikron) = p->hp0[15] * p->hp0[38];

                        //        omikron(art_id == 1 & sexf_id == 1 & risk_id == 3) = hp0(17); % stopping ART in FSW // not related to age!

            if IS_LIN(S->risk_id, 3) GET_LIN_ALL(S->omikron) = p->hp0[17];
          }

          if IS_LIN(S->sexm_id, 1) {
            if IS_LIN(S->age_id,  1) GET_LIN_ALL(S->omikron) = _h15_16 * p->hp0[34];
            else if IS_LIN(S->age_id,  2) GET_LIN_ALL(S->omikron) = _h15_16 * p->hp0[35];
            else if IS_LIN(S->age_id,  3) GET_LIN_ALL(S->omikron) = _h15_16 * p->hp0[36];
            else if IS_LIN(S->age_id,  4) GET_LIN_ALL(S->omikron) = _h15_16 * p->hp0[37];
            else if IS_LIN(S->age_id,  5) GET_LIN_ALL(S->omikron) = _h15_16 * p->hp0[38];
          }
        }
      }
    }

    //if n > tstep* (2017.5 - dp0(22))& n < tstep * (2018 - dp0(22)) + 1

    /////////////////////////////////////////// HIV sens
    if ((n > tstep* (2017.5 - p->dp0[31])) && (n < tstep * (2018 - p->dp0[31]) + 1)) {

      LOOP_LINEAR {
        if IS_LIN(S->art_id, 1) {

  //          if (IS(S->sexf_id, 1) || IS(S->sexm_id, 1)) {

  //      omikron(art_id == 1 & sexf_id == 1 & age_id == 1) = hp0(15) * hp0(34); % stopping ART in women
  //      omikron(art_id == 1 & sexf_id == 1 & age_id == 2) = hp0(15) * hp0(35); % stopping ART in women
  //      omikron(art_id == 1 & sexf_id == 1 & age_id == 3) = hp0(15) * hp0(36); % stopping ART in women
  //      omikron(art_id == 1 & sexf_id == 1 & age_id == 4) = hp0(15) * hp0(37); % stopping ART in women
  //      omikron(art_id == 1 & sexf_id == 1 & age_id == 5) = hp0(15) * hp0(38); % stopping ART in women
  //      omikron(art_id == 1 & sexm_id == 1 & age_id == 1) = hp0(15) * hp0(34); % stopping ART in men(faster)
  //      omikron(art_id == 1 & sexm_id == 1 & age_id == 2) = hp0(15) * hp0(35); % stopping ART in men(faster)
  //      omikron(art_id == 1 & sexm_id == 1 & age_id == 3) = hp0(15) * hp0(36); % stopping ART in men(faster)
  //      omikron(art_id == 1 & sexm_id == 1 & age_id == 4) = hp0(15) * hp0(37); % stopping ART in men(faster)
  //      omikron(art_id == 1 & sexm_id == 1 & age_id == 5) = hp0(15) * hp0(38); % stopping ART in men(faster)

          if IS_LIN(S->age_id, 1) GET_LIN_ALL(S->omikron) = p->hp0[15] * p->hp0[34];
          else if IS_LIN(S->age_id, 2) GET_LIN_ALL(S->omikron) = p->hp0[15] * p->hp0[35];
          else if IS_LIN(S->age_id, 3) GET_LIN_ALL(S->omikron) = p->hp0[15] * p->hp0[36];
          else if IS_LIN(S->age_id, 4) GET_LIN_ALL(S->omikron) = p->hp0[15] * p->hp0[37];
          else if IS_LIN(S->age_id, 5) GET_LIN_ALL(S->omikron) = p->hp0[15] * p->hp0[38];
          //      }

          //      omikron(art_id == 1 & sexf_id == 1 & risk_id == 3) = hp0(17); % stopping ART in FSW // not related to age!

          if IS_LIN(S->sexf_id, 1) {
            if IS_LIN(S->risk_id, 3) GET_LIN_ALL(S->omikron) = p->hp0[17];
          }
        }
      }
    }

    //    if n > tstep* (2018 - dp0(22))
    /////////////////////////////////////////// HIV sens
    if (n > tstep* (2018 - p->dp0[31])) {

      //      omikron(art_id == 1 & sexf_id == 1 & age_id == 1) = hp0(15); % stopping ART in women
      //      omikron(art_id == 1 & sexf_id == 1 & age_id == 2) = hp0(15); % stopping ART in women
      //      omikron(art_id == 1 & sexf_id == 1 & age_id == 3) = hp0(15); % stopping ART in women
      //      omikron(art_id == 1 & sexf_id == 1 & age_id == 4) = hp0(15); % stopping ART in women
      //      omikron(art_id == 1 & sexf_id == 1 & age_id == 5) = hp0(15); % stopping ART in women
      //      omikron(art_id == 1 & sexm_id == 1 & age_id == 1) = hp0(15); % stopping ART in men
      //      omikron(art_id == 1 & sexm_id == 1 & age_id == 2) = hp0(15); % stopping ART in men
      //      omikron(art_id == 1 & sexm_id == 1 & age_id == 3) = hp0(15); % stopping ART in men
      //      omikron(art_id == 1 & sexm_id == 1 & age_id == 4) = hp0(15); % stopping ART in men
      //      omikron(art_id == 1 & sexm_id == 1 & age_id == 5) = hp0(15); % stopping ART in men

      LOOP_LINEAR{
        if IS_LIN(S->art_id, 1) {
          if (IS_LIN(S->sexf_id, 1) || IS_LIN(S->sexm_id, 1)) {
            // Does this mean all 5 age groups? DO we need the ifs?
                 if IS_LIN(S->age_id,  1) GET_LIN_ALL(S->omikron) = p->hp0[15];
            else if IS_LIN(S->age_id,  2) GET_LIN_ALL(S->omikron) = p->hp0[15];
            else if IS_LIN(S->age_id,  3) GET_LIN_ALL(S->omikron) = p->hp0[15];
            else if IS_LIN(S->age_id,  4) GET_LIN_ALL(S->omikron) = p->hp0[15];
            else if IS_LIN(S->age_id,  5) GET_LIN_ALL(S->omikron) = p->hp0[15];
          }
        }

      //      omikron(art_id == 1 & sexf_id == 1 & risk_id == 3) = hp0(17); % stopping ART in FSW // not related to age!

        if (IS_LIN(S->sexf_id, 1) && IS_LIN(S->risk_id, 3)) {
          GET_LIN_ALL(S->omikron) = p->hp0[17];
        }
      }
        //    end
    }


    //    %% Time - varying ART impact on HPV transmission
    //      % ART has a larger impact(after 2010) on clearance, immunity, progression to
    //      % CIN2, and regression from CIN2 by reducing the impact of the HPV cofactor

    //      clrfastart = (pp0(5) + pp0(5) * tau_coff(n));% works when param between 0 - 1
    //      immslwart = 1 - (1 - tau_coff(n)) + pp0(16) * (1 - tau_coff(n));% works when param > 1
    //      regslwart = (pp0(23) + pp0(23) * tau_coff(n));
    //      prgfastart = 1 - (1 - tau_coff(n)) + pp0(48) * (1 - tau_coff(n));
    //      betaonart = 1 - (1 - tau_coff(n)) + pp0(4) * (1 - tau_coff(n));

    double clrfastart = p->pp0[5] + p->pp0[5] * tau_coff[n];
    double immslwart = 1.0 - (1.0 - tau_coff[n]) + p->pp0[16] * (1.0 - tau_coff[n]);
    double regslwart = p->pp0[23] + p->pp0[23] * tau_coff[n];
    double prgfastart = 1.0 - (1.0 - tau_coff[n]) + p->pp0[48] * (1.0 - tau_coff[n]);
    double betaonart = 1.0 - (1.0 - tau_coff[n]) + p->pp0[4] * (1.0 - tau_coff[n]);

    //    %% CC screening treatment success depending on baseline / intervention scenario
    //      if screenscen.setup == 0

            //txsuxx_cinhn = pp0(64);
    //        txsuxx_cinhp = pp0(65);

    //    txsuxx_cchn = pp0(66);
    //    txsuxx_cchp = pp0(67);
    //    end


    //      if (screenscen.setup == 1 | screenscen.setup == 2)
    //        if n < tstep * (2020 - dp0(22)) + 1
    //          txsuxx_cinhn = pp0(64);
    //          txsuxx_cinhp = pp0(65);
    //          txsuxx_cchn = pp0(66);
    //          txsuxx_cchp = pp0(67);
    //        end
    //        if n > tstep* (2020 - dp0(22))
    //          txsuxx_cinhn = screenscen.cintx(1);
    //          txsuxx_cinhp = screenscen.cintx(2);
    //          txsuxx_cchn = screenscen.cctx(1);
    //          txsuxx_cchp = screenscen.cctx(2);
    //        end
    //      end

    // DO THIS ? ((int)round(p->vaccscen_setup) > 0)

   // OLD if (p->screenscen_setup == 0) {
	if ((int)round(p->screenscen_setup) == 0) {
      p->txsuxx_cinhn = p->pp0[64];
      p->txsuxx_cinhp = p->pp0[65];
      p->txsuxx_cchn = p->pp0[66];
      p->txsuxx_cchp = p->pp0[67];

    }
    else if (((int)round(p->screenscen_setup) == 1) || ((int)round(p->screenscen_setup) == 2)) {
      if (n < tstep * (2020 - p->dp0[31]) + 1) {
        p->txsuxx_cinhn = p->pp0[64];
        p->txsuxx_cinhp = p->pp0[65];
        p->txsuxx_cchn = p->pp0[66];
        p->txsuxx_cchp = p->pp0[67];
      }
      if (n > tstep* (2020 - p->dp0[31])) {
        p->txsuxx_cinhn = p->screenscen_cintx[1];
        p->txsuxx_cinhp = p->screenscen_cintx[2];
        p->txsuxx_cchn = p->screenscen_cctx[1];
        p->txsuxx_cchp = p->screenscen_cctx[2];

       // double temp1 = p->screenscen_cintx[1];
       // double temp2 = p->screenscen_cintx[2];
       // double temp3 = p->screenscen_cctx[1];
       // double temp4 = p->screenscen_cctx[1];

       // printf("screenscen_cintx[1]", temp1);
       // printf("screenscen_cintx[2]", temp2);
       // printf("screenscen_cctx[1]", temp3);
       // printf("screenscen_cctx[2]", temp4);


      }
    }

    // The below is probably over-optimised...

    double _p10_clr = p->pp0[10] * clrfastart;
    double _p10_11_clr = _p10_clr * p->pp0[11];
    double _p10_12_clr = _p10_clr * p->pp0[12];
    double _p10_32_clr = _p10_clr * p->pp0[32];
    double _p10_32_11_clr = _p10_32_clr * p->pp0[11];
    double _p10_32_12_clr = _p10_32_clr * p->pp0[12];
    double _p13_imm = p->pp0[13] * immslwart;
    double _p13_14_imm = _p13_imm * p->pp0[14];
    double _p13_15_imm = _p13_imm * p->pp0[15];
    double _p20_regs = p->pp0[20] * regslwart;
    double _p20_p21_regs = p->pp0[20] * p->pp0[21] * regslwart;
    double _p20_p22_regs = p->pp0[20] * p->pp0[22] * regslwart;
    double _p20_regs_o24 = _p20_regs / p->pp0[24];
    double _p20_p21_regs_o24 = p->pp0[20] * p->pp0[21] * regslwart / p->pp0[24];
    double _p20_p22_regs_o24 = p->pp0[20] * p->pp0[22] * regslwart / p->pp0[24];
    double _p17_prg = p->pp0[17] * prgfastart;
    double _p17_18_prg = _p17_prg * p->pp0[18];
    double _p17_19_prg = _p17_prg * p->pp0[19];
    double _p1_beta = p->pp0[1] * betaonart;
    double _p1_2_beta = _p1_beta * p->pp0[2];
    double _p1_3_beta = _p1_beta * p->pp0[3];
    double _p_1_mcirc = p->hp0[1] * mcirc[n];
    double _p_1_5_mcirc = _p_1_mcirc * p->hp0[5];
    double _p_1_m6_mcirc = _p_1_mcirc * (1.0 - p->hp0[6]);
    double _p_3_mcirc = p->hp0[3] * mcirc[n];
    double _p_3_5_mcirc = _p_3_mcirc * p->hp0[5];
    double _p_3_m6_mcirc = _p_3_mcirc * (1.0 - p->hp0[6]);
    double _p_1_7_mcirc = _p_1_mcirc * p->hp0[7];
    double _p_1_5_7_mcirc = _p_1_5_mcirc * p->hp0[7];
    double _p_1_7_m6_mcirc = _p_1_7_mcirc * (1.0 - p->hp0[6]);
    double _th25_cin_29n = theta25to29hn[n] * p->txsuxx_cinhn;
    double _th25_cin_29p = theta25to29hp[n] * p->txsuxx_cinhp;
    double _th30_cin_34p = theta30to34hp[n] * p->txsuxx_cinhp;
    double _th35_cin_39n = theta35to39hn[n] * p->txsuxx_cinhn;
    double _th35_cin_39p = theta35to39hp[n] * p->txsuxx_cinhp;
    double _th40_cin_44p = theta40to44hp[n] * p->txsuxx_cinhp;
    double _th45_cin_49n = theta45to49hn[n] * p->txsuxx_cinhn;
    double _th45_cin_49p = theta45to49hp[n] * p->txsuxx_cinhp;

    double _th25_cc_29n = theta25to29hn[n] * p->txsuxx_cchn;
    double _th25_cc_29p = theta25to29hp[n] * p->txsuxx_cchp;
    double _th30_cc_34p = theta30to34hp[n] * p->txsuxx_cchp;
    double _th35_cc_39n = theta35to39hn[n] * p->txsuxx_cchn;
    double _th35_cc_39p = theta35to39hp[n] * p->txsuxx_cchp;
    double _th40_cc_44p = theta40to44hp[n] * p->txsuxx_cchp;
    double _th45_cc_49n = theta45to49hn[n] * p->txsuxx_cchn;
    double _th45_cc_49p = theta45to49hp[n] * p->txsuxx_cchp;

    LOOP_LINEAR{

      if IS_LIN(S->art_id, 1) {

        //    sigma1(hpv1_id == 1 & art_id == 1 & sexf_id == 1) = pp0(10) * clrfastart;
        //    sigma2(hpv2_id == 1 & art_id == 1 & sexf_id == 1) = pp0(10) * pp0(11) * clrfastart;
        //    sigma3(hpv3_id == 1 & art_id == 1 & sexf_id == 1) = pp0(10) * pp0(12) * clrfastart;

        if IS_LIN(S->sexf_id, 1) {
          if IS_LIN(S->hpv1_id, 1) GET_LIN_ALL(S->sigma1) = _p10_clr;
          if IS_LIN(S->hpv2_id, 1) GET_LIN_ALL(S->sigma2) = _p10_11_clr;
          if IS_LIN(S->hpv3_id, 1) GET_LIN_ALL(S->sigma3) = _p10_12_clr;
        }

        //    sigma1(hpv1_id == 1 & art_id == 1 & sexm_id == 1) = pp0(10) * pp0(32) * clrfastart;
        //    sigma2(hpv2_id == 1 & art_id == 1 & sexm_id == 1) = pp0(10) * pp0(32) * pp0(11) * clrfastart;
        //    sigma3(hpv3_id == 1 & art_id == 1 & sexm_id == 1) = pp0(10) * pp0(32) * pp0(12) * clrfastart;

        if IS_LIN(S->sexm_id, 1) {
          if IS_LIN(S->hpv1_id, 1) GET_LIN_ALL(S->sigma1) = _p10_32_clr;
          if IS_LIN(S->hpv2_id, 1) GET_LIN_ALL(S->sigma2) = _p10_32_11_clr;
          if IS_LIN(S->hpv3_id, 1) GET_LIN_ALL(S->sigma3) = _p10_32_12_clr;
        }

        //    delta1(hpv1R_id == 1 & art_id == 1) = pp0(13) * immslwart; %
        //    delta2(hpv2R_id == 1 & art_id == 1) = pp0(13) * pp0(14) * immslwart; %
        //    delta3(hpv3R_id == 1 & art_id == 1) = pp0(13) * pp0(15) * immslwart; %

        if IS_LIN(S->hpv1R_id, 1) GET_LIN_ALL(S->delta1) = _p13_imm;
        if IS_LIN(S->hpv2R_id, 1) GET_LIN_ALL(S->delta2) = _p13_14_imm;
        if IS_LIN(S->hpv3R_id, 1) GET_LIN_ALL(S->delta3) = _p13_15_imm;

        //    omega1(hpv1CIN_id == 1 & art_id == 1 & age_id < 4) = pp0(20) * regslwart;% regression from CIN2 / CIN3
        //    omega2(hpv2CIN_id == 1 & art_id == 1 & age_id < 4) = pp0(20)*pp0(21)*regslwart; % regression from CIN2/CIN3
        //    omega3(hpv3CIN_id == 1 & art_id == 1 & age_id < 4) = pp0(20)*pp0(22)*regslwart; % regression from CIN2/CIN3
        // CHANGED HERE SO THAT REDUCTION IN CIN2+ clearance starts at 50+
        if LT_LIN(S->age_id, 4) {
          if IS_LIN(S->hpv1CIN_id, 1) GET_LIN_ALL(S->omega1) = _p20_regs;
          if IS_LIN(S->hpv2CIN_id, 1) GET_LIN_ALL(S->omega2) = _p20_p21_regs;
          if IS_LIN(S->hpv3CIN_id, 1) GET_LIN_ALL(S->omega3) = _p20_p22_regs;
        }

        //    omega1(hpv1CIN_id == 1 & art_id == 1 & age_id > 3) = pp0(20) * regslwart / pp0(24);% regression from CIN2 / CIN3
        //    omega2(hpv2CIN_id == 1 & art_id == 1 & age_id > 3) = pp0(20)*pp0(21)*regslwart/pp0(24); % regression from CIN2/CIN3
        //    omega3(hpv3CIN_id == 1 & art_id == 1 & age_id > 3) = pp0(20)*pp0(22)*regslwart/pp0(24); % regression from CIN2/CIN3
        // CHANGED HERE SO THAT REDUCTION IN CIN2+ clearance starts at 50+
        if GT_LIN(S->age_id, 3) {
          if IS_LIN(S->hpv1CIN_id, 1) GET_LIN_ALL(S->omega1) = _p20_regs_o24;
          if IS_LIN(S->hpv2CIN_id, 1) GET_LIN_ALL(S->omega2) = _p20_p21_regs_o24;
          if IS_LIN(S->hpv3CIN_id, 1) GET_LIN_ALL(S->omega3) = _p20_p22_regs_o24;
        }

        //      psi1(hpv1I_id == 1 & art_id == 1) = pp0(17) * prgfastart;% progression to CIN2 / CIN3
        //      psi2(hpv2I_id == 1 & art_id == 1) = pp0(17) * pp0(18) * prgfastart;% progression to CIN2 / CIN3
        //      psi3(hpv3I_id == 1 & art_id == 1) = pp0(17) * pp0(19) * prgfastart;% progression to CIN2 / CIN3

        if IS_LIN(S->hpv1I_id, 1) GET_LIN_ALL(S->psi1) = _p17_prg;
        if IS_LIN(S->hpv2I_id, 1) GET_LIN_ALL(S->psi2) = _p17_18_prg;
        if IS_LIN(S->hpv3I_id, 1) GET_LIN_ALL(S->psi3) = _p17_19_prg;

        //      betahpv1(hpv1S_id == 1 & art_id == 1) = pp0(1) * betaonart;% vvtHPV acquisition_act prob
        //      betahpv2(hpv2S_id == 1 & art_id == 1) = pp0(1) * pp0(2) * betaonart;% vvtHPV acquisition_act prob
        //      betahpv3(hpv3S_id == 1 & art_id == 1) = pp0(1) * pp0(3) * betaonart;% vvtHPV acquisition_act prob

        if IS_LIN(S->hpv1S_id, 1) GET_LIN_ALL(S->betahpv1) = _p1_beta;
        if IS_LIN(S->hpv2S_id, 1) GET_LIN_ALL(S->betahpv2) = _p1_2_beta;
        if IS_LIN(S->hpv3S_id, 1) GET_LIN_ALL(S->betahpv3) = _p1_3_beta;
      }

    //      % so that beta cannot go over 1 when multiplied by RR
    //      betahpv1(betahpv1 > 1) = 1;
    //      betahpv2(betahpv2 > 1) = 1;
    //      betahpv3(betahpv3 > 1) = 1;

      if (GT_LIN(S->betahpv1, 1)) GET_LIN_ALL(S->betahpv1) = 1;
      if (GT_LIN(S->betahpv2, 1)) GET_LIN_ALL(S->betahpv2) = 1;
      if (GT_LIN(S->betahpv3, 1)) GET_LIN_ALL(S->betahpv3) = 1;

      if IS_LIN(S->sexf_id, 1) {

            //    %% Screeningand treatment for CIN2 + and CC
            //    theta1(hpv1CIN_id == 1 & age_id == 2 & sexf_id == 1 & hiv_id == 0) = theta15to24hn(n) * txsuxx_cinhn;
            //    theta1(hpv1CIN_id == 1 & age_id == 2 & sexf_id == 1 & hiv_id == 1) = theta15to24hp(n) * txsuxx_cinhp;
            //    theta1(hpv1CIN_id == 1 & age_id == 3 & sexf_id == 1 & hiv_id == 0) = theta25to34hn(n) * txsuxx_cinhn;
            //    theta1(hpv1CIN_id == 1 & age_id == 3 & sexf_id == 1 & hiv_id == 1) = theta25to34hp(n) * txsuxx_cinhp;
            //    theta1(hpv1CIN_id == 1 & age_id == 4 & sexf_id == 1 & hiv_id == 0) = theta35to49hn(n) * txsuxx_cinhn;
            //    theta1(hpv1CIN_id == 1 & age_id == 4 & sexf_id == 1 & hiv_id == 1) = theta35to49hp(n) * txsuxx_cinhp;

        if IS_LIN(S->hpv1CIN_id, 1) {

          if IS_LIN(S->age_id2, 4) {
            if IS_LIN(S->hiv_id, 0) GET_LIN_ALL(S->theta1) = _th25_cin_29n;
            else if IS_LIN(S->hiv_id, 1) GET_LIN_ALL(S->theta1) = _th25_cin_29p;

          } else if IS_LIN(S->age_id2, 5) {
			if IS_LIN(S->hiv_id, 1) GET_LIN_ALL(S->theta1) = _th30_cin_34p;

          } else if IS_LIN(S->age_id2, 6) {
            if IS_LIN(S->hiv_id, 0) GET_LIN_ALL(S->theta1) = _th35_cin_39n;
            else if IS_LIN(S->hiv_id, 1) GET_LIN_ALL(S->theta1) = _th35_cin_39p;

          } else if IS_LIN(S->age_id2, 7) {
			if IS_LIN(S->hiv_id, 1) GET_LIN_ALL(S->theta1) = _th40_cin_44p;

          } else if IS_LIN(S->age_id2, 8) {
            if IS_LIN(S->hiv_id, 0) GET_LIN_ALL(S->theta1) = _th45_cin_49n;
            else if IS_LIN(S->hiv_id, 1) GET_LIN_ALL(S->theta1) = _th45_cin_49p;
          }
        }

        //    theta2(hpv2CIN_id == 1 & age_id == 2 & sexf_id == 1 & hiv_id == 0) = theta15to24hn(n) * txsuxx_cinhn;
        //    theta2(hpv2CIN_id == 1 & age_id == 2 & sexf_id == 1 & hiv_id == 1) = theta15to24hp(n) * txsuxx_cinhp;
        //    theta2(hpv2CIN_id == 1 & age_id == 3 & sexf_id == 1 & hiv_id == 0) = theta25to34hn(n) * txsuxx_cinhn;
        //    theta2(hpv2CIN_id == 1 & age_id == 3 & sexf_id == 1 & hiv_id == 1) = theta25to34hp(n) * txsuxx_cinhp;
        //    theta2(hpv2CIN_id == 1 & age_id == 4 & sexf_id == 1 & hiv_id == 0) = theta35to49hn(n) * txsuxx_cinhn;
        //    theta2(hpv2CIN_id == 1 & age_id == 4 & sexf_id == 1 & hiv_id == 1) = theta35to49hp(n) * txsuxx_cinhp;

        if IS_LIN(S->hpv2CIN_id, 1) {
          if IS_LIN(S->age_id2, 4) {
            if IS_LIN(S->hiv_id, 0) GET_LIN_ALL(S->theta2) = _th25_cin_29n;
            else if IS_LIN(S->hiv_id, 1) GET_LIN_ALL(S->theta2) = _th25_cin_29p;

          } else if IS_LIN(S->age_id2, 5) {
			if IS_LIN(S->hiv_id, 1) GET_LIN_ALL(S->theta2) = _th30_cin_34p;

          }  else if IS_LIN(S->age_id2, 6) {
            if IS_LIN(S->hiv_id, 0) GET_LIN_ALL(S->theta2) = _th35_cin_39n;
            else if IS_LIN(S->hiv_id, 1) GET_LIN_ALL(S->theta2) = _th35_cin_39p;

          } else if IS_LIN(S->age_id2, 7) {
			if IS_LIN(S->hiv_id, 1) GET_LIN_ALL(S->theta2) = _th40_cin_44p;

          }  else if IS_LIN(S->age_id2, 8) {
            if IS_LIN(S->hiv_id, 0) GET_LIN_ALL(S->theta2) = _th45_cin_49n;
            else if IS_LIN(S->hiv_id, 1) GET_LIN_ALL(S->theta2) = _th45_cin_49p;
          }
        }

        //    theta3(hpv3CIN_id == 1 & age_id == 2 & sexf_id == 1 & hiv_id == 0) = theta15to24hn(n) * txsuxx_cinhn;
        //    theta3(hpv3CIN_id == 1 & age_id == 2 & sexf_id == 1 & hiv_id == 1) = theta15to24hp(n) * txsuxx_cinhp;
        //    theta3(hpv3CIN_id == 1 & age_id == 3 & sexf_id == 1 & hiv_id == 0) = theta25to34hn(n) * txsuxx_cinhn;
        //    theta3(hpv3CIN_id == 1 & age_id == 3 & sexf_id == 1 & hiv_id == 1) = theta25to34hp(n) * txsuxx_cinhp;
        //    theta3(hpv3CIN_id == 1 & age_id == 4 & sexf_id == 1 & hiv_id == 0) = theta35to49hn(n) * txsuxx_cinhn;
        //    theta3(hpv3CIN_id == 1 & age_id == 4 & sexf_id == 1 & hiv_id == 1) = theta35to49hp(n) * txsuxx_cinhp;

        if IS_LIN(S->hpv3CIN_id, 1) {
          if IS_LIN(S->age_id2, 4) {
            if IS_LIN(S->hiv_id, 0) GET_LIN_ALL(S->theta3) = _th25_cin_29n;
            else if IS_LIN(S->hiv_id, 1) GET_LIN_ALL(S->theta3) = _th25_cin_29p;

          } else if IS_LIN(S->age_id2, 5) {
			if IS_LIN(S->hiv_id, 1) GET_LIN_ALL(S->theta3) = _th30_cin_34p;

          }  else if IS_LIN(S->age_id2, 6) {
            if IS_LIN(S->hiv_id, 0) GET_LIN_ALL(S->theta3) = _th35_cin_39n;
            else if IS_LIN(S->hiv_id, 1) GET_LIN_ALL(S->theta3) = _th35_cin_39p;

          } else if IS_LIN(S->age_id2, 7) {
			if IS_LIN(S->hiv_id, 1) GET_LIN_ALL(S->theta3) = _th40_cin_44p;

          } else if IS_LIN(S->age_id2, 8) {
            if IS_LIN(S->hiv_id, 0) GET_LIN_ALL(S->theta3) = _th45_cin_49n;
            else if IS_LIN(S->hiv_id, 1) GET_LIN_ALL(S->theta3) = _th45_cin_49p;
          }
        }

        //    hups1(hpv1CC_id == 1 & age_id == 2 & sexf_id == 1 & hiv_id == 0) = theta15to24hn(n) * txsuxx_cchn;
        //    hups1(hpv1CC_id == 1 & age_id == 2 & sexf_id == 1 & hiv_id == 1) = theta15to24hp(n) * txsuxx_cchp;
        //    hups1(hpv1CC_id == 1 & age_id == 3 & sexf_id == 1 & hiv_id == 0) = theta25to34hn(n) * txsuxx_cchn;
        //    hups1(hpv1CC_id == 1 & age_id == 3 & sexf_id == 1 & hiv_id == 1) = theta25to34hp(n) * txsuxx_cchp;
        //    hups1(hpv1CC_id == 1 & age_id == 4 & sexf_id == 1 & hiv_id == 0) = theta35to49hn(n) * txsuxx_cchn;
        //    hups1(hpv1CC_id == 1 & age_id == 4 & sexf_id == 1 & hiv_id == 1) = theta35to49hp(n) * txsuxx_cchp;

        if IS_LIN(S->hpv1CC_id, 1) {
          if IS_LIN(S->age_id2, 4) {
            if IS_LIN(S->hiv_id, 0) GET_LIN_ALL(S->hups1) = _th25_cc_29n;
            else if IS_LIN(S->hiv_id, 1) GET_LIN_ALL(S->hups1) = _th25_cc_29p;

          } else if IS_LIN(S->age_id2, 5) {
			if IS_LIN(S->hiv_id, 1) GET_LIN_ALL(S->hups1) = _th30_cc_34p;

          } else if IS_LIN(S->age_id2, 6) {
            if IS_LIN(S->hiv_id, 0) GET_LIN_ALL(S->hups1) = _th35_cc_39n;
            else if IS_LIN(S->hiv_id, 1) GET_LIN_ALL(S->hups1) = _th35_cc_39p;

          } else if IS_LIN(S->age_id2, 7) {
         			if IS_LIN(S->hiv_id, 1) GET_LIN_ALL(S->hups1) = _th40_cc_44p;

          } else if IS_LIN(S->age_id2, 8) {
            if IS_LIN(S->hiv_id, 0) GET_LIN_ALL(S->hups1) = _th45_cc_49n;
            else if IS_LIN(S->hiv_id, 1) GET_LIN_ALL(S->hups1) = _th45_cc_49p;
          }
         }

        //    hups2(hpv2CC_id == 1 & age_id == 2 & sexf_id == 1 & hiv_id == 0) = theta15to24hn(n) * txsuxx_cchn;
        //    hups2(hpv2CC_id == 1 & age_id == 2 & sexf_id == 1 & hiv_id == 1) = theta15to24hp(n) * txsuxx_cchp;
        //    hups2(hpv2CC_id == 1 & age_id == 3 & sexf_id == 1 & hiv_id == 0) = theta25to34hn(n) * txsuxx_cchn;
        //    hups2(hpv2CC_id == 1 & age_id == 3 & sexf_id == 1 & hiv_id == 1) = theta25to34hp(n) * txsuxx_cchp;
        //    hups2(hpv2CC_id == 1 & age_id == 4 & sexf_id == 1 & hiv_id == 0) = theta35to49hn(n) * txsuxx_cchn;
        //    hups2(hpv2CC_id == 1 & age_id == 4 & sexf_id == 1 & hiv_id == 1) = theta35to49hp(n) * txsuxx_cchp;

        if IS_LIN(S->hpv2CC_id, 1) {
          if IS_LIN(S->age_id2, 4) {
            if IS_LIN(S->hiv_id, 0) GET_LIN_ALL(S->hups2) = _th25_cc_29n;
            else if IS_LIN(S->hiv_id, 1) GET_LIN_ALL(S->hups2) = _th25_cc_29p;

          } else if IS_LIN(S->age_id2, 5) {
			if IS_LIN(S->hiv_id, 1) GET_LIN_ALL(S->hups2) = _th30_cc_34p;

          } else if IS_LIN(S->age_id2, 6) {
            if IS_LIN(S->hiv_id, 0) GET_LIN_ALL(S->hups2) = _th35_cc_39n;
            else if IS_LIN(S->hiv_id, 1) GET_LIN_ALL(S->hups2) = _th35_cc_39p;

          } else if IS_LIN(S->age_id2, 7) {
			if IS_LIN(S->hiv_id, 1) GET_LIN_ALL(S->hups2) = _th40_cc_44p;

          } else if IS_LIN(S->age_id2, 8) {
            if IS_LIN(S->hiv_id, 0) GET_LIN_ALL(S->hups2) = _th45_cc_49n;
            else if IS_LIN(S->hiv_id, 1) GET_LIN_ALL(S->hups2) = _th45_cc_49p;
          }
        }

        //    hups3(hpv3CC_id == 1 & age_id == 2 & sexf_id == 1 & hiv_id == 0) = theta15to24hn(n) * txsuxx_cchn;
        //    hups3(hpv3CC_id == 1 & age_id == 2 & sexf_id == 1 & hiv_id == 1) = theta15to24hp(n) * txsuxx_cchp;
        //    hups3(hpv3CC_id == 1 & age_id == 3 & sexf_id == 1 & hiv_id == 0) = theta25to34hn(n) * txsuxx_cchn;
        //    hups3(hpv3CC_id == 1 & age_id == 3 & sexf_id == 1 & hiv_id == 1) = theta25to34hp(n) * txsuxx_cchp;
        //    hups3(hpv3CC_id == 1 & age_id == 4 & sexf_id == 1 & hiv_id == 0) = theta35to49hn(n) * txsuxx_cchn;
        //    hups3(hpv3CC_id == 1 & age_id == 4 & sexf_id == 1 & hiv_id == 1) = theta35to49hp(n) * txsuxx_cchp;

        if IS_LIN(S->hpv3CC_id, 1) {
          if IS_LIN(S->age_id2, 4) {
            if IS_LIN(S->hiv_id, 0) GET_LIN_ALL(S->hups3) = _th25_cc_29n;
            else if IS_LIN(S->hiv_id, 1) GET_LIN_ALL(S->hups3) = _th25_cc_29p;

          } else if IS_LIN(S->age_id2, 5) {
			if IS_LIN(S->hiv_id, 1) GET_LIN_ALL(S->hups3) = _th30_cc_34p;

          } else if IS_LIN(S->age_id2, 6) {
            if IS_LIN(S->hiv_id, 0) GET_LIN_ALL(S->hups3) = _th35_cc_39n;
            else if IS_LIN(S->hiv_id, 1) GET_LIN_ALL(S->hups3) = _th35_cc_39p;

          } else if IS_LIN(S->age_id2, 7) {
			if IS_LIN(S->hiv_id, 1) GET_LIN_ALL(S->hups3) = _th40_cc_44p;

          } else if IS_LIN(S->age_id2, 8) {
            if IS_LIN(S->hiv_id, 0) GET_LIN_ALL(S->hups3) = _th45_cc_49n;
            else if IS_LIN(S->hiv_id, 1) GET_LIN_ALL(S->hups3) = _th45_cc_49p;
          }
        }
      }


      if IS_LIN(S->hiv2_id, 1) {
        //  % ART
        //    tau(hiv2_id == 1) = tau_tv(n);% initiation of ART; this to become time - varying
        GET_LIN_ALL(S->tau) = tau_tv[n];
      }

      if (IS_LIN(S->hivS_id, 1) && IS_LIN(S->sexm_id, 1)) {

        //    betahivacute(hivS_id == 1 & sexm_id == 1) = hp0(1) * hp0(5) * mcirc(n);
        //    betahiv2(hivS_id == 1 & sexm_id == 1) = hp0(1) * mcirc(n);
        //    betahivart(hivS_id == 1 & sexm_id == 1) = hp0(1) * (1 - hp0(6)) * mcirc(n);
        //    betahivacute_ai(hivS_id == 1 & sexm_id == 1) = hp0(3) * hp0(5) * mcirc(n);
        //    betahiv2_ai(hivS_id == 1 & sexm_id == 1) = hp0(3) * mcirc(n);
        //    betahivart_ai(hivS_id == 1 & sexm_id == 1) = hp0(3) * (1 - hp0(6)) * mcirc(n);

        GET_LIN_ALL(S->betahivacute) = _p_1_5_mcirc;
        GET_LIN_ALL(S->betahiv2) = _p_1_mcirc;
        GET_LIN_ALL(S->betahivart) = _p_1_m6_mcirc;
        GET_LIN_ALL(S->betahivacute_ai) = _p_3_5_mcirc;
        GET_LIN_ALL(S->betahiv2_ai) = _p_3_mcirc;
        GET_LIN_ALL(S->betahivart_ai) = _p_3_m6_mcirc;

        if IS_LIN(S->hpv123_id, 1) {

          //    betahivacute(hivS_id == 1 & sexm_id == 1 & hpv123_id == 1) = hp0(1) * hp0(5) * mcirc(n) * hp0(7);% susceptibility
          //    betahiv2(hivS_id == 1 & sexm_id == 1 & hpv123_id == 1) = hp0(1) * mcirc(n) * hp0(7);
          //    betahivart(hivS_id == 1 & sexm_id == 1 & hpv123_id == 1) = hp0(1) * (1 - hp0(6)) * mcirc(n) * hp0(7);

          GET_LIN_ALL(S->betahivacute) = _p_1_5_7_mcirc;
          GET_LIN_ALL(S->betahiv2) = _p_1_7_mcirc;
          GET_LIN_ALL(S->betahivart) = _p_1_7_m6_mcirc;
        }
      }
    }

    // Break here; the next loop needs _ds, _dh1, _dh2, _dh3 separately, so do non-linear, non-parallel loop...

    LOOP_I LOOP_R LOOP_A LOOP_V1 {

      //    hups1(5, 5, 5, : , 1, : , : , : ) = hups1(5, 5, 5, :, 1, : , : , : ) / 3;
      //    hups2(5, 5, 5, : , 1, : , : , : ) = hups2(5, 5, 5, :, 1, : , : , : ) / 3;
      //    hups3(5, 5, 5, : , 1, : , : , : ) = hups3(5, 5, 5, :, 1, : , : , : ) / 3;

      GET(S->hups1, 5, 5, 5, _di, 1, _dr, _da, _dv1) /= 3.0;
      GET(S->hups2, 5, 5, 5, _di, 1, _dr, _da, _dv1) /= 3.0;
      GET(S->hups3, 5, 5, 5, _di, 1, _dr, _da, _dv1) /= 3.0;

      //    hups1(5, 5, :, : , 1, : , : , : ) = hups1(5, 5, :, : , 1, : , : , : ) / 2;
      //    hups2(5, 5, :, : , 1, : , : , : ) = hups2(5, 5, :, : , 1, : , : , : ) / 2;

      LOOP_H3 {
        GET(S->hups1, 5, 5, _dh3, _di, 1, _dr, _da, _dv1) /= 2.0;
        GET(S->hups2, 5, 5, _dh3, _di, 1, _dr, _da, _dv1) /= 2.0;
      }

      //    hups1(5, :, 5, : , 1, : , : , : ) = hups1(5, :, 5, : , 1, : , : , : ) / 2;
      //    hups3(5, :, 5, : , 1, : , : , : ) = hups3(5, :, 5, : , 1, : , : , : ) / 2;

      LOOP_H2  {
        GET(S->hups1, 5, _dh2, 5, _di, 1, _dr, _da, _dv1) /= 2.0;
        GET(S->hups3, 5, _dh2, 5, _di, 1, _dr, _da, _dv1) /= 2.0;
      }

      //    hups2(:, 5, 5, : , 1, : , : , : ) = hups2(:, 5, 5, : , 1, : , : , : ) / 2;
      //    hups3(:, 5, 5, : , 1, : , : , : ) = hups3(:, 5, 5, : , 1, : , : , : ) / 2;

      LOOP_H1 {
        GET(S->hups2, _dh1, 5, 5, _di, 1, _dr, _da, _dv1) /= 2.0;
        GET(S->hups3, _dh1, 5, 5, _di, 1, _dr, _da, _dv1) /= 2.0;
      }
    } // End Loop i,r,a,v1

//  % VMC - impacts susceptibility

//      % %FOI
//      % here 0 ^ 0 = 1

#define coneffhiv p->hp0[9]
//  % RR by act typeand age
//    commRRai = bp0(60);

#define commRRai p->bp0[60]

    double conuse_commtv_commRRai = conuse_commtv[n] * commRRai;
    double m_conuse_commtv_commRRai = 1.0 - conuse_commtv_commRRai;
    double m_conuse_commtv = 1.0 - conuse_commtv[n];
    double m_coneffhiv = 1.0 - coneffhiv;

#pragma omp parallel for schedule(static, 1)
    PAR_LOOP_LINEAR {

      // This turns out to be the most expensive bit of code, so will try and help the compiler as much as possible here.
      // commRRai is constant, and anything with [n] is constant for this iteration.

      double m_betahivacute = 1.0 - GET_LIN_ALL(S->betahivacute);
      double m_betahivacute_1mconeffhiv = 1.0 - (GET_LIN_ALL(S->betahivacute) * m_coneffhiv);
      double m_betahivacute_ai = 1.0 - GET_LIN_ALL(S->betahivacute_ai);
      double m_betahivacute_ai_1mconeffhiv = 1.0 - (GET_LIN_ALL(S->betahivacute_ai) * m_coneffhiv);

      // Right-hand-sides

      double conuse_tv_conage = conuse_tv[n] * GET_LIN_ALL(conage);
      double m_conuse_tv_conage = 1.0 - conuse_tv_conage;

      double nacts_mconuse_tv_conage = GET_LIN_ALL(S->nacts) * m_conuse_tv_conage;
      double nacts_conuse_tv_conage = GET_LIN_ALL(S->nacts) * conuse_tv_conage;

      double m_p_anal = (1.0 - GET_LIN_ALL(S->p_anal));
      double nacts_mconuse_tv_conage_m_p_anal = nacts_mconuse_tv_conage * m_p_anal;
      double nacts_conuse_tv_conage_m_p_anal = nacts_conuse_tv_conage * m_p_anal;
      double nacts_mconuse_tv_conage_p_anal = nacts_mconuse_tv_conage * GET_LIN_ALL(S->p_anal);
      double nacts_conuse_tv_conage_p_anal = nacts_conuse_tv_conage * GET_LIN_ALL(S->p_anal);

      //         betahivacute_esc = (1 - betahivacute)                        .^ (nacts  .*  (1 - (conuse_tv(n) .* conage)) .* (1 - p_anal)) .*...
      //                            (1 - (betahivacute    .*(1 - coneffhiv))) .^ (nacts  .*  conuse_tv(n)      .* conage    .* (1 - p_anal)) .* ...
      //                            (1 - betahivacute_ai)                     .^ (nacts  .*  (1 - (conuse_tv(n) .* conage)) .*      p_anal) .*...
      //                            (1 - (betahivacute_ai .*(1 - coneffhiv))) .^ (nacts  .*  conuse_tv(n)      .* conage    .*      p_anal)  ;% HIV escape probability

     GET_LIN_ALL(betahivacute_esc) =
        pow(m_betahivacute,                nacts_mconuse_tv_conage_m_p_anal) *
        pow(m_betahivacute_1mconeffhiv,    nacts_conuse_tv_conage_m_p_anal) *
        pow(m_betahivacute_ai,             nacts_mconuse_tv_conage_p_anal) *
        pow(m_betahivacute_ai_1mconeffhiv, nacts_conuse_tv_conage_p_anal);

     //    betahivacute_esc_comm = (1 - betahivacute)                        .^ (nacts_comm  .*  (1 - conuse_commtv(n))              .* (1 - pcomm_anal))    .*...
     //                            (1 - betahivacute    .*(1 - coneffhiv))   .^ (nacts_comm  .*  conuse_commtv(n)                    .* (1 - pcomm_anal)) .* ...
     //                            (1 - betahivacute_ai)                     .^ (nacts_comm  .*  (1 - (conuse_commtv(n).*commRRai))  .* pcomm_anal) .* ...
     //                            (1 - betahivacute_ai .*(1 - coneffhiv))   .^ (nacts_comm  .*  conuse_commtv(n) .* commRRai        .* pcomm_anal); % HIV escape probability


      double m_pcomm_anal = (1.0 - GET_LIN_ALL(S->pcomm_anal));
      double nacts_comm_m_conuse_commtv_m_pcomm_anal = GET_LIN_ALL(S->nacts_comm) * m_conuse_commtv * m_pcomm_anal;
      double nacts_comm_conuse_commtv_m_pcomm_anal = GET_LIN_ALL(S->nacts_comm) * conuse_commtv[n] * m_pcomm_anal;
      double nacts_comm_m_conuse_commtv_commRRai_pcomm_anal = GET_LIN_ALL(S->nacts_comm) * m_conuse_commtv_commRRai * GET_LIN_ALL(S->pcomm_anal);
      double nacts_comm_conuse_commtv_commRRai_pcomm_anal = GET_LIN_ALL(S->nacts_comm) * conuse_commtv_commRRai * GET_LIN_ALL(S->pcomm_anal);

      GET_LIN_ALL(betahivacute_esc_comm) =
        pow(m_betahivacute, nacts_comm_m_conuse_commtv_m_pcomm_anal) *
        pow(m_betahivacute_1mconeffhiv, nacts_comm_conuse_commtv_m_pcomm_anal) *
        pow(m_betahivacute_ai, nacts_comm_m_conuse_commtv_commRRai_pcomm_anal) *
        pow(m_betahivacute_ai_1mconeffhiv, nacts_comm_conuse_commtv_commRRai_pcomm_anal);

      //    betahiv2_esc =          (1 - betahiv2)                       .^ (nacts .* (1 - (conuse_tv(n) .* conage))        .* (1 - p_anal)) .*  ...
      //                            (1 - betahiv2   .*(1 - coneffhiv))   .^ (nacts .* conuse_tv(n) .* conage                .* (1 - p_anal)) .*...
      //                            (1 - betahiv2_ai)                    .^ (nacts .* (1 - (conuse_tv(n) .* conage))        .*      p_anal)  .*...
      //                            (1 - betahiv2_ai.*(1 - coneffhiv))   .^ (nacts .* conuse_tv(n) .* conage                .*      p_anal)  ;% HIV escape probability

      //      betahiv2_esc_comm =   (1 - betahiv2)                       .^ (nacts_comm .* (1 - conuse_commtv(n))           .* (1 - pcomm_anal))     .*...
      //                            (1 - betahiv2    .*(1 - coneffhiv))  .^ (nacts_comm .* conuse_commtv(n)                 .* (1 - pcomm_anal)) .*...
      //                            (1 - betahiv2_ai)                    .^ (nacts_comm .* (1 - conuse_commtv(n).*commRRai) .*      pcomm_anal)     .*...
      //                            (1 - betahiv2_ai .*(1 - coneffhiv))  .^ (nacts_comm .* conuse_commtv(n).*commRRai       .*      pcomm_anal);% HIV escape probability


      double m_betahiv2 = 1.0 - GET_LIN_ALL(S->betahiv2);
      double m_betahiv2_1mconeffhiv = 1.0 - (GET_LIN_ALL(S->betahiv2) * m_coneffhiv);
      double m_betahiv2_ai = 1.0 - GET_LIN_ALL(S->betahiv2_ai);
      double m_betahiv2_ai_1mconeffhiv = 1.0 - (GET_LIN_ALL(S->betahiv2_ai) * m_coneffhiv);

      GET_LIN_ALL(betahiv2_esc) =
        pow(m_betahiv2, nacts_mconuse_tv_conage_m_p_anal) *
        pow(m_betahiv2_1mconeffhiv, nacts_conuse_tv_conage_m_p_anal) *
        pow(m_betahiv2_ai, nacts_mconuse_tv_conage_p_anal) *
        pow(m_betahiv2_ai_1mconeffhiv, nacts_conuse_tv_conage_p_anal);

      GET_LIN_ALL(betahiv2_esc_comm) =
        pow(m_betahiv2, nacts_comm_m_conuse_commtv_m_pcomm_anal) *
        pow(m_betahiv2_1mconeffhiv, nacts_comm_conuse_commtv_m_pcomm_anal) *
        pow(m_betahiv2_ai, nacts_comm_m_conuse_commtv_commRRai_pcomm_anal) *
        pow(m_betahiv2_ai_1mconeffhiv, nacts_comm_conuse_commtv_commRRai_pcomm_anal);

     //      betahivart_esc = (1 - betahivart)                             .^  (nacts      .* (1 - (conuse_tv(n).*conage))         .* (1 - p_anal)) .*...
     //                       (1 - (betahivart   .*(1 - coneffhiv)))       .^  (nacts      .*       conuse_tv(n).*conage           .* (1 - p_anal)) .* ...
     //                       (1 - betahivart_ai)                          .^  (nacts      .* (1 - (conuse_tv(n).*conage))         .* p_anal).*...
     //                       (1 - (betahivart_ai.*(1 - coneffhiv)))       .^  (nacts      .*   conuse_tv(n)  .*  conage           .* p_anal);% HIV escape probability
     //      betahivart_esc_comm = (1 - betahivart)                        .^  (nacts_comm .* (1 - conuse_commtv(n))               .* (1 - pcomm_anal))   .*...
     //                            (1 - (betahivart   .*(1 - coneffhiv)))  .^  (nacts_comm .*   conuse_commtv(n)                   .* (1 - pcomm_anal)) .* ...
     //                            (1 - betahivart_ai)                     .^  (nacts_comm .* (1 - (conuse_commtv(n) .* commRRai)) .* pcomm_anal).*...
     //                            (1 - (betahivart_ai.*(1 - coneffhiv)))  .^  (nacts_comm .*   conuse_commtv(n) .* commRRai       .* pcomm_anal);% HIV escape probability

      double m_betahivart = 1.0 - GET_LIN_ALL(S->betahivart);
      double m_betahivart_1mconeffhiv = 1.0 - (GET_LIN_ALL(S->betahivart) * m_coneffhiv);
      double m_betahivart_ai = 1.0 - GET_LIN_ALL(S->betahivart_ai);
      double m_betahivart_ai_1mconeffhiv = 1.0 - (GET_LIN_ALL(S->betahivart_ai) * m_coneffhiv);

      GET_LIN_ALL(betahivart_esc) =
        pow(m_betahivart, nacts_mconuse_tv_conage_m_p_anal) *
        pow(m_betahivart_1mconeffhiv, nacts_conuse_tv_conage_m_p_anal) *
        pow(m_betahivart_ai, nacts_mconuse_tv_conage_p_anal) *
        pow(m_betahivart_ai_1mconeffhiv, nacts_conuse_tv_conage_p_anal);

      GET_LIN_ALL(betahivart_esc_comm) =
        pow(m_betahivart, nacts_comm_m_conuse_commtv_m_pcomm_anal) *
        pow(m_betahivart_1mconeffhiv, nacts_comm_conuse_commtv_m_pcomm_anal) *
        pow(m_betahivart_ai, nacts_comm_m_conuse_commtv_commRRai_pcomm_anal) *
        pow(m_betahivart_ai_1mconeffhiv, nacts_comm_conuse_commtv_commRRai_pcomm_anal);


      //  % here 1 - 1 = 0 (not a problem)
      //         betahivpp_acute = 1 - betahivacute_esc;% Pr of becoming infected in one partnership
      //         betahivpp_hiv2 = 1 - betahiv2_esc;
      //         betahivpp_art = 1 - betahivart_esc;
      //         betahivpp_acute_comm = 1 - betahivacute_esc_comm;% Pr of becoming infected in one partnership
      //         betahivpp_hiv2_comm = 1 - betahiv2_esc_comm;
      //         betahivpp_art_comm = 1 - betahivart_esc_comm;

      GET_LIN_ALL(betahivpp_acute) = 1.0 - GET_LIN_ALL(betahivacute_esc);
      GET_LIN_ALL(betahivpp_hiv2) = 1.0 - GET_LIN_ALL(betahiv2_esc);
      GET_LIN_ALL(betahivpp_art) = 1.0 - GET_LIN_ALL(betahivart_esc);
      GET_LIN_ALL(betahivpp_acute_comm) = 1.0 - GET_LIN_ALL(betahivacute_esc_comm);
      GET_LIN_ALL(betahivpp_hiv2_comm) = 1.0 - GET_LIN_ALL(betahiv2_esc_comm);
      GET_LIN_ALL(betahivpp_art_comm) = 1.0 - GET_LIN_ALL(betahivart_esc_comm);

    }
#pragma omp barrier
      // Unclear what purpose these assignments serve, since neither betahpv[1,2,3] nor betahpv1pp[1,2,3]pp, nor
      // betahpv[1,2,3]pp_comm get changed again.

      //      betahpv1pp = betahpv1;
      //      betahpv2pp = betahpv2;
      //      betahpv3pp = betahpv3;

#define betahpv1pp S->betahpv1
#define betahpv2pp S->betahpv2
#define betahpv3pp S->betahpv3

      //      betahpv1pp_comm = betahpv1./3; % N.B.lower per partner probability due to << acts pp
      //      betahpv2pp_comm = betahpv2./3;
      //      betahpv3pp_comm = betahpv3./3;

LOOP_LINEAR {
  GET_LIN_ALL(betahpv1pp_comm) = GET_LIN_ALL(S->betahpv1) / 3.0;
  GET_LIN_ALL(betahpv2pp_comm) = GET_LIN_ALL(S->betahpv2) / 3.0;
  GET_LIN_ALL(betahpv3pp_comm) = GET_LIN_ALL(S->betahpv3) / 3.0;
}

// Careful - the ordering of the reshape is non-trivial here...
    LOOP_S LOOP_R LOOP_A6 {
      const int index2 = GET_2D_S_RA6_INDEX(_ds, ((_da - 1) * _SIZE_R) + _dr);

      const double _tmp_popbyhiv1 = GET_4D_ISRA6(popbyhiv, 1, _ds, _dr, _da);
      const double _tmp_popbyhiv2 = GET_4D_ISRA6(popbyhiv, 2, _ds, _dr, _da);
      const double _tmp_popbyhiv3 = GET_4D_ISRA6(popbyhiv, 3, _ds, _dr, _da);
      const double _tmp_popbyhiv4 = GET_4D_ISRA6(popbyhiv, 4, _ds, _dr, _da);

      const double tot_popbyhiv = _tmp_popbyhiv1 + _tmp_popbyhiv2 + _tmp_popbyhiv3 + _tmp_popbyhiv4;

    //      prev_acute = reshape(popbyhiv(2, :, : , : ). / sum(popbyhiv(1:4, : , : , : )), s, r * a);
    //      prev_acute(isnan(prev_acute)) = 0;

      prev_acute[index2] = _tmp_popbyhiv2 / tot_popbyhiv;
      if (isnan(prev_acute[index2])) prev_acute[index2] = 0;

    //      prev_noart = reshape(popbyhiv(3, :, : , : ). / sum(popbyhiv(1:4, : , : , : )), s, r * a);
    //      prev_noart(isnan(prev_noart)) = 0;

      prev_noart[index2] = _tmp_popbyhiv3 / tot_popbyhiv;
      if (isnan(prev_noart[index2])) prev_noart[index2] = 0;


    //      prev_art = reshape(popbyhiv(4, :, : , : ). / sum(popbyhiv(1:4, : , : , : )), s, r * a);
    //      prev_art(isnan(prev_art)) = 0;

      prev_art[index2] = _tmp_popbyhiv4 / tot_popbyhiv;
      if (isnan(prev_art[index2])) prev_art[index2] = 0;

    //      prev_hpv1 = reshape(sum(popbyhpv1(([2 4 5]), :, : , : )). / sum(popbyhpv1(1:5, : , : , : )), s, r * a);
    //      prev_hpv1(isnan(prev_hpv1)) = 0;

      double _tmp_pophpvx_1 = GET_4D_H1SRA6(popbyhpv1, 1, _ds, _dr, _da);
      double _tmp_pophpvx_2 = GET_4D_H1SRA6(popbyhpv1, 2, _ds, _dr, _da);
      double _tmp_pophpvx_3 = GET_4D_H1SRA6(popbyhpv1, 3, _ds, _dr, _da);
      double _tmp_pophpvx_4 = GET_4D_H1SRA6(popbyhpv1, 4, _ds, _dr, _da);
      double _tmp_pophpvx_5 = GET_4D_H1SRA6(popbyhpv1, 5, _ds, _dr, _da);
      double _top = _tmp_pophpvx_2 + _tmp_pophpvx_4 + _tmp_pophpvx_5;
      double _bottom = _tmp_pophpvx_1 + _tmp_pophpvx_2 + _tmp_pophpvx_3 + _tmp_pophpvx_4 + _tmp_pophpvx_5;

      prev_hpv1[index2] = _top / _bottom;
      if (isnan(prev_hpv1[index2])) prev_hpv1[index2] = 0;

    //      prev_hpv2 = reshape(sum(popbyhpv2(([2 4 5]), :, : , : )). / sum(popbyhpv2(1:5, : , : , : )), s, r * a);
    //      prev_hpv2(isnan(prev_hpv2)) = 0;

      _tmp_pophpvx_1 = GET_4D_H2SRA6(popbyhpv2, 1, _ds, _dr, _da);
      _tmp_pophpvx_2 = GET_4D_H2SRA6(popbyhpv2, 2, _ds, _dr, _da);
      _tmp_pophpvx_3 = GET_4D_H2SRA6(popbyhpv2, 3, _ds, _dr, _da);
      _tmp_pophpvx_4 = GET_4D_H2SRA6(popbyhpv2, 4, _ds, _dr, _da);
      _tmp_pophpvx_5 = GET_4D_H2SRA6(popbyhpv2, 5, _ds, _dr, _da);
      _top = _tmp_pophpvx_2 + _tmp_pophpvx_4 + _tmp_pophpvx_5;
      _bottom = _tmp_pophpvx_1 + _tmp_pophpvx_2 + _tmp_pophpvx_3 + _tmp_pophpvx_4 + _tmp_pophpvx_5;

      prev_hpv2[index2] = _top / _bottom;
      if (isnan(prev_hpv2[index2])) prev_hpv2[index2] = 0;

    //      prev_hpv3 = reshape(sum(popbyhpv3(([2 4 5]), :, : , : )). / sum(popbyhpv3(1:5, : , : , : )), s, r * a);
    //      prev_hpv3(isnan(prev_hpv3)) = 0;

      _tmp_pophpvx_1 = GET_4D_H3SRA6(popbyhpv3, 1, _ds, _dr, _da);
      _tmp_pophpvx_2 = GET_4D_H3SRA6(popbyhpv3, 2, _ds, _dr, _da);
      _tmp_pophpvx_3 = GET_4D_H3SRA6(popbyhpv3, 3, _ds, _dr, _da);
      _tmp_pophpvx_4 = GET_4D_H3SRA6(popbyhpv3, 4, _ds, _dr, _da);
      _tmp_pophpvx_5 = GET_4D_H3SRA6(popbyhpv3, 5, _ds, _dr, _da);
      _top = _tmp_pophpvx_2 + _tmp_pophpvx_4 + _tmp_pophpvx_5;
      _bottom = _tmp_pophpvx_1 + _tmp_pophpvx_2 + _tmp_pophpvx_3 + _tmp_pophpvx_4 + _tmp_pophpvx_5;

      prev_hpv3[index2] = _top / _bottom;
      if (isnan(prev_hpv3[index2])) prev_hpv3[index2] = 0;

    }



    for (int _x = 1; _x <= _SIZE_R * _SIZE_A6; _x++) {
      const int pindex_1 = GET_2D_S_RA6_INDEX(1, _x);
      const int pindex_2 = GET_2D_S_RA6_INDEX(2, _x);

      for (int _y = 1; _y <= _SIZE_R * _SIZE_A6; _y++) {
        const int index = GET_2D_RA6RA6_INDEX(_y, _x);
        const double rhof_cf = rhof[index] * c_f[index];
        const double rhom_cm = rhom[index] * c_m[index];

        // tempfoih1 = rhof.*c_f.*reshape(prev_acute(2, :), 1, r * length(eaf(1, :)));
        // tempfoih2 = rhof.*c_f.*reshape(prev_noart(2, :), 1, r * length(eaf(1, :)));
        // tempfoih3 = rhof.*c_f.*reshape(prev_art(2, :), 1, r * length(eaf(1, :)));
        // tempmoih1 = rhom.*c_m.*reshape(prev_acute(1, :), 1, r * length(eaf(1, :)));
        // tempmoih2 = rhom.*c_m.*reshape(prev_noart(1, :), 1, r * length(eaf(1, :)));
        // tempmoih3 = rhom.*c_m.*reshape(prev_art(1, :), 1, r * length(eaf(1, :)));
        // tempfoihpv1 = rhof.*c_f.*reshape(prev_hpv1(2, :), 1, r * length(eaf(1, :)));
        // tempfoihpv2 = rhof.*c_f.*reshape(prev_hpv2(2, :), 1, r * length(eaf(1, :)));
        // tempfoihpv3 = rhof.*c_f.*reshape(prev_hpv3(2, :), 1, r * length(eaf(1, :)));
        // tempmoihpv1 = rhom.*c_m.*reshape(prev_hpv1(1, :), 1, r * length(eaf(1, :)));
        // tempmoihpv2 = rhom.*c_m.*reshape(prev_hpv2(1, :), 1, r * length(eaf(1, :)));
        // tempmoihpv3 = rhom.*c_m.*reshape(prev_hpv3(1, :), 1, r * length(eaf(1, :)));

        tempfoih1[index] = rhof_cf * prev_acute[pindex_2];
        tempfoih2[index] = rhof_cf * prev_noart[pindex_2];
        tempfoih3[index] = rhof_cf * prev_art[pindex_2];
        tempmoih1[index] = rhom_cm * prev_acute[pindex_1];
        tempmoih2[index] = rhom_cm * prev_noart[pindex_1];
        tempmoih3[index] = rhom_cm * prev_art[pindex_1];
        tempfoihpv1[index] = rhof_cf * prev_hpv1[pindex_2];
        tempfoihpv2[index] = rhof_cf * prev_hpv2[pindex_2];
        tempfoihpv3[index] = rhof_cf * prev_hpv3[pindex_2];
        tempmoihpv1[index] = rhom_cm * prev_hpv1[pindex_1];
        tempmoihpv2[index] = rhom_cm * prev_hpv2[pindex_1];
        tempmoihpv3[index] = rhom_cm * prev_hpv3[pindex_1];
      }
    }

        // tempfoih1_full = [tempfoih1(1:3, : ); repmat(tempfoih1(4:6, : ), 2, 1); repmat(tempfoih1(7:9, : ), 2, 1); repmat(tempfoih1(10:12, : ), 3, 1); repmat(tempfoih1(13:15, : ), 2, 1); repmat(tempfoih1(16:18, : ), 3, 1)];
        // tempfoih2_full = [tempfoih2(1:3, : ); repmat(tempfoih2(4:6, : ), 2, 1); repmat(tempfoih2(7:9, : ), 2, 1); repmat(tempfoih2(10:12, : ), 3, 1); repmat(tempfoih2(13:15, : ), 2, 1); repmat(tempfoih2(16:18, : ), 3, 1)];
        // tempfoih3_full = [tempfoih3(1:3, : ); repmat(tempfoih3(4:6, : ), 2, 1); repmat(tempfoih3(7:9, : ), 2, 1); repmat(tempfoih3(10:12, : ), 3, 1); repmat(tempfoih3(13:15, : ), 2, 1); repmat(tempfoih3(16:18, : ), 3, 1)];
        // tempmoih1_full = [tempmoih1(1:3, : ); repmat(tempmoih1(4:6, : ), 2, 1); repmat(tempmoih1(7:9, : ), 2, 1); repmat(tempmoih1(10:12, : ), 3, 1); repmat(tempmoih1(13:15, : ), 2, 1); repmat(tempmoih1(16:18, : ), 3, 1)];
        // tempmoih2_full = [tempmoih2(1:3, : ); repmat(tempmoih2(4:6, : ), 2, 1); repmat(tempmoih2(7:9, : ), 2, 1); repmat(tempmoih2(10:12, : ), 3, 1); repmat(tempmoih2(13:15, : ), 2, 1); repmat(tempmoih2(16:18, : ), 3, 1)];
        // tempmoih3_full = [tempmoih3(1:3, : ); repmat(tempmoih3(4:6, : ), 2, 1); repmat(tempmoih3(7:9, : ), 2, 1); repmat(tempmoih3(10:12, : ), 3, 1); repmat(tempmoih3(13:15, : ), 2, 1); repmat(tempmoih3(16:18, : ), 3, 1)];
        // tempfoihpv1_full = [tempfoihpv1(1:3, : ); repmat(tempfoihpv1(4:6, : ), 2, 1); repmat(tempfoihpv1(7:9, : ), 2, 1); repmat(tempfoihpv1(10:12, : ), 3, 1); repmat(tempfoihpv1(13:15, : ), 2, 1); repmat(tempfoihpv1(16:18, : ), 3, 1)];
        // tempfoihpv2_full = [tempfoihpv2(1:3, : ); repmat(tempfoihpv2(4:6, : ), 2, 1); repmat(tempfoihpv2(7:9, : ), 2, 1); repmat(tempfoihpv2(10:12, : ), 3, 1); repmat(tempfoihpv2(13:15, : ), 2, 1); repmat(tempfoihpv2(16:18, : ), 3, 1)];
        // tempfoihpv3_full = [tempfoihpv3(1:3, : ); repmat(tempfoihpv3(4:6, : ), 2, 1); repmat(tempfoihpv3(7:9, : ), 2, 1); repmat(tempfoihpv3(10:12, : ), 3, 1); repmat(tempfoihpv3(13:15, : ), 2, 1); repmat(tempfoihpv3(16:18, : ), 3, 1)];
        // tempmoihpv1_full = [tempmoihpv1(1:3, : ); repmat(tempmoihpv1(4:6, : ), 2, 1); repmat(tempmoihpv1(7:9, : ), 2, 1); repmat(tempmoihpv1(10:12, : ), 3, 1); repmat(tempmoihpv1(13:15, : ), 2, 1); repmat(tempmoihpv1(16:18, : ), 3, 1)];
        // tempmoihpv2_full = [tempmoihpv2(1:3, : ); repmat(tempmoihpv2(4:6, : ), 2, 1); repmat(tempmoihpv2(7:9, : ), 2, 1); repmat(tempmoihpv2(10:12, : ), 3, 1); repmat(tempmoihpv2(13:15, : ), 2, 1); repmat(tempmoihpv2(16:18, : ), 3, 1)];
        // tempmoihpv3_full = [tempmoihpv3(1:3, : ); repmat(tempmoihpv3(4:6, : ), 2, 1); repmat(tempmoihpv3(7:9, : ), 2, 1); repmat(tempmoihpv3(10:12, : ), 3, 1); repmat(tempmoihpv3(13:15, : ), 2, 1); repmat(tempmoihpv3(16:18, : ), 3, 1)];

      // See matlab.h for this - basically, it is copying chunks of 3 * 18 doubles, with chunk-repetitions (1,2,2,3,2,3) into a matrix.

    EXPAND_MATRIX_A6_A6(tempfoih1_full, tempfoih1);
    EXPAND_MATRIX_A6_A6(tempfoih2_full, tempfoih2);
    EXPAND_MATRIX_A6_A6(tempfoih3_full, tempfoih3);
    EXPAND_MATRIX_A6_A6(tempmoih1_full, tempmoih1);
    EXPAND_MATRIX_A6_A6(tempmoih2_full, tempmoih2);
    EXPAND_MATRIX_A6_A6(tempmoih3_full, tempmoih3);

    EXPAND_MATRIX_A6_A6(tempfoihpv1_full, tempfoihpv1);
    EXPAND_MATRIX_A6_A6(tempfoihpv2_full, tempfoihpv2);
    EXPAND_MATRIX_A6_A6(tempfoihpv3_full, tempfoihpv3);
    EXPAND_MATRIX_A6_A6(tempmoihpv1_full, tempmoihpv1);
    EXPAND_MATRIX_A6_A6(tempmoihpv2_full, tempmoihpv2);
    EXPAND_MATRIX_A6_A6(tempmoihpv3_full, tempmoihpv3);

    //  % comercial sex only between fswand clients

    // tempfoih1_comm   = reshape(pcr_comm_bal(1, :, : ), r * length(eaf(1, :)), 1).*sum(reshape(pshipsprop_age_com(2, :, : ), r * length(eaf(1, :)), 1).*reshape(prev_acute(2, :), r * length(eaf(1, :)), 1));
    // tempfoih2_comm   = reshape(pcr_comm_bal(1, :, : ), r * length(eaf(1, :)), 1).*sum(reshape(pshipsprop_age_com(2, :, : ), r * length(eaf(1, :)), 1).*reshape(prev_noart(2, :), r * length(eaf(1, :)), 1));
    // tempfoih3_comm   = reshape(pcr_comm_bal(1, :, : ), r * length(eaf(1, :)), 1).*sum(reshape(pshipsprop_age_com(2, :, : ), r * length(eaf(1, :)), 1).*reshape(prev_art(2, :), r * length(eaf(1, :)), 1));

    double sum = 0;
    LOOP_R LOOP_A6 sum += (GET_3D_SRA6(pshipsprop_age_com, 2, _dr, _da) * GET_2D_S_RA6(prev_acute, 2, ((_da - 1) * _SIZE_R) + _dr));
    LOOP_R LOOP_A6 tempfoih1_comm[((_da - 1) * _SIZE_R) + _dr] = GET_3D_SRA6(pcr_comm_bal, 1, _dr, _da) * sum;

    sum = 0;
    LOOP_R LOOP_A6 sum += (GET_3D_SRA6(pshipsprop_age_com, 2, _dr, _da) * GET_2D_S_RA6(prev_noart, 2, ((_da - 1) * _SIZE_R) + _dr));
    LOOP_R LOOP_A6 tempfoih2_comm[((_da - 1) * _SIZE_R) + _dr] = GET_3D_SRA6(pcr_comm_bal, 1, _dr, _da) * sum;

    sum = 0;
    LOOP_R LOOP_A6 sum += (GET_3D_SRA6(pshipsprop_age_com, 2, _dr, _da) * GET_2D_S_RA6(prev_art, 2, ((_da - 1) * _SIZE_R) + _dr));
    LOOP_R LOOP_A6 tempfoih3_comm[((_da - 1) * _SIZE_R) + _dr] = GET_3D_SRA6(pcr_comm_bal, 1, _dr, _da) * sum;

      // tempmoih1_comm   = reshape(pcr_comm_bal(2, :, : ), r * length(eaf(1, :)), 1).*sum(reshape(pshipsprop_age_com(1, :, : ), r * length(eaf(1, :)), 1).*reshape(prev_acute(1, :), r * length(eaf(1, :)), 1));
      // tempmoih2_comm   = reshape(pcr_comm_bal(2, :, : ), r * length(eaf(1, :)), 1).*sum(reshape(pshipsprop_age_com(1, :, : ), r * length(eaf(1, :)), 1).*reshape(prev_noart(1, :), r * length(eaf(1, :)), 1));
      // tempmoih3_comm   = reshape(pcr_comm_bal(2, :, : ), r * length(eaf(1, :)), 1).*sum(reshape(pshipsprop_age_com(1, :, : ), r * length(eaf(1, :)), 1).*reshape(prev_art(1, :), r * length(eaf(1, :)), 1));

    sum = 0;
    LOOP_R LOOP_A6 sum += (GET_3D_SRA6(pshipsprop_age_com, 1, _dr, _da) * GET_2D_S_RA6(prev_acute, 1, ((_da - 1) * _SIZE_R) + _dr));
    LOOP_R LOOP_A6 tempmoih1_comm[((_da - 1) * _SIZE_R) + _dr] = GET_3D_SRA6(pcr_comm_bal, 2, _dr, _da) * sum;

    sum = 0;
    LOOP_R LOOP_A6 sum += (GET_3D_SRA6(pshipsprop_age_com, 1, _dr, _da) * GET_2D_S_RA6(prev_noart, 1, ((_da - 1) * _SIZE_R) + _dr));
    LOOP_R LOOP_A6 tempmoih2_comm[((_da - 1) * _SIZE_R) + _dr] = GET_3D_SRA6(pcr_comm_bal, 2, _dr, _da) * sum;

    sum = 0;
    LOOP_R LOOP_A6 sum += (GET_3D_SRA6(pshipsprop_age_com, 1, _dr, _da) * GET_2D_S_RA6(prev_art, 1, ((_da - 1) * _SIZE_R) + _dr));
    LOOP_R LOOP_A6 tempmoih3_comm[((_da - 1) * _SIZE_R) + _dr] = GET_3D_SRA6(pcr_comm_bal, 2, _dr, _da) * sum;

    // tempfoihpv1_comm = reshape(pcr_comm_bal(1, :, : ), r * length(eaf(1, :)), 1).*
    //    sum(reshape(pshipsprop_age_com(2, :, : ), r * length(eaf(1, :)), 1).*
    //    reshape(prev_hpv1(2, :), r * length(eaf(1, :)), 1));

    sum = 0;
    LOOP_R LOOP_A6 sum += (GET_3D_SRA6(pshipsprop_age_com, 2, _dr, _da) *
                          GET_2D_S_RA6(prev_hpv1, 2, ((_da - 1) * _SIZE_R) + _dr));

    LOOP_R LOOP_A6 {

    tempfoihpv1_comm[((_da - 1) * _SIZE_R) + _dr] =
                          GET_3D_SRA6(pcr_comm_bal, 1, _dr, _da) * sum;
                          }

    // tempfoihpv2_comm = reshape(pcr_comm_bal(1, :, : ), r * length(eaf(1, :)), 1).*sum(reshape(pshipsprop_age_com(2, :, : ), r * length(eaf(1, :)), 1).*reshape(prev_hpv2(2, :), r * length(eaf(1, :)), 1));
    // tempfoihpv3_comm = reshape(pcr_comm_bal(1, :, : ), r * length(eaf(1, :)), 1).*sum(reshape(pshipsprop_age_com(2, :, : ), r * length(eaf(1, :)), 1).*reshape(prev_hpv3(2, :), r * length(eaf(1, :)), 1));

    sum = 0;
    LOOP_R LOOP_A6 sum += (GET_3D_SRA6(pshipsprop_age_com, 2, _dr, _da) * GET_2D_S_RA6(prev_hpv2, 2, ((_da - 1) * _SIZE_R) + _dr));
    LOOP_R LOOP_A6 tempfoihpv2_comm[((_da - 1) * _SIZE_R) + _dr] = GET_3D_SRA6(pcr_comm_bal, 1, _dr, _da) * sum;

    sum = 0;
    LOOP_R LOOP_A6 sum += (GET_3D_SRA6(pshipsprop_age_com, 2, _dr, _da) * GET_2D_S_RA6(prev_hpv3, 2, ((_da - 1) * _SIZE_R) + _dr));
    LOOP_R LOOP_A6 tempfoihpv3_comm[((_da - 1) * _SIZE_R) + _dr] = GET_3D_SRA6(pcr_comm_bal, 1, _dr, _da) * sum;

    // tempmoihpv1_comm = reshape(pcr_comm_bal(2, :, : ), r * length(eaf(1, :)), 1).*sum(reshape(pshipsprop_age_com(1, :, : ), r * length(eaf(1, :)), 1).*reshape(prev_hpv1(1, :), r * length(eaf(1, :)), 1));
    // tempmoihpv2_comm = reshape(pcr_comm_bal(2, :, : ), r * length(eaf(1, :)), 1).*sum(reshape(pshipsprop_age_com(1, :, : ), r * length(eaf(1, :)), 1).*reshape(prev_hpv2(1, :), r * length(eaf(1, :)), 1));
    // tempmoihpv3_comm = reshape(pcr_comm_bal(2, :, : ), r * length(eaf(1, :)), 1).*sum(reshape(pshipsprop_age_com(1, :, : ), r * length(eaf(1, :)), 1).*reshape(prev_hpv3(1, :), r * length(eaf(1, :)), 1));

    sum = 0;
    LOOP_R LOOP_A6 sum += (GET_3D_SRA6(pshipsprop_age_com, 1, _dr, _da) * GET_2D_S_RA6(prev_hpv1, 1, ((_da - 1) * _SIZE_R) + _dr));
    LOOP_R LOOP_A6 tempmoihpv1_comm[((_da - 1) * _SIZE_R) + _dr] = GET_3D_SRA6(pcr_comm_bal, 2, _dr, _da) * sum;

    sum = 0;
    LOOP_R LOOP_A6 sum += (GET_3D_SRA6(pshipsprop_age_com, 1, _dr, _da) * GET_2D_S_RA6(prev_hpv2, 1, ((_da - 1) * _SIZE_R) + _dr));
    LOOP_R LOOP_A6 tempmoihpv2_comm[((_da - 1) * _SIZE_R) + _dr] = GET_3D_SRA6(pcr_comm_bal, 2, _dr, _da) * sum;

    sum = 0;
    LOOP_R LOOP_A6 sum += (GET_3D_SRA6(pshipsprop_age_com, 1, _dr, _da) * GET_2D_S_RA6(prev_hpv3, 1, ((_da - 1) * _SIZE_R) + _dr));
    LOOP_R LOOP_A6 tempmoihpv3_comm[((_da - 1) * _SIZE_R) + _dr] = GET_3D_SRA6(pcr_comm_bal, 2, _dr, _da) * sum;

    // tempfoih1_comm_full = [tempfoih1_comm(1:3); repmat(tempfoih1_comm(4:6), 2, 1); repmat(tempfoih1_comm(7:9), 2, 1); repmat(tempfoih1_comm(10:12), 3, 1); repmat(tempfoih1_comm(13:15), 2, 1); repmat(tempfoih1_comm(16:18, : ), 3, 1)];
    // tempfoih2_comm_full = [tempfoih2_comm(1:3); repmat(tempfoih2_comm(4:6), 2, 1); repmat(tempfoih2_comm(7:9), 2, 1); repmat(tempfoih2_comm(10:12), 3, 1); repmat(tempfoih2_comm(13:15), 2, 1); repmat(tempfoih2_comm(16:18, : ), 3, 1)];
    // tempfoih3_comm_full = [tempfoih3_comm(1:3); repmat(tempfoih3_comm(4:6), 2, 1); repmat(tempfoih3_comm(7:9), 2, 1); repmat(tempfoih3_comm(10:12), 3, 1); repmat(tempfoih3_comm(13:15), 2, 1); repmat(tempfoih3_comm(16:18, : ), 3, 1)];
    // tempmoih1_comm_full = [tempmoih1_comm(1:3); repmat(tempmoih1_comm(4:6), 2, 1); repmat(tempmoih1_comm(7:9), 2, 1); repmat(tempmoih1_comm(10:12), 3, 1); repmat(tempmoih1_comm(13:15), 2, 1); repmat(tempmoih1_comm(16:18, : ), 3, 1)];
    // tempmoih2_comm_full = [tempmoih2_comm(1:3); repmat(tempmoih2_comm(4:6), 2, 1); repmat(tempmoih2_comm(7:9), 2, 1); repmat(tempmoih2_comm(10:12), 3, 1); repmat(tempmoih2_comm(13:15), 2, 1); repmat(tempmoih2_comm(16:18, : ), 3, 1)];
    // tempmoih3_comm_full = [tempmoih3_comm(1:3); repmat(tempmoih3_comm(4:6), 2, 1); repmat(tempmoih3_comm(7:9), 2, 1); repmat(tempmoih3_comm(10:12), 3, 1); repmat(tempmoih3_comm(13:15), 2, 1); repmat(tempmoih3_comm(16:18, : ), 3, 1)];
    // tempfoihpv1_comm_full = [tempfoihpv1_comm(1:3); repmat(tempfoihpv1_comm(4:6), 2, 1); repmat(tempfoihpv1_comm(7:9), 2, 1); repmat(tempfoihpv1_comm(10:12), 3, 1); repmat(tempfoihpv1_comm(13:15), 2, 1); repmat(tempfoihpv1_comm(16:18, : ), 3, 1)];
    // tempfoihpv2_comm_full = [tempfoihpv2_comm(1:3); repmat(tempfoihpv2_comm(4:6), 2, 1); repmat(tempfoihpv2_comm(7:9), 2, 1); repmat(tempfoihpv2_comm(10:12), 3, 1); repmat(tempfoihpv2_comm(13:15), 2, 1); repmat(tempfoihpv2_comm(16:18, : ), 3, 1)];
    // tempfoihpv3_comm_full = [tempfoihpv3_comm(1:3); repmat(tempfoihpv3_comm(4:6), 2, 1); repmat(tempfoihpv3_comm(7:9), 2, 1); repmat(tempfoihpv3_comm(10:12), 3, 1); repmat(tempfoihpv3_comm(13:15), 2, 1); repmat(tempfoihpv3_comm(16:18, : ), 3, 1)];
    // tempmoihpv1_comm_full = [tempmoihpv1_comm(1:3); repmat(tempmoihpv1_comm(4:6), 2, 1); repmat(tempmoihpv1_comm(7:9), 2, 1); repmat(tempmoihpv1_comm(10:12), 3, 1); repmat(tempmoihpv1_comm(13:15), 2, 1); repmat(tempmoihpv1_comm(16:18, : ), 3, 1)];
    // tempmoihpv2_comm_full = [tempmoihpv2_comm(1:3); repmat(tempmoihpv2_comm(4:6), 2, 1); repmat(tempmoihpv2_comm(7:9), 2, 1); repmat(tempmoihpv2_comm(10:12), 3, 1); repmat(tempmoihpv2_comm(13:15), 2, 1); repmat(tempmoihpv2_comm(16:18, : ), 3, 1)];
    // tempmoihpv3_comm_full = [tempmoihpv3_comm(1:3); repmat(tempmoihpv3_comm(4:6), 2, 1); repmat(tempmoihpv3_comm(7:9), 2, 1); repmat(tempmoihpv3_comm(10:12), 3, 1); repmat(tempmoihpv3_comm(13:15), 2, 1); repmat(tempmoihpv3_comm(16:18, : ), 3, 1)];

    EXPAND_MATRIX_A6_1(tempfoih1_comm_full, tempfoih1_comm);
    EXPAND_MATRIX_A6_1(tempfoih2_comm_full, tempfoih2_comm);
    EXPAND_MATRIX_A6_1(tempfoih3_comm_full, tempfoih3_comm);
    EXPAND_MATRIX_A6_1(tempmoih1_comm_full, tempmoih1_comm);
    EXPAND_MATRIX_A6_1(tempmoih2_comm_full, tempmoih2_comm);
    EXPAND_MATRIX_A6_1(tempmoih3_comm_full, tempmoih3_comm);
    EXPAND_MATRIX_A6_1(tempfoihpv1_comm_full, tempfoihpv1_comm);
    EXPAND_MATRIX_A6_1(tempfoihpv2_comm_full, tempfoihpv2_comm);
    EXPAND_MATRIX_A6_1(tempfoihpv3_comm_full, tempfoihpv3_comm);
    EXPAND_MATRIX_A6_1(tempmoihpv1_comm_full, tempmoihpv1_comm);
    EXPAND_MATRIX_A6_1(tempmoihpv2_comm_full, tempmoihpv2_comm);
    EXPAND_MATRIX_A6_1(tempmoihpv3_comm_full, tempmoihpv3_comm);

    //  % Final FOI

    //    for at = 1:a
    //      for rt = 1 : r

    int at_rt = _SIZE_A * _SIZE_R;

#pragma omp parallel for schedule(static, 1)
    for (int thread = 0; thread < p->threads; thread++)
      for (int loop_atrt = thread; loop_atrt < at_rt; loop_atrt += p->threads) {
        const int _da = 1 + (loop_atrt % _SIZE_A);
        const int _dr = 1 + (loop_atrt / _SIZE_A);
        const int index = ((_da - 1) * _SIZE_R) + _dr;
        const double sum2_tempfoih1_full_i = sum_second_dim(tempfoih1_full, index);
        const double sum2_tempfoih2_full_i = sum_second_dim(tempfoih2_full, index);
        const double sum2_tempfoih3_full_i = sum_second_dim(tempfoih3_full, index);
        const double sum2_tempmoih1_full_i = sum_second_dim(tempmoih1_full, index);
        const double sum2_tempmoih2_full_i = sum_second_dim(tempmoih2_full, index);
        const double sum2_tempmoih3_full_i = sum_second_dim(tempmoih3_full, index);
        const double sum2_tempfoihpv1_full_i = sum_second_dim(tempfoihpv1_full, index);
        const double sum2_tempfoihpv2_full_i = sum_second_dim(tempfoihpv2_full, index);
        const double sum2_tempfoihpv3_full_i = sum_second_dim(tempfoihpv3_full, index);
        const double sum2_tempmoihpv1_full_i = sum_second_dim(tempmoihpv1_full, index);
        const double sum2_tempmoihpv2_full_i = sum_second_dim(tempmoihpv2_full, index);
        const double sum2_tempmoihpv3_full_i = sum_second_dim(tempmoihpv3_full, index);

        //      foihiv(:, : , : , : , 1, rt, at, : ) = betahivpp_acute(:, : , : , : , 1, rt, at, : ).*sum(tempfoih1_full(at * r + rt - r, :), 2) + ...
        //                                              betahivpp_hiv2(:, : , : , : , 1, rt, at, : ).*sum(tempfoih2_full(at * r + rt - r, :), 2) + ...
        //                                               betahivpp_art(:, : , : , : , 1, rt, at, : ).*sum(tempfoih3_full(at * r + rt - r, :), 2) + ...
        //                                        betahivpp_acute_comm(:, : , : , : , 1, rt, at, : ).*tempfoih1_comm_full(at * r + rt - r, :) + ...
        //                                         betahivpp_hiv2_comm(:, : , : , : , 1, rt, at, : ).*tempfoih2_comm_full(at * r + rt - r, :) + ...
        //                                          betahivpp_art_comm(:, : , : , : , 1, rt, at, : ).*tempfoih3_comm_full(at * r + rt - r, :);


        LOOP_H1 LOOP_H2 LOOP_H3 LOOP_I LOOP_V1 {
           int _arrindex1 = GET_INDEX(_dh1,_dh2,_dh3,_di,1,_dr,_da,_dv1);
           int _arrindex2 = GET_INDEX(_dh1,_dh2,_dh3,_di,2,_dr,_da,_dv1);

                        foihiv[_arrindex1] =
               betahivpp_acute[_arrindex1] * sum2_tempfoih1_full_i +
                betahivpp_hiv2[_arrindex1] * sum2_tempfoih2_full_i +
                 betahivpp_art[_arrindex1] * sum2_tempfoih3_full_i +
          betahivpp_acute_comm[_arrindex1] * tempfoih1_comm_full[index] +
           betahivpp_hiv2_comm[_arrindex1] * tempfoih2_comm_full[index] +
            betahivpp_art_comm[_arrindex1] * tempfoih3_comm_full[index];

          //      foihiv(:, : , : , : , 2, rt, at, : ) = betahivpp_acute(:, : , : , : , 2, rt, at, : ).*sum(tempmoih1_full(at * r + rt - r, :), 2) + ...
          //                                              betahivpp_hiv2(:, : , : , : , 2, rt, at, : ).*sum(tempmoih2_full(at * r + rt - r, :), 2) + ...
          //                                               betahivpp_art(:, : , : , : , 2, rt, at, : ).*sum(tempmoih3_full(at * r + rt - r, :), 2) + ...
          //                                        betahivpp_acute_comm(:, : , : , : , 2, rt, at, : ).*tempmoih1_comm_full(at * r + rt - r, :) + ...
          //                                         betahivpp_hiv2_comm(:, : , : , : , 2, rt, at, : ).*tempmoih2_comm_full(at * r + rt - r, :) + ...
          //                                          betahivpp_art_comm(:, : , : , : , 2, rt, at, : ).*tempmoih3_comm_full(at * r + rt - r, :);

                        foihiv[_arrindex2] =
               betahivpp_acute[_arrindex2] * sum2_tempmoih1_full_i +
                betahivpp_hiv2[_arrindex2] * sum2_tempmoih2_full_i +
                 betahivpp_art[_arrindex2] * sum2_tempmoih3_full_i +
          betahivpp_acute_comm[_arrindex2] * tempmoih1_comm_full[index] +
           betahivpp_hiv2_comm[_arrindex2] * tempmoih2_comm_full[index] +
            betahivpp_art_comm[_arrindex2] * tempmoih3_comm_full[index];

          //      foihpv1(:, : , : , : , 1, rt, at, : ) = betahpv1pp(:, : , : , : , 1, rt, at, : ).*sum(tempfoihpv1_full(at * r + rt - r, :), 2) + ...
          //        betahpv1pp_comm(:, : , : , : , 1, rt, at, : ).*tempfoihpv1_comm_full(at * r + rt - r, :);
          //      foihpv2(:, : , : , : , 1, rt, at, : ) = betahpv2pp(:, : , : , : , 1, rt, at, : ).*sum(tempfoihpv2_full(at * r + rt - r, :), 2) + ...
          //        betahpv2pp_comm(:, : , : , : , 1, rt, at, : ).*tempfoihpv2_comm_full(at * r + rt - r, :);
          //      foihpv3(:, : , : , : , 1, rt, at, : ) = betahpv3pp(:, : , : , : , 1, rt, at, : ).*sum(tempfoihpv3_full(at * r + rt - r, :), 2) + ...
          //        betahpv3pp_comm(:, : , : , : , 1, rt, at, : ).*tempfoihpv3_comm_full(at * r + rt - r, :);
          //      foihpv1(:, : , : , : , 2, rt, at, : ) = betahpv1pp(:, : , : , : , 2, rt, at, : ).*sum(tempmoihpv1_full(at * r + rt - r, :), 2) + ...
          //        betahpv1pp_comm(:, : , : , : , 2, rt, at, : ).*tempmoihpv1_comm_full(at * r + rt - r, :);
          //      foihpv2(:, : , : , : , 2, rt, at, : ) = betahpv2pp(:, : , : , : , 2, rt, at, : ).*sum(tempmoihpv2_full(at * r + rt - r, :), 2) + ...
          //        betahpv2pp_comm(:, : , : , : , 2, rt, at, : ).*tempmoihpv2_comm_full(at * r + rt - r, :);
          //      foihpv3(:, : , : , : , 2, rt, at, : ) = betahpv3pp(:, : , : , : , 2, rt, at, : ).*sum(tempmoihpv3_full(at * r + rt - r, :), 2) + ...
          //        betahpv3pp_comm(:, : , : , : , 2, rt, at, : ).*tempmoihpv3_comm_full(at * r + rt - r, :);
          foihpv1[_arrindex1] = betahpv1pp[_arrindex1] * sum2_tempfoihpv1_full_i + betahpv1pp_comm[_arrindex1] * tempfoihpv1_comm_full[index];
          foihpv2[_arrindex1] = betahpv2pp[_arrindex1] * sum2_tempfoihpv2_full_i + betahpv2pp_comm[_arrindex1] * tempfoihpv2_comm_full[index];
          foihpv3[_arrindex1] = betahpv3pp[_arrindex1] * sum2_tempfoihpv3_full_i + betahpv3pp_comm[_arrindex1] * tempfoihpv3_comm_full[index];

          foihpv1[_arrindex2] = betahpv1pp[_arrindex2] * sum2_tempmoihpv1_full_i + betahpv1pp_comm[_arrindex2] * tempmoihpv1_comm_full[index];
          foihpv2[_arrindex2] = betahpv2pp[_arrindex2] * sum2_tempmoihpv2_full_i + betahpv2pp_comm[_arrindex2] * tempmoihpv2_comm_full[index];
          foihpv3[_arrindex2] = betahpv3pp[_arrindex2] * sum2_tempmoihpv3_full_i + betahpv3pp_comm[_arrindex2] * tempmoihpv3_comm_full[index];

        }
      }
#pragma omp barrier
#define phi p->pp0[43]
#define phi_hiv p->pp0[44]
#define phi2 p->pp0[45]

    LOOP_LINEAR

      //    foihpv3(vacc_id == 1) = (1 - phi2).*foihpv3(vacc_id == 1);
      //    foihpv1(vacc_id == 1 & hivS_id == 1) = (1 - phi).*foihpv1(vacc_id == 1 & hivS_id == 1);
      //    foihpv2(vacc_id == 1 & hivS_id == 1) = (1 - phi).*foihpv2(vacc_id == 1 & hivS_id == 1);

      if IS_LIN(S->vacc_id, 1) {
        GET_LIN_ALL(foihpv3) = (1.0 - phi2) * GET_LIN_ALL(foihpv3);

        if IS_LIN(S->hivS_id, 1) {
          GET_LIN_ALL(foihpv1) = (1.0 - phi) * GET_LIN_ALL(foihpv1);
          GET_LIN_ALL(foihpv2) = (1.0 - phi) * GET_LIN_ALL(foihpv2);
        }

        //    foihpv1(vacc_id == 1 & hiv_id == 1) = (1 - phi_hiv).*foihpv1(vacc_id == 1 & hiv_id == 1);
        //    foihpv2(vacc_id == 1 & hiv_id == 1) = (1 - phi_hiv).*foihpv2(vacc_id == 1 & hiv_id == 1);

        if IS_LIN(S->hiv_id, 1) {
          GET_LIN_ALL(foihpv1) = (1.0 - phi_hiv) * GET_LIN_ALL(foihpv1);
          GET_LIN_ALL(foihpv2) = (1.0 - phi_hiv) * GET_LIN_ALL(foihpv2);
        }
      }

    //      %% AGING
    //      %% DEATHS AND AGING OUT EVENTS
    //      % calculate how many are removed from the model(any reason) by sexand
    //      % risk group -- these are still in yearly measures

    //      % HIV related deaths

    //     agein1 = nu(:, : , : , 2, : , : , : , : ).*popn(:, : , : , 2, : , : , : , : )...
    //            + nu(:, : , : , 3, : , : , : , : ).*popn(:, : , : , 3, : , : , : , : )...
    //            + nu(:, : , : , 4, : , : , : , : ).*popn(:, : , : , 4, : , : , : , : );

    LOOP_H1 LOOP_H2 LOOP_H3 LOOP_S LOOP_R LOOP_A LOOP_V1
       GET(agein1, _dh1, _dh2, _dh3, 1, _ds, _dr, _da, _dv1) =
        GET(S->nu, _dh1, _dh2, _dh3, 2, _ds, _dr, _da, _dv1) * GET(popn, _dh1, _dh2, _dh3, 2, _ds, _dr, _da, _dv1) +
        GET(S->nu, _dh1, _dh2, _dh3, 3, _ds, _dr, _da, _dv1) * GET(popn, _dh1, _dh2, _dh3, 3, _ds, _dr, _da, _dv1) +
        GET(S->nu, _dh1, _dh2, _dh3, 4, _ds, _dr, _da, _dv1) * GET(popn, _dh1, _dh2, _dh3, 4, _ds, _dr, _da, _dv1);

    //    % Aging out
    //      agein2 = alpha(:, : , : , : , : , : , a, : ).*popn(:, : , : , : , : , : , a, : );
    LOOP_H1 LOOP_H2 LOOP_H3 LOOP_I LOOP_S LOOP_R LOOP_V1
        GET(agein2, _dh1, _dh2, _dh3, _di, _ds, _dr, 1, _dv1) =
      GET(S->alpha, _dh1, _dh2, _dh3, _di, _ds, _dr, _SIZE_A, _dv1) *
          GET(popn, _dh1, _dh2, _dh3, _di, _ds, _dr, _SIZE_A, _dv1);

    //    % Background death rate
    //      agein3 = mu(:, : , : , : , : , : , : , : ).*popn(:, : , : , : , : , : , : , : );

    LOOP_LINEAR GET_LIN_ALL(agein3) = GET_LIN_ALL(S->mu) * GET_LIN_ALL(popn);

    //      agein4 = (ups1(5, :, : , : , 1, : , : , : ) +
    //               hups1(5, :, : , : , 1, : , : , : )).*
    //                popn(5, :, : , : , 1, : , : , : );

    LOOP_H2 LOOP_H3 LOOP_I LOOP_R LOOP_A LOOP_V1
            GET(agein4, 1, _dh2, _dh3, _di, 1, _dr, _da, _dv1) =
          (GET(S->ups1, 5, _dh2, _dh3, _di, 1, _dr, _da, _dv1) +
          GET(S->hups1, 5, _dh2, _dh3, _di, 1, _dr, _da, _dv1)) *
              GET(popn, 5, _dh2, _dh3, _di, 1, _dr, _da, _dv1);

    //      agein5 = (ups2(:, 5, : , : , 1, : , : , : ) +
    //                hups2(:, 5, : , : , 1, : , : , : )).*popn(:, 5, : , : , 1, : , : , : );

    LOOP_H1 LOOP_H3 LOOP_I LOOP_R LOOP_A LOOP_V1
         GET(agein5, _dh1, 1, _dh3, _di, 1, _dr, _da, _dv1) =
       (GET(S->ups2, _dh1, 5, _dh3, _di, 1, _dr, _da, _dv1)
     + GET(S->hups2, _dh1, 5, _dh3, _di, 1, _dr, _da, _dv1))
         * GET(popn, _dh1, 5, _dh3, _di, 1, _dr, _da, _dv1);

    //      agein6 = (ups3(:, : , 5, : , 1, : , : , : ) +
    //               hups3(:, : , 5, : , 1, : , : , : )).*popn(:, : , 5, : , 1, : , : , : );

    LOOP_H1 LOOP_H2 LOOP_I LOOP_R LOOP_A LOOP_V1
        GET(agein6, _dh1, _dh2, 1, _di, 1, _dr, _da, _dv1) =
      (GET(S->ups3, _dh1, _dh2, 5, _di, 1, _dr, _da, _dv1)
    + GET(S->hups3, _dh1, _dh2, 5, _di, 1, _dr, _da, _dv1))
        * GET(popn, _dh1, _dh2, 5, _di, 1, _dr, _da, _dv1);

    // Break loops here, as we need agein1..6


//    ageinf(1) = sum(sum(sum(sum(sum(sum(agein1(:, : , : , : , 1, 1, : , : ))))))) ...
//              + sum(sum(sum(sum(sum(sum(agein2(:, : , : , : , 1, 1, : , : ))))))) ...
//              + sum(sum(sum(sum(sum(sum(agein3(:, : , : , : , 1, 1, : , : ))))))) ...
//              + sum(sum(sum(sum(sum(sum(agein4(:, : , : , : , 1, 1, : , : ))))))) ...
//              + sum(sum(sum(sum(sum(sum(agein5(:, : , : , : , 1, 1, : , : ))))))) ...
//              + sum(sum(sum(sum(sum(sum(agein6(:, : , : , : , 1, 1, : , : )))))));

//    ageinf(2) = sum(sum(sum(sum(sum(sum(agein1(:, : , : , : , 1, 2, : , : ))))))) ...
//              + sum(sum(sum(sum(sum(sum(agein2(:, : , : , : , 1, 2, : , : ))))))) ...
//              + sum(sum(sum(sum(sum(sum(agein3(:, : , : , : , 1, 2, : , : ))))))) ...
//              + sum(sum(sum(sum(sum(sum(agein4(:, : , : , : , 1, 2, : , : ))))))) ...
//              + sum(sum(sum(sum(sum(sum(agein5(:, : , : , : , 1, 2, : , : ))))))) ...
//              + sum(sum(sum(sum(sum(sum(agein6(:, : , : , : , 1, 2, : , : )))))));

//    ageinf(3) = sum(sum(sum(sum(sum(sum(agein1(:, : , : , : , 1, 3, : , : ))))))) ...
//              + sum(sum(sum(sum(sum(sum(agein2(:, : , : , : , 1, 3, : , : ))))))) ...
//              + sum(sum(sum(sum(sum(sum(agein3(:, : , : , : , 1, 3, : , : ))))))) ...
//              + sum(sum(sum(sum(sum(sum(agein4(:, : , : , : , 1, 3, : , : ))))))) ...
//              + sum(sum(sum(sum(sum(sum(agein5(:, : , : , : , 1, 3, : , : ))))))) ...
//              + sum(sum(sum(sum(sum(sum(agein6(:, : , : , : , 1, 3, : , : )))))));

//  ageinm(1) = sum(sum(sum(sum(sum(sum(agein1(:, : , : , : , 2, 1, : , : ))))))) ...
//            + sum(sum(sum(sum(sum(sum(agein2(:, : , : , : , 2, 1, : , : ))))))) ...
//            + sum(sum(sum(sum(sum(sum(agein3(:, : , : , : , 2, 1, : , : )))))));

//    ageinm(2) = sum(sum(sum(sum(sum(sum(agein1(:, : , : , : , 2, 2, : , : ))))))) ...
//              + sum(sum(sum(sum(sum(sum(agein2(:, : , : , : , 2, 2, : , : ))))))) ...
//              + sum(sum(sum(sum(sum(sum(agein3(:, : , : , : , 2, 2, : , : )))))));

//    ageinm(3) = sum(sum(sum(sum(sum(sum(agein1(:, : , : , : , 2, 3, : , : ))))))) ...
//              + sum(sum(sum(sum(sum(sum(agein2(:, : , : , : , 2, 3, : , : ))))))) ...
//              + sum(sum(sum(sum(sum(sum(agein3(:, : , : , : , 2, 3, : , : )))))));

    // Re-use the sums of first four indices as much as possible...

    for (int _r = 1; _r <= 3; _r++) { ageinf[_r] = 0; ageinm[_r] = 0; }

    sum_h1_h2_h3_i(agein1, tmpsum_h1_h2_h3_i);
    for (int _r = 1; _r <= 3; _r++) { ageinf[_r] += sum_ageinf(tmpsum_h1_h2_h3_i, 1, _r);
                                      ageinm[_r] += sum_ageinf(tmpsum_h1_h2_h3_i, 2, _r); }

    sum_h1_h2_h3_i(agein2, tmpsum_h1_h2_h3_i);
    for (int _r = 1; _r <= 3; _r++) { ageinf[_r] += sum_ageinf(tmpsum_h1_h2_h3_i, 1, _r);
                                      ageinm[_r] += sum_ageinf(tmpsum_h1_h2_h3_i, 2, _r); }

    sum_h1_h2_h3_i(agein3, tmpsum_h1_h2_h3_i);
    for (int _r = 1; _r <= 3; _r++) { ageinf[_r] += sum_ageinf(tmpsum_h1_h2_h3_i, 1, _r);
                                      ageinm[_r] += sum_ageinf(tmpsum_h1_h2_h3_i, 2, _r); }

    sum_h1_h2_h3_i(agein4, tmpsum_h1_h2_h3_i);
    for (int _r = 1; _r <= 3; _r++) ageinf[_r] += sum_ageinf(tmpsum_h1_h2_h3_i, 1, _r);

    sum_h1_h2_h3_i(agein5, tmpsum_h1_h2_h3_i);
    for (int _r = 1; _r <= 3; _r++) ageinf[_r] += sum_ageinf(tmpsum_h1_h2_h3_i, 1, _r);

    sum_h1_h2_h3_i(agein6, tmpsum_h1_h2_h3_i);
    for (int _r = 1; _r <= 3; _r++) ageinf[_r] += sum_ageinf(tmpsum_h1_h2_h3_i, 1, _r);

    //     %% Disease stages
    //     % Population matrices for calculating changes for each infection per time step

    //       [pdt1, pdt2, pdt3, pdt4, pdt5, pdt6, pdt7] = deal(zeros(h1, h2, h3, i, s, r, a, v1));

    // Move the above line to outside of the iteration loop...

    //     %% HIV - changes in one tstep
    LOOP_H1 LOOP_H2 LOOP_H3 LOOP_S LOOP_R LOOP_A LOOP_V1{

      //           pdt1(:, : , : , 1, : , : , : , : ) = ...
      //       (-foihiv(:, : , : , 1, : , : , : , : ).*popn(:, : , : , 1, : , : , : , : )).*dt; % HIV FOI

      int _index1 = GET_INDEX(_dh1,_dh2,_dh3,1,_ds,_dr,_da,_dv1);
      pdt1[_index1] = (-foihiv[_index1] * popn[_index1]) * dt;


      //       pdt1(:, : , : , 2, : , : , : , : ) = ...
      //       (-nu(:, : , : , 2, : , : , : , : ).*popn(:, : , : , 2, : , : , : , : )... % HIV - related death rate
      //         + foihiv(:, : , : , 1, : , : , : , : ).*popn(:, : , : , 1, : , : , : , : )... % HIV FOI
      //         - eta(:, : , : , 2, : , : , : , : ).*popn(:, : , : , 2, : , : , : , : )).*dt;% duration in acute phase

      int _index2 = GET_INDEX(_dh1, _dh2, _dh3, 2, _ds, _dr, _da, _dv1);
      pdt1[_index2] = (-S->nu[_index2] * popn[_index2] +
                       foihiv[_index1] * popn[_index1] -
                       S->eta[_index2] * popn[_index2]) * dt;

      //           pdt1(:, : , : , 3, : , : , : , : ) = ...
      //         + (-nu(:, : , : , 3, : , : , : , : ).*popn(:, : , : , 3, : , : , : , : )... % HIV - related death rate
      //          + eta(:, : , : , 2, : , : , : , : ).*popn(:, : , : , 2, : , : , : , : )... % acute phase
      //          - tau(:, : , : , 3, : , : , : , : ).*popn(:, : , : , 3, : , : , : , : )... % initiation of ART
      //      + omikron(:, : , : , 4, : , : , : , : ).*popn(:, : , : , 4, : , : , : , : )).*dt;% stopping ART

      int _index3 = GET_INDEX(_dh1, _dh2, _dh3, 3, _ds, _dr, _da, _dv1);
      int _index4 = GET_INDEX(_dh1, _dh2, _dh3, 4, _ds, _dr, _da, _dv1);

      pdt1[_index3] = (-S->nu[_index3] * popn[_index3]
                     + S->eta[_index2] * popn[_index2]
                     - S->tau[_index3] * popn[_index3]
                 + S->omikron[_index4] * popn[_index4]) * dt;

      //       pdt1(:, : , : , 4, : , : , : , : ) = ...
      //       + (-nu(:, : , : , 4, : , : , : , : ).*popn(:, : , : , 4, : , : , : , : )... % HIV - related death rate
      //         + tau(:, : , : , 3, : , : , : , : ).*popn(:, : , : , 3, : , : , : , : )... % initiation of ART
      //         - omikron(:, : , : , 4, : , : , : , : ).*popn(:, : , : , 4, : , : , : , : )).*dt;% stopping AR

      pdt1[_index4] = (-S->nu[_index4] * popn[_index4] +
                       S->tau[_index3] * popn[_index3] -
                   S->omikron[_index4] * popn[_index4]) * dt;

    }

    //       %% HPV 1

#define em p->pp0[7]
#define zeta p->pp0[8]
#define qu p->pp0[9]

    LOOP_H2 LOOP_H3 LOOP_I LOOP_S LOOP_R LOOP_A LOOP_V1 {

      //    pdt2(1, :, : , : , : , : , : , : ) = ... % %16 / 18 HPV FOI
      //                        (-foihpv1(1, :, : , : , : , : , : , : ).*popn(1, :, : , : , : , : , : , : )... % HPV1 FOI
      //           + (1 - em).*qu.*omega1(4, :, : , : , : , : , : , : ).*popn(4, :, : , : , : , : , : , : )... % Regression from CIN2 + / CIN3
      //         + (1 - em).*zeta.*theta1(4, :, : , : , : , : , : , : ).*popn(4, :, : , : , : , : , : , : )... % Screening and treatment
      //               + (1 - em).*sigma1(2, :, : , : , : , : , : , : ).*popn(2, :, : , : , : , : , : , : )... % Natural clearance from I state
      //                         + delta1(3, :, : , : , : , : , : , : ).*popn(3, :, : , : , : , : , : , : )).*dt;

      int _index1 = GET_INDEX(1, _dh2, _dh3, _di, _ds, _dr, _da, _dv1);
      int _index2 = GET_INDEX(2, _dh2, _dh3, _di, _ds, _dr, _da, _dv1);
      int _index3 = GET_INDEX(3, _dh2, _dh3, _di, _ds, _dr, _da, _dv1);
      int _index4 = GET_INDEX(4, _dh2, _dh3, _di, _ds, _dr, _da, _dv1);
      int _index5 = GET_INDEX(5, _dh2, _dh3, _di, _ds, _dr, _da, _dv1);

              pdt2[_index1] = (-foihpv1[_index1] * popn[_index1]
          + (1.0 - em) * qu * S->omega1[_index4] * popn[_index4]
        + (1.0 - em) * zeta * S->theta1[_index4] * popn[_index4]
               + (1.0 - em) * S->sigma1[_index2] * popn[_index2]
                            + S->delta1[_index3] * popn[_index3]) * dt;

      //     pdt2(2, :, : , : , : , : , : , : ) = ... %
      //                 (foihpv1(1, :, : , : , : , : , : , : ).*popn(1, :, : , : , : , : , : , : )... %
      //        +(1 - qu).*omega1(4, :, : , : , : , : , : , : ).*popn(4, :, : , : , : , : , : , : )... %
      //      +(1 - zeta).*theta1(4, :, : , : , : , : , : , : ).*popn(4, :, : , : , : , : , : , : )... %
      //                  -sigma1(2, :, : , : , : , : , : , : ).*popn(2, :, : , : , : , : , : , : )... %
      //                    -psi1(2, :, : , : , : , : , : , : ).*popn(2, :, : , : , : , : , : , : )).*dt;% Progression CIN2 + / CIN3

          pdt2[_index2] = (foihpv1[_index1] * popn[_index1]
          + (1.0 - qu) * S->omega1[_index4] * popn[_index4]
        + (1.0 - zeta) * S->theta1[_index4] * popn[_index4]
                       - S->sigma1[_index2] * popn[_index2]
                         - S->psi1[_index2] * popn[_index2]) * dt;

      //       pdt2(3, :, : , : , : , : , : , : ) = ... %
      //           (em.*qu.*omega1(4, :, : , : , : , : , : , : ).*popn(4, :, : , : , : , : , : , : )... %
      //         +em.*zeta.*theta1(4, :, : , : , : , : , : , : ).*popn(4, :, : , : , : , : , : , : )... %
      //               +em.*sigma1(2, :, : , : , : , : , : , : ).*popn(2, :, : , : , : , : , : , : )... %
      //                   -delta1(3, :, : , : , : , : , : , : ).*popn(3, :, : , : , : , : , : , : )).*dt; %

      pdt2[_index3] = (em * qu * S->omega1[_index4] * popn[_index4]
                   + em * zeta * S->theta1[_index4] * popn[_index4]
                          + em * S->sigma1[_index2] * popn[_index2]
                               - S->delta1[_index3] * popn[_index3]) * dt;

      //       pdt2(4, :, : , : , : , : , : , : ) = ...
      //            (psi1(2, :, : , : , : , : , : , : ).*popn(2, :, : , : , : , : , : , : )...
      //         - omega1(4, :, : , : , : , : , : , : ).*popn(4, :, : , : , : , : , : , : )... % Regression from CIN2 + / CIN3
      //         - theta1(4, :, : , : , : , : , : , : ).*popn(4, :, : , : , : , : , : , : )... % Screening and treatment
      //            - pi1(4, :, : , : , : , : , : , : ).*popn(4, :, : , : , : , : , : , : )).*dt;% Progression to CC

      pdt2[_index4] = (S->psi1[_index2] * popn[_index2]
                   - S->omega1[_index4] * popn[_index4]
                   - S->theta1[_index4] * popn[_index4]
                      - S->pi1[_index4] * popn[_index4]) * dt;


      //       pdt2(5, :, : , : , : , : , : , : ) = ...
      //           (pi1(4, :, : , : , : , : , : , : ).*popn(4, :, : , : , : , : , : , : )...
      //         - ups1(5, :, : , : , : , : , : , : ).*popn(5, :, : , : , : , : , : , : )... % CC deaths
      //        - hups1(5, :, : , : , : , : , : , : ).*popn(5, :, : , : , : , : , : , : )).*dt;% Removal from model due to hysterectomy

      pdt2[_index5] = (S->pi1[_index4] * popn[_index4]
                    - S->ups1[_index5] * popn[_index5]
                    - S->hups1[_index5] * popn[_index5]) * dt;
    }

    //       %% HPV 2

    LOOP_H1 LOOP_H3 LOOP_I LOOP_S LOOP_R LOOP_A LOOP_V1 {

      //          pdt3(:, 1, : , : , : , : , : , : ) = ... % 31 / 33 / 45 / 52 / 58 HPV FOI
      //                      + (-foihpv2(:, 1, : , : , : , : , : , : ).*popn(:, 1, : , : , : , : , : , : )... % HPV1 FOI
      //           + (1 - em).*qu.*omega2(:, 4, : , : , : , : , : , : ).*popn(:, 4, : , : , : , : , : , : )... % Regression from CIN2 + / CIN3
      //         + (1 - em).*zeta.*theta2(:, 4, : , : , : , : , : , : ).*popn(:, 4, : , : , : , : , : , : )... % Screening and treatment
      //               + (1 - em).*sigma2(:, 2, : , : , : , : , : , : ).*popn(:, 2, : , : , : , : , : , : )... % Natural clearance from I state
      //                         + delta2(:, 3, : , : , : , : , : , : ).*popn(:, 3, : , : , : , : , : , : )).*dt;

      int _index1 = GET_INDEX(_dh1, 1, _dh3, _di, _ds, _dr, _da, _dv1);
      int _index2 = GET_INDEX(_dh1, 2, _dh3, _di, _ds, _dr, _da, _dv1);
      int _index3 = GET_INDEX(_dh1, 3, _dh3, _di, _ds, _dr, _da, _dv1);
      int _index4 = GET_INDEX(_dh1, 4, _dh3, _di, _ds, _dr, _da, _dv1);
      int _index5 = GET_INDEX(_dh1, 5, _dh3, _di, _ds, _dr, _da, _dv1);

            pdt3[_index1] = (-foihpv2[_index1] * popn[_index1]
        + (1.0 - em) * qu * S->omega2[_index4] * popn[_index4]
      + (1.0 - em) * zeta * S->theta2[_index4] * popn[_index4]
             + (1.0 - em) * S->sigma2[_index2] * popn[_index2]
                          + S->delta2[_index3] * popn[_index3]) * dt;

      //     pdt3(:, 2, : , : , : , : , : , : ) = ... %
      //                +(+foihpv2(:, 1, : , : , : , : , : , : ).*popn(:, 1, : , : , : , : , : , : )... %
      //         +(1 - qu).*omega2(:, 4, : , : , : , : , : , : ).*popn(:, 4, : , : , : , : , : , : )... %
      //       +(1 - zeta).*theta2(:, 4, : , : , : , : , : , : ).*popn(:, 4, : , : , : , : , : , : )... %
      //                   -sigma2(:, 2, : , : , : , : , : , : ).*popn(:, 2, : , : , : , : , : , : )... %
      //                     -psi2(:, 2, : , : , : , : , : , : ).*popn(:, 2, : , : , : , : , : , : )).*dt;% Progression CIN2 + / CIN3

        pdt3[_index2] = (foihpv2[_index1] * popn[_index1]
        + (1.0 - qu) * S->omega2[_index4] * popn[_index4]
      + (1.0 - zeta) * S->theta2[_index4] * popn[_index4]
                     - S->sigma2[_index2] * popn[_index2]
                       - S->psi2[_index2] * popn[_index2]) * dt;

      //       pdt3(:, 3, : , : , : , : , : , : ) = ... %
      //       +(+em.*qu.*omega2(:, 4, : , : , : , : , : , : ).*popn(:, 4, : , : , : , : , : , : )... %
      //         +em.*zeta.*theta2(:, 4, : , : , : , : , : , : ).*popn(:, 4, : , : , : , : , : , : )... %
      //               +em.*sigma2(:, 2, : , : , : , : , : , : ).*popn(:, 2, : , : , : , : , : , : )... %
      //                   -delta2(:, 3, : , : , : , : , : , : ).*popn(:, 3, : , : , : , : , : , : )).*dt;%

      pdt3[_index3] = (em * qu * S->omega2[_index4] * popn[_index4]
                   + em * zeta * S->theta2[_index4] * popn[_index4]
                          + em * S->sigma2[_index2] * popn[_index2]
                               - S->delta2[_index3] * popn[_index3]) * dt;

      //       pdt3(:, 4, : , : , : , : , : , : ) = ...
      //       + (psi2(:, 2, : , : , : , : , : , : ).*popn(:, 2, : , : , : , : , : , : )...
      //         - omega2(:, 4, : , : , : , : , : , : ).*popn(:, 4, : , : , : , : , : , : )...
      //         - theta2(:, 4, : , : , : , : , : , : ).*popn(:, 4, : , : , : , : , : , : )...
      //         - pi2(:, 4, : , : , : , : , : , : ).*popn(:, 4, : , : , : , : , : , : )).*dt;% Progression to CC

      pdt3[_index4] = (S->psi2[_index2] * popn[_index2]
                   - S->omega2[_index4] * popn[_index4]
                   - S->theta2[_index4] * popn[_index4]
                      - S->pi2[_index4] * popn[_index4]) * dt;

      //       pdt3(:, 5, : , : , : , : , : , : ) = ...
//       + (+pi2(:, 4, : , : , : , : , : , : ).*popn(:, 4, : , : , : , : , : , : )...
//         - ups2(:, 5, : , : , : , : , : , : ).*popn(:, 5, : , : , : , : , : , : )... % CC deaths
//         - hups2(:, 5, : , : , : , : , : , : ).*popn(:, 5, : , : , : , : , : , : )).*dt;% Removal from model due to hysterectomy

      pdt3[_index5] = (S->pi2[_index4] * popn[_index4]
                     - S->ups2[_index5] * popn[_index5]
                     - S->hups2[_index5] * popn[_index5]) * dt;
    }

    //       %% HPV 3
    LOOP_H1 LOOP_H2 LOOP_I LOOP_S LOOP_R LOOP_A LOOP_V1{

      int _index1 = GET_INDEX(_dh1, _dh2, 1, _di, _ds, _dr, _da, _dv1);
      int _index2 = GET_INDEX(_dh1, _dh2, 2, _di, _ds, _dr, _da, _dv1);
      int _index3 = GET_INDEX(_dh1, _dh2, 3, _di, _ds, _dr, _da, _dv1);
      int _index4 = GET_INDEX(_dh1, _dh2, 4, _di, _ds, _dr, _da, _dv1);
      int _index5 = GET_INDEX(_dh1, _dh2, 5, _di, _ds, _dr, _da, _dv1);
      double m_em = 1.0 - em;

      // pdt4(:, : , 1, : , : , : , : , : ) = ... % nvt HR - HPV
      //                      + (-foihpv3(:, : , 1, : , : , : , : , : ).*popn(:, : , 1, : , : , : , : , : )... % HPV1 FOI
      //           + (1 - em).*qu.*omega3(:, : , 4, : , : , : , : , : ).*popn(:, : , 4, : , : , : , : , : )... % Regression from CIN2 + / CIN3
      //         + (1 - em).*zeta.*theta3(:, : , 4, : , : , : , : , : ).*popn(:, : , 4, : , : , : , : , : )... % Screening and treatment
      //               + (1 - em).*sigma3(:, : , 2, : , : , : , : , : ).*popn(:, : , 2, : , : , : , : , : )... % Natural clearance from I state
      //                         + delta3(:, : , 3, : , : , : , : , : ).*popn(:, : , 3, : , : , : , : , : )).*dt;

       pdt4[_index1] = (-foihpv3[_index1] * popn[_index1]
         + (m_em)*qu * S->omega3[_index4] * popn[_index4]
       + (m_em)*zeta * S->theta3[_index4] * popn[_index4]
              + (m_em)*S->sigma3[_index2] * popn[_index2]
                     + S->delta3[_index3] * popn[_index3]) * dt;

       //     pdt4(:, : , 2, : , : , : , : , : ) = ... %
       //                   +(+foihpv3(:, : , 1, : , : , : , : , : ).*popn(:, : , 1, : , : , : , : , : )... %
       //            +(1 - qu).*omega3(:, : , 4, : , : , : , : , : ).*popn(:, : , 4, : , : , : , : , : )...
       //         + (1 - zeta).*theta3(:, : , 4, : , : , : , : , : ).*popn(:, : , 4, : , : , : , : , : )... %
       //                      -sigma3(:, : , 2, : , : , : , : , : ).*popn(:, : , 2, : , : , : , : , : )... %
       //                        -psi3(:, : , 2, : , : , : , : , : ).*popn(:, : , 2, : , : , : , : , : )).*dt;% Progression CIN2 + / CIN3

         pdt4[_index2] = (foihpv3[_index1] * popn[_index1]
         + (1.0 - qu) * S->omega3[_index4] * popn[_index4]
       + (1.0 - zeta) * S->theta3[_index4] * popn[_index4]
                      - S->sigma3[_index2] * popn[_index2]
                        - S->psi3[_index2] * popn[_index2]) * dt;

         //                      pdt4(:, : , 3, : , : , : , : , : ) = ... %
         //       +(+em.*qu.*omega3(:, : , 4, : , : , : , : , : ).*popn(:, : , 4, : , : , : , : , : )... %
         //         +em.*zeta.*theta3(:, : , 4, : , : , : , : , : ).*popn(:, : , 4, : , : , : , : , : )... %
         //               +em.*sigma3(:, : , 2, : , : , : , : , : ).*popn(:, : , 2, : , : , : , : , : )... %
         //                   -delta3(:, : , 3, : , : , : , : , : ).*popn(:, : , 3, : , : , : , : , : )).*dt;%

        pdt4[_index3] = (em * qu * S->omega3[_index4] * popn[_index4]
                     + em * zeta * S->theta3[_index4] * popn[_index4]
                            + em * S->sigma3[_index2] * popn[_index2]
                                 - S->delta3[_index3] * popn[_index3]) * dt;

        //             pdt4(:, : , 4, : , : , : , : , : ) = ...
        //          + (psi3(:, : , 2, : , : , : , : , : ).*popn(:, : , 2, : , : , : , : , : )...
        //         - omega3(:, : , 4, : , : , : , : , : ).*popn(:, : , 4, : , : , : , : , : )...
        //         - theta3(:, : , 4, : , : , : , : , : ).*popn(:, : , 4, : , : , : , : , : )...
        //            - pi3(:, : , 4, : , : , : , : , : ).*popn(:, : , 4, : , : , : , : , : )).*dt;% Progression to CC

        pdt4[_index4] = (S->psi3[_index2] * popn[_index2]
                     - S->omega3[_index4] * popn[_index4]
                     - S->theta3[_index4] * popn[_index4]
                        - S->pi3[_index4] * popn[_index4]) * dt;

        //            pdt4(:, : , 5, : , : , : , : , : ) = ...
        //         + (+pi3(:, : , 4, : , : , : , : , : ).*popn(:, : , 4, : , : , : , : , : )...
        //          - ups3(:, : , 5, : , : , : , : , : ).*popn(:, : , 5, : , : , : , : , : )... % CC deaths
        //         - hups3(:, : , 5, : , : , : , : , : ).*popn(:, : , 5, : , : , : , : , : )).*dt;% Removal from model due to hysterectomy

        pdt4[_index5] = (S->pi3[_index4] * popn[_index4]
                      - S->ups3[_index5] * popn[_index5]
                     - S->hups3[_index5] * popn[_index5]) * dt;
    }

    //       %% AGING + background mortality

    LOOP_H1 LOOP_H2 LOOP_H3 LOOP_I LOOP_S LOOP_R LOOP_V1{

      int _index = GET_INDEX(_dh1, _dh2, _dh3, _di, _ds, _dr, 1, _dv1);
      int _index_prev;

      // %% AGING + background mortality
      // pdt5(:,:,:,:,:,:,1,:) =  + (-alpha(:,:,:,:,:,:,1,:).*popn(:,:,:,:,:,:,1,:)     - mu(:,:,:,:,:,:,1,:).*popn(:,:,:,:,:,:,1,:)).*dt;
      // pdt5(:,:,:,:,:,:,2,:) =  + (+alpha(:,:,:,:,:,:,1,:).*popn(:,:,:,:,:,:,1,:)  - alpha(:,:,:,:,:,:,2,:).*popn(:,:,:,:,:,:,2,:)- mu(:,:,:,:,:,:,2,:).*popn(:,:,:,:,:,:,2,:)).*dt;
      // pdt5(:,:,:,:,:,:,3,:) =  + (+alpha(:,:,:,:,:,:,2,:).*popn(:,:,:,:,:,:,2,:)  - alpha(:,:,:,:,:,:,3,:).*popn(:,:,:,:,:,:,3,:)- mu(:,:,:,:,:,:,3,:).*popn(:,:,:,:,:,:,3,:)).*dt;
      // pdt5(:,:,:,:,:,:,4,:) =  + (+alpha(:,:,:,:,:,:,3,:).*popn(:,:,:,:,:,:,3,:)  - alpha(:,:,:,:,:,:,4,:).*popn(:,:,:,:,:,:,4,:)- mu(:,:,:,:,:,:,4,:).*popn(:,:,:,:,:,:,4,:)).*dt;
      // pdt5(:,:,:,:,:,:, 5,:) = + (+alpha(:,:,:,:,:,:,4,:).*popn(:,:,:,:,:,:,4,:)  - alpha(:,:,:,:,:,:, 5,:).*popn(:,:,:,:,:,:,5,:) - mu(:,:,:,:,:,:,5,:).*popn(:,:,:,:,:,:,5,:)).*dt;
      // pdt5(:,:,:,:,:,:, 6,:) = + (+alpha(:,:,:,:,:,:,5,:).*popn(:,:,:,:,:,:,5,:)  - alpha(:,:,:,:,:,:, 6,:).*popn(:,:,:,:,:,:,6,:) - mu(:,:,:,:,:,:,6,:).*popn(:,:,:,:,:,:,6,:)).*dt;
      // pdt5(:,:,:,:,:,:, 7,:) = + (+alpha(:,:,:,:,:,:,6,:).*popn(:,:,:,:,:,:,6,:)  - alpha(:,:,:,:,:,:, 7,:).*popn(:,:,:,:,:,:,7,:) - mu(:,:,:,:,:,:,7,:).*popn(:,:,:,:,:,:,7,:)).*dt;
      // pdt5(:,:,:,:,:,:, 8,:) = + (+alpha(:,:,:,:,:,:,7,:).*popn(:,:,:,:,:,:,7,:)  - alpha(:,:,:,:,:,:, 8,:).*popn(:,:,:,:,:,:,8,:) - mu(:,:,:,:,:,:,8,:).*popn(:,:,:,:,:,:,8,:)).*dt;
      // pdt5(:,:,:,:,:,:, 9,:) = + (+alpha(:,:,:,:,:,:,8,:).*popn(:,:,:,:,:,:,8,:)  - alpha(:,:,:,:,:,:, 9,:).*popn(:,:,:,:,:,:,9,:) - mu(:,:,:,:,:,:,9,:).*popn(:,:,:,:,:,:,9,:)).*dt;
      // pdt5(:,:,:,:,:,:, 10,:) = +(+alpha(:,:,:,:,:,:,9,:).*popn(:,:,:,:,:,:,9,:)  - alpha(:,:,:,:,:,:, 10,:).*popn(:,:,:,:,:,:,10,:) - mu(:,:,:,:,:,:,10,:).*popn(:,:,:,:,:,:,10,:)).*dt;
      // pdt5(:,:,:,:,:,:, 11,:) = +(+alpha(:,:,:,:,:,:,10,:).*popn(:,:,:,:,:,:,10,:)- alpha(:,:,:,:,:,:, 11,:).*popn(:,:,:,:,:,:,11,:) - mu(:,:,:,:,:,:,11,:).*popn(:,:,:,:,:,:,11,:)).*dt;
      // pdt5(:,:,:,:,:,:, 12,:) = +(+alpha(:,:,:,:,:,:,11,:).*popn(:,:,:,:,:,:,11,:)- alpha(:,:,:,:,:,:, 12,:).*popn(:,:,:,:,:,:,12,:) - mu(:,:,:,:,:,:,12,:).*popn(:,:,:,:,:,:,12,:)).*dt;
      // pdt5(:,:,:,:,:,:, 13,:) = +(+alpha(:,:,:,:,:,:,12,:).*popn(:,:,:,:,:,:,12,:)- alpha(:,:,:,:,:,:, 13,:).*popn(:,:,:,:,:,:,13,:) - mu(:,:,:,:,:,:,13,:).*popn(:,:,:,:,:,:,13,:)).*dt;

      pdt5[_index] = (-S->alpha[_index] * popn[_index] - S->mu[_index] * popn[_index]) * dt;

      for (int _da = 2; _da <= _SIZE_A; _da++) {
         _index_prev = _index;
        _index = GET_INDEX(_dh1, _dh2, _dh3, _di, _ds, _dr, _da, _dv1);
        pdt5[_index] = (S->alpha[_index_prev] * popn[_index_prev] - S->alpha[_index] * popn[_index] - S->mu[_index] * popn[_index])* dt;
      }
    }
    //     %% POP GROWTH -- including the proportion of HIV infected people at start of 9 yos

    //     pdt6(1, 1, 1, 1, 1, 1, 1, 1) = (1 - (perin(n) + perin2(n))).*gr_tv(n).*ageinf(1).*dt;
    //     pdt6(1, 1, 1, 3, 1, 1, 1, 1) = perin(n).*gr_tv(n).*ageinf(1).*dt;
    //     pdt6(1, 1, 1, 4, 1, 1, 1, 1) = perin2(n).*gr_tv(n).*ageinf(1).*dt;

    SET(pdt6, 1, 1, 1, 1, 1, 1, 1, 1, (1 - (perin[n] + perin2[n]))* gr_tv[n] * ageinf[1] * dt);
    SET(pdt6, 1, 1, 1, 3, 1, 1, 1, 1, perin[n] * gr_tv[n] * ageinf[1] * dt);
    SET(pdt6, 1, 1, 1, 4, 1, 1, 1, 1, perin2[n] * gr_tv[n] * ageinf[1] * dt);

    //     pdt6(1, 1, 1, 1, 1, 2, 1, 1) = (1 - (perin(n) + perin2(n))).*gr_tv(n).*ageinf(2).*dt;
    //     pdt6(1, 1, 1, 3, 1, 2, 1, 1) = perin(n).*gr_tv(n).*ageinf(2).*dt;
    //     pdt6(1, 1, 1, 4, 1, 2, 1, 1) = perin2(n).*gr_tv(n).*ageinf(2).*dt;

    SET(pdt6, 1, 1, 1, 1, 1, 2, 1, 1, (1 - (perin[n] + perin2[n]))* gr_tv[n] * ageinf[2] * dt);
    SET(pdt6, 1, 1, 1, 3, 1, 2, 1, 1, perin[n] * gr_tv[n] * ageinf[2] * dt);
    SET(pdt6, 1, 1, 1, 4, 1, 2, 1, 1, perin2[n] * gr_tv[n] * ageinf[2] * dt);

    //     pdt6(1, 1, 1, 1, 1, 3, 1, 1) = (1 - (perin(n) + perin2(n))).*gr_tv(n).*ageinf(3).*dt;
    //     pdt6(1, 1, 1, 3, 1, 3, 1, 1) = perin(n).*gr_tv(n).*ageinf(3).*dt;
    //     pdt6(1, 1, 1, 4, 1, 3, 1, 1) = perin2(n).*gr_tv(n).*ageinf(3).*dt;

    SET(pdt6, 1, 1, 1, 1, 1, 3, 1, 1, (1 - (perin[n] + perin2[n]))* gr_tv[n] * ageinf[3] * dt);
    SET(pdt6, 1, 1, 1, 3, 1, 3, 1, 1, perin[n] * gr_tv[n] * ageinf[3] * dt);
    SET(pdt6, 1, 1, 1, 4, 1, 3, 1, 1, perin2[n] * gr_tv[n] * ageinf[3] * dt);

    //     pdt6(1, 1, 1, 1, 2, 1, 1, 1) = (1 - (perin(n) + perin2(n))).*gr_tv(n).*ageinm(1).*dt;
    //     pdt6(1, 1, 1, 3, 2, 1, 1, 1) = perin(n).*gr_tv(n).*ageinm(1).*dt;
    //     pdt6(1, 1, 1, 4, 2, 1, 1, 1) = perin2(n).*gr_tv(n).*ageinm(1).*dt;

    SET(pdt6, 1, 1, 1, 1, 2, 1, 1, 1, (1 - (perin[n] + perin2[n]))* gr_tv[n] * ageinm[1] * dt);
    SET(pdt6, 1, 1, 1, 3, 2, 1, 1, 1, perin[n] * gr_tv[n] * ageinm[1] * dt);
    SET(pdt6, 1, 1, 1, 4, 2, 1, 1, 1, perin2[n] * gr_tv[n] * ageinm[1] * dt);

    //     pdt6(1, 1, 1, 1, 2, 2, 1, 1) = (1 - (perin(n) + perin2(n))).*gr_tv(n).*ageinm(2).*dt;
    //     pdt6(1, 1, 1, 3, 2, 2, 1, 1) = perin(n).*gr_tv(n).*ageinm(2).*dt;
    //     pdt6(1, 1, 1, 4, 2, 2, 1, 1) = perin2(n).*gr_tv(n).*ageinm(2).*dt;

    SET(pdt6, 1, 1, 1, 1, 2, 2, 1, 1, (1 - (perin[n] + perin2[n]))* gr_tv[n] * ageinm[2] * dt);
    SET(pdt6, 1, 1, 1, 3, 2, 2, 1, 1, perin[n] * gr_tv[n] * ageinm[2] * dt);
    SET(pdt6, 1, 1, 1, 4, 2, 2, 1, 1, perin2[n] * gr_tv[n] * ageinm[2] * dt);

    //     pdt6(1, 1, 1, 1, 2, 3, 1, 1) = (1 - (perin(n) + perin2(n))).*gr_tv(n).*ageinm(3).*dt;
    //     pdt6(1, 1, 1, 3, 2, 3, 1, 1) = perin(n).*gr_tv(n).*ageinm(3).*dt;
    //     pdt6(1, 1, 1, 4, 2, 3, 1, 1) = perin2(n).*gr_tv(n).*ageinm(3).*dt;

    SET(pdt6, 1, 1, 1, 1, 2, 3, 1, 1, (1 - (perin[n] + perin2[n]))* gr_tv[n] * ageinm[3] * dt);
    SET(pdt6, 1, 1, 1, 3, 2, 3, 1, 1, perin[n] * gr_tv[n] * ageinm[3] * dt);
    SET(pdt6, 1, 1, 1, 4, 2, 3, 1, 1, perin2[n] * gr_tv[n] * ageinm[3] * dt);

    //     %% VACCINATION - here popdt as adding changes to previous
    LOOP_H1 LOOP_H2 LOOP_H3 LOOP_I LOOP_S LOOP_R LOOP_A{

      //       pdt7(:, : , : , : , : , : , : , 1) = ...
      //       (-vxx(:, : , : , : , : , : , : , 1).*popn(:, : , : , : , : , : , : , 1) ...
      //         + vxxwane(:, : , : , : , : , : , : , 2).*popn(:, : , : , : , : , : , : , 2)).*dt;

      int _index1 = GET_INDEX(_dh1, _dh2, _dh3, _di, _ds, _dr, _da, 1);
      int _index2 = GET_INDEX(_dh1, _dh2, _dh3, _di, _ds, _dr, _da, 2);


      pdt7[_index1] = (-vxx[_index1] * popn[_index1] + S->vxxwane[_index2] * popn[_index2]) * dt;

            //pdt7(:, : , : , : , : , : , : , 2) = ...
   //       (+vxx(:, : , : , : , : , : , : , 1).*popn(:, : , : , : , : , : , : , 1)...
   //         - vxxwane(:, : , : , : , : , : , : , 2).*popn(:, : , : , : , : , : , : , 2)).*dt;

      pdt7[_index2] = (vxx[_index1] * popn[_index1] - S->vxxwane[_index2] * popn[_index2]) * dt;

    }

    if (n > 1) {

      //hpv1cum_popt(n, :, : , : ) = hpv1cum_popt(n - 1, :, : , : ) + reshape(sum(sum(sum(sum(foihpv1(1, :, : , : , : , : , : , : ).*popn(1, :, : , : , : , : , : , : ).*dt, 2), 3), 4), 6), 1, s, a, v1);
      //hpv2cum_popt(n, :, : , : ) = hpv2cum_popt(n - 1, :, : , : ) + reshape(sum(sum(sum(sum(foihpv2(:, 1, : , : , : , : , : , : ).*popn(:, 1, : , : , : , : , : , : ).*dt, 1), 3), 4), 6), 1, s, a, v1);
      //hpv3cum_popt(n, :, : , : ) = hpv3cum_popt(n - 1, :, : , : ) + reshape(sum(sum(sum(sum(foihpv3(:, : , 1, : , : , : , : , : ).*popn(:, : , 1, : , : , : , : , : ).*dt, 1), 2), 4), 6), 1, s, a, v1);

      //hpv1hncum_popt(n, :, : , : ) = hpv1hncum_popt(n - 1, :, : , : ) + reshape(sum(sum(sum(sum(foihpv1(1, :, : , 1, : , : , : , : ).*popn(1, :, : , 1, : , : , : , : ).*dt, 2), 3), 4), 6), 1, s, a, v1);
      //hpv2hncum_popt(n, :, : , : ) = hpv2hncum_popt(n - 1, :, : , : ) + reshape(sum(sum(sum(sum(foihpv2(:, 1, : , 1, : , : , : , : ).*popn(:, 1, : , 1, : , : , : , : ).*dt, 1), 3), 4), 6), 1, s, a, v1);
      //hpv3hncum_popt(n, :, : , : ) = hpv3hncum_popt(n - 1, :, : , : ) + reshape(sum(sum(sum(sum(foihpv3(:, : , 1, 1, : , : , : , : ).*popn(:, : , 1, 1, : , : , : , : ).*dt, 1), 2), 4), 6), 1, s, a, v1);

      sum_prod_2346_578(spare_h2h3ir, foihpv1, popn, hpv1cum_popt, n, 1, _SIZE_I, dt);
      sum_prod_1346_578(spare_h1h3ir, foihpv2, popn, hpv2cum_popt, n, 1, _SIZE_I, dt);
      sum_prod_1246_578(spare_h1h2ir, foihpv3, popn, hpv3cum_popt, n, 1, _SIZE_I, dt);

      sum_prod_2346_578(spare_h2h3ir, foihpv1, popn, hpv1hncum_popt, n, 1, 1, dt);
      sum_prod_1346_578(spare_h1h3ir, foihpv2, popn, hpv2hncum_popt, n, 1, 1, dt);
      sum_prod_1246_578(spare_h1h2ir, foihpv3, popn, hpv3hncum_popt, n, 1, 1, dt);

      
      //% HIV incidence
      //hivinc_pop_cumt(n, :, : , : ) = hivinc_pop_cumt(n - 1, :, : , : ) + reshape(sum(sum(sum(sum(foihiv(:, : , : , 1, : , : , : , : ).*popn(:, : , : , 1, : , : , : , : ).*dt, 1), 2), 3), 6), 1, s, a, v1);

      sum_prod_h1_h2_h3(foihiv, popn, tmpsum_h1_h2_h3);
      LOOP_S LOOP_A LOOP_V1 {
        double sum=0;
        for (int _di=2; _di <= _SIZE_I; _di++) LOOP_R 
          GET(tmpsum_h1_h2_h3, 1, 1, 1, 1, _ds, _dr, _da, _dv1) += GET(tmpsum_h1_h2_h3, 1, 1, 1, _di, _ds, _dr, _da, _dv1);
        LOOP_R sum+=GET(tmpsum_h1_h2_h3, 1, 1, 1, 1, _ds, _dr, _da, _dv1);

        GET_4D_NSAV1(hivinc_pop_cumt, n, _ds, _da, _dv1) = GET_4D_NSAV1(hivinc_pop_cumt, n - 1, _ds, _da, _dv1) + (sum * dt);
      }

      LOOP_I LOOP_S LOOP_R LOOP_A LOOP_V1 {
        LOOP_H2 LOOP_H3 {

        // ccinc_pop(4, :, : , : , : , : , : , : ) = pi1(4, :, : , : , : , : , : , : ).*popn(4, :, : , : , : , : , : , : ).*dt;
        // cin_pop(2, :, : , : , : , : , : , : ) = psi1(2, :, : , : , : , : , : , : ).*popn(2, :, : , : , : , : , : , : ).*dt;

                SET(ccinc_pop, 4, _dh2, _dh3, _di, _ds, _dr, _da, _dv1,
                   GET(S->pi1, 4, _dh2, _dh3, _di, _ds, _dr, _da, _dv1)
                   * GET(popn, 4, _dh2, _dh3, _di, _ds, _dr, _da, _dv1) * dt);

                SET(cin_pop, 2, _dh2, _dh3, _di, _ds, _dr, _da, _dv1,
                GET(S->psi1, 2, _dh2, _dh3, _di, _ds, _dr, _da, _dv1)
                 * GET(popn, 2, _dh2, _dh3, _di, _ds, _dr, _da, _dv1) * dt);
        }

        LOOP_H1 LOOP_H3 {

        //    ccinc_pop(:, 4, : , : , : , : , : , : ) = pi2(:, 4, : , : , : , : , : , : ).*popn(:, 4, : , : , : , : , : , : ).*dt;
        //      cin_pop(:, 2, : , : , : , : , : , : ) = psi2(:, 2, : , : , : , : , : , : ).*popn(:, 2, : , : , : , : , : , : ).*dt;

                SET(ccinc_pop, _dh1, 4, _dh3, _di, _ds, _dr, _da, _dv1,
                   GET(S->pi2, _dh1, 4, _dh3, _di, _ds, _dr, _da, _dv1) *
                    GET(popn,  _dh1, 4, _dh3, _di, _ds, _dr, _da, _dv1) * dt);

                SET(cin_pop,   _dh1, 2, _dh3, _di, _ds, _dr, _da, _dv1,
                  GET(S->psi2, _dh1, 2, _dh3, _di, _ds, _dr, _da, _dv1) *
                    GET(popn,  _dh1, 2, _dh3, _di, _ds, _dr, _da, _dv1) * dt);
        }

        LOOP_H1 LOOP_H2 {


        //    ccinc_pop(:, : , 4, : , : , : , : , : ) = pi3(:, : , 4, : , : , : , : , : ).*popn(:, : , 4, : , : , : , : , : ).*dt;
        //      cin_pop(:, : , 2, : , : , : , : , : ) = psi3(:, : , 2, : , : , : , : , : ).*popn(:, : , 2, : , : , : , : , : ).*dt;

                SET(ccinc_pop,_dh1,_dh2,4,_di,_ds,_dr,_da,_dv1,
                  GET(S->pi3,_dh1,_dh2,4,_di,_ds,_dr,_da,_dv1) *
                    GET(popn, _dh1, _dh2, 4, _di, _ds, _dr, _da, _dv1) * dt);

                SET(cin_pop, _dh1, _dh2, 2, _di, _ds, _dr, _da, _dv1,
                  GET(S->psi3, _dh1, _dh2, 2, _di, _ds, _dr, _da, _dv1)*
                    GET(popn, _dh1, _dh2, 2, _di, _ds, _dr, _da, _dv1)* dt);
        }
      }
      //    % CC incidence

      sum_h1_h2(ccinc_pop, tmpsum_h1_h2);

      //  	% CIN2+ incidence // THIS NOW IN THE NEXT SECTION AFTER CC INCIDENCE (CALCULATED WRONG!)
      //  sum_h1_h2(cin_pop, tmpsum_h1_h2_h3);
      // sum_h1_h2_h3(cin_pop, tmpsum_h1_h2_h3);


      //                ccinc_nvtft_neg(n, :) = ccinc_nvtft_neg(n - 1, :) + reshape(sum(sum(sum(sum(sum(sum(ccinc_pop(:, : , 4, 1,     1, : , : , : ), 1), 2), 3), 4), 6), 8), 1, a);
      //                ccinc_nvtft_pos(n, :) = ccinc_nvtft_pos(n - 1, :) + reshape(sum(sum(sum(sum(sum(sum(ccinc_pop(:, : , 4, 2 : 3, 1, : , : , : ), 1), 2), 3), 4), 6), 8), 1, a);
      //                ccinc_nvtft_art(n, :) = ccinc_nvtft_art(n - 1, :) + reshape(sum(sum(sum(sum(sum(sum(ccinc_pop(:, : , 4, 4,     1, : , : , : ), 1), 2), 3), 4), 6), 8), 1, a);

      sum_12_4568_h3is(spare, tmpsum_h1_h2, 4, 4, 1, 1, 1, 1, sum1_a);

      int n_index = GET_2D_NA_INDEX(n, 1) - 1;
      int nm1_index = GET_2D_NA_INDEX(n - 1, 1) - 1;

      LOOP_A ccinc_nvtft_neg[n_index + _da] = ccinc_nvtft_neg[nm1_index + _da] + sum1_a[_da];
      sum_12_4568_h3is(spare, tmpsum_h1_h2, 4, 4, 2, 3, 1, 1, sum1_a);
      LOOP_A ccinc_nvtft_pos[n_index + _da] = ccinc_nvtft_pos[nm1_index + _da] + sum1_a[_da];
      sum_12_4568_h3is(spare, tmpsum_h1_h2, 4, 4, 4, 4, 1, 1, sum1_a);
      LOOP_A ccinc_nvtft_art[n_index + _da] = ccinc_nvtft_art[nm1_index + _da] + sum1_a[_da];

      sum_h3(tmpsum_h1_h2, tmpsum_h1_h2_h3);

      //      % CC inicdence nvt HPV
      //      % CIN2 + incidence (move this earlier to re-use cached h1+h2+h3


      	// HERE IS WHAT I DID WRONG IN THE MATLAB CODE
        //      ccinc_ft_neg(n, :) = ccinc_ft_neg(n - 1, :) + reshape(sum(sum(sum(sum(sum(sum(ccinc_pop(:, : , : , 1,     1, : , : , : ), 1), 2), 3), 4), 6), 8), 1, a);
        //      cin_ft_neg(n, :)   = cin_ft_neg(n - 1, :)   + reshape(sum(sum(sum(sum(sum(sum(ccinc_pop(:, : , : , 1,     1, : , : , : ), 1), 2), 3), 4), 6), 8), 1, a);
        //      ccinc_ft_pos(n, :) = ccinc_ft_pos(n - 1, :) + reshape(sum(sum(sum(sum(sum(sum(ccinc_pop(:, : , : , 2 : 3, 1, : , : , : ), 1), 2), 3), 4), 6), 8), 1, a);
        //      cin_ft_pos(n, :) =   cin_ft_pos(n - 1, :)   + reshape(sum(sum(sum(sum(sum(sum(ccinc_pop(:, : , : , 2 : 3, 1, : , : , : ), 1), 2), 3), 4), 6), 8), 1, a);
        //      ccinc_ft_art(n, :) = ccinc_ft_art(n - 1, :) + reshape(sum(sum(sum(sum(sum(sum(ccinc_pop(:, : , : , 4,     1, : , : , : ), 1), 2), 3), 4), 6), 8), 1, a);
        //      cin_ft_art(n, :) =   cin_ft_art(n - 1, :)   + reshape(sum(sum(sum(sum(sum(sum(ccinc_pop(:, : , : , 4, 1, : , : , : ), 1), 2), 3), 4), 6), 8), 1, a);

      sum_123_4568_is(spare, tmpsum_h1_h2_h3, 1, 1, 1, 1, sum1_a);
      LOOP_A ccinc_ft_neg[n_index + _da] = ccinc_ft_neg[nm1_index + _da] + sum1_a[_da];
    //  LOOP_A cin_ft_neg[n_index + _da] = cin_ft_neg[nm1_index + _da] + sum1_a[_da];

      sum_123_4568_is(spare, tmpsum_h1_h2_h3, 2, 3, 1, 1, sum1_a);
      LOOP_A ccinc_ft_pos[n_index + _da] = ccinc_ft_pos[nm1_index + _da] + sum1_a[_da];
   //   LOOP_A cin_ft_pos[n_index + _da] = cin_ft_pos[nm1_index + _da] + sum1_a[_da];

      sum_123_4568_is(spare, tmpsum_h1_h2_h3, 4, 4, 1, 1, sum1_a);
      LOOP_A ccinc_ft_art[n_index + _da] = ccinc_ft_art[nm1_index + _da] + sum1_a[_da];
   //   LOOP_A cin_ft_art[n_index + _da] = cin_ft_art[nm1_index + _da] + sum1_a[_da];

      sum_i(tmpsum_h1_h2_h3, tmpsum_h1_h2_h3_i);

      //      ccinc_ft(n, :) = ccinc_ft(n - 1, :) + reshape(sum(sum(sum(sum(sum(sum(ccinc_pop(:, : , : , : , 1, : , : , : ), 1), 2), 3), 4), 6), 8), 1, a);
      //        cin_ft(n, :) = cin_ft(n - 1, :)   + reshape(sum(sum(sum(sum(sum(sum(ccinc_pop(:, : , : , : , 1, : , : , : ), 1), 2), 3), 4), 6), 8), 1, a);

      sum_1234_68_s(spare, tmpsum_h1_h2_h3_i, 1, sum1_a);
      LOOP_A ccinc_ft[n_index + _da] = ccinc_ft[nm1_index + _da] + sum1_a[_da];
     // LOOP_A cin_ft[n_index + _da] = cin_ft[nm1_index + _da] + sum1_a[_da];

        //      ccinc_nvtft(n, :) = ccinc_nvtft(n - 1, :) +
        //reshape(sum(sum(sum(sum(sum(sum(ccinc_pop(:, : , 4, : , : , : , : , : ), 1), 2), 4), 5), 6), 8), 1, a);

      sum_12_4568_h3is(spare, tmpsum_h1_h2, 4, 4, 1, dims[DIM_i], 1, dims[DIM_s], sum1_a);

      LOOP_A ccinc_nvtft[n_index + _da] = ccinc_nvtft[nm1_index + _da] + sum1_a[_da];


      //////// NEW ////////

      //  	% CIN2+ incidence
             // sum_h1_h2(cin_pop, tmpsum_h1_h2_h3);
            sum_h1_h2_h3(cin_pop, tmpsum_h1_h2_h3);

            //      % CIN2 + incidence (move this earlier to re-use cached h1+h2+h3

              //      cin_ft_neg(n, :)   = cin_ft_neg(n - 1, :)   + reshape(sum(sum(sum(sum(sum(sum( cinc_pop(:, : , : , 1,     1, : , : , : ), 1), 2), 3), 4), 6), 8), 1, a);
              //      cin_ft_pos(n, :) =   cin_ft_pos(n - 1, :)   + reshape(sum(sum(sum(sum(sum(sum( cinc_pop(:, : , : , 2 : 3, 1, : , : , : ), 1), 2), 3), 4), 6), 8), 1, a);
              //      cin_ft_art(n, :) =   cin_ft_art(n - 1, :)   + reshape(sum(sum(sum(sum(sum(sum( cinc_pop(:, : , : , 4, 1, : , : , : ), 1), 2), 3), 4), 6), 8), 1, a);

            sum_123_4568_is(spare, tmpsum_h1_h2_h3, 1, 1, 1, 1, sum1_a);
            LOOP_A cin_ft_neg[n_index + _da] = cin_ft_neg[nm1_index + _da] + sum1_a[_da];

            sum_123_4568_is(spare, tmpsum_h1_h2_h3, 2, 3, 1, 1, sum1_a);
            LOOP_A cin_ft_pos[n_index + _da] = cin_ft_pos[nm1_index + _da] + sum1_a[_da];

            sum_123_4568_is(spare, tmpsum_h1_h2_h3, 4, 4, 1, 1, sum1_a);
            LOOP_A cin_ft_art[n_index + _da] = cin_ft_art[nm1_index + _da] + sum1_a[_da];

            sum_i(tmpsum_h1_h2_h3, tmpsum_h1_h2_h3_i);

            //        cin_ft(n, :) = cin_ft(n - 1, :)   + reshape(sum(sum(sum(sum(sum(sum( cinc_pop(:, : , : , : , 1, : , : , : ), 1), 2), 3), 4), 6), 8), 1, a);

            // i dont know if this is correct!?!
            sum_1234_68_s(spare, tmpsum_h1_h2_h3_i, 1, sum1_a);
            LOOP_A cin_ft[n_index + _da] = cin_ft[nm1_index + _da] + sum1_a[_da];


      ///////// END OF NEW /////

      //      % deaths

      //      hivdeaths_pop(:, : , : , 2, : , : , : , : ) = nu(:, : , : , 2, : , : , : , : ).*popn(:, : , : , 2, : , : , : , : ).*dt;
      //      hivdeaths_pop(:, : , : , 3, : , : , : , : ) = nu(:, : , : , 3, : , : , : , : ).*popn(:, : , : , 3, : , : , : , : ).*dt;
      //      hivdeaths_pop(:, : , : , 4, : , : , : , : ) = nu(:, : , : , 4, : , : , : , : ).*popn(:, : , : , 4, : , : , : , : ).*dt;

      LOOP_H1 LOOP_H2 LOOP_H3 LOOP_S LOOP_R LOOP_A LOOP_V1
      for (int _di = 2; _di <= 4; _di++)
        GET_ALL(hivdeaths_pop) = GET_ALL(S->nu) * GET_ALL(popn) * dt;

      sum_h1_h2_h3(hivdeaths_pop, tmpsum_h1_h2_h3);
      sum_i(tmpsum_h1_h2_h3, tmpsum_h1_h2_h3_i);


      //      hivd_ft(n, :) = hivd_ft(n - 1, :) + reshape(sum(sum(sum(sum(sum(sum(hivdeaths_pop(:, : , : , : , 1, : , : , : ), 1), 2), 3), 4), 6), 8), 1, a);
      //      hivd_mt(n, :) = hivd_mt(n - 1, :) + reshape(sum(sum(sum(sum(sum(sum(hivdeaths_pop(:, : , : , : , 2, : , : , : ), 1), 2), 3), 4), 6), 8), 1, a);

      sum_h1_h2_h3_i(hivdeaths_pop, tmpsum_h1_h2_h3_i);
      sum_1234_68_s(spare, tmpsum_h1_h2_h3_i, 1, sum1_a);
      sum_1234_68_s(spare, tmpsum_h1_h2_h3_i, 2, sum2_a);

      LOOP_A {
        hivd_ft[n_index + _da] = hivd_ft[nm1_index + _da] + sum1_a[_da];
        hivd_mt[n_index + _da] = hivd_mt[nm1_index + _da] + sum2_a[_da];
      }

      //      ccdeaths_pop(5, :, : , : , : , : , : , : ) = ups1(5, :, : , : , : , : , : , : ).*popn(5, :, : , : , : , : , : , : ).*dt;
      //      ccdeaths_pop(:, 5, : , : , : , : , : , : ) = ups2(:, 5, : , : , : , : , : , : ).*popn(:, 5, : , : , : , : , : , : ).*dt;
      //      ccdeaths_pop(:, : , 5, : , : , : , : , : ) = ups3(:, : , 5, : , : , : , : , : ).*popn(:, : , 5, : , : , : , : , : ).*dt;

      LOOP_I LOOP_S LOOP_R LOOP_A LOOP_V1 {
        LOOP_H2 LOOP_H3 {
          int _dh1 = 5;
          GET_ALL(ccdeaths_pop) = GET_ALL(S->ups1) * GET_ALL(popn) * dt;
        }
        LOOP_H1 {
          LOOP_H3 {
            int _dh2 = 5;
            GET_ALL(ccdeaths_pop) = GET_ALL(S->ups2) * GET_ALL(popn) * dt;
          }
          LOOP_H2 {
            int _dh3 = 5;
            GET_ALL(ccdeaths_pop) = GET_ALL(S->ups3) * GET_ALL(popn) * dt;
          }
        }
      }

      // See earlier for history update on ccd_ft, ccd_ft_neg, ccd_ft_pos, ccd_ft_art

      //                    ccd_ft(n, :) = ccd_ft(n - 1, :) + reshape(sum(sum(sum(sum(sum(sum(ccdeaths_pop(:, : , : , : , 1, : , : , : ), 1), 2), 3), 4), 6), 8), 1, a);

      sum_h1_h2_h3(ccdeaths_pop, tmpsum_h1_h2_h3);
      sum_i(tmpsum_h1_h2_h3, tmpsum_h1_h2_h3_i);

      sum_1234_68_s(spare, tmpsum_h1_h2_h3_i, 1, sum1_a);
      LOOP_A ccd_ft[n_index + _da] = ccd_ft[nm1_index + _da] + sum1_a[_da];


        //                    ccd_ft_neg(n, :) = ccd_ft_neg(n - 1, :) + reshape(sum(sum(sum(sum(sum(sum(ccdeaths_pop(:, : , : , 1,     1, : , : , : ), 1), 2), 3), 4), 6), 8), 1, a);
        //                    ccd_ft_pos(n, :) = ccd_ft_pos(n - 1, :) + reshape(sum(sum(sum(sum(sum(sum(ccdeaths_pop(:, : , : , 2 : 3, 1, : , : , : ), 1), 2), 3), 4), 6), 8), 1, a);
        //                    ccd_ft_art(n, :) = ccd_ft_art(n - 1, :) + reshape(sum(sum(sum(sum(sum(sum(ccdeaths_pop(:, : , : , 4,     1, : , : , : ), 1), 2), 3), 4), 6), 8), 1, a);

      sum_123_4568_is(spare, tmpsum_h1_h2_h3, 1, 1, 1, 1, sum1_a);
      LOOP_A ccd_ft_neg[n_index + _da] = ccd_ft_neg[nm1_index + _da] + sum1_a[_da];
      sum_123_4568_is(spare, tmpsum_h1_h2_h3, 2, 3, 1, 1, sum1_a);
      LOOP_A ccd_ft_pos[n_index + _da] = ccd_ft_pos[nm1_index + _da] + sum1_a[_da];
      sum_123_4568_is(spare, tmpsum_h1_h2_h3, 4, 4, 1, 1, sum1_a);
      LOOP_A ccd_ft_art[n_index + _da] = ccd_ft_art[nm1_index + _da] + sum1_a[_da];


          //      othdeaths_pop(:, : , : , : , : , : , 1, : ) = mu(:, : , : , : , : , : , 1, : ).*popn(:, : , : , : , : , : , 1, : ).*dt;
          //      othdeaths_pop(:, : , : , : , : , : , 2, : ) = mu(:, : , : , : , : , : , 2, : ).*popn(:, : , : , : , : , : , 2, : ).*dt;
          //      othdeaths_pop(:, : , : , : , : , : , 3, : ) = mu(:, : , : , : , : , : , 3, : ).*popn(:, : , : , : , : , : , 3, : ).*dt;
          //      othdeaths_pop(:, : , : , : , : , : , 4, : ) = mu(:, : , : , : , : , : , 4, : ).*popn(:, : , : , : , : , : , 4, : ).*dt;
          //      othdeaths_pop(:, : , : , : , : , : , 5, : ) = mu(:, : , : , : , : , : , 5, : ).*popn(:, : , : , : , : , : , 5, : ).*dt;

      LOOP_H1 LOOP_H2 LOOP_H3 LOOP_I LOOP_S LOOP_R LOOP_V1
        for (int _da = 1; _da <= 13; _da++)
          GET_ALL(othdeaths_pop) = GET_ALL(S->mu) * GET_ALL(popn) * dt;

      sum_h1_h2_h3_i(othdeaths_pop, tmpsum_h1_h2_h3_i);

      //                      othd_ft(n, :) = othd_ft(n - 1, :) + reshape(sum(sum(sum(sum(sum(sum(othdeaths_pop(:, : , : , : , 1, : , : , : ), 1), 2), 3), 4), 6), 8), 1, a);
      //                      othd_mt(n, :) = othd_mt(n - 1, :) + reshape(sum(sum(sum(sum(sum(sum(othdeaths_pop(:, : , : , : , 2, : , : , : ), 1), 2), 3), 4), 6), 8), 1, a);

      sum_1234_68_s(spare, tmpsum_h1_h2_h3_i, 1, sum1_a);
      LOOP_A othd_ft[n_index + _da] = othd_ft[nm1_index + _da] + sum1_a[_da];
      sum_1234_68_s(spare, tmpsum_h1_h2_h3_i, 2, sum1_a);
      LOOP_A othd_mt[n_index + _da] = othd_mt[nm1_index + _da] + sum1_a[_da];

      for (int _d1 = 1; _d1 <= 3; _d1++)
        for (int _d2 = 1; _d2 <= 2; _d2++)
          screencov_neg[GET_3D_N32_INDEX(n, _d1, _d2)] = screencov_neg[GET_3D_N32_INDEX(n - 1, _d1, _d2)];

     //      screencov(n, 1, 1) = screencov(n - 1, 1, 1) + theta15to24hn(n) * sum(sum(sum(sum(sum(sum(sum(sum(popn(:, : , : , 1, : , : , 2, : ), 1), 2), 3), 4), 5), 6), 7), 8) / 10;% only 25 screened, divide population size 10
     //      screencov(n, 2, 1) = screencov(n - 1, 2, 1) + theta25to34hn(n) * sum(sum(sum(sum(sum(sum(sum(sum(popn(:, : , : , 1, : , : , 3, : ), 1), 2), 3), 4), 5), 6), 7), 8) / 10;
     //      screencov(n, 3, 1) = screencov(n - 1, 3, 1) + theta35to49hn(n) * sum(sum(sum(sum(sum(sum(sum(sum(popn(:, : , : , 1, : , : , 4, : ), 1), 2), 3), 4), 5), 6), 7), 8) / 15;

      screencov_neg[GET_3D_N32_INDEX(n, 1, 1)] += theta25to29hn[n] * sum_all_isra_ch123(spare, popn_sum_h1_h2_h3, 1, 1, 1, 1, 1, dims[DIM_r], 4, 4) * dt;
      screencov_neg[GET_3D_N32_INDEX(n, 3, 1)] += theta35to39hn[n] * sum_all_isra_ch123(spare, popn_sum_h1_h2_h3, 1, 1, 1, 1, 1, dims[DIM_r], 6, 6) * dt ;
      screencov_neg[GET_3D_N32_INDEX(n, 2, 2)] += theta45to49hn[n] * sum_all_isra_ch123(spare, popn_sum_h1_h2_h3, 1, 1, 1, 1, 1, dims[DIM_r], 8, 8) * dt;

     //      screencov(n, 1, 2) = screencov(n - 1, 1, 2) + theta15to24hp(n) * sum(sum(sum(sum(sum(sum(sum(sum(popn(:, : , : , 2 : 4, : , : , 2, : ), 1), 2), 3), 4), 5), 6), 7), 8) / 10;
     //      screencov(n, 2, 2) = screencov(n - 1, 2, 2) + theta25to34hp(n) * sum(sum(sum(sum(sum(sum(sum(sum(popn(:, : , : , 2 : 4, : , : , 3, : ), 1), 2), 3), 4), 5), 6), 7), 8) / 10;
     //      screencov(n, 3, 2) = screencov(n - 1, 3, 2) + theta35to49hp(n) * sum(sum(sum(sum(sum(sum(sum(sum(popn(:, : , : , 2 : 4, : , : , 4, : ), 1), 2), 3), 4), 5), 6), 7), 8) / 15;

      for (int _d1 = 1; _d1 <= 3; _d1++)
        for (int _d2 = 1; _d2 <= 2; _d2++)
	     screencov_pos[GET_3D_N32_INDEX(n, _d1, _d2)] = screencov_pos[GET_3D_N32_INDEX(n - 1, _d1, _d2)];

       screencov_pos[GET_3D_N32_INDEX(n, 1, 1)] += (theta25to29hp[n]) * sum_all_isra_ch123(spare, popn_sum_h1_h2_h3, 2, 3, 1, 1, 1, dims[DIM_r], 4, 4) * dt ;
       screencov_pos[GET_3D_N32_INDEX(n, 2, 1)] += theta30to34hp[n] * sum_all_isra_ch123(spare, popn_sum_h1_h2_h3, 2, 3, 1, 1, 1, dims[DIM_r], 5, 5)  * dt;
       screencov_pos[GET_3D_N32_INDEX(n, 3, 1)] += (theta35to39hp[n]) * sum_all_isra_ch123(spare, popn_sum_h1_h2_h3, 2, 3, 1, 1, 1, dims[DIM_r], 6, 6)  * dt;
       screencov_pos[GET_3D_N32_INDEX(n, 1, 2)] += theta40to44hp[n] * sum_all_isra_ch123(spare, popn_sum_h1_h2_h3, 2, 3, 1, 1, 1, dims[DIM_r], 7, 7)  * dt ;
       screencov_pos[GET_3D_N32_INDEX(n, 2, 2)] += (theta45to49hp[n]) * sum_all_isra_ch123(spare, popn_sum_h1_h2_h3, 2, 3, 1, 1, 1, dims[DIM_r], 8, 8)  * dt ;

       for (int _d1 = 1; _d1 <= 3; _d1++)
         for (int _d2 = 1; _d2 <= 2; _d2++)
 		  screencov_art[GET_3D_N32_INDEX(n, _d1, _d2)] = screencov_art[GET_3D_N32_INDEX(n - 1, _d1, _d2)];

       screencov_art[GET_3D_N32_INDEX(n, 1, 1)] += (theta25to29hp[n])  * sum_all_isra_ch123(spare, popn_sum_h1_h2_h3, 4, 4, 1, 1, 1, dims[DIM_r], 4, 4)  * dt ;
       screencov_art[GET_3D_N32_INDEX(n, 2, 1)] += theta30to34hp[n] * sum_all_isra_ch123(spare, popn_sum_h1_h2_h3, 4, 4, 1, 1, 1, dims[DIM_r], 5, 5) * dt ;
       screencov_art[GET_3D_N32_INDEX(n, 3, 1)] += (theta35to39hp[n]) * sum_all_isra_ch123(spare, popn_sum_h1_h2_h3, 4, 4, 1, 1, 1, dims[DIM_r], 6, 6)  * dt;
       screencov_art[GET_3D_N32_INDEX(n, 1, 2)] += theta40to44hp[n] * sum_all_isra_ch123(spare, popn_sum_h1_h2_h3, 4, 4, 1, 1, 1, dims[DIM_r], 7, 7)  * dt ;
       screencov_art[GET_3D_N32_INDEX(n, 2, 2)] += (theta45to49hp[n]) * sum_all_isra_ch123(spare, popn_sum_h1_h2_h3, 4, 4, 1, 1, 1, dims[DIM_r], 8, 8)  * dt ;


        // end
    }

    //    %% Reset popn
    //    % popn is the "master" population matrix which is update with the changes described in pdt1 - pdt7
    //      popn = popn + pdt1 + pdt2 + pdt3 + pdt4 + pdt5 + pdt6 + pdt7;

    LOOP_LINEAR GET_LIN_ALL(popn) += GET_LIN_ALL(pdt1);
    LOOP_LINEAR GET_LIN_ALL(popn) += GET_LIN_ALL(pdt2);
    LOOP_LINEAR GET_LIN_ALL(popn) += GET_LIN_ALL(pdt3);
    LOOP_LINEAR GET_LIN_ALL(popn) += GET_LIN_ALL(pdt4);
    LOOP_LINEAR GET_LIN_ALL(popn) += GET_LIN_ALL(pdt5);
    LOOP_LINEAR GET_LIN_ALL(popn) += GET_LIN_ALL(pdt6);
    LOOP_LINEAR GET_LIN_ALL(popn) += GET_LIN_ALL(pdt7);

    sum_h1(popn, popn_sum_h1);
    sum_h2(popn_sum_h1, popn_sum_h1_h2);
    sum_h3(popn_sum_h1_h2, popn_sum_h1_h2_h3);
    sum_i(popn_sum_h1_h2_h3, popn_sum_h1_h2_h3_i);


    /*
    if (n == 30000) {
      FILE* fi = fopen("popn_out_mat.30000.bin", "rb");
      LOOP_ALL fread(&GET_ALL(pdt1), 8, 1, fi);
      LOOP_ALL fread(&GET_ALL(pdt2), 8, 1, fi);
      LOOP_ALL fread(&GET_ALL(pdt3), 8, 1, fi);
      LOOP_ALL fread(&GET_ALL(pdt4), 8, 1, fi);
      LOOP_ALL fread(&GET_ALL(pdt5), 8, 1, fi);
      LOOP_ALL fread(&GET_ALL(pdt6), 8, 1, fi);
      LOOP_ALL fread(&GET_ALL(pdt7), 8, 1, fi);
      LOOP_ALL fread(&GET_ALL(popn), 8, 1, fi);
      fclose(fi);
    }
    */
   /*
    if (n % 1000 == 0) {
    char ff[50];
    sprintf(ff, "popn_out_c.%d.bin", n);
    FILE* fo = fopen(ff, "wb");
    LOOP_ALL fwrite(&GET_ALL(pdt1), 8, 1, fo);
    LOOP_ALL fwrite(&GET_ALL(pdt2), 8, 1, fo);
    LOOP_ALL fwrite(&GET_ALL(pdt3), 8, 1, fo);
    LOOP_ALL fwrite(&GET_ALL(pdt4), 8, 1, fo);
    LOOP_ALL fwrite(&GET_ALL(pdt5), 8, 1, fo);
    LOOP_ALL fwrite(&GET_ALL(pdt6), 8, 1, fo);
    LOOP_ALL fwrite(&GET_ALL(pdt7), 8, 1, fo);
    LOOP_ALL fwrite(&GET_ALL(popn), 8, 1, fo);
    fclose(fo);

    }
   */



    //      %% Outputs
    //      % every year from HIV infection start
    //      htyrid = hyrid(gg);
    //      if n == htyrid

    int htyrid = (int)round(hyrid[gg]);
    if (n == htyrid) {

      //      hivsus_pop(gg, :, : , : , : ) = reshape(sum(sum(sum(sum(popn(:, : , : , 1, : , : , : , : ), 1), 2), 3), 6), 1, s, a, v1);

      sum_hivsus_pop(spare, popn_sum_h1_h2_h3, R->hivsus_pop[gg][1]);

      // hiv_15to49(gg, :) = reshape(sum(sum(sum(sum(sum(sum(sum(popn(:,:,:,2:4,:,:,2:8,:),1),2),3),4),6),7),8)./sum(sum(sum(sum(sum(sum(sum(popn(:,:,:,1:4,:,:,2:8,:),1),2),3),4),6),7),8),1,s);

      sum_123_4678_ia(spare, popn_sum_h1_h2_h3, 2, 4, 2, 8, sum1_s);
      sum_123_4678_ia(spare, popn_sum_h1_h2_h3, 1, 4, 2, 8, sum2_s);
      LOOP_S R->hiv_15to49[gg][_ds] = sum1_s[_ds] / sum2_s[_ds];

      // hiv_15to74(gg, :) = reshape(sum(sum(sum(sum(sum(sum(sum(popn(:,:,:,2:4,:,:,2:end,:),1),2),3),4),6),7),8) ./ sum(sum(sum(sum(sum(sum(sum(popn(:,:,:,1:4,:,:,2:end,:),1),2),3),4),6),7),8),1,s);

      sum_123_4678_ia(spare, popn_sum_h1_h2_h3, 2, 4, 2, _SIZE_A, sum1_s);
      sum_123_4678_ia(spare, popn_sum_h1_h2_h3, 1, 4, 2, _SIZE_A, sum2_s);
      LOOP_S R->hiv_15to74[gg][_ds] = sum1_s[_ds] / sum2_s[_ds];

      //       hiv_9to14(gg, :) = reshape(sum(sum(sum(sum(sum(sum(sum(popn(:, : , : , 2 : 4, : , : , 1, : ), 1), 2), 3), 4), 6), 7), 8) . /
      //                                  sum(sum(sum(sum(sum(sum(sum(popn(:, : , : , 1 : 4, : , : , 1, : ), 1), 2), 3), 4), 6), 7), 8), 1, s);

      sum_123_4678_ia(spare, popn_sum_h1_h2_h3, 2, 4, 1, 1, sum1_s);
      sum_123_4678_ia(spare, popn_sum_h1_h2_h3, 1, 4, 1, 1, sum2_s);
      LOOP_S R->hiv_9to14[gg][_ds] = sum1_s[_ds] / sum2_s[_ds];

      //      hivart_15to49(gg, :) = reshape(sum(sum(sum(sum(sum(sum(sum(popn(:, : , : , 4,     : , : , 2 : 8, : ), 1), 2), 3), 4), 6), 7), 8) . /
      //                                     sum(sum(sum(sum(sum(sum(sum(popn(:, : , : , 2 : 4, : , : , 2 : 8, : ), 1), 2), 3), 4), 6), 7), 8), 1, s);

      sum_123_4678_ia(spare, popn_sum_h1_h2_h3, 4, 4, 2, 8, sum1_s);
      sum_123_4678_ia(spare, popn_sum_h1_h2_h3, 2, 4, 2, 8, sum2_s);
      LOOP_S R->hivart_15to49[gg][_ds] = sum1_s[_ds] / sum2_s[_ds];

      //      hivart_15to74(gg, :) = reshape(sum(sum(sum(sum(sum(sum(sum(popn(:, : , : , 4,     : , : , 2 : end, : ), 1), 2), 3), 4), 6), 7), 8) . /
      //                                     sum(sum(sum(sum(sum(sum(sum(popn(:, : , : , 2 : 4, : , : , 2 : end, : ), 1), 2), 3), 4), 6), 7), 8), 1, s);

      sum_123_4678_ia(spare, popn_sum_h1_h2_h3, 4, 4, 2, _SIZE_A, sum1_s);
      sum_123_4678_ia(spare, popn_sum_h1_h2_h3, 2, 4, 2, _SIZE_A, sum2_s);
      LOOP_S R->hivart_15to74[gg][_ds] = sum1_s[_ds] / sum2_s[_ds];

      //      hivart_9to14(gg, :) = reshape(sum(sum(sum(sum(sum(sum(sum(popn(:, : , : , 4,     : , : , 1, : ), 1), 2), 3), 4), 6), 7), 8) . /
      //                                    sum(sum(sum(sum(sum(sum(sum(popn(:, : , : , 2 : 4, : , : , 1, : ), 1), 2), 3), 4), 6), 7), 8), 1, s);

      sum_123_4678_ia(spare, popn_sum_h1_h2_h3, 4, 4, 1, 1, sum1_s);
      sum_123_4678_ia(spare, popn_sum_h1_h2_h3, 2, 4, 1, 1, sum2_s);
      LOOP_S R->hivart_9to14[gg][_ds] = sum1_s[_ds] / sum2_s[_ds];

      //      art_15plus(gg, :) = reshape(sum(sum(sum(sum(sum(sum(sum(popn(:, : , : , 4,     : , : , 2 : end, : ), 1), 2), 3), 4), 6), 7), 8) . /
      //                                  sum(sum(sum(sum(sum(sum(sum(popn(:, : , : , 2 : 4, : , : , 2 : end, : ), 1), 2), 3), 4), 6), 7), 8), 1, s);

      sum_123_4678_ia(spare, popn_sum_h1_h2_h3, 4, 4, 2, _SIZE_A, sum1_s);
      sum_123_4678_ia(spare, popn_sum_h1_h2_h3, 2, 4, 2, _SIZE_A, sum2_s);
      LOOP_S R->art_15plus[gg][_ds] = sum1_s[_ds] / sum2_s[_ds];

      //art_age(gg,1)=reshape(sum(sum(sum(sum(sum(sum(sum(sum(popn(:,:,:,4,:,:,1,   :),1),2),3),4),5),6),7),8)./sum(sum(sum(sum(sum(sum(sum(sum(popn(:,:,:,2:4,:,:,1,   :),1),2),3),4),5),6),7),8),1,1);
      //art_age(gg,2)=reshape(sum(sum(sum(sum(sum(sum(sum(sum(popn(:,:,:,4,:,:,2:3, :),1),2),3),4),5),6),7),8)./sum(sum(sum(sum(sum(sum(sum(sum(popn(:,:,:,2:4,:,:,2:3, :),1),2),3),4),5),6),7),8),1,1);
      //art_age(gg,3)=reshape(sum(sum(sum(sum(sum(sum(sum(sum(popn(:,:,:,4,:,:,4:5, :),1),2),3),4),5),6),7),8)./sum(sum(sum(sum(sum(sum(sum(sum(popn(:,:,:,2:4,:,:,4:5, :),1),2),3),4),5),6),7),8),1,1);
      //art_age(gg,4)=reshape(sum(sum(sum(sum(sum(sum(sum(sum(popn(:,:,:,4,:,:,6:8, :),1),2),3),4),5),6),7),8)./sum(sum(sum(sum(sum(sum(sum(sum(popn(:,:,:,2:4,:,:,6:8, :),1),2),3),4),5),6),7),8),1,1);
      //art_age(gg,5)=reshape(sum(sum(sum(sum(sum(sum(sum(sum(popn(:,:,:,4,:,:,9:13,:),1),2),3),4),5),6),7),8)./sum(sum(sum(sum(sum(sum(sum(sum(popn(:,:,:,2:4,:,:,9:13,:),1),2),3),4),5),6),7),8),1,1);

      R->art_age[gg][1] = sum_all_isra_ch123(spare, popn_sum_h1_h2_h3, 4, 4, 1, _SIZE_S, 1, _SIZE_R, 1, 1) / sum_all_isra_ch123(spare, popn_sum_h1_h2_h3, 2, 4, 1, _SIZE_S, 1, _SIZE_R, 1, 1);
      R->art_age[gg][2] = sum_all_isra_ch123(spare, popn_sum_h1_h2_h3, 4, 4, 1, _SIZE_S, 1, _SIZE_R, 2, 3) / sum_all_isra_ch123(spare, popn_sum_h1_h2_h3, 2, 4, 1, _SIZE_S, 1, _SIZE_R, 2, 3);
      R->art_age[gg][3] = sum_all_isra_ch123(spare, popn_sum_h1_h2_h3, 4, 4, 1, _SIZE_S, 1, _SIZE_R, 4, 5) / sum_all_isra_ch123(spare, popn_sum_h1_h2_h3, 2, 4, 1, _SIZE_S, 1, _SIZE_R, 4, 5);
      R->art_age[gg][4] = sum_all_isra_ch123(spare, popn_sum_h1_h2_h3, 4, 4, 1, _SIZE_S, 1, _SIZE_R, 6, 8) / sum_all_isra_ch123(spare, popn_sum_h1_h2_h3, 2, 4, 1, _SIZE_S, 1, _SIZE_R, 6, 8);
      R->art_age[gg][5] = sum_all_isra_ch123(spare, popn_sum_h1_h2_h3, 4, 4, 1, _SIZE_S, 1, _SIZE_R, 9, 13) / sum_all_isra_ch123(spare, popn_sum_h1_h2_h3, 2, 4, 1, _SIZE_S, 1, _SIZE_R, 9, 13);

      // hiv_fage(gg, 1) = reshape(sum(sum(sum(sum(sum(sum(sum(sum(popn(:, : , : , 2 : 4, 1, : , 1, : ), 1), 2), 3), 4), 5), 6), 7), 8) . /     sum(sum(sum(sum(sum(sum(sum(sum(popn(:, : , : , 1 : 4, 1, : , 1, : ), 1), 2), 3), 4), 5), 6), 7), 8), 1, 1);
      // hiv_fage(gg, 2) = reshape(sum(sum(sum(sum(sum(sum(sum(sum(popn(:, : , : , 2 : 4, 1, : , 2 : 3, : ), 1), 2), 3), 4), 5), 6), 7), 8).  / sum(sum(sum(sum(sum(sum(sum(sum(popn(:, : , : , 1 : 4, 1, : , 2 : 3, : ), 1), 2), 3), 4), 5), 6), 7), 8), 1, 1);
      // hiv_fage(gg, 3) = reshape(sum(sum(sum(sum(sum(sum(sum(sum(popn(:, : , : , 2 : 4, 1, : , 4 : 5, : ), 1), 2), 3), 4), 5), 6), 7), 8).  / sum(sum(sum(sum(sum(sum(sum(sum(popn(:, : , : , 1 : 4, 1, : , 4 : 5, : ), 1), 2), 3), 4), 5), 6), 7), 8), 1, 1);
      // hiv_fage(gg, 4) = reshape(sum(sum(sum(sum(sum(sum(sum(sum(popn(:, : , : , 2 : 4, 1, : , 6 : 8, : ), 1), 2), 3), 4), 5), 6), 7), 8).  / sum(sum(sum(sum(sum(sum(sum(sum(popn(:, : , : , 1 : 4, 1, : , 6 : 8, : ), 1), 2), 3), 4), 5), 6), 7), 8), 1, 1);
      // hiv_fage(gg, 5) = reshape(sum(sum(sum(sum(sum(sum(sum(sum(popn(:, : , : , 2 : 4, 1, : , 9 : 13, : ), 1), 2), 3), 4), 5), 6), 7), 8). / sum(sum(sum(sum(sum(sum(sum(sum(popn(:, : , : , 1 : 4, 1, : , 9 : 13, : ), 1), 2), 3), 4), 5), 6), 7), 8), 1, 1);

      R->hiv_fage[gg][1] = sum_all_isra_ch123(spare, popn_sum_h1_h2_h3, 2, 4, 1, 1, 1, _SIZE_R, 1, 1) / sum_all_isra_ch123(spare, popn_sum_h1_h2_h3, 1, 4, 1, 1, 1, _SIZE_R, 1, 1);
      R->hiv_fage[gg][2] = sum_all_isra_ch123(spare, popn_sum_h1_h2_h3, 2, 4, 1, 1, 1, _SIZE_R, 2, 3) / sum_all_isra_ch123(spare, popn_sum_h1_h2_h3, 1, 4, 1, 1, 1, _SIZE_R, 2, 3);
      R->hiv_fage[gg][3] = sum_all_isra_ch123(spare, popn_sum_h1_h2_h3, 2, 4, 1, 1, 1, _SIZE_R, 4, 5) / sum_all_isra_ch123(spare, popn_sum_h1_h2_h3, 1, 4, 1, 1, 1, _SIZE_R, 4, 5);
      R->hiv_fage[gg][4] = sum_all_isra_ch123(spare, popn_sum_h1_h2_h3, 2, 4, 1, 1, 1, _SIZE_R, 6, 8) / sum_all_isra_ch123(spare, popn_sum_h1_h2_h3, 1, 4, 1, 1, 1, _SIZE_R, 6, 8);
      R->hiv_fage[gg][5] = sum_all_isra_ch123(spare, popn_sum_h1_h2_h3, 2, 4, 1, 1, 1, _SIZE_R, 9, 13) / sum_all_isra_ch123(spare, popn_sum_h1_h2_h3, 1, 4, 1, 1, 1, _SIZE_R, 9, 13);

      // hiv_mage(gg, 1) = reshape(sum(sum(sum(sum(sum(sum(sum(sum(popn(:, : , : , 2 : 4, 2, : , 1, : ), 1), 2), 3), 4), 5), 6), 7), 8) . / sum(sum(sum(sum(sum(sum(sum(sum(popn(:, : , : , 1 : 4, 2, : , 1, : ), 1), 2), 3), 4), 5), 6), 7), 8), 1, 1);
      // hiv_mage(gg, 2) = reshape(sum(sum(sum(sum(sum(sum(sum(sum(popn(:, : , : , 2 : 4, 2, : , 2 : 3, : ), 1), 2), 3), 4), 5), 6), 7), 8) . / sum(sum(sum(sum(sum(sum(sum(sum(popn(:, : , : , 1 : 4, 2, : , 2 : 3, : ), 1), 2), 3), 4), 5), 6), 7), 8), 1, 1);
      // hiv_mage(gg, 3) = reshape(sum(sum(sum(sum(sum(sum(sum(sum(popn(:, : , : , 2 : 4, 2, : , 4 : 5, : ), 1), 2), 3), 4), 5), 6), 7), 8) . / sum(sum(sum(sum(sum(sum(sum(sum(popn(:, : , : , 1 : 4, 2, : , 4 : 5, : ), 1), 2), 3), 4), 5), 6), 7), 8), 1, 1);
      // hiv_mage(gg, 4) = reshape(sum(sum(sum(sum(sum(sum(sum(sum(popn(:, : , : , 2 : 4, 2, : , 6 : 8, : ), 1), 2), 3), 4), 5), 6), 7), 8) . / sum(sum(sum(sum(sum(sum(sum(sum(popn(:, : , : , 1 : 4, 2, : , 6 : 8, : ), 1), 2), 3), 4), 5), 6), 7), 8), 1, 1);
      // hiv_mage(gg, 5) = reshape(sum(sum(sum(sum(sum(sum(sum(sum(popn(:, : , : , 2 : 4, 2, : , 9 : 13, : ), 1), 2), 3), 4), 5), 6), 7), 8) . / sum(sum(sum(sum(sum(sum(sum(sum(popn(:, : , : , 1 : 4, 2, : , 9 : 13, : ), 1), 2), 3), 4), 5), 6), 7), 8), 1, 1);

      R->hiv_mage[gg][1] = sum_all_isra_ch123(spare, popn_sum_h1_h2_h3, 2, 4, 2, 2, 1, _SIZE_R, 1, 1) / sum_all_isra_ch123(spare, popn_sum_h1_h2_h3, 1, 4, 2, 2, 1, _SIZE_R, 1, 1);
      R->hiv_mage[gg][2] = sum_all_isra_ch123(spare, popn_sum_h1_h2_h3, 2, 4, 2, 2, 1, _SIZE_R, 2, 3) / sum_all_isra_ch123(spare, popn_sum_h1_h2_h3, 1, 4, 2, 2, 1, _SIZE_R, 2, 3);
      R->hiv_mage[gg][3] = sum_all_isra_ch123(spare, popn_sum_h1_h2_h3, 2, 4, 2, 2, 1, _SIZE_R, 4, 5) / sum_all_isra_ch123(spare, popn_sum_h1_h2_h3, 1, 4, 2, 2, 1, _SIZE_R, 4, 5);
      R->hiv_mage[gg][4] = sum_all_isra_ch123(spare, popn_sum_h1_h2_h3, 2, 4, 2, 2, 1, _SIZE_R, 6, 8) / sum_all_isra_ch123(spare, popn_sum_h1_h2_h3, 1, 4, 2, 2, 1, _SIZE_R, 6, 8);
      R->hiv_mage[gg][5] = sum_all_isra_ch123(spare, popn_sum_h1_h2_h3, 2, 4, 2, 2, 1, _SIZE_R, 9, 13) / sum_all_isra_ch123(spare, popn_sum_h1_h2_h3, 1, 4, 2, 2, 1, _SIZE_R, 9, 13);

      // hiv_fsw(gg, :) = sum(sum(sum(sum(sum(sum(sum(sum(popn(:,:,:,3:4,1,3,2:11,:),1),2),3),4),5),6),7),8)./ sum(sum(sum(sum(sum(sum(sum(sum(popn(:,:,:,1:4,1,3,2:11,:),1),2),3),4),5),6),7),8); % 15 - 64 yos

      R->hiv_fsw[gg] = sum_all_isra_ch123(spare, popn_sum_h1_h2_h3, 2, 4, 1, 1, 3, 3, 2, 11) /
                       sum_all_isra_ch123(spare, popn_sum_h1_h2_h3, 1, 4, 1, 1, 3, 3, 2, 11);

      // hiv_mcli(gg, :) = sum(sum(sum(sum(sum(sum(sum(sum(popn(:, : , : , 2 : 4, 2, 3, 2 : 5, : ), 1), 2), 3), 4), 5), 6), 7), 8) . / sum(sum(sum(sum(sum(sum(sum(sum(popn(:, : , : , 1 : 4, 2, 3, 2 : 5, : ), 1), 2), 3), 4), 5), 6), 7), 8);
      // hiv_mcli(gg, :) = sum(sum(sum(sum(sum(sum(sum(sum(popn(:, : , : , 2 : 4, 2, 3, 2 : 11, : ), 1), 2), 3), 4), 5), 6), 7), 8) . / sum(sum(sum(sum(sum(sum(sum(sum(popn(:, : , : , 1 : 4, 2, 3, 2 : 11, : ), 1), 2), 3), 4), 5), 6), 7), 8); % 15 - 64 yos

      R->hiv_mcli[gg] = sum_all_isra_ch123(spare, popn_sum_h1_h2_h3, 2, 4, 2, 2, 3, 3, 2, 11) /
                        sum_all_isra_ch123(spare, popn_sum_h1_h2_h3, 1, 4, 2, 2, 3, 3, 2, 11);

      //      vacc_age_f(gg, :) = reshape(sum(sum(sum(sum(sum(sum(sum(popn(:, : , : , : , 1, : , : , 2), 1), 2), 3), 4), 5), 6), 8), 1, a);
      //      vacc_age_m(gg, :) = reshape(sum(sum(sum(sum(sum(sum(sum(popn(:, : , : , : , 2, : , : , 2), 1), 2), 3), 4), 5), 6), 8), 1, a);
      sum_1234_568_sv1(spare, popn_sum_h1_h2_h3_i, 1, 2, R->vacc_age_f[gg]);
      sum_1234_568_sv1(spare, popn_sum_h1_h2_h3_i, 2, 2, R->vacc_age_m[gg]);

        //      vacc_age_f_pos(gg, :) = reshape(sum(sum(sum(sum(sum(sum(sum(popn(:, : , : , 2 : 4, 1, : , : , 2), 1), 2), 3), 4), 5), 6), 8), 1, a);
        //      vacc_age_m_pos(gg, :) = reshape(sum(sum(sum(sum(sum(sum(sum(popn(:, : , : , 2 : 4, 2, : , : , 2), 1), 2), 3), 4), 5), 6), 8), 1, a);

      sum_123_4568_isv1(spare, popn_sum_h1_h2_h3, 2, 3, 1, 1, 2, 2, R->vacc_age_f_pos[gg]);
      sum_123_4568_isv1(spare, popn_sum_h1_h2_h3, 4, 4, 1, 1, 2, 2, R->vacc_age_f_art[gg]);
      sum_123_4568_isv1(spare, popn_sum_h1_h2_h3, 2, 4, 2, 2, 2, 2, R->vacc_age_m_pos[gg]);

        //      testtime(gg, :) = n;

      R->testtime[gg] = n;

      //      gg = gg + 1;
      gg++;

      //      % if gg == 5
      //      % keyboard
      //      % end

      //      end
    }

    //    tyrid = yrid(g);
    //    if n == tyrid

    int tyrid = (int)round(yrid[g]);
    if (n == tyrid) {

      //      psize(g, :) = reshape(sum(sum(sum(sum(sum(sum(sum(popn(:, : , : , : , : , : , : , : ), 1), 2), 3), 4), 6), 7), 8), 1, s);
      sum_all_1234_s(spare, popn_sum_h1_h2_h3_i, R->psize[g]);

      //    psize_age_f(g, :) =     reshape(sum(sum(sum(sum(sum(sum(sum(popn(:, : , : , :    , 1, : , : , : ), 1), 2), 3), 4), 5), 6), 8), 1, a);
      //    psize_age_m(g, :) =     reshape(sum(sum(sum(sum(sum(sum(sum(popn(:, : , : , : ,    2, : , : , : ), 1), 2), 3), 4), 5), 6), 8), 1, a);

      //    psize_age_f_neg(g, :) = reshape(sum(sum(sum(sum(sum(sum(sum(popn(:, : , : , 1,     1, : , : , : ), 1), 2), 3), 4), 5), 6), 8), 1, a);
      //    psize_age_f_pos(g, :) = reshape(sum(sum(sum(sum(sum(sum(sum(popn(:, : , : , 2 : 3, 1, : , : , : ), 1), 2), 3), 4), 5), 6), 8), 1, a);
      //    psize_age_f_art(g, :) = reshape(sum(sum(sum(sum(sum(sum(sum(popn(:, : , : , 4,     1, : , : , : ), 1), 2), 3), 4), 5), 6), 8), 1, a);
      //    psize_age_m_pos(g, :) = reshape(sum(sum(sum(sum(sum(sum(sum(popn(:, : , : , 2 : 3, 2, : , : , : ), 1), 2), 3), 4), 5), 6), 8), 1, a);
      //    psize_age_m_art(g, :) = reshape(sum(sum(sum(sum(sum(sum(sum(popn(:, : , : , 4,     2, : , : , : ), 1), 2), 3), 4), 5), 6), 8), 1, a);

      sum_1234_568_s(spare, popn_sum_h1_h2_h3_i, 1, R->psize_age_f[g]);
      sum_1234_568_s(spare, popn_sum_h1_h2_h3_i, 2, R->psize_age_m[g]);

      sum_123_4568_is(spare, popn_sum_h1_h2_h3, 1, 1, 1, 1, R->psize_age_f_neg[g]);
      sum_123_4568_is(spare, popn_sum_h1_h2_h3, 2, 3, 1, 1, R->psize_age_f_pos[g]);
      sum_123_4568_is(spare, popn_sum_h1_h2_h3, 4, 4, 1, 1, R->psize_age_f_art[g]);
      sum_123_4568_is(spare, popn_sum_h1_h2_h3, 1, 1, 2, 2, R->psize_age_m_neg[g]);
      sum_123_4568_is(spare, popn_sum_h1_h2_h3, 2, 3, 2, 2, R->psize_age_m_pos[g]);
      sum_123_4568_is(spare, popn_sum_h1_h2_h3, 4, 4, 2, 2, R->psize_age_m_art[g]);

        //      hpv1618_f(g, :, : ) = reshape(sum(sum(sum(sum(sum(popn([2 4 5], :, : , : , 1, : , : , : ), 1), 2), 3), 6), 8), 1, i, a);
        //      hpv1618_m(g, :, : ) = reshape(sum(sum(sum(sum(sum(popn([2 4 5], :, : , : , 2, : , : , : ), 1), 2), 3), 6), 8), 1, i, a);

      LOOP_H2 LOOP_H3 LOOP_I LOOP_R LOOP_A LOOP_V1 {
        R->hpv1618_f[g][_di][_da] += GET(popn, 2, _dh2, _dh3, _di, 1, _dr, _da, _dv1) +
                                     GET(popn, 4, _dh2, _dh3, _di, 1, _dr, _da, _dv1) +
                                     GET(popn, 5, _dh2, _dh3, _di, 1, _dr, _da, _dv1);
        R->hpv1618_m[g][_di][_da] += GET(popn, 2, _dh2, _dh3, _di, 2, _dr, _da, _dv1) +
                                     GET(popn, 4, _dh2, _dh3, _di, 2, _dr, _da, _dv1) +
                                     GET(popn, 5, _dh2, _dh3, _di, 2, _dr, _da, _dv1);
      }

        //      nvthpv_f(g, :, : ) = reshape(sum(sum(sum(sum(sum(popn(:, : , [2 4 5], : , 1, : , : , : ), 1), 2), 3), 6), 8), 1, i, a);
        //      nvthpv_m(g, :, : ) = reshape(sum(sum(sum(sum(sum(popn(:, : , [2 4 5], : , 2, : , : , : ), 1), 2), 3), 6), 8), 1, i, a);

      LOOP_I LOOP_R LOOP_A LOOP_V1 {
        R->nvthpv_f[g][_di][_da] += GET(popn_sum_h1_h2, 1, 1, 2, _di, 1, _dr, _da, _dv1)
                                  + GET(popn_sum_h1_h2, 1, 1, 4, _di, 1, _dr, _da, _dv1)
                                  + GET(popn_sum_h1_h2, 1, 1, 5, _di, 1, _dr, _da, _dv1);
        R->nvthpv_m[g][_di][_da] += GET(popn_sum_h1_h2, 1, 1, 2, _di, 2, _dr, _da, _dv1)
                                  + GET(popn_sum_h1_h2, 1, 1, 4, _di, 2, _dr, _da, _dv1)
                                  + GET(popn_sum_h1_h2, 1, 1, 5, _di, 2, _dr, _da, _dv1);
      }

        //      hpv9vt_f(g, 1, 1) = sum(popn(hpv12_id == 1 & sexf_id == 1 & hivS_id == 1 & age_id == 1));
        //      hpv9vt_f(g, 1, 2) = sum(popn(hpv12_id == 1 & sexf_id == 1 & hivS_id == 1 & age_id == 2));
        //      hpv9vt_f(g, 1, 3) = sum(popn(hpv12_id == 1 & sexf_id == 1 & hivS_id == 1 & age_id == 3));
        //      hpv9vt_f(g, 1, 4) = sum(popn(hpv12_id == 1 & sexf_id == 1 & hivS_id == 1 & age_id == 4));
        //      hpv9vt_f(g, 1, 5) = sum(popn(hpv12_id == 1 & sexf_id == 1 & hivS_id == 1 & age_id == 5));

        //      hpv9vt_f(g, 2, 1) = sum(popn(hpv12_id == 1 & sexf_id == 1 & noart_id == 1 & age_id == 1));
        //      hpv9vt_f(g, 2, 2) = sum(popn(hpv12_id == 1 & sexf_id == 1 & noart_id == 1 & age_id == 2));
        //      hpv9vt_f(g, 2, 3) = sum(popn(hpv12_id == 1 & sexf_id == 1 & noart_id == 1 & age_id == 3));
        //      hpv9vt_f(g, 2, 4) = sum(popn(hpv12_id == 1 & sexf_id == 1 & noart_id == 1 & age_id == 4));
        //      hpv9vt_f(g, 2, 5) = sum(popn(hpv12_id == 1 & sexf_id == 1 & noart_id == 1 & age_id == 5));

        //      hpv9vt_f(g, 3, 1) = sum(popn(hpv12_id == 1 & sexf_id == 1 & art_id == 1 & age_id == 1));
        //      hpv9vt_f(g, 3, 2) = sum(popn(hpv12_id == 1 & sexf_id == 1 & art_id == 1 & age_id == 2));
        //      hpv9vt_f(g, 3, 3) = sum(popn(hpv12_id == 1 & sexf_id == 1 & art_id == 1 & age_id == 3));
        //      hpv9vt_f(g, 3, 4) = sum(popn(hpv12_id == 1 & sexf_id == 1 & art_id == 1 & age_id == 4));
        //      hpv9vt_f(g, 3, 5) = sum(popn(hpv12_id == 1 & sexf_id == 1 & art_id == 1 & age_id == 5));

        //      hpv9vt_m(g, 1, 1) = sum(popn(hpv12_id == 1 & sexm_id == 1 & hivS_id == 1 & age_id == 1));
        //      hpv9vt_m(g, 1, 2) = sum(popn(hpv12_id == 1 & sexm_id == 1 & hivS_id == 1 & age_id == 2));
        //      hpv9vt_m(g, 1, 3) = sum(popn(hpv12_id == 1 & sexm_id == 1 & hivS_id == 1 & age_id == 3));
        //      hpv9vt_m(g, 1, 4) = sum(popn(hpv12_id == 1 & sexm_id == 1 & hivS_id == 1 & age_id == 4));
        //      hpv9vt_m(g, 1, 5) = sum(popn(hpv12_id == 1 & sexm_id == 1 & hivS_id == 1 & age_id == 5));

        //      hpv9vt_m(g, 2, 1) = sum(popn(hpv12_id == 1 & sexm_id == 1 & noart_id == 1 & age_id == 1));
        //      hpv9vt_m(g, 2, 2) = sum(popn(hpv12_id == 1 & sexm_id == 1 & noart_id == 1 & age_id == 2));
        //      hpv9vt_m(g, 2, 3) = sum(popn(hpv12_id == 1 & sexm_id == 1 & noart_id == 1 & age_id == 3));
        //      hpv9vt_m(g, 2, 4) = sum(popn(hpv12_id == 1 & sexm_id == 1 & noart_id == 1 & age_id == 4));
        //      hpv9vt_m(g, 2, 5) = sum(popn(hpv12_id == 1 & sexm_id == 1 & noart_id == 1 & age_id == 5));

        //      hpv9vt_m(g, 3, 1) = sum(popn(hpv12_id == 1 & sexm_id == 1 & art_id == 1 & age_id == 1));
        //      hpv9vt_m(g, 3, 2) = sum(popn(hpv12_id == 1 & sexm_id == 1 & art_id == 1 & age_id == 2));
        //      hpv9vt_m(g, 3, 3) = sum(popn(hpv12_id == 1 & sexm_id == 1 & art_id == 1 & age_id == 3));
        //      hpv9vt_m(g, 3, 4) = sum(popn(hpv12_id == 1 & sexm_id == 1 & art_id == 1 & age_id == 4));
        //      hpv9vt_m(g, 3, 5) = sum(popn(hpv12_id == 1 & sexm_id == 1 & art_id == 1 & age_id == 5));

        //      hpv_f(g, 1, 1) = sum(popn(hpv123_id == 1 & sexf_id == 1 & hivS_id == 1 & age_id == 1));
        //      hpv_f(g, 1, 2) = sum((hpv123_id == 1 & sexf_id == 1 & hivS_id == 1 & age_id == 2));
        //      hpv_f(g, 1, 3) = sum(popn(hpv123_id == 1 & sexf_id == 1 & hivS_id == 1 & age_id == 3));
        //      hpv_f(g, 1, 4) = sum(popn(hpv123_id == 1 & sexf_id == 1 & hivS_id == 1 & age_id == 4));
        //      hpv_f(g, 1, 5) = sum(popn(hpv123_id == 1 & sexf_id == 1 & hivS_id == 1 & age_id == 5));

        //      hpv_f(g, 2, 1) = sum(popn(hpv123_id == 1 & sexf_id == 1 & noart_id == 1 & age_id == 1));
        //      hpv_f(g, 2, 2) = sum(popn(hpv123_id == 1 & sexf_id == 1 & noart_id == 1 & age_id == 2));
        //      hpv_f(g, 2, 3) = sum(popn(hpv123_id == 1 & sexf_id == 1 & noart_id == 1 & age_id == 3));
        //      hpv_f(g, 2, 4) = sum(popn(hpv123_id == 1 & sexf_id == 1 & noart_id == 1 & age_id == 4));
        //      hpv_f(g, 2, 5) = sum(popn(hpv123_id == 1 & sexf_id == 1 & noart_id == 1 & age_id == 5));

        //      hpv_f(g, 3, 1) = sum(popn(hpv123_id == 1 & sexf_id == 1 & art_id == 1 & age_id == 1));
        //      hpv_f(g, 3, 2) = sum(popn(hpv123_id == 1 & sexf_id == 1 & art_id == 1 & age_id == 2));
        //      hpv_f(g, 3, 3) = sum(popn(hpv123_id == 1 & sexf_id == 1 & art_id == 1 & age_id == 3));
        //      hpv_f(g, 3, 4) = sum(popn(hpv123_id == 1 & sexf_id == 1 & art_id == 1 & age_id == 4));
        //      hpv_f(g, 3, 5) = sum(popn(hpv123_id == 1 & sexf_id    1 & art_id == 1 & age_id == 5));


        //      hpv_m(g, 1, 1) = sum(popn(hpv123_id == 1 & sexm_id == 1 & hivS_id == 1 & age_id == 1));
        //      hpv_m(g, 1, 2) = sum(popn(hpv123_id == 1 & sexm_id =grep = 1 & hivS_id == 1 & age_id == 2));
        //      hpv_m(g, 1, 3) = sum(popn(hpv123_id == 1 & sexm_id == 1 & hivS_id == 1 & age_id == 3));
        //      hpv_m(g, 1, 4) = sum(popn(hpv123_id == 1 & sexm_id == 1 & hivS_id == 1 & age_id == 4));
        //      hpv_m(g, 1, 5) = sum(popn(hpv123_id == 1 & sexm_id == 1 & hivS_id == 1 & age_id == 5));

        //      hpv_m(g, 2, 1) = sum(popn(hpv123_id == 1 & sexm_id == 1 & noart_id == 1 & age_id == 1));
        //      hpv_m(g, 2, 2) = sum(popn(hpv123_id == 1 & sexm_id == 1 & noart_id == 1 & age_id == 2));
        //      hpv_m(g, 2, 3) = sum(popn(hpv123_id == 1 & sexm_id == 1 & noart_id == 1 & age_id == 3));
        //      hpv_m(g, 2, 4) = sum(popn(hpv123_id == 1 & sexm_id == 1 & noart_id == 1 & age_id == 4));
        //      hpv_m(g, 2, 5) = sum(popn(hpv123_id == 1 & sexm_id == 1 & noart_id == 1 & age_id == 5));

        //      hpv_m(g, 3, 1) = sum(popn(hpv123_id == 1 & sexm_id == 1 & art_id == 1 & age_id == 1));
        //      hpv_m(g, 3, 2) = sum(popn(hpv123_id == 1 & sexm_id == 1 & art_id == 1 & age_id == 2));
        //      hpv_m(g, 3, 3) = sum(popn(hpv123_id == 1 & sexm_id == 1 & art_id == 1 & age_id == 3));
        //      hpv_m(g, 3, 4) = sum(popn(hpv123_id == 1 & sexm_id == 1 & art_id == 1 & age_id == 4));
        //      hpv_m(g, 3, 5) = sum(popn(hpv123_id == 1 & sexm_id == 1 & art_id == 1 & age_id == 5));

        //        cinpr_f(g, 1, 1) = sum(popn(CIN_id == 1 & sexf_id == 1 & hivS_id == 1 & age_id == 1));
        //        cinpr_f(g, 1, 2) = sum(popn(CIN_id == 1 & sexf_id == 1 & hivS_id == 1 & age_id == 2));
        //        cinpr_f(g, 1, 3) = sum(popn(CIN_id == 1 & sexf_id == 1 & hivS_id == 1 & age_id == 3));
        //        cinpr_f(g, 1, 4) = sum(popn(CIN_id == 1 & sexf_id == 1 & hivS_id == 1 & age_id == 4));
        //        cinpr_f(g, 1, 5) = sum(popn(CIN_id == 1 & sexf_id == 1 & hivS_id == 1 & age_id == 5));

        //        cinpr_f(g, 2, 1) = sum(popn(CIN_id == 1 & sexf_id == 1 & noart_id == 1 & age_id == 1));
        //        cinpr_f(g, 2, 2) = sum(popn(CIN_id == 1 & sexf_id == 1 & noart_id == 1 & age_id == 2));
        //        cinpr_f(g, 2, 3) = sum(popn(CIN_id == 1 & sexf_id == 1 & noart_id == 1 & age_id == 3));
        //        cinpr_f(g, 2, 4) = sum(popn(CIN_id == 1 & sexf_id == 1 & noart_id == 1 & age_id == 4));
        //        cinpr_f(g, 2, 5) = sum(popn(CIN_id == 1 & sexf_id == 1 & noart_id == 1 & age_id == 5));

        //        cinpr_f(g, 3, 1) = sum(popn(CIN_id == 1 & sexf_id == 1 & art_id == 1 & age_id == 1));
        //        cinpr_f(g, 3, 2) = sum(popn(CIN_id == 1 & sexf_id == 1 & art_id == 1 & age_id == 2));
        //        cinpr_f(g, 3, 3) = sum(popn(CIN_id == 1 & sexf_id == 1 & art_id == 1 & age_id == 3));
        //        cinpr_f(g, 3, 4) = sum(popn(CIN_id == 1 & sexf_id == 1 & art_id == 1 & age_id == 4));
        //        cinpr_f(g, 3, 5) = sum(popn(CIN_id == 1 & sexf_id == 1 & art_id == 1 & age_id == 5));

        // New ones:

        // ccpr_f(g, 1, 1) = sum(popn(CC_id == 1 & sexf_id == 1 & hivS_id == 1 & age_id == 1));
        // ccpr_f(g, 1, 2) = sum(popn(CC_id == 1 & sexf_id == 1 & hivS_id == 1 & age_id == 2));
        // ccpr_f(g, 1, 3) = sum(popn(CC_id == 1 & sexf_id == 1 & hivS_id == 1 & age_id == 3));
        // ccpr_f(g, 1, 4) = sum(popn(CC_id == 1 & sexf_id == 1 & hivS_id == 1 & age_id == 4));
        // ccpr_f(g, 1, 5) = sum(popn(CC_id == 1 & sexf_id == 1 & hivS_id == 1 & age_id == 5));

        // ccpr_f(g, 2, 1) = sum(popn(CC_id == 1 & sexf_id == 1 & noart_id == 1 & age_id == 1));
        // ccpr_f(g, 2, 2) = sum(popn(CC_id == 1 & sexf_id == 1 & noart_id == 1 & age_id == 2));
        // ccpr_f(g, 2, 3) = sum(popn(CC_id == 1 & sexf_id == 1 & noart_id == 1 & age_id == 3));
        // ccpr_f(g, 2, 4) = sum(popn(CC_id == 1 & sexf_id == 1 & noart_id == 1 & age_id == 4));
        // ccpr_f(g, 2, 5) = sum(popn(CC_id == 1 & sexf_id == 1 & noart_id == 1 & age_id == 5));

        // ccpr_f(g, 3, 1) = sum(popn(CC_id == 1 & sexf_id == 1 & art_id == 1 & age_id == 1));
        // ccpr_f(g, 3, 2) = sum(popn(CC_id == 1 & sexf_id == 1 & art_id == 1 & age_id == 2));
        // ccpr_f(g, 3, 3) = sum(popn(CC_id == 1 & sexf_id == 1 & art_id == 1 & age_id == 3));
        // ccpr_f(g, 3, 4) = sum(popn(CC_id == 1 & sexf_id == 1 & art_id == 1 & age_id == 4));
        // ccpr_f(g, 3, 5) = sum(popn(CC_id == 1 & sexf_id == 1 & art_id == 1 & age_id == 5));

        LOOP_ALL {
          if IS(S->hpv12_id, 1) {
            if IS(S->sexf_id, 1) {
              if IS(S->hivS_id, 1) R->hpv9vt_f[g][1][GET_ALL(S->age_id2)] += GET_ALL(popn);
              if IS(S->noart_id, 1) R->hpv9vt_f[g][2][GET_ALL(S->age_id2)] += GET_ALL(popn);
              if IS(S->art_id, 1) R->hpv9vt_f[g][3][GET_ALL(S->age_id2)] += GET_ALL(popn);
            }
            if IS(S->sexm_id, 1) {
              if IS(S->hivS_id, 1) R->hpv9vt_m[g][1][GET_ALL(S->age_id2)] += GET_ALL(popn);
              if IS(S->noart_id, 1) R->hpv9vt_m[g][2][GET_ALL(S->age_id2)] += GET_ALL(popn);
              if IS(S->art_id, 1) R->hpv9vt_m[g][3][GET_ALL(S->age_id2)] += GET_ALL(popn);
            }
          }
          if IS(S->hpv123_id, 1) {
            if IS(S->sexf_id, 1) {
              if IS(S->hivS_id, 1) R->hpv_f[g][1][GET_ALL(S->age_id2)] += GET_ALL(popn);
              if IS(S->noart_id, 1) R->hpv_f[g][2][GET_ALL(S->age_id2)] += GET_ALL(popn);
              if IS(S->art_id, 1) R->hpv_f[g][3][GET_ALL(S->age_id2)] += GET_ALL(popn);
            }
            if IS(S->sexm_id, 1) {
              if IS(S->hivS_id, 1) R->hpv_m[g][1][GET_ALL(S->age_id2)] += GET_ALL(popn);
              if IS(S->noart_id, 1) R->hpv_m[g][2][GET_ALL(S->age_id2)] += GET_ALL(popn);
              if IS(S->art_id, 1) R->hpv_m[g][3][GET_ALL(S->age_id2)] += GET_ALL(popn);
            }
          }
          if IS(S->CIN_id, 1) {
            if IS(S->sexf_id, 1) {
              if IS(S->hivS_id, 1) R->cinpr_f[g][1][GET_ALL(S->age_id2)] += GET_ALL(popn);
              if IS(S->noart_id, 1) R->cinpr_f[g][2][GET_ALL(S->age_id2)] += GET_ALL(popn);
              if IS(S->art_id, 1) R->cinpr_f[g][3][GET_ALL(S->age_id2)] += GET_ALL(popn);
            }
          }
          if (IS(S->CC_id, 1)) {
            if (IS(S->sexf_id, 1)) {
              if IS(S->hivS_id, 1) R->ccpr_f[g][1][GET_ALL(S->age_id2)] += GET_ALL(popn);
              if IS(S->noart_id, 1) R->ccpr_f[g][2][GET_ALL(S->age_id2)] += GET_ALL(popn);
              if IS(S->art_id, 1) R->ccpr_f[g][3][GET_ALL(S->age_id2)] += GET_ALL(popn);
            }
          }
        }

        //testtime1(g, :) = n;

      R->testtime1yr[g] = n;

      // g = g + 1;
      g++;
      //end
    }

    //    %%

    //      h5tyrid = h5yrid(ggg);
    //      if n == h5tyrid
    int h5tyrid = (int)round(h5yrid[ggg]);
    if (n == h5tyrid) {

      double _top1 = 0.0;
      double _bottom1 = 0.0;
      double _top2 = 0.0;
      double _bottom2 = 0.0;

      //        popbyrisk_f(ggg, :) = reshape(sum(sum(sum(sum(sum(sum(popn(:, : , : , : , 1, 1 : r, : , : ), 1), 2), 3), 4), 7), 8), 1, r) . /
      //                                  sum(sum(sum(sum(sum(sum(sum(popn(:, : , : , : , 1, :    , : , : ), 1), 2), 3), 4), 6), 7), 8);
      //      popbyrisk_m(ggg, :) = reshape(sum(sum(sum(sum(sum(sum(popn(:, : , : , : , 2, 1 : r, : , : ), 1), 2), 3), 4), 7), 8), 1, r) . /
      //                                sum(sum(sum(sum(sum(sum(sum(popn(:, : , : , : , 2, : , : , : ), 1), 2), 3), 4), 6), 7), 8);

      LOOP_R LOOP_A LOOP_V1{
        _bottom1 += GET(popn_sum_h1_h2_h3_i, 1, 1, 1, 1, 1, _dr, _da, _dv1);
        _bottom2 += GET(popn_sum_h1_h2_h3_i, 1, 1, 1, 1, 2, _dr, _da, _dv1);
      }

        LOOP_R {
          _top1 = 0.0;
          _top2 = 0.0;
          LOOP_A LOOP_V1 {
            _top1 += GET(popn_sum_h1_h2_h3_i, 1, 1, 1, 1, 1, _dr, _da, _dv1);
            _top2 += GET(popn_sum_h1_h2_h3_i, 1, 1, 1, 1, 2, _dr, _da, _dv1);
          }
          R->popbyrisk_f[ggg][_dr] = _top1 / _bottom1;
          R->popbyrisk_m[ggg][_dr] = _top2 / _bottom2;
      }

        //    ggg = ggg + 1;
        //    end
      ggg++;
    }

    //    if n == tstep * (1980 - dp0(22))
    //        pop1980 = popn;
    //    end

    if (n == tstep * (2000 - p->dp0[31])) LOOP_LINEAR GET_LIN_ALL(R->pop2000) = GET_LIN_ALL(popn);

    //      if n == tstep * (2018 - dp0(22))
    //          pop2018 = popn;
    //    end

    if (n == tstep * (2018 - p->dp0[31])) LOOP_LINEAR GET_LIN_ALL(R->pop2018) = GET_LIN_ALL(popn);

    //   if n == tstep * (2120 - dp0(22))
    //        pop2120 = popn;
    //    end

    if (n == tstep * (2120 - p->dp0[31])) LOOP_LINEAR GET_LIN_ALL(R->pop2120) = GET_LIN_ALL(popn);

  } // End of n-iterations loop
  // timer0 += clock();
  // printf("Done in time %d ms. Avg per step: %f\n", timer0, (float) timer0 / (float) niterations);
  // fflush(stdout);

  // hpv1cum_pop = hpv1cum_popt(yrid, :);
  // hpv2cum_pop = hpv2cum_popt(yrid, :);
  // hpv3cum_pop = hpv3cum_popt(yrid, :);
  // hpv1hncum_pop = hpv1hncum_popt(yrid, :);
  // hpv2hncum_pop = hpv2hncum_popt(yrid, :);
  // hpv3hncum_pop = hpv3hncum_popt(yrid, :);

  // Above is an implicit reshape from (n,s,a,v1) into (n,s*a*v1) - with
  // inner-loop s, then a, and last v1.

  for (int i = 1; i <= tend + 1; i++) {
    int dindex = 1;
    const int yindex = (int)round(yrid[i]);
    LOOP_V1 LOOP_A LOOP_S {
      const int index = GET_4D_NSAV1_INDEX(yindex, _ds, _da, _dv1);
      R->hpv1cum_pop[i][dindex] = hpv1cum_popt[index];
      R->hpv2cum_pop[i][dindex] = hpv2cum_popt[index];
      R->hpv3cum_pop[i][dindex] = hpv3cum_popt[index];
      R->hpv1hncum_pop[i][dindex] = hpv1hncum_popt[index];
      R->hpv2hncum_pop[i][dindex] = hpv2hncum_popt[index];
      R->hpv3hncum_pop[i][dindex++] = hpv3hncum_popt[index];
    }
  }




  // hivd_f = hivd_ft(yrid, :);
  // hivd_m = hivd_mt(yrid, :);
  // ccd_f = ccd_ft(yrid, :);
  // ccd_f_neg = ccd_ft_neg(yrid, :);
  // ccd_f_pos = ccd_ft_pos(yrid, :);
  // ccd_f_art = ccd_ft_art(yrid, :);
  // othd_f = othd_ft(yrid, :);
  // othd_m = othd_mt(yrid, :);
  // ccinc_f = ccinc_ft(yrid, :);
  // ccinc_f_neg = ccinc_ft_neg(yrid, :);
  // ccinc_f_pos = ccinc_ft_pos(yrid, :);
  // ccinc_f_art = ccinc_ft_art(yrid, :);
  // ccinc_nvtf = ccinc_nvtft(yrid, :);
  // ccinc_nvtf_neg = ccinc_nvtft_neg(yrid, :);
  // ccinc_nvtf_pos = ccinc_nvtft_pos(yrid, :);
  // ccinc_nvtf_art = ccinc_nvtft_art(yrid, :);
  // cinin_f = cin_ft(yrid, :);
  // cinin_f_neg = cin_ft_neg(yrid, :);
  // cinin_f_pos = cin_ft_pos(yrid, :);
  // cinin_f_art = cin_ft_art(yrid, :);

  int* yrid_a = new int[3]{ 0, tend + 1, _SIZE_A };
  R->hivd_f = CREATE_DMATRIX(yrid_a, R->hivd_f, 0)
  R->hivd_m = CREATE_DMATRIX(yrid_a, R->hivd_m, 0)
  R->ccd_f = CREATE_DMATRIX(yrid_a, R->ccd_f, 0)
  R->ccd_f_neg = CREATE_DMATRIX(yrid_a, R->ccd_f_neg, 0)
  R->ccd_f_pos = CREATE_DMATRIX(yrid_a, R->ccd_f_pos, 0)
  R->ccd_f_art = CREATE_DMATRIX(yrid_a, R->ccd_f_art, 0)
  R->othd_f = CREATE_DMATRIX(yrid_a, R->othd_f, 0)
  R->othd_m = CREATE_DMATRIX(yrid_a, R->othd_m, 0)
  R->ccinc_f = CREATE_DMATRIX(yrid_a, R->ccinc_f, 0)
  R->ccinc_f_neg = CREATE_DMATRIX(yrid_a, R->ccinc_f_neg, 0)
  R->ccinc_f_pos = CREATE_DMATRIX(yrid_a, R->ccinc_f_pos, 0)
  R->ccinc_f_art = CREATE_DMATRIX(yrid_a, R->ccinc_f_art, 0)
  R->ccinc_nvtf = CREATE_DMATRIX(yrid_a, R->ccinc_nvtf, 0)
  R->ccinc_nvtf_neg = CREATE_DMATRIX(yrid_a, R->ccinc_nvtf_neg, 0)
  R->ccinc_nvtf_pos = CREATE_DMATRIX(yrid_a, R->ccinc_nvtf_pos, 0)
  R->ccinc_nvtf_art = CREATE_DMATRIX(yrid_a, R->ccinc_nvtf_art, 0)
  R->cinin_f = CREATE_DMATRIX(yrid_a, R->cinin_f, 0)
  R->cinin_f_neg = CREATE_DMATRIX(yrid_a, R->cinin_f_neg, 0)
  R->cinin_f_pos = CREATE_DMATRIX(yrid_a, R->cinin_f_pos, 0)
  R->cinin_f_art = CREATE_DMATRIX(yrid_a, R->cinin_f_art, 0)

  for (int j = 1; j <= _SIZE_A; j++) {
    for (int i = 1; i <= tend + 1; i++) {
      int index = GET_2D_NA_INDEX((int)round(yrid[i]), j);
      R->hivd_f[i][j] = hivd_ft[index];
      R->hivd_m[i][j] = hivd_mt[index];
      R->ccd_f[i][j] = ccd_ft[index];
      R->ccd_f_neg[i][j] = ccd_ft_neg[index];
      R->ccd_f_pos[i][j] = ccd_ft_pos[index];
      R->ccd_f_art[i][j] = ccd_ft_art[index];
      R->othd_f[i][j] = othd_ft[index];
      R->othd_m[i][j] = othd_mt[index];
      R->ccinc_f[i][j] = ccinc_ft[index];
      R->ccinc_f_neg[i][j] = ccinc_ft_neg[index];
      R->ccinc_f_pos[i][j] = ccinc_ft_pos[index];
      R->ccinc_f_art[i][j] = ccinc_ft_art[index];
      R->ccinc_nvtf[i][j] = ccinc_nvtft[index];
      R->ccinc_nvtf_neg[i][j] = ccinc_nvtft_neg[index];
      R->ccinc_nvtf_pos[i][j] = ccinc_nvtft_pos[index];
      R->ccinc_nvtf_art[i][j] = ccinc_nvtft_art[index];
      R->cinin_f[i][j] = cin_ft[index];
      R->cinin_f_neg[i][j] = cin_ft_neg[index];
      R->cinin_f_pos[i][j] = cin_ft_pos[index];
      R->cinin_f_art[i][j] = cin_ft_art[index];
    }
  }

  // hivinc_pop_cum = hivinc_pop_cumt(yrid, :, : , : );

  int* yrid_sav1 = new int[5]{ 0, tend + 1, _SIZE_S, _SIZE_A, _SIZE_V1};
  R->hivinc_pop_cum = CREATE_DMATRIX4D(yrid_sav1, R->hivinc_pop_cum, 0.0)

    LOOP_S LOOP_A LOOP_V1
      for (int i = 1; i <= tend + 1; i++)
        R->hivinc_pop_cum[i][_ds][_da][_dv1] = hivinc_pop_cumt[GET_4D_NSAV1_INDEX((int)round(yrid[i]), _ds, _da, _dv1)];

  // screencov_yr = screencov(yrid, :, : );

  int* yrid_3_2 = new int[4]{ 0, tend + 1, 3, 2 };
  R->screencov_neg_yr = CREATE_DMATRIX3D(yrid_3_2, R->screencov_neg_yr, 0.0);
  for (int i = 1; i <= tend + 1; i++)
    for (int j = 1; j <= 3; j++)
      for (int k = 1; k <= 2; k++)
        R->screencov_neg_yr[i][j][k] = screencov_neg[GET_3D_N32_INDEX((int)round(yrid[i]), j, k)];

 // int* yrid_3_2 = new int[4]{ 0, tend + 1, 3, 2 };
  R->screencov_pos_yr = CREATE_DMATRIX3D(yrid_3_2, R->screencov_pos_yr, 0.0);
  for (int i = 1; i <= tend + 1; i++)
    for (int j = 1; j <= 3; j++)
      for (int k = 1; k <= 2; k++)
        R->screencov_pos_yr[i][j][k] = screencov_pos[GET_3D_N32_INDEX((int)round(yrid[i]), j, k)];

 // int* yrid_3_2 = new int[4]{ 0, tend + 1, 3, 2 };
  R->screencov_art_yr = CREATE_DMATRIX3D(yrid_3_2, R->screencov_art_yr, 0.0);
  for (int i = 1; i <= tend + 1; i++)
    for (int j = 1; j <= 3; j++)
      for (int k = 1; k <= 2; k++)
        R->screencov_art_yr[i][j][k] = screencov_art[GET_3D_N32_INDEX((int)round(yrid[i]), j, k)];



  return R;

}
