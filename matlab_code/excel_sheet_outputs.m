
 addpath('/Users/minttu/Dropbox/work files - continuous use/Human Papilloma Virus/Dynamic HPVx2+HIV+vacc model') 
 addpath('/Users/minttu/hpv_hiv_matlab' ) 
 addpath('/Users/minttu/matlab_figures' ) 

% 
sc0 = load('sc0-2021-05-14_17:24.mat');
sc00 = load('sc00-2021-05-15_16:03.mat');
sc0b = load('S0b-2021-05-15_06:48.mat');
sc0c = load('S0c-2021-05-15_23:34.mat');

sc1 = load('sc1-2021-05-14.mat');
sc2 = load('sc2-2021-05-14_17:22.mat');
sc3 = load('sc3-2021-05-14_17:43.mat');
sc4 = load('sc4-2021-05-14_18:25.mat');
sc5 = load('sc5-2021-05-14_17:54.mat');
sc6 = load('sc6-2021-05-14_17:33.mat');
sc7 = load('sc7-2021-05-14_18:26.mat');
sc7a = load('sc7a-2021-05-15_03:11.mat');
sc7b = load('sc7b-2021-05-15_02:00.mat');
sc8 = load('sc8-2021-05-14_18:26.mat');
sc9 = load('sc9-2021-05-14_17:37.mat');
sc10 = load('sc10-2021-05-14_17:27.mat');
sc11 = load('sc11-2021-05-14_18:18.mat');
%%%%%%%%%%%%%%%%%%%%%%
sc00 = load('sc00-2021-05-15_16:03.mat');

sc0c = load('S0c-2021-05-15_23:34.mat');


pop2015 = [ ...   
              325428
              311262
              295693
              287187
              291738
              299655
              272348
              247167
              240167
              226750
              201603
              171975
              150562
              113118
              82266
              64484
              42237
              23477
              9261
              2155 ];

%sc0.res = res;

yr1985=1985-1950+1;
yr1990=1990-1950+1;
yr2020=2020-1950+1;

% lighter
col1=0.7*[0.8  0.8  0.8];

dat1 = [1 .6 .6];
col2=0.6* [0 1 1];%dark cyan
col3=0.7* [0 0 1];%dark blue

col4=[0 0.5 0.5]*0.1;%[0 .5 .5];
col5=[0 0.5 0.5];

set(0,'DefaultFigureWindowStyle','docked')
cmap = colormap(lines);

%% baseline - Sc0
[sc0.cc_incid_all, sc0.cc_incid_neg, sc0.cc_incid_posall, sc0.cc_incid_art, sc0.cc_incid_pos] = ccincid(sc0.res); % pp

sc0.stall = strd2015_results(sc0.cc_incid_all, pop2015, 2);
sc0.stneg = strd2015_results(sc0.cc_incid_neg, pop2015, 2);
sc0.stpos = strd2015_results(sc0.cc_incid_posall, pop2015, 37);
sc0.start = strd2015_results(sc0.cc_incid_art, pop2015, 56);
sc0.stnart = strd2015_results(sc0.cc_incid_pos, pop2015, 37);

sc0.tab_ccall = res_table(sc0.stall(:, 69:end), sc0.cc_incid_all(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S0.xlsx"], "Pop(All) ICC");            
sc0.tab_ccneg = res_table(sc0.stneg(:, 69:end), sc0.cc_incid_neg(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S0.xlsx"], "HIV- (ICC)");
sc0.tab_ccpos = res_table(sc0.stpos(:, 34:end), sc0.cc_incid_posall(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S0.xlsx"], "HIV+ (ICC)");
sc0.tab_ccart = res_table(sc0.start(:, 15:end), sc0.cc_incid_art(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S0.xlsx"], "HIV+ ART (ICC)");
sc0.tab_ccnart = res_table(sc0.stnart(:, 34:end), sc0.cc_incid_pos(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S0.xlsx"], "HIV+ no ART (ICC)");

sc0.tab_pops = res_poptable(sc0.res, ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S0.xlsx"], 102);

sc0.tab_ccc = res_ccctable(sc0.res, ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S0.xlsx"]);

sc0.cum.stall = strd2015_numcum_results(sc0.cc_incid_all(:, 70:end, :), pop2015, ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S0.xlsx"], "Pop(All) CCC");
sc0.cum.stneg = strd2015_numcum_results(sc0.cc_incid_neg(:, 70:end, :), pop2015, ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S0.xlsx"], "HIV- (CCC)");
sc0.cum.stpos = strd2015_numcum_results(sc0.cc_incid_posall(:, 70:end, :), pop2015, ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S0.xlsx"], "HIV+ (CCC)" );
sc0.cum.start = strd2015_numcum_results(sc0.cc_incid_art(:, 70:end, :), pop2015, ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S0.xlsx"], "HIV+ ART (CCC)" );
sc0.cum.stnart = strd2015_numcum_results(sc0.cc_incid_pos(:, 70:end, :), pop2015, ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S0.xlsx"], "HIV+ no ART (CCC)" );

% 5% LUB

sc0.tab_ccall_p5 = res_table_prctile(sc0.stall(:, 70:end), sc0.cc_incid_all(:, 71:end, :), 5, ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)LUB_S0.xlsx"], "Pop(All) ICC");            
sc0.tab_ccneg_p5 = res_table_prctile(sc0.stneg(:, 70:end), sc0.cc_incid_neg(:, 71:end, :), 5,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)LUB_S0.xlsx"], "HIV- (ICC)");
sc0.tab_ccpos_p5 = res_table_prctile(sc0.stpos(:, 35:end), sc0.cc_incid_posall(:, 71:end, :), 5,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)LUB_S0.xlsx"], "HIV+ (ICC)");
sc0.tab_ccart_p5 = res_table_prctile(sc0.start(:, 16:end), sc0.cc_incid_art(:, 71:end, :), 5,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)LUB_S0.xlsx"], "HIV+ ART (ICC)");
sc0.tab_ccnart_p5 = res_table_prctile(sc0.stnart(:, 35:end), sc0.cc_incid_pos(:, 71:end, :), 5,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)LUB_S0.xlsx"], "HIV+ no ART (ICC)");

sc0.tab_pops_p5 = res_poptable_prctile(sc0.res, 5, ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)LUB_S0.xlsx"], 101);


% 95% LUB

sc0.tab_ccall_p95 = res_table_prctile(sc0.stall(:, 70:end), sc0.cc_incid_all(:, 71:end, :), 95, ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)HUB_S0.xlsx"], "Pop(All) ICC");            
sc0.tab_ccneg_p95 = res_table_prctile(sc0.stneg(:, 70:end), sc0.cc_incid_neg(:, 71:end, :), 95,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)HUB_S0.xlsx"], "HIV- (ICC)");
sc0.tab_ccpos_p95 = res_table_prctile(sc0.stpos(:, 35:end), sc0.cc_incid_posall(:, 71:end, :), 95,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)HUB_S0.xlsx"], "HIV+ (ICC)");
sc0.tab_ccart_p95 = res_table_prctile(sc0.start(:, 16:end), sc0.cc_incid_art(:, 71:end, :), 95,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)HUB_S0.xlsx"], "HIV+ ART (ICC)");
sc0.tab_ccnart_p95 = res_table_prctile(sc0.stnart(:, 35:end), sc0.cc_incid_pos(:, 71:end, :), 95,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)HUB_S0.xlsx"], "HIV+ no ART (ICC)");

sc0.tab_pops_p95 = res_poptable_prctile(sc0.res, 95, ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)HUB_S0.xlsx"], 101);

% hr-hpv women

[sc0.hrhpvf_all, sc0.hrhpvf_neg, sc0.hrhpvf_posall, sc0.hrhpvf_art, sc0.hrhpvf_pos, ...
 sc0.hrhpvf_all_num, sc0.hrhpvf_neg_num, sc0.hrhpvf_posall_num, sc0.hrhpvf_art_num, sc0.hrhpvf_pos_num] = hrhpvprev(sc0.res); % pp

sc0.stall_hrhpvf = strd2015_results_num(sc0.hrhpvf_all, pop2015, 2);
sc0.stneg_hrhpvf = strd2015_results_num(sc0.hrhpvf_neg, pop2015, 2);
sc0.stpos_hrhpvf = strd2015_results_num(sc0.hrhpvf_posall, pop2015, 37);
sc0.start_hrhpvf = strd2015_results_num(sc0.hrhpvf_art, pop2015, 56);
sc0.stnart_hrhpvf = strd2015_results_num(sc0.hrhpvf_pos, pop2015, 37);

% here weird spacing
sc0.tab_all_hrhpvf = res_table_prctile_num(sc0.stall_hrhpvf(:, 70:end), sc0.hrhpvf_all_num(:, 71:end, :), 50, ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_hrHPV_Prevalence-standardised-(2020-2120)_S0.xlsx"], "Pop(All) (hrHPV)");            
sc0.tab_neg_hrhpvf = res_table_prctile_num(sc0.stneg_hrhpvf(:, 70:end), sc0.hrhpvf_neg_num(:, 71:end, :), 50,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_hrHPV_Prevalence-standardised-(2020-2120)_S0.xlsx"], "HIV-  (hrHPV)");
sc0.tab_pos_hrhpvf = res_table_prctile_num(sc0.stpos_hrhpvf(:, 35:end), sc0.hrhpvf_posall_num(:, 71:end, :), 50,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_hrHPV_Prevalence-standardised-(2020-2120)_S0.xlsx"], "HIV+   (hrHPV)"); % this has 3 spaces
sc0.tab_art_hrhpvf = res_table_prctile_num(sc0.start_hrhpvf(:, 16:end), sc0.hrhpvf_art_num(:, 71:end, :), 50,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_hrHPV_Prevalence-standardised-(2020-2120)_S0.xlsx"], "HIV+ ART  (hrHPV)");
sc0.tab_nart_hrhpvf = res_table_prctile_num(sc0.stnart_hrhpvf(:, 35:end), sc0.hrhpvf_pos_num(:, 71:end, :), 50,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_hrHPV_Prevalence-standardised-(2020-2120)_S0.xlsx"], "HIV+ no ART  (hrHPV)");

sc0.tab_pops_hrhpvf = res_poptable_prctile(sc0.res, 50, ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_hrHPV_Prevalence-standardised-(2020-2120)_S0.xlsx"], 101);

% [sc0.vacc, sc0.screen] = coverageest(sc0.res);

% VT-hpv women

[sc0.hpv9vtf_all, sc0.hpv9vtf_neg, sc0.hpv9vtf_posall, sc0.hpv9vtf_art, sc0.hpv9vtf_pos, ...
 sc0.hpv9vtf_all_num, sc0.hpv9vtf_neg_num, sc0.hpv9vtf_posall_num, sc0.hpv9vtf_art_num, sc0.hpv9vtf_pos_num] = vthpvprev(sc0.res); % pp

sc0.stall_hpv9vtf = strd2015_results_num(sc0.hpv9vtf_all, pop2015, 2);
sc0.stneg_hpv9vtf = strd2015_results_num(sc0.hpv9vtf_neg, pop2015, 2);
sc0.stpos_hpv9vtf = strd2015_results_num(sc0.hpv9vtf_posall, pop2015, 37);
sc0.start_hpv9vtf = strd2015_results_num(sc0.hpv9vtf_art, pop2015, 56);
sc0.stnart_hpv9vtf = strd2015_results_num(sc0.hpv9vtf_pos, pop2015, 37);

% here weird spacing
sc0.tab_all_hpv9vtf = res_table_prctile_num(sc0.stall_hpv9vtf(:, 70:end), sc0.hpv9vtf_all_num(:, 71:end, :), 50, ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_VThrHPV_Prevalence-standardised-(2020-2120)_S0.xlsx"], "Pop(All) (vtHPV)");            
sc0.tab_neg_hpv9vtf = res_table_prctile_num(sc0.stneg_hpv9vtf(:, 70:end), sc0.hpv9vtf_neg_num(:, 71:end, :), 50,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_VThrHPV_Prevalence-standardised-(2020-2120)_S0.xlsx"], "HIV-  (vtHPV)");
sc0.tab_pos_hpv9vtf = res_table_prctile_num(sc0.stpos_hpv9vtf(:, 35:end), sc0.hpv9vtf_posall_num(:, 71:end, :), 50,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_VThrHPV_Prevalence-standardised-(2020-2120)_S0.xlsx"], "HIV+   (vtHPV)"); % this has 3 spaces
sc0.tab_art_hpv9vtf = res_table_prctile_num(sc0.start_hpv9vtf(:, 16:end), sc0.hpv9vtf_art_num(:, 71:end, :), 50,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_VThrHPV_Prevalence-standardised-(2020-2120)_S0.xlsx"], "HIV+ ART  (vtHPV)");
sc0.tab_nart_hpv9vtf = res_table_prctile_num(sc0.stnart_hpv9vtf(:, 35:end), sc0.hpv9vtf_pos_num(:, 71:end, :), 50,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_VThrHPV_Prevalence-standardised-(2020-2120)_S0.xlsx"], "HIV+ no ART  (vtHPV)");

sc0.tab_pops_hpv9vtf = res_poptable_prctile(sc0.res, 50, ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_VThrHPV_Prevalence-standardised-(2020-2120)_S0.xlsx"], 101);

% CIN2+ prevalence women

[sc0.cinprf_all, sc0.cinprf_neg, sc0.cinprf_posall, sc0.cinprf_art, sc0.cinprf_pos, ...
 sc0.cinprf_all_num, sc0.cinprf_neg_num, sc0.cinprf_posall_num, sc0.cinprf_art_num, sc0.cinprf_pos_num] = cin2prev(sc0.res); % pp

sc0.stall_cinprf = strd2015_results_num(sc0.cinprf_all, pop2015, 2);
sc0.stneg_cinprf = strd2015_results_num(sc0.cinprf_neg, pop2015, 2);
sc0.stpos_cinprf = strd2015_results_num(sc0.cinprf_posall, pop2015, 37);
sc0.start_cinprf = strd2015_results_num(sc0.cinprf_art, pop2015, 56);
sc0.stnart_cinprf = strd2015_results_num(sc0.cinprf_pos, pop2015, 37);

% here weird spacing
sc0.tab_all_cinprf = res_table_prctile_num(sc0.stall_cinprf(:, 70:end), sc0.cinprf_all_num(:, 71:end, :), 50, ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CIN2+_Prevalence-standardised-(2020-2120)_S0.xlsx"], "Pop(All) (CIN2+)");            
sc0.tab_neg_cinprf = res_table_prctile_num(sc0.stneg_cinprf(:, 70:end), sc0.cinprf_neg_num(:, 71:end, :), 50,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CIN2+_Prevalence-standardised-(2020-2120)_S0.xlsx"], "HIV- (CIN2+)"); % 1 space only
sc0.tab_pos_cinprf = res_table_prctile_num(sc0.stpos_cinprf(:, 35:end), sc0.cinprf_posall_num(:, 71:end, :), 50,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CIN2+_Prevalence-standardised-(2020-2120)_S0.xlsx"], "HIV+  (CIN2+)"); % this has 2 spaces
sc0.tab_art_cinprf = res_table_prctile_num(sc0.start_cinprf(:, 16:end), sc0.cinprf_art_num(:, 71:end, :), 50,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CIN2+_Prevalence-standardised-(2020-2120)_S0.xlsx"], "HIV+ ART  (CIN2+)");
sc0.tab_nart_cinprf = res_table_prctile_num(sc0.stnart_cinprf(:, 35:end), sc0.cinprf_pos_num(:, 71:end, :), 50,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CIN2+_Prevalence-standardised-(2020-2120)_S0.xlsx"], "HIV+ no ART  (CIN2+)");

sc0.tab_pops_cinprf = res_poptable_prctile(sc0.res, 50, ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CIN2+_Prevalence-standardised-(2020-2120)_S0.xlsx"], 101);

%%%%%%
% pre-impact

sc0.stpos2 = strd2015_results(sc0.cc_incid_posall, pop2015, 2);
sc0.start2 = strd2015_results(sc0.cc_incid_art, pop2015, 2);
sc0.stnart2 = strd2015_results(sc0.cc_incid_pos, pop2015, 2);

sc0.stpos2(isnan(sc0.stpos2))=0;
sc0.start2(isnan(sc0.start2))=0;
sc0.stnart2(isnan(sc0.stnart2))=0;

sc0.tab_ccall = res_pretable(sc0.stall(:, yr1990-1:end), sc0.cc_incid_all(:, yr1990:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Pre-Impact_CC_IncidenceRates-standardised-(Before_2020)_S0f.xlsx"], "Pop(All) ICC");            
sc0.tab_ccneg = res_pretable(sc0.stneg(:, yr1990-1:end), sc0.cc_incid_neg(:, yr1990:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Pre-Impact_CC_IncidenceRates-standardised-(Before_2020)_S0f.xlsx"], "HIV- (ICC)");
            
sc0.tab_ccpos = res_pretable(sc0.stpos2(:, yr1990-1:end), sc0.cc_incid_posall(:, yr1990:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Pre-Impact_CC_IncidenceRates-standardised-(Before_2020)_S0f.xlsx"], "HIV+ (ICC)");
sc0.tab_ccart = res_pretable(sc0.start2(:, yr1990-1:end), sc0.cc_incid_art(:, yr1990:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Pre-Impact_CC_IncidenceRates-standardised-(Before_2020)_S0f.xlsx"], "HIV+ ART (ICC)");
sc0.tab_ccnart = res_pretable(sc0.stnart2(:, yr1990-1:end), sc0.cc_incid_pos(:, yr1990:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Pre-Impact_CC_IncidenceRates-standardised-(Before_2020)_S0f.xlsx"], "HIV+ no ART (ICC)");

sc0.tab_pops = res_poptable(sc0.res, ["/Users/minttu/who_results_sheets/ICL_HSPH_Pre-Impact_CC_IncidenceRates-standardised-(Before_2020)_S0f.xlsx"], 132);

sc0.tab_coverage = res_covtable(sc0.res, ["/Users/minttu/who_results_sheets/ICL_HSPH_Pre-Impact_CC_IncidenceRates-standardised-(Before_2020)_S0f.xlsx"], 101+30);


%[sc0.vacc, sc0.vacc_cov_all, sc0.screen] = coverageest_nums(sc0.res);



%% baseline - Sc00
[sc00.cc_incid_all, sc00.cc_incid_neg, sc00.cc_incid_posall, sc00.cc_incid_art, sc00.cc_incid_pos] = ccincid(sc00.res); % pp

sc00.stall = strd2015_results(sc00.cc_incid_all, pop2015, 2);
sc00.stneg = strd2015_results(sc00.cc_incid_neg, pop2015, 2);
sc00.stpos = strd2015_results(sc00.cc_incid_posall, pop2015, 37);
sc00.start = strd2015_results(sc00.cc_incid_art, pop2015, 56);
sc00.stnart = strd2015_results(sc00.cc_incid_pos, pop2015, 37);

sc00.tab_ccall = res_table(sc00.stall(:, 69:end), sc00.cc_incid_all(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S00.xlsx"], "Pop(All) ICC");            
sc00.tab_ccneg = res_table(sc00.stneg(:, 69:end), sc00.cc_incid_neg(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S00.xlsx"], "HIV- (ICC)");
sc00.tab_ccpos = res_table(sc00.stpos(:, 34:end), sc00.cc_incid_posall(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S00.xlsx"], "HIV+ (ICC)");
sc00.tab_ccart = res_table(sc00.start(:, 15:end), sc00.cc_incid_art(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S00.xlsx"], "HIV+ ART (ICC)");
sc00.tab_ccnart = res_table(sc00.stnart(:, 34:end), sc00.cc_incid_pos(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S00.xlsx"], "HIV+ no ART (ICC)");

sc00.tab_pops = res_poptable(sc00.res, ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S00.xlsx"], 102);

sc00.tab_ccc = res_ccctable(sc00.res, ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S00.xlsx"]);


sc00.cum.stall = strd2015_numcum_results(sc00.cc_incid_all(:, 70:end, :), pop2015,  ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S00.xlsx"], "Pop(All) CCC");
sc00.cum.stneg = strd2015_numcum_results(sc00.cc_incid_neg(:, 70:end, :), pop2015,  ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S00.xlsx"], "HIV- (CCC)");
sc00.cum.stpos = strd2015_numcum_results(sc00.cc_incid_posall(:, 70:end, :), pop2015,  ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S00.xlsx"], "HIV+ (CCC)" );
sc00.cum.start = strd2015_numcum_results(sc00.cc_incid_art(:, 70:end, :), pop2015,  ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S00.xlsx"], "HIV+ ART (CCC)" );
sc00.cum.stnart = strd2015_numcum_results(sc00.cc_incid_pos(:, 70:end, :), pop2015,  ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S00.xlsx"], "HIV+ no ART (CCC)" );


%[sc00.vacc, sc00.screen] = coverageest(sc00.res);
%% baseline 2 w/ bivalent - Sc0b

[sc0b.cc_incid_all, sc0b.cc_incid_neg, sc0b.cc_incid_posall, sc0b.cc_incid_art, sc0b.cc_incid_pos] = ccincid(sc0b.res); % pp

sc0b.stall = strd2015_results(sc0b.cc_incid_all, pop2015, 2);
sc0b.stneg = strd2015_results(sc0b.cc_incid_neg, pop2015, 2);
sc0b.stpos = strd2015_results(sc0b.cc_incid_posall, pop2015, 37);
sc0b.start = strd2015_results(sc0b.cc_incid_art, pop2015, 56);
sc0b.stnart = strd2015_results(sc0b.cc_incid_pos, pop2015, 37);

sc0b.tab_ccall = res_table(sc0b.stall(:, 69:end), sc0b.cc_incid_all(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S0b.xlsx"], "Pop(All) ICC");            
sc0b.tab_ccneg = res_table(sc0b.stneg(:, 69:end), sc0b.cc_incid_neg(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S0b.xlsx"], "HIV- (ICC)");
sc0b.tab_ccpos = res_table(sc0b.stpos(:, 34:end), sc0b.cc_incid_posall(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S0b.xlsx"], "HIV+ (ICC)");
sc0b.tab_ccart = res_table(sc0b.start(:, 15:end), sc0b.cc_incid_art(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S0b.xlsx"], "HIV+ ART (ICC)");
sc0b.tab_ccnart = res_table(sc0b.stnart(:, 34:end), sc0b.cc_incid_pos(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S0b.xlsx"], "HIV+ no ART (ICC)");

sc0b.tab_pops = res_poptable(sc0b.res, ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S0b.xlsx"], 102);

sc0b.tab_ccc = res_ccctable(sc0b.res, ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S0b.xlsx"]);

sc0b.cum.stall = strd2015_numcum_results(sc0b.cc_incid_all(:, 70:end, :), pop2015,  ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S0b.xlsx"], "Pop(All) CCC");
sc0b.cum.stneg = strd2015_numcum_results(sc0b.cc_incid_neg(:, 70:end, :), pop2015, ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S0b.xlsx"], "HIV- (CCC)");
sc0b.cum.stpos = strd2015_numcum_results(sc0b.cc_incid_posall(:, 70:end, :), pop2015, ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S0b.xlsx"], "HIV+ (CCC)" );
sc0b.cum.start = strd2015_numcum_results(sc0b.cc_incid_art(:, 70:end, :), pop2015, ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S0b.xlsx"], "HIV+ ART (CCC)" );
sc0b.cum.stnart = strd2015_numcum_results(sc0b.cc_incid_pos(:, 70:end, :), pop2015, ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S0b.xlsx"], "HIV+ no ART (CCC)" );


%[sc0b.vacc, sc0b.vaccpos, sc0b.screen] = coverageest(sc0b.res);

% %% baseline 2 w/ bivalent - Sc0b
% [sc0b2.cc_incid_all, sc0b2.cc_incid_neg, sc0b2.cc_incid_posall, sc0b2.cc_incid_art, sc0b2.cc_incid_pos] = ccincid(sc0b2.res); % pp
% 
% sc0b2.stall = strd2015_results(sc0b2.cc_incid_all, pop2015, 2);
% sc0b2.stneg = strd2015_results(sc0b2.cc_incid_neg, pop2015, 2);
% sc0b2.stpos = strd2015_results(sc0b2.cc_incid_posall, pop2015, 37);
% sc0b2.start = strd2015_results(sc0b2.cc_incid_art, pop2015, 56);
% sc0b2.stnart = strd2015_results(sc0b2.cc_incid_pos, pop2015, 37);
% 
% [sc0b2.vacc, sc0b2.screen] = coverageest(sc0b2.res);

%% baseline 2 w/ bivalent - Sc0c


% for i=2:22
%     
% sc0c.res2{i-1} = sc0c.res{i};
% 
% end

[sc0c.cc_incid_all, sc0c.cc_incid_neg, sc0c.cc_incid_posall, sc0c.cc_incid_art, sc0c.cc_incid_pos] = ccincid(sc0c.res); % pp

sc0c.stall = strd2015_results(sc0c.cc_incid_all, pop2015, 2);
sc0c.stneg = strd2015_results(sc0c.cc_incid_neg, pop2015, 2);
sc0c.stpos = strd2015_results(sc0c.cc_incid_posall, pop2015, 37);
sc0c.start = strd2015_results(sc0c.cc_incid_art, pop2015, 56);
sc0c.stnart = strd2015_results(sc0c.cc_incid_pos, pop2015, 37);

sc0c.tab_ccall = res_table(sc0c.stall(:, 69:end), sc0c.cc_incid_all(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S0c.xlsx"], "Pop(All) ICC");            
sc0c.tab_ccneg = res_table(sc0c.stneg(:, 69:end), sc0c.cc_incid_neg(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S0c.xlsx"], "HIV- (ICC)");
sc0c.tab_ccpos = res_table(sc0c.stpos(:, 34:end), sc0c.cc_incid_posall(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S0c.xlsx"], "HIV+ (ICC)");
sc0c.tab_ccart = res_table(sc0c.start(:, 15:end), sc0c.cc_incid_art(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S0c.xlsx"], "HIV+ ART (ICC)");
sc0c.tab_ccnart = res_table(sc0c.stnart(:, 34:end), sc0c.cc_incid_pos(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S0c.xlsx"], "HIV+ no ART (ICC)");

sc0c.tab_pops = res_poptable(sc0c.res, ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S0c.xlsx"], 102);

sc0c.tab_ccc = res_ccctable(sc0c.res, ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S0c.xlsx"]);

sc0c.cum.stall = strd2015_numcum_results(sc0c.cc_incid_all(:, 70:end, :), pop2015, ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S0c.xlsx"], "Pop(All) CCC");
sc0c.cum.stneg = strd2015_numcum_results(sc0c.cc_incid_neg(:, 70:end, :), pop2015, ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S0c.xlsx"], "HIV- (CCC)");
sc0c.cum.stpos = strd2015_numcum_results(sc0c.cc_incid_posall(:, 70:end, :), pop2015, ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S0c.xlsx"], "HIV+ (CCC)" );
sc0c.cum.start = strd2015_numcum_results(sc0c.cc_incid_art(:, 70:end, :), pop2015, ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S0c.xlsx"], "HIV+ ART (CCC)" );
sc0c.cum.stnart = strd2015_numcum_results(sc0c.cc_incid_pos(:, 70:end, :), pop2015, ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S0c.xlsx"], "HIV+ no ART (CCC)" );


%[sc0c.vacc, sc0c.vaccpos, sc0c.screen] = coverageest(sc0c.res);

%% Sc1
[sc1.cc_incid_all, sc1.cc_incid_neg, sc1.cc_incid_posall, sc1.cc_incid_art, sc1.cc_incid_pos] = ccincid(sc1.res); % pp

sc1.stall = strd2015_results(sc1.cc_incid_all, pop2015, 2);
sc1.stneg = strd2015_results(sc1.cc_incid_neg, pop2015, 2);
sc1.stpos = strd2015_results(sc1.cc_incid_posall, pop2015, 37);
sc1.start = strd2015_results(sc1.cc_incid_art, pop2015, 56);
sc1.stnart = strd2015_results(sc1.cc_incid_pos, pop2015, 37);

sc1.tab_ccall = res_table(sc1.stall(:, 69:end), sc1.cc_incid_all(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S1.xlsx"], "Pop(All) ICC");            
sc1.tab_ccneg = res_table(sc1.stneg(:, 69:end), sc1.cc_incid_neg(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S1.xlsx"], "HIV- (ICC)");
sc1.tab_ccpos = res_table(sc1.stpos(:, 34:end), sc1.cc_incid_posall(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S1.xlsx"], "HIV+ (ICC)");
sc1.tab_ccart = res_table(sc1.start(:, 15:end), sc1.cc_incid_art(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S1.xlsx"], "HIV+ ART (ICC)");
sc1.tab_ccnart = res_table(sc1.stnart(:, 34:end), sc1.cc_incid_pos(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S1.xlsx"], "HIV+ no ART (ICC)");

sc1.tab_pops = res_poptable(sc1.res, ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S1.xlsx"], 102);

sc1.tab_ccc = res_ccctable(sc1.res, ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S1.xlsx"]);

sc1.cum.stall = strd2015_numcum_results(sc1.cc_incid_all(:, 70:end, :), pop2015,  ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S1.xlsx"], "Pop(All) CCC");
sc1.cum.stneg = strd2015_numcum_results(sc1.cc_incid_neg(:, 70:end, :), pop2015,  ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S1.xlsx"], "HIV- (CCC)");
sc1.cum.stpos = strd2015_numcum_results(sc1.cc_incid_posall(:, 70:end, :), pop2015,  ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S1.xlsx"], "HIV+ (CCC)" );
sc1.cum.start = strd2015_numcum_results(sc1.cc_incid_art(:, 70:end, :), pop2015,  ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S1.xlsx"], "HIV+ ART (CCC)" );
sc1.cum.stnart = strd2015_numcum_results(sc1.cc_incid_pos(:, 70:end, :), pop2015,  ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S1.xlsx"], "HIV+ no ART (CCC)" );

%[sc1.vacc, sc1.screen] = coverageest(sc1.res);

%%%%%%%

for i=1:length(sc0.stneg(:,1))
    for t=1:170
         sc1.prreduc(i,t) = 1 - sc1.stall(i,t) / sc0.stall(i,t);
    end
end

sc1.summ(:,1) = 2019:2120;
sc1.summ(:,2) = prctile(sc1.prreduc(:,69:end), 50);
sc1.summ(:,3) = prctile(sc1.prreduc(:,69:end), 5);
sc1.summ(:,4) = prctile(sc1.prreduc(:,69:end), 95);

tabnew2 = array2table(sc1.summ);
writetable(tabnew2, "/Users/minttu/who_results_sheets/ICL_HSPH_PercentReduct_S1.xlsx", 'WriteVariableNames',0, 'Sheet', "All")


for i=1:length(sc0.stpos(:,1))
    for t=1:135
         sc1.prreduc2(i,t) = 1 - sc1.stpos(i,t) / sc0.stpos(i,t);
    end
end

sc1.summ2(:,1) = 2019:2120;
sc1.summ2(:,2) = prctile(sc1.prreduc2(:,34:end), 50);
sc1.summ2(:,3) = prctile(sc1.prreduc2(:,34:end), 5);
sc1.summ2(:,4) = prctile(sc1.prreduc2(:,34:end), 95);

tabnew3 = array2table(sc1.summ2);
writetable(tabnew3, "/Users/minttu/who_results_sheets/ICL_HSPH_PercentReduct_S1.xlsx", 'WriteVariableNames',0, 'Sheet', "HIV+")



%% Sc2
[sc2.cc_incid_all, sc2.cc_incid_neg, sc2.cc_incid_posall, sc2.cc_incid_art, sc2.cc_incid_pos] = ccincid(sc2.res); % pp

sc2.stall = strd2015_results(sc2.cc_incid_all, pop2015, 2);
sc2.stneg = strd2015_results(sc2.cc_incid_neg, pop2015, 2);
sc2.stpos = strd2015_results(sc2.cc_incid_posall, pop2015, 37);
sc2.start = strd2015_results(sc2.cc_incid_art, pop2015, 56);
sc2.stnart = strd2015_results(sc2.cc_incid_pos, pop2015, 37);

sc2.tab_ccall = res_table(sc2.stall(:, 69:end), sc2.cc_incid_all(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S2.xlsx"], "Pop(All) ICC");            
sc2.tab_ccneg = res_table(sc2.stneg(:, 69:end), sc2.cc_incid_neg(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S2.xlsx"], "HIV- (ICC)");
sc2.tab_ccpos = res_table(sc2.stpos(:, 34:end), sc2.cc_incid_posall(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S2.xlsx"], "HIV+ (ICC)");
sc2.tab_ccart = res_table(sc2.start(:, 15:end), sc2.cc_incid_art(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S2.xlsx"], "HIV+ ART (ICC)");
sc2.tab_ccnart = res_table(sc2.stnart(:, 34:end), sc2.cc_incid_pos(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S2.xlsx"], "HIV+ no ART (ICC)");

sc2.tab_pops = res_poptable(sc2.res, ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S2.xlsx"], 102);

sc2.tab_ccc = res_ccctable(sc2.res, ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S2.xlsx"]);

sc2.cum.stall = strd2015_numcum_results(sc2.cc_incid_all(:, 70:end, :), pop2015, ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S2.xlsx"], "Pop(All) CCC");
sc2.cum.stneg = strd2015_numcum_results(sc2.cc_incid_neg(:, 70:end, :), pop2015, ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S2.xlsx"], "HIV- (CCC)");
sc2.cum.stpos = strd2015_numcum_results(sc2.cc_incid_posall(:, 70:end, :), pop2015,  ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S2.xlsx"], "HIV+ (CCC)" );
sc2.cum.start = strd2015_numcum_results(sc2.cc_incid_art(:, 70:end, :), pop2015,  ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S2.xlsx"], "HIV+ ART (CCC)" );
sc2.cum.stnart = strd2015_numcum_results(sc2.cc_incid_pos(:, 70:end, :), pop2015,  ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S2.xlsx"], "HIV+ no ART (CCC)" );

sc2.tab_coverage = res_covtable(sc2.res, ["/Users/minttu/who_results_sheets/ICL_HSPH_Coverage-(2020-2120)_S2.xlsx"], 101);

% 5% LUB

sc2.tab_ccall_p5 = res_table_prctile(sc2.stall(:, 70:end), sc2.cc_incid_all(:, 71:end, :), 5, ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)LUB_S2.xlsx"], "Pop(All) ICC");            
sc2.tab_ccneg_p5 = res_table_prctile(sc2.stneg(:, 70:end), sc2.cc_incid_neg(:, 71:end, :), 5,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)LUB_S2.xlsx"], "HIV- (ICC)");
sc2.tab_ccpos_p5 = res_table_prctile(sc2.stpos(:, 35:end), sc2.cc_incid_posall(:, 71:end, :), 5,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)LUB_S2.xlsx"], "HIV+ (ICC)");
sc2.tab_ccart_p5 = res_table_prctile(sc2.start(:, 16:end), sc2.cc_incid_art(:, 71:end, :), 5,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)LUB_S2.xlsx"], "HIV+ ART (ICC)");
sc2.tab_ccnart_p5 = res_table_prctile(sc2.stnart(:, 35:end), sc2.cc_incid_pos(:, 71:end, :), 5,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)LUB_S2.xlsx"], "HIV+ no ART (ICC)");

sc2.tab_pops_p5 = res_poptable_prctile(sc2.res, 5, ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)LUB_S2.xlsx"], 101);


% 95% LUB

sc2.tab_ccall_p95 = res_table_prctile(sc2.stall(:, 70:end), sc2.cc_incid_all(:, 71:end, :), 95, ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)HUB_S2.xlsx"], "Pop(All) ICC");            
sc2.tab_ccneg_p95 = res_table_prctile(sc2.stneg(:, 70:end), sc2.cc_incid_neg(:, 71:end, :), 95,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)HUB_S2.xlsx"], "HIV- (ICC)");
sc2.tab_ccpos_p95 = res_table_prctile(sc2.stpos(:, 35:end), sc2.cc_incid_posall(:, 71:end, :), 95,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)HUB_S2.xlsx"], "HIV+ (ICC)");
sc2.tab_ccart_p95 = res_table_prctile(sc2.start(:, 16:end), sc2.cc_incid_art(:, 71:end, :), 95,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)HUB_S2.xlsx"], "HIV+ ART (ICC)");
sc2.tab_ccnart_p95 = res_table_prctile(sc2.stnart(:, 35:end), sc2.cc_incid_pos(:, 71:end, :), 95,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)HUB_S2.xlsx"], "HIV+ no ART (ICC)");

sc2.tab_pops_p95 = res_poptable_prctile(sc2.res, 95, ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)HUB_S2.xlsx"], 101);

% hr HPV women

[sc2.hrhpvf_all, sc2.hrhpvf_neg, sc2.hrhpvf_posall, sc2.hrhpvf_art, sc2.hrhpvf_pos, ...
 sc2.hrhpvf_all_num, sc2.hrhpvf_neg_num, sc2.hrhpvf_posall_num, sc2.hrhpvf_art_num, sc2.hrhpvf_pos_num] = hrhpvprev(sc2.res); % pp

sc2.stall_hrhpvf = strd2015_results_num(sc2.hrhpvf_all, pop2015, 2);
sc2.stneg_hrhpvf = strd2015_results_num(sc2.hrhpvf_neg, pop2015, 2);
sc2.stpos_hrhpvf = strd2015_results_num(sc2.hrhpvf_posall, pop2015, 37);
sc2.start_hrhpvf = strd2015_results_num(sc2.hrhpvf_art, pop2015, 56);
sc2.stnart_hrhpvf = strd2015_results_num(sc2.hrhpvf_pos, pop2015, 37);

% here weird spacing
sc2.tab_all_hrhpvf = res_table_prctile_num(sc2.stall_hrhpvf(:, 70:end), sc2.hrhpvf_all_num(:, 71:end, :), 50, ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_hrHPV_Prevalence-standardised-(2020-2120)_S2.xlsx"], "Pop(All) (hrHPV)");            
sc2.tab_neg_hrhpvf = res_table_prctile_num(sc2.stneg_hrhpvf(:, 70:end), sc2.hrhpvf_neg_num(:, 71:end, :), 50,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_hrHPV_Prevalence-standardised-(2020-2120)_S2.xlsx"], "HIV-  (hrHPV)");
sc2.tab_pos_hrhpvf = res_table_prctile_num(sc2.stpos_hrhpvf(:, 35:end), sc2.hrhpvf_posall_num(:, 71:end, :), 50,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_hrHPV_Prevalence-standardised-(2020-2120)_S2.xlsx"], "HIV+   (hrHPV)"); % this has 3 spaces
sc2.tab_art_hrhpvf = res_table_prctile_num(sc2.start_hrhpvf(:, 16:end), sc2.hrhpvf_art_num(:, 71:end, :), 50,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_hrHPV_Prevalence-standardised-(2020-2120)_S2.xlsx"], "HIV+ ART  (hrHPV)");
sc2.tab_nart_hrhpvf = res_table_prctile_num(sc2.stnart_hrhpvf(:, 35:end), sc2.hrhpvf_pos_num(:, 71:end, :), 50,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_hrHPV_Prevalence-standardised-(2020-2120)_S2.xlsx"], "HIV+ no ART  (hrHPV)");

sc2.tab_pops_hrhpvf = res_poptable_prctile(sc2.res, 50, ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_hrHPV_Prevalence-standardised-(2020-2120)_S2.xlsx"], 101);


% VT-hpv women

[sc2.hpv9vtf_all, sc2.hpv9vtf_neg, sc2.hpv9vtf_posall, sc2.hpv9vtf_art, sc2.hpv9vtf_pos, ...
 sc2.hpv9vtf_all_num, sc2.hpv9vtf_neg_num, sc2.hpv9vtf_posall_num, sc2.hpv9vtf_art_num, sc2.hpv9vtf_pos_num] = vthpvprev(sc2.res); % pp

sc2.stall_hpv9vtf = strd2015_results_num(sc2.hpv9vtf_all, pop2015, 2);
sc2.stneg_hpv9vtf = strd2015_results_num(sc2.hpv9vtf_neg, pop2015, 2);
sc2.stpos_hpv9vtf = strd2015_results_num(sc2.hpv9vtf_posall, pop2015, 37);
sc2.start_hpv9vtf = strd2015_results_num(sc2.hpv9vtf_art, pop2015, 56);
sc2.stnart_hpv9vtf = strd2015_results_num(sc2.hpv9vtf_pos, pop2015, 37);

% here weird spacing
sc2.tab_all_hpv9vtf = res_table_prctile_num(sc2.stall_hpv9vtf(:, 70:end), sc2.hpv9vtf_all_num(:, 71:end, :), 50, ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_VThrHPV_Prevalence-standardised-(2020-2120)_S2.xlsx"], "Pop(All) (vtHPV)");            
sc2.tab_neg_hpv9vtf = res_table_prctile_num(sc2.stneg_hpv9vtf(:, 70:end), sc2.hpv9vtf_neg_num(:, 71:end, :), 50,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_VThrHPV_Prevalence-standardised-(2020-2120)_S2.xlsx"], "HIV-  (vtHPV)");
sc2.tab_pos_hpv9vtf = res_table_prctile_num(sc2.stpos_hpv9vtf(:, 35:end), sc2.hpv9vtf_posall_num(:, 71:end, :), 50,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_VThrHPV_Prevalence-standardised-(2020-2120)_S2.xlsx"], "HIV+   (vtHPV)"); % this has 3 spaces
sc2.tab_art_hpv9vtf = res_table_prctile_num(sc2.start_hpv9vtf(:, 16:end), sc2.hpv9vtf_art_num(:, 71:end, :), 50,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_VThrHPV_Prevalence-standardised-(2020-2120)_S2.xlsx"], "HIV+ ART  (vtHPV)");
sc2.tab_nart_hpv9vtf = res_table_prctile_num(sc2.stnart_hpv9vtf(:, 35:end), sc2.hpv9vtf_pos_num(:, 71:end, :), 50,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_VThrHPV_Prevalence-standardised-(2020-2120)_S2.xlsx"], "HIV+ no ART  (vtHPV)");

sc2.tab_pops_hpv9vtf = res_poptable_prctile(sc2.res, 50, ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_VThrHPV_Prevalence-standardised-(2020-2120)_S2.xlsx"], 101);


% CIN2+ prevalence 

[sc2.cinprf_all, sc2.cinprf_neg, sc2.cinprf_posall, sc2.cinprf_art, sc2.cinprf_pos, ...
 sc2.cinprf_all_num, sc2.cinprf_neg_num, sc2.cinprf_posall_num, sc2.cinprf_art_num, sc2.cinprf_pos_num] = cin2prev(sc2.res); % pp

sc2.stall_cinprf = strd2015_results_num(sc2.cinprf_all, pop2015, 2);
sc2.stneg_cinprf = strd2015_results_num(sc2.cinprf_neg, pop2015, 2);
sc2.stpos_cinprf = strd2015_results_num(sc2.cinprf_posall, pop2015, 37);
sc2.start_cinprf = strd2015_results_num(sc2.cinprf_art, pop2015, 56);
sc2.stnart_cinprf = strd2015_results_num(sc2.cinprf_pos, pop2015, 37);

% here weird spacing
sc2.tab_all_cinprf = res_table_prctile_num(sc2.stall_cinprf(:, 70:end), sc2.cinprf_all_num(:, 71:end, :), 50, ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CIN2+_Prevalence-standardised-(2020-2120)_S2.xlsx"], "Pop(All) (CIN2+)");            
sc2.tab_neg_cinprf = res_table_prctile_num(sc2.stneg_cinprf(:, 70:end), sc2.cinprf_neg_num(:, 71:end, :), 50,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CIN2+_Prevalence-standardised-(2020-2120)_S2.xlsx"], "HIV- (CIN2+)"); % 1 space only
sc2.tab_pos_cinprf = res_table_prctile_num(sc2.stpos_cinprf(:, 35:end), sc2.cinprf_posall_num(:, 71:end, :), 50,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CIN2+_Prevalence-standardised-(2020-2120)_S2.xlsx"], "HIV+  (CIN2+)"); % this has 2 spaces
sc2.tab_art_cinprf = res_table_prctile_num(sc2.start_cinprf(:, 16:end), sc2.cinprf_art_num(:, 71:end, :), 50,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CIN2+_Prevalence-standardised-(2020-2120)_S2.xlsx"], "HIV+ ART  (CIN2+)");
sc2.tab_nart_cinprf = res_table_prctile_num(sc2.stnart_cinprf(:, 35:end), sc2.cinprf_pos_num(:, 71:end, :), 50,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CIN2+_Prevalence-standardised-(2020-2120)_S2.xlsx"], "HIV+ no ART  (CIN2+)");

sc2.tab_pops_cinprf = res_poptable_prctile(sc2.res, 50, ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CIN2+_Prevalence-standardised-(2020-2120)_S2.xlsx"], 101);

%%%%%%%

for i=1:length(sc0.stneg(:,1))
    for t=1:170
         sc2.prreduc(i,t) = 1 - sc2.stall(i,t) / sc0.stall(i,t);
    end
end

% alternative way takes reduction of the median
% sc2.summ2(:,1) = 2019:2120;
% sc2.summ2(:,2) = 1 - median(sc2.stall(:,69:end)) ./ median(sc0.stall(:,69:end));
% sc2.summ2(:,3) = 1 - prctile(sc2.stall(:,69:end), 5) ./ prctile(sc0.stall(:,69:end), 5);
% sc2.summ2(:,4) = 1 - prctile(sc2.stall(:,69:end), 95) ./ prctile(sc0.stall(:,69:end), 95);

sc2.summ(:,1) = 2019:2120;
sc2.summ(:,2) = prctile(sc2.prreduc(:,69:end), 50);
sc2.summ(:,3) = prctile(sc2.prreduc(:,69:end), 5);
sc2.summ(:,4) = prctile(sc2.prreduc(:,69:end), 95);

tabnew2 = array2table(sc2.summ);
writetable(tabnew2, "/Users/minttu/who_results_sheets/ICL_HSPH_PercentReduct_S2.xlsx", 'WriteVariableNames',0, 'Sheet', "All")


for i=1:length(sc0.stpos(:,1))
    for t=1:135
         sc2.prreduc2(i,t) = 1 - sc2.stpos(i,t) / sc0.stpos(i,t);
    end
end

sc2.summ2(:,1) = 2019:2120;
sc2.summ2(:,2) = prctile(sc2.prreduc2(:,34:end), 50);
sc2.summ2(:,3) = prctile(sc2.prreduc2(:,34:end), 5);
sc2.summ2(:,4) = prctile(sc2.prreduc2(:,34:end), 95);

tabnew3 = array2table(sc2.summ2);
writetable(tabnew3, "/Users/minttu/who_results_sheets/ICL_HSPH_PercentReduct_S2.xlsx", 'WriteVariableNames',0, 'Sheet', "HIV+")



%% Sc3
[sc3.cc_incid_all, sc3.cc_incid_neg, sc3.cc_incid_posall, sc3.cc_incid_art, sc3.cc_incid_pos] = ccincid(sc3.res); % pp

sc3.stall = strd2015_results(sc3.cc_incid_all, pop2015, 2);
sc3.stneg = strd2015_results(sc3.cc_incid_neg, pop2015, 2);
sc3.stpos = strd2015_results(sc3.cc_incid_posall, pop2015, 37);
sc3.start = strd2015_results(sc3.cc_incid_art, pop2015, 56);
sc3.stnart = strd2015_results(sc3.cc_incid_pos, pop2015, 37);

% sc3.tab_ccall = res_table(sc3.stall(:, 70:end), sc3.cc_incid_all(:, 71:end, :));
% sc3.tab_ccneg = res_table(sc3.stneg(:, 70:end), sc3.cc_incid_neg(:, 71:end, :));
% sc3.tab_ccpos = res_table(sc3.stpos(:, 35:end), sc3.cc_incid_posall(:, 71:end, :));
% sc3.tab_ccart = res_table(sc3.start(:, 16:end), sc3.cc_incid_art(:, 71:end, :));
% sc3.tab_ccnart = res_table(sc3.stnart(:, 35:end),sc3.cc_incid_pos(:, 71:end, :));

sc3.tab_ccall = res_table(sc3.stall(:, 69:end), sc3.cc_incid_all(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S3.xlsx"], "Pop(All) ICC");            
sc3.tab_ccneg = res_table(sc3.stneg(:, 69:end), sc3.cc_incid_neg(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S3.xlsx"], "HIV- (ICC)");
sc3.tab_ccpos = res_table(sc3.stpos(:, 34:end), sc3.cc_incid_posall(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S3.xlsx"], "HIV+ (ICC)");
sc3.tab_ccart = res_table(sc3.start(:, 15:end), sc3.cc_incid_art(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S3.xlsx"], "HIV+ ART (ICC)");
sc3.tab_ccnart = res_table(sc3.stnart(:, 34:end), sc3.cc_incid_pos(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S3.xlsx"], "HIV+ no ART (ICC)");

sc3.tab_pops = res_poptable(sc3.res, ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S3.xlsx"], 102);

sc3.tab_ccc = res_ccctable(sc3.res, ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S3.xlsx"]);

sc3.cum.stall = strd2015_numcum_results(sc3.cc_incid_all(:, 70:end, :), pop2015, ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S3.xlsx"], "Pop(All) CCC");
sc3.cum.stneg = strd2015_numcum_results(sc3.cc_incid_neg(:, 70:end, :), pop2015, ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S3.xlsx"], "HIV- (CCC)");
sc3.cum.stpos = strd2015_numcum_results(sc3.cc_incid_posall(:, 70:end, :), pop2015, ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S3.xlsx"], "HIV+ (CCC)" );
sc3.cum.start = strd2015_numcum_results(sc3.cc_incid_art(:, 70:end, :), pop2015, ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S3.xlsx"], "HIV+ ART (CCC)" );
sc3.cum.stnart = strd2015_numcum_results(sc3.cc_incid_pos(:, 70:end, :), pop2015, ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S3.xlsx"], "HIV+ no ART (CCC)" );

%%%%%%

for i=1:length(sc3.stall(:,1))
    for t=1:170
         sc3.prreduc(i,t) = 1 - sc3.stall(i,t) / sc0.stall(i,t);
    end
end

sc3.summ(:,1) = 2019:2120;
sc3.summ(:,2) = prctile(sc3.prreduc(:,69:end), 50);
sc3.summ(:,3) = prctile(sc3.prreduc(:,69:end), 5);
sc3.summ(:,4) = prctile(sc3.prreduc(:,69:end), 95);

tabnew2 = array2table(sc3.summ);
writetable(tabnew2, "/Users/minttu/who_results_sheets/ICL_HSPH_PercentReduct_S3.xlsx", 'WriteVariableNames',0, 'Sheet', "All")

for i=1:length(sc0.stpos(:,1))
    for t=1:135
         sc3.prreduc2(i,t) = 1 - sc3.stpos(i,t) / sc0.stpos(i,t);
    end
end

sc3.summ2(:,1) = 2019:2120;
sc3.summ2(:,2) = prctile(sc3.prreduc2(:,34:end), 50);
sc3.summ2(:,3) = prctile(sc3.prreduc2(:,34:end), 5);
sc3.summ2(:,4) = prctile(sc3.prreduc2(:,34:end), 95);

tabnew3 = array2table(sc3.summ2);
writetable(tabnew3, "/Users/minttu/who_results_sheets/ICL_HSPH_PercentReduct_S3.xlsx", 'WriteVariableNames',0, 'Sheet', "HIV+")



%% Sc4
[sc4.cc_incid_all, sc4.cc_incid_neg, sc4.cc_incid_posall, sc4.cc_incid_art, sc4.cc_incid_pos] = ccincid(sc4.res); % pp

sc4.stall = strd2015_results(sc4.cc_incid_all, pop2015, 2);
sc4.stneg = strd2015_results(sc4.cc_incid_neg, pop2015, 2);
sc4.stpos = strd2015_results(sc4.cc_incid_posall, pop2015, 37);
sc4.start = strd2015_results(sc4.cc_incid_art, pop2015, 56);
sc4.stnart = strd2015_results(sc4.cc_incid_pos, pop2015, 37);


sc4.tab_ccall = res_table(sc4.stall(:, 69:end), sc4.cc_incid_all(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S4.xlsx"], "Pop(All) ICC");            
sc4.tab_ccneg = res_table(sc4.stneg(:, 69:end), sc4.cc_incid_neg(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S4.xlsx"], "HIV- (ICC)");
sc4.tab_ccpos = res_table(sc4.stpos(:, 34:end), sc4.cc_incid_posall(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S4.xlsx"], "HIV+ (ICC)");
sc4.tab_ccart = res_table(sc4.start(:, 15:end), sc4.cc_incid_art(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S4.xlsx"], "HIV+ ART (ICC)");
sc4.tab_ccnart = res_table(sc4.stnart(:, 34:end), sc4.cc_incid_pos(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S4.xlsx"], "HIV+ no ART (ICC)");

sc4.tab_pops = res_poptable(sc4.res, ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S4.xlsx"], 102);

sc4.tab_ccc = res_ccctable(sc4.res, ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S4.xlsx"]);

sc4.cum.stall = strd2015_numcum_results(sc4.cc_incid_all(:, 70:end, :), pop2015, ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S4.xlsx"], "Pop(All) CCC");
sc4.cum.stneg = strd2015_numcum_results(sc4.cc_incid_neg(:, 70:end, :), pop2015, ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S4.xlsx"], "HIV- (CCC)");
sc4.cum.stpos = strd2015_numcum_results(sc4.cc_incid_posall(:, 70:end, :), pop2015, ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S4.xlsx"], "HIV+ (CCC)" );
sc4.cum.start = strd2015_numcum_results(sc4.cc_incid_art(:, 70:end, :), pop2015, ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S4.xlsx"], "HIV+ ART (CCC)" );
sc4.cum.stnart = strd2015_numcum_results(sc4.cc_incid_pos(:, 70:end, :), pop2015, ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S4.xlsx"], "HIV+ no ART (CCC)" );

sc4.tab_coverage = res_covtable(sc4.res, ["/Users/minttu/who_results_sheets/ICL_HSPH_Coverage-(2020-2120)_S4.xlsx"], 101);

% 5% LUB

sc4.tab_ccall_p5 = res_table_prctile(sc4.stall(:, 70:end), sc4.cc_incid_all(:, 71:end, :), 5, ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)LUB_S4.xlsx"], "Pop(All) ICC");            
sc4.tab_ccneg_p5 = res_table_prctile(sc4.stneg(:, 70:end), sc4.cc_incid_neg(:, 71:end, :), 5,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)LUB_S4.xlsx"], "HIV- (ICC)");
sc4.tab_ccpos_p5 = res_table_prctile(sc4.stpos(:, 35:end), sc4.cc_incid_posall(:, 71:end, :), 5,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)LUB_S4.xlsx"], "HIV+ (ICC)");
sc4.tab_ccart_p5 = res_table_prctile(sc4.start(:, 16:end), sc4.cc_incid_art(:, 71:end, :), 5,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)LUB_S4.xlsx"], "HIV+ ART (ICC)");
sc4.tab_ccnart_p5 = res_table_prctile(sc4.stnart(:, 35:end), sc4.cc_incid_pos(:, 71:end, :), 5,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)LUB_S4.xlsx"], "HIV+ no ART (ICC)");

sc4.tab_pops_p5 = res_poptable_prctile(sc4.res, 5, ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)LUB_S4.xlsx"], 101);


% 95% LUB

sc4.tab_ccall_p95 = res_table_prctile(sc4.stall(:, 70:end), sc4.cc_incid_all(:, 71:end, :), 95, ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)HUB_S4.xlsx"], "Pop(All) ICC");            
sc4.tab_ccneg_p95 = res_table_prctile(sc4.stneg(:, 70:end), sc4.cc_incid_neg(:, 71:end, :), 95,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)HUB_S4.xlsx"], "HIV- (ICC)");
sc4.tab_ccpos_p95 = res_table_prctile(sc4.stpos(:, 35:end), sc4.cc_incid_posall(:, 71:end, :), 95,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)HUB_S4.xlsx"], "HIV+ (ICC)");
sc4.tab_ccart_p95 = res_table_prctile(sc4.start(:, 16:end), sc4.cc_incid_art(:, 71:end, :), 95,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)HUB_S4.xlsx"], "HIV+ ART (ICC)");
sc4.tab_ccnart_p95 = res_table_prctile(sc4.stnart(:, 35:end), sc4.cc_incid_pos(:, 71:end, :), 95,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)HUB_S4.xlsx"], "HIV+ no ART (ICC)");

sc4.tab_pops_p95 = res_poptable_prctile(sc4.res, 95, ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)HUB_S4.xlsx"], 101);

%%%%%

for i=1:length(sc4.stall(:,1))
    for t=1:170
         sc4.prreduc(i,t) = 1 - sc4.stall(i,t) / sc0.stall(i,t);
    end
end

sc4.summ(:,1) = 2019:2120;
sc4.summ(:,2) = prctile(sc4.prreduc(:,69:end), 50);
sc4.summ(:,3) = prctile(sc4.prreduc(:,69:end), 5);
sc4.summ(:,4) = prctile(sc4.prreduc(:,69:end), 95);

tabnew2 = array2table(sc4.summ);
writetable(tabnew2, "/Users/minttu/who_results_sheets/ICL_HSPH_PercentReduct_S4.xlsx", 'WriteVariableNames',0, 'Sheet', "All")

for i=1:length(sc0.stpos(:,1))
    for t=1:135
         sc4.prreduc2(i,t) = 1 - sc4.stpos(i,t) / sc0.stpos(i,t);
    end
end

sc4.summ2(:,1) = 2019:2120;
sc4.summ2(:,2) = prctile(sc4.prreduc2(:,34:end), 50);
sc4.summ2(:,3) = prctile(sc4.prreduc2(:,34:end), 5);
sc4.summ2(:,4) = prctile(sc4.prreduc2(:,34:end), 95);

tabnew3 = array2table(sc4.summ2);
writetable(tabnew3, "/Users/minttu/who_results_sheets/ICL_HSPH_PercentReduct_S4.xlsx", 'WriteVariableNames',0, 'Sheet', "HIV+")


%% Sc5
[sc5.cc_incid_all, sc5.cc_incid_neg, sc5.cc_incid_posall, sc5.cc_incid_art, sc5.cc_incid_pos] = ccincid(sc5.res); % pp

sc5.stall = strd2015_results(sc5.cc_incid_all, pop2015, 2);
sc5.stneg = strd2015_results(sc5.cc_incid_neg, pop2015, 2);
sc5.stpos = strd2015_results(sc5.cc_incid_posall, pop2015, 37);
sc5.start = strd2015_results(sc5.cc_incid_art, pop2015, 56);
sc5.stnart = strd2015_results(sc5.cc_incid_pos, pop2015, 37);

sc5.tab_ccall = res_table(sc5.stall(:, 69:end), sc5.cc_incid_all(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S5.xlsx"], "Pop(All) ICC");            
sc5.tab_ccneg = res_table(sc5.stneg(:, 69:end), sc5.cc_incid_neg(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S5.xlsx"], "HIV- (ICC)");
sc5.tab_ccpos = res_table(sc5.stpos(:, 34:end), sc5.cc_incid_posall(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S5.xlsx"], "HIV+ (ICC)");
sc5.tab_ccart = res_table(sc5.start(:, 15:end), sc5.cc_incid_art(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S5.xlsx"], "HIV+ ART (ICC)");
sc5.tab_ccnart = res_table(sc5.stnart(:, 34:end), sc5.cc_incid_pos(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S5.xlsx"], "HIV+ no ART (ICC)");

sc5.tab_pops = res_poptable(sc5.res, ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S5.xlsx"], 102);

sc5.tab_ccc = res_ccctable(sc5.res, ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S5.xlsx"]);

sc5.cum.stall = strd2015_numcum_results(sc5.cc_incid_all(:, 70:end, :), pop2015, ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S5.xlsx"], "Pop(All) CCC");
sc5.cum.stneg = strd2015_numcum_results(sc5.cc_incid_neg(:, 70:end, :), pop2015, ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S5.xlsx"], "HIV- (CCC)");
sc5.cum.stpos = strd2015_numcum_results(sc5.cc_incid_posall(:, 70:end, :), pop2015, ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S5.xlsx"], "HIV+ (CCC)" );
sc5.cum.start = strd2015_numcum_results(sc5.cc_incid_art(:, 70:end, :), pop2015, ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S5.xlsx"], "HIV+ ART (CCC)" );
sc5.cum.stnart = strd2015_numcum_results(sc5.cc_incid_pos(:, 70:end, :), pop2015, ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S5.xlsx"], "HIV+ no ART (CCC)" );


%%%%%%%

for i=1:length(sc0.stneg(:,1))
    for t=1:170
         sc5.prreduc(i,t) = 1 - sc5.stall(i,t) / sc0.stall(i,t);
    end
end

sc5.summ(:,1) = 2019:2120;
sc5.summ(:,2) = prctile(sc5.prreduc(:,69:end), 50);
sc5.summ(:,3) = prctile(sc5.prreduc(:,69:end), 5);
sc5.summ(:,4) = prctile(sc5.prreduc(:,69:end), 95);

tabnew2 = array2table(sc5.summ);
writetable(tabnew2, "/Users/minttu/who_results_sheets/ICL_HSPH_PercentReduct_S5.xlsx", 'WriteVariableNames',0, 'Sheet', "All")


for i=1:length(sc0.stpos(:,1))
    for t=1:135
         sc5.prreduc2(i,t) = 1 - sc5.stpos(i,t) / sc0.stpos(i,t);
    end
end

sc5.summ2(:,1) = 2019:2120;
sc5.summ2(:,2) = prctile(sc5.prreduc2(:,34:end), 50);
sc5.summ2(:,3) = prctile(sc5.prreduc2(:,34:end), 5);
sc5.summ2(:,4) = prctile(sc5.prreduc2(:,34:end), 95);

tabnew3 = array2table(sc5.summ2);
writetable(tabnew3, "/Users/minttu/who_results_sheets/ICL_HSPH_PercentReduct_S5.xlsx", 'WriteVariableNames',0, 'Sheet', "HIV+")


%% Sc6
[sc6.cc_incid_all, sc6.cc_incid_neg, sc6.cc_incid_posall, sc6.cc_incid_art, sc6.cc_incid_pos] = ccincid(sc6.res); % pp

sc6.stall = strd2015_results(sc6.cc_incid_all, pop2015, 2);
sc6.stneg = strd2015_results(sc6.cc_incid_neg, pop2015, 2);
sc6.stpos = strd2015_results(sc6.cc_incid_posall, pop2015, 37);
sc6.start = strd2015_results(sc6.cc_incid_art, pop2015, 56);
sc6.stnart = strd2015_results(sc6.cc_incid_pos, pop2015, 37);

sc6.tab_ccall = res_table(sc6.stall(:, 69:end), sc6.cc_incid_all(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S6.xlsx"], "Pop(All) ICC");            
sc6.tab_ccneg = res_table(sc6.stneg(:, 69:end), sc6.cc_incid_neg(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S6.xlsx"], "HIV- (ICC)");
sc6.tab_ccpos = res_table(sc6.stpos(:, 34:end), sc6.cc_incid_posall(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S6.xlsx"], "HIV+ (ICC)");
sc6.tab_ccart = res_table(sc6.start(:, 15:end), sc6.cc_incid_art(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S6.xlsx"], "HIV+ ART (ICC)");
sc6.tab_ccnart = res_table(sc6.stnart(:, 34:end), sc6.cc_incid_pos(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S6.xlsx"], "HIV+ no ART (ICC)");

sc6.tab_pops = res_poptable(sc6.res, ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S6.xlsx"], 102);

sc6.tab_ccc = res_ccctable(sc6.res, ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S6.xlsx"]);

sc6.cum.stall = strd2015_numcum_results(sc6.cc_incid_all(:, 70:end, :), pop2015, ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S6.xlsx"], "Pop(All) CCC");
sc6.cum.stneg = strd2015_numcum_results(sc6.cc_incid_neg(:, 70:end, :), pop2015, ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S6.xlsx"], "HIV- (CCC)");
sc6.cum.stpos = strd2015_numcum_results(sc6.cc_incid_posall(:, 70:end, :), pop2015, ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S6.xlsx"], "HIV+ (CCC)" );
sc6.cum.start = strd2015_numcum_results(sc6.cc_incid_art(:, 70:end, :), pop2015, ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S6.xlsx"], "HIV+ ART (CCC)" );
sc6.cum.stnart = strd2015_numcum_results(sc6.cc_incid_pos(:, 70:end, :), pop2015, ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S6.xlsx"], "HIV+ no ART (CCC)" );

sc6.tab_coverage = res_covtable(sc6.res, ["/Users/minttu/who_results_sheets/ICL_HSPH_Coverage-(2020-2120)_S6.xlsx"], 101);

% 5% LUB

sc6.tab_ccall_p5 = res_table_prctile(sc6.stall(:, 70:end), sc6.cc_incid_all(:, 71:end, :), 5, ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)LUB_S6.xlsx"], "Pop(All) ICC");            
sc6.tab_ccneg_p5 = res_table_prctile(sc6.stneg(:, 70:end), sc6.cc_incid_neg(:, 71:end, :), 5,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)LUB_S6.xlsx"], "HIV- (ICC)");
sc6.tab_ccpos_p5 = res_table_prctile(sc6.stpos(:, 35:end), sc6.cc_incid_posall(:, 71:end, :), 5,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)LUB_S6.xlsx"], "HIV+ (ICC)");
sc6.tab_ccart_p5 = res_table_prctile(sc6.start(:, 16:end), sc6.cc_incid_art(:, 71:end, :), 5,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)LUB_S6.xlsx"], "HIV+ ART (ICC)");
sc6.tab_ccnart_p5 = res_table_prctile(sc6.stnart(:, 35:end), sc6.cc_incid_pos(:, 71:end, :), 5,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)LUB_S6.xlsx"], "HIV+ no ART (ICC)");

sc6.tab_pops_p5 = res_poptable_prctile(sc6.res, 5, ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)LUB_S6.xlsx"], 101);


% 95% LUB

sc6.tab_ccall_p95 = res_table_prctile(sc6.stall(:, 70:end), sc6.cc_incid_all(:, 71:end, :), 95, ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)HUB_S6.xlsx"], "Pop(All) ICC");            
sc6.tab_ccneg_p95 = res_table_prctile(sc6.stneg(:, 70:end), sc6.cc_incid_neg(:, 71:end, :), 95,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)HUB_S6.xlsx"], "HIV- (ICC)");
sc6.tab_ccpos_p95 = res_table_prctile(sc6.stpos(:, 35:end), sc6.cc_incid_posall(:, 71:end, :), 95,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)HUB_S6.xlsx"], "HIV+ (ICC)");
sc6.tab_ccart_p95 = res_table_prctile(sc6.start(:, 16:end), sc6.cc_incid_art(:, 71:end, :), 95,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)HUB_S6.xlsx"], "HIV+ ART (ICC)");
sc6.tab_ccnart_p95 = res_table_prctile(sc6.stnart(:, 35:end), sc6.cc_incid_pos(:, 71:end, :), 95,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)HUB_S6.xlsx"], "HIV+ no ART (ICC)");

sc6.tab_pops_p95 = res_poptable_prctile(sc6.res, 95, ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)HUB_S6.xlsx"], 101);

% hr HPV women
[sc6.hrhpvf_all, sc6.hrhpvf_neg, sc6.hrhpvf_posall, sc6.hrhpvf_art, sc6.hrhpvf_pos, ...
 sc6.hrhpvf_all_num, sc6.hrhpvf_neg_num, sc6.hrhpvf_posall_num, sc6.hrhpvf_art_num, sc6.hrhpvf_pos_num] = hrhpvprev(sc6.res); % pp

sc6.stall_hrhpvf = strd2015_results_num(sc6.hrhpvf_all, pop2015, 2);
sc6.stneg_hrhpvf = strd2015_results_num(sc6.hrhpvf_neg, pop2015, 2);
sc6.stpos_hrhpvf = strd2015_results_num(sc6.hrhpvf_posall, pop2015, 37);
sc6.start_hrhpvf = strd2015_results_num(sc6.hrhpvf_art, pop2015, 56);
sc6.stnart_hrhpvf = strd2015_results_num(sc6.hrhpvf_pos, pop2015, 37);

% here weird spacing
sc6.tab_all_hrhpvf = res_table_prctile_num(sc6.stall_hrhpvf(:, 70:end), sc6.hrhpvf_all_num(:, 71:end, :), 50, ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_hrHPV_Prevalence-standardised-(2020-2120)_S6.xlsx"], "Pop(All) (hrHPV)");            
sc6.tab_neg_hrhpvf = res_table_prctile_num(sc6.stneg_hrhpvf(:, 70:end), sc6.hrhpvf_neg_num(:, 71:end, :), 50,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_hrHPV_Prevalence-standardised-(2020-2120)_S6.xlsx"], "HIV-  (hrHPV)");
sc6.tab_pos_hrhpvf = res_table_prctile_num(sc6.stpos_hrhpvf(:, 35:end), sc6.hrhpvf_posall_num(:, 71:end, :), 50,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_hrHPV_Prevalence-standardised-(2020-2120)_S6.xlsx"], "HIV+   (hrHPV)"); % this has 3 spaces
sc6.tab_art_hrhpvf = res_table_prctile_num(sc6.start_hrhpvf(:, 16:end), sc6.hrhpvf_art_num(:, 71:end, :), 50,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_hrHPV_Prevalence-standardised-(2020-2120)_S6.xlsx"], "HIV+ ART  (hrHPV)");
sc6.tab_nart_hrhpvf = res_table_prctile_num(sc6.stnart_hrhpvf(:, 35:end), sc6.hrhpvf_pos_num(:, 71:end, :), 50,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_hrHPV_Prevalence-standardised-(2020-2120)_S6.xlsx"], "HIV+ no ART  (hrHPV)");

sc6.tab_pops_hrhpvf = res_poptable_prctile(sc6.res, 50, ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_hrHPV_Prevalence-standardised-(2020-2120)_S6.xlsx"], 101);

% VT-hpv women

[sc6.hpv9vtf_all, sc6.hpv9vtf_neg, sc6.hpv9vtf_posall, sc6.hpv9vtf_art, sc6.hpv9vtf_pos, ...
 sc6.hpv9vtf_all_num, sc6.hpv9vtf_neg_num, sc6.hpv9vtf_posall_num, sc6.hpv9vtf_art_num, sc6.hpv9vtf_pos_num] = vthpvprev(sc6.res); % pp

sc6.stall_hpv9vtf = strd2015_results_num(sc6.hpv9vtf_all, pop2015, 2);
sc6.stneg_hpv9vtf = strd2015_results_num(sc6.hpv9vtf_neg, pop2015, 2);
sc6.stpos_hpv9vtf = strd2015_results_num(sc6.hpv9vtf_posall, pop2015, 37);
sc6.start_hpv9vtf = strd2015_results_num(sc6.hpv9vtf_art, pop2015, 56);
sc6.stnart_hpv9vtf = strd2015_results_num(sc6.hpv9vtf_pos, pop2015, 37);

% here weird spacing
sc6.tab_all_hpv9vtf = res_table_prctile_num(sc6.stall_hpv9vtf(:, 70:end), sc6.hpv9vtf_all_num(:, 71:end, :), 50, ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_VThrHPV_Prevalence-standardised-(2020-2120)_S6.xlsx"], "Pop(All) (vtHPV)");            
sc6.tab_neg_hpv9vtf = res_table_prctile_num(sc6.stneg_hpv9vtf(:, 70:end), sc6.hpv9vtf_neg_num(:, 71:end, :), 50,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_VThrHPV_Prevalence-standardised-(2020-2120)_S6.xlsx"], "HIV-  (vtHPV)");
sc6.tab_pos_hpv9vtf = res_table_prctile_num(sc6.stpos_hpv9vtf(:, 35:end), sc6.hpv9vtf_posall_num(:, 71:end, :), 50,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_VThrHPV_Prevalence-standardised-(2020-2120)_S6.xlsx"], "HIV+   (vtHPV)"); % this has 3 spaces
sc6.tab_art_hpv9vtf = res_table_prctile_num(sc6.start_hpv9vtf(:, 16:end), sc6.hpv9vtf_art_num(:, 71:end, :), 50,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_VThrHPV_Prevalence-standardised-(2020-2120)_S6.xlsx"], "HIV+ ART  (vtHPV)");
sc6.tab_nart_hpv9vtf = res_table_prctile_num(sc6.stnart_hpv9vtf(:, 35:end), sc6.hpv9vtf_pos_num(:, 71:end, :), 50,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_VThrHPV_Prevalence-standardised-(2020-2120)_S6.xlsx"], "HIV+ no ART  (vtHPV)");

sc6.tab_pops_hpv9vtf = res_poptable_prctile(sc6.res, 50, ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_VThrHPV_Prevalence-standardised-(2020-2120)_S6.xlsx"], 101);

% CIN2+ prevalence women
[sc6.cinprf_all, sc6.cinprf_neg, sc6.cinprf_posall, sc6.cinprf_art, sc6.cinprf_pos, ...
 sc6.cinprf_all_num, sc6.cinprf_neg_num, sc6.cinprf_posall_num, sc6.cinprf_art_num, sc6.cinprf_pos_num] = cin2prev(sc6.res); % pp

sc6.stall_cinprf = strd2015_results_num(sc6.cinprf_all, pop2015, 2);
sc6.stneg_cinprf = strd2015_results_num(sc6.cinprf_neg, pop2015, 2);
sc6.stpos_cinprf = strd2015_results_num(sc6.cinprf_posall, pop2015, 37);
sc6.start_cinprf = strd2015_results_num(sc6.cinprf_art, pop2015, 56);
sc6.stnart_cinprf = strd2015_results_num(sc6.cinprf_pos, pop2015, 37);

% here weird spacing
sc6.tab_all_cinprf = res_table_prctile_num(sc6.stall_cinprf(:, 70:end), sc6.cinprf_all_num(:, 71:end, :), 50, ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CIN2+_Prevalence-standardised-(2020-2120)_S6.xlsx"], "Pop(All) (CIN2+)");            
sc6.tab_neg_cinprf = res_table_prctile_num(sc6.stneg_cinprf(:, 70:end), sc6.cinprf_neg_num(:, 71:end, :), 50,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CIN2+_Prevalence-standardised-(2020-2120)_S6.xlsx"], "HIV- (CIN2+)"); % 1 space only
sc6.tab_pos_cinprf = res_table_prctile_num(sc6.stpos_cinprf(:, 35:end), sc6.cinprf_posall_num(:, 71:end, :), 50,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CIN2+_Prevalence-standardised-(2020-2120)_S6.xlsx"], "HIV+  (CIN2+)"); % this has 2 spaces
sc6.tab_art_cinprf = res_table_prctile_num(sc6.start_cinprf(:, 16:end), sc6.cinprf_art_num(:, 71:end, :), 50,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CIN2+_Prevalence-standardised-(2020-2120)_S6.xlsx"], "HIV+ ART  (CIN2+)");
sc6.tab_nart_cinprf = res_table_prctile_num(sc6.stnart_cinprf(:, 35:end), sc6.cinprf_pos_num(:, 71:end, :), 50,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CIN2+_Prevalence-standardised-(2020-2120)_S6.xlsx"], "HIV+ no ART  (CIN2+)");

sc6.tab_pops_cinprf = res_poptable_prctile(sc6.res, 50, ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CIN2+_Prevalence-standardised-(2020-2120)_S6.xlsx"], 101);


%%%%%%%

for i=1:length(sc0.stneg(:,1))
    for t=1:170
         sc6.prreduc(i,t) = 1 - sc6.stall(i,t) / sc0.stall(i,t);
    end
end

sc6.summ(:,1) = 2019:2120;
sc6.summ(:,2) = prctile(sc6.prreduc(:,69:end), 50);
sc6.summ(:,3) = prctile(sc6.prreduc(:,69:end), 5);
sc6.summ(:,4) = prctile(sc6.prreduc(:,69:end), 95);

tabnew2 = array2table(sc6.summ);
writetable(tabnew2, "/Users/minttu/who_results_sheets/ICL_HSPH_PercentReduct_S6.xlsx", 'WriteVariableNames',0, 'Sheet', "All")


for i=1:length(sc0.stpos(:,1))
    for t=1:135
         sc6.prreduc2(i,t) = 1 - sc6.stpos(i,t) / sc0.stpos(i,t);
    end
end

sc6.summ2(:,1) = 2019:2120;
sc6.summ2(:,2) = prctile(sc6.prreduc2(:,34:end), 50);
sc6.summ2(:,3) = prctile(sc6.prreduc2(:,34:end), 5);
sc6.summ2(:,4) = prctile(sc6.prreduc2(:,34:end), 95);

tabnew3 = array2table(sc6.summ2);
writetable(tabnew3, "/Users/minttu/who_results_sheets/ICL_HSPH_PercentReduct_S6.xlsx", 'WriteVariableNames',0, 'Sheet', "HIV+")


%% Sc7
[sc7.cc_incid_all, sc7.cc_incid_neg, sc7.cc_incid_posall, sc7.cc_incid_art, sc7.cc_incid_pos] = ccincid(sc7.res); % pp

sc7.stall = strd2015_results(sc7.cc_incid_all, pop2015, 2);
sc7.stneg = strd2015_results(sc7.cc_incid_neg, pop2015, 2);
sc7.stpos = strd2015_results(sc7.cc_incid_posall, pop2015, 37);
sc7.start = strd2015_results(sc7.cc_incid_art, pop2015, 56);
sc7.stnart = strd2015_results(sc7.cc_incid_pos, pop2015, 37);
 

sc7.tab_ccall = res_table(sc7.stall(:, 69:end), sc7.cc_incid_all(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S7.xlsx"], "Pop(All) ICC");            
sc7.tab_ccneg = res_table(sc7.stneg(:, 69:end), sc7.cc_incid_neg(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S7.xlsx"], "HIV- (ICC)");
sc7.tab_ccpos = res_table(sc7.stpos(:, 34:end), sc7.cc_incid_posall(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S7.xlsx"], "HIV+ (ICC)");
sc7.tab_ccart = res_table(sc7.start(:, 15:end), sc7.cc_incid_art(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S7.xlsx"], "HIV+ ART (ICC)");
sc7.tab_ccnart = res_table(sc7.stnart(:, 34:end), sc7.cc_incid_pos(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S7.xlsx"], "HIV+ no ART (ICC)");

sc7.tab_pops = res_poptable(sc7.res, ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S7.xlsx"], 102);

sc7.tab_ccc = res_ccctable(sc7.res, ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S7.xlsx"]);

sc7.cum.stall = strd2015_numcum_results(sc7.cc_incid_all(:, 70:end, :), pop2015,["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S7.xlsx"], "Pop(All) CCC");
sc7.cum.stneg = strd2015_numcum_results(sc7.cc_incid_neg(:, 70:end, :), pop2015,["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S7.xlsx"], "HIV- (CCC)");
sc7.cum.stpos = strd2015_numcum_results(sc7.cc_incid_posall(:, 70:end, :), pop2015,["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S7.xlsx"], "HIV+ (CCC)" );
sc7.cum.start = strd2015_numcum_results(sc7.cc_incid_art(:, 70:end, :), pop2015,["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S7.xlsx"], "HIV+ ART (CCC)" );
sc7.cum.stnart = strd2015_numcum_results(sc7.cc_incid_pos(:, 70:end, :), pop2015, ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S7.xlsx"], "HIV+ no ART (CCC)" );


% [sc7.vacc, sc7.vaccpos, sc7.screen] = coverageest(sc7.res);
% [sc7.vaccnum, sc7.screennum] = coverageest_nums(sc7.res);


%% Sc7a
[sc7a.cc_incid_all, sc7a.cc_incid_neg, sc7a.cc_incid_posall, sc7a.cc_incid_art, sc7a.cc_incid_pos] = ccincid(sc7a.res); % pp

sc7a.stall = strd2015_results(sc7a.cc_incid_all, pop2015, 2);
sc7a.stneg = strd2015_results(sc7a.cc_incid_neg, pop2015, 2);
sc7a.stpos = strd2015_results(sc7a.cc_incid_posall, pop2015, 37);
sc7a.start = strd2015_results(sc7a.cc_incid_art, pop2015, 56);
sc7a.stnart = strd2015_results(sc7a.cc_incid_pos, pop2015, 37);

sc7a.tab_ccall = res_table(sc7a.stall(:, 69:end), sc7a.cc_incid_all(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S7a.xlsx"], "Pop(All) ICC");            
sc7a.tab_ccneg = res_table(sc7a.stneg(:, 69:end), sc7a.cc_incid_neg(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S7a.xlsx"], "HIV- (ICC)");
sc7a.tab_ccpos = res_table(sc7a.stpos(:, 34:end), sc7a.cc_incid_posall(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S7a.xlsx"], "HIV+ (ICC)");
sc7a.tab_ccart = res_table(sc7a.start(:, 15:end), sc7a.cc_incid_art(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S7a.xlsx"], "HIV+ ART (ICC)");
sc7a.tab_ccnart = res_table(sc7a.stnart(:, 34:end), sc7a.cc_incid_pos(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S7a.xlsx"], "HIV+ no ART (ICC)");

sc7a.tab_pops = res_poptable(sc7a.res, ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S7a.xlsx"], 102);

sc7a.tab_ccc = res_ccctable(sc7a.res, ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S7a.xlsx"]);

sc7a.cum.stall = strd2015_numcum_results(sc7a.cc_incid_all(:, 70:end, :), pop2015,["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S7a.xlsx"], "Pop(All) CCC");
sc7a.cum.stneg = strd2015_numcum_results(sc7a.cc_incid_neg(:, 70:end, :), pop2015,["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S7a.xlsx"], "HIV- (CCC)");
sc7a.cum.stpos = strd2015_numcum_results(sc7a.cc_incid_posall(:, 70:end, :), pop2015,["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S7a.xlsx"], "HIV+ (CCC)" );
sc7a.cum.start = strd2015_numcum_results(sc7a.cc_incid_art(:, 70:end, :), pop2015,["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S7a.xlsx"], "HIV+ ART (CCC)" );
sc7a.cum.stnart = strd2015_numcum_results(sc7a.cc_incid_pos(:, 70:end, :), pop2015, ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S7a.xlsx"], "HIV+ no ART (CCC)" );



%% Sc7b
[sc7b.cc_incid_all, sc7b.cc_incid_neg, sc7b.cc_incid_posall, sc7b.cc_incid_art, sc7b.cc_incid_pos] = ccincid(sc7b.res); % pp

sc7b.stall = strd2015_results(sc7b.cc_incid_all, pop2015, 2);
sc7b.stneg = strd2015_results(sc7b.cc_incid_neg, pop2015, 2);
sc7b.stpos = strd2015_results(sc7b.cc_incid_posall, pop2015, 37);
sc7b.start = strd2015_results(sc7b.cc_incid_art, pop2015, 56);
sc7b.stnart = strd2015_results(sc7b.cc_incid_pos, pop2015, 37);

sc7b.tab_ccall = res_table(sc7b.stall(:, 69:end), sc7b.cc_incid_all(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S7b.xlsx"], "Pop(All) ICC");            
sc7b.tab_ccneg = res_table(sc7b.stneg(:, 69:end), sc7b.cc_incid_neg(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S7b.xlsx"], "HIV- (ICC)");
sc7b.tab_ccpos = res_table(sc7b.stpos(:, 34:end), sc7b.cc_incid_posall(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S7b.xlsx"], "HIV+ (ICC)");
sc7b.tab_ccart = res_table(sc7b.start(:, 15:end), sc7b.cc_incid_art(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S7b.xlsx"], "HIV+ ART (ICC)");
sc7b.tab_ccnart = res_table(sc7b.stnart(:, 34:end), sc7b.cc_incid_pos(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S7b.xlsx"], "HIV+ no ART (ICC)");

sc7b.tab_pops = res_poptable(sc7b.res, ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S7b.xlsx"], 102);

sc7b.tab_ccc = res_ccctable(sc7b.res, ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S7b.xlsx"]);

sc7b.cum.stall = strd2015_numcum_results(sc7b.cc_incid_all(:, 70:end, :), pop2015,["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S7b.xlsx"], "Pop(All) CCC");
sc7b.cum.stneg = strd2015_numcum_results(sc7b.cc_incid_neg(:, 70:end, :), pop2015,["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S7b.xlsx"], "HIV- (CCC)");
sc7b.cum.stpos = strd2015_numcum_results(sc7b.cc_incid_posall(:, 70:end, :), pop2015,["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S7b.xlsx"], "HIV+ (CCC)" );
sc7b.cum.start = strd2015_numcum_results(sc7b.cc_incid_art(:, 70:end, :), pop2015,["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S7b.xlsx"], "HIV+ ART (CCC)" );
sc7b.cum.stnart = strd2015_numcum_results(sc7b.cc_incid_pos(:, 70:end, :), pop2015, ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S7b.xlsx"], "HIV+ no ART (CCC)" );


%% Sc8
[sc8.cc_incid_all, sc8.cc_incid_neg, sc8.cc_incid_posall, sc8.cc_incid_art, sc8.cc_incid_pos] = ccincid(sc8.res); % pp

sc8.stall = strd2015_results(sc8.cc_incid_all, pop2015, 2);
sc8.stneg = strd2015_results(sc8.cc_incid_neg, pop2015, 2);
sc8.stpos = strd2015_results(sc8.cc_incid_posall, pop2015, 37);
sc8.start = strd2015_results(sc8.cc_incid_art, pop2015, 56);
sc8.stnart = strd2015_results(sc8.cc_incid_pos, pop2015, 37);

sc8.tab_ccall = res_table(sc8.stall(:, 69:end), sc8.cc_incid_all(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S8.xlsx"], "Pop(All) ICC");            
sc8.tab_ccneg = res_table(sc8.stneg(:, 69:end), sc8.cc_incid_neg(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S8.xlsx"], "HIV- (ICC)");
sc8.tab_ccpos = res_table(sc8.stpos(:, 34:end), sc8.cc_incid_posall(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S8.xlsx"], "HIV+ (ICC)");
sc8.tab_ccart = res_table(sc8.start(:, 15:end), sc8.cc_incid_art(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S8.xlsx"], "HIV+ ART (ICC)");
sc8.tab_ccnart = res_table(sc8.stnart(:, 34:end), sc8.cc_incid_pos(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S8.xlsx"], "HIV+ no ART (ICC)");

sc8.tab_pops = res_poptable(sc8.res, ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S8.xlsx"], 102);

sc8.tab_ccc = res_ccctable(sc8.res, ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S8.xlsx"]);

sc8.cum.stall = strd2015_numcum_results(sc8.cc_incid_all(:, 70:end, :), pop2015,  ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S8.xlsx"], "Pop(All) CCC");
sc8.cum.stneg = strd2015_numcum_results(sc8.cc_incid_neg(:, 70:end, :), pop2015, ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S8.xlsx"], "HIV- (CCC)");
sc8.cum.stpos = strd2015_numcum_results(sc8.cc_incid_posall(:, 70:end, :), pop2015, ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S8.xlsx"], "HIV+ (CCC)" );
sc8.cum.start = strd2015_numcum_results(sc8.cc_incid_art(:, 70:end, :), pop2015, ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S8.xlsx"], "HIV+ ART (CCC)" );
sc8.cum.stnart = strd2015_numcum_results(sc8.cc_incid_pos(:, 70:end, :), pop2015, ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S8.xlsx"], "HIV+ no ART (CCC)" );


% [sc8.vacc, sc8.vaccpos, sc8.screen] = coverageest(sc8.res);
% [sc8.vaccnum, sc8.screennum] = coverageest_nums(sc8.res);

sc8.tab_coverage = res_covtable(sc8.res, ["/Users/minttu/who_results_sheets/ICL_HSPH_Coverage-(2020-2120)_S8.xlsx"], 101);

% LL & UL

% 5% LUB

sc8.tab_ccall_p5 = res_table_prctile(sc8.stall(:, 70:end), sc8.cc_incid_all(:, 71:end, :), 5, ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)LUB_S8.xlsx"], "Pop(All) ICC");            
sc8.tab_ccneg_p5 = res_table_prctile(sc8.stneg(:, 70:end), sc8.cc_incid_neg(:, 71:end, :), 5,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)LUB_S8.xlsx"], "HIV- (ICC)");
sc8.tab_ccpos_p5 = res_table_prctile(sc8.stpos(:, 35:end), sc8.cc_incid_posall(:, 71:end, :), 5,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)LUB_S8.xlsx"], "HIV+ (ICC)");
sc8.tab_ccart_p5 = res_table_prctile(sc8.start(:, 16:end), sc8.cc_incid_art(:, 71:end, :), 5,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)LUB_S8.xlsx"], "HIV+ ART (ICC)");
sc8.tab_ccnart_p5 = res_table_prctile(sc8.stnart(:, 35:end), sc8.cc_incid_pos(:, 71:end, :), 5,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)LUB_S8.xlsx"], "HIV+ no ART (ICC)");

sc8.tab_pops_p5 = res_poptable_prctile(sc8.res, 5, ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)LUB_S8.xlsx"], 101);


% 95% LUB

sc8.tab_ccall_p95 = res_table_prctile(sc8.stall(:, 70:end), sc8.cc_incid_all(:, 71:end, :), 95, ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)HUB_S8.xlsx"], "Pop(All) ICC");            
sc8.tab_ccneg_p95 = res_table_prctile(sc8.stneg(:, 70:end), sc8.cc_incid_neg(:, 71:end, :), 95,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)HUB_S8.xlsx"], "HIV- (ICC)");
sc8.tab_ccpos_p95 = res_table_prctile(sc8.stpos(:, 35:end), sc8.cc_incid_posall(:, 71:end, :), 95,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)HUB_S8.xlsx"], "HIV+ (ICC)");
sc8.tab_ccart_p95 = res_table_prctile(sc8.start(:, 16:end), sc8.cc_incid_art(:, 71:end, :), 95,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)HUB_S8.xlsx"], "HIV+ ART (ICC)");
sc8.tab_ccnart_p95 = res_table_prctile(sc8.stnart(:, 35:end), sc8.cc_incid_pos(:, 71:end, :), 95,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)HUB_S8.xlsx"], "HIV+ no ART (ICC)");

sc8.tab_pops_p95 = res_poptable_prctile(sc8.res, 95, ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)HUB_S8.xlsx"], 101);


%%%%%%%

for i=1:length(sc0.stneg(:,1))
    for t=1:170
         sc8.prreduc(i,t) = 1 - sc8.stall(i,t) / sc0.stall(i,t);
    end
end

sc8.summ(:,1) = 2019:2120;
sc8.summ(:,2) = prctile(sc8.prreduc(:,69:end), 50);
sc8.summ(:,3) = prctile(sc8.prreduc(:,69:end), 5);
sc8.summ(:,4) = prctile(sc8.prreduc(:,69:end), 95);

tabnew2 = array2table(sc8.summ);
writetable(tabnew2, "/Users/minttu/who_results_sheets/ICL_HSPH_PercentReduct_S8.xlsx", 'WriteVariableNames',0, 'Sheet', "All")


for i=1:length(sc0.stpos(:,1))
    for t=1:135
         sc8.prreduc2(i,t) = 1 - sc8.stpos(i,t) / sc0.stpos(i,t);
    end
end

sc8.summ2(:,1) = 2019:2120;
sc8.summ2(:,2) = prctile(sc8.prreduc2(:,34:end), 50);
sc8.summ2(:,3) = prctile(sc8.prreduc2(:,34:end), 5);
sc8.summ2(:,4) = prctile(sc8.prreduc2(:,34:end), 95);

tabnew3 = array2table(sc8.summ2);
writetable(tabnew3, "/Users/minttu/who_results_sheets/ICL_HSPH_PercentReduct_S8.xlsx", 'WriteVariableNames',0, 'Sheet', "HIV+")




%% Sc9
[sc9.cc_incid_all, sc9.cc_incid_neg, sc9.cc_incid_posall, sc9.cc_incid_art, sc9.cc_incid_pos] = ccincid(sc9.res); % pp

sc9.stall = strd2015_results(sc9.cc_incid_all, pop2015, 2);
sc9.stneg = strd2015_results(sc9.cc_incid_neg, pop2015, 2);
sc9.stpos = strd2015_results(sc9.cc_incid_posall, pop2015, 37);
sc9.start = strd2015_results(sc9.cc_incid_art, pop2015, 56);
sc9.stnart = strd2015_results(sc9.cc_incid_pos, pop2015, 37);

% [sc9.vacc, sc9.screen] = coverageest(sc9.res);

sc9.tab_ccall = res_table(sc9.stall(:, 69:end), sc9.cc_incid_all(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S9.xlsx"], "Pop(All) ICC");            
sc9.tab_ccneg = res_table(sc9.stneg(:, 69:end), sc9.cc_incid_neg(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S9.xlsx"], "HIV- (ICC)");
sc9.tab_ccpos = res_table(sc9.stpos(:, 34:end), sc9.cc_incid_posall(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S9.xlsx"], "HIV+ (ICC)");
sc9.tab_ccart = res_table(sc9.start(:, 15:end), sc9.cc_incid_art(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S9.xlsx"], "HIV+ ART (ICC)");
sc9.tab_ccnart = res_table(sc9.stnart(:, 34:end), sc9.cc_incid_pos(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S9.xlsx"], "HIV+ no ART (ICC)");

sc9.tab_pops = res_poptable(sc9.res, ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S9.xlsx"], 102);

sc9.tab_ccc = res_ccctable(sc9.res, ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S9.xlsx"]);

sc9.cum.stall = strd2015_numcum_results(sc9.cc_incid_all(:, 70:end, :), pop2015, ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S9.xlsx"], "Pop(All) CCC");
sc9.cum.stneg = strd2015_numcum_results(sc9.cc_incid_neg(:, 70:end, :), pop2015, ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S9.xlsx"], "HIV- (CCC)");
sc9.cum.stpos = strd2015_numcum_results(sc9.cc_incid_posall(:, 70:end, :), pop2015, ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S9.xlsx"], "HIV+ (CCC)" );
sc9.cum.start = strd2015_numcum_results(sc9.cc_incid_art(:, 70:end, :), pop2015, ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S9.xlsx"], "HIV+ ART (CCC)" );
sc9.cum.stnart = strd2015_numcum_results(sc9.cc_incid_pos(:, 70:end, :), pop2015, ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S9.xlsx"], "HIV+ no ART (CCC)" );

% UL & LL


% 5% LUB

sc9.tab_ccall_p5 = res_table_prctile(sc9.stall(:, 70:end), sc9.cc_incid_all(:, 71:end, :), 5, ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)LUB_S9.xlsx"], "Pop(All) ICC");            
sc9.tab_ccneg_p5 = res_table_prctile(sc9.stneg(:, 70:end), sc9.cc_incid_neg(:, 71:end, :), 5,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)LUB_S9.xlsx"], "HIV- (ICC)");
sc9.tab_ccpos_p5 = res_table_prctile(sc9.stpos(:, 35:end), sc9.cc_incid_posall(:, 71:end, :), 5,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)LUB_S9.xlsx"], "HIV+ (ICC)");
sc9.tab_ccart_p5 = res_table_prctile(sc9.start(:, 16:end), sc9.cc_incid_art(:, 71:end, :), 5,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)LUB_S9.xlsx"], "HIV+ ART (ICC)");
sc9.tab_ccnart_p5 = res_table_prctile(sc9.stnart(:, 35:end), sc9.cc_incid_pos(:, 71:end, :), 5,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)LUB_S9.xlsx"], "HIV+ no ART (ICC)");

sc9.tab_pops_p5 = res_poptable_prctile(sc9.res, 5, ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)LUB_S9.xlsx"], 101);


% 95% LUB

sc9.tab_ccall_p95 = res_table_prctile(sc9.stall(:, 70:end), sc9.cc_incid_all(:, 71:end, :), 95, ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)HUB_S9.xlsx"], "Pop(All) ICC");            
sc9.tab_ccneg_p95 = res_table_prctile(sc9.stneg(:, 70:end), sc9.cc_incid_neg(:, 71:end, :), 95,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)HUB_S9.xlsx"], "HIV- (ICC)");
sc9.tab_ccpos_p95 = res_table_prctile(sc9.stpos(:, 35:end), sc9.cc_incid_posall(:, 71:end, :), 95,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)HUB_S9.xlsx"], "HIV+ (ICC)");
sc9.tab_ccart_p95 = res_table_prctile(sc9.start(:, 16:end), sc9.cc_incid_art(:, 71:end, :), 95,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)HUB_S9.xlsx"], "HIV+ ART (ICC)");
sc9.tab_ccnart_p95 = res_table_prctile(sc9.stnart(:, 35:end), sc9.cc_incid_pos(:, 71:end, :), 95,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)HUB_S9.xlsx"], "HIV+ no ART (ICC)");

sc9.tab_pops_p95 = res_poptable_prctile(sc9.res, 95, ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)HUB_S9.xlsx"], 101);

%%%%%

for i=1:length(sc0.stneg(:,1))
    for t=1:170
         sc9.prreduc(i,t) = 1 - sc9.stall(i,t) / sc0.stall(i,t);
    end
end

sc9.summ(:,1) = 2019:2120;
sc9.summ(:,2) = prctile(sc9.prreduc(:,69:end), 50);
sc9.summ(:,3) = prctile(sc9.prreduc(:,69:end), 5);
sc9.summ(:,4) = prctile(sc9.prreduc(:,69:end), 95);

tabnew2 = array2table(sc9.summ);
writetable(tabnew2, "/Users/minttu/who_results_sheets/ICL_HSPH_PercentReduct_S9.xlsx", 'WriteVariableNames',0, 'Sheet', "All")


for i=1:length(sc0.stpos(:,1))
    for t=1:135
         sc9.prreduc2(i,t) = 1 - sc9.stpos(i,t) / sc0.stpos(i,t);
    end
end

sc9.summ2(:,1) = 2019:2120;
sc9.summ2(:,2) = prctile(sc9.prreduc2(:,34:end), 50);
sc9.summ2(:,3) = prctile(sc9.prreduc2(:,34:end), 5);
sc9.summ2(:,4) = prctile(sc9.prreduc2(:,34:end), 95);

tabnew3 = array2table(sc9.summ2);
writetable(tabnew3, "/Users/minttu/who_results_sheets/ICL_HSPH_PercentReduct_S9.xlsx", 'WriteVariableNames',0, 'Sheet', "HIV+")



%% Sc10
[sc10.cc_incid_all, sc10.cc_incid_neg, sc10.cc_incid_posall, sc10.cc_incid_art, sc10.cc_incid_pos] = ccincid(sc10.res); % pp

sc10.stall = strd2015_results(sc10.cc_incid_all, pop2015, 2);
sc10.stneg = strd2015_results(sc10.cc_incid_neg, pop2015, 2);
sc10.stpos = strd2015_results(sc10.cc_incid_posall, pop2015, 37);
sc10.start = strd2015_results(sc10.cc_incid_art, pop2015, 56);
sc10.stnart = strd2015_results(sc10.cc_incid_pos, pop2015, 37);

sc10.tab_ccall = res_table(sc10.stall(:, 69:end), sc10.cc_incid_all(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S10.xlsx"], "Pop(All) ICC");            
sc10.tab_ccneg = res_table(sc10.stneg(:, 69:end), sc10.cc_incid_neg(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S10.xlsx"], "HIV- (ICC)");
sc10.tab_ccpos = res_table(sc10.stpos(:, 34:end), sc10.cc_incid_posall(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S10.xlsx"], "HIV+ (ICC)");
sc10.tab_ccart = res_table(sc10.start(:, 15:end), sc10.cc_incid_art(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S10.xlsx"], "HIV+ ART (ICC)");
sc10.tab_ccnart = res_table(sc10.stnart(:, 34:end), sc10.cc_incid_pos(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S10.xlsx"], "HIV+ no ART (ICC)");

sc10.tab_pops = res_poptable(sc10.res, ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S10.xlsx"], 102);

sc10.tab_ccc = res_ccctable(sc10.res, ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S10.xlsx"]);

sc10.cum.stall = strd2015_numcum_results(sc10.cc_incid_all(:, 70:end, :), pop2015, ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S10.xlsx"], "Pop(All) CCC");
sc10.cum.stneg = strd2015_numcum_results(sc10.cc_incid_neg(:, 70:end, :), pop2015, ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S10.xlsx"], "HIV- (CCC)");
sc10.cum.stpos = strd2015_numcum_results(sc10.cc_incid_posall(:, 70:end, :), pop2015, ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S10.xlsx"], "HIV+ (CCC)" );
sc10.cum.start = strd2015_numcum_results(sc10.cc_incid_art(:, 70:end, :), pop2015, ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S10.xlsx"], "HIV+ ART (CCC)" );
sc10.cum.stnart = strd2015_numcum_results(sc10.cc_incid_pos(:, 70:end, :), pop2015, ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S10.xlsx"], "HIV+ no ART (CCC)" );


%% Sc11
[sc11.cc_incid_all, sc11.cc_incid_neg, sc11.cc_incid_posall, sc11.cc_incid_art, sc11.cc_incid_pos] = ccincid(sc11.res); % pp

sc11.stall = strd2015_results(sc11.cc_incid_all, pop2015, 2);
sc11.stneg = strd2015_results(sc11.cc_incid_neg, pop2015, 2);
sc11.stpos = strd2015_results(sc11.cc_incid_posall, pop2015, 37);
sc11.start = strd2015_results(sc11.cc_incid_art, pop2015, 56);
sc11.stnart = strd2015_results(sc11.cc_incid_pos, pop2015, 37);

sc11.tab_ccall = res_table(sc11.stall(:, 69:end), sc11.cc_incid_all(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S11.xlsx"], "Pop(All) ICC");            
sc11.tab_ccneg = res_table(sc11.stneg(:, 69:end), sc11.cc_incid_neg(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S11.xlsx"], "HIV- (ICC)");
sc11.tab_ccpos = res_table(sc11.stpos(:, 34:end), sc11.cc_incid_posall(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S11.xlsx"], "HIV+ (ICC)");
sc11.tab_ccart = res_table(sc11.start(:, 15:end), sc11.cc_incid_art(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S11.xlsx"], "HIV+ ART (ICC)");
sc11.tab_ccnart = res_table(sc11.stnart(:, 34:end), sc11.cc_incid_pos(:, 70:end, :), ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S11.xlsx"], "HIV+ no ART (ICC)");

sc11.tab_pops = res_poptable(sc11.res, ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)_S11.xlsx"], 102);

sc11.tab_ccc = res_ccctable(sc11.res, ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S11.xlsx"]);

sc11.cum.stall = strd2015_numcum_results(sc11.cc_incid_all(:, 70:end, :), pop2015, ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S11.xlsx"], "Pop(All) CCC");
sc11.cum.stneg = strd2015_numcum_results(sc11.cc_incid_neg(:, 70:end, :), pop2015, ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S11.xlsx"], "HIV- (CCC)");
sc11.cum.stpos = strd2015_numcum_results(sc11.cc_incid_posall(:, 70:end, :), pop2015, ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S11.xlsx"], "HIV+ (CCC)" );
sc11.cum.start = strd2015_numcum_results(sc11.cc_incid_art(:, 70:end, :), pop2015, ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S11.xlsx"], "HIV+ ART (CCC)" );
sc11.cum.stnart = strd2015_numcum_results(sc11.cc_incid_pos(:, 70:end, :), pop2015, ["/Users/minttu/who_results_sheets/ICL_HSPH_CumulativeImpact_CC-standardised-(2020-2120)_S11.xlsx"], "HIV+ no ART (CCC)" );

% LL & UL

% 5% LUB

sc11.tab_ccall_p5 = res_table_prctile(sc11.stall(:, 70:end), sc11.cc_incid_all(:, 71:end, :), 5, ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)LUB_S11.xlsx"], "Pop(All) ICC");            
sc11.tab_ccneg_p5 = res_table_prctile(sc11.stneg(:, 70:end), sc11.cc_incid_neg(:, 71:end, :), 5,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)LUB_S11.xlsx"], "HIV- (ICC)");
sc11.tab_ccpos_p5 = res_table_prctile(sc11.stpos(:, 35:end), sc11.cc_incid_posall(:, 71:end, :), 5,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)LUB_S11.xlsx"], "HIV+ (ICC)");
sc11.tab_ccart_p5 = res_table_prctile(sc11.start(:, 16:end), sc11.cc_incid_art(:, 71:end, :), 5,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)LUB_S11.xlsx"], "HIV+ ART (ICC)");
sc11.tab_ccnart_p5 = res_table_prctile(sc11.stnart(:, 35:end), sc11.cc_incid_pos(:, 71:end, :), 5,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)LUB_S11.xlsx"], "HIV+ no ART (ICC)");

sc11.tab_pops_p5 = res_poptable_prctile(sc11.res, 5, ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)LUB_S11.xlsx"], 101);


% 95% LUB

sc11.tab_ccall_p95 = res_table_prctile(sc11.stall(:, 70:end), sc11.cc_incid_all(:, 71:end, :), 95, ...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)HUB_S11.xlsx"], "Pop(All) ICC");            
sc11.tab_ccneg_p95 = res_table_prctile(sc11.stneg(:, 70:end), sc11.cc_incid_neg(:, 71:end, :), 95,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)HUB_S11.xlsx"], "HIV- (ICC)");
sc11.tab_ccpos_p95 = res_table_prctile(sc11.stpos(:, 35:end), sc11.cc_incid_posall(:, 71:end, :), 95,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)HUB_S11.xlsx"], "HIV+ (ICC)");
sc11.tab_ccart_p95 = res_table_prctile(sc11.start(:, 16:end), sc11.cc_incid_art(:, 71:end, :), 95,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)HUB_S11.xlsx"], "HIV+ ART (ICC)");
sc11.tab_ccnart_p95 = res_table_prctile(sc11.stnart(:, 35:end), sc11.cc_incid_pos(:, 71:end, :), 95,...
                ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)HUB_S11.xlsx"], "HIV+ no ART (ICC)");

sc11.tab_pops_p95 = res_poptable_prctile(sc11.res, 95, ["/Users/minttu/who_results_sheets/ICL_HSPH_Impact_CC_IncidenceRates-standardised-(2020-2120)HUB_S11.xlsx"], 101);

%%%%%%%

for i=1:length(sc0.stneg(:,1))
    for t=1:170
         sc11.prreduc(i,t) = 1 - sc11.stall(i,t) / sc0.stall(i,t);
    end
end

sc11.summ(:,1) = 2019:2120;
sc11.summ(:,2) = prctile(sc11.prreduc(:,69:end), 50);
sc11.summ(:,3) = prctile(sc11.prreduc(:,69:end), 5);
sc11.summ(:,4) = prctile(sc11.prreduc(:,69:end), 95);

tabnew2 = array2table(sc11.summ);
writetable(tabnew2, "/Users/minttu/who_results_sheets/ICL_HSPH_PercentReduct_S11.xlsx", 'WriteVariableNames',0, 'Sheet', "All")


for i=1:length(sc0.stpos(:,1))
    for t=1:135
         sc11.prreduc2(i,t) = 1 - sc11.stpos(i,t) / sc0.stpos(i,t);
    end
end

sc11.summ2(:,1) = 2019:2120;
sc11.summ2(:,2) = prctile(sc11.prreduc2(:,34:end), 50);
sc11.summ2(:,3) = prctile(sc11.prreduc2(:,34:end), 5);
sc11.summ2(:,4) = prctile(sc11.prreduc2(:,34:end), 95);

tabnew3 = array2table(sc11.summ2);
writetable(tabnew3, "/Users/minttu/who_results_sheets/ICL_HSPH_PercentReduct_S11.xlsx", 'WriteVariableNames',0, 'Sheet', "HIV+")


%%

% baseline plot

figure
for i=1:length(sc0.stneg(:,1))
hold on
     plot(1951:2120, sc0.stneg(i,:), 'Color', col1, 'LineWidth', mruns)
     
  %  plot(1951:2120, sc0.stall(i,:), 'Color', col3, 'LineWidth', mruns)
%      
  %   plot(1986:2120, sc0.stpos(i,:), 'Color', col2, 'LineWidth', mruns)
%     % plot(1950:2120, 100000.*sc0.cc_incid_neg, 'Color', col1, 'LineWidth', mruns)
hold off
end

hold on

  
     plot(1951:2120, median(sc0.stall), 'Color', col2, 'LineWidth', 2)
     plot(1951:2120, median(sc0.stneg), 'Color', col3, 'LineWidth', 2)
     plot(1986:2120, median(sc0.stpos), 'Color', col1, 'LineWidth', 2)
     
     
     
     plot(1951:2120, median(sc0.stall), 'Color', col2, 'LineWidth', 2,  'LineStyle', '--')
     plot(1951:2120, median(sc0.stneg), 'Color', col3, 'LineWidth', 2, 'LineStyle', '--')
     plot(1986:2120, median(sc0.stpos), 'Color', col1, 'LineWidth', 2,  'LineStyle', '--')
     
    
     

yline(4)
yline(10, '-.')

   legend('Sc0.all', 'Sc0.negative', 'Sc0.positive')%, 'Sc1', 'Sc2')% 'Sc3 (Sc1+1)', 'Sc5 (Sc1+1)', 'Sc2', 'Sc4 (Sc2+1)',)
   

 axis([2020 2120 0 120])
 title({'Standardized. CC incidence by HIV'}, 'FontSize', 11); set(gca,'FontSize',11); ...
    ylabel({'per 100,000'}, 'FontSize', 11); xlabel({'Year'}, 'FontSize', 11); set(gca,'FontSize',11)

 
 filename2 = ['/Users/minttu/matlab_figures/CCI-Sc0-1-new' ,'-', datestr(now,'yyyy-mm-dd'),'.png'] ; 
saveas(gcf, filename2)


% baseline plot 2


figure
% for i=1:length(sc0.stneg(:,1))
% hold on
%      plot(2005:2120, sc0.start(i,:), 'Color', cmap(5,:), 'LineWidth', 0.01)
%      
%      plot(1986:2120, sc0.stnart(i,:), 'Color', cmap(6,:), 'LineWidth', 0.01)
%     % plot(1950:2120, 100000.*sc0.cc_incid_neg, 'Color', col1, 'LineWidth', mruns)
% %hold off
% end

hold on
     plot(2005:2120, median(sc0.start), 'Color', 'k', 'LineWidth', 2)
     plot(1986:2120, median(sc0.stnart), 'Color', 'b', 'LineWidth', 2, 'LineStyle', '--')
     
yline(4)
yline(10, '-.')

   legend('HIV on ART', 'HIV not on ART')%, 'Sc1', 'Sc2')% 'Sc3 (Sc1+1)', 'Sc5 (Sc1+1)', 'Sc2', 'Sc4 (Sc2+1)',)


 title({'Standardized. CC incidence by ART'}, 'FontSize', 11); set(gca,'FontSize',11); ...
    ylabel({'per 100,000'}, 'FontSize', 11); xlabel({'Year'}, 'FontSize', 11); set(gca,'FontSize',11)

  axis([2020 2120 0 230])
  
 filename2 = ['/Users/minttu/matlab_figures/CCI-Sc0-ART' ,'-', datestr(now,'yyyy-mm-dd'),'.png'] ; 
saveas(gcf, filename2)



% ART ratio of CCI 

art_ratio = sc0.start ./  sc0.stnart(:, 20:end);

figure
for i=1:length(sc0.stneg(:,1))
hold on

plot(2005:2120, art_ratio(i,:), 'Color', 'k')
end

 plot(2005:2120, median(art_ratio), 'Color', 'b', 'LineWidth', 3)
 
 
 title({'Ratio, CC incidence (age-standardized) on ART / not on ART'}, 'FontSize', 11); set(gca,'FontSize',11); ...
    ylabel({'per 100,000'}, 'FontSize', 11); xlabel({'Year'}, 'FontSize', 11); set(gca,'FontSize',11)

  axis([2000 2121 0.2 1.2])
  
 filename2 = ['/Users/minttu/matlab_figures/CCI-ARTratio' ,'-', datestr(now,'yyyy-mm-dd'),'.png'] ; 
saveas(gcf, filename2)
%%
% scenario plots ALL

figure


hold on
     plot(1951:2120, median(sc0.stall), 'Color', cmap(1,:), 'LineWidth', 2)
      plot(1951:2120, median(sc00.stall), 'Color', cmap(9,:), 'LineWidth', 2)
     
     plot(1951:2120, median(sc0b.stall),'Color',cmap(2,:), 'LineWidth', 1)
     plot(1951:2120, median(sc0c.stall),'Color',cmap(2,:), 'LineWidth', 1)
     
     plot(1951:2120, median(sc1.stall), 'Color', cmap(3,:), 'LineWidth', 2)
     plot(1951:2120, median(sc3.stall), 'Color', cmap(3,:), 'LineWidth', 1, 'LineStyle', ':')
     plot(1951:2120, median(sc5.stall), 'Color', cmap(3,:), 'LineWidth', 1, 'LineStyle', '-.')
           
      plot(1951:2120, median(sc2.stall), 'Color', cmap(4,:), 'LineWidth', 2)
     plot(1951:2120, median(sc4.stall), 'Color', cmap(4,:), 'LineWidth', 1, 'LineStyle', ':')
     plot(1951:2120, median(sc6.stall), 'Color', cmap(4,:), 'LineWidth', 1, 'LineStyle', '-.')
     
     plot(1951:2120, median(sc7.stall), 'Color', cmap(5,:), 'LineWidth', 1)
     plot(1951:2120, median(sc7a.stall), 'Color', cmap(5,:), 'LineWidth', 1)
     plot(1951:2120, median(sc7b.stall), 'Color', cmap(5,:), 'LineWidth', 1)
     plot(1951:2120, median(sc8.stall), 'Color', cmap(6,:), 'LineWidth', 2)
     plot(1951:2120, median(sc9.stall), 'Color', cmap(5,:), 'LineWidth', 1)
     plot(1951:2120, median(sc10.stall), 'Color', cmap(5,:), 'LineWidth', 1)
     plot(1951:2120, median(sc11.stall), 'Color', cmap(5,:), 'LineWidth', 1)
  
  
  
     
      axis([2020 2120 0 80])
      
      yline(4)
      yline(10, '-.')
      
      legend('Sc0', 'Sc0b', 'S1', 'Sc2' , 'Sc7'  , 'Sc8')%, 'Sc1', 'Sc2')% 'Sc3 (Sc1+1)', 'Sc5 (Sc1+1)', 'Sc2', 'Sc4 (Sc2+1)',)
      

 

 title({'CC incidence, standardised'}, 'FontSize', 11); set(gca,'FontSize',11); ...
    ylabel({'per 100,000'}, 'FontSize', 11); xlabel({'Year'}, 'FontSize', 11); set(gca,'FontSize',11)

 filename2 = ['/Users/minttu/matlab_figures/Sc0_to_Sc6' ,'-', datestr(now,'yyyy-mm-dd'),'.png'] ; 
saveas(gcf, filename2)

%%
figure
% for i=1:length(sc0.stneg(:,1))
% hold on
%      plot(2005:2120, sc0.start(i,:), 'Color', cmap(1,:), 'LineWidth', mruns)
%      
%      plot(1986:2120, sc0.stnart(i,:), 'Color', cmap(2,:), 'LineWidth', mruns)
%     % plot(1950:2120, 100000.*sc0.cc_incid_neg, 'Color', col1, 'LineWidth', mruns)
% %hold off
% end

hold on
     plot(1951:2120, median(sc0.stall), 'Color', cmap(1,:), 'LineWidth', 1)
        plot(1951:2120, median(sc00.stall), 'Color', cmap(2,:), 'LineWidth', 2)
     
     plot(1951:2120, median(sc0b.stall),'Color',cmap(2,:), 'LineWidth', 1)
     plot(1951:2120, median(sc0c.stall),'Color',cmap(2,:), 'LineWidth', 1)
     
     plot(1951:2120, median(sc1.stall), 'Color', cmap(3,:), 'LineWidth', 2)
     plot(1951:2120, median(sc3.stall), 'Color', cmap(3,:), 'LineWidth', 1, 'LineStyle', ':')
     plot(1951:2120, median(sc5.stall), 'Color', cmap(3,:), 'LineWidth', 1, 'LineStyle', '-.')
           
     plot(1951:2120, median(sc2.stall), 'Color', cmap(4,:), 'LineWidth', 1)
     plot(1951:2120, median(sc4.stall), 'Color', cmap(4,:), 'LineWidth', 1, 'LineStyle', ':')
     plot(1951:2120, median(sc6.stall), 'Color', cmap(4,:), 'LineWidth', 1, 'LineStyle', '-.')
     
     plot(1951:2120, median(sc7.stall), 'Color', cmap(5,:), 'LineWidth', 1)
     plot(1951:2120, median(sc8.stall), 'Color', cmap(6,:), 'LineWidth', 1)
     plot(1951:2120, median(sc9.stall), 'Color', cmap(7,:), 'LineWidth', 1)
     plot(1951:2120, median(sc10.stall), 'Color', cmap(10,:), 'LineWidth', 1)
     plot(1951:2120, median(sc11.stall), 'Color', cmap(12,:), 'LineWidth', 1)
  
  
  
     
      axis([2016 2120 0 50])
      
      yline(4)
      yline(10, '-.')
      
      legend('Sc0', 'Sc2' , 'Sc7' , 'Sc8', 'Sc9' , 'Sc10' , 'Sc11')%, 'Sc1', 'Sc2')% 'Sc3 (Sc1+1)', 'Sc5 (Sc1+1)', 'Sc2', 'Sc4 (Sc2+1)',)
      
      % legend('Sc0', 'Sc00')

 

 title({'CC incidence, standardised'}, 'FontSize', 11); set(gca,'FontSize',11); ...
    ylabel({'per 100,000'}, 'FontSize', 11); xlabel({'Year'}, 'FontSize', 11); set(gca,'FontSize',11)

 filename2 = ['/Users/minttu/matlab_figures/Sc0_to_Sc11-new' ,'-', datestr(now,'yyyy-mm-dd'),'.png'] ; 
saveas(gcf, filename2)

%%
figure
% for i=1:length(sc0.stneg(:,1))
% hold on
%      plot(2005:2120, sc0.start(i,:), 'Color', cmap(1,:), 'LineWidth', mruns)
%      
%      plot(1986:2120, sc0.stnart(i,:), 'Color', cmap(2,:), 'LineWidth', mruns)
%     % plot(1950:2120, 100000.*sc0.cc_incid_neg, 'Color', col1, 'LineWidth', mruns)
% %hold off
% end

hold on
     plot(1951:2120, median(sc0.stall), 'Color', cmap(1,:), 'LineWidth', 2)
     
     plot(1951:2120, median(sc7.stall), 'Color', cmap(6,:), 'LineWidth', 1)   
     plot(1951:2120, median(sc8.stall), 'Color', cmap(7,:), 'LineWidth', 1)
     plot(1951:2120, median(sc9.stall), 'Color', cmap(8,:), 'LineWidth', 1, 'LineStyle', '-.')
     plot(1951:2120, median(sc10.stall), 'Color', cmap(11,:), 'LineWidth', 1, 'LineStyle', ':')
     plot(1951:2120, median(sc11.stall), 'Color', cmap(16,:), 'LineWidth', 1, 'LineStyle', ':')
  
  
  
     
      axis([2016 2120 0 35])
      
      yline(4)
      yline(10, '-.')
      
      legend('Sc0', 'S7', 'Sc8' , 'Sc9'  , 'Sc10' , 'Sc11')%, 'Sc1', 'Sc2')% 'Sc3 (Sc1+1)', 'Sc5 (Sc1+1)', 'Sc2', 'Sc4 (Sc2+1)',)
      

 

 title({'CC incidence, standardised'}, 'FontSize', 11); set(gca,'FontSize',11); ...
    ylabel({'per 100,000'}, 'FontSize', 11); xlabel({'Year'}, 'FontSize', 11); set(gca,'FontSize',11)

 filename2 = ['/Users/minttu/matlab_figures/Sc7_to_Sc11' ,'-', datestr(now,'yyyy-mm-dd'),'.png'] ; 
saveas(gcf, filename2)
