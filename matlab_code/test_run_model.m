

addpath('/Users/minttu/det_hiv_hpv/matlab_code/') 

 
%% C code set up

vaccscen.setup = 0;
screenscen.setup = 0;

% this is here to generate a unique c_output for each model run
% as multiple will be run in the same environment

% this is for cluster
run_name = string(randi([1 1000000]));

% calibrated parameters, divided into behavior, demography, hiv, and hpv
filename = ['calibparams-3.mat'];
load(filename)

%% HERE NEED TO CHANGE WHEN IN CLUSTER
% Setting up the format in which parameters are read when executing the
% model

pop1950 = load('pop1950.mat');

c_command = strcat('c_code', filesep, 'hpv_hiv'); %% N.B. This path needs to lead to the executable 

c_threads = '12'; % this needs to be the same number as -n on the slurm
c_params = 'matlab_code/input/params_age13_0.bin';
c_pop = 'matlab_code/input/pop_13ages.bin';
c_scenario = strcat("matlab_code/input/scenario_", "0", ".bin");
c_output = strcat("outputc_", run_name, ".bin");

write_scenario_c(vaccscen, screenscen, c_scenario);
write_popfile_c(pop1950.pop1950, c_pop);

for ir=1:length(dpg(1,:))
          
    dp = dpg(:, ir);
    bp = bpg(:, ir);
    hp = hpg(:, ir);
    pp = ppg(:, ir);
    
    write_params_c(bp, dp, hp, pp, c_params);

    %% Run the C code - 9 mins on my comp
   % tic
   % if it gives an error, check the file names lead to right locations for
   % these
    system(strcat(c_command," ", c_params," ",c_pop," ", c_scenario," ", c_output," ",c_threads));
   % toc/60
    
    ir

    res{ir} = read_c_output(c_output);
    

end

     filename = ['test-', datestr(now,'yyyy-mm-dd_HH:MM'), '.mat'];  
      save(filename, 'res', 'dpg', 'bpg', 'hpg', 'ppg')
