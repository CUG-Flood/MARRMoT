% Copyright (C) 2022 Dongdong Kong
%% 1. Prepare data
% Load the data
data_raw = readtable('input_2009_2012.csv');
input.time = datenum(data_raw.date, 'yyyy-mm-ddTHH:MM:SSZ');
input.streamflow = data_raw.R; %  input_data.Q_after / 3600 ; 
n = length(input.time);

%% Create a climatology data input structure. 
% NOTE: the names of all structure fields are hard-coded in each model
% file. These should not be changed.
input_climatology.precip  = data_raw.prcp_OBS; % Daily data: P rate  [mm/d]
input_climatology.temp    = data_raw.t2m; % Daily data: mean T  [degree C]
input_climatology.pet     = data_raw.PET; % Daily data: Ep rate [mm/d]
input_climatology.delta_t = 1/24;                         % time step size of the inputs: 1 [d]

%% Time periods for calibration.
warmup = 24*14; % 3162 input1 2006（157） 2007（3006）
% cal_idx = (warmup+1):(warmup+365*2);
cal_idx = (warmup + 1) : round(n*0.7);
eval_idx = max(cal_idx):n;

Q_obs = input.streamflow;

%% 5. Calibrate the model
% `optim_fun`: function to use for optimisation (must have same structure as fminsearch)
get_models_all;

n = length(models_all);
res = cell(n, 1);
for i = 1:n
    model = models_all{i};
    fprintf("[%02d]: running %s...\n", i, model)

    % model = 'm_10_susannah2_6p_2s';
    % model = 'h_28_xinanjiang_12p_4s';
    [m, mopts, optim_opts] = init_option(model);
    m.input_climate = input_climatology;
    
    tic
    [par_opt, of_cal, stopflag, output] =  m.calibrate(Q_obs, cal_idx, ...
          'my_cmaes',...      
          mopts.par_ini,...   
          optim_opts,...      % options to optim_fun
          mopts.of_name, ... 
          1,1,...             % should the OF be inversed?   Should I display details about the calibration?
          mopts.weights);     % additional arguments to of_name
    toc
    Q_sim = m.get_streamflow([],[],par_opt);  
    of_eval = feval(mopts.of_name, Q_obs, Q_sim, eval_idx, mopts.weights);
    
    res{i} = struct('name', model, 'of_cal', of_cal, 'of_eval', of_eval, ...
        'Q_sim', Q_sim);
end


%% 6. Evaluate the calibrated parameters on unseen data
% Run the model with calibrated parameters, get only the streamflow

%% 7. Visualise the results
plot_runoff(input.time, Q_obs, Q_sim, warmup, cal_idx, eval_idx)
% writefig('Plot_ODE.pdf', 12, 6)
