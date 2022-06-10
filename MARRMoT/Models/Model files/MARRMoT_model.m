classdef MARRMoT_model < MARRMoT_model0
    % Superclass for all MARRMoT models

    % Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
    % This file is part of the Modular Assessment of Rainfall-Runoff Models
    % Toolbox (MARRMoT).
    % MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
    % WARRANTY. See <https://www.gnu.org/licenses/> for details.

    properties
        % attribute to store whether we are running MATLAB or Octave
        % isOctave % 1 if we're on Octave, 0 if MATLAB
        % % static attributes, set for each models in the model definition
        % numStores % number of model stores
        % numFluxes % number of model fluxes
        % numParams % number of model parameters
        % parRanges % default parameter ranges
        % JacobPattern % pattern of the Jacobian matrix of model store ODEs
        % StoreNames % Names for the stores
        % FluxNames % Names for the fluxes
        % FluxGroups % Grouping of fluxes (useful for water balance and output)
        % StoreSigns % Signs to give to stores (-1 is a deficit store), assumes all 1 if not given
        % % attributes set at the beginning of the simulation
        % % directly by the user
        % theta % Set of parameters
        % delta_t % time step
        % S0 % initial store values
        % input_climate % vector of input climate
        % solver_opts % options for numerical solving of ODEs
        % % automatically, based on parameter set
        % store_min % store minimum values
        % store_max % store maximum values
        % % attributes created and updated automatically throughout a
        % % simulation
        % t % current timestep
        % fluxes % vector of all fluxes
        % stores % vector of all stores
        % uhs % unit hydrographs and still-to-flow fluxes
        % solver_data % step-by-step info of solver used and residuals
        % status % 0 = model created, 1 = simulation ended

        % add two paramter
        use_sceua
        use_ode
    end

    methods

        function obj = init_(obj)
            % min and max of stores
            obj.store_min = zeros(obj.numStores, 1);
            obj.store_max = inf(obj.numStores, 1);

            % empty vectors of fluxes and stores
            t_end = size(obj.input_climate, 1);
            obj.stores = zeros(t_end, obj.numStores);
            obj.fluxes = zeros(t_end, obj.numFluxes);

            % empty struct with the solver data
            obj.solver_data.resnorm = zeros(t_end, 1);
            obj.solver_data.solver = zeros(t_end, 1);
            if (~obj.isOctave); obj.solver_data.solver = categorical(obj.solver_data.solver); end;
            obj.solver_data.iter = zeros(t_end, 1);

            % model specific initialisation
            obj.init();
        end

        function [] = run(obj, input_climate, S0, theta, solver_opts)

            if nargin > 4 && ~isempty(solver_opts)
                obj.solver_opts = solver_opts;
            end

            if nargin > 3 && ~isempty(theta)
                obj.theta = theta;
            end

            if nargin > 2 && ~isempty(S0)
                obj.S0 = S0;
            end

            if nargin > 1 && ~isempty(input_climate)
                obj.input_climate = input_climate;
            end

            % run INIT_ method, this will calculate all auxiliary parameters
            % and set up routing vectors and store limits
            obj.init_();

            t_end = size(obj.input_climate, 1);

            resnorm = -999;
            solver = -999;
            iter = 999;

            for t = 1:t_end
                obj.t = t;

                if t == 1
                    Sold = obj.S0(:);
                else
                    Sold = obj.stores(t - 1, :)';
                end

                %% rm the ODE solve part;
                % Sold: previous time
                if obj.use_ode
                    [Snew, resnorm, solver, iter] = obj.solve_stores(Sold);
                    obj.solver_data.resnorm(t) = resnorm;
                    obj.solver_data.solver(t) = solver;
                    obj.solver_data.iter(t) = iter;
                    [dS, f] = obj.model_fun(Snew);
                else
                    %% 一个循环的关系
                    % dS: W, S, Si, Sg
                    % [dS, f] = obj.model_fun(Snew);
                    % Snew = Sold + dS' * obj.delta_t;
                    [dS, f] = obj.model_fun(Sold);
                end

                obj.fluxes(t, :) = f * obj.delta_t;
                obj.stores(t, :) = Sold + dS' * obj.delta_t;
                obj.step();
            end

            obj.status = 1;
        end

        % CALIBRATE uses the chosen algorithm to find the optimal parameter
        % set, given model inputs, objective function and observed streamflow.
        % the function chosen in algorithm should have the same inputs and
        % outputs as MATLAB's fminsearch.

        %% Arguments
        % ---
        %
        %% Return
        % ---
        % `par_opt`  : optimal parameter set
        % `of_cal`   : value of objective function at par_opt
        % `stopflag` : flag indicating reason the algorithm stopped
        % `output`   : see fminsearch for detail
        function [par_opt, of_cal, stopflag, output] = calibrate(obj, ...
            Q_obs, ... % observed streamflow
                cal_idx, ... % timesteps to use for model calibration
                optim_fun, ... % function to use for optimisation (must have same structure as fminsearch)
                par_ini, ... % initial parameter estimates
                optim_opts, ... % options to optim_fun
                of_name, ... % name of objective function to use
                inverse_flag, ... % should the OF be inversed?
                display, ... % should I display information about the calibration?
                varargin) % additional arguments to the objective function

            if isempty(obj.input_climate) || isempty(obj.delta_t) || ...
                    isempty(obj.S0) || isempty(obj.solver_opts)
                error(['input_climate, delta_t, S0 and solver_opts ' ...
                        'attributes must be specified before calling ' ...
                    'calibrate.']);
            end
            % obj.init_();
            
            % if the list of timesteps to use for calibration is empty,
            % use all steps
            if isempty(cal_idx)
                cal_idx = 1:length(Q_obs);
            end

            % use display by default
            if isempty(display)
                display = true;
            end

            % use the data from the start to the last value of cal_idx to
            % run the simulation
            if islogical(cal_idx); cal_idx = find(cal_idx); end
            input_climate_all = obj.input_climate;
            obj.input_climate = input_climate_all(1:max(cal_idx), :);
            Q_obs = Q_obs(1:max(cal_idx));

            % if the initial parameter set isn't set,  start from mean
            % values of parameter range
            if isempty(par_ini)
                par_ini = mean(obj.parRanges, 2);
            end

            % helper function to calculate fitness given a set of
            % parameters
            function fitness = fitness_fun(par)
                Q_sim = obj.get_streamflow([], [], par);
                fitness = (-1)^inverse_flag * feval(of_name, Q_obs, Q_sim, cal_idx, varargin{:});
            end

            % display some useful things for the user to make sure they
            % used the right settings
            if display
                disp('---')
                disp(['Starting calibration of model ' class(obj) '.'])
                disp(['Simulation will run for timesteps 1-' num2str(max(cal_idx)) '.'])

                % this is a bit ugly, but it formats the list of cal_idx in
                % a pretty and concise way
                cal_idx = sort(cal_idx);
                i = 1;
                previous = cal_idx(i);
                cal_idx_str = num2str(previous);

                while i < numel(cal_idx)
                    i = i + 1;

                    if cal_idx(i) - previous == 1
                        i = find(diff(cal_idx(i:end)) ~= 1, 1) + i - 1;
                        if isempty(i); i = numel(cal_idx); end
                        previous = cal_idx(i);
                        cal_idx_str = append(cal_idx_str, '-', num2str(previous));
                    else
                        previous = cal_idx(i);
                        cal_idx_str = append(cal_idx_str, ', ', num2str(previous));
                    end

                end

                disp(['Objective function ' of_name ' will be calculated in time steps ' cal_idx_str '.'])
                disp(['The optimiser ' optim_fun ' will be used to optimise the objective function.'])
                disp(['Options passed to the optimiser:'])
                disp(optim_opts)
                disp('All other options are left to their default,')
                disp('check the source code of the optimiser to find these default values.')
                disp('---')
            end
            
            %% Update by Dongdong Kong
            obj.use_sceua = false;
            obj.use_ode = false;
            
            %% MarrMot default option
            % obj.use_sceua = false;
            % obj.use_ode = true;
            
            % 自带的优化算法效率已经比较高了
            % rng(-1);
            rand('seed', -1);
            if obj.use_sceua
                bl = obj.parRanges(:, 1);
                bu = obj.parRanges(:, 2);
                [par_opt, of_cal, stopflag, output] = sceua(@fitness_fun, par_ini, bl, bu, optim_opts.MaxFunEvals);
            else
                % see fminsearch for detail
                [par_opt, of_cal, stopflag, output] = feval(optim_fun, @fitness_fun, par_ini, optim_opts);
            end

            % if of_cal was inverted, invert it back before returning
            of_cal = (-1)^inverse_flag * of_cal;

            % reset the whole input climate as it was before the
            % calibration
            obj.input_climate = input_climate_all;
        end

    end

end
