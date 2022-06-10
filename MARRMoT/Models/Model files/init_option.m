function [m, mopts, defopts] = init_option(model)
    % model = 'm_10_susannah2_6p_2s';
    m = feval(model);
    parRanges = m.parRanges; % Parameter ranges
    numParams = m.numParams; % Number of parameters
    % numStores = m.numStores; % Number of stores
    m.S0 = zeros(m.numStores, 1); % Initial storages (see note in paragraph 5 on model warm-up)
    
    % initial parameter set
    mopts.par_ini = mean(parRanges, 2); % same as default value
    % Choose the objective function
    mopts.of_name = 'of_KGE'; % This function is provided as part of MARRMoT. See ./MARRMoT/Functions/Objective functions
    mopts.weights = [1, 1, 1]; % Weights for the three KGE components

    %% 3. Define the solver settings
    % Create a solver settings data input structure.
    % NOTE: the names of all structure fields are hard-coded in each model
    % file. These should not be changed.
    input_solver_opts.resnorm_tolerance = 0.1; % Root-finding convergence tolerance;
    % users have reported differences in simulation accuracy (KGE scores) during calibration between Matlab and Octave for a given tolerance.
    % In certain cases, Octave seems to require tigther tolerances to obtain the same KGE scores as Matlab does.
    input_solver_opts.resnorm_maxiter = 6; % Maximum number of re-runs
    m.solver_opts = input_solver_opts;

    %% 4. Define calibration settings
    % Settings for 'my_cmaes'
    % the opts struct is made up of two fields, the names are hardcoded, so
    % they cannot be changed:
    %    .sigma0:     initial value of sigma
    %    .cmaes_opts: struct of options for cmaes, see cmaes documentation
    %                 or type cmaes to see list of options and default values

    % starting sigma
    defopts.insigma = .3 * (parRanges(:, 2) - parRanges(:, 1)); % starting sigma (this is default, could have left it blank)
    defopts.LBounds = parRanges(:, 1); % lower bounds of parameters
    defopts.UBounds = parRanges(:, 2); % upper bounds of parameters
    defopts.PopSize = 4 + floor(3 * log(numParams)); % population size (default)
    defopts.TolX = 1e-6 * min(defopts.insigma); % stopping criterion on changes to parameters
    defopts.TolFun = 1e-4; % stopping criterion on changes to fitness function
    defopts.TolHistFun = 1e-5; % stopping criterion on changes to fitness function
    % defopts.SaveFilename = 'wf_ex_4_cmaesvars.mat'; % output file of cmaes variables
    % defopts.LogFilenamePrefix = 'wf_ex_4_'; % prefix for cmaes log-files
    
    defopts.MaxFunEvals = 1500;
    defopts.MaxIter = 200;
    
    % note that saving to ".mat" file of CMA-ES output is disabled for Octave
    % (lines 1795 and 1839 of my_cmaes.m) due to a bug that otherwise crashes the calibration.

    % Other useful options
    % change to true to run in parallel on a pool of CPUs (e.g. on a cluster)
    defopts.EvalParallel = false;
    % uncomment to restart r-times until the condition
    % r = 2;
    % defopts.Restarts  = r;
    % defopts.RestartIf = 'fmin > -.8'; % OF is inverted, so this restarts
    %                                        unless the OF (KGE here) > 0.8
    % first set up the model
    % m.input_climate = input_climatology;
    %m.delta_t       = input_climatology.delta_t;         % unnecessary if input_climate already contains .delta_t
    % debugging options
    %defopts.MaxIter = 5;                                                    % just do 5 iterations, to check if it works
    %defopts.Seed = 1234;                                                    % for reproducibility
end
