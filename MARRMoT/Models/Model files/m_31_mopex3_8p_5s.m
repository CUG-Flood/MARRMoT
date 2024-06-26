classdef m_31_mopex3_8p_5s < MARRMoT_model
% Class for hydrologic conceptual model: MOPEX-3

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Model reference
% Ye, S., Yaeger, M., Coopersmith, E., Cheng, L., & Sivapalan, M. (2012). 
% Exploring the physical controls of regional patterns of flow duration 
% curves - Part 2: Role of seasonality, the regime curve, and associated 
% process controls. Hydrology and Earth System Sciences, 16(11), 4447�4465.
% http://doi.org/10.5194/hess-16-4447-2012

    properties
        % model-specific attributes
    end
    methods
        
        % creator method
        function obj = m_31_mopex3_8p_5s()
            obj.numStores = 5;                                             % number of model stores
            obj.numFluxes = 11;                                            % number of model fluxes
            obj.numParams = 8;
            
            obj.JacobPattern  = [1,0,0,0,0;
                                 1,1,0,0,0;
                                 0,1,1,0,0;
                                 0,1,1,1,0;
                                 0,0,1,0,1];                               % Jacobian matrix of model store ODEs
                             
            obj.parRanges = [-3  , 3;       % tcrit, Snowfall & snowmelt temperature [oC]
                             0   , 20;      % ddf, Degree-day factor for snowmelt [mm/oC/d]
                             1   , 2000;    % Sb1, Maximum soil moisture storage [mm]
                             0   , 1 ;      % tw, Groundwater leakage time [d-1]
                             0   , 1 ;      % tu, Slow flow routing response time [d-1]
                             0.05, 0.95;    % se, Root zone storage capacity as fraction of Sb2 [-]
                             1   , 2000;    % Sb2, Root zone storage capacity [mm]
                             0   , 1];      % tc, Mean residence time [d-1]
            
            obj.StoreNames = {"S1", "S2" "S3" "S4" "S5"};                   % Names for the stores
            obj.FluxNames  = {"ps", "pr",  "qn",  "et1", "q1f",...
                              "qw", "et2", "q2f", "q2u", "qf", "qs"};      % Names for the fluxes
            
            obj.FluxGroups.Ea = [4 7];                                     % Index or indices of fluxes to add to Actual ET
            obj.FluxGroups.Q  = [10 11];                                   % Index or indices of fluxes to add to Streamflow

        end
        
        % INITialisation function
        function obj = init(obj)          
        end
        
        % MODEL_FUN are the model governing equations in state-space formulation
        function [dS, fluxes] = model_fun(obj, S)
            % parameters
            theta   = obj.theta;
            tcrit = theta(1);     % Snowfall & snowmelt temperature [oC]
            ddf   = theta(2);     % Degree-day factor for snowmelt [mm/oC/d]
            s2max = theta(3);     % Maximum soil moisture storage [mm]
            tw    = theta(4);     % Groundwater leakage time [d-1]
            tu    = theta(5);     % Slow flow routing response time [d-1]
            se    = theta(6);     % Root zone storage capacity as fraction of s3max [-]
            s3max = theta(7);     % Maximum groundwater storage [mm]
            tc    = theta(8);     % Mean residence time [d-1]
            
            % delta_t
            delta_t = obj.delta_t;
            
            % stores
            S1 = S(1);
            S2 = S(2);
            S3 = S(3);
            S4 = S(4);
            S5 = S(5);
            
            % climate input
            t = obj.t;                             % this time step
            climate_in = obj.input_climate(t,:);   % climate at this step
            P  = climate_in(1);
            Ep = climate_in(2);
            T  = climate_in(3);
            
            % fluxes functions
            flux_ps   = snowfall_1(P,T,tcrit);
            flux_pr   = rainfall_1(P,T,tcrit);
            flux_qn   = melt_1(ddf,tcrit,T,S1,delta_t);
            flux_et1  = evap_7(S2,s2max,Ep,delta_t);
            flux_q1f  = saturation_1(flux_pr+flux_qn,S2,s2max);
            flux_qw   = recharge_3(tw,S2);
            flux_et2  = evap_7(S3,se*s3max,Ep,delta_t);
            flux_q2f  = saturation_1(flux_qw,S3,s3max);
            flux_q2u  = baseflow_1(tu,S3);
            flux_qf   = baseflow_1(tc,S4);
            flux_qs   = baseflow_1(tc,S5);

            % stores ODEs
            dS1 = flux_ps - flux_qn;
            dS2 = flux_pr + flux_qn - flux_et1 - flux_q1f - flux_qw;
            dS3 = flux_qw - flux_et2 - flux_q2f - flux_q2u;
            dS4 = flux_q1f + flux_q2f - flux_qf;
            dS5 = flux_q2u - flux_qs;
            
            
            % outputs
            dS = [dS1 dS2 dS3 dS4 dS5];
            fluxes = [flux_ps, flux_pr, flux_qn, flux_et1,...
                      flux_q1f, flux_qw, flux_et2, flux_q2f,...
                      flux_q2u, flux_qf, flux_qs];
        end
        
        % STEP runs at the end of every timestep
        function obj = step(obj)
        end
    end
end