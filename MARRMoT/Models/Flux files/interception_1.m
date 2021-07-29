function [out] = interception_1(In,S,Smax,varargin)
%interception_1 Creates function for store overflow: uses logistic smoother.
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Flux function
% ------------------
% Description:  Interception excess when maximum capacity is reached
% Constraints:  -
% @(Inputs):    In   - incoming flux [mm/d]
%               S    - current storage [mm]
%               Smax - maximum storage [mm]
%               varargin(1) - smoothing variable r (default 0.01)
%               varargin(2) - smoothing variable e (default 5.00)

if size(varargin,2) == 0
    out = In.*(1-smoothThreshold_storage_logistic(S,Smax));
elseif size(varargin,2) == 1
    out = In.*(1-smoothThreshold_storage_logistic(S,Smax,varargin(1)));
elseif size(varargin,2) == 2
    out = In.*(1-smoothThreshold_storage_logistic(S,Smax,varargin(1),varargin(2)));    
end

end