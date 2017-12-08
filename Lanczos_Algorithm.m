function [ alpha, beta_new, q_new ] = Lanczos_Algorithm( q_old, q_old_old, hess_f, beta_old, varargin )
% Iteration k of the Lanczos Algorithm.
%   Input - 
%       q_old: the last lanczos vector
%       q_old_old: the one before that
%       hess_f: current hessian matrix of the objective
%       beta_old: last beta value
% 
%   Output - 
%       alpha: the nex diagonal entry
%       beta_new: the new off diagonal entry
%       q_new: the new lanczos vector

%% Input handeling

if nargin == 4
    Lanczos_Tol = 1e-10;
else 
    options = varargin{1};
    if isfield(options, 'Lanczos_Tol')
        Lanczos_Tol = options.Lanczos_Tol;
    else
        Lanczos_Tol = 1e-10;
    end
end

%% Carry out one Lanczos Iteration

% compute the new quantities
w = hess_f*q_old;
alpha = q_old'*w;
w = w-beta_old*q_old_old- alpha*q_old;
beta_new = norm(w);

% check for a breakdown
if beta_new < Lanczos_Tol
    q_new = NaN;
else
    q_new = w/beta_new;
end

end

