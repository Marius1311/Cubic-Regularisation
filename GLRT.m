function [ s, stat ] = GLRT(grad_f, hess_f, sigma, varargin)
% Lanczos Type Method for finding the global minimum of the regularised
% local quadratic model m(s) in a subspace S_j
%
% Input -
%   f: objective function
%   CurrentX: current iterate
%   grad_f: current gradient
%   hess_f: current hessian
%
% Optional Input - 
%   theta: algorithm parameter
%   sigma0: initial regularisation
%   outputLevel: 0, 1 or 2
%
% Output -
%   s, global minimiser over a subspace
%   stat: algortihm exit status

%% Process the Input

if nargin == 3
    % use standart settings
    theta = 1e-4; outputLevel = 2;
else
    options = varargin{1};
    [theta, outputLevel] = GLRT_Input(options);
end
    
%% Initialise

j = 1; % Loop counter
n = length(grad_f); % problem dimension
gamma0 = norm(grad_f);
Lanczos_Tol = 1e-10;
Lanczos_Options = struct('Lanczos_Tol', Lanczos_Tol);
Newton_Options = struct('Epsilon_Newton', 1e-6, 'IterMax_Newton', n*20, 'OutputLevel', outputLevel);

% Initialise the Algortihm 
q_old_old = 0; q_old = grad_f/gamma0; beta_old = 0;
[ alpha, beta_new, q_new ] = Lanczos_Algorithm( q_old, q_old_old, hess_f, beta_old);
V = [q_old, q_new]; T = [alpha; beta_new];

% Compute the first multiplier
lambda = - alpha/2 + sqrt( (alpha/2)^2 + sigma*gamma0);
u = - lambda/sigma;

% Define a function. Shall accept row vectors as input
% m = @(s) s*grad_f + 1/2*s*hess_f*s' + 1/3*sigma*norm(s)^3;
% grad_m = @(s) (hess_f + sigma*norm(s) * eye(n))*s + grad_f;

%% Give out information
if outputLevel > 0
    fprintf('\n \t \t GLRT: sigma = %1.2f. \n ', sigma);
    fprintf( [ '  \n \t \t   iter       lambda    ', ...
        '  beta_k * |u_k^(k)|    theta*||u||^2   \t ||u|| \n' ] );
    
    fprintf( ' \t \t %6.0f %13.3e \t %12.3e \t \t \t %8.3e \t \t   %8.3e \n', ...
        j, lambda, beta_new * abs(u(end)), theta*norm(u)^2, norm(u) );
    
end

%% Inner Loop over expanding subspaces

%while norm(gradM) > min([theta, sqrt(norm(grad_f))])*norm(grad_f) && j <= n
while beta_new * abs(u(end)) > theta*norm(u)^2 && j <= n
    %% Create tridiagonal T and basis Q
    
    j = j+1;
    
    %[ Q, T ] = Lanczos( j, n, Q, T, CurrentX, hess_f);
    
    % Move values around a little
    q_old_old = q_old;
    q_old = q_new;
    beta_old = beta_new;
    
    % Create new quantities
    [ alpha, beta_new, q_new ] = Lanczos_Algorithm( q_old, q_old_old, hess_f, beta_old, Lanczos_Options);
    
    % adjust the matrices
    T = [[T; zeros(1, j-1)], [zeros(j-2, 1); beta_old; alpha; beta_new]];
    
    % Check for a breakdown in the Lanczos algortihm
    if beta_new < Lanczos_Tol
        disp('there was a breakdown');
    else
        V = [V, q_new];
    end
    
    
    %% Solve for u in the current subspace
    
    % We need to do an efficient cholesky factorisation and we need to
    % estimate the smallest eigenvalue of T.
    
    % Compute u by a Newton Type Method
    [ u, lambda, phi ] = AdaptedNewton(n, T(1:j, 1:j), sigma, gamma0, Newton_Options) ;
      
    % Does this u really minimise the local model over the current
    % subspace?
     % Create the identity
     
%     if j ~= n
%         I = eye(j+1);
%     else
%         I = eye(j);
%     end
    I = eye(j);
    error = norm((T(1:j, 1:j) + lambda*I)*u + gamma0*I(:, 1));
    
    % Also check positive semi definiteness:
    smallestEV = eigs(T(1:j, 1:j) + lambda*I, 1, 'sa');
    
    % Result: u does indeed minimise the local model globaly in the current
    % subspace. Next, we can try to check whether the Matrices Q and T do
    % as they should:
    MatrixError = norm(V(:, 1:j)'*hess_f*V(:, 1:j) - T(1:j, 1:j));
    %disp(MatrixError);
    
    %% Update the model gradient and the loop counter
    
    % gradM = grad_f + hess_f*s + lambda * s; % Gradient of the local model  
    
    
    %% Display information
    
    if outputLevel > 0
        fprintf( [ '  \n \t \t   iter       lambda    ', ...
            '  beta_k * |u_k^(k)|    theta*||u||^2   phi \t ||u|| \t \t Error \t \t smallestEV \t MatErr \n' ] );
        
        fprintf( ' \t \t %6.0f %13.3e %12.3e \t \t \t \t %8.3e %8.3e \t %1.2e \t %1.4e \t %1.3f \t \t %1.3f \n', ...
            j-1, lambda, beta_new * abs(u(end)), theta*norm(u)^2, phi, norm(u), error, smallestEV, MatrixError );
    end
    % We will stop the iteration as soon as we either drive the norm of the gradient under the
    % relative stopping criterion or we hit the iteration limit n, which is the problem dimension           
    
end

% Recover s from here:
s = V(:, 1:j)*u;

if j >= n && beta_new * abs(u(end)) > theta*norm(u)^2
    stat = 0; % Failed to converge while trying to minimise the local model
else
    stat = 1; % everything is fine
end

end

