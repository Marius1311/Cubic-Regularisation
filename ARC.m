function [x, k, iterates, RES, gradients] = ARC( f, g, H, x0, varargin )
% Adaptive Regularisation Using Cubics for Optimisation
%
% Input:
%   f: the objective function
%   g: a function handle to the gradient, returning column vectors
%   x0: initial guess
%
% Optional Input: 
%   options: a struct containing any of the following fields:
%       epsilon:  target precision
%       IterMax:   max number of iterations
%       sigma0: initial regularisation parameter
%       sigmaMin: minimu value for the regularisation parameter
%       theta: relative stopping criterion model minimisation
%       eta1: smallest decrease ratio we are happy with
%       eta2: definition very sucessful step
%       gamma1: decrease in the regularisation parameter
%       gamma2: small increase in the regularisation paramter
%       gamma3: big increase in the regularisation parameter 
%       outputLevel: 0 for nop output, 1 for a little output, 2 for full
%           output
%
% Output:
%       x: solution
%       k: final number of iterations
%       gradients: vector that contains the 2-norm of the succesive gradients
%       RES: final gradient 2-norm
%       iterates: matrix where each row corresponds to one iterate

%% Parameters

if nargin == 4
    % use standart settings
    epsilon = 1.e-5; IterMax = 100; sigma0 = 1; theta = 1e-4; eta1 = 0.1; eta2 = 0.9; outputLevel = 2;
else 
    % use whatever is in the options struct
    options = varargin{1};
   [epsilon, IterMax, sigma0, theta, eta1, eta2, outputLevel] = ARC_Input(options); 
end

%% Initialise

k = 0; % Initialise outer loop count
n = length(x0); % gives the dimension of the objective function
iterates = zeros(IterMax, n); % An array to store the iterates
CurrentX = x0; % use initial guess
iterates(1, :) = CurrentX;
gradients = zeros(IterMax, 1); % We will store the magnitude of the gradient in here, not the actual gradients
grad_f = g(CurrentX);
hess_f = H(CurrentX);
gradients(1, 1) = norm(grad_f, 2); % Check norm of gradient at initial guess
RES = gradients(1, 1);
sigma = sigma0;
options_GLRT = struct('theta', theta, 'outputLevel', outputLevel);

%% Initial Information
if outputLevel > 0
    fprintf('------------------------------------------------\n');
    fprintf('ARC for smooth non-convex optimisation.');
    fprintf( [ '  \n  iter       f    ', ...
        '        ||g||    sigma    ratio \n' ] );
    
    fprintf( '%6.0f %13.6e %12.6e %8.2e %8.2e \n', ...
        0, f(CurrentX), RES, sigma, 0 );
    fprintf('------------------------------------------------\n');
end

%% Outer Loop

while k <= IterMax && RES > epsilon
    if outputLevel > 0
        fprintf('Starting iteration %1.0f \n', k+1);
    end
    
    %% Step 1: Find search direction s
    [s, stat] = GLRT(grad_f, hess_f, sigma, options_GLRT); % Calculates a step that gives sufficient decrease
    s = s';
    % We should check whether this really minimiser the local model in the
    % current subspace:
    
    if stat == 0
        break;
    end
%     disp(CurrentX);
%     disp(grad_f);
%     disp(hess_f);
%     disp(s);
    m = f(CurrentX) + s*grad_f + 0.5*s*hess_f*s' + 1/3* sigma * norm(s, 2)^3; % calculates the value of the model
    
    %% Step 2: Check whether this was a sucessful step
     rho = (  f(CurrentX ) - f(CurrentX + s) ) / ...
         ( f(CurrentX ) - m);  
     
     %% Step 3: Update x
     
     if rho >= eta1
         CurrentX = CurrentX + s;
     end
         
    %% Step 4: Save new iterate and new gradient norm
    iterates(k+1, :) = CurrentX;
    grad_f = g(CurrentX);
    hess_f = H(CurrentX);
    gradients(k+1, 1) = norm(grad_f, 2); 
    RES = gradients(k+1, 1); % Calculate residual
    k = k+1;
    
    %% Step 5: Update sigma
     if rho > eta2
         sigma = max([min([sigma, norm(grad_f, 2)]), 1e-16]);
     end
     if rho < eta1
         sigma = sigma * 2;
     end
    
    %% Step 6: Display information
    
    if outputLevel > 0
    fprintf( [ '  \n  iter       f    ', ...
                         '        ||g||    sigma    ratio \n' ] ); 
                     
    fprintf( '%6.0f %13.6e %12.6e %8.2e %8.2e \n', ...
                      k, f(CurrentX), RES, sigma, rho );
    fprintf('------------------------------------------------\n');
    end
    
end


%% Return

if stat == 0   
fprintf('------------------------------------------------\n');
fprintf('Could not minimize the local model. Failed in');
    fprintf( [ '  \n  iter       f    ', ...
                         '        ||g||  \n' ] ); 
                     
    fprintf( '%6.0f %13.6e %12.6e %8.2e %8.2e', ...
                      k, f(CurrentX), RES );
    fprintf('\n \n');
    disp('Current x value = ');
    disp(CurrentX);
fprintf(' ------------------------------------------------\n');
fprintf(' ------------------------------------------------\n');
else
fprintf('------------------------------------------------\n');
fprintf('Sucess! Iteration Converged.');
    fprintf( [ '  \n  iter       f    ', ...
                         '        ||g||  \n' ] ); 
                     
    fprintf( '%6.0f %13.6e %12.6e %8.2e %8.2e \n', ...
                      k, f(CurrentX), RES );
fprintf(' \n ------------------------------------------------\n');
fprintf(' ------------------------------------------------\n');
end

x = CurrentX;


end

