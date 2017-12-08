function [ u, lambda, phi ] = AdaptedNewton(n,  T, sigma, gamma0, varargin )                      
% Adapted Newton
%   A Newton-Type method for rootfinding of a non-linear function
%
%   Input -
%   n: Problem dimension
%   T: Tridiagonal reduction of the Hessian
%   sigma: Regularisation parameter
%   gamma0: Norm of the gradient
% 
%   Optional Input - 
%   Epsilon_Newton: Tolerance for the algorithm
%   IterMax_Newton: Iteration Limit
%   OutputLevel: 0, 1 or 2
% 
%   Output - 
%   u: Solution in the current subspace
%   lambda: Corresponding Multiplier
%   phi: Corresponding function value

%% Process the Input

if nargin == 4
    % use standart settings
    Epsilon_Newton = 1e-6; IterMax_Newton = n*20; outputLevel = 2;
else
    % use custom settings
    options = varargin{1};
    
    if isfield(options, 'Epsilon_Newton')
        Epsilon_Newton = options.Epsilon_Newton;
    else
        Epsilon_Newton = 1e-6;
    end
    
    if isfield(options, 'IterMax_Newton')
        IterMax_Newton = options.IterMax_Newton;
    else
        IterMax_Newton = n*20;
    end
    
    if isfield(options, 'OutputLevel')
        outputLevel = options.OutputLevel;
    else
        outputLevel = 2;
    end
end
    
%% Initialise

    % Create the identity
    [j, ~] = size(T);
    I = eye(j);
    
%     if j ~= n
%         I = eye(j+1);
%     else
%         I = eye(j);
%     end
     
      % Initialise lambda. This is not yet perfect, as it could be more
      % efficient by making use of the previous value of lambda, by testing
      % whether or not it works. This however requires an efficient way to
      % get the nex factorisation from the old one.
      
      eSmallest = eigs(T, 1, 'sa');
      if eSmallest > 0
          epsLambda = (eSmallest + 1)*sqrt(1e-16);
          lambda = epsLambda;
      else    
      epsLambda = (-eSmallest + 1)*sqrt(1e-16);
      lambda =-eSmallest + epsLambda;
      end
      
      flag = 0;
      i = 1;
         
     while flag == 0   
     try
         L = chol(T + lambda * I,'lower'); 
         flag = 1;
     catch
         lambda = lambda + i*epsLambda;
         flag = 0;
         i = i*2;
     end
     end
     
     % We actually also need to check whether phi for this lambda is
     % smaller than zero. Only in this case we can global convergence
          
    % Initialise value of phi:
    phi = 1;
    
    % Loop counter
    l = 1;
    
    %% Display information
     
    if outputLevel > 0
        fprintf('\n \t \t \t \t Newton: lambda0 = %1.3f \n ', lambda);
    end
    
%% Main Iteration
    
    while abs(phi) > Epsilon_Newton && l <= IterMax_Newton
        
        if l ~= 1
            % Step 1: Factorise T + lambda I
            L = chol(T + lambda * I,'lower');
        end
    
    % Step 2: Solve L L^T s = -g
    u = (L*L')\(-gamma0*I(:, 1));
    
    % Step 3: Solve Lw = s;
    w = L\u;
    
    % Step 4: Compute the newton correction
    DeltaLambda = lambda* ( norm(u)-lambda/sigma ) / ( norm(u)+lambda/sigma*(lambda*norm(w)^2/norm(u)^2) );
    
    % Step 5: Update lambda
     lambda = lambda + DeltaLambda;
    
    % Compute value of phi:
    phi = 1/norm(u) - sigma/lambda;
    
   %% Display information
    
   if outputLevel > 1
       fprintf(  ' \n \t \t \t \t  iter      lambda    DeltaLambda       phi      \n' );
       
       fprintf( ' \t \t \t \t %6.0f %13.3e %12.3e %12.3e \n', ...
           l, lambda, DeltaLambda, phi );
   end
                  
    %%
    
    % Update loop counter
    l = l+1;
    
    end

end

