%% This is a basic implementation of the ARC optimisation algorithm.

%% State the problem, initialise

% Create an objective function, gradient handle and hessian handle as
% symbolic objects
syms  a b c;
f_sym = symfun(a^4 + b^2 + c^6, [ a b c]);
g_sym = gradient(f_sym);
H_sym = hessian(f_sym);

% Convert to matlab function handles
f = matlabFunction(f_sym,'Vars',{[ a b c]});
g = matlabFunction(g_sym,'Vars',{[ a b c]});
H = matlabFunction(H_sym,'Vars',{[ a b c]});

% Initial guess
x0 = [1 2 3];    

%% If the input argument is two dimensional:

if length(x0) == 2
xmin = -5;
xmax = 5;
ymin= -5;
ymax = 5;
px = 100; % Points in either direction
py = 100;      
figure('Name', 'Surface plot of f');
plotF( f, xmin, xmax, ymin, ymax, px, py );
end

%% Call the optimiser

% There is a problem in the lanczos algorithm
outputLevel = 2;
options = struct('outputLevel', outputLevel);
[MIN, iterates, gradients, k, RES] = ARC( f, g, H,  x0, options );



