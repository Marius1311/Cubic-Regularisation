function [ y ] = evalF(f, X )
% Evaluates an objective function f at a point X. 
% X can either be a column vecor or a matrix composed of column vectors
%
% Input:
% f: function handle to an objective function. SHALL ACCEPT ROW VECTORS
% X: Matrix of points to evaluate. Each column corresponds to one point,
% each row to one of the components
%
% Output:
% y: row vector of function values, corresponding to columns in X

% Initialise
[a, b] = size(X);

% Reshape, in case a single point shall be evaluated
if a == 1 || b == 1
   X = X(:); 
   b = 1;
end

y = zeros(1, b);

% Fill y
for j = 1:b
   y(1, j) = f(X(:, j)'); 
end

end

