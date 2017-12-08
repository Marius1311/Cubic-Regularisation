function [  ] = plotF( f, xmin, xmax, ymin, ymax, px, py, varargin )
%This simply plots f in case there are only two input arguments

v = linspace(xmin, xmax, px);
w = linspace(ymin, ymax, py);
[X, Y] = meshgrid(v, w);
P = reshape(X, [px * py, 1]);
Q = reshape(Y, [px * py, 1]);
RawZ = evalF(f, [P'; Q']);
Z = reshape(RawZ', [py, px]);

if nargin == 7
    mesh(X, Y, Z), xlabel('x'), ylabel('y'), zlabel('z'), title('Surface plot of the objective function f');
else
    x = varargin{1};
    mesh(X, Y, Z), xlabel('x'), ylabel('y'), zlabel('z'), title('Surface plot of the objective function f');
    hold on;
    plot3(x(:, 1), x(:, 2), evalF(f, x'), 'ro-', 'MarkerSize', 3, 'LineWidth', 2 );
    hold off;

end

