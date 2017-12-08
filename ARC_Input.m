function [epsilon, IterMax, sigma0, theta, eta1, eta2, outputLevel] = ARC_Input( options )
%ARC_Input
%   parses the input arguments for main ARC algorithm


if isfield(options, 'epsilon')
    epsilon = options.epsilon;
else
    epsilon = 1e-5;
end

if isfield(options, 'IterMax')
    IterMax = options.IterMax;
else
    IterMax = 100;
end

if isfield(options, 'sigma0')
    sigma0 = options.sigma0;
else
    sigma0 = 1;
end

if isfield(options, 'theta')
    theta = options.theta;
else
    theta = 1e-4;
end

if isfield(options, 'eta1')
    eta1 = options.eta1;
else
    eta1 = 0.1;
end

if isfield(options, 'eta2')
    eta2 = options.eta2;
else
    eta2 = 0.9;
end

if isfield(options, 'outputLevel')
    outputLevel = options.outputLevel;
else
    outputLevel = 2;
end

end

