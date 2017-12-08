function [theta, sigma0, outputLevel] = GLRT_Input( options )
% GLRT_Input
%   Input handeling for GLRT

if isfield(options, 'theta')
    theta = options.theta;
else
    theta = 1e-4;
end

if isfield(options, 'sigma0')
    sigma0 = options.sigma0;
else
    sigma0 = 1;
end

if isfield(options, 'outputLevel')
    outputLevel = options.outputLevel;
else
    outputLevel = 2;
end

end

