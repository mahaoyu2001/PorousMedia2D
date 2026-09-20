function [coefficient, correlation, correlation2] = cal_stress_coef(sigma11, sigma22)
%CAL_STRESS_COEF Substitute the remote stress into symbolic coefficients.
%   The bundled MAT file stores the analytical expressions used in the
%   paper. This function requires MATLAB Symbolic Math Toolbox.

functionDir = fileparts(mfilename('fullpath'));
dataFile = fullfile(functionDir, 'stress_stats_symbolic.mat');
if ~isfile(dataFile)
    error('Analytical coefficient file not found: %s', dataFile);
end

data = load(dataFile, 'k_n', 'corr', 'corr2', 'Sigma_1', 'Sigma_2');
coefficient = double(subs(data.k_n, ...
    {data.Sigma_1, data.Sigma_2}, {sigma11, sigma22}));
correlation = double(subs(data.corr, ...
    {data.Sigma_1, data.Sigma_2}, {sigma11, sigma22}));
correlation2 = double(subs(data.corr2, ...
    {data.Sigma_1, data.Sigma_2}, {sigma11, sigma22}));
end
