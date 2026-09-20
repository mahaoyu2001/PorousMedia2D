function [mu, variance, skewnessValue, kurtosisValue] = ...
    cal_stress_statistics(Sigma, porosity)
%CAL_STRESS_STATISTICS Evaluate the analytical stress moments.
%   SIGMA is [Sigma11; Sigma22; Sigma12]. POROSITY may be a scalar or row
%   vector. The returned rows correspond to the three in-plane components.

validateattributes(Sigma, {'numeric'}, {'vector', 'numel', 3, 'finite'});
validateattributes(porosity, {'numeric'}, ...
    {'real', 'nonnegative', '<', 1, 'finite'});

[coefficient, ~, ~] = cal_stress_coef(Sigma(1), Sigma(2));
k2 = coefficient(:, 2);
k3 = coefficient(:, 3);
k4 = coefficient(:, 4);

mu = Sigma(:) ./ (1 - porosity);
variance = k2 .* porosity ./ (1 - porosity).^3;
skewnessValue = k3 .* sqrt((1 - porosity) ./ porosity);
kurtosisValue = k4 ./ k2.^2 .* (1 - porosity) ./ porosity ...
    + 3 * (1 - porosity);
end
