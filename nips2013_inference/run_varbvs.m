function [alpha_rnd, mu_rnd, s_sq_rnd] = run_varbvs(X, Y, sigma_n, sigma_s, alpha, options)
%RUN_VARBVS
%
% centers the data and runs varbvs "inner loop" with specified 
% hyperparameters.
%

hyper_sigma_n_guess = sigma_n^2;
hyper_sigma_s_guess = sigma_s.^2 ./ sigma_n.^2;
hyper_logodds = logit(alpha);

% for now. in theory this function could take hyperparameter ranges.
assert(length(hyper_sigma_n_guess) == 1);
assert(length(hyper_sigma_s_guess) == 1);
assert(length(hyper_logodds) == 1);

if nargin < 6
    options = struct('verbose', true);
end

% CENTER Y and X:
[X_center, Y_center] = center_data(X, Y);

% RUN VARBVS
[~, alpha_rnd, mu_rnd, s_sq_rnd] = varbvs(X_center,Y_center,hyper_sigma_n_guess,hyper_sigma_s_guess,hyper_logodds,options);
    