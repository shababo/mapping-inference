function [alpha, mu, s_sq] = run_varbvs_general(X, Y, sigma_n, sigma_s, alpha, eta, options)
%RUN_VARBVS
%
% centers the data and runs varbvs "inner loop" with specified 
% hyperparameters.
%

hyper_sigma_n_guess = sigma_n^2;
hyper_sigma_s_guess = sigma_s.^2;
hyper_logodds = logit(alpha);


if ~exist('options')
        options = struct('verbose', true);
end

    if ~isfield(options,'center')
        options.center = 0;
    end

% CENTER Y and X:
if options.center
    [X_center, Y_center] = center_data(X, Y);
else
    X_center = X;
    Y_center = Y;
end

% RUN VARBVS
[~, alpha, mu, s_sq] = varbvs_general(X_center,Y_center,hyper_sigma_n_guess,hyper_sigma_s_guess,hyper_logodds,eta,options);
