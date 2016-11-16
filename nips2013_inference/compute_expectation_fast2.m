function [E, values] = compute_expectation_fast2(num_samples, obj_function, alpha, mu, s_squared, sigma_n, sigma_s, hyperparam_alpha, X_t_plus_1, y_t,eta)
%COMPUTE_EXPECTATION(obj_function, alpha, mu, s_squared, sigma_n, sigma_s, hyperparam_alpha, X_t_plus_1, y_t)
%
% brute-force compute expectation of @obj_function over the distribution
% p(y_{t+1}|X_{1:t+1},y_{1:t})
%

K = size(X_t_plus_1,2);

% set varbvs hyperparameters
hyper_sigma_n_guess = sigma_n^2;
hyper_sigma_s_guess = sigma_s.^2 ./ sigma_n.^2;
hyper_logodds = logit(hyperparam_alpha);
   
% for now. in theory this function could take hyperparameter ranges.
assert(length(hyper_sigma_n_guess) == 1);
assert(length(hyper_sigma_s_guess) == 1);
assert(length(hyper_logodds) == 1);

options = struct('verbose', false);



% store 
values = zeros(num_samples, 1);

% we dont need to cmpute this everytime
% [X_center, dummy] = center_data(X_t_plus_1, zeros(length(y_t)+1,1));
% XtX = X_center'*X_center;
% XtX_alpha_mu = bsxfun(@times,bsxfun(@times,XtX,alpha),mu);


% Generate data from the fitted model:


for s = 1:num_samples
    % Draw the response status based on the probability 
    
    X_status = X_t_plus_1>rand(size(X_t_plus_1,1),K);
    
    % ancestral sample y_next:

    % 1. sample gamma
    gamma = rand(K,1) < alpha;

    % 2. sample weights
    w = normrnd(mu,sqrt(s_squared)).*gamma;

    % sample y_next
    % TODO note that this depends on sigma_n! when sigma_n is unknown, we
    % will have to revisit.
    y_next = X_status(end,:)*w  + normrnd(0,sigma_n);

    % update y vector
    y_test = [y_t; y_next];
    
    % center y
    %[dummy, y_center] = center_data(0, y_test);
    
    
    % can we get away with single updates!?
    
%     Xty = X_center'*y_center;
%     sig_update = hyper_sigma_n_guess./(diag(XtX) + 1/hyper_sigma_s_guess);
%     mu_update = (1/hyper_sigma_n_guess) * sig_update .* (Xty - sum(XtX_alpha_mu,2) + diag(XtX_alpha_mu));
%     ratio = hyperparam_alpha/(1 - hyperparam_alpha) * sqrt(sig_update) ./ (sigma_n*sigma_s) .* exp(0.5*mu_update.^2./sig_update);
%     alpha_update = ratio./(1 + ratio);
    
    % below is in case single updates are too slow, we can instead just
    % update THIS neuron
%     sig_update = hyper_sigma_n_guess/(X_center(:,k)'*X_center(:,k) + 1/hyper_sigma_s_guess);
%     mu_update = sig_update / hyper_sigma_n_guess * (X_center(:,k)'*y_center - ...
%         sum(X_center'*X_center(:,k).*alpha.*mu) + X_center(:,k)'*X_center(:,k)*alpha(k)*mu(k));
%     ratio = hyperparam_alpha/(1 - hyperparam_alpha) * sqrt(sig_update) / (sigma_n*sigma_s) * exp(0.5*mu_update^2/sig_update);
%     alpha_update = ratio/(1+ratio);

	% run vbs
    [~, alpha_update, mu_update, sig_update] = ...
        varbvs_general(X_status, y_test, hyper_sigma_n_guess, hyper_sigma_s_guess, hyper_logodds, eta,options);

    % store output
    values(s) = obj_function(alpha_update, mu_update, sig_update);

end

% take expectation
E = mean(values);
