
function [log_likelihood] = glm_lif_loglikelihood_oneloc(fit_params, responses, covs, params)

% repsonses is N x T
% covs is T x 3 x N
% fit_params is 3 x 1
% params is a struct
% params.dt
% params.link - link function
% params.dlink - derivative of link

[N, T] = size(responses);
dt = 1;% params.dt;

total_filter_out = zeros(size(responses));


for i = 1:N
    
    total_filter_out(i,:) = (fit_params*squeeze(covs(i,:,:)));
    
end

lambda = dt*params.invlink(total_filter_out);
d_lambda = dt*params.dinvlink(total_filter_out)*1e5;
assignin('base','d_lambda',d_lambda)

log_likelihood = -sum(sum(total_filter_out.*responses - lambda));
 

% grad_ll = zeros(size(fit_params));
% for i = 1:length(fit_params)
%     
%     this_cov = squeeze(covs(:,i,:));
% %     size(responses)
% %     size(lambda)
% %     size(d_lambda)
% %     size(this_cov)
%     grad_ll(i) = ...
%         -sum(sum(responses./lambda.*d_lambda.*this_cov - d_lambda.*this_cov));
% 
% end
% 
