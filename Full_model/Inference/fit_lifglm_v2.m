function [stats_conv] = fit_lifglm_v2(responses,stims,in_params)
% responses is an N x T binary matrix of spike times where we have N trials
% stims is an N x T matrix of the stimulus timeseries which is scaled by power - if using a shape to
% fit the lif-glm then these should be scaled by the shape as well.
% stims_ind is an N x 1 vector of location indices for each stimulus, this is if you are fitting
% each location independently, if using ashape set all to 1
% in_params is struct with the field g which is 1/tau for this neuron

% in_params.link_type
% in_params.fit_init
% in_params.g

% [num_trials, num_samps] = size(responses);

[cov_resting, cov_reset, cov_current]  = gconv_v3(stims',responses',in_params.g);
% covs = permute(covs,[3 2 1]);
%assignin('base','covs_tmp',covs)

% params.dt = 1/20000;

% params.link = @(mu) log(exp(mu)-1);  %link = @(mu) mu + log(1-exp(-mu));
% params.dlink = @(mu) exp(mu)./(exp(mu)-1);
% params.dinvlink = @(resp) exp(resp)./(1 + exp(resp));
% params.invlink = @(resp) log(1 + exp(resp));

params.invlink = @invlink_test;
params.dlink = @derlink_test;
params.link = @link_test;
params.dinvlink = @derlink_test;

% params.link = @(mu) log(mu);  %link = @(mu) mu + log(1-exp(-mu));
% params.dlink = @(mu) 1./mu;
% params.dinvlink = @(resp) exp(resp);
% params.invlink = @(resp) exp(resp);
F = {params.link, params.dlink, params.invlink};
%  = @(x) exp(x);
%  = @(x) exp(x);

covs_1trial = zeros(length(responses(:)),3); 
covs_1trial(:,1) = cov_reset(:);
covs_1trial(:,2) = cov_current(:);
covs_1trial(:,3) = cov_resting(:); % not fitting resting 
r_temp = responses';
% this sets v_th at 15

%[betahat_glm,dev,stats_conv]=glmfit(covs_1trial(:,1:2),r_temp(:),...
%    'poisson','link',F,'constant','off','offset',-15*ones(length(responses(:)),1));



[betahat_glm,dev,stats_conv]=glmfit(covs_1trial(:,2),r_temp(:),...
    'poisson','link',F,'constant','off','offset',-4e3*covs_1trial(:,1)-15*ones(length(responses(:)),1));


% this fits v_th
% [betahat_glm,dev,stats_conv]=glmfit(covs_1trial(:,[1 2 3:end]),responses(:),...
%     'poisson','link',F,'constant','off');%

stats_conv.dev = dev;
% betahat_glm = []; stats_conv = [];
% [betahat_fmin, ll_fmin] = fmincon(obj_fun,x0,[],[],[],[],lb,ub,[],options);
betahat_fmin = []; ll_fmin = 0;