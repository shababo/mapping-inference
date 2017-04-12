function betahat = fit_lifglm_oneloc(responses,stims,stims_ind,in_params)

% in_params.link_type
% in_params.fit_init
% in_params.g

% [num_trials, num_samps] = size(responses);

covs = gconv(stims',stims_ind,responses',in_params.g);
covs = permute(covs,[3 2 1]);
params.dt = 1/20000;
params.link = @(x) exp(x);
params.dlink = @(x) exp(x);
x0 = [0 0 .0001*ones(1,max(stims_ind))];
obj_fun = @(x) glm_lif_loglikelihood_oneloc(x, responses, covs, params);
options = optimoptions('fmincon','Algorithm','interior-point',...
                                'Display','iter','GradObj','on',...
                                'Diagnostics','on',...
                                'UseParallel',true); 
                            
ub = Inf*ones(size(x0));                            
lb = -Inf*ones(size(x0));
lb(3:end) = 0;


betahat = fmincon(obj_fun,x0,[],[],[],[],lb,ub,[],options);

