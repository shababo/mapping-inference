addpath(genpath('../../psc-detection'),genpath('../../mapping-inference'),genpath('../../mapping-core'));
%% Simulate sums of a few Bernoulli random variables
n_cell=10;
gamma_truth=zeros(n_cell,1);
gamma_truth = (rand([n_cell 1])<0.5).*(0.5+0.5*rand([n_cell 1]));
% gamma_truth = (0.5+0.5*rand([n_cell 1]));
background_rate = 0.03;
%% Generate design matrix
n_trial= 2000;
designs = zeros(n_trial,n_cell);
outputs=zeros(n_trial,1);
for i_trial = 1:n_trial
    designs(i_trial,randsample(n_cell,3)) = 1;
    outputs(i_trial)=sum( rand(n_cell,1)< (designs(i_trial,:)'.*gamma_truth));
    outputs(i_trial)=outputs(i_trial)+ (rand(1,1)< background_rate);
end

%% Set the parameters 
% Initialize the variational family (spike-and-slab with beta distribution)
variational_params=struct([]);
for i_cell = 1:n_cell
    variational_params(i_cell).pi = 0.01;
    variational_params(i_cell).p_logit = log(variational_params(i_cell).pi/(1-variational_params(i_cell).pi));
    variational_params(i_cell).log_alpha = 0; %log_alpha
    variational_params(i_cell).log_beta = 0;%log_alpha
end

prioir_params.pi0= 0.7*ones(n_cell,1);
prioir_params.alpha0= ones(n_cell,1);
prioir_params.beta0 = ones(n_cell,1);
C_threshold = 0.01;maxit=1000;
S=200;epsilon=0.01;eta_logit=0;eta_beta=0.01;
%%
[parameter_history,change_history] = fit_working_model_vi(...
    designs,outputs, background_rate,...
    variational_params,prioir_params,C_threshold,...
S,epsilon,eta_logit,eta_beta,maxit);
%% Calculate the variational means 
last_iter = size(parameter_history.pi,2);

mean_gamma= (1-parameter_history.pi(:,last_iter)).*...
     (C_threshold+ (1-C_threshold)./(1+parameter_history.beta(:,last_iter)./parameter_history.alpha(:,last_iter)));
   
mean_gamma(gamma_truth==0)
mean_gamma(gamma_truth>0)
gamma_truth(gamma_truth>0)


