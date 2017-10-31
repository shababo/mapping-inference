%% Parameters in the design stage
% Design parameters

% Initialize the five cell groups
undefined_cells= cell(0); undefined_cells{1}=ones(n_cell_this_plane,1);%A
disconnected_cells= cell(0); disconnected_cells{1}=zeros(n_cell_this_plane,1);%B
dead_cells= cell(0); dead_cells{1}=zeros(n_cell_this_plane,1);%D
connected_cells= cell(0); connected_cells{1}=zeros(n_cell_this_plane,1);%C
alive_cells= cell(0);alive_cells{1}=zeros(n_cell_this_plane,1);%E

% Initialize the trial records 
mpp_undefined=struct([]);mpp_connected=struct([]);

%%
related_cell_list=[target_cell_list.primary target_cell_list.secondary];
n_related_cell=length(related_cell_list);
%% Initialize the posterior/variational distributions
switch prior_info_type
    case 1 
        var_gain_prior=0;
        gain_bias=0;
        disp('Use good prior')
    case 2 %   
        var_gain_prior=0; 
        gain_bias=0;
        disp('Use uninformative prior')
    case 3 %
        var_gain_prior=0;
        gain_bias=0.005;
        disp('Biased prior')
    otherwise
end

prior_pi0=0.5;var_pi_ini=0.5;
var_alpha_initial=0;var_beta_initial=1;
perturbed_gain=max(0.0055,gain_truth(related_cell_list)+gain_bias*(2*(rand([n_related_cell, 1])-0.5) ));
var_alpha_gain_initial=log((perturbed_gain  - gain_bound.low)./(gain_bound.up-perturbed_gain));
var_beta_gain_initial=var_gain_prior; % uncertainty of the prior

prior_params=struct;
temp=num2cell(log(prior_pi0/(1-prior_pi0))*ones(n_related_cell,1));[prior_params(1:n_related_cell).p_logit]=temp{:};
temp=num2cell(var_alpha_initial*ones(n_related_cell,1));[prior_params(1:n_related_cell).alpha]=temp{:};
temp=num2cell(var_beta_initial*ones(n_related_cell,1));[prior_params(1:n_related_cell).beta]=temp{:};
temp=num2cell(var_alpha_gain_initial);[prior_params(1:n_related_cell).alpha_gain]=temp{:};
temp=num2cell(var_beta_gain_initial*ones(n_related_cell,1));[prior_params(1:n_related_cell).beta_gain]=temp{:};
%%
% The prioir info serves as the first variational distributions 
variational_params_path=struct;
variational_params_path=prior_params;

% Initialize storage for the fitted parameters in the experiment
[parameter_path] = calculate_posterior(prior_params,gamma_bound,gain_bound,quantiles_prob,true);
% 
% parameter_path=struct([]);
% temp=num2cell(0.5*ones(length(related_cell_list),1));[parameter_path(1,1:n_related_cell).gamma_mean]=temp{:};
% temp=num2cell(zeros(length(related_cell_list),1));[parameter_path(1,1:n_related_cell).gamma_upper_quantile]=temp{:};
% temp=num2cell(zeros(length(related_cell_list),1));[parameter_path(1,1:n_related_cell).gamma_lower_quantile]=temp{:};
% temp=num2cell(zeros(length(related_cell_list),1));[parameter_path(1,1:n_related_cell).gain_mean]=temp{:};
% temp=num2cell(zeros(length(related_cell_list),1));[parameter_path(1:n_related_cell).alpha_gain]=temp{:};
% temp=num2cell(zeros(length(related_cell_list),1));[parameter_path(1,1:n_related_cell).gain_var]=temp{:};
